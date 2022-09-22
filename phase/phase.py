import os
import logging
import itertools
from collections import OrderedDict

import numpy as np
import pandas as pd
import pysam
from scipy.stats import fisher_exact

from .sample import Sample
from .mutation import Mutation
from .mutation import ClonalStructure
from .writer import *


# make the ID variables float for handling these with 'nan'
PHASE_ID = 0.0
GERMLINE_ID = 0.0

MNP_TYPE = {1:'SNP', 2:'DNP', 3:'TNP'}


def check_input_args(args):
    """
    Check input parameters
    """
    # check required arguments
    if args.ref_fasta == None:
        raise Exception('provide reference FASTA')
    if not os.path.exists(args.ref_fasta):
        raise FileNotFoundError('reference FASTA not found: '+args.ref_fasta)
    
    if args.input_maf == None:
        raise Exception('provide input MAF')
    if not os.path.exists(args.input_maf):
        raise FileNotFoundError('input MAF not found: '+args.input_maf)
    
    if args.tumor_name == None:
        raise Exception('provide tumor sample name in the input MAF')

    if args.tumor_bam == None:
        raise Exception('provide tumor BAM')
    if not os.path.exists(args.tumor_bam):
        raise FileNotFoundError('tumor BAM not found: '+args.tumor_bam)
    
    # check optional arguments
    if args.tumor_vcf:
        if not os.path.exists(args.tumor_vcf):
            raise FileNotFoundError('tumor VCF not found')

    if args.normal_name or args.normal_bam or args.normal_vcf:
        if args.normal_name == None:
            raise Exception('provide normal sample name')
        if not (args.normal_bam and args.normal_vcf):
            raise Exception('provide both normal VCF and BAM')
        if args.normal_bam == None:
            raise Exception('normal VCF provided, but cannot find normal BAM')
        if args.normal_vcf == None:
            raise Exception('normal BAM provided, but cannot find normal VCF')
        if not os.path.exists(args.normal_vcf):
            raise FileNotFoundError('normal VCF not found: '+args.normal_vcf)
        if not os.path.exists(args.normal_bam):
            raise FileNotFoundError('normal BAM not found: '+args.normal_bam)

    if args.phylogic_ccfs or args.phylogic_tree:
        if args.phylogic_ccfs == None:
            raise Exception('phylogic tree provided, but cannot find the CCF file')
        if args.phylogic_tree == None:
            raise Exception('phylogic CCF provided, but cannot find the tree file')
        if not os.path.exists(args.phylogic_ccfs):
            raise FileNotFoundError('phylogic CCF not found: '+args.phylogic_ccfs)
        if not os.path.exists(args.phylogic_tree):
            raise FileNotFoundError('phylogic tree not found: '+args.phylogic_tree)
        if args.phylogic_index <= 0:
            raise ValueError('invalid tree index; provide a value larger than 0')

    if args.max_mnp_size <= 3:
        raise ValueError('invalid maximum MNP; provide a value larger than 3')
    if args.max_phase_radius <= args.max_mnp_size:
        raise ValueError('invalid max phasing radius; provide a value larger than max MNP length')
    if args.vaf_skew_pval <= 0 or 1 <= args.vaf_skew_pval:
        raise ValueError('invalid VAF skew p-value; a value nees to be 0 < p < 1')
    if args.min_map_quality < 0:
        raise ValueError('invalid min mapping quality; provide a value larger than or equal to 0')
    if args.min_base_quality < 0:
        raise ValueError('invalid min base quality; provide a value larger than or equal to 0')
    if args.min_coverage < 0:
        raise ValueError('invalid min coverage cutoff; provide a value larger than 0')
    if args.min_phased_altc < 0:
        raise ValueError('invalid min phased alt counts; provide a value larger than 0')
    
    
def group_mutations(df, max_dist, group_id, min_cov, excl_indel=False):
    """
    Group mutations within a certain distance, but exclude sites
    below a coverage cutoff from assigning group IDs
    Grouped mutations will have the same numerical IDs (start=1)
    """
    df[group_id] = np.full(len(df.index), np.nan)
    df.sort_values('Start_position', ignore_index=True, inplace=True)
    curr_id = 0
    prev_id = 0
    for i, dist in list(df['Start_position'].diff().items())[1:]:
        is_indel = sum(df.loc[i-1:i]['Variant_Type'].isin(['DEL', 'INS']))
        cov = min(df.loc[i-1:i]['t_alt_count'] + df.loc[i-1:i]['t_ref_count'])
        if (excl_indel and is_indel) or (cov < min_cov):
            curr_id += 1
            continue
        if dist <= max_dist:
            curr_id = prev_id+1 if curr_id > prev_id else curr_id
            df.at[i-1, group_id] = curr_id
            df.at[i, group_id] = curr_id
            prev_id = curr_id
        else:
            curr_id += 1
    return df


def get_read_names(x, g, bam, min_base_q, min_map_q):
    r = {'alt': set(), 'ref': set()}
    if bam == None:
        return r
    ref = x['Reference_Allele']
    alt = x['Tumor_Seq_Allele']
    mut_type = x['Variant_Type']
    chrom = x['Chromosome']
    pos = x['Start_position']-(0 if mut_type == 'INS' else 1)
    end = x['Start_position']+(1 if mut_type == 'INS' else len(alt)-1)
    for pileupcolumn in bam.pileup(chrom, pos, end,
                                   truncate=True,
                                   max_depth=1000000,
                                   min_base_quality=min_base_q,
                                   min_mapping_quality=min_map_q):
        for read in pileupcolumn.pileups:
            # filter reads with the following flags
            if read.alignment.is_duplicate or \
               read.alignment.is_qcfail or \
               read.alignment.is_secondary or \
               read.alignment.is_unmapped or \
               read.alignment.is_supplementary:
                continue
            rg = None
            base = None
            if read.query_position:
                base = read.alignment.seq[read.query_position]

            if mut_type in ['INS', 'DEL']:
                blocks = read.alignment.get_blocks()
                if len(blocks) > 1:
                    prev_block_end = blocks[0][1]
                    dist_to_mut = abs(prev_block_end - pos)
                    for block_start, block_end in blocks[1:]:
                        ref_space = block_start - prev_block_end
                        dist_to_mut = min(dist_to_mut, abs(block_start - pos))
                        # insertion
                        if ref_space == 0 and \
                           mut_type == 'INS' and \
                           pos == prev_block_end:
                            rg = 'alt'
                            break
                        # deletion
                        elif ref_space == len(ref) and \
                             mut_type == 'DEL' and \
                             prev_block_end <= pos <= block_start:
                            rg = 'alt'
                            break
                        prev_block_end = block_end
                    if rg == None and \
                       dist_to_mut > 1 and \
                       base == g.fetch(chrom, pos, end):
                        rg = 'ref'
                else:
                    if base == g.fetch(chrom, pos, end):
                        rg = 'ref'
            else: # point mutation
                if base == None:
                    continue
                offset = pileupcolumn.reference_pos - pos
                alt_base = alt[offset]
                ref_base = ref[offset]
                if base == alt_base:
                    rg = 'alt'
                elif base == ref_base:
                    rg = 'ref'
            try:
                r[rg].add(read.alignment.query_name)
            except KeyError:
                continue
    return r


def vaf_from_reads(r):
    alt = len(r['alt'])
    ref = len(r['ref'])
    try:
        vaf = alt/float(alt + ref)
        return alt, ref, vaf
    except ZeroDivisionError:
        return 0, 0, 0


def calc_phase_prob(r1, r2, min_phased_altc):
    r1_alt, r1_ref, r1_vaf = vaf_from_reads(r1)
    r2_alt, r2_ref, r2_vaf = vaf_from_reads(r2)
    # check if VAFs of the two mutations are the same
    vaf_odds, vaf_pval = fisher_exact([[r1_alt, r1_ref],
                                       [r2_alt, r2_ref]])
    # raw alt and ref counts from a mutation with smaller VAF
    if r1_vaf > r2_vaf:
        u_alt = r2_alt
        u_ref = r2_ref
    else:
        u_alt = r1_alt
        u_ref = r1_ref
    # phased alt and ref counts
    p_alt = len(r1['alt'].intersection(r2['alt']))
    p_ref = len(r1['ref'].intersection(r2['ref']))
    phase_odds, phase_pval = fisher_exact([[u_alt, u_ref],
                                           [p_alt, p_ref]])
    if p_alt < min_phased_altc:
        phase_pval = 0.0
    if p_alt > min_phased_altc * 4:
        phase_pval = 1.0
    return (phase_pval, vaf_pval)


def cluster_mutations(pos, p, cutoff):
    # p-value is insignificant when variants are phased
    if len(pos) == 1:
        return set([(pos[0],)])
    gi = dict()
    for x, y in itertools.combinations(pos, 2):
        hx, hy = (x,), (y,)
        if p[(x, y)] > cutoff:
            hx = hy = (x, y)
        gi[x] = gi.get(x,set()).union(hx)
        gi[y] = gi.get(y,set()).union(hy)
    gi = set([tuple(sorted(h)) for h in gi.values()])
    return gi


def chain_mnps(cl, pos, df, mnp_dict, g, chrom, t_reads, n_reads):
    start, end = list(df.loc[[pos[0],pos[-1]],'Start_position'])
    posr = list(df.loc[pos,'Start_position'] - start)
    maf_idx = set()
    ref = list(g.fetch(chrom, start-1, end))
    alt = np.array(ref, dtype=str)
    t_alt_count, t_ref_count = [], []
    n_alt_count, n_ref_count = [], []
    for i in range(len(pos)):
        maf_idx = maf_idx.union(df.loc[pos[i], 'maf_idx'])
        # this can handle existing DNP, TNP, ONP
        alt_base = df.loc[pos[i], 'Tumor_Seq_Allele']
        alt[posr[i]:posr[i]+len(alt_base)] = list(alt_base)
        t_alt_count.append(t_reads[pos[i]]['alt'])
        t_ref_count.append(t_reads[pos[i]]['ref'])
        n_alt_count.append(n_reads[pos[i]]['alt'])
        n_ref_count.append(n_reads[pos[i]]['ref'])
    row = OrderedDict(df.loc[pos[0]])
    row['maf_idx'] = maf_idx
    row['ClonalStructure'] = set([cl])
    row['Reference_Allele'] = ''.join(ref)
    row['Tumor_Seq_Allele'] = ''.join(alt)
    row['Variant_Type'] = MNP_TYPE.get(len(alt), 'MNP')
    key = tuple([row[k] for k in list(row.keys())[:5]])
    if key in mnp_dict:
        mnp_dict[key]['ClonalStructure'].add(cl)
        return mnp_dict
    if row['Variant_Type'] != 'SNP':
        row['t_alt_count'] = len(set.intersection(*t_alt_count))
        row['t_ref_count'] = len(set.intersection(*t_ref_count))
        if 'n_alt_count' in df:
            row['n_alt_count'] = len(set.intersection(*n_alt_count))
            row['n_ref_count'] = len(set.intersection(*n_ref_count))
    mnp_dict.setdefault(key, row)
    return mnp_dict


def phase_mutations(df, chrom,
                    t, n, genome,
                    clonal_struct,
                    vaf_skew_pval,
                    min_phased_altc,
                    min_base_q,
                    min_map_q):
    """
    Look at reads supporting candidate mutations to determine phasing
    """
    global PHASE_ID
    df = df.reset_index(drop=True)
    g = pysam.FastaFile(genome)
    tbam = pysam.AlignmentFile(t.bam, 'rb')
    nbam = pysam.AlignmentFile(n.bam, 'rb') if n.bam else None
    # get alt and ref reads from tumor and normal samples
    t_reads = []
    n_reads = []
    for i, row in df.iterrows():
        t_reads.append(get_read_names(row, g, tbam, min_base_q, min_map_q))
        n_reads.append(get_read_names(row, g, nbam, min_base_q, min_map_q))
    # compute p-values: [0] phased reads, [1] allele fractions
    unp = 0
    phase_p = dict()
    vaf_p = dict()
    for x, y in itertools.combinations(df.index, 2):
        phase_p[(x,y)], vaf_p[(x,y)] = calc_phase_prob(t_reads[x],
                                                       t_reads[y],
                                                       min_phased_altc)
        if phase_p[(x,y)] <= vaf_skew_pval:
            unp += 1
    # return the original df if not phased
    if unp == len(phase_p):
        df[['PhaseID', 'MNP']] = np.nan
        return df
    
    mutcl = pd.Series([np.nan]*len(df.index))
    if clonal_struct:
        mutcl = df['Cluster_Assignment']
        # switch to VAF estimate if undefined clone clusters exist
        if clonal_struct.excluded_clone.intersection(set(mutcl.unique())):
            mutcl = pd.Series([np.nan]*len(df.index))
    # estimate clonal structure from VAF if the information is not available
    if sum(mutcl.isna()) > 0:
        vaf_tiers = cluster_mutations(df.index, vaf_p, vaf_skew_pval)
        vaf_means = dict()
        for vt in vaf_tiers:
            alt = df.loc[list(vt),'t_alt_count']
            ref = df.loc[list(vt),'t_ref_count']
            vaf_means.setdefault(vt, np.mean(alt/(alt+ref)))
        rank = sorted(vaf_means.items(),
                      key=lambda kv: kv[1], reverse=True)
        input_tree = []
        for i, vt in enumerate(rank):
            mutcl[list(vt[0])] = str(i+1)
            input_tree.append(str(i+1) + '-' + str(i+2))
        input_tree = ','.join(input_tree)
        clonal_struct = ClonalStructure(input_tree, 0, set(mutcl.unique()))
    # identify phased mutations per clone
    clones = dict()
    for cl in mutcl.unique():
        pos = mutcl[mutcl == cl].index
        for tier in cluster_mutations(pos, phase_p, vaf_skew_pval):
            clones.setdefault(tier, [cl, set(tier), set(tier)])
    # generate phased mutations
    for tx, ty in itertools.combinations(clones.keys(), 2):
        cx = clones[tx][0]
        cy = clones[ty][0]
        if cx == cy:
            continue
        if not clonal_struct.is_related(cx, cy):
            continue
        tier = cluster_mutations(list(tx)+list(ty),
                                 phase_p, vaf_skew_pval)
        if not tier.difference(set([tx, ty])):
            continue
        txy = set(list(tx) + list(ty))
        if clonal_struct.tree[cx].depth > clonal_struct.tree[cy].depth:
            clones[tx][1] = clones[tx][1].union(txy)
        else:
            clones[ty][1] = clones[ty][1].union(txy)
        clones[tx][2] = clones[tx][2].union(txy)
        clones[ty][2] = clones[ty][2].union(txy)

    # assign phase ID and resolve clonal structure
    df['ClonalStructure'] = [set([i]) for i in mutcl]
    df['PhaseID'] = np.nan
    blocks = set()
    for cl, tier, block in clones.values():
        if len(tier) == 1:
            continue
        df.loc[list(tier), 'ClonalStructure'].apply(lambda x: x.add(cl))
        if tier == block:
            PHASE_ID += 1
            df.loc[list(block), 'PhaseID'] = PHASE_ID

    is_unphased = df['PhaseID'].isna()
    df.at[is_unphased, 'MNP'] = np.nan
    df.at[is_unphased, 'ClonalStructure'] = [set() for i in range(sum(is_unphased))]
    if sum(df['MNP'].isna()) == len(df.index):
        return df
    
    # chain multi-nucleotide polymorphisms
    index = []
    mnp_dict = dict()
    for mnp_id in df['MNP'].dropna().unique():
        mnp_i = df['MNP'] == mnp_id
        mutcl = set.union(*df[mnp_i]['ClonalStructure'])
        for cl in mutcl:
            pos = df.index[mnp_i & df['ClonalStructure'].apply(lambda x: cl in x)]
            mnp_dict = chain_mnps(cl, pos, df, mnp_dict, g, chrom, t_reads, n_reads)
            index.extend(pos)
    index = list(set(index))
    df = df.drop(index)
    for key in mnp_dict:
        df = df.append(mnp_dict[key], ignore_index=True)
    if len(df.index) == 1:
        # reset phase ID
        PHASE_ID -= 1
        df['PhaseID'] = np.nan
    return df.reset_index(drop=True)


def integrate_germline(df, chrom, t, n, genome, max_dist, min_base_q, min_map_q):
    global GERMLINE_ID
    germline = []
    gm_size = len(germline)
    g = pysam.FastaFile(genome)
    vcf = pysam.VariantFile(n.vcf, 'r')
    tbam = pysam.AlignmentFile(t.bam, 'rb')
    for i, row in df.iterrows():
        GERMLINE_ID += 1
        alt = row['Tumor_Seq_Allele']
        mut_type = row['Variant_Type']
        pos = row['Start_position']-(0 if mut_type == 'INS' else 1)
        end = row['Start_position']+(1 if mut_type == 'INS' else len(alt)-1)
        
        var_list = [v for v in vcf.fetch(chrom, pos-max_dist, end+max_dist)]
        if len(var_list) == 0:
            GERMLINE_ID -= 1
            continue
        
        t_reads = get_read_names(row, g, tbam, min_base_q, min_map_q)['alt']
        for v in var_list:
            if not 'PASS' in v.filter:
                continue
            allele_depth = v.samples[n.vcf_sm]['AD']
            for j, gm_alt in enumerate(v.alts):
                # skip spanning deletion
                if gm_alt == '*':
                    continue
                if len(v.ref) == len(gm_alt):
                    gm_mut_type = MNP_TYPE.get(len(alt), 'MNP')
                elif len(v.ref) > len(gm_alt):
                    gm_mut_type = 'DEL'
                elif len(v.ref) < len(gm_alt):
                    gm_mut_type = 'INS'
                gm = {
                    'Chromosome': v.chrom,
                    'Start_position': v.pos,
                    'Reference_Allele': v.ref,
                    'Tumor_Seq_Allele': gm_alt,
                    'Variant_Type': gm_mut_type,
                    't_ref_count': 0,
                    't_alt_count': 0,
                    'n_ref_count': allele_depth[0],
                    'n_alt_count': allele_depth[j+1],
                    'genotype': '/'.join(map(str,v.samples[n.vcf_sm]['GT'])),
                    'GermlineID': GERMLINE_ID
                }
                n_reads = get_read_names(gm, g, tbam, min_base_q, min_map_q)['alt']
                if set.intersection(*[t_reads, n_reads]):
                    if gm not in germline:
                        germline.append(gm)
        if len(germline) > gm_size:
            df.at[i, 'GermlineID'] = GERMLINE_ID
            gm_size = len(germline)
        else:
            GERMLINE_ID -= 1
    return df, pd.DataFrame(germline)


def run(args):    
    check_input_args(args)

    logging.basicConfig(level=logging.INFO,
                        format='phase | %(asctime)s [%(levelname)s] %(message)s',
                        datefmt='%Y-%m-%d %H:%M:%S')
    
    t = Sample(args.tumor_bam, args.tumor_vcf)
    n = Sample(args.normal_bam, args.normal_vcf)

    logging.info('loading mutations from input MAF')
    m = Mutation(args.input_maf,
                 args.tumor_name,
                 t.bam_sm,
                 args.phylogic_ccfs,
                 args.phylogic_tree,
                 args.phylogic_index)

    if 'n_alt_count' in m.muts and n.bam == None:
        logging.warn('normal BAM not supplied; normal allele depths for MNP will be zero')
    
    smuts = []
    gmuts = []
    m.muts['maf_idx'] = [set([i]) for i in m.muts.index]
    
    # use the original chromosome order (i.e. sort = False)
    for chrom, df in m.muts.groupby('Chromosome', sort=False):
        logging.info('processing chromosome {}'.format(chrom))
        
        df = group_mutations(df, args.max_phase_radius,
                             'PhaseID', args.min_coverage)
        df = group_mutations(df, args.max_mnp_size,
                             'MNP', args.min_coverage,
                             excl_indel=True)
        df['ClonalStructure'] = [set() for i in df.index]
        
        dfs = [df[df['PhaseID'].isna()]]
        if len(dfs[0].index) < len(df.index):
            for i, x in df.groupby('PhaseID', dropna=True, sort=False):
                dfs.append(
                    phase_mutations(
                        x, chrom, t, n,
                        args.ref_fasta,
                        m.clonal_struct,
                        args.vaf_skew_pval,
                        args.min_phased_altc,
                        args.min_base_quality,
                        args.min_map_quality
                    )
                )
        dfs = pd.concat(dfs)
        dfs.sort_values('Start_position', ignore_index=True, inplace=True)
        
        dfs['ClonalStructure'] = dfs['ClonalStructure'].apply(
            lambda x: ','.join(sorted(x)) if len(x) else np.nan
        )
        dfs.drop(columns='MNP', inplace=True)

        dfs['GermlineID'] = np.nan
        if n.vcf:
            dfs, dfg = integrate_germline(dfs, chrom, t, n,
                                          args.ref_fasta,
                                          args.max_phase_radius,
                                          args.min_base_quality,
                                          args.min_map_quality)
            gmuts.append(dfg)

        smuts.append(dfs)

    # extract phased mutations
    smuts = pd.concat(smuts).reset_index(drop=True)
    pi = ~smuts['ClonalStructure'].isna() | ~smuts['GermlineID'].isna()

    # write unphased mutations as a subset of the original MAF
    logging.info('writing {}.unphased.maf'.format(args.tumor_name))
    write_subset_maf(args.tumor_name+'.unphased.maf',
                     [x.pop() for x in smuts[~pi]['maf_idx']],
                     m.maf)
    
    # take minimum CCF for chained MNPs
    cols_to_keep = ['ClonalStructure','PhaseID','GermlineID']
    is_mnp = pi & smuts['Variant_Type'].isin(['DNP','TNP','MNP'])
    if 'ccf_hat' in smuts and sum(is_mnp) > 0:
        smuts.at[is_mnp,'ccf_hat'] = smuts[is_mnp]['maf_idx'].apply(
            lambda x: m.muts.loc[list(x), 'ccf_hat'].min()
        )
        cols_to_keep.insert(0,'ccf_hat')
    
    # write somatic phased mutations
    if sum(pi) > 0:
        logging.info('found {} somatic phased mutations'.format(sum(pi)))
        smuts = smuts[pi].reset_index(drop=True)
        count = list(smuts.columns[smuts.columns.str.endswith('_count')])
        smuts[count] = smuts[count].astype(int)
        logging.info('writing {}.phased.vcf'.format(args.tumor_name))
        if t.vcf:
            write_phase_vcf(args.tumor_name+'.phased.vcf',
                            smuts, t, cols_to_keep)
        else:
            write_phase_vcf_from_scratch(args.tumor_name+'.phased.vcf',
                                         t.bam_sm, n.bam_sm, smuts)
        logging.info('writing {}.phased.maflite.tsv'.format(args.tumor_name))
        write_phase_maflite(args.tumor_name+'.phased.maflite.tsv',
                            smuts, cols_to_keep)

        # write germline variants phased with somatic variants
        if gmuts:
            gmuts = pd.concat(gmuts).reset_index(drop=True)
            logging.info('found {} phased germline mutations'.format(len(gmuts.index)))
            count = list(gmuts.columns[gmuts.columns.str.endswith('_count')])
            gmuts[count] = gmuts[count].astype(int)
            logging.info('writing {}.phased.vcf'.format(args.normal_name))
            write_phase_vcf(args.normal_name+'.phased.vcf',
                            gmuts, n, ['GermlineID'])
            logging.info('writing {}.phased.maflite.tsv'.format(args.normal_name))
            write_phase_maflite(args.normal_name+'.phased.maflite.tsv',
                                gmuts,
                                ['genotype','GermlineID'])
    logging.info('all done!')
