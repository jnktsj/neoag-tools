import os
import copy
import re
import logging
import itertools
import argparse
import sys

import numpy as np
import pandas as pd

from .seq import Seq
from .seq import standard_code
from .seq import create_codon_table
from .io import write_neoorf, write_peptide, write_fasta, read_rsem_gene, read_maf, neoorf_header, common_header, neoorf_mut_class, peptide_header, coding_mut_class
from .gtf import Annotation


def validate_path(name: str, path: str | None, is_required=True):
    """
    Checks that the path is present.

    Args:
        name (str): The name of the object pointed to by the path.
        path (str): A path to a necessary file.
    """
    if path is None:
        if is_required:
            raise Exception(f'Provide {name}')
    else:
        if not os.path.exists(path):
            raise FileNotFoundError(f'{name} not found: {path}')


def check_input_args(args: argparse.Namespace):
    """
    Check input parameters
    """
    # check required arguments
    validate_path('reference FASTA', args.ref_fasta)
    validate_path('input MAF', args.input_maf)
    validate_path('gene annotation GTF', args.gtf)
    validate_path('codon table', args.codon_table, is_required=False)
    validate_path('gene expression matrix', args.gene_expr, is_required=False)

    if args.tumor_name is None:
        raise Exception('provide tumor sample name in the input MAF')
    

    # check optional arguments
    if args.normal_name or args.germline_maf:
        if args.normal_name is None:
            raise Exception('provide normal sample name')
        if args.germline_maf is None:
            raise Exception('cannot find germline MAF')
        if not os.path.exists(args.germline_maf):
            raise FileNotFoundError('germline MAF not found: '+args.germline_maf)
        
    if args.flank_peplen < 0:
        raise ValueError('invalid peptide upstream/downstream length; provide a value larger than or equal to 0')
    
    if args.peptide_length == '':
        raise ValueError('provide peptide length(s)')
    for p in args.peptide_length.split(','):
        try:
            pep_len = int(p)
            if pep_len not in [8, 9, 10, 11]:
                raise ValueError('invalid peptide length: {} ({})'.format(p, args.peptide_length))
        except ValueError:
            raise ValueError('invalid peptide length: {} ({})'.format(p, args.peptide_length))


def mutate_sequence(txid, seq, pos, muts: pd.DataFrame):
    mut_seq = copy.deepcopy(seq)
    idx_seq = [[-1 for x in s] for s in mut_seq]
    for x, m in muts.iterrows():
        ref = m['Reference_Allele']
        alt = m['Tumor_Seq_Allele']
        for index in range(len(pos)):
            if pos[index].start <= m['Start_position'] and m['Start_position'] <= pos[index].end:
                start = m['Start_position'] - pos[index].start
                break
        else:
            if m['Variant_Classification'] in coding_mut_class:
                logging.warning('skipping {}:{}; position not in {}'.format(
                    m['Chromosome'], m['Start_position'], txid
                ))
            continue
        if m['Variant_Type'] == 'DEL':
            end = start + len(ref)
            mut_seq[index][start:end] = [''] * len(ref)
            idx_seq[index][start:end] = [x] * len(ref)
        elif m['Variant_Type'] == 'INS':
            mut_seq[index][start] += alt
            idx_seq[index][start] = x
        elif len(ref) == len(alt):
            end = start + len(alt)
            for i, mb in zip(range(start, end), list(alt)):
                if mut_seq[index][i] == '':
                    continue
                mut_seq[index][i] = mb + mut_seq[index][i][1:]
                idx_seq[index][i] = x
        else:
            logging.warning('unknown variant type: ' + m['Variant_Type'])
    return mut_seq, idx_seq


def get_mutation_notation(wt, mt, idx):
    mt_str = ''
    mt_pos = []
    mt_idx = dict()
    for i, (wt_base,mt_base) in enumerate(zip(wt, mt)):
        if mt_base == '':
            mt_str += 'd'*len(wt_base)
            if idx[i] not in mt_pos:
                mt_pos.append(idx[i])
        elif len(mt_base) > len(wt_base):
            mt_str += ' '*len(wt_base) + 'i'*(len(mt_base)-len(wt_base))
            mt_pos.append(idx[i])
        elif mt_base != wt_base:
            mt_str += 'm'*len(mt_base)
            mt_pos.append(idx[i])
        else:
            mt_str += ' '*len(wt_base)
            
    # move sequential deletion characters to deletion edges
    mt_str = re.sub(' d+ ', 'dd', mt_str)

    for i, m in enumerate(re.finditer('[mdi]+', mt_str)):
        mt_idx[mt_pos[i]] = m.span()
    return mt_str, mt_idx


def get_peptide_notation(wt: str, mt: str) -> str:
    """
    Generates a string graphically showing the differences between a wild-type sequence and a mutated sequence.

    Args:
        wt (str): A wild-type sequence.
        mt (str): A mutated sequence.

    Returns:
        str: A string with ' ' where the WT and mutated sequences match, and '*' where they don't.
    """
    return ''.join('*' if w == m else ' ' for w, m in zip(wt, mt))
    aa_str = ""
    for i in range(min(len(wt), len(mt))):
        if wt[i] != mt[i]:
            aa_str += "*"
        else:
            aa_str += " "
    return aa_str


def translate_mutation(fo_fasta, fo_peptide, fo_neoorf, tumor_name, smuts,
                       normal_name, gmuts, gtf, gene_tpm, flank_length,
                       peptide_lengths, genome, codon_table):
    for cl in smuts:
        for txid in smuts[cl]['tx']:
            tx = gtf.transcripts[txid]
            if tx.is_coding == False:
                logging.warning('skipping {}; CDS not defined'.format(txid))
                continue
                
            seq, pos, exon_str = tx.raw_cds_with_downstream(genome)

            # get coding mutation index
            vc_i = smuts[cl]['mut']['Variant_Classification'].isin(coding_mut_class)
            at_i = smuts[cl]['mut']['Annotation_Transcript'] == txid
            mx = list(smuts[cl]['mut'][vc_i & at_i].index)
                
            # append germline mutations
            gx = []
            if not gmuts.empty:
                gv_i = gmuts['Variant_Classification'].isin(coding_mut_class)
                gt_i = gmuts['Annotation_Transcript'] == txid
                gx = list(gmuts[gv_i & gt_i].index)
                seq, _ = mutate_sequence(txid, seq, pos, gmuts)

            # append somatic mutations
            mseq, midx = mutate_sequence(txid, seq, pos, smuts[cl]['mut'])
                
            # remove exon structure
            seq = list(itertools.chain(*seq))
            mseq = list(itertools.chain(*mseq))
            midx = list(itertools.chain(*midx))
                
            wt_seq = Seq(''.join(seq))
            mt_seq = Seq(''.join(mseq))
                
            if tx.pos.strand == "-":
                wt_seq.reverse_complement()
                mt_seq.reverse_complement()
                exon_str = exon_str[::-1]
                midx = midx[::-1]
                mx = mx[::-1]
                gx = gx[::-1]
                # note: these sequences are not reverse complement, just reverse.
                #       but it's sufficient for mutation notation.
                seq = seq[::-1]
                mseq = mseq[::-1]

            wt_aa, _ = wt_seq.translate(codon_table=codon_table, padchar=' ')
            mt_aa, _ = mt_seq.translate(codon_table=codon_table, padchar=' ')

            mt_str, mt_idx = get_mutation_notation(seq, mseq, midx)
            aa_str = get_peptide_notation(wt_aa, mt_aa)
                
            write_fasta(fo_fasta, tumor_name, normal_name,
                        smuts[cl]['mut'].loc[mx], gmuts.loc[gx], txid,
                        wt_seq, mt_seq, aa_str, mt_str, exon_str, cl)
            
            tpm = gene_tpm.get(txid.split('.')[0], 'nan')
            write_peptide(fo_peptide, smuts[cl]['mut'].loc[mx],
                          cl, wt_seq, mt_seq, aa_str, mt_idx,
                          tumor_name, txid, tx.gene.id, tpm,
                          flank_length, peptide_lengths)

            nf_i = smuts[cl]['mut'].loc[mx, 'Variant_Classification'].isin(neoorf_mut_class)
            mx = list(smuts[cl]['mut'][vc_i & at_i & nf_i].index)
            if len(mx) > 0:
                write_neoorf(fo_neoorf, smuts[cl]['mut'].loc[mx],
                             cl, wt_seq, mt_seq, aa_str, mt_idx,
                             tumor_name, txid, tx.gene.id, tpm)


def run(args):
    check_input_args(args)

    logging.basicConfig(level=logging.INFO,
                        format='epitope | %(asctime)s [%(levelname)s] %(message)s',
                        datefmt='%Y-%m-%d %H:%M:%S')
    
    # load mutations
    logging.info('loading mutations from input MAF')
    dfs = read_maf(args.input_maf, args.tumor_name, 'Tumor_Sample_Barcode')
    dfg = read_maf(args.germline_maf, args.normal_name, 'Matched_Norm_Sample_Barcode')
    
    # check phasing and clonal structure columns
    if args.ignore_phasing:
        dfs['PhaseID'] = np.nan
    else:
        if 'PhaseID' not in dfs:
            raise Exception('cannot find "PhaseID" in: '+args.input_maf)
    if args.ignore_clones:
        dfs['ClonalStructure'] = np.nan
    else:
        if 'ClonalStructure' not in dfs:
            raise Exception('cannot find "ClonalStructure" in: '+args.input_maf)
    if not dfg.empty:
        if 'GermlineID' not in dfs:
            raise Exception('cannot find "GermlineID" in: '+args.input_maf)
        if 'GermlineID' not in dfg:
            raise Exception('cannot find "GermlineID" in: '+args.germline_maf)

    # exit if there are no coding transcripts with return code 0
    row = dfs['Variant_Classification'].isin(coding_mut_class)
    if sum(row) == 0:
        logging.warning('no coding mutations found! exiting...')
        sys.exit(0)
    transcripts = list(dfs.loc[row, 'Annotation_Transcript'].unique())

    # load transcript structures from GTF
    logging.info('loading {} transcripts with mutations from GTF'.format(len(transcripts)))
    gtf = Annotation(args.gtf, transcript_list=transcripts)
    logging.info('finished loading {} transcripts from GTF'.format(len(gtf.transcripts)))

    # load codon table if supplied
    codon_table = standard_code
    if args.codon_table:
        logging.info('custom codon table supplied: {}'.format(args.codon_table))
        codon_table = create_codon_table(args.codon_table)

    gene_tpm = dict()
    # load RSEM gene matrix if supplied
    if args.gene_expr:
        logging.info('loading gene expression data')
        gene_tpm = read_rsem_gene(args.gene_expr, transcript_list=transcripts)

    peptide_lengths = list(map(int,args.peptide_length.split(',')))
        
    # output FASTA
    fo_fasta = open(args.tumor_name + '.muts.fa', 'w')

    # output peptide file
    fo_peptide = open(args.tumor_name + '.peptide.tsv', 'w')
    fo_peptide.write('\t'.join(peptide_header + common_header)+'\n')

    # output neoorf file
    fo_neoorf = open(args.tumor_name + '.neoorf.tsv', 'w')
    fo_neoorf.write('\t'.join(neoorf_header + common_header)+'\n')

    logging.info('start translating mutations')
    analyzed = []
    for i, mut in dfs.loc[row].iterrows():
        if i in analyzed:
            continue
        
        index = [i]
        if pd.isna(mut['PhaseID']) == False:
            index.extend(list(dfs[dfs['PhaseID']==mut['PhaseID']].index))
            index = sorted(set(index))
        x = dfs.loc[index].reset_index(drop=True).copy()
        analyzed.extend(index)

        # subset somatic mutations
        smuts = dict()
        clones = list(x['ClonalStructure'].dropna())
        if len(clones) > 0:
            for cl in sorted(set.union(*clones)):
                cs = x['ClonalStructure'].apply(lambda clstr: cl in clstr)
                xi = x.loc[cs].copy()
                # avoid a clone composed of all non-coding mutations
                is_coding_muts = xi['Variant_Classification'].isin(coding_mut_class)
                if sum(is_coding_muts) == 0:
                    continue
                xi = xi.sort_values('Start_position', ignore_index=True)
                tx = set(list(xi.loc[is_coding_muts, 'Annotation_Transcript']))
                smuts.setdefault(cl, {'mut':xi, 'tx':tx})
        else:
            is_coding_muts = x['Variant_Classification'].isin(coding_mut_class)
            tx = set(list(x.loc[is_coding_muts, 'Annotation_Transcript']))
            smuts.setdefault('0', {'mut':x, 'tx':tx})

        # subset germline mutations if supplied
        gmuts = pd.DataFrame()
        if not dfg.empty:
            gid = list(x['GermlineID'].dropna())
            if len(gid) > 0:
                gmuts = dfg[dfg['GermlineID'].isin(set(gid))].copy()

        translate_mutation(fo_fasta, fo_peptide, fo_neoorf, args.tumor_name, smuts,
                           args.normal_name, gmuts, gtf, gene_tpm, args.flank_peplen,
                           peptide_lengths, args.ref_fasta, codon_table)
    # close file handles and end
    fo_fasta.close()
    fo_peptide.close()
    fo_neoorf.close()
    logging.info('all done!')
