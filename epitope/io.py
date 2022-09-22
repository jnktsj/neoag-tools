import os
import re
import gzip
from collections import OrderedDict

import numpy as np
import pandas as pd

from .seq import unpad_peptide


# coding mutations
# excluding 'Splice_Site' from amino acid translation
coding_mut_class = [
    'Missense_Mutation',
    'Frame_Shift_Ins',
    'Frame_Shift_Del',
    'Nonstop_Mutation',
    'In_Frame_Ins',
    'In_Frame_Del'
]
neoorf_mut_class = [
    'Frame_Shift_Ins',
    'Frame_Shift_Del',
    'Nonstop_Mutation'
]

# maf columns
required_cols = OrderedDict({
    'Hugo_Symbol'           : str,
    'Chromosome'            : str,
    'Start_position'        : int,
    'End_position'          : int,
    'Reference_Allele'      : str,
    'Tumor_Seq_Allele'      : str,
    'Variant_Classification': str,
    'Variant_Type'          : str,
    'Genome_Change'         : str,
    'Annotation_Transcript' : str,
    'cDNA_Change'           : str,
    'Codon_Change'          : str,
    'Protein_Change'        : str,
})
optional_cols = OrderedDict({
    'ClonalStructure': np.dtype('object'),
    'PhaseID'        : float,
    'GermlineID'     : float
})

# peptide and neoorf file columns
peptide_header = ['Hugo_Symbol',
                  'pep_wt',
                  'pep',
                  'ctex_up',
                  'ctex_dn',
                  'pep_start',
                  'pep_end']

neoorf_header = ['Hugo_Symbol',
                 'upstream',
                 'neoORF',
                 'neoORF_start']

common_header = ['transcript_id',
                 'gene_id',
                 'TPM',
                 'Tumor_Sample_Barcode',
                 'clone',
                 'Chromosome',
                 'Start_position',
                 'End_position',
                 'Variant_Classification',
                 'Variant_Type',
                 'Genome_Change',
                 'cDNA_Change',
                 'Codon_Change',
                 'Protein_Change']


def read_maf(maf, name, name_col):
    if maf == None:
        return pd.DataFrame()

    # reading MAF in a primivie way to handle commend characters '#'
    # in the middle of the line which pd.read_csv couldn't handle
    header = []
    record = []
    for line in open(maf):
        if line.startswith('#'):
            continue
        line = line.rstrip('\n').split('\t')
        if line[0] == 'Hugo_Symbol':
            header = line
        else:
            record.append({h:x for h,x in zip(header,line)})
    df = pd.DataFrame(data=record, dtype=str)

    tid = set(df[name_col])
    if name not in tid:
        raise Exception(name + " not in " + os.path.basename(maf))
    if len(tid) > 1:
        row = df[name_col] == name
        df[row].reset_index(drop=True, inplace=True)
    tid = name

    # rename MAF header
    df.rename(columns={'Tumor_Seq_Allele2':'Tumor_Seq_Allele'}, inplace=True)
    if 'Start_Position' in df:
        df.rename(columns={'Start_Position':'Start_position'}, inplace=True)
    if 'End_Position' in df:
        df.rename(columns={'End_Position':'End_position'}, inplace=True)

    # format MAF columns
    if 'PhaseID' in df:
        df.loc[df['PhaseID'].isin(['nan','']), 'PhaseID'] = np.nan
    if 'GermlineID' in df:
        df.loc[df['GermlineID'].isin(['nan','']), 'GermlineID'] = np.nan
        
    # subset MAF columns
    maf_cols = required_cols.copy()
    for oc in optional_cols:
        if oc in df:
            maf_cols[oc] = optional_cols[oc]
    try:
        df = df[maf_cols.keys()].astype(maf_cols)
    except Exception as e:
            raise KeyError(str(e).replace('index', 'MAF'))

    # process clone information if supplied
    if 'ClonalStructure' in df:
        df.loc[df['ClonalStructure'].isin(['nan','']), 'ClonalStructure'] = np.nan
        i = ~df['ClonalStructure'].isna()
        make_set = lambda x: set([xi.strip() for xi in re.sub('\[|\]','',x).split(',') if xi.strip().isnumeric()])
        df.loc[i, 'ClonalStructure'] = df.loc[i, 'ClonalStructure'].apply(make_set)
        
    return df


def read_rsem_gene(rsem_gene, transcript_list=[]):
    """
    Returns Transcript ID (without versioning) and TPM
    from RSEM gene expression matrix
    """
    # strip off transcript ID versions
    if transcript_list != []:
        transcript_list = set([txid.split('.')[0] for txid in transcript_list])

    gene_tpm = dict()
    try:
        if rsem_gene.endswith('.gz'):
            rsem_gene_fi = gzip.open(rsem_gene, 'rt')
        else:
            rsem_gene_fi = open(rsem_gene, 'r')
    except IOError:
        raise Exception('cannot read: {}'.format(rsem_gene))

    header = []
    for row in rsem_gene_fi:
        # skip comment lines if any
        if row.startswith('#'):
            continue
        row = row.rstrip('\n').split('\t')
        if row[0] == 'gene_id':
            continue
        
        transcript_ids = set([txid.split('.')[0] for txid in row[1].split(',')])
        tpm = row[5]
        
        if transcript_list != []:
            overlap = transcript_list.intersection(transcript_ids)
            if len(overlap) == 0:
                continue
            transcript_ids = overlap

        for txid in transcript_ids:
            gene_tpm.setdefault(txid, tpm)

    return gene_tpm


def get_fasta_header(m, name, txid, clone=None):
    gene = m['Hugo_Symbol'].unique()[0]
    nt = ';'.join(m['cDNA_Change'].unique())
    aa = ';'.join(m['Protein_Change'].unique())
    if clone != None:
        return '|'.join([name, txid, gene, nt, aa, 'clone='+clone])
    else:
        return '|'.join([name, txid, gene, nt, aa])


def write_fasta(fo, tumor_name, normal_name, muts, gmuts,
                txid, wt, mt, aa_str, mt_str, exon_str, clone):
    """
    Write peptide and cDNA fasta.

    Each line per record represents:
      1. FASTA entry header
      2. wild type protein sequence
      3. amino acid concordance string
      4. mutated protein sequence
      5. mutated cDNA sequence
      6. cDNA concordance string
      7. wild type cDNA sequence
      8. exon boundary string

    FASTA entry header is formatted:
      > {somatic_mutation_field} ; {germline_mutation_field}
 
    Germline mutation filed will be only written if there are
    phased coding germline mutations.
    """
    tumor_field = get_fasta_header(muts, tumor_name, txid, clone)
    normal_field = ''
    if not gmuts.empty:
        normal_field = get_fasta_header(gmuts, normal_name, txid)
    header = '>' + tumor_field + ';' + normal_field
    fo.write(header+'\n')
    fo.write(wt.aa+'\n')
    fo.write(aa_str+'\n')
    fo.write(mt.aa+'\n')
    fo.write(mt.seq+'\n')
    fo.write(mt_str+'\n')
    fo.write(wt.seq+'\n')
    fo.write(exon_str+'\n')


def write_peptide(fo, smuts, clone, wt, mt, mstr, idx,
                  name, txid, gnid, tpm, flen, pep_lens):
    flank_length = flen+max(pep_lens)-1
    for index, m in smuts.iterrows():
        start, end = idx[index]
        start = start - start%3 + 1
        end = (end-1) - (end-1)%3 + 1
        
        # adjust reading frame for shifted amino acid mismatches
        for i in range(start,len(mstr),3):
            if mstr[i] == '*':
                start = i
                break
        for i in range(end,start-3,-3):
            if mstr[i] == '*':
                end = i
        
        if m['Variant_Classification'].startswith('Frame_Shift'):
            wt_aa = unpad_peptide(wt.aa[start:])
            mt_aa = unpad_peptide(mt.aa[start:])
            # make wild type sequence same length as mutated sequence
            wt_aa += '-'*(len(mt_aa)-len(wt_aa))
            wt_aa = wt_aa[:len(mt_aa)]
            mt_aa_dn = ''
            wt_aa_dn = ''
        elif m['Variant_Classification'].startswith('Nonstop_Mutation'):
            mt_aa = unpad_peptide(mt.aa[start:])
            wt_aa = '-'*len(mt_aa)
            mt_aa_dn = ''
            wt_aa_dn = ''            
        else:
            wt_aa = unpad_peptide(wt.aa[start:end+1])
            mt_aa = unpad_peptide(mt.aa[start:end+1])
            mt_aa_dn = unpad_peptide(mt.aa[(end+3):])
            wt_aa_dn = unpad_peptide(wt.aa[(end+3):])
            
            # trim downstream sequence with specified flanking peptide length
            mt_aa_dn = mt_aa_dn[:flank_length] + '-'*(flank_length-len(mt_aa_dn))
            wt_aa_dn = wt_aa_dn[:flank_length] + '-'*(flank_length-len(wt_aa_dn))

        mt_aa_up = unpad_peptide(mt.aa[:start])
        wt_aa_up = unpad_peptide(wt.aa[:start])

        # mutated protein sequence coordinates
        aa_start = len(mt_aa_up) + 1
        aa_stop = aa_start + len(mt_aa) - 1

        # trim upstream sequence with specified flanking peptide length
        mt_aa_up = '-'*(flank_length-len(mt_aa_up)) + mt_aa_up[-flank_length:]
        wt_aa_up = '-'*(flank_length-len(wt_aa_up)) + wt_aa_up[-flank_length:]
        
        start = len(mt_aa_up)
        stop = start + len(mt_aa)
        
        wt_aa = wt_aa_up + wt_aa + wt_aa_dn
        mt_aa = mt_aa_up + mt_aa + mt_aa_dn
        for plen in pep_lens:
            for i in range(start-plen+1, stop):
                pep_wt = wt_aa[i:i+plen]
                pep = mt_aa[i:i+plen]
                if len(pep) < plen:
                    continue
                elif '-' in pep:
                    continue
                elif pep == pep_wt:
                    continue
                ctex_up = mt_aa[i-flen:i]
                ctex_dn = mt_aa[i+plen:i+plen+flen]
                ctex_dn += '-'*(flen-len(ctex_dn))
                fo.write('\t'.join([
                    m['Hugo_Symbol'],
                    pep_wt, pep, ctex_up, ctex_dn,
                    str(aa_start),
                    str(aa_stop),
                    txid, gnid, tpm, name, clone,
                    m['Chromosome'],
                    str(m['Start_position']),
                    str(m['End_position']),
                    m['Variant_Classification'],
                    m['Variant_Type'],
                    m['Genome_Change'],
                    m['cDNA_Change'],
                    m['Codon_Change'],
                    m['Protein_Change'],
                ])+'\n')


def write_neoorf(fo, smuts, clone, wt, mt, mstr, idx, name, txid, gnid, tpm):
    for index, m in smuts.iterrows():
        start, end = idx[index]
        start = start - start%3 + 1
        end = (end-1) - (end-1)%3 + 1

        if m['Variant_Classification'].startswith('Frame_Shift'):
            # adjust reading frame for shifted amino acid mismatches
            for i in range(start,len(mstr),3):
                if mstr[i] == '*':
                    start = i
                    break
        neoorf = unpad_peptide(mt.aa[start:])
        upstream = unpad_peptide(mt.aa[:start])
        neoorf_start = len(upstream) + 1

        fo.write('\t'.join([
            m['Hugo_Symbol'],
            upstream, neoorf, str(neoorf_start),
            txid, gnid, tpm, name, clone,
            m['Chromosome'],
            str(m['Start_position']),
            str(m['End_position']),
            m['Variant_Classification'],
            m['Variant_Type'],
            m['Genome_Change'],
            m['cDNA_Change'],
            m['Codon_Change'],
            m['Protein_Change']
        ])+'\n')
