import pysam
import numpy as np

def write_subset_maf(outname, index, maf):
    """
    write out MAF file based on supplied index values
    """
    fo = open(outname, 'w')
    i = -1
    for line in open(maf, 'r'):
        if line.startswith('#'):
            continue
        if line.startswith('Hugo'):
            fo.write(line)
            continue
        i += 1
        if i in index:
            fo.write(line)
    fo.close()

    
def write_phase_vcf(outname, df, sample, extra_cols=[]):
    """
    write out VCF file based on supplied index values
    the function reuses the original VCF information
    """
    rvcf = pysam.VariantFile(sample.vcf, 'r')
    t_name = sample.vcf_sm
    n_name = set(rvcf.header.samples).difference(set([t_name]))
    n_name = n_name.pop() if n_name else None
    is_somatic = False
    if 'ClonalStructure' in df:
        is_somatic = True
        rvcf.header.info.add('PhaseID',
                             number=1,
                             type='String',
                             description='Block ID for phased somatic variants')
        rvcf.header.info.add('MNPVariantType',
                             number=1,
                             type='String',
                             description='Oligo-nucleotide polymorphism (DNP,TNP,MNP)')
        rvcf.header.info.add('ClonalStructure',
                             number=1,
                             type='String',
                             description='Possible clones having the mutation')
        rvcf.header.add_meta('tumor_sample', value=t_name)
        if n_name:
            rvcf.header.add_meta('normal_sample', value=n_name)
        if 'ccf_hat' in df:
            rvcf.header.info.add('ccf_hat',
                                 number=1,
                                 type='Float',
                                 description='Cancer cell fraction of the mutation')
    if 'GermlineID' in df:
        rvcf.header.info.add('GermlineID',
                             number=1,
                             type='String',
                             description='Block ID for phased germline and somatic variants')
    ovcf = pysam.VariantFile(outname, 'w', header=rvcf.header)                             
    for _, row in df.iterrows():
        alt = row['Tumor_Seq_Allele']
        mut_type = row['Variant_Type']
        chrom = row['Chromosome']
        pos = row['Start_position']-(0 if mut_type == 'INS' else 1)
        end = row['Start_position']+(1 if mut_type == 'INS' else len(alt)-1)
        for v in rvcf.fetch(chrom, pos, end):
            if row['Start_position'] != v.pos:
                continue
            for c in extra_cols:
                v.info[c] = str(row[c])
            if is_somatic:
                v.ref = row['Reference_Allele']
                v.alts = (alt,)
                refc = row['t_ref_count']
                altc = row['t_alt_count']
                v.samples[t_name]['AD'] = (refc, altc)
                v.samples[t_name]['DP'] = refc + altc
                v.samples[t_name]['AF'] = altc/float(refc+altc)
                v.samples[t_name]['GT'] = (0, 1)
                if n_name and 'n_ref_count' in df:
                    refc = row['n_ref_count']
                    altc = row['n_alt_count']
                    v.samples[n_name]['AD'] = (refc, altc)
                    v.samples[n_name]['DP'] = refc + altc
                    v.samples[n_name]['AF'] = altc/max(float(refc+altc), 1)
                    v.samples[n_name]['GT'] = (0, 0)
                if len(row['maf_idx']) > 1:
                    v.info['MNPVariantType'] = row['Variant_Type']
            ovcf.write(v)
    ovcf.close()
    rvcf.close()


def write_phase_vcf_from_scratch(outname, t_name, n_name, df):
    """
    write out VCF file from scratch
    if the original VCF is not supplied
    """
    info_cols = ['PhaseID', 'ClonalStructure']
    
    fo = open(outname, 'w')

    # write header
    fo.write('##fileformat=VCFv4.2\n')
    fo.write('##source=neoag-tools\n')
    fo.write('##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths of ref and alt alleles">\n')
    fo.write('##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Total read depth">\n')
    fo.write('##FORMAT=<ID=AF,Number=A,Type=Float,Description="Allelic fractions of alt alleles">\n')
    fo.write('##INFO=<ID=PhaseID,Number=1,Type=String,Description="Block ID for phased somatic variants">\n')
    fo.write('##INFO=<ID=MNPVariantType,Number=1,Type=String,Description="Oligo-nucleotide polymorphism (DNP,TNP,MNP)">\n')
    fo.write('##INFO=<ID=ClonalStructure,Number=1,Type=String,Description="Possible clones having the mutation">\n')
    if 'GermlineID' in df:
        fo.write('##INFO=<ID=GermlineID,Number=1,Type=String,Description="Block ID for phased germline and somatic variants">\n')
        info_cols.append('GermlineID')
    if 'ccf_hat' in df:
        fo.write('##INFO=<ID=ccf_hat,Number=1,Type=Float,Description="Cancer cell fraction of the mutation">\n')
        info_cols.append('ccf_hat')
        
    # write sample names
    fo.write('##tumor_sample={}\n'.format(t_name))
    if n_name and 'n_alt_count' in df:
        fo.write('##normal_sample={}\n'.format(n_name))
    fo.write('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t'+'\t'.join([n_name, t_name]).lstrip('\t')+'\n')
    
    # write mutation records
    for _, row in df.iterrows():
        info_str = ';'.join(['{}={}'.format(c, str(row[c])) for c in info_cols])
        if len(row['maf_idx']) > 1:
            info_str += ';MNPVariantType={}'.format(row['Variant_Type'])
        refc = row['t_ref_count']
        altc = row['t_alt_count']
        tumor_format = [ str(refc)+','+str(altc),
                         str(refc+altc),
                         str(altc/float(refc+altc)) ]
        v = [ row['Chromosome'],
              row['Start_position'],
              '.',
              row['Reference_Allele'],
              row['Tumor_Seq_Allele'],
              '.',
              'PASS',
              info_str,
              'AD:DP:AF',
              ';'.join(tumor_format) ]
        if n_name and 'n_ref_count' in df:
            refc = row['n_ref_count']
            altc = row['n_alt_count']
            normal_format = [ str(refc)+','+str(altc),
                              str(refc+altc),
                              str(altc/float(refc+altc))]
            v.insert(len(v)-1, ';'.join(normal_format))
        fo.write('\t'.join(list(map(str,v)))+'\n')
    fo.close()


def write_phase_maflite(outname, df, extra_cols=[]):
    """
    write out MAFLITE file for Oncotator input
    """
    h = ['chr', 'start', 'end', 'ref_allele', 'alt_allele',
         't_ref_count', 't_alt_count']
    if 'n_ref_count' in df:
        h.extend(['n_ref_count', 'n_alt_count'])
    h.extend(extra_cols)
    
    # insert end cooridnates
    offset = df['Tumor_Seq_Allele'].str.len()-1
    offset[df['Variant_Type'] == 'INS'] = 1
    df['end'] = df['Start_position'] + offset

    x = df.rename(columns={'Chromosome':'chr',
                           'Start_position': 'start',
                           'Reference_Allele': 'ref_allele',
                           'Tumor_Seq_Allele': 'alt_allele'})
    x[h].to_csv(outname, index=False, sep='\t', na_rep='nan')
