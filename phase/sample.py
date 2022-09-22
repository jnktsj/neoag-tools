import os

import numpy as np
import pysam

class Sample:
    """
    bam   : bam file name
    bam_sm: sample name in BAM
    vcf   : vcf file name
    vcf_sm: sample name in VCF
    """
    def __init__(self, bam, vcf):
        self.bam = None
        self.bam_sm = None
        self.vcf = None
        self.vcf_sm = None
        
        if bam:
            self.bam = pysam.AlignmentFile(bam, 'rb')
            bam_sm = set([rg['SM'] for rg in self.bam.header['RG']])
            if len(bam_sm) > 1:
                raise ValueError('Multiple sample names in {}: {}'.format(
                    os.path.basename(bam), ' '.join(bam_sm)
                ))
            self.bam_sm = bam_sm.pop()
            self.bam = bam

        if vcf:
            self.vcf = pysam.VariantFile(vcf, 'r')
            vcf_sm = set(self.vcf.header.samples)
            if len(vcf_sm) > 2:
                raise ValueError('Multiple sample names in {}: {}'.format(
                    os.path.basename(vcf), ' '.join(vcf_sm)
                ))
            if self.bam_sm in vcf_sm:
                self.vcf_sm = self.bam_sm
            elif 'TUMOR' in vcf_sm:
                self.vcf_sm = 'TUMOR'
            elif len(vcf_sm) == 1:
                self.vcf_sm = vcf_sm.pop()
            else:
                raise ValueError('Cannot read sample name in {}: {}'.format(
                    os.path.basename(vcf), ' '.join(bam_sm)))
            if not self.bam_sm and self.vcf_sm != 'TUMOR':
                self.bam_sm = self.vcf_sm
            self.vcf = vcf
