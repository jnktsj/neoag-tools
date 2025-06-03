import os
from typing import Optional

import numpy as np
import pysam

class Sample:
    """
    bam   : bam file name
    bam_sm: sample name in BAM
    vcf   : vcf file name
    vcf_sm: sample name in VCF
    """
    bam: str
    """The path to the BAM file."""
    bam_sm: str
    """The BAM sample name."""
    vcf: Optional[str]
    """The path to the sample VCF file."""
    vcf_sm: str
    """The sample name in the VCF file."""

    def __init__(self, bam_path: str, vcf_path: Optional[str]):
        self.bam = bam_path
        self.vcf = vcf_path

        # Getting the sample name from the BAM file.
        bam = pysam.AlignmentFile(self.bam, 'rb')
        bam_sample_names = set([rg['SM'] for rg in bam.header['RG']])  # type: ignore
        if len(bam_sample_names) > 1:
            raise ValueError('Multiple sample names in {}: {}'.format(
                os.path.basename(self.bam), ' '.join(bam_sample_names)
            ))
        self.bam_sm = bam_sample_names.pop()

        if vcf_path is not None:
            # Getting the sample name from the VCF.
            vcf = pysam.VariantFile(vcf_path, 'r')
            vcf_sample_names = set(vcf.header.samples)

            if self.bam_sm in vcf_sample_names:
                self.vcf_sm = self.bam_sm
            if len(vcf_sample_names) > 2:
                raise ValueError('Multiple sample names in {}: {}'.format(
                    os.path.basename(vcf_path), ' '.join(vcf_sample_names)
                ))
            elif 'TUMOR' in vcf_sample_names:
                self.vcf_sm = 'TUMOR'
            elif len(vcf_sample_names) == 1:
                self.vcf_sm = vcf_sample_names.pop()
            else:
                raise ValueError('Cannot read sample name in {}: {}'.format(
                    os.path.basename(vcf_path), ' '.join(vcf_sample_names)))
