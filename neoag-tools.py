#!/usr/bin/env/python3
# author: Junko Tsuji

import os
import sys
import logging
import argparse

sys.path.append(os.path.dirname(os.path.realpath(__file__)) + "/")

import phase.phase
import epitope.epitope


def build_parser():
    parser = argparse.ArgumentParser(description='Run neoantigen generation tools')
    
    # common parameters
    base_parser = argparse.ArgumentParser(add_help=False)
    base_parser.add_argument('-r', '--ref-fasta',
                             type=str,
                             action='store',
                             help='Reference fasta, required')
    base_parser.add_argument('-m', '--input-maf',
                             type=str,
                             action='store',
                             help='Input MAF, required')
    base_parser.add_argument('-t', '--tumor-name',
                             type=str,
                             action='store',
                             help='Tumor sample name (it has to match with "Tumor_Sample_Barcode" in MAF), required')
    base_parser.add_argument('-n', '--normal-name',
                             type=str,
                             action='store',
                             help='Normal sample name')

    # tool specific parser
    subparsers = parser.add_subparsers(title='neoag-tools',
                                       description='Choose a tool to run',
                                       help='')
    
    phasing = subparsers.add_parser('phase',
                                    help='Generate phased somatic (and germline) mutations',
                                    parents=[base_parser])
    # change 'optional arguments' to 'input arguments'
    phasing._action_groups[1].title = 'input arguments'

    phasing.add_argument('--tumor-bam',
                         type=str,
                         action='store',
                         help='Tumor BAM, required')
    phasing.add_argument('--normal-bam',
                         type=str,
                         action='store',
                         help='Normal BAM, optional')
    phasing.add_argument('--tumor-vcf',
                         type=str,
                         action='store',
                         help='Somatic variant VCF, optional')
    phasing.add_argument('--normal-vcf',
                         type=str,
                         action='store',
                         help='Germline variant VCF, optional')
    phasing.add_argument('--phylogic-ccfs',
                         type=str,
                         action='store',
                         help='Clustered somatic mutations \'mut_ccfs\' from PhylogicNDT, optional')
    phasing.add_argument('--phylogic-tree',
                         type=str,
                         action='store',
                         help='Clonal tree information \'build_tree_posteriors\' from PhylogicNDT, optional')
    phasing.add_argument('--phylogic-index',
                         type=int,
                         action='store',
                         default=1,
                         help='N-th tree in PhylogicNDT for clonal structure (default=%(default)s)')
    phasing.add_argument('--vaf-skew-pval',
                         type=float,
                         action='store',
                         default=0.05,
                         help='P-values to evaluate if a pair of mutation VAFs is skewed (default=%(default)s)')
    phasing.add_argument('--max-mnp-size',
                         type=int,
                         action='store',
                         default=10,
                         help='Maximum length of phased somatic mutations to chain as MNP/ONP (default=%(default)s)')
    phasing.add_argument('--max-phase-radius',
                         type=int,
                         action='store',
                         default=120,
                         help='Maximum radius from somatic mutations to incorporate phased germline mutations (default=%(default)s)')
    phasing.add_argument('--min-coverage',
                         type=int,
                         action='store',
                         default=10,
                         help='Minimum coverage to investigate phasing (default=%(default)s)')
    phasing.add_argument('--min-phased-altc',
                         type=int,
                         action='store',
                         default=3,
                         help='Minimum alt counts of phased mutation (default=%(default)s)')
    phasing.add_argument('--min-base-quality',
                         type=int,
                         action='store',
                         default=5,
                         help='Minimum base quality to investigate phasing (default=%(default)s)')
    phasing.add_argument('--min-map-quality',
                         type=int,
                         action='store',
                         default=5,
                         help='Minimum mapping quality to investigate phasing (default=%(default)s)')
    phasing.set_defaults(func=phase.phase.run)


    antigen = subparsers.add_parser('epitope',
                                    help='Construct neoantigens with supplied mutations',
                                    parents=[base_parser])
    # change 'optional arguments' to 'input arguments'
    antigen._action_groups[1].title = 'input arguments'
    antigen.add_argument('--gtf',
                         type=str,
                         action='store',
                         help='GENCODE gene annotation GTF, required')
    antigen.add_argument('--germline-maf',
                         type=str,
                         action='store',
                         help='Phased germline MAF, optional')
    antigen.add_argument('--gene-expr',
                         type=str,
                         action='store',
                         help='RSEM gene expression matrix, optional')
    antigen.add_argument('--codon-table',
                         type=str,
                         action='store',
                         help='Codon table to use, see "epitope/seq.py" for details (default="standard genetic code")')
    antigen.add_argument('--peptide-length',
                         type=str,
                         action='store',
                         default='8,9,10,11',
                         help='Comma-separated list of peptide lengths. Lengths 8, 9, 10, and 11 are supported (default=%(default)s)')
    antigen.add_argument('--flank-peplen',
                         type=int,
                         action='store',
                         default=30,
                         help='Upstream and downstream peptide lengths to report (default=%(default)s)')
    antigen.add_argument('--ignore-clones',
                         action='store_true',
                         help='Ignore mutation clonal structure during peptide translation (default=off)')
    antigen.add_argument('--ignore-phasing',
                         action='store_true',
                         help='Ignore phasing information during peptide translation (default=off)')
    antigen.set_defaults(func=epitope.epitope.run)
    
    if len(sys.argv) < 2:
        parser.print_help(sys.stderr)
        sys.exit()

    print_help = '-h' in sys.argv or '--help' in sys.argv
    if 'phase' in sys.argv:
        if len(sys.argv) < 3 or print_help:
            phasing.print_help(sys.stderr)
            sys.exit()

    if 'epitope' in sys.argv:
        if len(sys.argv) < 3 or print_help:
            antigen.print_help(sys.stderr)
            sys.exit()
        
    return parser.parse_args()

    
if __name__ == '__main__':
    args = build_parser()
    args.func(args)
