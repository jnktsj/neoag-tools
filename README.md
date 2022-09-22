# `neoag-tools`: Modules for neoantigen analysis

This repo provides two main modules for neoantigen analysis.

## Requirements (version tested)
 - Python3 (3.9.7)
 - NumPy (1.21.4)
 - SciPy (1.8.0)
 - pandas (1.3.4)
 - Pysam (0.18.0)

`requiements.txt` includes the exact versions tested, but
the modules should be most likely functional with the latest
versions of the above packages.

## Phase mode

This tool connects separately reported somatic and germline
mutations if these exist in cis (i.e. phased mutations).

### Input

#### Required input

The tool can run without a tumor VCF. If a tumor VCF is supplied,
the output VCFs will inherit the headers in the supplied VCF.
If the VCF is not available, the tool will generate custom VCF headers.
See `[Example 1]`.

 - Reference genome (required)
 - Tumor MAF (required)
 - Tumor sample name (has to match with 'Tumor_Sample_Barcode' in MAF and 'Sample_ID' in PhylogicNDT; required)
 - Tumor BAM (required)
 - Tumor VCF (optional)

#### Germline input

**Highly recommended to use germline information for more accurate phasing analysis.**
Some germline mutations are adjacent to somatic mutations in cis, and
these impact amino acid changes. For obtaining accurate patient/sample-specific
neopeptides and neo-ORFs, it is important to take into account such
cis-germline mutations. See `[Example 2]`.

 - Normal sample name
 - Normal BAM
 - Normal VCF

#### Clonal information input

Due to tumor heterogeneity, there is a case that subclonal mutations
exist in cis of clonal mutations. Although this observation is rare,
the tool can add tumor clonal structure information to phased somatic
mutations by using PhylogicNDT results. See `[Example 3]`.

 - Phylogic CCFs (`*.mut_ccfs.txt`)
 - Phylogic tree (`*_build_tree_posteriors.tsv`)
 - Phylogic index (integer value of a specific tree, e.g. 1)


### Output

#### <tumor_name>.unphased.maf

MAF file with unphased somatic mutations. If no phased mutations are
found, this MAF contents will be identical to the input MAF.

#### <tumor_name>.phased.vcf

VCF file with phased somatic mutations. This file will be generated
only if phased somatic mutations exist. If a tumor VCF is supplied,
the tool will use the header information of the VCF. As the downstream
analysis, this file can be used as an input for an annotation tool, e.g.
Funcotator.

#### <tumor_name>.phased.maflite.tsv

MAFLITE file with phased somatic mutations. The file will be generated
only if phased somatic mutations exist. As the downstream analysis,
this file can be used as an input for an annotation tool, e.g.
Oncotator.


#### <normal_name>.phased.vcf

VCF file with germline mutations phased with somatic mutations.
`GermlineID` field in both somatic and germline VCFs links
between phased somatic and germline mutations. This file
will be generated only if normal BAM/VCF is supplied, and
such phased germline mutations exist.

#### <normal_name>.phased.maflite.tsv
 
MAFLITE file with germline mutations phased with somatic mutations.
`GermlineID` field in both somatic and germline MAFLITEs links
between phased somatic and germline mutations. This file
will be generated only if normal BAM/VCF is supplied, and
such phased germline mutations exist.


### Examples


```shell
# [Example 1] Basic run with tumor only
python3 neoag-tools.py phase \
  -r genome.fa \
  -m somatic_muts.maf \
  -t tumor_name \
  --tumor-bam tumor.bam \
  --tumor-vcf tumor.vcf

# [Example 2] Add phased germline mutations (recommended)
python3 neoag-tools.py phase \
  -r genome.fa \
  -m somatic_muts.maf \
  -t tumor_name \
  --tumor-bam tumor.bam \
  --tumor-vcf tumor.vcf \
  -n normal_name \
  --normal-bam normal.bam \
  --normal-vcf normal.vcf

# [Example 3] Add clonal structure to phased mutations
python3 neoag-tools.py phase \
  -r genome.fa \
  -m somatic_muts.maf \
  -t tumor_name \
  --tumor-bam tumor.bam \
  --tumor-vcf tumor.vcf \
  -n normal_name \
  --normal-bam normal.bam \
  --normal-vcf normal.vcf \
  --phylogic-ccfs phylogic.mut_ccfs.txt \
  --phylogic-tree phylogic_build_tree_posteriors.tsv \
  --phylogic-index 1
```


## Epitope mode

This tool translates phased germline and somatic mutations
processed with the phase mode into peptide sequences.
The tool will read clonal structure information if available.

### Input

#### Required input

  - Reference genome
  - Tumor MAF
  - Tumor sample name (has to match with 'Tumor_Sample_Barcode' in MAF)
  - GENCODE GTF gene annotation

See `[Example 1]`.

#### Germline input

**Highly recommended to use germline information for more accurate translation.**

 - Normal sample name
 - Normal MAF

See `[Example 2]`.

#### Gene expression

RSEM gene expression matrix can be supplied to attach TPM information.

#### Peptide lengths

Only lengths 8, 9, 10, and 11 are supported.

#### Codon table

If an input sample needs to use non-standard genetic code, a custom
codon table can be supplied. Download the corresponding code from
NCBI's genetic code website
([link](https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi)).
To find which genetic code you need to use, search your species
through NCBI's taxonomy database
([link](https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi)).

Once you find the genetic code of interest, copy and paste the
following format from the database and create a text file.

```
  AAs  = FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG
Starts = ---M------**--*----M---------------M----------------------------
Base1  = TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG
Base2  = TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG
Base3  = TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG
```

### Examples

```shell
# [Example 1] Basic run with tumor only
python3 neoag-tools.py epitope \
  -r genome.fa \
  -m somatic_muts.maf \
  -t tumor_name \
  --gtf gencode.gtf


# [Example 2] Add phased germline mutations (recommended)
python3 neoag-tools.py epitope \
  -r genome.fa \
  -m somatic_muts.maf \
  -t tumor_name \
  --gtf gencode.gtf
  -n normal_name \
  --germline-maf germline_muts.maf \
  --gene-expr rsem.gene.results.gz
```
