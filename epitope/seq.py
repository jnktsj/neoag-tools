from typing import Any, Optional, Tuple
import pandas as pd

# standard code (transl_table=1)
standard_code = {
    'AAA': 'K', 'AAC': 'N', 'AAG': 'K', 'AAT': 'N', 
    'ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACT': 'T', 
    'AGA': 'R', 'AGC': 'S', 'AGG': 'R', 'AGT': 'S', 
    'ATA': 'I', 'ATC': 'I', 'ATG': 'M', 'ATT': 'I', 
    'CAA': 'Q', 'CAC': 'H', 'CAG': 'Q', 'CAT': 'H', 
    'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCT': 'P', 
    'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGT': 'R', 
    'CTA': 'L', 'CTC': 'L', 'CTG': 'L', 'CTT': 'L', 
    'GAA': 'E', 'GAC': 'D', 'GAG': 'E', 'GAT': 'D', 
    'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCT': 'A', 
    'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGT': 'G', 
    'GTA': 'V', 'GTC': 'V', 'GTG': 'V', 'GTT': 'V', 
    'TAA': '*', 'TAC': 'Y', 'TAG': '*', 'TAT': 'Y', 
    'TCA': 'S', 'TCC': 'S', 'TCG': 'S', 'TCT': 'S', 
    'TGA': '*', 'TGC': 'C', 'TGG': 'W', 'TGT': 'C', 
    'TTA': 'L', 'TTC': 'F', 'TTG': 'L', 'TTT': 'F'
  }

# vertebrate mitochondrial code (transl_table=2)
mitochondrial_code = {
    'AAA': 'K', 'AAC': 'N', 'AAG': 'K', 'AAT': 'N', 
    'ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACT': 'T', 
    'AGA': '*', 'AGC': 'S', 'AGG': '*', 'AGT': 'S', 
    'ATA': 'M', 'ATC': 'I', 'ATG': 'M', 'ATT': 'I', 
    'CAA': 'Q', 'CAC': 'H', 'CAG': 'Q', 'CAT': 'H', 
    'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCT': 'P', 
    'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGT': 'R', 
    'CTA': 'L', 'CTC': 'L', 'CTG': 'L', 'CTT': 'L', 
    'GAA': 'E', 'GAC': 'D', 'GAG': 'E', 'GAT': 'D', 
    'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCT': 'A', 
    'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGT': 'G', 
    'GTA': 'V', 'GTC': 'V', 'GTG': 'V', 'GTT': 'V', 
    'TAA': '*', 'TAC': 'Y', 'TAG': '*', 'TAT': 'Y', 
    'TCA': 'S', 'TCC': 'S', 'TCG': 'S', 'TCT': 'S', 
    'TGA': 'W', 'TGC': 'C', 'TGG': 'W', 'TGT': 'C', 
    'TTA': 'L', 'TTC': 'F', 'TTG': 'L', 'TTT': 'F'
  }

def create_codon_table(file_path):
    """
    The above pre-compiled codon table is standard genetic code.
    If your organism or organelle (e.g. mitochondria, plastid) uses
    a different genetic code, download the equivalent code from
    NCBI's genetic code website:
        https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi
    To find which genetic code you need to use, search your species
    through NCBI's taxonomy database:
        https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi

    Once you find the genetic code of interest, copy and paste the
    following format from the database and create a text file.

      AAs  = FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG
    Starts = ---M------**--*----M---------------M----------------------------
    Base1  = TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG
    Base2  = TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG
    Base3  = TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG
    """
    aa = []
    codon= []
    for c in open(file_path, 'r'):
        c = c.strip().replace(' ', '').split('=')
        if c[0].startswith('Starts'):
            continue
        elif c[0].startswith('AAs'):
            aa = list(c[1])
        elif c[0].startswith('Base'):
            codon.append(c[1])
            
    codon_table = dict()
    for i, c in enumerate(zip(*codon)):
        codon_table.setdefault(''.join(c), aa[i])
        
    return codon_table


complement = str.maketrans("ACGTRYKMBDHV", "TGCAYRMKVHDB")
class Seq:
    """A class holding a nucleotide sequence.
    """
    seq: str
    """The nucleotide sequence."""
    aa: Optional[str]
    """An amino acid sequence corresponding to the translated peptide sequence."""

    def __init__(self, seq: str):
        self.seq = seq.upper()
        self.aa = None
        
    def reverse_complement(self):
        self.seq = self.seq[::-1].translate(complement)
        return self.seq

    def translate(self, start: Optional[int]=None, end: Optional[int]=None, codon_table=standard_code, padchar='') -> Tuple[str, int]:
        """
        Translates the sequence from a nucleotide sequence to an amino acid sequence.

        Args:
            start (int, optional): The starting nucleotide. Defaults to 0.
            end (int, optional): The final nucleotide. Defaults to the end of the sequence.
            codon_table (int, optional): A codon table used to translate nucleotides to amino acids. Defaults to standard_code.
            padchar (str, optional): _description_. Defaults to ''.

        Returns:
            A tuple, containing the translated sequence and the number of untranslated nucleotides remaining.
        """
        if start is None:
            start = 0
        if end is None:
            end = len(self.seq)
        peptide = ''
        i = start
        for i in range(start, end, 3):
            aa = '*'
            c = self.seq[i:min(i+3, end)]
            if len(c) != 3:
                break
            try:
                aa = codon_table[c]
            except KeyError:
                raise Exception('unknown codon code: {}'.format(c))
            peptide += padchar + aa + padchar
            if aa == '*':
                break
        self.aa = peptide
        return peptide, end-min(i+3,end)


def unpad_peptide(aa):
    return aa.replace(' ', '').rstrip('*')

