import os
from pathlib import Path
import re
import gzip
from collections import OrderedDict
from typing import Any, Dict, List, Optional, TextIO, Tuple, Union

import numpy as np
import pandas as pd

from .seq import Seq, unpad_peptide


# coding mutations
# excluding 'Splice_Site' from amino acid translation
coding_mut_class = [
    "Missense_Mutation",
    "Frame_Shift_Ins",
    "Frame_Shift_Del",
    "Nonstop_Mutation",
    "In_Frame_Ins",
    "In_Frame_Del",
]
neoorf_mut_class = ["Frame_Shift_Ins", "Frame_Shift_Del", "Nonstop_Mutation"]

# maf columns
required_cols = OrderedDict(
    {
        "Hugo_Symbol": str,
        "Chromosome": str,
        "Start_position": int,
        "End_position": int,
        "Reference_Allele": str,
        "Tumor_Seq_Allele": str,
        "Variant_Classification": str,
        "Variant_Type": str,
        "Genome_Change": str,
        "Annotation_Transcript": str,
        "cDNA_Change": str,
        "Codon_Change": str,
        "Protein_Change": str,
    }
)
optional_cols = OrderedDict(
    {
        "ClonalStructure": np.dtype("object"),
        "PhaseID": float,
        "GermlineID": float,
    }
)

# peptide and neoorf file columns
peptide_header = [
    "Hugo_Symbol",
    "pep_wt",
    "pep",
    "ctex_up",
    "ctex_dn",
    "pep_start",
    "pep_end",
]

neoorf_header = ["Hugo_Symbol", "upstream", "neoORF", "neoORF_start"]

common_header = [
    "transcript_id",
    "gene_id",
    "TPM",
    "Tumor_Sample_Barcode",
    "clone",
    "Chromosome",
    "Start_position",
    "End_position",
    "Variant_Classification",
    "Variant_Type",
    "Genome_Change",
    "cDNA_Change",
    "Codon_Change",
    "Protein_Change",
]


def read_maf(
    maf: Union[str, Path, None], name: str, name_col: str
) -> pd.DataFrame:
    """
    Reads a MAF file into a pandas DataFrame.
    The MAF files handled both use # as a way to mark comments and as a character inside fields,
    which would confuse pd.read_csv.

    Args:
        maf (Union[str, Path, None]): The path to the MAF file.
        name (str): The name of the sample.
        name_col (str): The name of the column where names are stored.

    Returns:
        pd.DataFrame: A DataFrame holding the MAF data.
    """
    if maf is None:
        return pd.DataFrame()

    # reading MAF in a primitive way to handle command characters '#'
    # in the middle of the line which pd.read_csv couldn't handle
    header = []
    record = []
    for line in open(maf):
        if line.startswith("#"):
            continue
        line = line.rstrip("\n").split("\t")
        if line[0] == "Hugo_Symbol":
            header = line
        else:
            record.append({h: x for h, x in zip(header, line)})
    df = pd.DataFrame(data=record, dtype=str)

    tid = set(df[name_col])
    if name not in tid:
        raise Exception(name + " not in " + os.path.basename(maf))
    if len(tid) > 1:
        row = df[name_col] == name
        df[row].reset_index(drop=True, inplace=True)
    tid = name

    # rename MAF header
    df.rename(columns={"Tumor_Seq_Allele2": "Tumor_Seq_Allele"}, inplace=True)
    if "Start_Position" in df:
        df.rename(columns={"Start_Position": "Start_position"}, inplace=True)
    if "End_Position" in df:
        df.rename(columns={"End_Position": "End_position"}, inplace=True)

    # format MAF columns
    if "PhaseID" in df:
        df.loc[df["PhaseID"].isin(["nan", ""]), "PhaseID"] = np.nan
    if "GermlineID" in df:
        df.loc[df["GermlineID"].isin(["nan", ""]), "GermlineID"] = np.nan

    # subset MAF columns
    maf_cols: OrderedDict[str, Any] = required_cols.copy()
    for oc in optional_cols:
        if oc in df:
            maf_cols[oc] = optional_cols[oc]
    try:
        df = df[maf_cols.keys()].astype(maf_cols)
    except Exception as e:
        raise KeyError(str(e).replace("index", "MAF"))

    # process clone information if supplied
    if "ClonalStructure" in df:
        df.loc[df["ClonalStructure"].isin(["nan", ""]), "ClonalStructure"] = (
            np.nan
        )

        i = ~df["ClonalStructure"].isna()

        def make_set(x):
            """Converts a string of the form "[a, b, c]" to a set {a, b, c}."""
            return set(
                [
                    xi.strip()
                    for xi in re.sub(r"\[|\]", "", x).split(",")
                    if xi.strip().isnumeric()
                ]
            )

        df.loc[i, "ClonalStructure"] = df.loc[i, "ClonalStructure"].apply(
            make_set
        )

    return df


def read_rsem_gene(rsem_gene: str, transcript_list=[]):
    """
    Returns Transcript ID (without versioning) and TPM
    from RSEM gene expression matrix
    """
    # strip off transcript ID versions
    if transcript_list != []:
        transcript_list = set([txid.split(".")[0] for txid in transcript_list])

    gene_tpm = dict()
    try:
        if rsem_gene.endswith(".gz"):
            rsem_gene_fi = gzip.open(rsem_gene, "rt")
        else:
            rsem_gene_fi = open(rsem_gene, "r")
    except IOError:
        raise Exception("cannot read: {}".format(rsem_gene))

    for row in rsem_gene_fi:
        # skip comment lines if any
        if row.startswith("#"):
            continue
        row = row.rstrip("\n").split("\t")
        if row[0] == "gene_id":
            continue

        transcript_ids = set([txid.split(".")[0] for txid in row[1].split(",")])
        tpm = row[5]

        if transcript_list != []:
            overlap = transcript_list.intersection(transcript_ids)
            if len(overlap) == 0:
                continue
            transcript_ids = overlap

        for txid in transcript_ids:
            gene_tpm.setdefault(txid, tpm)

    return gene_tpm


def get_fasta_header(
    m: pd.DataFrame, name: str, txid: str, clone: Optional[str] = None
):
    gene = m["Hugo_Symbol"].unique()[0]
    nt = ";".join(m["cDNA_Change"].unique())
    aa = ";".join(m["Protein_Change"].unique())
    if clone is not None:
        return "|".join([name, txid, gene, nt, aa, "clone=" + clone])
    else:
        return "|".join([name, txid, gene, nt, aa])


def write_fasta(
    fo: TextIO,
    tumor_name: str,
    normal_name: str,
    muts: pd.DataFrame,
    gmuts: pd.DataFrame,
    txid: str,
    wt: Seq,
    mt: Seq,
    aa_str: str,
    mt_str: str,
    exon_str: str,
    clone: str,
):
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
    assert wt.aa is not None, (
        "Wild-type sequence was not translated before writing to the FASTA output!"
    )
    assert mt.aa is not None, (
        "Mutated sequence was not translated before writing to the FASTA output!"
    )

    tumor_field = get_fasta_header(muts, tumor_name, txid, clone)
    normal_field = ""
    if not gmuts.empty:
        normal_field = get_fasta_header(gmuts, normal_name, txid)
    header = ">" + tumor_field + ";" + normal_field
    fo.write(header + "\n")
    fo.write(wt.aa + "\n")
    fo.write(aa_str + "\n")
    fo.write(mt.aa + "\n")
    fo.write(mt.seq + "\n")
    fo.write(mt_str + "\n")
    fo.write(wt.seq + "\n")
    fo.write(exon_str + "\n")


def write_peptide(
    fo: TextIO,
    smuts: pd.DataFrame,
    clone: str,
    wt: Seq,
    mt: Seq,
    mstr: str,
    idx: Dict[int, Tuple[int, int]],
    name: str,
    txid: str,
    gnid: str,
    tpm: str,
    flen: int,
    pep_lens: List[int],
):
    """
    Write the given peptide to the output file.

    Args:
        fo (TextIO): An open pipe, writing to the output file.
        smuts (pd.DataFrame): A dataframe containing a set of mutations.
        clone (str): The clone from the clonal structure of the mutations.
        wt (Seq): The wild-type sequence.
        mt (Seq): The mutated sequence.
        mstr (str): The representation of the differences between the wild-type and the mutation, produced by `get_peptide_notation`.
        idx (Dict[int, Tuple[int, int]]): A map from some indices to spans of mutated regions in the mutated sequence. These are produced by `get_mutation_notation`.
            The index only really has to be correct for the mutated sequence, and the wild-type peptide is neither really well defined nor handled properly:
            If there are two deletions in the same peptide, the first will induce a shift in the alignment of the next two.
        name (str): The tumor name.
        txid (str): ???
        gnid (str): The ID of the gene.
        tpm (str): The expression level of the transcript in TPM?
        flen (int): The length of the flanking regions.
        pep_lens (List[int]): The lengths of peptides to take from the transcript.
    """
    assert wt.aa is not None
    assert mt.aa is not None

    flank_length = flen + max(pep_lens) - 1

    index: int
    for index, m in smuts.iterrows():  # type: ignore
        start, end = idx[index]
        start = start - start % 3
        end = (end + 2) - (end + 2) % 3

        indices = range(start, min(len(mstr), end), 3)
        # The mutated range can theoretically be in a non-coding sequence, in which case we need to ignore it as it has no transcript.
        if len(indices) == 0:
            print(
                f"Skipped: Mutation is in a non-coding sequence. {start=} {end=}"
            )
            continue

        # adjust reading frame for shifted amino acid mismatches
        for i in indices:
            if mstr[i] == "*":
                start = i
                break

        for i in reversed(indices):
            if mstr[i] == "*":
                end = i

        if m["Variant_Classification"].startswith("Frame_Shift"):
            wt_aa = unpad_peptide(wt.aa[start:])
            mt_aa = unpad_peptide(mt.aa[start:])
            # make wild type sequence same length as mutated sequence
            wt_aa += "-" * (len(mt_aa) - len(wt_aa))
            wt_aa = wt_aa[: len(mt_aa)]
            mt_aa_dn = ""
            wt_aa_dn = ""
        elif m["Variant_Classification"].startswith("Nonstop_Mutation"):
            mt_aa = unpad_peptide(mt.aa[start:])
            wt_aa = "-" * len(mt_aa)
            mt_aa_dn = ""
            wt_aa_dn = ""
        else:
            wt_aa = unpad_peptide(wt.aa[start:end])
            mt_aa = unpad_peptide(mt.aa[start:end])
            mt_aa_dn = unpad_peptide(mt.aa[end:])
            wt_aa_dn = unpad_peptide(wt.aa[end:])

            # trim downstream sequence with specified flanking peptide length
            mt_aa_dn = mt_aa_dn[:flank_length] + "-" * (
                flank_length - len(mt_aa_dn)
            )
            wt_aa_dn = wt_aa_dn[:flank_length] + "-" * (
                flank_length - len(wt_aa_dn)
            )

        mt_aa_up = unpad_peptide(mt.aa[:start])
        wt_aa_up = unpad_peptide(wt.aa[:start])

        # mutated protein sequence coordinates
        aa_start = len(mt_aa_up) + 1
        aa_stop = aa_start + len(mt_aa) - 1

        # trim upstream sequence with specified flanking peptide length
        mt_aa_up = (
            "-" * (flank_length - len(mt_aa_up)) + mt_aa_up[-flank_length:]
        )
        wt_aa_up = (
            "-" * (flank_length - len(wt_aa_up)) + wt_aa_up[-flank_length:]
        )

        start = len(mt_aa_up)
        stop = start + len(mt_aa)

        wt_aa = wt_aa_up + wt_aa + wt_aa_dn
        mt_aa = mt_aa_up + mt_aa + mt_aa_dn

        for plen in pep_lens:
            for i in range(start - plen + 1, stop):
                pep_wt = wt_aa[i : i + plen]
                pep = mt_aa[i : i + plen]
                if len(pep) < plen:
                    continue
                elif "-" in pep:
                    continue
                elif pep == pep_wt:
                    continue
                ctex_up = mt_aa[i - flen : i]
                ctex_dn = mt_aa[i + plen : i + plen + flen]
                ctex_dn += "-" * (flen - len(ctex_dn))
                fo.write(
                    "\t".join(
                        [
                            m["Hugo_Symbol"],
                            pep_wt,
                            pep,
                            ctex_up,
                            ctex_dn,
                            str(aa_start),
                            str(aa_stop),
                            txid,
                            gnid,
                            tpm,
                            name,
                            clone,
                            m["Chromosome"],
                            str(m["Start_position"]),
                            str(m["End_position"]),
                            m["Variant_Classification"],
                            m["Variant_Type"],
                            m["Genome_Change"],
                            m["cDNA_Change"],
                            m["Codon_Change"],
                            m["Protein_Change"],
                        ]
                    )
                    + "\n"
                )


def write_neoorf(
    fo: TextIO,
    smuts: pd.DataFrame,
    clone: str,
    wt: Seq,
    mt: Seq,
    mstr: str,
    idx: Dict[int, Tuple[int, int]],
    name: str,
    txid: str,
    gnid: str,
    tpm: str,
):
    index: int
    for index, m in smuts.iterrows():  # type: ignore
        start, end = idx[index]
        start = start - start % 3 + 1
        end = (end - 1) - (end - 1) % 3 + 1

        if m["Variant_Classification"].startswith("Frame_Shift"):
            # adjust reading frame for shifted amino acid mismatches
            for i in range(start, len(mstr), 3):
                if mstr[i] == "*":
                    start = i
                    break

        assert mt.aa is not None
        neoorf = unpad_peptide(mt.aa[start:])
        upstream = unpad_peptide(mt.aa[:start])
        neoorf_start = len(upstream) + 1

        fo.write(
            "\t".join(
                [
                    m["Hugo_Symbol"],
                    upstream,
                    neoorf,
                    str(neoorf_start),
                    txid,
                    gnid,
                    tpm,
                    name,
                    clone,
                    m["Chromosome"],
                    str(m["Start_position"]),
                    str(m["End_position"]),
                    m["Variant_Classification"],
                    m["Variant_Type"],
                    m["Genome_Change"],
                    m["cDNA_Change"],
                    m["Codon_Change"],
                    m["Protein_Change"],
                ]
            )
            + "\n"
        )
