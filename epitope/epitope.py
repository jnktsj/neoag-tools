from collections import defaultdict
import os
import copy
from pathlib import Path
import pickle
import re
import logging
import itertools
import argparse
import sys
import pysam
from typing import (
    Any,
    Callable,
    Dict,
    Hashable,
    Iterable,
    List,
    Optional,
    TextIO,
    Tuple,
    TypeVar,
)

import numpy as np
import pandas as pd

from .seq import Mutation, Seq, standard_code, create_codon_table
from .io import (
    write_neoorf,
    write_peptide,
    write_fasta,
    read_rsem_gene,
    read_maf,
    neoorf_header,
    common_header,
    neoorf_mut_class,
    peptide_header,
    coding_mut_class,
)
from .gtf import Annotation, Transcript, Locus

VALID_PEPTIDE_LENGTHS = [8, 9, 10, 11]


def validate_path(name: str, path: Optional[str], is_required=True):
    """
    Checks that the path is present.

    Args:
        name (str): The name of the object pointed to by the path.
        path (str): A path to a necessary file.
    """
    if path is None:
        if is_required:
            raise Exception(f"Provide {name}")
    else:
        if not os.path.exists(path):
            raise FileNotFoundError(f"{name} not found: {path}")


def check_input_args(args: argparse.Namespace):
    """
    Check input parameters
    """
    # check required arguments
    validate_path("reference FASTA", args.ref_fasta)
    validate_path("input MAF", args.input_maf)
    validate_path("gene annotation GTF", args.gtf)
    validate_path("codon table", args.codon_table, is_required=False)
    validate_path("gene expression matrix", args.gene_expr, is_required=False)

    if tumor_name is None:
        raise Exception("provide tumor sample name in the input MAF")

    # check optional arguments
    if normal_name or args.germline_maf:
        if normal_name is None:
            raise Exception("provide normal sample name")
        validate_path("germline MAF", args.germline_maf)

    if args.flank_peplen < 0:
        raise ValueError(
            "invalid peptide upstream/downstream length; provide a value larger than or equal to 0"
        )

    if args.peptide_length == "":
        raise ValueError("provide peptide length(s)")
    for p in args.peptide_length.split(","):
        try:
            pep_len = int(p)
            if pep_len not in VALID_PEPTIDE_LENGTHS:
                raise ValueError(
                    f"invalid peptide length: {p} ({args.peptide_length})"
                )
        except ValueError:
            raise ValueError(
                f"invalid peptide length: {p} ({args.peptide_length})"
            )


CONTIG_COL = "Chromosome"
POS_COL = "Start_position"
REF_COL = "Reference_Allele"
ALT_COL = "Tumor_Seq_Allele"
TRANSCRIPT_ID_COL = "Annotation_Transcript"
PHASE_ID_COL = "PhaseID"
CLONAL_STRUCTURE_COL = "ClonalStructure"
GERMLINE_ID_COL = "GermlineID"


def parse_mutations(mutations_df: pd.DataFrame) -> List[Mutation]:
    """Parses the given dataframe containing mutation data to a list of mutation objects.

    Args:
        mutations_df (pd.DataFrame): A pandas DataFrame containing the mutation data.

    Returns:
        list[Mutation]: A list of the mutations.
    """

    res = []
    for _, row in mutations_df.iterrows():
        phase_id = None if pd.isna(row[PHASE_ID_COL]) else row[PHASE_ID_COL]
        clonal_structure = (
            None
            if pd.isna(row[CLONAL_STRUCTURE_COL])
            else row[CLONAL_STRUCTURE_COL]
        )
        is_germline = pd.isna(row[GERMLINE_ID_COL])

        res.append(
            Mutation(
                row[CONTIG_COL],
                row[POS_COL],
                row[REF_COL],
                row[ALT_COL],
                row[TRANSCRIPT_ID_COL],
                phase_id,
                clonal_structure,
                is_germline,
            )
        )

    return res


def group_mutations(muts: Iterable[Mutation]) -> List[List[Mutation]]:
    """
    Groups the mutations by their phase ID.

    Args:
        muts (Iterable[Mutation]): An iterable containing mutations.

    Returns:
        List[List[Mutation]]: A list of lists of mutations, each list containing mutations with the same phase ID. Mutations with a None
    """
    res = []
    by_phase_id = defaultdict(list)

    for mut in muts:
        if mut.phase_id is None:
            res.append([mut])
        else:
            by_phase_id[mut.phase_id].append(mut)

    res.extend(by_phase_id.values())
    return res


def mutate_sequence(
    txid: str, seq: List[List[str]], pos: List[Locus], muts: pd.DataFrame
) -> Tuple[List[List[str]], List[List[int]]]:
    """
    Applies mutations to the given sequence.

    Args:
        txid (str): The trancript ID.
        seq (List[List[str]]): The sequence, represented as a list of the coding regions and the UTR.
        pos (List[Locus]): The loci of the coding regions and the UTRs.
        muts (pd.DataFrame): A dataframe containing the mutation data.

    Returns:
        Tuple: A tuple, containing the mutated sequence in the same format as the input sequence, and the mutation IDs, represented as a list of lists of ints,
        where each int represents the mutation ID that affected the corresponding bases.
    """
    mut_seq = copy.deepcopy(seq)
    idx_seq = [[-1 for _ in s] for s in mut_seq]

    x: int
    for x, m in muts.iterrows():  # type: ignore
        ref = m["Reference_Allele"]
        alt = m["Tumor_Seq_Allele"]
        # Finding `index`, the index of the locus where the mutation lies, and `start`, the start position of the mutation in the sequence.
        for index in range(len(pos)):
            if (
                pos[index].start <= m["Start_position"]
                and m["Start_position"] <= pos[index].end
            ):
                start = m["Start_position"] - pos[index].start
                break
        else:
            if m["Variant_Classification"] in coding_mut_class:
                logging.warning(
                    "skipping {}:{}; position not in {}".format(
                        m["Chromosome"], m["Start_position"], txid
                    )
                )
            continue
        if m["Variant_Type"] == "DEL":
            end = start + len(ref)
            mut_seq[index][start:end] = [""] * len(ref)
            idx_seq[index][start:end] = [x] * len(ref)
        elif m["Variant_Type"] == "INS":
            mut_seq[index][start] += alt
            idx_seq[index][start] = x
        elif len(ref) == len(alt):
            end = start + len(alt)
            for i, mb in zip(range(start, end), list(alt)):
                if mut_seq[index][i] == "":
                    continue
                mut_seq[index][i] = mb + mut_seq[index][i][1:]
                idx_seq[index][i] = x
        else:
            logging.warning("unknown variant type: " + m["Variant_Type"])
    return mut_seq, idx_seq


def get_mutation_notation(
    wt: List[str], mt: List[str], idx: List[int]
) -> Tuple[str, Dict[int, Tuple[int, int]]]:
    """
    Gets a string describing the differences between the WT and the mutated sequences.
    Insertions are mapped as 'i', deletions as 'd', mutations as 'm', and matches as ' '.

    Args:
        wt (List[str]): The wild-type sequence, represented as aligned chunks.
        mt (List[str]): The mutated sequence, represented as aligned chunks.
        idx (List[int]): A list of the mutation IDs of every position, or -1 if the position is not mutated.
    Returns:
        Tuple: A tuple containing:
          * The mutation representation
          * A map from the mutation IDs to the mutated spans.

    The spanning regions for each mutation are defined as the minimal bounding sets for the mutation.
    This is meant to handle multi-nucleotide substitutions, which can produce non-contiguous mutated regions but are sometimes defined as a single mutation
    by the upstream code.
    """
    mt_str = ""
    mt_pos = []
    mt_idx = {}

    assert len(wt) == len(mt) == len(idx)
    for wt_base, mt_base, mut_idx in zip(wt, mt, idx):
        # Handling unmutated regions.
        if wt_base == mt_base:
            mt_str += " " * len(wt_base)
            mt_pos += [-1] * len(wt_base)
            continue

        mt_pos += [mut_idx] * len(wt_base)
        # Deletions
        if mt_base == "":
            mt_idx.setdefault(mut_idx, [len(mt_str), len(mt_str)])[1] = len(
                mt_str
            ) + len(wt_base)
            mt_str += "d" * len(wt_base)
        # Insertions
        elif len(mt_base) > len(wt_base):
            mt_idx.setdefault(mut_idx, [len(mt_str), len(mt_str)])[1] = len(
                mt_str
            ) + len(mt_base)
            mt_str += " " * len(wt_base) + "i" * (len(mt_base) - len(wt_base))
        else:
            assert len(wt_base) == len(mt_base), (
                'Deletions should be represented as matches to ""'
            )
            mt_idx.setdefault(mut_idx, [len(mt_str), len(mt_str)])[1] = len(
                mt_str
            ) + len(mt_base)
            mt_str += "m" * len(mt_base)

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
    return "".join("*" if w != m else " " for w, m in zip(wt, mt))
    aa_str = ""
    for i in range(min(len(wt), len(mt))):
        if wt[i] != mt[i]:
            aa_str += "*"
        else:
            aa_str += " "
    return aa_str


MUT_IDX = 0


def translate_mutation(
    fo_fasta: TextIO,
    fo_peptide: TextIO,
    fo_neoorf: TextIO,
    tumor_name: str,
    smuts: Dict[str, Dict[str, Any]],
    normal_name: str,
    gmuts: pd.DataFrame,
    gtf: Annotation,
    gene_tpm: Dict[str, str],
    flank_length: int,
    peptide_lengths: List[int],
    genome: str,
    codon_table: Dict[str, str],
):
    global MUT_IDX
    for cl in smuts:
        for txid in smuts[cl]["tx"]:
            MUT_IDX += 1
            if txid not in gtf.transcripts:
                logging.warning(
                    f"Skipping the transcript {txid}; not in the GTF file"
                )
                continue
            tx: Transcript = gtf.transcripts[txid]
            if not tx.is_coding:
                logging.warning(
                    f"Skipping the transcript {txid}; CDS not defined"
                )
                continue

            seq, pos, exon_str = tx.raw_cds_with_downstream(genome)
            # get coding mutation index
            vc_i = smuts[cl]["mut"]["Variant_Classification"].isin(
                coding_mut_class
            )
            at_i = smuts[cl]["mut"]["Annotation_Transcript"] == txid
            mx = list(smuts[cl]["mut"][vc_i & at_i].index)

            # append germline mutations
            gx = []
            if not gmuts.empty:
                gv_i = gmuts["Variant_Classification"].isin(coding_mut_class)
                gt_i = gmuts["Annotation_Transcript"] == txid
                gx = list(gmuts[gv_i & gt_i].index)
                seq, _ = mutate_sequence(txid, seq, pos, gmuts)

            # append somatic mutations
            mseq, midx = mutate_sequence(txid, seq, pos, smuts[cl]["mut"])

            # remove exon structure
            seq = list(itertools.chain(*seq))
            mseq = list(itertools.chain(*mseq))
            midx = list(itertools.chain(*midx))

            wt_seq = Seq("".join(seq))
            mt_seq = Seq("".join(mseq))

            if tx.pos.strand == "-":
                wt_seq = wt_seq.reverse_complement()
                mt_seq = mt_seq.reverse_complement()
                exon_str = exon_str[::-1]
                midx = midx[::-1]
                mx = mx[::-1]
                gx = gx[::-1]
                # note: these sequences are not reverse complement, just reverse.
                #       but it's sufficient for mutation notation.
                seq = seq[::-1]
                mseq = mseq[::-1]

            wt_aa, _ = wt_seq.translate(codon_table=codon_table, padchar=" ")
            mt_aa, _ = mt_seq.translate(codon_table=codon_table, padchar=" ")

            mt_str, mt_idx = get_mutation_notation(seq, mseq, midx)
            aa_str = get_peptide_notation(wt_aa, mt_aa)
            # print(f'Writing mutation {MUT_IDX}: {wt_aa} {mt_aa}')
            write_fasta(
                fo_fasta,
                tumor_name,
                normal_name,
                smuts[cl]["mut"].loc[mx],
                gmuts.loc[gx],
                txid,
                wt_seq,
                mt_seq,
                aa_str,
                mt_str,
                exon_str,
                cl,
            )

            tpm = gene_tpm.get(txid.split(".")[0], "nan")

            write_peptide(
                fo_peptide,
                smuts[cl]["mut"].loc[mx],
                cl,
                wt_seq,
                mt_seq,
                aa_str,
                mt_idx,
                tumor_name,
                txid,
                tx.gene.id,
                tpm,
                flank_length,
                peptide_lengths,
            )

            nf_i = (
                smuts[cl]["mut"]
                .loc[mx, "Variant_Classification"]
                .isin(neoorf_mut_class)
            )
            mx = list(smuts[cl]["mut"][vc_i & at_i & nf_i].index)
            if len(mx) > 0:
                write_neoorf(
                    fo_neoorf,
                    smuts[cl]["mut"].loc[mx],
                    cl,
                    wt_seq,
                    mt_seq,
                    aa_str,
                    mt_idx,
                    tumor_name,
                    txid,
                    tx.gene.id,
                    tpm,
                )


# def mutate_sequence(transcript: Transcript, mutations: List[Mutation], reference: pysam.FastaFile) -> Tuple[str, str, List[int]]:
#     """
#     Computes the mutated nucleotide sequence.

#     Args:
#         seq (Seq): The sequence to mutate.
#         mutations (List[Mutation]): The list of mutations to apply.
#         reference (pysam.FastaFile): A fasta file containing the reference genome.

#     Returns:

#     A tuple containing:
#     * The unmutated sequence.
#     * The mutated sequence.
#     * A list of integers, mapping each nucleotide in the muated sequence to a nucleotide in the unmutated sequence.
#     """
#     wt = []
#     mt = []
#     map_indices = []

#     for exon_idx in transcript.get_exon_order():
#         exon = transcript.exons[exon_idx]
#         exon.cds_pos


def run(args):
    check_input_args(args)

    logging.basicConfig(
        level=logging.INFO,
        format="epitope | %(asctime)s [%(levelname)s] %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )

    # load mutations
    logging.info("loading mutations from input MAF")

    tumor_name = args.tumor_name
    normal_name = args.normal_name

    mutations = read_maf(args.input_maf, tumor_name, "Tumor_Sample_Barcode")
    germline_mutations = read_maf(
        args.germline_maf, normal_name, "Matched_Norm_Sample_Barcode"
    )

    # If the names are comma-separated lists, take the first element as the chosen name.
    # (Otherwise, they will remain the same)
    tumor_name = tumor_name.split(",")[0]
    normal_name = normal_name.split(",")[0]

    # Ignoring the phasing if we are instructed to.
    if args.ignore_phasing:
        mutations[PHASE_ID_COL] = np.nan
        germline_mutations[PHASE_ID_COL] = np.nan
    else:
        if PHASE_ID_COL not in mutations:
            raise Exception('cannot find "PhaseID" in: ' + args.input_maf)
        if PHASE_ID_COL not in germline_mutations:
            germline_mutations[PHASE_ID_COL] = np.nan

    # Ignoring the clonal structure if we are instructed to.
    if args.ignore_clones:
        mutations[CLONAL_STRUCTURE_COL] = np.nan
        germline_mutations[CLONAL_STRUCTURE_COL] = np.nan
    else:
        if CLONAL_STRUCTURE_COL not in mutations:
            raise Exception(
                'cannot find "ClonalStructure" in: ' + args.input_maf
            )
        if CLONAL_STRUCTURE_COL not in germline_mutations:
            germline_mutations[CLONAL_STRUCTURE_COL] = np.nan

    # Making sure we have the germline annotation.
    if not germline_mutations.empty:
        if "GermlineID" not in mutations:
            raise Exception('cannot find "GermlineID" in: ' + args.input_maf)
        if "GermlineID" not in germline_mutations:
            raise Exception('cannot find "GermlineID" in: ' + args.germline_maf)

    transcripts = list(mutations["Annotation_Transcript"].unique())
    # load transcript structures from GTF
    logging.info(
        "loading {} transcripts with mutations from GTF".format(
            len(transcripts)
        )
    )

    DEBUG_CACHE = Path("cache.pickle")
    if DEBUG_CACHE.exists():
        gtf = pickle.load(DEBUG_CACHE.open("rb"))
    else:
        gtf = Annotation(args.gtf, transcript_list=transcripts)
        pickle.dump(gtf, DEBUG_CACHE.open("wb"))

    logging.info(
        "finished loading {} transcripts from GTF".format(len(gtf.transcripts))
    )

    # exit if there are no coding transcripts with return code 0
    is_coding = mutations["Variant_Classification"].isin(coding_mut_class)
    if sum(is_coding) == 0:
        logging.warning("no coding mutations found! exiting...")
        sys.exit(0)

    coding_mutations = mutations.loc[is_coding]
    logging.info(f"Translating {len(coding_mutations)} coding mutations.")

    # load codon table if supplied
    codon_table = standard_code
    if args.codon_table:
        logging.info("custom codon table supplied: {}".format(args.codon_table))
        codon_table = create_codon_table(args.codon_table)

    gene_tpm = dict()
    # load RSEM gene matrix if supplied
    if args.gene_expr:
        logging.info("loading gene expression data")
        gene_tpm = read_rsem_gene(args.gene_expr, transcript_list=transcripts)

    peptide_lengths = list(map(int, args.peptide_length.split(",")))

    # Parsing the mutations.
    tumor_muts = parse_mutations(mutations)
    germline_muts = parse_mutations(germline_mutations)

    # output FASTA
    fo_fasta = open(tumor_name + ".muts.fa", "w")

    # output peptide file
    fo_peptide = open(tumor_name + ".peptide.tsv", "w")
    fo_peptide.write("\t".join(peptide_header + common_header) + "\n")

    # output neoorf file
    fo_neoorf = open(tumor_name + ".neoorf.tsv", "w")
    fo_neoorf.write("\t".join(neoorf_header + common_header) + "\n")

    logging.info("start translating mutations")
    analyzed = set()

    for i, mut in mutations.loc[is_coding].iterrows():
        index: List[Any] = [i]
        if not pd.isna(mut["PhaseID"]):
            curr_phase_id = mut["PhaseID"]  # noqa: F841  # This is used two lines below with a pandas query.
            index.extend(
                list(coding_mutations.query("PhaseID == @curr_phase_id").index)
            )
            index = sorted(set(index))
            analyzed.update(index)

        phased_coding_mutations = (
            coding_mutations.loc[index].reset_index(drop=True).copy()
        )

        # subset somatic mutations
        smuts = dict()
        clones = list(phased_coding_mutations["ClonalStructure"].dropna())
        if len(clones) > 0:
            for cl in sorted(set.union(*clones)):
                cs = phased_coding_mutations["ClonalStructure"].apply(
                    lambda clstr: cl in clstr
                )
                xi = phased_coding_mutations.loc[cs].copy()
                # avoid a clone composed of all non-coding mutations
                xi = xi.sort_values("Start_position", ignore_index=True)
                is_coding_muts = xi["Variant_Classification"].isin(
                    coding_mut_class
                )
                if sum(is_coding_muts) == 0:
                    print("Skipped because non-coding")
                    continue
                tx = set(list(xi.loc[is_coding_muts, "Annotation_Transcript"]))
                smuts.setdefault(cl, {"mut": xi, "tx": tx})
        else:
            is_coding_muts = phased_coding_mutations[
                "Variant_Classification"
            ].isin(coding_mut_class)
            tx = set(
                list(
                    phased_coding_mutations.loc[
                        is_coding_muts, "Annotation_Transcript"
                    ]
                )
            )
            smuts.setdefault("0", {"mut": phased_coding_mutations, "tx": tx})

        # subset germline mutations if supplied
        gmuts = pd.DataFrame()
        if not germline_mutations.empty:
            gid = list(phased_coding_mutations["GermlineID"].dropna())
            if len(gid) > 0:
                gmuts = germline_mutations[
                    germline_mutations["GermlineID"].isin(set(gid))
                ].copy()

        translate_mutation(
            fo_fasta,
            fo_peptide,
            fo_neoorf,
            tumor_name,
            smuts,
            normal_name,
            gmuts,
            gtf,
            gene_tpm,
            args.flank_peplen,
            peptide_lengths,
            args.ref_fasta,
            codon_table,
        )
    # close file handles and end
    fo_fasta.close()
    fo_peptide.close()
    fo_neoorf.close()
    logging.info("all done!")
