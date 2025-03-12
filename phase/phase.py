import os
import logging
import itertools
from collections import OrderedDict
from pathlib import Path
from typing import Dict, List, Optional, Set, Tuple, Union

import numpy as np
import pandas as pd
import pysam
from scipy.stats import fisher_exact

from .sample import Sample
from .mutation import Mutation
from .mutation import ClonalStructure
from .writer import (
    write_phase_vcf,
    write_phase_vcf_from_scratch,
    write_phase_maflite,
    write_subset_maf,
)


# make the ID variables float for handling these with 'nan'
PHASE_ID = 0.0
GERMLINE_ID = 0.0

MNP_TYPE = {1: "SNP", 2: "DNP", 3: "TNP"}


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


def check_input_args(args):
    """
    Check input parameters
    """
    # check required arguments
    validate_path("reference FASTA", args.ref_fasta)
    validate_path("input MAF", args.input_maf)
    validate_path("tumor BAM", args.tumor_bam)
    validate_path("tumor VCF", args.tumor_vcf, is_required=False)

    if args.tumor_name is None:
        raise Exception("provide tumor sample name in the input MAF")

    if args.normal_name or args.normal_bam or args.normal_vcf:
        if args.normal_name is None:
            raise Exception("provide normal sample name")

        validate_path("normal VCF", args.normal_vcf)
        validate_path("normal BAM", args.normal_bam)

    if args.phylogic_ccfs or args.phylogic_tree:
        validate_path("phylogic CCF", args.phylogic_ccfs)
        validate_path("phylogic tree", args.phylogic_tree)

        if args.phylogic_index <= 0:
            raise ValueError(
                "invalid tree index; provide a value larger than 0"
            )

    if args.max_mnp_size <= 3:
        raise ValueError("invalid maximum MNP; provide a value larger than 3")
    if args.max_phase_radius <= args.max_mnp_size:
        raise ValueError(
            "invalid max phasing radius; provide a value larger than max MNP length"
        )
    if args.vaf_skew_pval <= 0 or 1 <= args.vaf_skew_pval:
        raise ValueError(
            "invalid VAF skew p-value; a value nees to be 0 < p < 1"
        )
    if args.min_map_quality < 0:
        raise ValueError(
            "invalid min mapping quality; provide a value larger than or equal to 0"
        )
    if args.min_base_quality < 0:
        raise ValueError(
            "invalid min base quality; provide a value larger than or equal to 0"
        )
    if args.min_coverage < 0:
        raise ValueError(
            "invalid min coverage cutoff; provide a value larger than 0"
        )
    if args.min_phased_altc < 0:
        raise ValueError(
            "invalid min phased alt counts; provide a value larger than 0"
        )


def group_mutations(
    df: pd.DataFrame,
    max_dist: int,
    group_id: str,
    min_cov: int,
    excl_indel=False,
) -> pd.DataFrame:
    """
    Groups mutations within a certain distance, but excludes sites
    below a coverage cutoff from assigning group IDs
    Grouped mutations will have the same numerical IDs (start=1)

    The dataframe will be edited to have a "group ID" column. The group ID is an int/float that is the same for mutations in the same group.

    Args:
        df (pd.DataFrame): The data frame containing the mutations. Must include the columns ["Start_position", "Variant_Type", "t_alt_count", "t_ref_count"].
        max_dist (int): The maximum distance between mutations in the same group.
        group_id (str): The name of the group ID column that is created.
        min_cov (int): The minimal coverage of a mutation that can be included in a group.
        excl_indel (bool, optional): Whether to also group indel mutations, or only substitutions. Defaults to False.
    Returns:
        pd.DataFrame: The mutation dataframe, together with the column specifying the mutation grouping.
    """
    df[group_id] = np.full(len(df.index), np.nan)
    df.sort_values("Start_position", ignore_index=True, inplace=True)
    curr_id = 0
    prev_id = 0

    i: int
    dist: int
    for i, dist in list(df["Start_position"].diff().items())[1:]:  # type: ignore
        is_indel = sum(df.loc[i - 1 : i]["Variant_Type"].isin(["DEL", "INS"]))
        cov = min(
            df.loc[i - 1 : i]["t_alt_count"] + df.loc[i - 1 : i]["t_ref_count"]
        )
        if (excl_indel and is_indel) or (cov < min_cov):
            curr_id += 1
            continue
        if dist <= max_dist:
            curr_id = prev_id + 1 if curr_id > prev_id else curr_id
            df.loc[i - 1, group_id] = curr_id
            df.loc[i, group_id] = curr_id
            prev_id = curr_id
        else:
            curr_id += 1
    return df


def get_read_names(
    x: pd.Series,
    g: pysam.FastaFile,
    bam: pysam.AlignmentFile,
    min_base_q,
    min_map_q,
) -> Dict[str, Set[str]]:
    """Gets the names of the reads supporting the given mutations (?).

    Args:
        x (pd.Series): The row in the mutation dataframe describing the given mutation. Should contain the fields: ["Reference_Allele", "Tumor_Seq_Allele", "Variant_Type", "Chromosome", "Start_position"].
        g (pysam.FastaFile): The FASTA file of the human genome.
        bam (pysam.AlignmentFile): The BAM file containing the reads.
        min_base_q (_type_): The minimal base quality allowed in the reads.
        min_map_q (_type_): The minimal alignment quality allowed in the reads.

    Returns:
        Dict[str, Set[str]]: A dictionary of the form {"ref": [list of reference read names], "alt": [list of alternative read names]}
    """
    r = {"alt": set(), "ref": set()}
    if bam is None:
        return r
    ref = x["Reference_Allele"]
    alt = x["Tumor_Seq_Allele"]
    mut_type = x["Variant_Type"]
    chrom = x["Chromosome"]
    pos = x["Start_position"] - (0 if mut_type == "INS" else 1)
    end = x["Start_position"] + (1 if mut_type == "INS" else len(alt) - 1)
    for pileupcolumn in bam.pileup(
        chrom,
        pos,
        end,
        truncate=True,
        max_depth=1000000,
        min_base_quality=min_base_q,
        min_mapping_quality=min_map_q,
    ):
        for read in pileupcolumn.pileups:
            # filter reads with the following flags
            if (
                read.alignment.is_duplicate
                or read.alignment.is_qcfail
                or read.alignment.is_secondary
                or read.alignment.is_unmapped
                or read.alignment.is_supplementary
            ):
                continue
            rg = None
            base = None
            if read.query_position:
                base = read.alignment.seq[read.query_position]

            if mut_type in ["INS", "DEL"]:
                blocks = read.alignment.get_blocks()
                if len(blocks) > 1:
                    prev_block_end = blocks[0][1]
                    dist_to_mut = abs(prev_block_end - pos)
                    for block_start, block_end in blocks[1:]:
                        ref_space = block_start - prev_block_end
                        dist_to_mut = min(dist_to_mut, abs(block_start - pos))
                        # insertion
                        if (
                            ref_space == 0
                            and mut_type == "INS"
                            and pos == prev_block_end
                        ):
                            rg = "alt"
                            break
                        # deletion
                        elif (
                            ref_space == len(ref)
                            and mut_type == "DEL"
                            and prev_block_end <= pos <= block_start
                        ):
                            rg = "alt"
                            break
                        prev_block_end = block_end
                    if (
                        rg is None
                        and dist_to_mut > 1
                        and base == g.fetch(chrom, pos, end)
                    ):
                        rg = "ref"
                else:
                    if base == g.fetch(chrom, pos, end):
                        rg = "ref"
            else:  # point mutation
                if base is None:
                    continue
                offset = pileupcolumn.reference_pos - pos
                alt_base = alt[offset]
                ref_base = ref[offset]
                if base == alt_base:
                    rg = "alt"
                elif base == ref_base:
                    rg = "ref"
            if rg in r:
                r[rg].add(read.alignment.query_name)
    return r


def vaf_from_reads(r: Dict[str, Set[str]]) -> Tuple[int, int, float]:
    """Computes the variant allele frequency of the mutation.

    Args:
        r (Dict[str, Set[str]]): A dictionary of the form {"ref": [list of reference read names], "alt": [list of alternative read names]}

    Returns:
        Tuple[int, int, float]: A tuple containing the number of alternative reads, the number of reference reads, and the VAF.
    """
    alt = len(r["alt"])
    ref = len(r["ref"])
    try:
        vaf = alt / float(alt + ref)
        return alt, ref, vaf
    except ZeroDivisionError:
        return 0, 0, 0


def calc_phase_prob(
    r1: Dict[str, Set[str]], r2: Dict[str, Set[str]], min_phased_altc: float
) -> Tuple[float, float]:
    """Computes the P-Values of the VAFs of the mutations being the same, and of the mutations sharing phasing.

    Args:
        r1 (Dict[str, Set[str]]): A dictionary of the form {"ref": [list of reference read names], "alt": [list of alternative read names]}
        r2 (Dict[str, Set[str]]): A dictionary of the form {"ref": [list of reference read names], "alt": [list of alternative read names]}
        min_phased_altc (float): ???

    Returns:
        Tuple[float, float]: A tuple containing the P-values for the mutations having the same VAF and the mutations being phased together.
    """
    r1_alt, r1_ref, r1_vaf = vaf_from_reads(r1)
    r2_alt, r2_ref, r2_vaf = vaf_from_reads(r2)
    # check if VAFs of the two mutations are the same
    vaf_odds, vaf_pval = fisher_exact([[r1_alt, r1_ref], [r2_alt, r2_ref]])
    # raw alt and ref counts from a mutation with smaller VAF
    if r1_vaf > r2_vaf:
        u_alt = r2_alt
        u_ref = r2_ref
    else:
        u_alt = r1_alt
        u_ref = r1_ref
    # phased alt and ref counts
    p_alt = len(r1["alt"].intersection(r2["alt"]))
    p_ref = len(r1["ref"].intersection(r2["ref"]))
    phase_odds, phase_pval = fisher_exact([[u_alt, u_ref], [p_alt, p_ref]])
    if p_alt < min_phased_altc:
        phase_pval = 0.0
    if p_alt > min_phased_altc * 4:
        phase_pval = 1.0
    return (phase_pval, vaf_pval)  # type: ignore


def cluster_mutations(
    pos: Union[List[int], pd.Index],
    p: Dict[Tuple[int, int], float],
    cutoff: float,
) -> Set[Tuple[int, ...]]:
    """
    Clusters the mutations according to the given p-values.
    Given `pos`, a set of positions that serve as the nodes of the graph,
    and `p`, a dictionary mapping pairs of nodes to the p-values,
    computes for every node the set of neighbors whose edges have a p-value smaller than the cutoff.

    TODO: Ask someone about this clustering algorithm. It doesn't seem to do what one would expect.
    Clustering with large p-values is a bad idea in general, and could lead to a collapse when there is a mutation with very little data in general.
    This clustering really returns non-disjoint clusters.

    Args:
        pos (List[int]): The list of positions.
        p (Dict[Tuple[int, int], float]): A dictionary mapping pairs of positions to the p-value of them sharing a cluster.
        cutoff (float): The cutoff for an edge to be significant.

    Returns:
        Set[Tuple[int, ...]]: A set of tuples of mutations belonging in the same cluster.
    """
    # p-value is insignificant when variants are phased
    if len(pos) == 1:
        return set([(pos[0],)])
    gi: Dict[int, Set[int]] = dict()
    for pair in itertools.combinations(pos, 2):
        x, y = tuple(sorted(pair))
        hx, hy = (x,), (y,)
        if p[(x, y)] > cutoff:
            hx = hy = (x, y)
        gi[x] = gi.get(x, set()).union(hx)
        gi[y] = gi.get(y, set()).union(hy)

    return set([tuple(sorted(h)) for h in gi.values()])


def chain_mnps(cl, pos, df, mnp_dict, g, chrom, t_reads, n_reads):
    start, end = list(df.loc[[pos[0], pos[-1]], "Start_position"])
    posr = list(df.loc[pos, "Start_position"] - start)
    maf_idx = set()
    ref = list(g.fetch(chrom, start - 1, end))
    alt = np.array(ref, dtype=str)
    t_alt_count, t_ref_count = [], []
    n_alt_count, n_ref_count = [], []
    for i in range(len(pos)):
        maf_idx = maf_idx.union(df.loc[pos[i], "maf_idx"])
        # this can handle existing DNP, TNP, ONP
        alt_base = df.loc[pos[i], "Tumor_Seq_Allele"]
        alt[posr[i] : posr[i] + len(alt_base)] = list(alt_base)
        t_alt_count.append(t_reads[pos[i]]["alt"])
        t_ref_count.append(t_reads[pos[i]]["ref"])
        n_alt_count.append(n_reads[pos[i]]["alt"])
        n_ref_count.append(n_reads[pos[i]]["ref"])
    row = OrderedDict(df.loc[pos[0]])
    row["maf_idx"] = maf_idx
    row["ClonalStructure"] = set([cl])
    row["Reference_Allele"] = "".join(ref)
    row["Tumor_Seq_Allele"] = "".join(alt)
    row["Variant_Type"] = MNP_TYPE.get(len(alt), "MNP")
    key = tuple([row[k] for k in list(row.keys())[:5]])
    if key in mnp_dict:
        mnp_dict[key]["ClonalStructure"].add(cl)
        return mnp_dict
    if row["Variant_Type"] != "SNP":
        row["t_alt_count"] = len(set.intersection(*t_alt_count))
        row["t_ref_count"] = len(set.intersection(*t_ref_count))
        if "n_alt_count" in df:
            row["n_alt_count"] = len(set.intersection(*n_alt_count))
            row["n_ref_count"] = len(set.intersection(*n_ref_count))
    mnp_dict.setdefault(key, row)
    return mnp_dict


def phase_mutations(
    df: pd.DataFrame,
    chrom: str,
    t: Sample,
    n: Sample,
    genome: Union[str, Path],
    clonal_struct: Optional[ClonalStructure],
    vaf_skew_pval,
    min_phased_altc,
    min_base_q,
    min_map_q,
) -> pd.DataFrame:
    """
    Look at reads supporting candidate mutations to determine phasing

    When determining phasing, we check every pair of mutations, and run Fisher's exact test to determine if they have the same VAF.
    """
    global PHASE_ID
    df = df.reset_index(drop=True)
    g = pysam.FastaFile(genome)
    tbam = pysam.AlignmentFile(t.bam, "rb")
    nbam = pysam.AlignmentFile(n.bam, "rb") if n.bam else None
    # get alt and ref reads from tumor and normal samples
    t_reads: List[Dict[str, Set[str]]] = []
    n_reads: List[Dict[str, Set[str]]] = []
    for i, row in df.iterrows():
        t_reads.append(get_read_names(row, g, tbam, min_base_q, min_map_q))
        n_reads.append(get_read_names(row, g, nbam, min_base_q, min_map_q))
    # compute p-values: [0] phased reads, [1] allele fractions
    unp = 0
    phase_p = dict()
    vaf_p = dict()
    for pair in itertools.combinations(df.index, 2):
        x, y = tuple(sorted(pair))
        phase_p[(x, y)], vaf_p[(x, y)] = calc_phase_prob(
            t_reads[x], t_reads[y], min_phased_altc
        )
        if phase_p[(x, y)] <= vaf_skew_pval:
            unp += 1
    # return the original df if not phased
    if unp == len(phase_p):
        df[["PhaseID", "MNP"]] = np.nan
        return df

    mutcl = pd.Series([np.nan] * len(df.index))
    if clonal_struct is not None:
        mutcl = df["Cluster_Assignment"]
        # switch to VAF estimate if undefined clone clusters exist
        if clonal_struct.excluded_clone.intersection(set(mutcl.unique())):
            mutcl = pd.Series([np.nan] * len(df.index))

    # estimate clonal structure from VAF if the information is not available
    # We don't run phylogic here, and the estimated tree is just a line.
    if sum(mutcl.isna()) > 0:
        vaf_tiers = cluster_mutations(df.index, vaf_p, vaf_skew_pval)
        vaf_means = dict()
        for vt in vaf_tiers:
            alt = df.loc[list(vt), "t_alt_count"]
            ref = df.loc[list(vt), "t_ref_count"]
            vaf_means.setdefault(vt, np.mean(alt / (alt + ref)))
            # TODO: Is this the correct calculation for the cluster VAF?
            # TODO: Are the clusters what we want them to be?

        rank = sorted(vaf_means.items(), key=lambda kv: kv[1], reverse=True)
        input_tree = []
        for i, vt in enumerate(rank):
            mutcl[list(vt[0])] = str(i + 1)
            input_tree.append(str(i + 1) + "-" + str(i + 2))
        input_tree = ",".join(input_tree)
        clonal_struct = ClonalStructure(input_tree, 0, set(mutcl.unique()))

    # identify phased mutations per clone
    clones = dict()
    for cl in mutcl.unique():
        pos = mutcl[mutcl == cl].index
        for tier in cluster_mutations(pos, phase_p, vaf_skew_pval):
            clones.setdefault(tier, [cl, set(tier), set(tier)])
    # generate phased mutations
    for tx, ty in itertools.combinations(clones.keys(), 2):
        cx = clones[tx][0]
        cy = clones[ty][0]
        if cx == cy:
            continue
        if not clonal_struct.is_related(cx, cy):
            continue
        tier = cluster_mutations(list(tx) + list(ty), phase_p, vaf_skew_pval)
        if not tier.difference(set([tx, ty])):
            continue
        txy = set(list(tx) + list(ty))
        if clonal_struct.tree[cx].depth > clonal_struct.tree[cy].depth:
            clones[tx][1] = clones[tx][1].union(txy)
        else:
            clones[ty][1] = clones[ty][1].union(txy)
        clones[tx][2] = clones[tx][2].union(txy)
        clones[ty][2] = clones[ty][2].union(txy)

    # assign phase ID and resolve clonal structure
    df["ClonalStructure"] = [set([i]) for i in mutcl]
    df["PhaseID"] = np.nan
    blocks = set()
    for cl, tier, block in clones.values():
        if len(tier) == 1:
            continue
        df.loc[list(tier), "ClonalStructure"].apply(lambda x: x.add(cl))
        if tier == block:
            PHASE_ID += 1
            df.loc[list(block), "PhaseID"] = PHASE_ID

    is_unphased = df["PhaseID"].isna()
    df.loc[is_unphased, "MNP"] = np.nan
    df.loc[is_unphased, "ClonalStructure"] = [
        set() for i in range(sum(is_unphased))
    ]
    if sum(df["MNP"].isna()) == len(df.index):
        return df

    # chain multi-nucleotide polymorphisms
    index = []
    mnp_dict = dict()
    for mnp_id in df["MNP"].dropna().unique():
        mnp_i = df["MNP"] == mnp_id
        mutcl = set.union(*df[mnp_i]["ClonalStructure"])
        for cl in mutcl:
            pos = df.index[
                mnp_i & df["ClonalStructure"].apply(lambda x: cl in x)
            ]
            mnp_dict = chain_mnps(
                cl, pos, df, mnp_dict, g, chrom, t_reads, n_reads
            )
            index.extend(pos)
    index = list(set(index))
    df = df.drop(index)
    for key in mnp_dict:
        df = df.append(mnp_dict[key], ignore_index=True)
    if len(df.index) == 1:
        # reset phase ID
        PHASE_ID -= 1
        df["PhaseID"] = np.nan
    return df.reset_index(drop=True)


def integrate_germline(
    df: pd.DataFrame,
    chrom: str,
    t: Sample,
    n: Sample,
    genome: Union[str, Path],
    max_dist: int,
    min_base_q,
    min_map_q,
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    global GERMLINE_ID
    germline = []
    gm_size = len(germline)
    g = pysam.FastaFile(genome)
    vcf = pysam.VariantFile(n.vcf, "r")
    tbam = pysam.AlignmentFile(t.bam, "rb")
    i: int
    for i, row in df.iterrows():  # type: ignore
        GERMLINE_ID += 1
        alt = row["Tumor_Seq_Allele"]
        mut_type = row["Variant_Type"]
        pos = row["Start_position"] - (0 if mut_type == "INS" else 1)
        end = row["Start_position"] + (1 if mut_type == "INS" else len(alt) - 1)

        var_list = [v for v in vcf.fetch(chrom, pos - max_dist, end + max_dist)]
        if len(var_list) == 0:
            GERMLINE_ID -= 1
            continue

        t_reads = get_read_names(row, g, tbam, min_base_q, min_map_q)["alt"]
        for v in var_list:
            if "PASS" not in v.filter:
                continue
            allele_depth = v.samples[n.vcf_sm]["AD"]
            for j, gm_alt in enumerate(v.alts):
                # skip spanning deletion
                if gm_alt == "*":
                    continue
                if len(v.ref) == len(gm_alt):
                    gm_mut_type = MNP_TYPE.get(len(alt), "MNP")
                elif len(v.ref) > len(gm_alt):
                    gm_mut_type = "DEL"
                else:
                    gm_mut_type = "INS"

                gm = {
                    "Chromosome": v.chrom,
                    "Start_position": v.pos,
                    "Reference_Allele": v.ref,
                    "Tumor_Seq_Allele": gm_alt,
                    "Variant_Type": gm_mut_type,
                    "t_ref_count": 0,
                    "t_alt_count": 0,
                    "n_ref_count": allele_depth[0],
                    "n_alt_count": allele_depth[j + 1],
                    "genotype": "/".join(map(str, v.samples[n.vcf_sm]["GT"])),
                    "GermlineID": GERMLINE_ID,
                }
                n_reads = get_read_names(gm, g, tbam, min_base_q, min_map_q)[
                    "alt"
                ]
                if set.intersection(*[t_reads, n_reads]):
                    if gm not in germline:
                        germline.append(gm)
        if len(germline) > gm_size:
            df.loc[i, "GermlineID"] = GERMLINE_ID
            gm_size = len(germline)
        else:
            GERMLINE_ID -= 1
    return df, pd.DataFrame(germline)


def run(args):
    check_input_args(args)

    logging.basicConfig(
        level=logging.INFO,
        format="phase | %(asctime)s [%(levelname)s] %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )

    t = Sample(args.tumor_bam, args.tumor_vcf)
    n = Sample(args.normal_bam, args.normal_vcf)

    logging.info("loading mutations from input MAF")
    m = Mutation(
        args.input_maf,
        args.tumor_name,
        t.bam_sm,
        args.phylogic_ccfs,
        args.phylogic_tree,
        args.phylogic_index,
    )

    if "n_alt_count" in m.muts and n.bam is None:
        logging.warn(
            "normal BAM not supplied; normal allele depths for MNP will be zero"
        )

    smuts = []
    gmuts = []
    m.muts["maf_idx"] = [set([i]) for i in m.muts.index]

    # use the original chromosome order (i.e. sort = False)
    chrom: str
    for chrom, df in m.muts.groupby("Chromosome", sort=False):  # type: ignore
        logging.info("processing chromosome {}".format(chrom))
        df = group_mutations(
            df, args.max_phase_radius, "PhaseID", args.min_coverage
        )
        df = group_mutations(
            df, args.max_mnp_size, "MNP", args.min_coverage, excl_indel=True
        )
        df["ClonalStructure"] = [set() for i in df.index]

        dfs = [df[df["PhaseID"].isna()]]
        if len(dfs[0].index) < len(df.index):
            for i, x in df.groupby("PhaseID", dropna=True, sort=False):
                dfs.append(
                    phase_mutations(
                        x,
                        chrom,
                        t,
                        n,
                        args.ref_fasta,
                        m.clonal_struct,
                        args.vaf_skew_pval,
                        args.min_phased_altc,
                        args.min_base_quality,
                        args.min_map_quality,
                    )
                )
        dfs = pd.concat(dfs)
        dfs.sort_values("Start_position", ignore_index=True, inplace=True)

        dfs["ClonalStructure"] = dfs["ClonalStructure"].apply(
            lambda x: ",".join(sorted(x)) if len(x) else np.nan
        )
        dfs.drop(columns="MNP", inplace=True)

        dfs["GermlineID"] = np.nan
        if n.vcf:
            dfs, dfg = integrate_germline(
                dfs,
                chrom,
                t,
                n,
                args.ref_fasta,
                args.max_phase_radius,
                args.min_base_quality,
                args.min_map_quality,
            )

            if len(dfs) > 0:
                dfs.drop_duplicates(
                    subset=[
                        "Chromosome",
                        "Start_position",
                        "Reference_Allele",
                        "Tumor_Seq_Allele",
                    ],
                    inplace=True,
                )
                dfs.sort_values(
                    "Start_position", ignore_index=True, inplace=True
                )
            if len(dfg) > 0:
                dfg.drop_duplicates(
                    subset=[
                        "Chromosome",
                        "Start_position",
                        "Reference_Allele",
                        "Tumor_Seq_Allele",
                    ],
                    inplace=True,
                )
                dfg.sort_values(
                    "Start_position", ignore_index=True, inplace=True
                )
            gmuts.append(dfg)

        smuts.append(dfs)

    # extract phased mutations
    smuts = pd.concat(smuts).reset_index(drop=True)
    pi = ~smuts["ClonalStructure"].isna() | ~smuts["GermlineID"].isna()

    # write unphased mutations as a subset of the original MAF
    logging.info("writing {}.unphased.maf".format(args.tumor_name))
    write_subset_maf(
        args.tumor_name + ".unphased.maf",
        [x.pop() for x in smuts[~pi]["maf_idx"]],
        m.maf,
    )

    # take minimum CCF for chained MNPs
    cols_to_keep = ["ClonalStructure", "PhaseID", "GermlineID"]
    is_mnp = pi & smuts["Variant_Type"].isin(["DNP", "TNP", "MNP"])
    if "ccf_hat" in smuts and sum(is_mnp) > 0:
        smuts.loc[is_mnp, "ccf_hat"] = smuts[is_mnp]["maf_idx"].apply(
            lambda x: m.muts.loc[list(x), "ccf_hat"].min()
        )
        cols_to_keep.insert(0, "ccf_hat")

    # write somatic phased mutations
    if sum(pi) > 0:
        logging.info("found {} somatic phased mutations".format(sum(pi)))
        smuts = smuts[pi].reset_index(drop=True)
        count = list(smuts.columns[smuts.columns.str.endswith("_count")])
        smuts[count] = smuts[count].fillna(value=0).astype(int)
        logging.info("writing {}.phased.vcf".format(args.tumor_name))
        if t.vcf:
            write_phase_vcf(
                args.tumor_name + ".phased.vcf", smuts, t, cols_to_keep
            )
        else:
            write_phase_vcf_from_scratch(
                args.tumor_name + ".phased.vcf", t.bam_sm, n.bam_sm, smuts
            )
        logging.info("writing {}.phased.maflite.tsv".format(args.tumor_name))
        write_phase_maflite(
            args.tumor_name + ".phased.maflite.tsv", smuts, cols_to_keep
        )

        # write germline variants phased with somatic variants
        if gmuts:
            gmuts = pd.concat(gmuts).reset_index(drop=True)
            logging.info(
                "found {} phased germline mutations".format(len(gmuts.index))
            )
            count = list(gmuts.columns[gmuts.columns.str.endswith("_count")])
            gmuts[count] = gmuts[count].astype(int)
            logging.info("writing {}.phased.vcf".format(args.normal_name))
            write_phase_vcf(
                args.normal_name + ".phased.vcf", gmuts, n, ["GermlineID"]
            )
            logging.info(
                "writing {}.phased.maflite.tsv".format(args.normal_name)
            )
            write_phase_maflite(
                args.normal_name + ".phased.maflite.tsv",
                gmuts,
                ["genotype", "GermlineID"],
            )
    logging.info("all done!")
