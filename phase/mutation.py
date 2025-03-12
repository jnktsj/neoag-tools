import os
from collections import OrderedDict
from pathlib import Path
from typing import Dict, List, Optional, Set, Union

import pandas as pd


required_maf_cols = OrderedDict(
    {
        "Chromosome": str,
        "Start_position": int,
        "Reference_Allele": str,
        "Tumor_Seq_Allele": str,
        "Variant_Type": str,
        "t_alt_count": float,
        "t_ref_count": float,
        "ref_context": str,
    }
)

optional_maf_cols = OrderedDict(
    {"n_alt_count": float, "n_ref_count": float, "ccf_hat": float}
)

phylogic_ccf_cols = OrderedDict(
    {
        "Chromosome": str,
        "Start_position": int,
        "Reference_Allele": str,
        "Tumor_Seq_Allele": str,
        "Variant_Type": str,
        "Cluster_Assignment": str,
    }
)


class Tree:
    """
    A node in a phylogenetic tree.
    """

    name: str
    """The name of the node."""
    parent: "Optional[str]"
    """The parent node of the node."""
    children: "List[str]"
    """All child nodes of the node."""
    depth: int
    """The depth of the node in the tree - the distance from the root."""

    def __init__(
        self,
        name: str,
        parent: "Optional[str]" = None,
        children: "Optional[List[str]]" = None,
    ):
        """Initializes a new node in the phylogenetic itree.

        Args:
            name (str): The name of the node.
            parent (Tree, optional): The name of the parent of the node. Defaults to None.
            children (List[Tree], optional): The child nodes of the node. Defaults to None.
        """
        self.name = name
        self.parent = parent
        self.children = children if children is not None else []
        self.depth = -1

    def __repr__(self):
        return f"node={self.name}@[{self.depth}]: parent={self.parent}, children=[{', '.join(str(child) for child in self.children)}]"


CLONAL_CLUSTER = "1"
"""The name of the clonal cluster."""


def set_depths(tree: Dict[str, Tree], clonal_cluster=CLONAL_CLUSTER):
    """Sets the depths of all nodes in the tree.

    Args:
        tree (Dict[str, Tree]): A dictionary from node names to the node data.
    """
    depth = 0
    tree[clonal_cluster].depth = depth
    c = tree[clonal_cluster].children
    while c != []:
        depth += 1
        next_c = []
        for x in c:
            tree[x].depth = depth
            next_c.extend(tree[x].children)
        c = next_c


class ClonalStructure:
    tree: Dict[str, Tree]
    excluded_clone: Set[str]

    def __init__(self, input_tree: str, index: int, cluster: Set[str]):
        """Initializes the clonal structure.

        Args:
            input_tree (str): Either a tree in the format "parent1-child1,parent2-child2,...", or a path to a file containing trees in this format.
            index (int): The index of the tree in the file.
            cluster (_type_): _description_
        """
        self.excluded_clone = set()

        tr = input_tree
        if os.path.exists(input_tree):
            # Getting the tree from a file.
            tr = [row.strip().split("\t")[-1] for row in open(input_tree, "r")]
            if len(tr) < index or index < 1:
                raise ValueError("{}-th tree does not exist".format(index))
            tr = tr[index]

        self.tree = dict()
        for pair in tr.split(","):
            if pair.startswith("None"):
                continue
            parent, child = pair.split("-")
            self.tree.setdefault(
                parent, Tree(parent, None, [])
            ).children.append(child)
            self.tree.setdefault(child, Tree(child, parent, []))

        set_depths(self.tree)
        self.excluded_clone = cluster.difference(set(self.tree.keys()))

    def get_offsprings(self, n: str) -> List[str]:
        """Gets all nodes in the subtree rooted by the given node, excluding the node."""
        children = []
        c = self.tree[n].children
        while c != []:
            children.extend(c)
            next_c = []
            for x in c:
                next_c.extend(self.tree[x].children)
            c = next_c
        return children

    def is_related(self, x: str, y: str) -> bool:
        """Checks if one of the nodes is the ancestor of the other.

        Args:
            x (str): The name of a node in the phylogenetic tree.
            y (str): The name of a node in the phylogenetic tree.

        Returns:
            bool: True iff x is the ancestor of y or y is the ancestor of x.
        """
        return (
            (x == y)
            or (x in self.get_offsprings(y))
            or (y in self.get_offsprings(x))
        )


class Mutation:
    maf: str | Path
    """The path to the mutation MAF file."""
    muts: pd.DataFrame
    """The dataframe holding the mutations."""
    tid: str
    """The tumor ID."""
    num: int
    """The number of mutations stored by the object."""
    clonal_struct: Optional[ClonalStructure]
    """The subclonal structure of the tumor sample."""

    def __init__(
        self,
        maf: Union[str, Path],
        name: str,
        smid: str,
        phylogic_ccf: Union[str, Path, None],
        phylogic_tree: Optional[str],
        tree_index: int,
    ):
        """A class holding the mutation data.

        Args:
            maf (Union[str, Path]): The path to the mutation MAF file.
            name (str): The name of the sample (?).
            smid (str): The ID of the sample.
            phylogic_ccf (Union[str, Path, None]): A path to a file mapping every mutation to a single CCF.
            phylogic_tree (Optional[str]): A string describing a tree or a path to a tree file.
            tree_index (int): The index of the tree in the tree file.
        """
        self.maf = maf
        self.muts = pd.read_csv(maf, sep="\t", comment="#", dtype=str)
        tid = set(self.muts["Tumor_Sample_Barcode"])
        self.num = len(self.muts.index)
        self.clonal_struct = None

        if name not in tid:
            raise Exception(name + " not in " + os.path.basename(maf))

        if len(tid) > 1:
            row = self.muts["Tumor_Sample_Barcode"] == name
            self.muts[row].reset_index(drop=True, inplace=True)
            self.num = sum(row)
        self.tid = name

        # rename MAF header
        self.muts.rename(
            columns={"Tumor_Seq_Allele2": "Tumor_Seq_Allele"}, inplace=True
        )
        if "Start_Position" in self.muts:
            self.muts.rename(
                columns={"Start_Position": "Start_position"}, inplace=True
            )

        # gather columns from MAF
        maf_cols = required_maf_cols.copy()
        for oc in optional_maf_cols:
            if oc in self.muts.columns:
                maf_cols[oc] = optional_maf_cols[oc]
        try:
            self.muts = self.muts[maf_cols.keys()].astype(maf_cols)
        except Exception as e:
            raise KeyError(str(e).replace("index", "MAF"))

        # set Phylogic CCF and tree information
        if phylogic_ccf and phylogic_tree:
            df = pd.read_csv(phylogic_ccf, sep="\t", comment="#", dtype=str)
            t1 = df["Sample_ID"] == self.tid
            t2 = df["Sample_ID"] == smid
            if sum(t1) == 0 and sum(t2) == 0:
                raise Exception(
                    "cannot find {} in Phylogic CCF".format(
                        " or ".join({self.tid, smid})
                    )
                )
            elif abs(sum(t1) - self.num) > abs(sum(t2) - self.num):
                df = df[t2]
            else:
                df = df[t1]
            df.reset_index(drop=True, inplace=True)
            df = df[phylogic_ccf_cols.keys()].astype(phylogic_ccf_cols)
            self.muts = pd.merge(
                self.muts, df, how="left", on=list(required_maf_cols.keys())[:5]
            )
            self.clonal_struct = ClonalStructure(
                phylogic_tree, tree_index, set(self.muts["Cluster_Assignment"])
            )
