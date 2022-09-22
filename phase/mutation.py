import os
from collections import OrderedDict

import numpy as np
import pandas as pd


required_maf_cols = OrderedDict({ 'Chromosome'      : str,
                                  'Start_position'  : int,
                                  'Reference_Allele': str,
                                  'Tumor_Seq_Allele': str,
                                  'Variant_Type'    : str,
                                  't_alt_count'     : float,
                                  't_ref_count'     : float })

optional_maf_cols = OrderedDict({ 'n_alt_count': float,
                                  'n_ref_count': float,
                                  'ccf_hat'    : float })

phylogic_ccf_cols = OrderedDict({ 'Chromosome'        : str,
                                  'Start_position'    : int,
                                  'Reference_Allele'  : str,
                                  'Tumor_Seq_Allele'  : str,
                                  'Variant_Type'      : str,
                                  'Cluster_Assignment': str })

class Tree:
    def __init__(self, name, parent=None, children=[]):
        self.name = name
        self.parent = parent
        self.children = children
        self.depth = -1
    def __repr__(self):
        fmt = 'node={}@[{}]: parent={}, children=[{}]'
        return fmt.format(self.name, self.depth, self.parent,
                          ', '.join(map(str, self.children)))


CLONAL_CLUSTER = '1'
class ClonalStructure:
    def __init__(self, input_tree, index, cluster):
        self.tree = None
        self.excluded_clone = set()
        
        tr = input_tree
        if os.path.exists(input_tree):
            tr = [row.strip().split('\t')[-1] for row in open(input_tree, 'r')]
            if len(tr) < index or index < 1:
                raise ValueError('{}-th tree does not exist'.format(index))
            tr = tr[index]
        tree = dict()
        for pair in tr.split(','):
            if pair.startswith('None'):
                continue
            parent, child = pair.split('-')
            tree.setdefault(parent, Tree(parent, None, [])).children.append(child)
            tree.setdefault(child, Tree(child, parent, []))
            
        self.tree = self.set_depths(tree)
        self.excluded_clone = cluster.difference(set(self.tree.keys()))
        
    def set_depths(self, tree):
        depth = 0
        tree[CLONAL_CLUSTER].depth = depth
        c = tree[CLONAL_CLUSTER].children
        while c != []:
            depth += 1
            next_c = []
            for x in c:
                tree[x].depth = depth
                next_c.extend(tree[x].children)
            c = next_c
        return tree
    
    def get_offsprings(self, n):
        children = []
        c = self.tree[n].children
        while c != []:
            children.extend(c)
            next_c = []
            for x in c:
                next_c.extend(self.tree[x].children)
            c = next_c
        return children
    
    def is_related(self, x, y):
        if x == y:
            return True
        elif y in self.get_offsprings(x):
            return True
        elif x in self.get_offsprings(y):
            return True
        else:
            return False

        
class Mutation:
    def __init__(self, maf, name, smid, phylogic_ccf, phylogic_tree, tree_index):
        self.maf = maf
        self.muts = pd.read_csv(maf, sep='\t', comment='#', dtype=str)
        self.tid = set(self.muts['Tumor_Sample_Barcode'])
        self.num = len(self.muts.index)
        self.clonal_struct = None

        if name not in self.tid:
            raise Exception(name + " not in " + os.path.basename(maf))
        
        if len(self.tid) > 1:
            row = self.muts['Tumor_Sample_Barcode'] == name
            self.muts[row].reset_index(drop=True, inplace=True)
            self.num = sum(row)
        self.tid = name
        
        # rename MAF header
        self.muts.rename(columns={'Tumor_Seq_Allele2':'Tumor_Seq_Allele'}, inplace=True)
        if 'Start_Position' in self.muts:
            self.muts.rename(columns={'Start_Position':'Start_position'}, inplace=True)

        # gather columns from MAF
        maf_cols = required_maf_cols.copy()
        for oc in optional_maf_cols:
            if oc in self.muts.columns:
                maf_cols[oc] = optional_maf_cols[oc]
        try:
            self.muts = self.muts[maf_cols.keys()].astype(maf_cols)
        except Exception as e:
            raise KeyError(str(e).replace('index', 'MAF'))

        # set Phylogic CCF and tree information
        if phylogic_ccf and phylogic_tree:
            df = pd.read_csv(phylogic_ccf, sep='\t', comment='#', dtype=str)
            t1 = df['Sample_ID'] == self.tid
            t2 = df['Sample_ID'] == smid
            if sum(t1) == 0 and sum(t2) == 0:
                raise Exception('cannot find {} in Phylogic CCF'.format(
                    ' or '.join({self.tid, smid})))
            elif abs(sum(t1)-self.num) > abs(sum(t2)-self.num):
                df = df[t2]
            else:
                df = df[t1]
            df.reset_index(drop=True, inplace=True)
            df = df[phylogic_ccf_cols.keys()].astype(phylogic_ccf_cols)
            self.muts = pd.merge(self.muts, df, how='left',
                                 on=list(required_maf_cols.keys())[:5])
            self.clonal_struct = ClonalStructure(phylogic_tree, tree_index,
                                                 set(self.muts['Cluster_Assignment']))
