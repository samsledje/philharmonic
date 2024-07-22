#!/usr/bin/env python


###################################################################
## Authors:  Rohit Singh rsingh@alum.mit.edu
## License: MIT
###################################################################

import pandas as pd
import numpy as np
import scipy, sklearn, os, sys, string, fileinput, glob, re, math, itertools, gc, csv
from functools import partial, reduce
import  copy, multiprocessing, traceback, logging, pickle, traceback
import scipy.stats, sklearn.decomposition, sklearn.preprocessing, sklearn.covariance, sklearn.datasets
from scipy.stats import describe
import os.path
import scipy.sparse
from scipy.sparse import csr_matrix, csc_matrix
from sklearn.preprocessing import normalize
from collections import defaultdict

from base_config import *


def parse_go_graph(go_graph_file):

    go2children = defaultdict(list)
    go2desc = {}
    
    def process_block(b):
        if not ('id' in b and 'name' in b and 'namespace' in b): return
        
        go2desc[b['id']] = (b['namespace'], b['name'])
        if 'parent' in b:
            for c in b['parent']:
                go2children[c].append(b['id'])        
    
    block = {}
    for line in open(go_graph_file):
        if line.startswith("[Term]"):
            if block:
                #dbg_print("Flag 231.20 ", block["id"])
                process_block(block)
                block = {}
        if line.startswith("id: "):
            block["id"] = line.split()[1]
            
        if line.startswith("is_a: "):
            if "parent" not in block:
                block["parent"] = []
            block["parent"].append(line.split()[1])
            
        if line.startswith("name:"):
            block["name"] = line[6:].strip()
            
        if line.startswith("namespace:"):
            block["namespace"] = line[11:].strip()

    if block: process_block(block)

    dbg_print("Flag 231.40 ", len(go2children), len(go2desc))
    return go2children, go2desc


        

def subset_go_graph(go_graph_file, go_included_terms):
    dbg_print("Flag 341.01 ")
    go2children, go2desc = parse_go_graph(go_graph_file)
    dbg_print("Flag 341.05 ", len(go2children), len(go2desc))
    
    visited_go_terms = set()
    
    def dfs_visit(go2children, goid):
        for c in go2children[goid]:
            dfs_visit(go2children, c)
        visited_go_terms.add(goid)

    for c1 in go_included_terms:
        dfs_visit(go2children, c1)

    l = []
    for c in visited_go_terms:
        l.append((c, go2desc[c][0], go2desc[c][1], ";".join(go2children[c])))
        
    return sorted(l)



def run():
    go_file = PROJ_DIR + "/data/raw/go/go.obo"
    dbg_print ("Flag 469.01 ")
    
    selected_go_terms = []
    go_markedup_file = PROJ_DIR + "/data/raw/go/go_level2_marked-up.csv"
    for r in csv.reader(open(go_markedup_file)):
        if r[2]!="useflag" and int(r[2])==1: 
            selected_go_terms.append(r[0].strip())

    dbg_print ("Flag 469.10 ", selected_go_terms)
    
    l = subset_go_graph(go_file, selected_go_terms)

    dbg_print ("Flag 469.30 ", len(l))
    csvout = csv.writer(sys.stdout)
    csvout.writerow("GOID,namespace,name,child_terms".split(","))
    csvout.writerows(l)

    
##########################################################################

if __name__ == "__main__":
    run()

    
    
