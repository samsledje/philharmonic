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




def identify_allowed_pfam_domains():
    df_go = pd.read_csv(PROJ_DIR + "/data/processed/selected_GO_terms.csv").rename(columns={"GOID":"GO"})

    godomain_miner_files = glob.glob(PROJ_DIR + "/data/raw/go/pfam_go??_most_specific.txt")
    df_gopfam_map = pd.concat([ pd.read_csv(f, delimiter=";") for f in godomain_miner_files]).reset_index(drop=True)
    dbg_print("Flag 754.01 ", df_go.shape, df_gopfam_map.shape)

    dfx = pd.merge(df_go[["GO","namespace","name"]], df_gopfam_map[["GO","PFAM"]])
    dfx = dfx.sort_values("PFAM").drop_duplicates().reset_index()
    
    dbg_print("Flag 754.10 ", dfx.shape) #, dfx.head().values)
    
    allowed_pfam_domains = set(dfx["PFAM"].values)

    pfam2go_map = defaultdict(list)
    for i in range(dfx.shape[0]):
        pfam2go_map[ dfx["PFAM"].iat[i]].append("%s|%s|%s" % ( dfx["GO"].iat[i], dfx["namespace"].iat[i], dfx["name"].iat[i]))

    fh = open("/tmp/t.domainminer_allowed_pfam","w")
    for k in allowed_pfam_domains: print(k,file=fh)
    fh2 = open("/tmp/t.domainminer_pfam2go_map","w")
    for k,v in pfam2go_map.items(): print(k,v,file=fh2)
    
    return allowed_pfam_domains, pfam2go_map




def collate_hmmscan_results(species):
    domtblout_files = glob.glob(PROJ_DIR + "/data/processed/{0}_hmmscan_chunk-*-*.domtblout".format(species))
    dbg_print("Flag 563.01 ", domtblout_files)
            
    pfam2info = defaultdict(list)
    for f in domtblout_files:
        for line in open(f,'r'):
            try:
                _ = int(line[70:76])
            except:
                continue
            #if (species+"_") not in line[30:65]: continue
            pfamid = line[21:32].strip()
            pdamid = line[38:59].strip() #using "pdam" as a proxy name for generic species
            pdamlen = line[70:76].strip()
            #dbg_print("Flag 563.20 ", (pfamid, pdamid, pdamlen, line[:80]))
            try:
                if 50 < int(pdamlen) <800:
                    pfam2info[pfamid.split(".")[0]].append((pfamid, pdamid))
            except:
                print (line)
                raise
            
    fh = open("/tmp/t.collate_hmmscan","w")
    for k,v in pfam2info.items(): print(k,v,file=fh)
    
    return pfam2info




def select_species_proteins(allowed_pfam_domains, pfam2go_map, pfam2hmmscan, manual_annot_proteins):
    dbg_print("Flag 842.10 ", len(allowed_pfam_domains), len(pfam2go_map), len(pfam2hmmscan), len(manual_annot_proteins))
    #using "pdam" as a proxy name for generic species
    
    d = {}; n=m=0
    for pfamid_base, L in pfam2hmmscan.items():
        n += 1
        if pfamid_base not in allowed_pfam_domains: continue
        m += 1
        for pfamid, pdamid in L:
            if pdamid not in d:
                special = ''
                for c, pdamid_list in manual_annot_proteins.items():
                    if pdamid in pdamid_list:
                        special = c

                d[pdamid] = (special, set([pfamid]), set(pfam2go_map[pfamid_base]))
            else:
                d[pdamid][1].add(pfamid)
                d[pdamid][2].update(pfam2go_map[pfamid_base])
        
    for c, pdamid_list in manual_annot_proteins.items():
        for pdamid in pdamid_list:
            if pdamid in d: continue
            d[pdamid] = (c, set(), set())
        
    dbg_print("Flag 842.50 ", len(d), m, n)
    
    l = []
    for k,v in d.items(): 
        l.append( (k, v[0], ";".join(v[1]), ";".join(v[2])))

    hdr = "species_id,manual_annot,pfam_list,GO_list".split(",")

    return hdr, l




def write_candidate_pair_list(df, nPointPairs, outfile2):
    N = df.shape[0]
    l0 = list(range(N))

    dbg_print("Flag 569.10 ", df.shape, nPointPairs, outfile2)
    j_u = np.random.choice(l0, 10*int(nPointPairs), True, (df["wt"]/df["wt"].sum()).values)
    j_v = np.random.choice(l0, 10*int(nPointPairs), True, (df["wt"]/df["wt"].sum()).values)
    valid = j_u < j_v  #get rid of potential duplicates (x1,x2) and (x2,x1) as well as (x1,x1)
    i_u = j_u[valid]
    i_v = j_v[valid]

    V0 = df["protein"].values
    df2 = pd.DataFrame( zip(V0[i_u], V0[i_v]), columns=["protein_A","protein_B"])
    df2 = df2.drop_duplicates().reset_index(drop=True)
    df2 = df2.head(nPointPairs)


    dbg_print("Flag 569.30 ", len(j_u), sum(valid), len(i_u), df2.shape)
    df2.to_csv(outfile2, index=False)
    return df2



def write_protein_list(outfile1, species):
    dbg_print("Flag 469.01 ")

    allowed_pfam_domains, pfam2go_map = identify_allowed_pfam_domains()
    dbg_print("Flag 469.10 ", len(allowed_pfam_domains), len(pfam2go_map))
    
    pfam2hmmscan = collate_hmmscan_results(species)
    dbg_print("Flag 469.20 ", len(pfam2hmmscan))

    manual_annot_proteins = {}
    if species=="pdam":
        receptor_proteins = [l.strip() for l in open(PROJ_DIR + "/data/raw/Judith_receptor_proteins_20200612.txt",'r') if l.strip()!=""]
        woundhealing_proteins =  [l.strip() for l in open(PROJ_DIR + "/data/raw/Judith_woundhealing_proteins_20200612.txt",'r') if l.strip()!=""]
        manual_annot_proteins["Receptor"] = receptor_proteins
        manual_annot_proteins["Wound-heal"] = woundhealing_proteins
        
        dbg_print("Flag 469.30 ", len(receptor_proteins), len(woundhealing_proteins))
        
    elif species in ["cowV2"]:
        rumen_proteins = [l.strip() for l in open("/afs/csail.mit.edu/u/s/samsl/Work/DSCRIPT_Validation_RECOMB/Cow_CellSystems/RumenSpecificProteins.txt",'r') if l.strip()!=""]
        manual_annot_proteins["Rumen"] = rumen_proteins

        
        dbg_print("Flag 469.31 ", len(rumen_proteins))

    hdr, l = select_species_proteins(allowed_pfam_domains, pfam2go_map, pfam2hmmscan, manual_annot_proteins)

    dbg_print ("Flag 469.50 ", hdr, len(l))
    csvout = csv.writer(open(outfile1,'w'))
    csvout.writerow(hdr)
    csvout.writerows(l)
    
    return hdr, l


    
def run(args, extra_args):
    
    dbg_print("Flag 401.01 ")

    outfile1 = "{0}/proteins_list_{1}_{2}.csv".format(args.outdir, args.species, args.outsfx)
    outfile2 = "{0}/candidate_pairs_{1}_{2}.csv".format(args.outdir, args.species, args.outsfx)
    
    hdrs, l = write_protein_list(outfile1, args.species)
    
    l1 = []
    for p, ma, _, _ in l:
        l1.append((p, 1 if not ma else args.manual_annot_wt))
        
    df = pd.DataFrame(l1, columns=["protein","wt"]).sort_values("protein").reset_index(drop=True)

    write_candidate_pair_list(df, int(args.paircount), outfile2)
    


    
##########################################################################

if __name__ == "__main__":
    sys.path.append(os.path.join(sys.path[0],PROJ_DIR+'/src'))


    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("--species", help="species to analyze. Data should be in PROJ_DIR/data/proessed . Current choices are pdam , hvul, tpseud, symbc1, or cow ", type=str, default='pdam') #choices=['pdam','hvul','tpseud','symbc1', 'cow'], default='pdam')
    parser.add_argument("--paircount", help="number of protein pairs to sample", type=int, default=10000000)
    parser.add_argument("--outdir", help="output directory (can set to '.')", type=str, default=PROJ_DIR+"/data/processed/")
    parser.add_argument("--outsfx", help="suffix to use when producing output files")
    parser.add_argument("--manual_annot_wt", help="relative weight of proteins with manual annotation", type=float, default=100)
    parser.add_argument("--extra", help="put this as the LAST option and arbitrary space-separated key=val pairs after that", type=str, nargs='*')


    args = parser.parse_args()
    extra_args = dict([a.split("=") for a in args.extra]) if args.extra else {}


    run(args, extra_args)

    
    
