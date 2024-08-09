# python src/philharmonic_candidates.py -o {output.candidates} --go_list {input.go_shortlist} --protein_list {input.protein_shortlist} {input.sequences}

import argparse
import pandas as pd
import numpy as np
from Bio import SeqIO
from collections import defaultdict

from utils import parse_GO_map
from itertools import combinations


def parse_go_graph(go_graph_file):
    go2children = defaultdict(list)
    go2desc = {}

    def process_block(b):
        if not ("id" in b and "name" in b and "namespace" in b):
            return

        go2desc[b["id"]] = (b["namespace"], b["name"])
        if "parent" in b:
            for c in b["parent"]:
                go2children[c].append(b["id"])

    block = {}
    for line in open(go_graph_file):
        if line.startswith("[Term]"):
            if block:
                # dbg_print("Flag 231.20 ", block["id"])
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

    if block:
        process_block(block)

    # dbg_print("Flag 231.40 ", len(go2children), len(go2desc))
    return go2children, go2desc


def subset_go_graph(go_graph_file, go_included_terms):
    # dbg_print("Flag 341.01 ")
    go2children, go2desc = parse_go_graph(go_graph_file)
    # dbg_print("Flag 341.05 ", len(go2children), len(go2desc))

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


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Generate a list of candidate PPIs")
    parser.add_argument("--sequences", required=True, type=str, help="Sequences file")
    parser.add_argument("-o", "--output", required=True, type=str, help="Output file")
    parser.add_argument("--go_map", type=str, required=True, help="GO map file")
    parser.add_argument(
        "--paircount",
        type=int,
        help="Number of protein pairs to sample",
        default=10000000,
    )
    parser.add_argument("--go_database", type=str, help="GO database file")
    parser.add_argument("--go_filter", type=str, help="GO terms to filter by")
    parser.add_argument("--go_list", type=str, help="GO shortlist to include")
    parser.add_argument(
        "--manual_annot_wt",
        type=float,
        help="Relative weight of proteins with manual annotation",
        default=100,
    )
    parser.add_argument("--protein_list", type=str, help="Protein shortlist")

    args = parser.parse_args()

    sequences = SeqIO.to_dict(SeqIO.parse(args.sequences, "fasta"))
    protein_names = list(sequences.keys())
    candidates = list(combinations(sequences.keys(), 2))

    # Create list of filtered GO terms
    go_filter = pd.read_csv(args.go_filter, sep=",")
    allowed_go = list(go_filter[go_filter["useflag"] == 1]["GOID"].values)
    shortlisted_go = [line.strip() for line in open(args.go_list, "r")]
    allowed_go.extend(shortlisted_go)

    # Get children of allowed GO terms
    allowed_go = [i[0] for i in subset_go_graph(args.go_database, allowed_go)]
    allowed_go = set(allowed_go)

    # Filter proteins by GO terms
    go_map = parse_GO_map(args.go_map)
    allowed_proteins = []
    for protein, go_terms in go_map.items():
        if go_terms is not None and any(gt in allowed_go for gt in go_terms):
            allowed_proteins.append(protein)
    shortlisted_proteins = [line.strip() for line in open(args.protein_list, "r")]
    allowed_proteins.extend(shortlisted_proteins)
    allowed_proteins = set(allowed_proteins)

    # Set weight for sampling
    protein_weight_df = []
    for ap in allowed_proteins:
        protein_weight_df.append(
            (ap, args.manual_annot_wt if ap in shortlisted_proteins else 1)
        )
    protein_weight_df = pd.DataFrame(protein_weight_df, columns=["protein", "wt"])

    l0 = list(range(len(protein_weight_df)))
    j_u = np.random.choice(
        l0,
        int(args.paircount),
        True,
        (protein_weight_df["wt"] / protein_weight_df["wt"].sum()).values,
    )
    j_v = np.random.choice(
        l0,
        int(args.paircount),
        True,
        (protein_weight_df["wt"] / protein_weight_df["wt"].sum()).values,
    )
    valid = (
        j_u < j_v
    )  # get rid of potential duplicates (x1,x2) and (x2,x1) as well as (x1,x1)
    i_u = j_u[valid]
    i_v = j_v[valid]

    V0 = protein_weight_df["protein"].values
    candidates = pd.DataFrame(zip(V0[i_u], V0[i_v]), columns=["protein_A", "protein_B"])
    candidates = candidates.drop_duplicates().reset_index(drop=True)

    candidates.to_csv(args.output, index=False, header=None, sep="\t")

# def write_candidate_pair_list(df, nPointPairs, outfile2):
#     N = df.shape[0]
#     l0 = list(range(N))

#     dbg_print("Flag 569.10 ", df.shape, nPointPairs, outfile2)
#     j_u = np.random.choice(l0, 10*int(nPointPairs), True, (df["wt"]/df["wt"].sum()).values)
#     j_v = np.random.choice(l0, 10*int(nPointPairs), True, (df["wt"]/df["wt"].sum()).values)
#     valid = j_u < j_v  #get rid of potential duplicates (x1,x2) and (x2,x1) as well as (x1,x1)
#     i_u = j_u[valid]
#     i_v = j_v[valid]

#     V0 = df["protein"].values
#     df2 = pd.DataFrame( zip(V0[i_u], V0[i_v]), columns=["protein_A","protein_B"])
#     df2 = df2.drop_duplicates().reset_index(drop=True)
#     df2 = df2.head(nPointPairs)


#     dbg_print("Flag 569.30 ", len(j_u), sum(valid), len(i_u), df2.shape)
#     df2.to_csv(outfile2, index=False)
#     return df2


# def write_protein_list(outfile1, species):
#     dbg_print("Flag 469.01 ")

#     allowed_pfam_domains, pfam2go_map = identify_allowed_pfam_domains()
#     dbg_print("Flag 469.10 ", len(allowed_pfam_domains), len(pfam2go_map))

#     pfam2hmmscan = collate_hmmscan_results(species)
#     dbg_print("Flag 469.20 ", len(pfam2hmmscan))

#     manual_annot_proteins = {}
#     if species=="pdam":
#         receptor_proteins = [l.strip() for l in open(PROJ_DIR + "/data/raw/Judith_receptor_proteins_20200612.txt",'r') if l.strip()!=""]
#         woundhealing_proteins =  [l.strip() for l in open(PROJ_DIR + "/data/raw/Judith_woundhealing_proteins_20200612.txt",'r') if l.strip()!=""]
#         manual_annot_proteins["Receptor"] = receptor_proteins
#         manual_annot_proteins["Wound-heal"] = woundhealing_proteins

#         dbg_print("Flag 469.30 ", len(receptor_proteins), len(woundhealing_proteins))

#     elif species in ["cowV2"]:
#         rumen_proteins = [l.strip() for l in open("/afs/csail.mit.edu/u/s/samsl/Work/DSCRIPT_Validation_RECOMB/Cow_CellSystems/RumenSpecificProteins.txt",'r') if l.strip()!=""]
#         manual_annot_proteins["Rumen"] = rumen_proteins


#         dbg_print("Flag 469.31 ", len(rumen_proteins))

#     hdr, l = select_species_proteins(allowed_pfam_domains, pfam2go_map, pfam2hmmscan, manual_annot_proteins)

#     dbg_print ("Flag 469.50 ", hdr, len(l))
#     csvout = csv.writer(open(outfile1,'w'))
#     csvout.writerow(hdr)
#     csvout.writerows(l)

#     return hdr, l


# def run(args, extra_args):

#     dbg_print("Flag 401.01 ")

#     outfile1 = "{0}/proteins_list_{1}_{2}.csv".format(args.outdir, args.species, args.outsfx)
#     outfile2 = "{0}/candidate_pairs_{1}_{2}.csv".format(args.outdir, args.species, args.outsfx)

#     hdrs, l = write_protein_list(outfile1, args.species)

#     l1 = []
#     for p, ma, _, _ in l:
#         l1.append((p, 1 if not ma else args.manual_annot_wt))

#     df = pd.DataFrame(l1, columns=["protein","wt"]).sort_values("protein").reset_index(drop=True)

#     write_candidate_pair_list(df, int(args.paircount), outfile2)


# ##########################################################################

# if __name__ == "__main__":

#     parser = argparse.ArgumentParser()
#     parser.add_argument("--species", help="species to analyze. Data should be in PROJ_DIR/data/proessed . Current choices are pdam , hvul, tpseud, symbc1, or cow ", type=str, default='pdam') #choices=['pdam','hvul','tpseud','symbc1', 'cow'], default='pdam')
#     parser.add_argument("--paircount", help="number of protein pairs to sample", type=int, default=10000000)
#     parser.add_argument("--outdir", help="output directory (can set to '.')", type=str, default=PROJ_DIR+"/data/processed/")
#     parser.add_argument("--outsfx", help="suffix to use when producing output files")
#     parser.add_argument("--manual_annot_wt", help="relative weight of proteins with manual annotation", type=float, default=100)
#     parser.add_argument("--extra", help="put this as the LAST option and arbitrary space-separated key=val pairs after that", type=str, nargs='*')


#     args = parser.parse_args()
#     extra_args = dict([a.split("=") for a in args.extra]) if args.extra else {}


#     run(args, extra_args)
