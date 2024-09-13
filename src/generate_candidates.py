# python src/philharmonic_candidates.py -o {output.candidates} --go_list {input.go_shortlist} --protein_list {input.protein_shortlist} {input.sequences}

import argparse
import random
import itertools
import scipy
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

    go_subset = []
    for c in visited_go_terms:
        go_subset.append((c, go2desc[c][0], go2desc[c][1], ";".join(go2children[c])))

    return sorted(go_subset)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Generate a list of candidate PPIs")
    parser.add_argument("--sequences", required=True, type=str, help="Sequences file")
    parser.add_argument("-o", "--output", required=True, type=str, help="Output file")
    parser.add_argument("--seq_out", type=str, help="Optional output: report which sequences were kept after GO filtering")
    parser.add_argument(
        "--paircount",
        type=int,
        help="Number of protein pairs to sample",
        default=-1,
    )
    parser.add_argument("--go_filter", type=str, help="GO terms by which to filter")
    parser.add_argument("--go_map", type=str, help="GO map file. Required if --go_filter is set")
    parser.add_argument("--go_database", type=str, help="GO database file. Required if --go_filter is set")
    parser.add_argument("--seed", type=int, help="Random seed", default=42)

    args = parser.parse_args()

    sequences = SeqIO.to_dict(SeqIO.parse(args.sequences, "fasta"))
    protein_names = list(sequences.keys())

    # Create list of filtered GO terms
    allowed_go = []
    allowed_proteins = protein_names
    if args.go_filter is not None:
        with open(args.go_filter, "r") as f:
            allowed_go = [l.strip() for l in f]

        # Get children of allowed GO terms
        allowed_go = [i[0] for i in subset_go_graph(args.go_database, allowed_go)]
        allowed_go = set(allowed_go)

        # Filter proteins by GO terms
        go_map = parse_GO_map(args.go_map)
        allowed_proteins = []
        for protein, go_terms in go_map.items():
            if go_terms is not None and any(gt in allowed_go for gt in go_terms):
                allowed_proteins.append(protein)
        allowed_proteins = set(allowed_proteins)

    if args.seq_out is not None:
        with open(args.seq_out, "w") as f:
            for p in allowed_proteins:
                f.write(f">{p}\n{sequences[p].seq}\n")
    
    # Create DataFrame
    g = np.random.Generator(np.random.PCG64(args.seed))
    protein_pairs = list(combinations(allowed_proteins, 2))

    if args.paircount < 0:
        sampled_pairs = protein_pairs
    elif args.paircount > scipy.special.comb(len(allowed_proteins), 2):
        raise ValueError("Pair count exceeds number of possible pairs")
    else:
        sampled_pairs = g.choice(protein_pairs, args.paircount, replace=False)
    
    candidates = pd.DataFrame(sampled_pairs, columns=["protein_A", "protein_B"])
    candidates = candidates.drop_duplicates().reset_index(drop=True)

    candidates.to_csv(args.output, index=False, header=None, sep="\t")