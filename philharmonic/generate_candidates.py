# python src/philharmonic_candidates.py -o {output.candidates} --go_list {input.go_shortlist} --protein_list {input.protein_shortlist} {input.sequences}

import typer
import scipy
import pandas as pd
import numpy as np
from Bio import SeqIO
from collections import defaultdict
from pathlib import Path
from itertools import combinations
from loguru import logger

from .utils import parse_GO_map

app = typer.Typer()

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

@app.command()
def main(
    sequences: Path = typer.Option(..., "--sequences", help="Sequences file"),
    output: Path = typer.Option(..., "-o", "--output", help="Output file"),
    seq_out: Path = typer.Option(None, "--seq_out", help="Optional output: report which sequences were kept after GO filtering"),
    paircount: int = typer.Option(-1, "--paircount", help="Number of protein pairs to sample"),
    go_filter: str = typer.Option(None, "--go_filter", help="GO terms by which to filter"),
    go_map: str = typer.Option(..., "--go_map", help="GO map file. Required if --go_filter is set"),
    go_database: str = typer.Option(..., "--go_database", help="GO database file. Required if --go_filter is set"),
    seed: int = typer.Option(42, "--seed", help="Random seed"),
):
    """Generate candidate protein pairs from a sequences file"""
    sequences = SeqIO.to_dict(SeqIO.parse(sequences, "fasta"))
    protein_names = list(sequences.keys())

    # Create list of filtered GO terms
    allowed_go = []
    allowed_proteins = protein_names
    if go_filter is not None:
        with open(go_filter, "r") as f:
            allowed_go = [l.strip() for l in f]

        # Get children of allowed GO terms
        allowed_go = [i[0] for i in subset_go_graph(go_database, allowed_go)]
        allowed_go = set(allowed_go)

        # Filter proteins by GO terms
        go_map = parse_GO_map(go_map)
        allowed_proteins = []
        for protein, go_terms in go_map.items():
            if go_terms is not None and any(gt in allowed_go for gt in go_terms):
                allowed_proteins.append(protein)
        allowed_proteins = set(allowed_proteins)

    if seq_out is not None:
        with open(seq_out, "w") as f:
            for p in allowed_proteins:
                f.write(f">{p}\n{sequences[p].seq}\n")
    
    # Create DataFrame
    g = np.random.Generator(np.random.PCG64(seed))
    protein_pairs = list(combinations(allowed_proteins, 2))

    if args.paircount < 0:
        sampled_pairs = protein_pairs
    elif paircount > scipy.special.comb(len(allowed_proteins), 2):
        raise ValueError("Pair count exceeds number of possible pairs")
    else:
        sampled_pairs = g.choice(protein_pairs, paircount, replace=False)
    
    candidates = pd.DataFrame(sampled_pairs, columns=["protein_A", "protein_B"])
    candidates = candidates.drop_duplicates().reset_index(drop=True)

    candidates.to_csv(output, index=False, header=None, sep="\t")



if __name__ == "__main__":
    app()

