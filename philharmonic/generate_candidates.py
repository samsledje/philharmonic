# python src/philharmonic_candidates.py -o {output.candidates} --go_list {input.go_shortlist} --protein_list {input.protein_shortlist} {input.sequences}

import typer
import scipy
import pandas as pd
import numpy as np
from Bio import SeqIO
from pathlib import Path
from itertools import combinations

from .utils import filter_proteins_GO

app = typer.Typer()

@app.command()
def main(
    sequences: Path = typer.Option(..., "--sequences", help="Sequences file"),
    output: Path = typer.Option(..., "-o", "--output", help="Output file"),
    seq_out: Path = typer.Option(
        None,
        "--seq_out",
        help="Optional output: report which sequences were kept after GO filtering",
    ),
    paircount: int = typer.Option(
        -1, "--paircount", help="Number of protein pairs to sample"
    ),
    go_filter: str = typer.Option(
        None, "--go_filter", help="GO terms by which to filter"
    ),
    go_map: str = typer.Option(
        ..., "--go_map", help="GO map file. Required if --go_filter is set"
    ),
    go_database: str = typer.Option(
        ..., "--go_database", help="GO database file. Required if --go_filter is set"
    ),
    seed: int = typer.Option(42, "--seed", help="Random seed"),
):
    """Generate candidate protein pairs from a sequences file"""
    sequences = SeqIO.to_dict(SeqIO.parse(sequences, "fasta"))
    protein_names = list(sequences.keys())

    # Create list of filtered GO terms
    allowed_proteins = filter_proteins_GO(protein_names, go_map, go_database, go_filter)

    if seq_out is not None:
        with open(seq_out, "w") as f:
            for p in allowed_proteins:
                f.write(f">{p}\n{sequences[p].seq}\n")

    # Create DataFrame
    g = np.random.Generator(np.random.PCG64(seed))
    protein_pairs = list(combinations(allowed_proteins, 2))

    if paircount < 0:
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
