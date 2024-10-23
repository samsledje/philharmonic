# python src/philharmonic_candidates.py -o {output.candidates} --go_list {input.go_shortlist} --protein_list {input.protein_shortlist} {input.sequences}

from itertools import combinations
from pathlib import Path

import numpy as np
import pandas as pd
import scipy
import typer
from Bio import SeqIO
from loguru import logger

from .utils import filter_proteins_GO

app = typer.Typer()


@app.command()
def main(
    sequence_file_path: Path = typer.Option(..., "--sequences", help="Sequences file"),
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
        ..., "--go_filter", help="GO terms by which to filter"
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
    sequences = SeqIO.to_dict(SeqIO.parse(sequence_file_path, "fasta"))
    protein_names = list(sequences.keys())

    # Create list of filtered GO terms
    allowed_proteins = filter_proteins_GO(protein_names, go_filter, go_map, go_database)
    logger.info(
        f"Generating candidates from {len(allowed_proteins)} proteins after GO filtering"
    )

    if seq_out is not None:
        with open(seq_out, "w") as f:
            for p in allowed_proteins:
                f.write(f">{p}\n{sequences[p].seq}\n")

    # Create DataFrame
    rng = np.random.Generator(np.random.PCG64(seed))
    all_pairs = list(combinations(allowed_proteins, 2))

    if paircount > scipy.special.comb(len(allowed_proteins), 2):
        raise ValueError("Pair count exceeds number of possible pairs")
    elif paircount < 0:
        logger.info(f"Using all {len(all_pairs)} possible pairs")
        sampled_pairs = np.asarray(all_pairs)
    else:
        logger.info(f"Sampling {paircount} pairs")
        sampled_pairs = rng.choice(all_pairs, paircount, replace=False)

    logger.info(f"Writing {len(sampled_pairs)} pairs to {output}")
    candidates = pd.DataFrame(sampled_pairs, columns=["protein_A", "protein_B"])
    candidates = candidates.drop_duplicates().reset_index(drop=True)

    candidates.to_csv(str(output), sep="\t", index=False, header=False)


if __name__ == "__main__":
    app()
