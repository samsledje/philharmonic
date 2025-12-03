# %% Imports
import re
from multiprocessing import cpu_count
from pathlib import Path

import pandas as pd
import typer
from biotite.sequence.io import fasta
from datasets import Dataset
from loguru import logger
from tqdm import tqdm

from philharmonic.utils import load_cluster_json


def process_attribute(row, key, preface="", suffix=""):
    id_list = re.findall(f"{key}=([^;]+)", row["attributes"])
    return {
        "id": f"{preface}{id_list[0]}{suffix}" if len(id_list) else None,
        "chromosome": row["seqid"],
        "start": row["start"],
        "end": row["end"],
        "strand": row["strand"],
    }


def main(
    json_file: Path = typer.Argument(..., help="Input JSON file with cluster data"),
    go_file: Path = typer.Argument(..., help="GO mapping CSV file"),
    fasta_file: Path = typer.Argument(..., help="Protein sequences FASTA file"),
    gff_file: Path = typer.Argument(..., help="GFF annotation file"),
    outfile: Path = typer.Argument(..., help="Output CSV file"),
    # GFF parsing options
    gff_feature_type: str = typer.Option("CDS", help="GFF feature type to extract"),
    gff_attribute_key: str = typer.Option(
        "protein_id", help="GFF attribute key for protein ID"
    ),
    gff_id_preface: str = typer.Option(" ", help="Prefix to add to protein IDs"),
    # Taxonomy options
    species: str = typer.Option(..., help="Species name"),
    genus: str = typer.Option(..., help="Genus name"),
    family: str = typer.Option(..., help="Family name"),
    order: str = typer.Option(..., help="Order name"),
    class_name: str = typer.Option(..., "--class", help="Class name"),
    phylum: str = typer.Option(..., help="Phylum name"),
    domain: str = typer.Option(..., help="Domain name"),
    # Processing options
    num_proc: int = typer.Option(
        None,
        help="Number of processes for parallel processing (default: CPU count - 4)",
    ),
):
    """
    Convert Philharmonic cluster data to per-protein functional annotations.

    Combines cluster information, GO terms, Pfam domains, sequences, and genomic locations
    into a single CSV file.
    """
    if num_proc is None:
        num_proc = max(1, cpu_count() - 4)

    logger.info(f"Loading Philharmonic data from JSON: {json_file}")
    philharmonic_data = load_cluster_json(json_file)

    logger.info(f"Loading GO mapping from CSV: {go_file}")
    go_map = pd.read_csv(go_file)

    logger.info(f"Loading sequences from FASTA: {fasta_file}")
    fasta_fh = fasta.FastaFile.read(str(fasta_file))
    sequence_dict = {
        seq_id.split()[0]: seq for seq_id, seq in fasta.get_sequences(fasta_fh).items()
    }

    logger.info(f"Loading GFF data from file: {gff_file}")
    gff_df = pd.read_csv(
        gff_file,
        sep="\t",
        comment="#",
        header=None,
    ).dropna()
    gff_df.columns = [
        "seqid",
        "source",
        "type",
        "start",
        "end",
        "score",
        "strand",
        "phase",
        "attributes",
    ]

    # Clean GFF data
    gff_df = gff_df[gff_df["type"].isin([gff_feature_type])]
    gff_df["seqid"] = gff_df["seqid"].astype(str)
    gff_df["score"] = pd.to_numeric(gff_df["score"], errors="coerce")
    gff_df["score"] = gff_df["score"].fillna(0)

    # Apply the mapping function in parallel
    gff_dataset = Dataset.from_pandas(gff_df)

    location_map = {}
    processed_dataset = gff_dataset.map(
        process_attribute,
        num_proc=num_proc,
        desc="Processing attributes in parallel",
        fn_kwargs={"key": gff_attribute_key, "preface": gff_id_preface},
    )

    # Convert back to a dictionary for location_map
    for row in tqdm(processed_dataset):
        location_map[row["id"]] = {
            "chromosome": row["chromosome"],
            "start": row["start"],
            "end": row["end"],
            "strand": row["strand"],
        }

    # %% Prepare output data structure
    out_columns = [
        "id",
        "sequence",
        "length",
        "species",
        "genus",
        "family",
        "order",
        "class",
        "phylum",
        "domain",
        "chromosome",
        "start",
        "end",
        "strand",
        "cluster_id",
        "cluster_name",
        "cluster_description",
        "pfam_domains",  # semicolon separated
        "go_terms",  # semicolon separated
        "direct_interactions",  # semicolon separated
        "cluster_interactions",  # semicolon separated
    ]

    protein_data: dict[str, dict] = {}

    # %% Process GO and Pfam data
    logger.info("Processing GO and Pfam data")
    for _, row in tqdm(go_map.iterrows(), total=len(go_map)):
        protein_id = row["prot_id"]
        protein_data.setdefault(protein_id, {})
        protein_data[protein_id]["pfam_domains"] = row["pfam_list"]
        protein_data[protein_id]["go_terms"] = row["GO_list"]

    logger.info("Processing protein sequences")
    for protein_id, sequence in tqdm(sequence_dict.items(), total=len(sequence_dict)):
        protein_data.setdefault(protein_id, {})
        protein_data[protein_id]["sequence"] = str(sequence)
        protein_data[protein_id]["length"] = len(sequence)

    logger.info("Adding genomic location data")
    for protein_id in tqdm(protein_data.keys(), total=len(protein_data)):
        loc_info = location_map.get(
            protein_id,
            {"chromosome": None, "start": None, "end": None, "strand": None},
        )
        protein_data[protein_id]["chromosome"] = loc_info["chromosome"]
        protein_data[protein_id]["start"] = loc_info["start"]
        protein_data[protein_id]["end"] = loc_info["end"]
        protein_data[protein_id]["strand"] = loc_info["strand"]

    logger.info("Processing Philharmonic cluster data")
    for cluster_id, cluster in tqdm(
        philharmonic_data.items(), total=len(philharmonic_data)
    ):
        for protein_id in cluster["members"]:
            protein_data.setdefault(protein_id, {})
            protein_data[protein_id]["cluster_id"] = str(cluster_id)
            protein_data[protein_id]["cluster_name"] = cluster["llm_name"]
            protein_data[protein_id]["cluster_description"] = cluster["llm_explanation"]

            direct_interactions = []
            for source, target, _ in cluster["graph"]:
                if source == protein_id:
                    direct_interactions.append(target)
                elif target == protein_id:
                    direct_interactions.append(source)
            if len(direct_interactions):
                protein_data[protein_id]["direct_interactions"] = ";".join(
                    direct_interactions
                )

            protein_data[protein_id]["cluster_interactions"] = ";".join(
                [i for i in cluster["members"] if i != protein_id]
            )

    logger.info("Adding taxonomy information")
    for protein_id, protein_info in tqdm(protein_data.items(), total=len(protein_data)):
        protein_data[protein_id]["species"] = species
        protein_data[protein_id]["genus"] = genus
        protein_data[protein_id]["family"] = family
        protein_data[protein_id]["order"] = order
        protein_data[protein_id]["class"] = class_name
        protein_data[protein_id]["phylum"] = phylum
        protein_data[protein_id]["domain"] = domain

    logger.info("Converting protein data to DataFrame and saving as CSV")
    df = pd.DataFrame.from_dict(protein_data, orient="index")
    df.index.name = "id"
    df = df.reset_index()

    # Ensure all required columns exist, filling with NaN if missing
    for col in out_columns:
        if col not in df.columns:
            df[col] = pd.NA

    df = df[out_columns]
    df.to_csv(outfile, index=False)

    logger.info(f"Processed {len(df)} proteins.")
    logger.info(
        f"Includes {len(df[~df.cluster_id.isna()])} proteins with cluster information."
    )
    logger.info(
        f"Includes {len(df[~df.pfam_domains.isna()])} proteins with Pfam domain information."
    )
    logger.info(f"Includes {len(df[~df.sequence.isna()])} proteins with sequences.")
    logger.info(f"Includes {len(df[~df.start.isna()])} proteins with genome location.")
    logger.info(f"Data saved to {outfile}")


if __name__ == "__main__":
    typer.run(main)
