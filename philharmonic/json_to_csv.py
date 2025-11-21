# %% Imports
import re
from multiprocessing import cpu_count

import hydra
import pandas as pd
from biotite.sequence.io import fasta
from datasets import Dataset
from loguru import logger
from omegaconf import DictConfig
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


# %% Main function with Hydra configuration
@hydra.main(config_path="configs")
def main(cfg: DictConfig):
    logger.info(f"Loading Philharmonic data from JSON: {cfg.JSON_FILE}")
    philharmonic_data = load_cluster_json(cfg.JSON_FILE)

    logger.info(f"Loading GO mapping from CSV: {cfg.GO_FILE}")
    go_map = pd.read_csv(cfg.GO_FILE)

    logger.info(f"Loading sequences from FASTA: {cfg.FASTA_FILE}")
    fasta_fh = fasta.FastaFile.read(cfg.FASTA_FILE)
    sequence_dict = {
        seq_id.split()[0]: seq for seq_id, seq in fasta.get_sequences(fasta_fh).items()
    }

    logger.info(f"Loading GFF data from file: {cfg.GFF_FILE}")
    gff_df = pd.read_csv(
        cfg.GFF_FILE,
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
    gff_df = gff_df[gff_df["type"].isin([cfg.GFF_FEATURE_TYPE])]
    gff_df["seqid"] = gff_df["seqid"].astype(str)
    gff_df["score"] = pd.to_numeric(gff_df["score"], errors="coerce")
    gff_df["score"] = gff_df["score"].fillna(0)

    # Apply the mapping function in parallel
    gff_dataset = Dataset.from_pandas(gff_df)

    location_map = {}
    processed_dataset = gff_dataset.map(
        process_attribute,
        num_proc=cpu_count()
        - 4,  # Adjust the number of processes based on your CPU cores
        desc="Processing attributes in parallel",
        fn_kwargs={"key": cfg.GFF_ATTRIBUTE_KEY, "preface": cfg.GFF_ID_PREFACE},
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
            for source, target, confidence in cluster["graph"]:
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
        protein_data[protein_id]["species"] = cfg.SPECIES
        protein_data[protein_id]["genus"] = cfg.GENUS
        protein_data[protein_id]["family"] = cfg.FAMILY
        protein_data[protein_id]["order"] = cfg.ORDER
        protein_data[protein_id]["class"] = cfg.CLASS
        protein_data[protein_id]["phylum"] = cfg.PHYLUM
        protein_data[protein_id]["domain"] = cfg.DOMAIN

    logger.info("Converting protein data to DataFrame and saving as CSV")
    df = pd.DataFrame.from_dict(protein_data, orient="index")
    df.index.name = "id"
    df = df.reset_index()
    df = df[out_columns]
    df.to_csv(cfg.OUTFILE, index=False)

    logger.info(f"Processed {len(df)} proteins.")
    logger.info(
        f"Includes {len(df[~df.cluster_id.isna()])} proteins with cluster information."
    )
    logger.info(
        f"Includes {len(df[~df.pfam_domains.isna()])} proteins with Pfam domain information."
    )
    logger.info(f"Includes {len(df[~df.sequence.isna()])} proteins with sequences.")
    logger.info(f"Includes {len(df[~df.start.isna()])} proteins with genome location.")
    logger.info(f"Data saved to {cfg.OUTFILE}")


if __name__ == "__main__":
    main()
