# python src/add_cluster_functions.py -o {output.clusters_functional} -cfp {input.clusters} --go_map {input.go_map}

import typer
from pathlib import Path
import json
from loguru import logger

from .utils import parse_GO_map, load_cluster_json


def add_GO_function(cluster, go_map):
    """
    Keep track of how many proteins in the cluster have a given GO term
    """
    go_terms = {}
    for protein in cluster["members"]:
        if protein in go_map:
            for gid in set(go_map[protein]):
                go_terms[gid] = go_terms.setdefault(gid, 0) + 1
    return go_terms

app = typer.Typer()

@app.command()
def main(
    output: Path = typer.Option(..., "-o", "--output", help="Output file"),
    clusters_file_path: Path = typer.Option(..., "-cfp", "--clusters-file-path", help="Clusters file path"),
    go_map: Path = typer.Option(..., "--go-map", help="GO map file"),
):
    """Add GO terms to clusters"""
    # Load clusters
    clusters = load_cluster_json(clusters_file_path)

    # Add GO Annotations
    go_map_data = parse_GO_map(go_map)

    logger.info("Adding GO Annotations...")
    for k, clust in clusters.items():
        clusters[k]["GO_terms"] = add_GO_function(clust, go_map_data)

    with open(output, "w") as f:
        json.dump(clusters, f, indent=4)

if __name__ == "__main__":
    app()