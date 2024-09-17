# python src/add_cluster_functions.py -o {output.clusters_functional} -cfp {input.clusters} --go_map {input.go_map}
import argparse
import json

from utils import parse_GO_map, load_cluster_json, log


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


if __name__ == "__main__":
    parser = argparse.ArgumentParser("Add GO terms to clusters")
    parser.add_argument("-o", "--output", required=True, type=str, help="Output file")
    parser.add_argument(
        "-cfp",
        "--clusters_file_path",
        required=True,
        type=str,
        help="Clusters file path",
    )
    parser.add_argument("--go_map", required=True, type=str, help="GO map file")

    args = parser.parse_args()

    # Load clusters
    clusters = load_cluster_json(args.clusters_file_path)

    # Add GO Annotations
    go_map = parse_GO_map(args.go_map)

    log("Adding GO Annotations...")
    for k, clust in clusters.items():
        clusters[k]["GO_terms"] = add_GO_function(clust, go_map)

    with open(args.output, "w") as f:
        json.dump(clusters, f, indent=4)
