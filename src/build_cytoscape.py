# "python src/philharmonic_cytoscape.py -s {input.styles} -o {output} {input.clusters}

import argparse
import pandas
import py4cytoscape as p4c
import networkx as nx
from pathlib import Path
from utils import load_cluster_json

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Create a Cytoscape style file from a cluster file')
    parser.add_argument('-s', '--styles', type=str, help='Cytoscape style file')
    parser.add_argument('-o', '--output', type=str, help='Output file')
    parser.add_argument('--name', type=str, help='Name of the network')
    parser.add_argument('-cfp', '--cluster_file_path', type=str, help='Cluster file path')
    parser.add_argument('-nfp', '--network_file_path', type=str, help='Node file path')

    args = parser.parse_args()

    netx = nx.read_weighted_edgelist(args.network_file_path)
    clusters = load_cluster_json(args.cluster_file_path)

    p4c.create_network_from_networkx(netx, title=f"{args.name} Main", collection=args.name)

    cluster_group_map = {}
    for k, clust in clusters.items():
        # p4c.create_network_from_networkx(netx.subgraph(clust["members"]), title=f"{args.name} Cluster {k}", collection=args.name)
        p4c.create_group(nodes=clust["members"], group_name=f"Cluster {k}", nodes_by_col="name")

    p4c.collapse_group('all')

    p4c.save_session(args.output)
    # Path(args.output).touch()