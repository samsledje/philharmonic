# "python src/philharmonic_cytoscape.py -s {input.styles} -o {output} {input.clusters}

import typer
import py4cytoscape as p4c
import networkx as nx
from .utils import load_cluster_json

app = typer.Typer()


@app.command()
def main(
    styles: str, output: str, cluster_file_path: str, network_file_path: str, name: str
):
    """Build a Cytoscape network from a cluster file and a network file"""
    netx = nx.read_weighted_edgelist(network_file_path)
    clusters = load_cluster_json(cluster_file_path)

    p4c.create_network_from_networkx(netx, title=f"{name} Main", collection=name)

    # cluster_group_map = {}
    for k, clust in clusters.items():
        # p4c.create_network_from_networkx(netx.subgraph(clust["members"]), title=f"{args.name} Cluster {k}", collection=args.name)
        p4c.create_group(
            nodes=clust["members"], group_name=f"Cluster {k}", nodes_by_col="name"
        )

    p4c.collapse_group("all")
    p4c.save_session(output)


if __name__ == "__main__":
    app()
