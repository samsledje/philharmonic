from collections import Counter

import networkx as nx
import pandas as pd
import typer
from loguru import logger
from tqdm import tqdm
from typing_extensions import Annotated

from .utils import (
    add_GO_function,
    clean_top_terms,
    load_cluster_json,
    nx_graph_cluster,
    parse_GO_database,
    parse_GO_map,
)

app = typer.Typer()


@app.command()
def main(
    output: str = typer.Option(..., "-o", "--output", help="Output file"),
    coc_functions: str = typer.Option(
        ..., "-coc", "--coc_functions", help="Function of clusters of clusters"
    ),
    cluster_file_path: str = typer.Option(
        ..., "-cfp", "--cluster_file_path", help="Cluster file path"
    ),
    network_file_path: str = typer.Option(
        ..., "-nfp", "--network_file_path", help="Network file path"
    ),
    go_map_path: str = typer.Option(..., "--go_map", help="GO map file"),
    go_db_path: str = typer.Option(..., "--go_db", help="GO database file"),
    n_crossing_edges: int = typer.Option(
        10, "--n_crossing_edges", help="Clusters with n crossing edges are connected"
    ),
    include_recipe: Annotated[
        bool, typer.Option(..., "--include_recipe", help="Include recipe")
    ] = False,
    recipe_metric: str = typer.Option(
        "degree", "--recipe_metric", help="Recipe metric"
    ),
    recipe_cthresh: str = typer.Option(
        "0.75", "--recipe_cthresh", help="Recipe threshold"
    ),
):
    """Build a graph of clusters"""
    # Load clusters
    cluster_dict = load_cluster_json(cluster_file_path)

    # Load network
    full_G = nx.read_weighted_edgelist(network_file_path)

    # Add GO Annotations
    go_map = parse_GO_map(go_map_path)
    go_db = parse_GO_database(go_db_path)

    for clust in cluster_dict.values():
        clust["GO_terms"] = add_GO_function(clust, go_map, go_db=go_db)
        for gt in clust["GO_terms"].keys():
            assert gt in go_db.keys()

    cluster_top_terms = {
        k: clean_top_terms(c, go_db, return_counts=False)
        for k, c in cluster_dict.items()
    }
    counter = Counter(cluster_top_terms.values())
    logger.info(counter)

    def test_cluster_connection(c1_graph, c2_graph, parent_graph):
        return len(list(nx.edge_boundary(parent_graph, c1_graph.nodes, c2_graph.nodes)))

    def edge_filter(x1, x2, n_crossing):
        return clusG[x1][x2].get("weight") > n_crossing

    clusG = nx.Graph()
    for ckey, cdict in cluster_dict.items():
        clusG.add_node(ckey, size=len(cdict["members"]))

    cgraphs = {
        k: nx_graph_cluster(
            v,
            use_recipe_nodes=include_recipe,
            full_G=full_G,
            recipe_metric=recipe_metric,
            recipe_cthresh=recipe_cthresh,
        )
        for k, v in cluster_dict.items()
    }
    for n in tqdm(clusG.nodes):
        for n2 in clusG.nodes:
            if n != n2:
                n_neighbors = test_cluster_connection(cgraphs[n], cgraphs[n2], full_G)
                clusG.add_edge(n, n2, weight=n_neighbors)

    subG = nx.subgraph_view(
        clusG, filter_edge=lambda x1, x2: edge_filter(x1, x2, n_crossing_edges)
    )
    logger.info(
        f"Cluster graph has {subG.number_of_nodes()} nodes and {subG.number_of_edges()} edges."
    )

    edge_table = nx.to_pandas_edgelist(subG)

    # Nodes
    size_table = pd.DataFrame(
        [(k, len(cluster_dict[k]["members"])) for k in subG.nodes],
        columns=["key", "size"],
    )
    fn_table = pd.DataFrame(cluster_top_terms.items(), columns=["key", "go_fn"])
    node_table = pd.merge(size_table, fn_table)

    n_singletons = len(subG.nodes) - len(
        set(edge_table["source"].values).union(edge_table["target"].values)
    )
    logger.info(
        f"There are {n_singletons} singleton clusters that don't have more than {n_crossing_edges} edges to any other cluster"
    )

    edge_table.to_csv(output, index=False, sep="\t")
    node_table.to_csv(coc_functions, index=False, sep="\t")


if __name__ == "__main__":
    app()
