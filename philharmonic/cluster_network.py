# python src/cluster_network.py --network_file {input.network} --dsd_file {input.distances} --output {output.clusters} --min_cluster_size {config[clustering][min_cluster_size]} --cluster_divisor {config[clustering][cluster_divisor]} --init_k {config[clustering][init_k]}

import json
import typing as T
from queue import PriorityQueue

import networkx as nx
import numpy as np
import pandas as pd
import typer
from loguru import logger
from sklearn.cluster import SpectralClustering

from .utils import hash_cluster

app = typer.Typer()


def RBF(D, sigma=None):
    """
    Convert distance matrix D into similarity matrix S using Radial Basis Function (RBF) Kernel
    RBF(x,x') = exp( -((x - x')**2 / 2sigma**@))
    """
    sigma = sigma or np.sqrt(np.max(D))
    return np.exp(-1 * (np.square(D) / (2 * sigma**2)))


def dsd_to_similarity(D, sparsity_threshold=1e-5):
    sim = RBF(D)
    sim[sim < sparsity_threshold] = 0
    return sim


def extract_clusters(
    SC: SpectralClustering, og_clust: list, og_csize: int, verbosity=2
):
    subCs = []
    for label in set(SC.labels_):  # go through all labels
        # make a list of all nodes in the cluster assigned to given label
        subCs.append([og_clust[i] for i in range(og_csize) if SC.labels_[i] == label])

    if verbosity > 0:
        to_print = f"Cluster of size {og_csize} -> split into {len(subCs)} clusters "

        if verbosity > 1:
            to_print += "of sizes: "
            for subC in subCs:
                to_print += f"{len(subC)}/"
        logger.info(to_print[:-1])

    # return list of tuples of priority of cluster and list of nodes in cluster
    return [(1 / len(subC), subC) for subC in subCs]


def read_network_subset(network_file, protein_names, output_stats=True):
    """
    Read in the network and filter edges based on DSD confidence threshold
    """
    fullG = nx.read_weighted_edgelist(network_file)
    logger.info("Selecting subset proteins...")
    G = fullG.subgraph(protein_names)

    if output_stats:
        # print a table of network statistics
        label = ["Nodes", "Edges", "Degree (Med)", "Degree (Avg)", "Sparsity"]
        degrees = [i[1] for i in list(G.degree())]
        value = [
            len(G.nodes),
            len(G.edges),
            np.median(degrees),
            np.mean(degrees),
            len(G.edges()) / len(G) ** 2,
        ]
        stats = pd.DataFrame([label, value]).T
        stats.columns = ["", "Value"]
        stats = stats.set_index("")
        logger.info(stats)
        # save a histogram of protein degree

    return G


@app.command()
def main(
    init_k: int = typer.Option(500, "--init_k", help="Initial number of clusters"),
    min_cluster_size: int = typer.Option(
        3, "--min_cluster_size", help="Minimum cluster size"
    ),
    cluster_divisor: int = typer.Option(
        20, "--cluster_divisor", help="Cluster divisor"
    ),
    sparsity: float = typer.Option(1e-5, "--sparsity", help="Sparsity threshold"),
    random_seed: int = typer.Option(42, "--random_seed", help="Random seed"),
    output: str = typer.Option(..., "-o", "--output", help="Output file"),
    dsd_file: str = typer.Option(..., "--dsd_file", help="Distances file"),
    network_file: str = typer.Option(..., "--network_file", help="Network file"),
):
    """Cluster a network using a DSD file"""
    logger.info(f"Reading DSD file: {dsd_file}")
    dsd_df = pd.read_csv(dsd_file, sep="\t", index_col=0, header=0)
    protein_names = [str(i) for i in dsd_df.index]
    dsd = dsd_df.values

    logger.info(f"Reading network file: {network_file}")
    G = read_network_subset(network_file, protein_names)

    logger.info(
        f"Converting DSD to similarity matrix with sparsity threshold: {sparsity}"
    )
    sim = dsd_to_similarity(dsd, sparsity_threshold=sparsity)

    logger.info(f"Fitting {init_k} spectral clusters...")
    SC = SpectralClustering(
        n_clusters=init_k,
        assign_labels="discretize",
        random_state=random_seed,
        affinity="precomputed",
    )
    SC.fit(sim)

    clusts = [
        [j for j in range(len(SC.labels_)) if SC.labels_[j] == i]
        for i in range(max(SC.labels_) + 1)
        if i in SC.labels_
    ]
    clusts.sort(key=lambda x: len(x), reverse=True)

    logger.info(f"Initial Clustering: {len(clusts)} clusters")
    clustQ: PriorityQueue[T.Tuple[float, T.List[int]]] = PriorityQueue()
    for c in clusts:
        clustQ.put((1 / len(c), c))

    logger.info(f"Splitting large clusters into {cluster_divisor} clusters...")
    while True:
        priority, c = clustQ.get()
        csize = int(1 / priority)

        n_clusters = int(
            np.round(csize / cluster_divisor)
        )  # NOTE: this changed so that we are left with clusters of max size 30

        if n_clusters < 2:
            clustQ.put((priority, c))
            break

        SC2 = SpectralClustering(
            n_clusters=n_clusters,
            assign_labels="discretize",
            random_state=random_seed,
            affinity="precomputed",
        )
        SC2.fit(sim[c, :][:, c])

        subClusts = extract_clusters(SC2, c, csize, verbosity=2)
        for subClust in subClusts:
            clustQ.put(subClust)

    logger.info(f"Removing small clusters (<{min_cluster_size})...")
    filtered_clusters = []
    while not clustQ.empty():
        _, c = clustQ.get()
        filtered_clusters.append(c)
    filtered_clusters = [i for i in filtered_clusters if len(i) >= min_cluster_size]
    logger.info(f"Final Clustering: {len(filtered_clusters)} clusters")

    clustsNames = [[protein_names[i] for i in cl] for cl in filtered_clusters]

    clustsDict = {}
    for cl in clustsNames:
        uqid = hash_cluster(cl)
        members = cl
        subgraph = G.subgraph(cl)
        clustsDict[uqid] = {
            "members": members,
            "graph": [
                (i[0], i[1], subgraph.get_edge_data(i[0], i[1])["weight"])
                for i in subgraph.edges()
            ],
        }

    logger.info(f"Writing clusters to: {output}")
    with open(output, "w+") as f:
        json.dump(clustsDict, f, indent=4)


if __name__ == "__main__":
    app()
