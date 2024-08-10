# python src/cluster_network.py --network_file {input.network} --dsd_file {input.distances} --output {output.clusters} --min_cluster_size {config[clustering][min_cluster_size]} --cluster_divisor {config[clustering][cluster_divisor]} --init_k {config[clustering][init_k]}

import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import networkx as nx
import argparse
import json
from sklearn.cluster import SpectralClustering
from queue import PriorityQueue

from utils import log, hash_cluster


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
        log(to_print[:-1])

    # return list of tuples of priority of cluster and list of nodes in cluster
    return [(1 / len(subC), subC) for subC in subCs]


def writeClusters(outfile, clusts):
    with open(outfile, "w+") as f:
        json.dump(clusts, f, indent=4)


def read_network_subset(network_file, protein_names, output_stats=True):
    """
    Read in the network and filter edges based on DSD confidence threshold
    """
    fullG = nx.read_weighted_edgelist(network_file)
    log("Selecting subset proteins...")
    G = fullG.subgraph(protein_names)
    degrees = [i[1] for i in list(G.degree())]

    if output_stats:
        # print a table of network statistics
        label = ["Nodes", "Edges", "Degree (Med)", "Degree (Avg)", "Sparsity"]
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
        log(stats)
        # save a histogram of protein degree

    return G


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Cluster PPI network")
    parser.add_argument(
        "--init_k", type=int, default=500, help="Initial number of clusters"
    )
    parser.add_argument(
        "--min_cluster_size", type=int, default=3, help="Minimum cluster size"
    )
    parser.add_argument(
        "--cluster_divisor", type=int, default=20, help="Cluster divisor"
    )
    parser.add_argument(
        "--sparsity", type=float, default=1e-5, help="Sparsity threshold"
    )
    parser.add_argument("--random_seed", type=int, default=42, help="Random seed")
    parser.add_argument("-o", "--output", type=str, required=True, help="Output file")
    parser.add_argument("--dsd_file", type=str, required=True, help="Distances file")
    parser.add_argument("--network_file", type=str, required=True, help="Network file")

    args = parser.parse_args()

    log(f"Reading DSD file: {args.dsd_file}")
    dsd_df = pd.read_csv(args.dsd_file, sep="\t", index_col=0, header=0)
    protein_names = [str(i) for i in dsd_df.index]
    dsd = dsd_df.values

    log(f"Reading network file: {args.network_file}")
    G = read_network_subset(args.network_file, protein_names)

    log(f"Converting DSD to similarity matrix with sparsity threshold: {args.sparsity}")
    sim = dsd_to_similarity(dsd, sparsity_threshold=args.sparsity)

    log(f"Fitting {args.init_k} spectral clusters...")
    SC = SpectralClustering(
        n_clusters=args.init_k,
        assign_labels="discretize",
        random_state=args.random_seed,
        affinity="precomputed",
    )
    SC.fit(sim)

    clusts = [
        [j for j in range(len(SC.labels_)) if SC.labels_[j] == i]
        for i in range(max(SC.labels_) + 1)
        if i in SC.labels_
    ]
    clusts.sort(key=lambda x: len(x), reverse=True)

    log(f"Initial Clustering: {len(clusts)} clusters")
    clustQ = PriorityQueue()
    for c in clusts:
        clustQ.put((1 / len(c), c))

    log(f"Splitting large clusters into {args.cluster_divisor} clusters...")
    while True:
        priority, c = clustQ.get()
        csize = int(1 / priority)

        n_clusters = int(
            np.round(csize / args.cluster_divisor)
        )  # NOTE: this changed so that we are left with clusters of max size 30

        if n_clusters < 2:
            clustQ.put((priority, c))
            break

        SC2 = SpectralClustering(
            n_clusters=n_clusters,
            assign_labels="discretize",
            random_state=args.random_seed,
            affinity="precomputed",
        )
        SC2.fit(sim[c, :][:, c])

        subClusts = extract_clusters(SC2, c, csize, verbosity=2)
        for subClust in subClusts:
            clustQ.put(subClust)

    log(f"Removing small clusters (<{args.min_cluster_size})...")
    filteredClusters = []
    while not clustQ.empty():
        wght, c = clustQ.get()
        filteredClusters.append(c)
    filteredClusters = [i for i in filteredClusters if len(i) >= args.min_cluster_size]
    log(f"Final Clustering: {len(filteredClusters)} clusters")

    clustsNames = [[protein_names[i] for i in cl] for cl in filteredClusters]

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

    # clustsDict = {hash_cluster(cl): {"members": cl} for cl in clustsNames}

    log(f"Writing clusters to: {args.output}")
    writeClusters(args.output, clustsDict)
