import sys
import argparse
import logging
import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
import pandas as pd
import pickle as pk
import json
import hashlib
import collections
from tqdm.notebook import tqdm
from sklearn.cluster import KMeans
from sklearn.cluster import SpectralClustering

from philharmonic_cluster import (
    Cluster,
    readClusterObjects,
    writeClusters,
    cluster_from_json,
)

from add_cluster_functions import (
    GO,
    extract_GO_id_from_list,
    clean_GO_map,
    read_GO_obo,
    GO_search,
    protein_search,
    triangle_search,
    node_search,
)

import numpy as np
import pandas as pd


class GO:
    def __init__(self, ID, features):
        self.ID = ID
        self.D = features
        self.name = features["name"]

    def __repr__(self):
        return "{} - <{}>".format(self.ID, self.name)

    def __eq__(self, other):
        return self.ID == other.ID

    def __hash__(self):
        return hash(self.ID)


def extract_GO_id_from_list(l):
    if isinstance(l, list):
        return [i.split("|")[0] for i in l]
    else:
        return None


def clean_GO_map(f):
    seqDb = pd.read_csv(f, sep=",")
    seqDb.columns = ["seq", "manual_annot", "pfam_list", "GO_list"]
    seqDb["GO_str"] = seqDb["GO_list"]
    seqDb["GO_list"] = seqDb["GO_str"].str.split(";")

    def extract_GO_id_from_list(l):
        if isinstance(l, list):
            return [i.split("|")[0] for i in l]
        else:
            return None

    seqDb["GO_ids"] = seqDb["GO_list"].apply(extract_GO_id_from_list)
    seq2GO = seqDb[["seq", "GO_ids"]]
    seq2GO.columns = ["seq", "GO_ids"]
    return seq2GO


def read_GO_obo(infile):
    terms = {}
    with open(infile, "r") as f:
        for line in f:
            tDict = {}
            line = line.strip()
            if line == "[Term]":
                line = f.readline().strip().split(": ")
                while not line == [""]:
                    tDict[line[0]] = "".join(line[1:])
                    line = f.readline().strip().split(": ")
                for k, v in tDict.items():
                    k = k.strip()
                    v = v.strip()
                    tDict[k] = v
                terms[tDict["id"]] = GO(tDict["id"], tDict)
    return terms


def GO_search(clusters, GO_term, GO_OBJECTS, N=20):
    if isinstance(GO_term, str):
        GO_term = GO_OBJECTS[GO_term]
    return [c for c in clusters if GO_term in [i[0] for i in c.get_top_terms(N)]]


def protein_search(clusters, protein_list):
    if isinstance(protein_list, str):
        protein_list = [protein_list]
    plist = [
        c for c in clusters if np.array([p in c.proteins for p in protein_list]).all()
    ]
    if len(plist) == 0:
        return None
    elif len(plist) == 1:
        return plist[0]
    else:
        return plist


def triangle_search(clusters, min_triangles=0, max_triangles=np.infty):
    return [
        c
        for c in clusters
        if c.triangles() >= min_triangles and c.triangles() <= max_triangles
    ]


def node_search(clusters, min_nodes=0, max_nodes=np.infty):
    return [c for c in clusters if len(c) >= min_nodes and len(c) <= max_nodes]


def RBF(D, sigma=None):
    """
    Convert distance matrix D into similarity matrix S using Radial Basis Function (RBF) Kernel
    RBF(x,x') = exp( -((x - x')**2 / 2sigma**@))
    """
    sigma = sigma or np.sqrt(np.max(D))
    return np.exp(-1 * (np.square(D) / (2 * sigma**2)))


if __name__ == "__main__":
    parser = argparse.ArgumentParser("PHILHARMONIC network clustering")

    parser.add_argument("--network-id", type=str, required=True)
    parser.add_argument("--network-file", type=str, required=True)
    parser.add_argument("--dsd-file", type=str, required=True)
    parser.add_argument("--go-map", type=str, required=True)

    parser.add_argument("--log-file", type=str, default="philharmonic.log")
    parser.add_argument("--edge-weight-thresh", default=0.5, type=float)
    parser.add_argument("--sparse-similarity-thresh", default=1e-5, type=float)
    parser.add_argument("--k-clusts", default=500, type=int)
    parser.add_argument("--min-cluster-size", default=3, type=int)
    parser.add_argument("--max-cluster-size", default=100, type=int)

    parser.add_argument("--random-state", default=6191998, type=int)
    parser.add_argument("--go-db", default="go.obo", type=str)

    args = parser.parse_args()

    network_id = args.network_id
    network_file = args.network_file
    dsd_file = args.dsd_file
    go_map = args.go_map
    go_db = args.go_db
    log_file = args.log_file

    edge_weight_thresh = args.edge_weight_thresh
    k_clusts = args.k_clusts
    min_cluster_size = args.min_cluster_size
    max_cluster_size = args.max_cluster_size
    sparse_sim_thresh = args.sparse_similarity_thresh
    random_state = args.random_state

    logg = logging.getLogger(network_id)
    logg.setLevel(logging.DEBUG)
    formatter = logging.Formatter(
        "%(asctime)s - %(name)s - %(levelname)s - %(message)s"
    )

    streamHandler = logging.StreamHandler(sys.stdout)
    streamHandler.setLevel(logging.DEBUG)
    streamHandler.setFormatter(formatter)
    logg.addHandler(streamHandler)

    if log_file is not None:
        fileHandler = logging.FileHandler(log_file)
        fileHandler.setFormatter(formatter)
        logg.addHandler(fileHandler)

    logg.info(args)

    # Read Network and DSD File
    logg.info("Reading DSD File...")
    df = pd.read_csv(dsd_file, sep="\t", index_col=0, header=0)
    protein_names = [str(i) for i in df.index]
    DSD = df.values

    fullG = nx.read_weighted_edgelist(network_file)
    logg.info("Selecting DSD connected component...")
    G = fullG.subgraph(protein_names)
    logg.info(
        "Filtering edges with confidence threshold {}...".format(edge_weight_thresh)
    )
    wG = nx.Graph()
    for u, v, d in tqdm(G.edges.data()):
        if d["weight"] >= edge_weight_thresh:
            wG.add_edge(u, v, weight=d["weight"])
    del G
    G = wG
    A = nx.to_numpy_array(G, nodelist=protein_names)
    degrees = [i[1] for i in list(G.degree())]

    # Print Network Statistics
    label = ["Nodes", "Edges", "Degree (Med)", "Degree (Avg)", "Sparsity"]
    value = [
        len(G.nodes),
        len(G.edges),
        np.median(degrees),
        np.mean(degrees),
        len(G.edges()) / len(G) ** 2,
    ]
    df = pd.DataFrame([label, value]).T
    df.columns = ["", network_id]
    df = df.set_index("")

    logg.info(df.to_markdown())

    # DSD --> Similarity
    logg.info("Computing similarity scores...")
    simDSD = RBF(DSD)

    logg.info("Sparsifying similarity scores...")
    simRav = simDSD.ravel()
    simRav[simRav < sparse_sim_thresh] = 0
    simRav = simRav.reshape(simDSD.shape)
    simDSD = simRav

    # Spectral Clustering
    logg.info("Fitting spectral clusters...")
    SC = SpectralClustering(
        n_clusters=k_clusts,
        assign_labels="discretize",
        random_state=random_state,
        affinity="precomputed",
    )
    SC.fit(simDSD)

    # Filter and Split Clusters
    clusts = [
        [j for j in range(len(SC.labels_)) if SC.labels_[j] == i]
        for i in range(max(SC.labels_) + 1)
        if i in SC.labels_
    ]
    clusts.sort(key=lambda x: len(x), reverse=True)
    logg.info(f"Initial clustering results in {len(clusts)} clusters")

    from queue import PriorityQueue

    clustQ = PriorityQueue()

    for c in clusts:
        clustQ.put((1 / len(c), c))

    logg.info(f"Splitting large clusters (>{max_cluster_size})...")
    while True:
        c = clustQ.get()
        csize = int(1 / c[0])
        if csize <= max_cluster_size:
            break
        SC2 = SpectralClustering(
            n_clusters=2,
            assign_labels="discretize",
            random_state=random_state,
            affinity="precomputed",
        )
        SC2.fit(simMatrix[c[1], :][:, c[1]])
        subC_0 = [c[1][i] for i in range(csize) if SC2.labels_[i] == 0]
        subC_1 = [c[1][i] for i in range(csize) if SC2.labels_[i] == 1]
        logg.info("{} -> {}/{}".format(csize, len(subC_0), len(subC_1)))
        clustQ.put((1 / len(subC_0), subC_0))
        clustQ.put((1 / len(subC_1), subC_1))

    logg.info(f"Removing small clusters (<{min_cluster_size})...")
    filteredClusters = []
    while not clustQ.empty():
        filteredClusters.append(clustQ.get()[1])
    filteredClusters = [i for i in filteredClusters if len(i) >= min_cluster_size]

    clustsNames = [[protein_names[i] for i in cl] for cl in filteredClusters]
    clustsNames.sort(key=lambda x: len(x), reverse=True)

    # Write and Load Clusters
    writeClusters(f"{network_id}.clusters.csv", clustsNames)
    clusters = readClusterObjects(f"{network_id}.clusters.csv")
    logg.info(f"Final number of clusters: {len(clusters)}")

    for i in clusters:
        i.set_graph(G)

    # Add GO Annotations
    goMap = clean_GO_map(go_map)
    goDB = read_GO_obo(go_db)

    logg.info("Adding GO Annotations...")
    for clust in tqdm(clusters):
        clust.add_GO_terms(goMap, goDB)

    clusters.sort(key=lambda x: len(x), reverse=True)

    # Save Clusters
    logg.info("Saving clusters...")
    pk.dump(clusters, open(f"{network_id}.clusters.pk", "wb"))
    clusterDictJson = {hash(i): i.to_dict() for i in clusters}
    with open(f"{network_id}.clusters.json", "w+") as f:
        f.write(json.dumps(clusterDictJson))
