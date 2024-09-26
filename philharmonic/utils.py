import hashlib
import json
import numpy as np
import pandas as pd
import networkx as nx
from matplotlib import pyplot as plt
import seaborn as sns
from Bio import SeqIO


class Cluster:
    def __init__(self, cluster_dict):
        for k, v in cluster_dict.items():
            setattr(self, k, v)
        assert hasattr(self, "members"), "Cluster must have 'members' attribute"

        if hasattr(self, "graph"):
            self.G = nx.Graph()
            for edge in self.graph:
                self.G.add_edge(edge[0], edge[1], weight=edge[2])

    def __len__(self):
        return len(self.members)

    def __repr__(self):
        reprStr = "Cluster of {} [{}] (hash {})".format(
            len(self), self._member_list(), hash(self)
        )
        if hasattr(self, "recipe"):
            pass
        if hasattr(self, "G"):
            reprStr += "\nTriangles: {}\nMax Degree: {}".format(
                self.triangles(), max(self.G.degree(), key=lambda x: x[1])[1]
            )
        if hasattr(self, "GO_terms"):
            reprStr += "\nTop Terms:\n\t{}".format(
                "\n\t".join(
                    # ['{} ({})'.format(i[0], i[1]) for i in self.get_top_terms(5)]
                    [f"{i[0]} - <{self.go_db[i[0]]}> ({i[1]})" for i in self.get_top_terms(5)]
                )
            )
        return reprStr

    def __hash__(self):
        return hash_cluster(self.members)

    def _member_list(self):
        if len(self.members) < 3:
            return ", ".join(self.members)
        else:
            return "{}, ...".format(", ".join(self.members[:3]))

    def add_GO_terms(self, go_map, go_db):
        self.GO_terms = {}
        self.go_db = go_db.copy()
        for prot in self.members:
            goIds = go_map.get(prot, None)
            if goIds is None or len(goIds) == 0:
                continue
            for gid in goIds:
                if gid in go_db:
                    goCount = self.GO_terms.setdefault(gid, 0)
                    self.GO_terms[gid] = goCount + 1

    def set_graph(self, G):
        self.G = G.subgraph(self.members)

    def draw_graph(self):
        if not hasattr(self, "G"):
            raise ValueError("Run .set_graph() method on this cluster first")
        G = self.G
        nx.draw_kamada_kawai(G, with_labels=True, node_size=600, font_size=8)

def add_GO_function(cluster, go_map, go_db = None):
    """
    Keep track of how many proteins in the cluster have a given GO term
    """
    go_terms = {}
    for protein in cluster["members"]:
        if protein in go_map:
            for gid in set(go_map[protein]):
                if go_db is not None:
                    if gid not in go_db:
                        continue
                go_terms[gid] = go_terms.setdefault(gid, 0) + 1
    return go_terms

def triangles(graph: nx.Graph):
    return int(sum([i for i in nx.triangles(graph).values()]) / 3)

def get_top_terms(cluster, N=5, go_map=None):
    if (go_map is not None) and ("GO_terms" not in cluster):
        cluster["GO_terms"] = add_GO_function(cluster, go_map)
    term_dict = cluster["GO_terms"]
    if N == -1:
        N = len(term_dict)
    return sorted(term_dict.items(), key=lambda x: x[1], reverse=True)[:N] 

def print_cluster(clust, go_database, n_terms=5, return_str=False):
    description_string = ""

    if "llm_name" in clust:
        description_string += f"Cluster Name: {clust.llm_name}\n"

    members = clust["members"]
    short_mem_string = ", ".join(members[:3])
    description_string += (
        f"Cluster of {len(members)} proteins [{short_mem_string}, ...] (hash {hash_cluster(members)})\n"
    )

    if "recipe" in clust:
        recipe_dict = clust["recipe"]
        recipe_metrics = list(recipe_dict.keys())
        for rm in recipe_metrics:
            for deg in recipe_dict[rm].keys():
                nadded = len(recipe_dict[rm][deg])
                description_string += f"{nadded} proteins re-added by ReCIPE ({rm}, {deg})\n"

    if "graph" in clust:
        G = nx.Graph()
        for edge in clust["graph"]:
            G.add_edge(edge[0], edge[1], weight=edge[2])

        description_string += f"Edges: {len(G.edges())}\n"
        description_string += f"Triangles: {triangles(G)}\n"
        description_string += f"Max Degree: {0 if not len(G.edges()) else max(G.degree(), key=lambda x: x[1])[1]}\n"

    if "GO_terms" in clust:
        top_terms = get_top_terms(clust, n_terms)
        description_string += "Top Terms:\n"
        for gid, freq in top_terms:
            try:
                go_name = go_database[gid]
            except KeyError:
                go_name = "Unknown"
            description_string += f"\t\t{gid} - <{go_name}> ({freq})\n"

    if "llm_explanation" in clust:
        llm_desc = clust["llm_explanation"]
        description_string += f"LLM Explanation: {llm_desc}\n"
        # description_string += f"LLM Confidence: {llm_desc}\n"

    if return_str:
        return description_string
    else:
        print(description_string)


def hash_cluster(protein_list):
    return int(hashlib.md5("".join(sorted(protein_list)).encode()).hexdigest(), 16) % (
        2**61 - 1
    )


def nx_graph_cluster(cluster):
    G = nx.Graph()
    for edge in cluster["graph"]:
        G.add_edge(edge[0], edge[1], weight=edge[2])
    return G


def load_cluster_json(infile):
    with open(infile, "r") as f:
        clusters = json.load(f)
    return clusters


def parse_GO_database(infile):
    terms = {}
    with open(infile, "r") as f:
        for line in f:
            line = line.strip()
            if line == "[Term]":
                term_ids = []
                line = f.readline().strip()
                while not line == "":
                    if line.startswith("id:") or line.startswith("alt_id:"):
                        term_ids.append(line.split(": ")[1])
                    elif line.startswith("name:"):
                        main_name = line.split(": ")[1]
                    line = f.readline().strip()
                for tid in term_ids:
                    terms[tid] = main_name
    return terms


def parse_GO_map(f):
    seqDb = pd.read_csv(f, sep=",")
    seqDb.columns = ["seq", "manual_annot", "pfam_list", "GO_list"]
    seqDb["GO_str"] = seqDb["GO_list"]
    seqDb["GO_list"] = seqDb["GO_str"].str.split(";")

    def extract_GO_id_from_list(go_list):
        if isinstance(go_list, list):
            return [i.split("|")[0] for i in go_list]
        else:
            return None

    seqDb["GO_ids"] = seqDb["GO_list"].apply(extract_GO_id_from_list)
    seq2GO = seqDb[["seq", "GO_ids"]]
    seq2GO.columns = ["seq", "GO_ids"]
    go_map = {}
    for _, r in seq2GO.iterrows():
        if r.GO_ids is not None:
            go_map[r.seq] = r.GO_ids
    return go_map

def clean_top_terms(c, go_db, return_counts=True, n_filter=3):
    csize = len(c)
    tt = get_top_terms(c, 1)
    if len(tt):
        if tt[0][1] < n_filter:
            pass
        elif return_counts:
            return (go_db[tt[0][0]], tt[0][1], csize)
        else:
            return go_db[tt[0][0]]
    if return_counts:
        return ('', None, csize)
    else:
        return ''


def plot_cluster(cluster, full_graph, name="Graph", node_labels=False, savefig=None):
    # From https://networkx.org/documentation/stable/auto_examples/drawing/plot_degree.html
    G = nx.subgraph(full_graph, cluster["members"])
    degree_sequence = sorted((d for n, d in G.degree()), reverse=True)
    
    fig = plt.figure("Degree of a random graph", figsize=(8, 8))
    # Create a gridspec for adding subplots of different sizes
    axgrid = fig.add_gridspec(5, 4)

    ax0 = fig.add_subplot(axgrid[0:3, :])
    Gcc = G.subgraph(sorted(nx.connected_components(G), key=len, reverse=True)[0])
    pos = nx.spring_layout(Gcc, seed=10396953)
    nx.draw_networkx_nodes(Gcc, pos, ax=ax0, node_size=20)
    nx.draw_networkx_edges(Gcc, pos, ax=ax0, alpha=0.4)
    if node_labels:
        nx.draw_networkx_labels(Gcc, pos, ax=ax0)
    ax0.set_axis_off()

    ax1 = fig.add_subplot(axgrid[3:, :2])
    ax1.plot(degree_sequence, "b-", marker="o")
    ax1.set_title("Degree Rank Plot")
    ax1.set_ylabel("Degree")
    ax1.set_xlabel("Rank")

    ax2 = fig.add_subplot(axgrid[3:, 2:])
    ax2.bar(*np.unique(degree_sequence, return_counts=True))
    ax2.set_title("Degree histogram")
    ax2.set_xlabel("Degree")
    ax2.set_ylabel("# of Nodes")

    fig.tight_layout()
    sns.despine()
    plt.suptitle(f"{name} ({len(G)} nodes / {len(G.edges())} edges)")

    if savefig:
        plt.savefig(savefig,dpi=300,bbox_inches="tight")
    plt.show()


def write_cluster_fasta(cluster_file, sequence_file, directory = ".", prefix = "cluster"):
    cluster_dict = load_cluster_json(cluster_file)
    seq_dict = SeqIO.to_dict(SeqIO.parse(sequence_file, "fasta"))
    
    file_names = []
    for k, clust in cluster_dict.items():
        fname = f"{directory}/{prefix}_{k}.fasta"
        file_names.append(fname)
        with open(fname, "w+") as f:
            for p in clust["members"]:
                f.write(f">{p}\n")
                f.write(f"{seq_dict[k]}\n")
    return None
