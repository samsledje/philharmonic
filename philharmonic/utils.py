import hashlib
import json
import typing as T
from collections import defaultdict
from pathlib import Path

import networkx as nx
import numpy as np
import pandas as pd
import seaborn as sns
from Bio import SeqIO
from loguru import logger
from matplotlib import pyplot as plt
from matplotlib.colors import LinearSegmentedColormap


def add_GO_function(
    cluster: T.Dict, go_map: T.Dict, go_db: T.Optional[T.Dict] = None
) -> T.Dict[str, int]:
    """
    Keep track of how many proteins in the cluster have a given GO term
    """
    go_terms: dict[str, int] = dict()
    for protein in cluster["members"]:
        if protein in go_map:
            for gid in set(go_map[protein]):
                if go_db is not None:
                    if gid not in go_db:
                        continue
                go_terms[gid] = go_terms.setdefault(gid, 0) + 1
    return go_terms


def calculate_graph_triangles(graph: nx.Graph) -> int:
    return int(sum([i for i in nx.triangles(graph).values()]) / 3)


def get_cluster_top_terms(
    cluster: T.Dict, N: int = 10, go_map: T.Optional[T.Dict] = None
) -> T.List[T.Tuple[str, int]]:
    if (go_map is not None) and ("GO_terms" not in cluster):
        cluster["GO_terms"] = add_GO_function(cluster, go_map)
    term_dict = cluster["GO_terms"]
    if N == -1:
        N = len(term_dict)
    return sorted(term_dict.items(), key=lambda x: x[1], reverse=True)[:N]


def print_cluster(
    cluster: T.Dict,
    go_database: T.Dict,
    n_terms: int = 10,
    return_str: bool = False,
) -> T.Optional[str]:
    description_string = ""

    if "llm_name" in cluster:
        description_string += f"Cluster Name: {cluster['llm_name']}\n"

    members = cluster["members"]
    short_mem_string = ", ".join(members[:3])
    description_string += f"Cluster of {len(members)} proteins [{short_mem_string}, ...] (hash {hash_cluster(members)})\n"

    if "recipe" in cluster:
        recipe_dict = cluster["recipe"]
        recipe_metrics = list(recipe_dict.keys())
        for rm in recipe_metrics:
            for deg in recipe_dict[rm].keys():
                nadded = len(recipe_dict[rm][deg])
                description_string += (
                    f"{nadded} proteins re-added by ReCIPE ({rm}, {deg})\n"
                )

    if "graph" in cluster:
        G = nx.Graph()
        for edge in cluster["graph"]:
            G.add_edge(edge[0], edge[1], weight=edge[2])

        description_string += f"Edges: {len(G.edges())}\n"
        description_string += f"Triangles: {calculate_graph_triangles(G)}\n"
        description_string += f"Max Degree: {0 if not len(G.edges()) else max(G.degree(), key=lambda x: x[1])[1]}\n"

    if "GO_terms" in cluster:
        top_terms = get_cluster_top_terms(cluster, n_terms)
        description_string += "Top Terms:\n"
        for gid, freq in top_terms:
            try:
                go_name = go_database[gid]
            except KeyError:
                go_name = "Unknown"
            description_string += f"\t\t{gid} - <{go_name}> ({freq})\n"

    if "llm_explanation" in cluster:
        llm_desc = cluster["llm_explanation"]
        llm_confidence = cluster["llm_confidence"]
        description_string += f"LLM Explanation: {llm_desc}\n"
        description_string += f"LLM Confidence: {llm_confidence}\n"

    if return_str:
        return description_string
    else:
        print(description_string)
        return None


def hash_cluster(protein_list: T.List[str]) -> int:
    return int(hashlib.md5("".join(sorted(protein_list)).encode()).hexdigest(), 16) % (
        2**61 - 1
    )


def nx_graph_cluster(
    cluster: T.Dict,
    full_G: T.Optional[nx.Graph] = None,
    use_recipe_nodes: bool = False,
    recipe_metric: str = "degree",
    recipe_cthresh: str = "0.75",
) -> nx.Graph:
    if use_recipe_nodes:
        if full_G is None:
            logger.error("No full graph provided")
            raise ValueError("No full graph provided")
        recipe_prots = cluster["recipe"][recipe_metric][recipe_cthresh]
        if not isinstance(recipe_prots, list):
            recipe_prots = list(recipe_prots)
        base_prots = cluster["members"]
        clustG = full_G.subgraph(base_prots + recipe_prots)
    else:
        clustG = nx.Graph()
        for edge in cluster["graph"]:
            clustG.add_edge(edge[0], edge[1], weight=edge[2])

    return clustG


def load_cluster_json(infile: T.Union[str, Path]) -> T.Dict:
    with open(infile, "r") as f:
        clusters = json.load(f)
    return clusters


def parse_GO_graph(go_graph_file: T.Union[str, Path]) -> T.Tuple[T.Dict, T.Dict]:
    go2children = defaultdict(list)
    go2desc = dict()

    def process_block(b):
        if not ("id" in b and "name" in b and "namespace" in b):
            return

        go2desc[b["id"]] = (b["namespace"], b["name"])
        if "parent" in b:
            for c in b["parent"]:
                go2children[c].append(b["id"])

    block: T.Dict[str, T.Any] = dict()
    for line in open(go_graph_file):
        if line.startswith("[Term]"):
            if block:
                process_block(block)
                block = {}
        if line.startswith("id: "):
            block["id"] = line.split()[1]

        if line.startswith("is_a: "):
            if "parent" not in block:
                block["parent"] = []
            block["parent"].append(line.split()[1])

        if line.startswith("name:"):
            block["name"] = line[6:].strip()

        if line.startswith("namespace:"):
            block["namespace"] = line[11:].strip()

    if block:
        process_block(block)

    return go2children, go2desc


def subset_GO_graph(
    go_graph_file: T.Union[str, Path], go_included_terms: T.List[str]
) -> T.List:
    go2children, go2desc = parse_GO_graph(go_graph_file)

    visited_go_terms = set()

    def dfs_visit(go2children, goid):
        for c in go2children[goid]:
            dfs_visit(go2children, c)
        visited_go_terms.add(goid)

    for c1 in go_included_terms:
        dfs_visit(go2children, c1)

    go_subset = []
    for c in visited_go_terms:
        go_subset.append((c, go2desc[c][0], go2desc[c][1], ";".join(go2children[c])))

    return sorted(go_subset)


def parse_GO_database(infile: T.Union[str, Path]) -> T.Dict:
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


def parse_GO_map(file_path: T.Union[str, Path]) -> T.Dict[str, T.List[str]]:
    seqDb = pd.read_csv(file_path, sep=",")
    seqDb.columns = pd.Index(["seq", "manual_annot", "pfam_list", "GO_list"])
    seqDb["GO_str"] = seqDb["GO_list"]
    seqDb["GO_list"] = seqDb["GO_str"].str.split(";")

    def extract_GO_id_from_list(go_list):
        if isinstance(go_list, list):
            return [i.split("|")[0] for i in go_list]
        else:
            return None

    seqDb["GO_ids"] = seqDb["GO_list"].apply(extract_GO_id_from_list)
    seq2GO = seqDb[["seq", "GO_ids"]]
    seq2GO.columns = pd.Index(["seq", "GO_ids"])
    go_map: T.Dict[str, T.List[str]] = dict()
    for _, r in seq2GO.iterrows():
        if r.GO_ids is not None:
            go_map[r.seq] = r.GO_ids
    return go_map


def filter_proteins_GO(
    proteins: T.List[str],
    go_filter_f: T.Union[str, Path],
    go_map_f: T.Union[str, Path],
    go_database_f: T.Union[str, Path],
) -> T.Set[str]:
    with open(go_filter_f, "r") as f:
        # Get children of allowed GO terms
        allowed_go_initial = [line.strip() for line in f]
        allowed_go = set(
            [i[0] for i in subset_GO_graph(go_database_f, allowed_go_initial)]
        )

        # Filter proteins by GO terms
        go_map = parse_GO_map(go_map_f)
        allowed_proteins = []
        for protein, go_terms in go_map.items():
            if go_terms is not None and any(gt in allowed_go for gt in go_terms):
                allowed_proteins.append(protein)

    return set(allowed_proteins).intersection(proteins)


def clean_top_terms(
    clust: T.Dict,
    go_db: T.Dict[str, str],
    return_counts: bool = False,
    n_filter: int = 3,
):
    csize = len(clust["members"])
    tt = get_cluster_top_terms(clust, 1)
    if len(tt):
        if tt[0][1] < n_filter:
            pass
        elif return_counts:
            return (go_db[tt[0][0]], tt[0][1], csize)
        else:
            return go_db[tt[0][0]]
    if return_counts:
        return ("", None, csize)
    else:
        return ""


def get_node_colors(
    cluster: T.Dict,
    recipe_metric: str = "degree",
    recipe_cthresh: str = "0.75",
    base_color: str = "blue",
    recipe_color: str = "red",
) -> T.Dict[str, str]:
    colors: T.Dict[str, str] = dict()
    for k in cluster["members"]:
        colors[k] = base_color
    for k in cluster["recipe"][recipe_metric][recipe_cthresh]:
        colors[k] = recipe_color
    return colors


def plot_degree(
    G: nx.Graph,
    name: str = "Graph",
    node_colors: T.Optional[T.Dict[str, str]] = None,
    node_labels: T.Optional[T.List[str]] = None,
    savefig: T.Optional[T.List[str]] = None,
) -> None:
    # From https://networkx.org/documentation/stable/auto_examples/drawing/plot_degree.html
    degree_sequence = sorted((d for n, d in G.degree()), reverse=True)
    fig = plt.figure(name, figsize=(8, 8))
    # Create a gridspec for adding subplots of different sizes
    axgrid = fig.add_gridspec(5, 4)

    ax0 = fig.add_subplot(axgrid[0:3, :])
    # Gcc = G.subgraph(sorted(nx.connected_components(G), key=len, reverse=True)[0])
    Gcc = G
    pos = nx.spring_layout(Gcc, seed=10396953)
    if node_colors is not None:
        nx.draw_networkx_nodes(
            Gcc,
            pos,
            ax=ax0,
            node_size=20,
            node_color=[node_colors[n] for n in G.nodes()],
        )
    else:
        nx.draw_networkx_nodes(Gcc, pos, ax=ax0, node_size=20)
    nx.draw_networkx_edges(Gcc, pos, ax=ax0, alpha=0.4)
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
        plt.savefig(savefig, dpi=300, bbox_inches="tight")
    plt.show()


def plot_cluster(
    cluster: T.Dict,
    full_graph: nx.Graph,
    name: str = "Graph",
    node_labels: T.Optional[T.List[str]] = None,
    use_recipe: bool = True,
    recipe_metric: str = "degree",
    recipe_cthresh: str = "0.75",
    savefig: T.Optional[T.List[str]] = None,
) -> None:
    # From https://networkx.org/documentation/stable/auto_examples/drawing/plot_degree.html

    G = nx_graph_cluster(cluster, full_G=full_graph, use_recipe_nodes=use_recipe)
    node_colors = get_node_colors(
        cluster, recipe_metric=recipe_metric, recipe_cthresh=recipe_cthresh
    )
    plot_degree(
        G, name=name, node_colors=node_colors, node_labels=node_labels, savefig=savefig
    )


def write_cluster_fasta(
    cluster_file: T.Union[str, Path],
    sequence_file: T.Union[str, Path],
    directory: T.Union[str, Path] = ".",
    prefix: str = "cluster",
) -> None:
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


def write_cluster_cytoscape(
    cluster: T.Dict,
    full_G: nx.Graph,
    outfile: T.Union[str, Path] = Path("cytoscape_input.txt"),
    with_recipe: bool = True,
    recipe_metric: str = "degree",
    recipe_cthresh: str = "0.75",
) -> None:
    if not isinstance(outfile, Path):
        outfile = Path(outfile)

    G = nx_graph_cluster(
        cluster,
        full_G=full_G,
        use_recipe_nodes=with_recipe,
        recipe_metric=recipe_metric,
        recipe_cthresh=recipe_cthresh,
    )

    with open(outfile, "w") as f:
        f.write("Source\tTarget\tWeight\n")
        for edge in G.edges(data=True):
            f.write(f"{edge[0]}\t{edge[1]}\t{edge[2]['weight']}\n")

    with open(outfile.with_suffix(".nodes.txt"), "w") as f:
        f.write("Node\tType\n")
        for node in G.nodes():
            if node in cluster["recipe"][recipe_metric][recipe_cthresh]:
                f.write(f"{node}\tRecipe\n")
            else:
                f.write(f"{node}\tMember\n")

    return None


def create_rainbow_colorbar(
    vmin: T.Union[int, float] = 0,
    vmax: T.Union[int, float] = 100,
    size: T.Tuple[int, int] = (1, 1),
    step: int = 25,
    label: str = "pLDDT",
    savefig: T.Optional[T.Union[str, Path]] = None,
) -> plt.Figure:
    # Create figure and axes
    fig, ax = plt.subplots(figsize=size)

    # Create rainbow color map
    colors = ["red", "orange", "yellow", "green", "blue"]
    n_bins = 100
    cmap = LinearSegmentedColormap.from_list("rainbow", colors, N=n_bins)

    # Create color bar
    norm = plt.Normalize(vmin, vmax)
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])

    cbar = plt.colorbar(sm, cax=ax, orientation="vertical")
    cbar.set_ticks([float(i) for i in np.arange(vmin, vmax + step, step)])
    cbar.set_ticklabels([str(i) for i in np.arange(vmin, vmax + step, step)])

    # Add label to the left side
    ax.yaxis.set_label_position("left")
    ax.yaxis.set_ticks_position("right")
    ax.set_ylabel(label, rotation=90, va="bottom", labelpad=5)

    # Remove frame
    ax.set_frame_on(False)

    # Save figure
    if savefig:
        plt.savefig(savefig, bbox_inches="tight", dpi=300)

    return fig
