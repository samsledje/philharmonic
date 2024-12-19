from itertools import combinations
from pathlib import Path

import os
import argparse
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from loguru import logger
from scipy.spatial.distance import jaccard
from scipy.stats import ttest_ind
from tqdm import tqdm

from philharmonic.utils import (
    add_GO_function,
    load_cluster_json,
    parse_GO_database,
    parse_GO_map,
)

parser = argparse.ArgumentParser("Functional permutation test")
parser.add_argument(
    "--results_dir",
    type=str,
    help="Directory containing results",
    default="results",
)
parser.add_argument(
    "--run_name",
    type=str,
    help="Name of the run",
    default="run",
)
parser.add_argument(
    "--use_go_slim",
    action="store_true",
    help="Use GO slim",
)

if __name__ == "__main__":

    args = parser.parse_args()

    results_dir = Path(args.results_dir)
    run_name = args.run_name
    img_dir = results_dir / "img"
    os.makedirs(img_dir, exist_ok=True)

    GO_OBO_PATH = results_dir / "go.obo"
    GO_SLIM_PATH = results_dir / "goslim_generic.obo"
    GO_MAP_PATH = results_dir / f"{run_name}_GO_map.csv"

    CLUSTER_FILE_PATH = results_dir / f"{run_name}_clusters.json"
    IMG_DIR = results_dir / "img"

    clusters = load_cluster_json(CLUSTER_FILE_PATH)
    go_map = parse_GO_map(GO_MAP_PATH)

    if args.use_go_slim:
        go_db = parse_GO_database(GO_SLIM_PATH)
    else:
        go_db = parse_GO_database(GO_OBO_PATH)

    for clust in clusters.values():
        clust["GO_terms"] = add_GO_function(clust, go_map, go_db=go_db)
        for gt in clust["GO_terms"].keys():
            assert gt in go_db.keys()

    go_assigned = set()
    for clust in tqdm(clusters.values(), desc="Collecting GO terms"):
        for m in clust["members"]:
            go_assigned.update(go_map.get(m, []))
    go_assigned = sorted(list(go_assigned.intersection(go_db.keys())))

    logger.info(f"{len(go_assigned)} GO terms assigned")

    proteins_in_clusters = []
    for clust in tqdm(clusters.values(), desc="Collecting proteins"):
        proteins_in_clusters.extend(clust["members"])
        proteins_in_clusters.extend(list(clust["recipe"]["degree"]["0.75"]))
    logger.info(f"{len(proteins_in_clusters)} proteins in clusters")

    def protein_GO_bit_vector(
        protein_id, go_map, full_go_list, id_col="seq", go_col="GO_ids"
    ):
        go_bv = np.zeros(len(full_go_list))
        prot_go = go_map.get(protein_id)
        if prot_go is not None:
            for gid in prot_go:
                if gid in full_go_list:
                    go_bv[full_go_list.index(gid)] = 1
        return go_bv


    # Compute GO bit vectors for each protein
    protein_GO_bvs = {}
    for pid in tqdm(proteins_in_clusters, desc="Computing bit vectors"):
        protein_GO_bvs[pid] = protein_GO_bit_vector(pid, go_map, go_assigned)

    cluster_jaccards = {}

    for k, clust in tqdm(clusters.items(), desc="Computing Jaccard Similarity"):
        cjaccard = []
        for p1, p2 in combinations(
            clust["members"] + list(clust["recipe"]["degree"]["0.75"]), 2
        ):
            jc = 1 - jaccard(protein_GO_bvs[p1], protein_GO_bvs[p2])
            cjaccard.append(jc)
        cluster_jaccards[k] = np.array(cjaccard)

    rng = np.random.default_rng(seed=42)
    shuffled_bit_vectors = {
        k: v
        for k, v in zip(
            protein_GO_bvs.keys(),
            rng.permutation(list(protein_GO_bvs.values())),
            strict=False,
        )
    }

    cluster_jaccards_perm = {}

    for k, clust in tqdm(clusters.items(), desc="Computing permuted Jaccard Similarity"):
        cjaccard = []
        for p1, p2 in combinations(
            clust["members"] + list(clust["recipe"]["degree"]["0.75"]), 2
        ):
            jc = 1 - jaccard(shuffled_bit_vectors[p1], shuffled_bit_vectors[p2])
            cjaccard.append(jc)
        cluster_jaccards_perm[k] = np.array(cjaccard)

    phil_mean = [np.mean(i) for i in cluster_jaccards.values()]
    permute_mean = [np.mean(i) for i in cluster_jaccards_perm.values()]
    coherence_df = (
        pd.DataFrame(
            {
                "cluster": list(cluster_jaccards.keys()),
                "PHILHARMONIC": phil_mean,
                "Random Clustering": permute_mean,
            }
        )
        .melt("cluster")
        .rename(
            {"variable": "Clustering Method", "value": "Mean Jaccard Similarity"}, axis=1
        )
    )
    coherence_df.head()

    from scipy import interpolate
    def find_histogram_intersection(data1, data2, bins=50, return_curves=False):
        """
        Find the intersection point of two histograms.

        Parameters:
        data1: array-like, first dataset
        data2: array-like, second dataset
        bins: int, number of bins for histogram
        return_curves: bool, whether to return the interpolated curves

        Returns:
        float: x-value where histograms intersect
        """
        # Create histograms
        hist1, bin_edges = np.histogram(data1, bins=bins, density=True)
        hist2, _ = np.histogram(data2, bins=bin_edges, density=True)

        # Get bin centers
        bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2

        # Create interpolation functions
        f1 = interpolate.interp1d(
            bin_centers, hist1, kind="linear", fill_value="extrapolate"
        )
        f2 = interpolate.interp1d(
            bin_centers, hist2, kind="linear", fill_value="extrapolate"
        )

        # Find intersection points
        x_range = np.linspace(bin_centers[0], bin_centers[-1], 1000)
        y1 = f1(x_range)
        y2 = f2(x_range)

        # Find where the difference changes sign
        diff = y1 - y2
        intersection_indices = np.where(np.diff(np.signbit(diff)))[0]

        if len(intersection_indices) == 0:
            raise ValueError("No intersection found")

        # Get x values at intersections
        intersection_points = x_range[intersection_indices]

        if return_curves:
            return max(intersection_points), f1, f2, x_range

        return max(intersection_points)

    sns.set_palette("colorblind")
    sns.set_theme(style="white", palette="pastel", font_scale=1)

    # Create figure and gridspec
    fig = plt.figure(figsize=(10, 8))
    gs = fig.add_gridspec(2, 1, height_ratios=[3, 2], hspace=0.05)

    # Create top subplot for histogram
    ax0 = fig.add_subplot(gs[0])
    ax1 = fig.add_subplot(gs[1])

    sns.histplot(
        data=coherence_df,
        x="Mean Jaccard Similarity",
        hue="Clustering Method",
        alpha=0.3,
        bins=np.arange(0, 1.05, 0.01),
        kde=True,
        palette=["blue", "red"],
        common_norm=False,
        ec="white",
        ax=ax0,
    )

    sns.boxplot(
        data=coherence_df,
        x="Mean Jaccard Similarity",
        hue="Clustering Method",
        palette=["blue", "red"],
        ax=ax1,
    )

    ax0.set_xlabel("")  # Remove x-label from top plot
    ax0.set_xticklabels([])  # Remove x-ticks from top plot

    ax0.set_xlim(ax1.get_xlim())  # Align the x-axis of both subplots
    ax1.get_legend().remove()  # Remove legend from bottom plot

    tstat, p = ttest_ind(phil_mean, permute_mean, alternative="greater")
    go_type = "GO Slim" if args.use_go_slim else "All GO Terms"
    ax0.set_title(f"Functional Enrichment: PHILHARMONIC, {go_type} (p={p:.3})")

    intersection_x = find_histogram_intersection(phil_mean, permute_mean, bins=50)
    logger.info(f"Distributions cross at {intersection_x:.3f}")
    ax0.axvline(intersection_x, linestyle="--", color="black")

    # Show the plot
    sns.despine()
    finame = (
        f"{run_name}_function_enrichment_GOfull.png"
        if not args.use_go_slim
        else f"{run_name}_function_enrichment_GOslim.png"
    )
    plt.savefig(IMG_DIR / finame, bbox_inches="tight", dpi=300)
    # plt.show()

    with open(results_dir / f"{run_name}_stats.txt", "w+") as f:
        f.write(f"# statistics: {p:.3e} {tstat:.3f}\n")
        for _, r in coherence_df.iterrows():
            f.write(f"{r['cluster']}\t{r['Clustering Method']}\t{r['Mean Jaccard Similarity']}\n")
