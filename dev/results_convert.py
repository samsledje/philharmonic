# results are in a dict. load this
# python results_convert.py -icf pdam_dscript_distances.DSD1.clusters.csv -rrf ../ReCIPE/results/_3-100_top_100_linear_0.25_qualifying_proteins.json --sep , --outfile pdam_results.csv
import argparse
import json
import collections

def main(args):
    initial_clusters_file = args.init_clusters_file
    recipe_results_file = args.recipe_results_file
    outfile = args.outfile

    qualifying_proteins_by_metric = {}

    with open(initial_clusters_file, "r") as f:
        initial_clusters = f.read().split("\n")

    with open(recipe_results_file, "rb") as f:
        qualifying_proteins_by_metric = json.load(f)

    # print(f"qualifying: {qualifying_proteins_by_metric.keys()}")
    # write output to clusters
    # for metric in qualifying_proteins_by_metric.keys():
    #     with open(metric + "_" + outfile, "w") as f:
        
    metric = args.metric
    percent_conec = args.percent_conec


    # IGNORE PROTEINS ADDED MORE THAN 15 TIMES
    prot_counts = collections.defaultdict(int)
    for k in qualifying_proteins_by_metric[metric][percent_conec].keys():
        for prot in qualifying_proteins_by_metric[metric][percent_conec][k]:
            prot_counts[prot] += 1

    readdition_threshold = args.readdition_threshold
    addition_proteins = set([k for k in prot_counts.keys() if prot_counts[k] <= readdition_threshold])

    for (cluster, prots) in qualifying_proteins_by_metric[metric][percent_conec].items():
        if len(prots) > 0:
            qualifying_proteins_by_metric[metric][percent_conec][cluster] = [prot for prot in prots if prot in addition_proteins]



    print(f"writing to {outfile}")
    with open(outfile, "w") as f:
        for i in range(0, len(initial_clusters)):
            if (len(initial_clusters[i]) == 0):
                continue

            f.write(initial_clusters[i])
            
            if str(i) in qualifying_proteins_by_metric[metric][percent_conec].keys():
                print(f"{len(qualifying_proteins_by_metric[metric][percent_conec][str(i)])} proteins added to cluster {i} for metric {metric} and percent_conec {percent_conec}")
                f.write(args.sep)
                f.write(",".join(qualifying_proteins_by_metric[metric][percent_conec][str(i)]))
            f.write("\n")

        

def get_args():
    """
    """
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--init_clusters_file", "-icf",
        required=True,
        help = "path and name for the initial clusters file, before modification by recipe", 
        type = str, 
    )
    parser.add_argument(
        "--recipe_results_file", "-rrf",
        required=True,
        help = "path and name for the recipe results file", 
        type = str, 
    )
    
    parser.add_argument( # TODO: currently does not allow for multiple metrics
        "--outfile", 
        required=True,
        help = "name for the output file of clusters", 
        type = str, 
    )
    # parser.add_argument( # TODO: which metric to use, which percent_conec to use
    parser.add_argument(
        "--sep",
        required=False,
        help = "separator for the input file (tab in between old and new proteins)",
        type = str,
        default = "\t"
    )

    parser.add_argument(
        "--metric",
        required=False,
        help = "metric to use for qualifying proteins",
        type = str,
        default = "degree"
    )
    parser.add_argument(
        "--percent_conec",
        required=False,
        help = "percent_conec to use for qualifying proteins",
        type = str,
        default = "0.75"
    )

    parser.add_argument(
        "--readdition-threshold",
        required=False,
        help="if a protein is added back more than this many times, ignore it",
        type=int,
        default=15
    )
  
    return parser.parse_args()


if __name__ == "__main__":
    main(get_args())
