#!/Users/charlottev/micromamba/envs/Research/bin/python


# This script is used to add proteins to clusters to reconnect the clusters, based on the structure of the network.
# run with python recipe.py --network ../to-run-recipe-on/20240202_apall_hits_v_all_dscript_out.positive.tsv --cluster-filepath ../to-run-recipe-on/protein_to_cluster.json --lr .1 --max 100 -cthresh 0.75
# the result is a .json file with the added proteins for each cluster, for each metric and connectivity threshold
# the results are structured in a dict with the metric as the key, and within the metric, keys for each connectivity threshold, and within each connectivity threshold, keys for each cluster name. The corresponding values are the added proteins for each cluster.

import time
import json
import sys
import os
sys.path.append(os.getcwd())
import argparse

from math import sqrt
from math import floor
from math import ceil
from collections import defaultdict


import numpy as np
import pandas as pd

import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)
warnings.simplefilter(action='ignore', category=RuntimeWarning)


# importing the classes recipe uses
from matrix_class import ProteinMatrix # ppi matrix 
from cluster_class import AllClusters # dictionary to hold all clusters (in form number of cluster : list of proteins in that cluster)
from degreelist_class import DegreeList # creates a list of all proteins in order of their degree
from matrix_class import SubMatrix

# helper functions for setting up program
from recipe_utils import initialize_matrix_clusters_degreelist

def compute_qualifying_proteins(
    matrix,
    degreelist,
    clusters,
    metric = 'degree',
    sort_ascending = False,
    connectivity_threshold = 1,
    degree_cutoff = 5000,
    linear_ratio = None,
    max_added_proteins = None,
    lower_bound = 3,
    upper_bound = 8
):
    proteins = matrix.get_list_of_proteins()
    degree_dict = dict(degreelist.sorted_protein_degree_dict)
    matrix_df = matrix.get_matrix()
    all_proteins_to_add = {}

    for cluster_num in clusters.filter_clusters_by_size(lower_bound, upper_bound).keys():
        # initialise list of proteins to add
        added_proteins = []
        protein_to_add = "initialise"
        # get all the proteins associated to a cluster
        cluster_proteins = clusters.get_cluster_proteins(cluster_num)
        # get the list of potential proteins to add to cluster 
        potential_proteins = list(filter(lambda prot: prot not in cluster_proteins and degree_dict[prot] < degree_cutoff, proteins))
        # TODO: is there a better way of doing this??
        submatrix = SubMatrix(cluster_proteins, matrix)
        components_and_labels = submatrix.get_num_components_and_labels()
        num_components = components_and_labels[0]
        # current ratio of clusters to proteins
        num_proteins = len(cluster_proteins)
        percent_connectivity = 1 - (num_components - 1) / (num_proteins - 1)
        # loop through all the proteins and add proteins based on score
        while protein_to_add and percent_connectivity < connectivity_threshold:
            qualifying_proteins = {}
            connection_sitch = None
            # get sqrt of number of components in the subgraph
            if linear_ratio:
                connection_sitch = floor(linear_ratio * len(np.unique(components_and_labels[1])))
            else:
                connection_sitch = floor(sqrt(len(np.unique(components_and_labels[1]))))
            min_connection = connection_sitch if connection_sitch > 1 else 2

            for protein in potential_proteins:
                a = protein not in matrix_df
                if a:
                    with open("mismatched_proteins.txt", 'a+') as f:
                        f.write(f"{protein}\n")
                else:
                    protein_degree = degree_dict[protein]
                    if protein_degree >= min_connection:
                        # create component dictionary
                        protein_component_dictionary = dict(zip(submatrix.get_matrix().index, components_and_labels[1]))
                        # swap the values so the component number is the key
                        component_dictionary = defaultdict(list) 
                        for key, val in protein_component_dictionary.items():
                            component_dictionary[val].append(key)
                        # get number of connected components
                        num_components_protein_connects = 0
                        for component_number in range(num_components):
                            if next((prot for prot in component_dictionary[component_number] if prot in matrix_df and matrix_df[prot][protein]), None):
                                num_components_protein_connects = num_components_protein_connects + 1
                        # if connection, greater than cutoff, consider for re-addition
                        if num_components_protein_connects >= min_connection:
                            qualifying_proteins[protein] = {
                                'components_connected': num_components_protein_connects,
                                'degree': protein_degree,
                                'score': num_components_protein_connects * (1 / protein_degree)
                            }
            if max_added_proteins:
                sorted_qualifying_proteins = sorted(qualifying_proteins.items(), key = lambda x: x[1][metric], reverse=not sort_ascending)
                for _ in range(max_added_proteins):
                    if len(sorted_qualifying_proteins):
                        added_proteins.append(sorted_qualifying_proteins.pop()[0])
                    else:
                        break
                protein_to_add = None
                percent_connectivity = 1
            else:
                protein_to_add = sorted(qualifying_proteins.items(), key = lambda x: x[1][metric], reverse=sort_ascending)[0][0] if qualifying_proteins else None
                if protein_to_add:
                    potential_proteins.remove(protein_to_add)
                    added_proteins.append(protein_to_add)
                    percent_connectivity = 1 - (num_components - 1) / (num_proteins - 1)
                    # get number of components in original cluster
                    submatrix = SubMatrix(cluster_proteins + added_proteins, matrix)
                    components_and_labels = submatrix.get_num_components_and_labels()
                    num_components = components_and_labels[0]

        if len(added_proteins):
            all_proteins_to_add[cluster_num] = added_proteins
    return all_proteins_to_add

def main(args):
    max_added_proteins = args.max # 3
    size = f"{args.lb}-{args.ub}"
    lr = args.lr

    matrix, clusters, degreelist = initialize_matrix_clusters_degreelist(args.network_filepath, args.cluster_filepath, csv_clusters_have_labels=args.clusters_labeled)
    qualifying_proteins_by_metric = {} # Dict of structure: {metric: {connectivity_threshold: {cluster_id: [qualifying_proteins]}}}
    
    connectivity_thresholds = []
    if args.connectivity_threshold == -1.0:
        connectivity_thresholds = [0.1, 0.25, 0.5, 0.75, 1.0] 
    else:
        connectivity_thresholds = [args.connectivity_threshold]
    print("connectivity_thresholds", connectivity_thresholds, time.ctime())
    
    metrics = {'degree': False, 'components_connected': True, 'score': True}
    for metric in metrics.items():
        # skip all metrics not specified
        if args.metric != "all" and metric[0] not in args.metric: continue 

        print("starting metric:", metric[0])
        qualifying_proteins_at_threshold = {}

        for connectivity_threshold in connectivity_thresholds:
            print("   -starting threshold", connectivity_threshold, time.ctime())
            adjusted_lb = args.lb if args.lb == 3 else args.lb + 1
            qualifying_proteins = compute_qualifying_proteins(
                matrix,
                degreelist,
                clusters,
                metric = metric[0],
                sort_ascending = metric[1],
                connectivity_threshold = connectivity_threshold,
                degree_cutoff = 5000,
                linear_ratio = lr,
                max_added_proteins = max_added_proteins,
                lower_bound=adjusted_lb,
                upper_bound=args.ub
            )

            if qualifying_proteins:
                qualifying_proteins_at_threshold[connectivity_threshold] = qualifying_proteins
                avg_proteins_added = sum([len(proteins) for proteins in qualifying_proteins.values()]) / len(qualifying_proteins)
                
                # NOTMETRIC = f"lr: {lr}" if lr else "sqrt" # TODO figure out how to do LR and SQRT Qualifiers
                print(f"at threshold {connectivity_threshold}, and {metric}, {len(qualifying_proteins)} clusters have an average of {avg_proteins_added} proteins added")
        if qualifying_proteins_at_threshold:
            qualifying_proteins_by_metric[metric[0]] = qualifying_proteins_at_threshold
            print(f"qualifying_proteins_by_metric[{metric[0]}]" , "TODO: remove")
            
    # output results
    if not args.modify_clusters: # print dict of all clusters
        
        result_file = args.outfile
        with open(result_file, "wb") as f:
            f.write(json.dumps(qualifying_proteins_by_metric).encode('utf-8'))
        print(f"results saved to {result_file}")
    # else: # print modified clusters in the same format as the original clusters
    #     print("HELP!! TODO: if cluster number exists, retain that otherwise use i to map clusters to waht they were added")
    #     result_file = f"{args.outdir}/{args.new_clusters_outfile}"
    #     with open(result_file, "w") as f:
    #         if args.clusters_labeled:

        

def get_args():
    """
    """
    parser = argparse.ArgumentParser()
    # TODO: Charlotte: add a default value for the network filepath
    # TODO: Charlotte: add an option for structure of the network filepath
    # TODO: Charlotte: add a default value for the cluster filepath
    # Location Arguments: (where data is, and where to save it)
    parser.add_argument(
        "-cfp", "--cluster-filepath", 
        required=True,
        help = "Cluster filepath", 
        type = str
    )
    parser.add_argument (
        "-cfl", "--clusters-labeled", 
        required=False,
        help = "If a CSV file of clusters is passed, clusters have labels. Default: False",
        type = bool,
        default = False
    )
    parser.add_argument(
        "-nfp", "--network-filepath", 
        help = "Network filepath", 
        required=True,
        type = str
    )
    parser.add_argument("--outfile", help = "Output file", type = str, required=True)

    # Arguments for output file options
    parser.add_argument (
        "--modify-clusters", 
        help = "Format of the output file. default is false, meaning the dict (which maps added proteins to clusters, and retains all param options) is printed. if set to true, the modified clusters are printed. Default: False", 
        type = bool, 
        required = False,
        default = False
    )
    # parser.add_argument (
    #     "--new-clusters-outfile", 
    #     required=False,
    #     help = "name for the output file of clusters with qualifying proteins added", 
    #     type = str, 
    #     default = "updated_clusters.csv"
    # )

    # Arguments for which clusters to run ReCIPE on
    # Based on number of proteins in cluster
    parser.add_argument( 
        "--lb", 
        required=False,
        help = "Lower bound (inclusive) for cluster size. Default: 3", 
        type = int,
        default = 3,
    )
    parser.add_argument(
        "--ub", 
        required=False,
        help = "Upper bound (exclusive) for cluster size. Default: 100", 
        type = int,
        default=100,
    )
    # Arguments for ReCIPE
    parser.add_argument(
        "--lr", 
        required=True,
        help = "Linear ratio (if not using sqrt). Default = None", 
        type = float,
        default = None,  

    )
    parser.add_argument(
        "--connectivity-threshold", "-cthresh",
        required=False,
        help = 'Connectivity threshold to add proteins until. Default = -1.0 (yields connectivity thresholds [0.1, 0.25, 0.5, 0.75, 1.0]) (if only a single option is desired, 0.75 is recommended)',
        default = -1.0,
        type=float
    )
    parser.add_argument(
        "--metric", "-wm",
        required=False,
        help = "Which metric(s) to use to rank proteins to be added back. Default: all. Options: all, degree, components_connected, score",
        type = str,
        default = "all"
    )

    parser.add_argument(
        "--max", 
        required=False,
        help = "Max number of proteins to add to a cluster. Default = 20", 
        type=int,
        default = 20 # TODO: SHOULD BE NONE TO ALLOW FOR CONNECTIVITY THRESHOLDS
                      # TODO: Do connectivity parameter also 
    )
    # parser.add_argument("--ic", help = "Spectral parameter", type = int)
    return parser.parse_args()

if __name__ == "__main__":
    main(get_args())

