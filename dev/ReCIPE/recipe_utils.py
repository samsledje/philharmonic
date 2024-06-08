"""
Author: Charlotte Versavel
Date:   July 2022
Last touched: Oct 2022

                            recipe_utils.py

"""

from matrix_class import ProteinMatrix
from matrix_class import SubMatrix
from cluster_class import AllClusters
from degreelist_class import DegreeList

import numpy as np
import pandas as pd
import func_e.vocabs.all as vocabs
from func_e.FUNC_E import FUNC_E 

from json import load
from math import sqrt
from math import ceil

"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
 * * * * * * * * * * * * * * * FUNCTIONS * * * * * * * * * * * * * * *
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

# def initialize_clusters(clusters_filepath: str):
#     """
#     TODO
#     a file that has interactions of the form protein1 TAB protein2 TAB interaction
#     a file that containing a dictionary, for each protein in a cluster, the protein's name is linked to its cluster number
#     """
#     clusters_dict = {}
#     # convert actual cluster file to a dictionary!!
#     with open(clusters_filepath, "r") as cluster_dict_file:
#         clusters_dict = load(cluster_dict_file)
#     clusters = AllClusters(protein_to_cluster_dict=clusters_dict)
#     # NOTE:  the above is commented out because the file format has changed from a .json to a csv file

#     return clusters

# def initialize_matrix_degreelist(interactions_filepath: str):
#     """
#     TODO
#     a file that has interactions of the form protein1 TAB protein2 TAB interaction
#     a file that containing a dictionary, for each protein in a cluster, the protein's name is linked to its cluster number
#     """
    
#     matrix = ProteinMatrix(interactions_filepath)
#     # clusters = AllClusters(csv_filename=clusters_filepath)
#     degreelist = DegreeList(matrix)

#     return matrix, degreelist

def initialize_matrix_clusters_degreelist(interactions_filepath: str, clusters_filepath: str, csv_clusters_have_labels: bool = False):
    """
    TODO
    interactions_filepath: a file that has interactions of the form protein1 TAB protein2 TAB interaction
    clusters_filepath: either a csv of clusters (one each line), or a file that containing a dictionary with protein keys and cluster identifier values
    """
    
    clusters: AllClusters = None
    if clusters_filepath.endswith(".csv"):
        clusters = AllClusters(csv_filename=clusters_filepath, csv_clusters_have_labels=csv_clusters_have_labels)
        
    else: # file is a .json dict mapping 
        clusters_dict = {}
        with open(clusters_filepath, "r") as cluster_dict_file:
            clusters_dict = load(cluster_dict_file)
        clusters = AllClusters(protein_to_cluster_dict=clusters_dict)


    matrix = ProteinMatrix(interactions_filepath)
    degreelist = DegreeList(matrix)

    return matrix, clusters, degreelist


def print_protein_background_to_file(matrix: ProteinMatrix, filename: str = "background_proteinlist.txt") -> None:
    """ TODO """
    with open(filename, 'w') as file:
        for protein in matrix.get_list_of_proteins():
            file.write(f"{protein}\n")


def find_clusters_and_proteins_together(matrix: ProteinMatrix, clusters: AllClusters, degreelist: DegreeList, cluster_ratio: float = 1, cluster_constant: int = 0, protein_ratio: float = .05, protein_constant: int = 2, min_components_that_protein_connects: int = 2, max_degree: int = 500, use_sqrt:bool = False, find_clusters_that_are_MORE_connected: bool = False) -> list() and dict():
    """
    function is a version of find_clusters_that_match_criteria, that, once it finds the cluster, finds corresponding proteins at the same time so that the submatrix doesn't need to be reconstructed

    Parameters: 
        matrix - a ProteinMatrix of all protein interactions
        clusters - an AllClusters containing proteins grouped into clusters
        cluster_ratio and cluster_constant - used together to determine which clusters qualify, with the output of the function being cluster_ratio * input + cluster_constant
        when cluster ratio is 1, all clusters will qualify
        TODO: remaining parameters
    Purpose:    determines clusters that are mostly highly connected, then 
                determines which proteins that, when added to the cluster, will 
                increase it's connectedness
    Returns:    a list containing the numbers of the clusters that qualify, and 
                a dictionary linking each cluster, to a list of the qualifying 
                proteins
    TODO: a cluster will qualify if -> IT IS HIGHLY UNCONNECTED
    """
    
    cluster_nums_that_qualify = list()
    qualifying_proteins_dict = dict()

    for cluster_num in clusters.get_all_clusters():
        # create a submatrix out of the proteins in the cluster
        submatrix = SubMatrix(clusters.get_cluster_proteins(cluster_num), matrix)
        num_components, labels = submatrix.get_num_components_and_labels()
        # print(f"num components is {num_components}. num proteins is {len(submatrix.get_list_of_proteins())}")
        if (num_components >= (cluster_ratio * len(submatrix.get_list_of_proteins()) + cluster_constant)) != find_clusters_that_are_MORE_connected:

            # add cluster to list showing that it qualifies, 
            cluster_nums_that_qualify.append(cluster_num)


######## min components has to do with sqrt of the number components

            ### PLEASE NOTE: THIS FUNCTION CALL FINDS PROTEINS BASED ON THE NUMBER OF COMPONENTS IN A CLUSTER ###
            qualifying_proteins = qualifying_proteins_using_num_components(cluster_num, submatrix, clusters, degreelist, ratio=protein_ratio, constant=protein_constant, max_degree=max_degree, use_sqrt=use_sqrt, min_components_that_protein_connects=min_components_that_protein_connects)

########
            if qualifying_proteins: # not empty
                qualifying_proteins_dict[cluster_num] = qualifying_proteins
            
    return cluster_nums_that_qualify, qualifying_proteins_dict


def qualifying_proteins_using_submatrix(cluster_num: int, submatrix: SubMatrix, clusters: AllClusters, degreelist: DegreeList, ratio: float = .5, constant: int = 0, min_components_that_protein_connects: int = 2, max_degree: int = 500, use_sqrt:bool = False) -> list():
    """
    TODO : a revised version of the find_proteins_that_match_criteria function that takes in a submatrix as a parameter, and therefore doesn't need to construct one. 
    TODO: will need to check max degree of the protein somewhere (maybe in the find clusters and proteins together function)
    """
    if use_sqrt:
        if ratio == 0: # can use 0 as a value to avoid thinking about the ratio
            ratio == 1
        min_components_that_protein_connects = max(int(constant + ratio * sqrt(len(clusters.get_cluster_proteins(cluster_num)))), min_components_that_protein_connects)
    else:
        min_components_that_protein_connects = max(constant + ratio * len(clusters.get_cluster_proteins(cluster_num)), min_components_that_protein_connects)
        
    num_components, labels = submatrix.get_num_components_and_labels()

    ### POPULATE COMPONENT DICTIONARY ###
    component_dictionary = dict() # protein : component_num
    j = 0
    for array in [(np.array(submatrix.get_list_of_proteins())[np.nonzero(labels == i)]) for i in range(num_components)]:
        for protein in array:
            component_dictionary[protein] = j
        j += 1
    
    ## FIND CONNECTED PROTEINS AND DETERMINE IF THEY QUALIFY 
    qualifying_proteins = list()

    for protein in (degreelist.get_list_of_proteins_sorted_by_degree()):   

        if (degreelist.get_degree_of_protein(protein, max_degree=max_degree) <= max_degree):
            num_edges, which_proteins = degreelist.determine_num_edges_to_cluster(protein, clusters.get_cluster_proteins(cluster_num), also_return_which_proteins=True)
                    
            if (num_edges >= min_components_that_protein_connects):
                set_of_components_that_protein_connects = degreelist.which_components_of_a_cluster_would_a_protein_connect(protein, clusters.get_cluster_proteins(cluster_num), component_dictionary, connected_proteins_within_cluster=which_proteins)

                if len(set_of_components_that_protein_connects) >= min_components_that_protein_connects:
                    qualifying_proteins.append(protein)

    return qualifying_proteins


#############################################################################




#############################################################################


def qualifying_proteins_using_num_components(cluster_num: int, submatrix: SubMatrix, clusters: AllClusters, degreelist: DegreeList, ratio: float = .5, constant: int = 0, max_degree: int = 500, use_sqrt:bool = False,     min_components_that_protein_connects = 2) -> list():
    """
    TODO : a revised version of the qualifying_proteins_using_submatrix fxn that incorperates the number of components in a cluster when determining if a protein qualifies to reconnect it.


    """
    num_components, labels = submatrix.get_num_components_and_labels()


    if use_sqrt:
        min_components_that_protein_connects = max(min_components_that_protein_connects, ceil(constant + ratio * sqrt(num_components)))
        # min_components_that_protein_connects = ceil(constant + ratio * sqrt(num_components))

    else:
        min_components_that_protein_connects = max(min_components_that_protein_connects, ceil(constant + ratio * (num_components)))
        
    ### POPULATE COMPONENT DICTIONARY ###
    component_dictionary = dict() # protein : component_num
    j = 0
    for array in [(np.array(submatrix.get_list_of_proteins())[np.nonzero(labels == i)]) for i in range(num_components)]:
        for protein in array:
            component_dictionary[protein] = j
        j += 1
    
    ## FIND CONNECTED PROTEINS AND DETERMINE IF THEY QUALIFY 
    qualifying_proteins = list()

    for protein in (degreelist.get_list_of_proteins_sorted_by_degree()):   

        if (degreelist.get_degree_of_protein(protein, max_degree=max_degree) <= max_degree):
            num_edges, which_proteins = degreelist.determine_num_edges_to_cluster(protein, clusters.get_cluster_proteins(cluster_num), also_return_which_proteins=True)
                    
            if (num_edges >= min_components_that_protein_connects):
                set_of_components_that_protein_connects = degreelist.which_components_of_a_cluster_would_a_protein_connect(protein, clusters.get_cluster_proteins(cluster_num), component_dictionary, connected_proteins_within_cluster=which_proteins)

                if len(set_of_components_that_protein_connects) >= min_components_that_protein_connects:
                    qualifying_proteins.append(protein)
                

    return qualifying_proteins








# def pick_ratio(num_clusters: int):
#     """
#     will determine an approximate ratio to start with based on the total number of clusters
#     """
#     if (num_clusters > 1000):
#         return .5
#     elif num_clusters > 500:
#         return .7
#     elif num_clusters > 200:
#         return .9
#     elif num_clusters > 100:
#         return .925
#     elif num_clusters > 50:
#         return .995
#     else: 
#         return 1

def print_querylist_of_clusters_to_file(clusters: AllClusters, clusters_to_print: list(), query_filepath: str = "querylist.txt", proteins_to_add: dict() = dict()):
    """
    clusters_to_print -> specify a list of which clusters to print
    TODO
    proteins_to_add -> dictionary containing key: clusternum and value: list of proteins to add to that cluster
    """
    output_file = open(query_filepath, 'w')

    for cluster_num in clusters_to_print:
        for protein in clusters.get_cluster_proteins(cluster_num):
            output_file.write(f"{protein}\tcluster_{cluster_num}\n")
    
    if proteins_to_add: # dict not empty -> dict contains 
        for cluster_num in proteins_to_add: 
            for protein in proteins_to_add[cluster_num]: 
                output_file.write(f"{protein}\tcluster_{cluster_num}\n")

    output_file.close()


def get_initialized_fe(background_filepath: str, terms2features_filepath: str, termlist: pd.DataFrame() = vocabs.getTerms(['GO']), ecut: float = 0.01) -> FUNC_E():
    """TODO"""
    fe = FUNC_E()

    fe.importFiles({
        'background': background_filepath, 
        'terms2features': terms2features_filepath })
    fe.setTerms(termlist)
    fe.setEnrichmentSettings({'ecut': ecut})

    # now all that is left to do is upload the querylist using fe.importFiles({'query': querylist }), and running it, using fe.run(cluster=False)
    return fe


def create_term_mapping_list(go_terms_filepath: str, term_mapping_filepath: str = 'term_mapping.txt'):
    """
    the original file (go_terms_filepath) is in form GOTERM tab PROTEIN, while the term mapping file (term_mapping_filepath) is printed in form PROTEIN tab GOTERM. if a protein has multiple they appear on seperate lines
    """
    with open(term_mapping_filepath, 'w') as file:
        with open(go_terms_filepath, 'r') as go_annotation_file:
            for _ in range(1): # first line of file has column titles, and should be skipped
                next(go_annotation_file)
            for line in go_annotation_file:
                terms = line.split()
                file.write(f"{terms[1]}\t{terms[0]}\n")


            


########

def get_cluster_connectivity (
    matrix:ProteinMatrix,
    degreelist:DegreeList,
    clusters:AllClusters,
    added_proteins:dict={},
    percentages:bool=True,
    sort_it:bool=False,


):
    """
    returns a dictionary of cluster_num : percent_connectivity
    note: uses SubMatrix from matrix class

    can specify if you want sorted.
    if added_proteins is specified (not empty), then it will add those proteins to the cluster before calculating connectivity

    """
    proteins = matrix.get_list_of_proteins()
    degree_dict = dict(degreelist.sorted_protein_degree_dict)
    matrix_df = matrix.get_matrix()
    cluster_connectivity = {}

    for cluster_num in clusters.get_all_cluster_labels():
        # get all the proteins associated to a cluster  
        cluster_proteins = clusters.get_cluster_proteins(cluster_num)
        
        # added_cluster_proteins is empty in the case that none have been added, or if added proteins was not specified
        added_cluster_proteins = [] if not added_proteins or cluster_num not in added_proteins else added_proteins[cluster_num]
        # get the list of potential proteins to add to cluster 
        submatrix = SubMatrix(list(set(cluster_proteins + added_cluster_proteins)), matrix)
        components_and_labels = submatrix.get_num_components_and_labels()
        num_components = components_and_labels[0]
        # current ratio of clusters to proteins
        if (percentages):
            num_proteins = len(cluster_proteins)
            # The following is a linear relationship that guarantees that a cluster with 1 component is 100% connected, and a cluster with 0 edges between proteins is 0% connected  
            percent_connectivity = 1 - (num_components - 1) / (num_proteins - 1)
            cluster_connectivity[cluster_num] = percent_connectivity
        else:
            cluster_connectivity[cluster_num] = num_components

    if sort_it:
        sorted_cluster_connectivity:dict = {k: v for k, v in sorted(cluster_connectivity.items(), key=lambda item: item[1], reverse=False)}
        return sorted_cluster_connectivity
    
    return cluster_connectivity


def top_n_proteins(
        qualifying_proteins:dict,
        n:int # max # of proteins to return
):
    n_proteins = dict()
    for key in qualifying_proteins:
        n_proteins[key] = qualifying_proteins[key][0:n]
    
    return n_proteins


def calculate_connectivity_difference(connectivity_before:dict, connectivity_after:dict, sort_it:bool=True):
    """
    Calculates the difference between the connectivity of each cluster before and after the addition of a new protein
    """
    difference = {}
    for cluster in connectivity_before:
        difference[cluster] = connectivity_after[cluster] - connectivity_before[cluster]
    
    
    if sort_it:
        sorted_difference:dict= {k: v for k, v in sorted(difference.items(), key=lambda item: item[1], reverse=True)}
        return sorted_difference
    
    
    return difference
