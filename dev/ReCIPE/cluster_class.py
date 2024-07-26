"""
Author: Charlotte Versavel
Date:   June 2022
Last Edit: Nov 2022

                             cluster_class.py

Purpose: a class to store the protein clusters and allow for access of a 
         specific cluster. 
         Also, allows 

"""

import pandas as pd 
import numpy as np
from collections import defaultdict

from sklearn import cluster

class AllClusters:

    """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
    * * * * * * * * * * * * * MEMBER VARIABLES * * * * * * * * * * * * * *  
    """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

    # clusters = defaultdict(lambda: []) # a dict of relation {cluster_num : list_of_proteins_in_cluster}
    
    """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
    * * * * * * * * * * * * * * INITIALIZERS * * * * * * * * * * * * * * *  
    """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

    def __init__(self, csv_filename: str = "", csv_clusters_have_labels: bool = False, protein_to_cluster_dict: dict = {}) -> None:
        """  
        Parameters: csv_filename is the name of a csv file containing several 
                    clusters of proteins 
                        in the form [cluster_identifier]  Protein1    Protein2 ... (note that cluster_identifiers may not be present, and if not, they will be assigned in order of appearance)
                    protein_to_cluster_dict is a dictionary with the form { protein : cluster_identifier }
        Purpose:    to populate several single clusters with data from a CSV 
                    file, or from a dictionary
        Returns:    n/a
        """

        self.clusters = defaultdict(lambda: [])

        if csv_filename != "":
            try:
                with open(csv_filename, "r") as data:
                    
                    for i, line in enumerate(data):
                        list_of_proteins = line.strip().split(",")
                        cluster_id = None
                        if csv_clusters_have_labels:
                            cluster_id = list_of_proteins.pop(0)
                        else:
                            cluster_id = i
                       
                        self.clusters[cluster_id] = list_of_proteins

            except FileNotFoundError:
                print(f"ERROR! file: {csv_filename} not found.")
        
        elif protein_to_cluster_dict: # dictionary not empty
            for protein in protein_to_cluster_dict.keys():
                self.add_protein_to_cluster(protein, int(protein_to_cluster_dict[protein]))
        
        else: # no filename or dictionary passed in
            print(f"ERROR! please specify a [csv_filename] or a [protein_to_cluster_dict] to initialize the clusters.")
            



    def __repr__(self): 
        """             
        Purpose:    Overloaded Print function - Prints a message indicating how to print clusters
        Returns:    a new message to print
        """
        return f"AllClusters has {len(self.clusters)} clusters (use the print_all method to see them)"

    """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
    * * * * * * * * * * * * * * * SETTERS * * * * * * * * * * * * * * * * *  
    """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

    def add_protein_to_cluster(self, protein:str, cluster_num) -> None:
        """             
        Parameters: 
            -   protein is the protein to add to a specified cluster
            -   cluster_num is the num of the cluster to add a protein to
        Purpose:    to add a protein to a cluster
        Returns:    n/a
        """
        self.clusters[cluster_num].append(protein)
        # print(f"appended cluster {cluster_num}: {self.clusters[cluster_num]}")

    def sort_dictionary(self) -> None:
        """             
        Purpose:    to sort the dictionary by number of proteins in each cluster
        Returns:    n/a
        """
        sorted_clusters = dict(sorted(self.clusters.items(), key=lambda x: len(x[1])))
        self.clusters = sorted_clusters
        # print(f"appended cluster {cluster_num}: {self.clusters[cluster_num]}")


    
    """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
    * * * * * * * * * * * * * * * GETTERS * * * * * * * * * * * * * * * * *  
    """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

    def get_cluster_proteins(self, cluster_number) -> list:
        """             
        Parameters: cluster_number is the number of the cluster to get
        Purpose:    to get the list of proteins from a cluster
        Returns:    the list of proteins in the cluster
        """

        return self.clusters[cluster_number]

    def get_num_clusters(self) -> int:
        """
        Purpose:    to access the number of clusters
        Returns:    the number of clusters
        """
        return len(self.clusters)

    def get_all_cluster_labels(self) -> list():
        """
        Purpose:    to access all labels (cluster nums)
        Returns:    the labels of the clusters
        """
        return self.clusters.keys()

    def get_all_clusters(self) -> dict():
        """
        Purpose:    to access all of the clusters
        Returns:    all clusters in format {cluster_num: [list_of_proteins]}
        """
        return dict(self.clusters)


    def print_all(self) -> None:
        """             
        Purpose:    to print all the clusters in the dictionary
        Returns:    n/a
        """
        print(self.clusters.keys())
        
        for cluster_num in self.clusters.keys():
            print(f"Cluster {cluster_num}: {self.get_cluster_proteins(cluster_num)}")
    
    
    def filter_clusters_by_size(self, min_size, max_size):
        """             
        Purpose:    to retrieve a dictionary that only contains clusters within a certain size range
        Returns:    dictionary
        """
        filtered_dict = {key: value for key, value in self.clusters.items() if min_size <= len(value) <= max_size}
        return filtered_dict



