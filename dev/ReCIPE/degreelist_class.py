"""
Author: Charlotte Versavel
Date:   June 2022
Last Edit: Oct 2022

                             degreelist_class.py

Purpose: a class TODO
TODO: figure out how to use (oct 2022)


"""

import pandas as pd 
import numpy as np

from matrix_class import ProteinMatrix




class DegreeList:
    
    """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
    * * * * * * * * * * * * * MEMBER VARIABLES * * * * * * * * * * * * * *  
    """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
    protein_matrix : ProteinMatrix = ProteinMatrix

    sorted_protein_degree_dict = dict()


    """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
    * * * * * * * * * * * * * * INITIALIZERS * * * * * * * * * * * * * * *  
    """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
    def __init__(self, matrix : ProteinMatrix) -> None:
        """            
        Parameters: matrix is populated with proteins and their interaction 
                    weights
        Purpose:    to take in a proteinMatrix (or submatrix) and create a 
                    sorted dictionary of protein:degree for all proteins in the 
                    matrix.
        Returns:    n/a
        """
        self.protein_matrix = matrix

        protein_degree_dict = {name:matrix.find_degree(name) for name in matrix.get_list_of_proteins()}

        self.sorted_protein_degree_dict = sorted(protein_degree_dict.items(), key=lambda x: x[1], reverse=True)


    def __repr__(self): 
        """             
        Purpose:    to override the print function for this class to print the 
                    sorted dictionary when called
        Returns:    a string of the dictionary
        """
        
        return str(self.sorted_protein_degree_dict)

    """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
    * * * * * * * * * * * * * * * GETTERS * * * * * * * * * * * * * * * * *  
    """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
    def get_degree_list(self) -> list():
        """             
        Purpose:    to allow access to the sorted degree list
        Returns:    a list of tuples of (protein, degree)
        """
        return self.sorted_protein_degree_dict


    def get_list_of_proteins_sorted_by_degree(self) -> list():
        """
        TODO
        """
        reverse_list_of_proteins = []

        for protein_degree_pair in reversed(self.sorted_protein_degree_dict):
            reverse_list_of_proteins.append(protein_degree_pair[0])

        return reverse_list_of_proteins

    def get_protein_at_index(self, index : int, degree = False) -> str or tuple:
        """             
        Parameters: index is the index of the protein in the sorted list
                    degree is a boolean that determines if the degree is returned as well
        Purpose:    to return the protein at the specified index
        Returns:    the protein at the specified index, or if degree is True, a 
                    tuple of (protein, degree)
        """
        if not degree:
            return self.sorted_protein_degree_dict[index][0]
        else:
            return self.sorted_protein_degree_dict[index]
    

    # def get_degree_of_protein(self, protein: str) -> int:
    #     """
    #     Parameters: 
    #         -   protein is the name of a protein whose degree should be determined
    #     Purpose:    to determine how well-connected a given protein is
    #     Returns:    protein's degree
    #     """
    #     try:
    #         degree = self.sorted_protein_degree_dict[protein]
    #         return degree
    #     except KeyError:
    #         print(f"get_degree_of_protein({protein}) caused an error because {protein} is not in the degreelist dictionary. returning -1")
    #         return -1

    
    
    def determine_num_edges_to_cluster(self, protein : str, cluster_list : list(), max_edges_until_return : int = -1, also_return_which_proteins: bool = False) -> int and list:
        """             
        Parameters: protein is a single protein in the matrix
                    cluster_list, a list of proteins in a cluster
                    max_edges_until_return allows the function to stop counting edges once a certain target is reached
                    also_return_which_proteins if set to true, will return which proteins have have edges to the given protein. should not be set to true when max_edges is set to a specific value
        Purpose:    to determine the number of edges between the protein and the proteins in the cluster
        Returns:    the number of edges, and identify_connections: a list of bools in the same order as proteins in the cluster, with true if a the protein is connected, false otherwise
        """

        for already_in_cluster in cluster_list:
            if protein == already_in_cluster:
                if also_return_which_proteins:
                    return 0, []
                return 0
        
        num_edges = 0
        which_proteins = list() 

        if max_edges_until_return == -1: # max_edges_until_return has been left unspecified
            i = 0
            for cluster_protein in cluster_list:
                if (self.protein_matrix).has_edge(protein, cluster_protein):
                    num_edges += 1
                    which_proteins.append(cluster_protein)
                i += 1
        else: # max_edges_until_return has been specified
            for cluster_protein in cluster_list:
                # print(f"about to call has edge 1")
                if (self.protein_matrix).has_edge(protein, cluster_protein):
                    num_edges += 1
                if num_edges >= max_edges_until_return:
                    if also_return_which_proteins:
                        return 0, []
                    return num_edges
        
        # print(f"num_edges: {num_edges}, which_proteins: {which_proteins}")

        if (also_return_which_proteins):
            return num_edges, which_proteins
        return num_edges
        


    def create_list_of_proteins_connected_to_cluster(self, list_of_proteins: np.array, cluster_list: np.array, max_list_length: int = -1, min_num_connections: int = 3) -> list:
        """             
        Parameters: cluster_list is a list of proteins in a cluster
                    max_list_length is an upper bound for the number of proteins to return in a list. If None, all proteins with at least min_num_connections connections are added to the list
                    min_num_connections is the minimum number of connections a protein must have to be added to the list and considered 'connected' to the cluster
        Purpose:    to create a list of proteins that are connected to the cluster
        Returns:    a list of proteins that are connected to the cluster
        """
        
        qualifying_proteins = []

        if max_list_length == -1: # max_list_length left unspecified
            for protein in list_of_proteins:
                num_edges = self.determine_num_edges_to_cluster(protein, cluster_list, max_edges_until_return=min_num_connections)

                if (num_edges >= min_num_connections):
                        qualifying_proteins.append(protein)
        else: # max_list_length has been specified    
            for protein in list_of_proteins:
                num_edges = self.determine_num_edges_to_cluster(protein, cluster_list, max_edges_until_return=min_num_connections)
                
                if (num_edges >= min_num_connections):
                        qualifying_proteins.append(protein)
                        if (len(qualifying_proteins) >= max_list_length):
                            return qualifying_proteins
            
        
        return qualifying_proteins
        

    def which_components_of_a_cluster_would_a_protein_connect(self, protein: str, cluster_of_proteins, protein_to_component_dict: dict(), connected_proteins_within_cluster:list = []) -> set:
        """
        Parameters: 
            -   cluster is a cluster containing a group of proteins that are 
                in some way related
            -   TODO
        Purpose:    to determine if a protein could re-attach different 
                    components of a cluster
        Returns:    a set of the components that would be connected by the 
                    given protein
        """
        which_components_were_connected = set()

        connected_proteins = connected_proteins_within_cluster
        if not connected_proteins:
            # obtain which cluster proteins the given protein is connected to
            num_edges, connected_proteins = self.determine_num_edges_to_cluster(protein, cluster_of_proteins, also_return_which_proteins = True)

        
        # for each protein that was connected, store it's component label
        # print(f"connected proteins: {connected_proteins}. labels: {protein_to_component_dict}")
        for protein in connected_proteins:
            try:
                which_components_were_connected.add(protein_to_component_dict[protein])
            except KeyError:
                pass
                # print(f"adding the following label {cluster_component_labels[i]}")

        return which_components_were_connected # type: set


    def determine_if_a_protein_will_connect_a_cluster(self, protein: str, cluster_of_proteins: list(), protein_to_component_dict: dict(), min_num_connections: int = 2) -> bool: 
        """
        Parameters: 
            -   cluster is a cluster containing a group of proteins that are 
                in some way related
            -   TODO -> dict
            -   min_num_connections is the minimum number of 
                components that need to be connected for a protein to be 
                considered successful in connecting a cluster
        Purpose:    to determine if the components that a protein could 
                    reconnect are satisfactory.
        Returns:    true if a protein connects >= 2 elements in a cluster
        """
        set_of_components = self.which_components_of_a_cluster_would_a_protein_connect(protein, cluster_of_proteins, protein_to_component_dict)

        if ((len(set_of_components)) >= min_num_connections):
            return True

        return False


    def get_degree_of_protein(self, protein: str, max_degree: int = 501):
        """ TODO """
        try:
            return dict(self.sorted_protein_degree_dict)[protein]
        except:
            print(f"{protein} not in degreelist dictionary")
            return (max_degree + 1)
            
