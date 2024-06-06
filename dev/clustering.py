# Philharmonic Clustering Pipeline
# Last modified: June 5, 2024

'''
Run with: python3 clustering.py --net-name SymbC1 --interaction-file ./data/SymbC1_predictions_positive.tsv --dsd-file ./data/SymbC1_dscript_distances.DSD1.tsv --results-dir results/ 
clustering: 
  num_init_clusts: 500
  min_clust_size: 3
  clust_divisor: 20
'''

import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import pickle as pk
import argparse
import json
import hashlib
import collections
from tqdm.notebook import tqdm
from sklearn.cluster import KMeans
from sklearn.cluster import SpectralClustering
from sklearn.metrics import silhouette_score
from networkx.linalg.graphmatrix import adjacency_matrix
from queue import PriorityQueue




def get_args():
    """
    Required Arguments:
        --net-name: Name of the network
        --interaction-file: Path to the interaction file
        --dsd-file: Path to the DSD file containing a matrix of distances between nodes


    Optional Arguments:
        --base-dir: Base directory for the project (default: ./)
        --edge-weight-thresh: Threshold for edge weights (default: 0.5)
        --rand-seed: Random seed for reproducibility (default: 6191998)
        --min-clust-size: Minimum cluster size (default: 3)
        --go-db-file: Path to the GO database file (default: go.obo)
        --n-init-clusters: Number of clusters for initial spectral clustering (default: len(G.nodes) // 20)
        --cluster-divisor: Number of clusters to split a too-big cluster into (default: 20)
        --results-dir: Directory to output experiment results into (default: ./)

    """
    parser = argparse.ArgumentParser()
    ## Directories and Network Name:
    parser.add_argument("--base-dir", type=str, default="./", help="Base directory for the project", required=False)
    parser.add_argument("--results-dir", type=str, default="./", help="Directory to output experiment results into (ie results)", required=False)
    parser.add_argument("--net-name", type=str, help="Name of the network (to label output)", required=True)

    ## Input Files:
    parser.add_argument("--go-db-file", type=str, help="Path to the GO database file", required=False, default="go.obo")
    parser.add_argument("--interaction-file", type=str, help="Path to the interaction file", required=True)
    parser.add_argument("--dsd-file", type=str, help="Path to the DSD file", required=True)

    ## Spectral Clustering Args:
    parser.add_argument("--n-init-clusters", type=int, help="Number of clusters for initial spectral clustering", required=False, default=None)
    parser.add_argument("--edge-weight-thresh", type=float, help="Threshold for edge weights", required=False, default=0.5)
    parser.add_argument("--rand-seed", type=int, help="Random seed for reproducibility", required=False, default=6191998)
    
    parser.add_argument("--min-clust-size", type=int, help="Minimum cluster size", required=False, default=3)
    parser.add_argument("--cluster-divisor", type=int, help="Number of clusters to recursively split a too-big cluster into", required=False, default=20)

    return parser.parse_args()



def read_network(DSDfile, interactionFile, edge_weight_thresh=0.5, net_name="Network", output_stats=True, results_dir='./'):
    '''
    Read in the network and filter edges based on DSD confidence threshold
    '''
    print(f'Reading DSD File: {DSDfile}...')
    dsd_df = pd.read_csv(DSDfile, sep='\t', index_col=0, header=0)
    protein_names = [str(i) for i in dsd_df.index]
    DSD = dsd_df.values

    fullG = nx.read_weighted_edgelist(interactionFile)
    print('Selecting DSD connected component...')
    G = fullG.subgraph(protein_names)
    print('Filtering edges with confidence threshold {}...'.format(edge_weight_thresh))
    wG = nx.Graph()
    for (u,v,d) in tqdm(G.edges.data()):
        if d['weight'] >= edge_weight_thresh:
            wG.add_edge(u,v,weight=d['weight'])
    del G
    G = wG 
    A = nx.to_numpy_array(G, nodelist=protein_names)
    degrees = [i[1] for i in list(G.degree())]

    if output_stats:
        # print a table of network statistics
        label = ['Nodes','Edges','Degree (Med)','Degree (Avg)','Sparsity']
        value = [len(G.nodes), len(G.edges), np.median(degrees), np.mean(degrees), len(G.edges()) / len(G)**2]
        stats = pd.DataFrame([label,value]).T
        stats.columns = ['',net_name]
        stats = stats.set_index('')
        print(stats)
        # save a histogram of protein degree 
        display_degree_dist(G, degrees, net_name, results_dir)

    return G, DSD, protein_names

def display_degree_dist(G, degrees, net_name, results_dir):
    '''
    Helper function to read_network that displays the degree distribution of the initial network
    '''
    degreeDist = {}
    for i in degrees:
        n = degreeDist.setdefault(i,0)
        degreeDist[i] = n + 1

    plt.xlabel('Degree')
    plt.ylabel('Proportion of Nodes')  # we already handled the x-label with ax1
    plt.title('Node Degree Distribution')
    plt.scatter(degreeDist.keys(), [i/len(G) for i in degreeDist.values()])
    print('Saving Degree Distribution Plot to: ', results_dir + net_name + '.degree_dist.png')
    plt.savefig(results_dir + net_name + '_degree_dist.png')

def RBF(D, sigma=None):
    """
    Convert distance matrix D into similarity matrix S using Radial Basis Function (RBF) Kernel
    RBF(x,x') = exp( -((x - x')**2 / 2sigma**@))
    """
    sigma = sigma or np.sqrt(np.max(D))
    return np.exp(-1 * (np.square(D) / (2 * sigma**2))) 

def compute_similarity_score(DSD):
  print('Computing similarity scores...')
  simDSD = RBF(DSD)
  print('Sparsifying similarity scores...')
  sparse_sim_thresh = 1e-5
  simRav = simDSD.ravel()
  simRav[simRav < sparse_sim_thresh] = 0
  simRav = simRav.reshape(simDSD.shape)
  simDSD = simRav
  return simDSD


def print_silhouette_stats(G, simDSD, SC): 
    print('Silhouette Score: ', silhouette_score(simDSD,SC.labels_))
    AM = (adjacency_matrix(G).todense() != 0).astype(int)
    print('Adjacency Martrix Silhouette Score: ', silhouette_score(AM,SC.labels_)) # TODO: IDK WHAT THIS IS



def extract_clusters(SC: SpectralClustering, og_clust:list, og_csize:int, verbosity=0): # TODO: convert Verbosity to logging (options for verbosity are 0, 1, 2)
    subCs = []
    for label in set(SC.labels_): # go through all labels
        # make a list of all nodes in the cluster assigned to given label
        subCs.append([og_clust[i] for i in range(og_csize) if SC.labels_[i] == label])
    
    if verbosity > 0:
        to_print = f"Cluster of size {og_csize} -> split into {len(subCs)} clusters "
        
        if verbosity > 1:
            to_print += "of sizes: "
            for subC in subCs:
                to_print += f"{len(subC)}/"        
        print(to_print[:-1])

    # return list of tuples of priority of cluster and list of nodes in cluster
    return [(1/len(subC), subC) for subC in subCs]


def cluster_network(G, simDSD, min_size, cluster_div, kclusts, rand_seed, results_dir, net_name, protein_names):
    
    # INITIAL CLUSTERING:
    print(f'Fitting {kclusts} spectral clusters...')
    SC = SpectralClustering(n_clusters=kclusts, assign_labels="discretize", random_state=rand_seed, affinity='precomputed')
    SC.fit(simDSD)

    initial_SC_pickle_file = f"{results_dir}/{net_name}_initial_SC.pkl"
    with open(initial_SC_pickle_file,'wb') as f:
      print('Saving initial spectral clustering to: ', initial_SC_pickle_file)
      pk.dump(SC, f)

    print_silhouette_stats(G, simDSD, SC)

    with open(initial_SC_pickle_file,'rb') as f:
      print('Loading initial spectral clustering from: ', initial_SC_pickle_file)
      clusteringToUse = pk.load(f)

    simMatrix = simDSD

    clusts = [[j for j in range(len(clusteringToUse.labels_)) if clusteringToUse.labels_[j] == i] for i in range(max(clusteringToUse.labels_)+1) if i in clusteringToUse.labels_]
    clusts.sort(key = lambda x: len(x), reverse=True)


    clustQ = PriorityQueue()
    for c in clusts:
        clustQ.put((1/len(c), c))
        
    print(f'Splitting large clusters...')
    while True:
        priority, c = clustQ.get()
        csize = int(1/priority)

        n_clusters = int(np.round(csize / cluster_div)) # NOTE: this changed so that we are left with clusters of max size 30

        if n_clusters < 2:
            clustQ.put((priority, c))
            break
      
        SC2 = SpectralClustering(n_clusters=n_clusters, assign_labels="discretize", random_state=rand_seed, affinity='precomputed')
        SC2.fit(simMatrix[c,:][:,c])
        
        subClusts = extract_clusters(SC2, c, csize, verbosity=2)
        for subClust in subClusts:
            clustQ.put(subClust)

    print(f'Removing small clusters (<{min_size})...')
    filteredClusters = []
    while not clustQ.empty():
        wght, c = clustQ.get()
        filteredClusters.append(c)
    filteredClusters = [i for i in filteredClusters if len(i) >= min_size]
    
    print(f'Final Clustering: {len(filteredClusters)} clusters')

    spectral_clustering_results(filteredClusters, results_dir)

    clustsNames = [[protein_names[i] for i in cl] for cl in filteredClusters]
    clustsNames.sort(key = lambda x: len(x), reverse=True)

    initial_clusters_output_file = f"{results_dir}/{net_name}_clusters.csv"
    print('Writing clusters to: ', initial_clusters_output_file)
    writeClusters(initial_clusters_output_file, clustsNames)

    return filteredClusters

def writeClusters(outfile, clusts):
    with open(outfile, 'w+') as f:
        for cl in clusts:
            f.write('{}\n'.format(','.join([str(i) for i in cl])))

def spectral_clustering_results(filteredClusters, results_dir):
    sizes = [len(i) for i in filteredClusters]
    plt.hist(sizes, bins=range(min(sizes), max(sizes)+1))
    plt.xlabel('Cluster Size')
    plt.ylabel('Number of Clusters')
    plt.title('Cluster Size Distribution')
    plt.savefig(f"{results_dir}/clust_size_dist.png")

def main(args):

    sns.set_theme(style='white')

    # DATA PARAMS
    BASE_DIR = args.base_dir
    RES_DIR = f'{BASE_DIR}/{args.results_dir}/'
    DATA_DIR = f'{BASE_DIR}/data/' # TODO: make this a command line argument {args.data_dir}
    interactionFile = f'{BASE_DIR}/{args.interaction_file}'
    DSDfile = f'{BASE_DIR}/{args.dsd_file}' 
    # go_db_file =f'{BASE_DIR}/{args.go_db_file}' 

    min_size = args.min_clust_size
    cluster_div = args.cluster_divisor

    
    print(RES_DIR)
    G, DSD, protein_names = read_network(DSDfile, interactionFile, edge_weight_thresh=args.edge_weight_thresh, net_name=args.net_name, output_stats=True, results_dir=RES_DIR)
    kclusts = args.n_init_clusters 
    if kclusts is None: kclusts = len(G.nodes) // 20

    DSD = compute_similarity_score(DSD)
    _ = cluster_network(G, DSD, min_size, cluster_div, kclusts, args.rand_seed, RES_DIR, args.net_name, protein_names)
    
    
    
    # if kclusts is None: kclusts = len(G.nodes) // 20




if __name__ == "__main__":
    main(get_args())

