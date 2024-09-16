#python src/build_cluster_graph.py -o {output.graph} -cfp {input.clusters} -nfp {input.network}
import argparse
import networkx as nx
import pandas as pd
from collections import Counter
from tqdm import tqdm

from utils import Cluster, parse_GO_map, parse_GO_database, load_cluster_json, clean_top_terms, log


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Create a cluster-of-clusters network to visualize overall organization"
    )
    parser.add_argument("-o", "--output", type=str, help="Output file")
    parser.add_argument("--coc_functions", type=str, help="Function of clusters of clusters")
    parser.add_argument(
        "-cfp", "--cluster_file_path", type=str, help="Cluster file path"
    )
    parser.add_argument("-nfp", "--network_file_path", type=str, help="Network file path")
    parser.add_argument("--go_map", required=True, type=str, help="GO map file")
    parser.add_argument("--go_db", required=True, type=str, help="GO database file")
    parser.add_argument("--n_crossing_edges", type=int, default=10, help="Clusters with n crossing edges are connected")
    
    args = parser.parse_args()

    # Load clusters
    cluster_dict = load_cluster_json(args.cluster_file_path)
    clusters = [Cluster(i) for i in cluster_dict.values()]

    # Load network
    full_G = nx.read_weighted_edgelist(args.network_file_path)
    
    # Add GO Annotations
    go_map = parse_GO_map(args.go_map)
    go_db = parse_GO_database(args.go_db)

    for clust in clusters:
        clust.add_GO_terms(go_map, go_db)
        for gt in clust.GO_terms.keys(): assert gt in go_db.keys()
            
    cluster_top_terms = {hash(c): clean_top_terms(c, go_db, return_counts=False) for c in clusters}

    counter = Counter(cluster_top_terms.values())
    print(counter)

    def test_cluster_connection(c1_graph, c2_graph, parent_graph):
        return len(list(nx.edge_boundary(parent_graph, c1_graph.nodes, c2_graph.nodes)))
    
    def edge_filter(x1,x2,n_crossing):
        return clusG[x1][x2].get("weight") > n_crossing
    
    clusG = nx.Graph()

    for c in clusters:
        clusG.add_node(c, size=len(c))
    
    for n in tqdm(clusG.nodes):
        for n2 in clusG.nodes:
            # thresh = (len(n) + len(n2) // 4)
            if n != n2:
                n_neighbors = test_cluster_connection(n.G, n2.G, full_G)
                clusG.add_edge(n, n2, weight=n_neighbors)
                    
    clusG = nx.relabel_nodes(clusG, {i: hash(i) for i in clusG.nodes})
    subG = nx.subgraph_view(clusG, filter_edge=lambda x1, x2: edge_filter(x1, x2, args.n_crossing_edges))

    log(f"Cluster graph has {subG.number_of_nodes()} nodes and {subG.number_of_edges()} edges.")

    edge_table = nx.to_pandas_edgelist(subG)
    
    # Nodes
    size_table = pd.DataFrame(nx.get_node_attributes(subG, "size").items(),columns=["key","size"])
    fn_table = pd.DataFrame(cluster_top_terms.items(),columns=["key","go_fn"])
    node_table = pd.merge(size_table, fn_table)

    n_singletons = len(subG.nodes) - len(set(edge_table['source'].values).union(edge_table['target'].values))
    log(f"There are {n_singletons} singleton clusters that don't have more than {args.n_crossing_edges} edges to any other cluster")

    edge_table.to_csv(args.output,index=False)
    node_table.to_csv(args.coc_functions,index=False)

    

    