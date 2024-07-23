import hashlib
import json
import re
import pandas as pd
import networkx as nx

class Cluster:
    def __init__(self,cluster_dict):
        for k,v in cluster_dict.items():
            setattr(self,k,v)
        assert hasattr(self,'members'), "Cluster must have 'members' attribute"

        if hasattr(self, "graph"):
            self.G = nx.Graph()
            for edge in self.graph:
                self.G.add_edge(edge[0], edge[1], weight=edge[2])

    def __len__(self):
        return len(self.members)
    
    def __repr__(self):
        reprStr = "Cluster of {} [{}] (hash {})".format(len(self), self._member_list(), hash(self))
        if hasattr(self, 'recipe'):
            pass
        if hasattr(self, 'G'):
            reprStr += "\nTriangles: {}\nMax Degree: {}".format(self.triangles(), max(self.G.degree(), key=lambda x: x[1])[1])
        if hasattr(self, 'GO_terms'):
            reprStr += "\nTop Terms:\n\t{}".format('\n\t'.join(
                    # ['{} ({})'.format(i[0], i[1]) for i in self.get_top_terms(5)]
                    ['{}'.format(i) for i in self.get_top_terms(5)]
            ))
        return reprStr

    def __hash__(self):
        return hash_cluster(self.members)
    
    def _member_list(self):
        if len(self.members) < 3:
            return ", ".join(self.members)
        else:
            return "{}, ...".format(", ".join(self.members[:3]))

    # def to_dict(self):
    #     D = {}
    #     D['id'] = hash(self)
    #     D['proteins'] = []
    #     for p in self.members:
    #         pD = {}
    #         pD['name'] = p
    #         if hasattr(self, 'GO_DB'):
    #             pD['go'] = self.GO_DB[self.GO_DB['seq'] == p]['GO_ids'].values[0]
    #         D['proteins'].append(pD)
    #     if hasattr(self, 'GO_DB'):
    #         D['go'] = sorted([{"id": i.ID, "desc": i.name, "freq": self.GO_terms[i]} for i in self.GO_terms], key = lambda x: x['freq'], reverse=True)
    #     if hasattr(self,'G'):
    #         D['graph'] = list(self.G.edges())
    #     return D

    # def to_json(self):
    #     return json.dumps(self.to_dict())

    def add_GO_terms(self, go_map, go_db):
        self.GO_terms = {}
        self.go_db = go_db.copy()
        for prot in self.members:
            goIds = go_map[go_map['seq'] == prot]['GO_ids'].values[0]
            if goIds is None or len(goIds) == 0:
                continue
            for gid in goIds:
                goCount = self.GO_terms.setdefault(gid,0)
                self.GO_terms[gid] = goCount + 1

    # def get_proteins_by_GO(self, GO_id):
    #     return [p for p in self.members if GO_id in prot_go_db.loc[p,'GO_ids']]

    # def get_GO_by_protein(self, protein):
    #     assert protein in self.members, "{} not in cluster".format(protein)
    #     return [gt for gt in coi.GO_terms if gt.ID in prot_go_db.loc[protein,'GO_ids']]

    def get_top_terms(self,N):
        if not hasattr(self, 'GO_terms'):
            raise NotImplementedError("GO Terms have not been added yet.")
        GOlist = list(self.GO_terms.keys())
        if N == -1:
            N = len(GOlist)
        sortedList = sorted(GOlist,key=lambda x: self.GO_terms[x],reverse=True)[:N]
        return list(zip(sortedList, [self.GO_terms[i] for i in sortedList]))

    def set_graph(self,G):
        self.G = G.subgraph(self.members)

    def triangles(self):
        return int(sum([i for i in nx.triangles(self.G).values()]) / 3)

    # def draw_degree_histogram(self,draw_graph=True):
    #     if not hasattr(self,'G'):
    #         raise ValueError('Run .set_graph() method on this cluster first')
    #     G = self.G
    #     degree_sequence = sorted([d for n, d in G.degree()], reverse=True)  # degree sequence
    #     degreeCount = collections.Counter(degree_sequence)
    #     deg, cnt = zip(*degreeCount.items())

    #     fig, ax = plt.subplots()
    #     plt.bar(deg, cnt, width=0.80, color='b')

    #     plt.title("Degree Histogram")
    #     plt.ylabel("Count")
    #     plt.xlabel("Degree")
    #     ax.set_xticks([d + 0.4 for d in deg])
    #     ax.set_xticklabels(deg)

    #     # draw graph in inset
    #     if draw_graph:
    #         plt.axes([0.4, 0.4, 0.5, 0.5])
    #         pos = nx.spring_layout(G, k=0.15,iterations=10)
    #         plt.axis('off')
    #         nx.draw_networkx_nodes(G, pos, node_size=20)
    #         nx.draw_networkx_edges(G, pos, alpha=0.4)
    #     plt.show()

    def draw_graph(self):
        if not hasattr(self,'G'):
            raise ValueError('Run .set_graph() method on this cluster first')
        G = self.G
        nx.draw_kamada_kawai(G, with_labels=True,node_size=600, font_size=8)
    

def log(x):
    print(x)

def hash_cluster(protein_list):
    return int(hashlib.md5(''.join(sorted(protein_list)).encode()).hexdigest(), 16) % (2**61-1)

def load_cluster_json(infile, return_objects = False):
    with open(infile,'r') as f:
        clusters = json.load(f)
    if return_objects:
        return {k: Cluster(v) for k,v in clusters.items()}
    else:
        return clusters

def parse_GO_database(infile):
    terms = {}
    with open(infile,'r') as f:
        for line in f:
            line = line.strip()
            if line == "[Term]":
                term_ids = []
                line = f.readline().strip()
                while not line == '':
                    if line.startswith('id:') or line.startswith('alt_id:'):
                        term_ids.append(line.split(': ')[1])
                    elif line.startswith('name:'):
                        main_name = line.split(': ')[1]
                    line = f.readline().strip()
                for tid in term_ids:
                    terms[tid] = main_name
    return terms


def parse_GO_map(f):
    seqDb = pd.read_csv(f,sep=',')
    seqDb.columns = ['seq','manual_annot','pfam_list','GO_list']
    seqDb['GO_str'] = seqDb['GO_list']
    seqDb['GO_list'] = seqDb['GO_str'].str.split(';')
    def extract_GO_id_from_list(l):
        if isinstance(l,list):
            return [i.split('|')[0] for i in l]
        else:
            return None
    seqDb['GO_ids'] = seqDb['GO_list'].apply(extract_GO_id_from_list)
    seq2GO = seqDb[['seq','GO_ids']]
    seq2GO.columns = ['seq','GO_ids']
    go_map = {}
    for _, r in seq2GO.iterrows():
        if r.GO_ids is not None:
            go_map[r.seq] = r.GO_ids
    return go_map

def repr_cluster(c):
    c = Cluster(c)
    return str(c)
    