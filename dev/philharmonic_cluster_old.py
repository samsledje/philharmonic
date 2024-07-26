import numpy as np
import pandas as pd
import networkx as nx
import hashlib
import collections
import matplotlib.pyplot as plt

from add_cluster_functions import GO

class Cluster:
    def __init__(self,proteins):
        self.proteins = proteins

    def __len__(self):
        return len(self.proteins)

    def __repr__(self):
        reprStr = "Cluster of {} [{},{},...] (hash {})".format(len(self), self.proteins[0], self.proteins[1], hash(self))
        if hasattr(self, 'G'):
            reprStr += "\nTriangles: {}\nMax Degree: {}".format(self.triangles(), max(self.G.degree(), key=lambda x: x[1])[1])
        if hasattr(self, 'GO_terms'):
            reprStr += "\nTop Terms:\n\t{}".format('\n\t'.join(
                    ['{} ({})'.format(i[0], i[1]) for i in self.get_top_terms(5)]
            ))
        return reprStr

    def __hash__(self):
        return int(hashlib.md5(''.join(self.proteins).encode()).hexdigest(), 16)

    def to_dict(self):
        D = {}
        D['id'] = hash(self)
        D['proteins'] = []
        for p in self.proteins:
            pD = {}
            pD['name'] = p
            if hasattr(self, 'GO_DB'):
                pD['go'] = self.GO_DB[self.GO_DB['seq'] == p]['GO_ids'].values[0]
            D['proteins'].append(pD)
        if hasattr(self, 'GO_DB'):
            D['go'] = sorted([{"id": i.ID, "desc": i.name, "freq": self.GO_terms[i]} for i in self.GO_terms], key = lambda x: x['freq'], reverse=True)
        if hasattr(self,'G'):
            D['graph'] = list(self.G.edges())
        return D

    def to_json(self):
        return json.dumps(self.to_dict())

    def add_GO_terms(self, go_db, GO_OBJECTS):
        self.GO_terms = {}
        self.GO_DB = go_db
        for prot in self.proteins:
            goIds = go_db[go_db['seq'] == prot]['GO_ids'].values[0]
            if goIds is None or len(goIds) == 0:
                continue
            for gid in goIds:
                try:
                    goObj = GO_OBJECTS[gid]
                except KeyError:
                    GO_OBJECTS[gid] = GO(gid,{'id':gid,'name':gid})
                    goObj = GO_OBJECTS[gid]
                goCount = self.GO_terms.setdefault(goObj,0)
                self.GO_terms[goObj] = goCount + 1

    def get_proteins_by_GO(self, GO_id):
        return [p for p in self.proteins if GO_id in prot_go_db.loc[p,'GO_ids']]

    def get_GO_by_protein(self, protein):
        assert protein in self.proteins, "{} not in cluster".format(protein)
        return [gt for gt in coi.GO_terms if gt.ID in prot_go_db.loc[protein,'GO_ids']]

    def get_top_terms(self,N):
        if not hasattr(self, 'GO_terms'):
            raise NotImplementedError("GO Terms have not been added yet.")
        GOlist = list(self.GO_terms.keys())
        if N == -1:
            N = len(GOlist)
        sortedList = sorted(GOlist,key=lambda x: self.GO_terms[x],reverse=True)[:N]
        return list(zip(sortedList, [self.GO_terms[i] for i in sortedList]))

    def set_graph(self,G):
        self.G = G.subgraph(self.proteins)

    def triangles(self):
        return sum([i for i in nx.triangles(self.G).values()]) / 3

    def draw_degree_histogram(self,draw_graph=True):
        if not hasattr(self,'G'):
            raise ValueError('Run .set_graph() method on this cluster first')
        G = self.G
        degree_sequence = sorted([d for n, d in G.degree()], reverse=True)  # degree sequence
        degreeCount = collections.Counter(degree_sequence)
        deg, cnt = zip(*degreeCount.items())

        fig, ax = plt.subplots()
        plt.bar(deg, cnt, width=0.80, color='b')

        plt.title("Degree Histogram")
        plt.ylabel("Count")
        plt.xlabel("Degree")
        ax.set_xticks([d + 0.4 for d in deg])
        ax.set_xticklabels(deg)

        # draw graph in inset
        if draw_graph:
            plt.axes([0.4, 0.4, 0.5, 0.5])
            pos = nx.spring_layout(G, k=0.15,iterations=10)
            plt.axis('off')
            nx.draw_networkx_nodes(G, pos, node_size=20)
            nx.draw_networkx_edges(G, pos, alpha=0.4)
        plt.show()

    def draw_graph(self):
        if not hasattr(self,'G'):
            raise ValueError('Run .set_graph() method on this cluster first')
        G = self.G
        nx.draw_kamada_kawai(G, with_labels=True,node_size=600, font_size=8)
    
def readClusterObjects(infile,sep=','):
    clusts = []
    with open(infile,'r') as f:
        for line in f:
            clusts.append(Cluster(line.strip().split(sep)))
    return clusts

def writeClusters(outfile, clusts):
    with open(outfile,'w+') as f:
        for cl in clusts:
            f.write('{}\n'.format(','.join([str(i) for i in cl])))

def cluster_from_json(jsonString, GO_OBJECTS):
        clust = Cluster([])
        D = json.loads(jsonString)
        clust.proteins = [i['name'] for i in D['proteins']]
        clust.GO_terms = {}
        for goDict in D['go']:
            gid = goDict['id']
            gdesc = goDict['desc']
            try:
                goObj = GO_OBJECTS[gid]
            except KeyError:
                GO_OBJECTS[gid] = GO(gid,{'id':gid,'name':gdesc})
                goObj = GO_OBJECTS[gid]
            clust.GO_terms[goObj] = goDict['freq']
        try:
            edgeList = D['graph']
            G = nx.Graph()
            for e in edgeList:
                G.add_edge(*e)
            clust.G = G
        except KeyError:
            pass
        return clust