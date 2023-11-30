import numpy as np
import pandas as pd


class GO:
    def __init__(self, ID, features):
        self.ID = ID
        self.D = features
        self.name = features['name']

    def __repr__(self):
        return '{} - <{}>'.format(self.ID, self.name)

    def __eq__(self, other):
        return self.ID == other.ID

    def __hash__(self):
        return hash(self.ID)

def extract_GO_id_from_list(l):
    if isinstance(l,list):
        return [i.split('|')[0] for i in l]
    else:
        return None
    
def clean_GO_map(f):
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
    return seq2GO

def read_GO_obo(infile):
    terms = {}
    with open(infile,'r') as f:
        for line in f:
            tDict = {}
            line = line.strip()
            if line == "[Term]":
                line = f.readline().strip().split(': ')
                while not line == ['']:
                    tDict[line[0]] = ''.join(line[1:])
                    line = f.readline().strip().split(': ')
                for k,v in tDict.items():
                    k = k.strip()
                    v = v.strip()
                    tDict[k] = v
                terms[tDict['id']] = GO(tDict['id'], tDict)
    return terms

def GO_search(clusters, GO_term, GO_OBJECTS, N=20):
    if isinstance(GO_term,str):
        GO_term = GO_OBJECTS[GO_term]
    return [c for c in clusters if GO_term in [i[0] for i in c.get_top_terms(N)]]

def protein_search(clusters, protein_list):
    if isinstance(protein_list, str):
        protein_list = [protein_list]
    plist = [c for c in clusters if np.array([p in c.proteins for p in protein_list]).all()]
    if len(plist) == 0:
        return None
    elif len(plist) == 1:
        return plist[0]
    else:
        return plist
    
def triangle_search(clusters, min_triangles=0, max_triangles=np.infty):
    return [c for c in clusters if c.triangles() >= min_triangles and c.triangles() <= max_triangles]

def node_search(clusters, min_nodes=0, max_nodes=np.infty):
    return [c for c in clusters if len(c) >= min_nodes and len(c) <= max_nodes]