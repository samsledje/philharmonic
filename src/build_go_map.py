import argparse
import pandas as pd
from collections import defaultdict
from Bio import SearchIO

if __name__ == "__main__":
    parser = argparse.ArgumentParser('Add GO terms to clusters')
    parser.add_argument('-o', '--output', required=True, type=str, help='Output file')
    parser.add_argument('--hhtblout', required=True, type=str, help='hhsearch tblout file')
    parser.add_argument('--pfam_go_files', type=str, nargs='+')

    args = parser.parse_args()

    print(args.output)
    print(args.hhtblout)

    hits = defaultdict(list)
    with open(args.hhtblout) as handle:
        for queryresult in SearchIO.parse(handle, 'hmmer3-tab'):
            for hit in queryresult.hits:
                prot_id = hit.query_id
                pfam_id = hit.accession.split('.')[0]
                hits[prot_id].append(pfam_id)

    pfam_go_map = pd.concat([pd.read_csv(f, sep=';', header=0) for f in args.pfam_go_files], ignore_index=True)

    protein_go_map = []
    for prot_id, pfam_hits in hits.items():
        go_terms = pfam_go_map[pfam_go_map['PFAM'].isin(pfam_hits)]['GO'].tolist()
        protein_go_map.append((prot_id, '', ";".join(pfam_hits), ";".join(go_terms)))
    protein_go_map = pd.DataFrame(protein_go_map, columns=['prot_id', 'manual annot', 'pfam_list', 'GO_list'])

    protein_go_map.to_csv(args.output, index=False)