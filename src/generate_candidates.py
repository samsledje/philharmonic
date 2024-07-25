# python src/philharmonic_candidates.py -o {output.candidates} --go_list {input.go_shortlist} --protein_list {input.protein_shortlist} {input.sequences}

import argparse
import numpy as np
from Bio import SeqIO

from utils import parse_GO_map
from itertools import combinations

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Generate a list of candidate PPIs')
    parser.add_argument('--sequences', required=True, type=str, help='Sequences file')
    parser.add_argument('-o', '--output', required=True, type=str, help='Output file')
    parser.add_argument('--paircount', type=int, help='Number of protein pairs to sample', default=10000000)
    parser.add_argument('--go_map', type=str, help='GO map file')
    parser.add_argument('--go_list', type=str, help='GO shortlist')
    parser.add_argument('--go_annot_wt', type=float, help='Relative weight of proteins with GO annotation', default=100)
    parser.add_argument('--protein_list', type=str, help='Protein shortlist')

    args = parser.parse_args()

    sequences = SeqIO.to_dict(SeqIO.parse(args.sequences, 'fasta'))
    protein_names = list(sequences.keys())
    candidates = list(combinations(sequences.keys(), 2))

    shortlisted_proteins = [line.strip() for line in open(args.protein_list,"r")]
    shortlisted_go = [line.strip() for line in open(args.go_list,"r")]
    go_map = parse_GO_map(args.go_map)



    with open(args.output, 'w') as f:
        for candidate in candidates:
            f.write(f'{candidate[0]}\t{candidate[1]}\n')

# def write_candidate_pair_list(df, nPointPairs, outfile2):
#     N = df.shape[0]
#     l0 = list(range(N))

#     dbg_print("Flag 569.10 ", df.shape, nPointPairs, outfile2)
#     j_u = np.random.choice(l0, 10*int(nPointPairs), True, (df["wt"]/df["wt"].sum()).values)
#     j_v = np.random.choice(l0, 10*int(nPointPairs), True, (df["wt"]/df["wt"].sum()).values)
#     valid = j_u < j_v  #get rid of potential duplicates (x1,x2) and (x2,x1) as well as (x1,x1)
#     i_u = j_u[valid]
#     i_v = j_v[valid]

#     V0 = df["protein"].values
#     df2 = pd.DataFrame( zip(V0[i_u], V0[i_v]), columns=["protein_A","protein_B"])
#     df2 = df2.drop_duplicates().reset_index(drop=True)
#     df2 = df2.head(nPointPairs)


#     dbg_print("Flag 569.30 ", len(j_u), sum(valid), len(i_u), df2.shape)
#     df2.to_csv(outfile2, index=False)
#     return df2



# def write_protein_list(outfile1, species):
#     dbg_print("Flag 469.01 ")

#     allowed_pfam_domains, pfam2go_map = identify_allowed_pfam_domains()
#     dbg_print("Flag 469.10 ", len(allowed_pfam_domains), len(pfam2go_map))
    
#     pfam2hmmscan = collate_hmmscan_results(species)
#     dbg_print("Flag 469.20 ", len(pfam2hmmscan))

#     manual_annot_proteins = {}
#     if species=="pdam":
#         receptor_proteins = [l.strip() for l in open(PROJ_DIR + "/data/raw/Judith_receptor_proteins_20200612.txt",'r') if l.strip()!=""]
#         woundhealing_proteins =  [l.strip() for l in open(PROJ_DIR + "/data/raw/Judith_woundhealing_proteins_20200612.txt",'r') if l.strip()!=""]
#         manual_annot_proteins["Receptor"] = receptor_proteins
#         manual_annot_proteins["Wound-heal"] = woundhealing_proteins
        
#         dbg_print("Flag 469.30 ", len(receptor_proteins), len(woundhealing_proteins))
        
#     elif species in ["cowV2"]:
#         rumen_proteins = [l.strip() for l in open("/afs/csail.mit.edu/u/s/samsl/Work/DSCRIPT_Validation_RECOMB/Cow_CellSystems/RumenSpecificProteins.txt",'r') if l.strip()!=""]
#         manual_annot_proteins["Rumen"] = rumen_proteins

        
#         dbg_print("Flag 469.31 ", len(rumen_proteins))

#     hdr, l = select_species_proteins(allowed_pfam_domains, pfam2go_map, pfam2hmmscan, manual_annot_proteins)

#     dbg_print ("Flag 469.50 ", hdr, len(l))
#     csvout = csv.writer(open(outfile1,'w'))
#     csvout.writerow(hdr)
#     csvout.writerows(l)
    
#     return hdr, l


    
# def run(args, extra_args):
    
#     dbg_print("Flag 401.01 ")

#     outfile1 = "{0}/proteins_list_{1}_{2}.csv".format(args.outdir, args.species, args.outsfx)
#     outfile2 = "{0}/candidate_pairs_{1}_{2}.csv".format(args.outdir, args.species, args.outsfx)
    
#     hdrs, l = write_protein_list(outfile1, args.species)
    
#     l1 = []
#     for p, ma, _, _ in l:
#         l1.append((p, 1 if not ma else args.manual_annot_wt))
        
#     df = pd.DataFrame(l1, columns=["protein","wt"]).sort_values("protein").reset_index(drop=True)

#     write_candidate_pair_list(df, int(args.paircount), outfile2)
    


    
# ##########################################################################

# if __name__ == "__main__":

#     parser = argparse.ArgumentParser()
#     parser.add_argument("--species", help="species to analyze. Data should be in PROJ_DIR/data/proessed . Current choices are pdam , hvul, tpseud, symbc1, or cow ", type=str, default='pdam') #choices=['pdam','hvul','tpseud','symbc1', 'cow'], default='pdam')
#     parser.add_argument("--paircount", help="number of protein pairs to sample", type=int, default=10000000)
#     parser.add_argument("--outdir", help="output directory (can set to '.')", type=str, default=PROJ_DIR+"/data/processed/")
#     parser.add_argument("--outsfx", help="suffix to use when producing output files")
#     parser.add_argument("--manual_annot_wt", help="relative weight of proteins with manual annotation", type=float, default=100)
#     parser.add_argument("--extra", help="put this as the LAST option and arbitrary space-separated key=val pairs after that", type=str, nargs='*')


#     args = parser.parse_args()
#     extra_args = dict([a.split("=") for a in args.extra]) if args.extra else {}


#     run(args, extra_args)

    
    
