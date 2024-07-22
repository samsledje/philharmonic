# philharmonic
Pipeline for decoding the functional networks of non-model organisms

Required Inputs: protein sequences
Optional Inputs: shortlist proteins, shortlist GO terms

## Steps

### Download necessary files
- go.obo https://current.geneontology.org/ontology/go.obo
- pfam go associations
    - http://godm.loria.fr/data/pfam_gomf_most_specific.txt
    - http://godm.loria.fr/data/pfam_gobp_most_specific.txt
    - http://godm.loria.fr/data/pfam_gocc_most_specific.txt

### Working backwards

- download go database
    - input: None
    - output: go database, pfam go associations
    - exec: wgets
- add pfam to sequences
    - input: protein sequences
    - output: protein-pfam associations
    - exec: hmmscan
- add functions to sequences
    - input: protein-pfam associations, pfam-go associations
    - output: protein-go associations
    - exec: python philharmonic_go.py, 
- generate candidate pairs
    - input: protein sequences, 
- predict network
    - input: candidate pairs
    - output: predicted network
    - exec: dscript
- compute dsd distances
    - input: predicted network
    - output: dsd matrix
    - exec: fastdsd
- cluster dsd matrix
    - input: dsd matrix
    - output: clusters
    - exec: python philharmonic_clustering.py
- add recipe proteins back to clusters
    - input: clusters
    - output: clusters with recipe proteins
    - exec: recipe clusters.json
- add pfam domains and go terms to clusters
    - input: clusters(after recipe), protein-pfam associations
    - other requires: pfam go associations
    - output: functionally annotated clusters
    - exec: python philharmonic_functional_annotation.py
- name and describe clusters with langchain
    - input: functionally annotated clusters
    - output: human readable clusters
    - exec: python langchain.py
- create cytoscape file
    - input: functionally annotated clusters, human readable clusters
    - other requires: styles.xml
    - output: cytoscape session
    - exec: python philharmonic_cytoscape.py

## Environment

```
mamba create -f environment.yml
mamba activate philharmonic
```
