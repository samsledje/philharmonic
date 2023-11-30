# philharmonic
Pipeline for decoding the functional networks of non-model organisms

Requirements:
    Pfam database

Required Inputs: protein sequences
Optional Inputs: shortlist proteins, shortlist GO terms

Steps:
- (optional: filter candidates with experimental requirements)
- generate list of candidate interactions
- predict interactions with dscript
- compute dsd distances
- cluster dsd matrix
- reconnect clusters with recipe
- functional annotation of clusters
