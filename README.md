# philharmonic
Pipeline for decoding the functional networks of non-model organisms

Required Inputs: protein sequences
Optional Inputs: shortlist proteins, shortlist GO terms

Steps:
- download go database
- filter go database to terms we care about
- map go terms to pfam domains (godomain miner)
- map pfam domains to protein sequences (hmmscan)
- (optional: filter candidates with experimental requirements)
- generate list of candidate interactions
- predict interactions with dscript
- compute dsd distances
- cluster dsd matrix
- reconnect clusters with recipe
- functional annotation of clusters
