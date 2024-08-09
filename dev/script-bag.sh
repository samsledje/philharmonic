

#download GO data
function download_GO_data {
    DIR=/afs/csail.mit.edu/u/r/rsingh/work/corals/data/raw/go
    mkdir -p $DIR
    cd $DIR

    wget https://current.geneontology.org/ontology/go.obo
    wget http://godm.loria.fr/data/pfam_gomf_most_specific.txt
    wget http://godm.loria.fr/data/pfam_gobp_most_specific.txt
    wget http://godm.loria.fr/data/pfam_gocc_most_specific.txt
}


#get the list of top and 2nd-level GO terms
function process_GO_toplevel {
    DIR=/afs/csail.mit.edu/u/r/rsingh/work/corals/data/raw/go
    cd $DIR

    grep '^is_a' go.obo | egrep 'molecular_function|biological_process|cellular_component' | sort | uniq |
	awk 'BEGIN {print "GOID,desc"} {print $2 "," $4}' | perl -pe 's/, /,/' > go_level1.csv
    
    sed 1,28d go.obo | awk '/^id:/ {id=$2;} /^name:/ {$1=""; nm=$0} /^is_a: .*(biological_process|molecular_function|cellular_component)/ {print id "," nm; id=nm=""}' |\
	sort | uniq | awk 'BEGIN {print "GOID,desc"} {print} '| perl -pe 's/, /,/'  > go_level2.csv
    
    cat go_level2.csv | awk -F, 'NR==1 {print "GOID,desc,useflag"} NR>1 { if ($2 ~ /(biological regulation)|(cellular anatomical entity)|(cellular process)|(developmental process)|(growth)|(intracellular)|(localization)|(metabolic process)|(other organism part)|(protein-containing complex)|(reproduction)|(reproductive process)|(transcription regulator activity)|(translation regulator activity)/) {print $1 "," $2 "," 0} else {print $1 "," $2 "," 1}} ' > go_level2_marked-up.csv

}


function run_pfam_on_pdam {
    DIR=/afs/csail.mit.edu/u/r/rsingh/work/corals/data/
    cd $DIR
    SEQDIR=${DIR}/raw/seqs
    ODIR=${DIR}/processed
    
    for x in 0:7000 7000:14000 14000:21000 21000:28000
    do
	a=$(echo $x | cut -d: -f1); b=$(echo $x | cut -d: -f2);
	sfx=chunk-${a}-${b}
	outf=${ODIR}/pdam_proteins_${sfx}.fasta.gz
	
	zcat ${SEQDIR}/pdam_proteins.fasta.gz | awk '/^>/ {seqcnt++} seqcnt>='$a' && seqcnt < '$b' {print}' | gzip -c > $outf

	hmmoutf=${ODIR}/pdam_hmmscan_${sfx}
	/scratch1/rsingh/tools/bin/hmmscan -o ${hmmoutf}.out --tblout ${hmmoutf}.tblout --domtblout ${hmmoutf}.domtblout --acc --noali --notextw --cut_ga ${DIR}/raw/pfam/Pfam-A.hmm $outf &
	
    done
    wait

}



function run_hmmscan_on_hvul {
    DIR=/afs/csail.mit.edu/u/r/rsingh/work/corals/data/
    cd $DIR
    SEQDIR=${DIR}/raw/seqs
    ODIR=${DIR}/processed
    
    for x in 0:4200 4200:8400 8400:12600 12600:19000
    do
	a=$(echo $x | cut -d: -f1); b=$(echo $x | cut -d: -f2);
	sfx=chunk-${a}-${b}
	outf=${ODIR}/hvul_proteins_${sfx}.fasta.gz
	
	cat ${SEQDIR}/hvul_proteins.fasta | awk '/^>/ {seqcnt++} seqcnt>='$a' && seqcnt < '$b' {print}' | gzip -c > $outf

	hmmoutf=${ODIR}/hvul_hmmscan_${sfx}
	/scratch1/rsingh/tools/bin/hmmscan -o ${hmmoutf}.out --tblout ${hmmoutf}.tblout --domtblout ${hmmoutf}.domtblout --acc --noali --notextw --cut_ga ${DIR}/raw/pfam/Pfam-A.hmm $outf &
	
    done
    wait

}


function extract_fasta_subset() {
    orig_fasta=$1
    seqid_list_csv=$2 
    seqid_colidx=${3:-1} #assumes there's a header and that seqid is in column 1
    
    zcat -f --stdout $orig_fasta |\
    awk 'BEGIN {nm=""} \
         /^>/ {if (nm!="") {seq[nm]=block;} nm=substr($1,2); block=$0} \
         !/^>/ {block = block "\n" $0} \
         END { seq[nm]=block; k='$seqid_colidx'; while ("sed 1d '$seqid_list_csv' " | getline) { split($0,a,","); if (a[k] in seq) {print seq[a[k]]} }}' 
}



function run_hmmscan_on_tpseud {
    DIR=/afs/csail.mit.edu/u/r/rsingh/work/corals/data/
    cd $DIR
    SEQDIR=${DIR}/raw/seqs
    ODIR=${DIR}/processed
    
    for x in 0:4500 4500:10000
    do
	a=$(echo $x | cut -d: -f1); b=$(echo $x | cut -d: -f2);
	sfx=chunk-${a}-${b}
	outf=${ODIR}/tpseud_proteins_${sfx}.fasta.gz
	
	cat ${SEQDIR}/t_pseudonana.fa | awk '/^>/ {seqcnt++} seqcnt>='$a' && seqcnt < '$b' {print}' | gzip -c > $outf

	hmmoutf=${ODIR}/tpseud_hmmscan_${sfx}
	/scratch1/rsingh/tools/bin/hmmscan -o ${hmmoutf}.out --tblout ${hmmoutf}.tblout --domtblout ${hmmoutf}.domtblout --acc --noali --notextw --cut_ga ${DIR}/raw/pfam/Pfam-A.hmm $outf &
	
    done
    wait

}



function run_hmmscan_on_symbc1 {
    DIR=/afs/csail.mit.edu/u/r/rsingh/work/corals/data/
    cd $DIR
    SEQDIR=${DIR}/raw/seqs
    ODIR=${DIR}/processed
    
    for x in 0:7000 7000:14000 14000:21000 21000:30000
    do
	a=$(echo $x | cut -d: -f1); b=$(echo $x | cut -d: -f2);
	sfx=chunk-${a}-${b}
	outf=${ODIR}/symbc1_proteins_${sfx}.fasta.gz
	
	cat ${SEQDIR}/symbc1.fa | awk '/^>/ {seqcnt++} seqcnt>='$a' && seqcnt < '$b' {print}' | gzip -c > $outf

	hmmoutf=${ODIR}/symbc1_hmmscan_${sfx}
	/scratch1/rsingh/tools/bin/hmmscan -o ${hmmoutf}.out --tblout ${hmmoutf}.tblout --domtblout ${hmmoutf}.domtblout --acc --noali --notextw --cut_ga ${DIR}/raw/pfam/Pfam-A.hmm $outf &
	
    done
    wait

}

#/afs/csail.mit.edu/u/r/rsingh/work/corals/data/raw/Bos_Taursus_longest_isoforms_protein_coding.fa

function run_hmmscan_on_cow {
    DIR=/afs/csail.mit.edu/u/r/rsingh/work/corals/data/
    cd $DIR
    SEQDIR=${DIR}/raw/
    ODIR=${DIR}/processed
    
    for x in 0:5000 5000:10000 10000:15000 15000:30000
    do
	a=$(echo $x | cut -d: -f1); b=$(echo $x | cut -d: -f2);
	sfx=chunk-${a}-${b}
	outf=${ODIR}/cow_proteins_${sfx}.fasta.gz
	
	cat ${SEQDIR}/Bos_Taursus_longest_isoforms_protein_coding.fa | awk '/^>/ {seqcnt++} seqcnt>='$a' && seqcnt < '$b' {print}' | gzip -c > $outf

	hmmoutf=${ODIR}/cow_hmmscan_${sfx}
	/scratch1/rsingh/tools/bin/hmmscan -o ${hmmoutf}.out --tblout ${hmmoutf}.tblout --domtblout ${hmmoutf}.domtblout --acc --noali --notextw --cut_ga ${DIR}/raw/pfam/Pfam-A.hmm $outf &
	
    done
    wait

}



# ### Python code for cross-species interaction list
#
# d1 = pd.read_csv("proteins_list_tpseud__16Mpairs_20200830.csv")
# spat = re.compile(r'^.*((transmembrane)|(ion channel)|(ion transport)|(transporter activity)).*$')
# d2 = pd.read_csv("proteins_list_symbc1_50Mpairs_20200830.csv")
# idx1 = d1.GO_list.apply(lambda s: spat.search(s) is not None).values
# idx2 = d2.GO_list.apply(lambda s: spat.search(s) is not None).values
# d1a = d1.iloc[idx1,:].copy().reset_index(drop=True)
# d2a = d2.iloc[idx2,:].copy().reset_index(drop=True)
#
# nPointPairs=5000000
# j_u = np.random.choice(range(d1a.shape[0]), int(nPointPairs))
# j_v = np.random.choice(range(d2.shape[0]), int(nPointPairs))
# df3a = pd.DataFrame( zip(d1a.species_id.values[j_u], d2.species_id.values[j_v]), columns=["protein_A","protein_B"])
#
# nPointPairs=7500000
# j_u = np.random.choice(range(d1.shape[0]), int(nPointPairs))
# j_v = np.random.choice(range(d2a.shape[0]), int(nPointPairs))
# df3b = pd.DataFrame( zip(d1.species_id.values[j_u], d2a.species_id.values[j_v]), columns=["protein_A","protein_B"])
#
# df3 = pd.concat([df3a, df3b])
# df3c = df3.drop_duplicates().reset_index(drop=True)
# df3c.to_csv("/afs/csail.mit.edu/u/r/rsingh/work/corals/data/processed/candidate_pairs_X_tpseud-symbc1_10Mpairs_20200830.csv", index=False)
#
# ####


# ### Python code for reading Cow genes
#
# id2fa = {}
# for line in open("/afs/csail.mit.edu/u/r/rsingh/work/corals/data/raw/GCF_002263795.1_ARS-UCD1.2_protein.faa", "rt"):
#     if line[0]==">":
#         if newseq!="": id2fa[refseq]=newseq
#         refseq=line[1:].split()[0]
#         newseq=line
#     else:
#         newseq += line
#     id2fa[refseq] = newseq

# df_gff = pd.read_csv("/afs/csail.mit.edu/u/r/rsingh/work/corals/data/raw/ARS-UCD1.2_RefSeq_all_proteincoding.gff3", sep="\t", comment="#", header=None)

# from collections import defaultdict
# gene2proteins = defaultdict(list); valid_protein_ids = set()
# for i in range(df_gff.shape[0]):
#     if df_gff.iat[i,2]!="mRNA": continue
#     s = df_gff.iat[i,8]
#     l = [a for a in s.split(';') if a.startswith("Dbxref=")]
#     if not l: continue
#     l2 = [a for a in l[0].split(',') if a.startswith("RefSeq_Prot:")]
#     if not l2: continue
#     pc = l2[0].split(':')[1]
#     l3 = [a for a in s.split(';') if a.startswith("symbol_ncbi=")]
#     if not l3: continue
#     gc = l3[0].split("=")[1]
#     valid_protein_ids.add(pc)
#     gene2proteins[gc].append(pc)

# gene2fa = {}
# for g, plist in gene2proteins.items():
#     l2 = [(f_protlen(id2fa[p]), p) for p in plist if p in id2fa]
#     if not l2: continue
#     _, p = max(l2)
#     hdr = id2fa[p].split("\n")[0]
#     hdr += " symbol_ncbi={}".format(g)
#     s2 = "\n".join([hdr] + id2fa[p].split("\n")[1:])
#     gene2fa[g] = s2

# list(gene2fa.keys())
# gene2fa['RNASE4']
# outfh=open('/afs/csail.mit.edu/u/r/rsingh/work/corals/data/raw/Bos_Taursus_longest_isoforms_protein_coding.fa','wt')
# for fa in gene2fa.values()
# for fa in gene2fa.values():
#     outfh.write(fa)
# outfh.close()

    

cd honig-lab/ecoli-recomb-retest
dscript embed --seq ecoli.fasta -o ecoli_EMBEDS.h5 -d 0
