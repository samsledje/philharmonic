from pathlib import Path
import os

# Define base paths for cleaner syntax
work_dir = Path(config['work_dir'])
run_name = config['run_name']

rule PHILHARMONIC:
    input:
        network = work_dir / f"{run_name}_network.positive.tsv",
        clusters = work_dir / f"{run_name}_clusters.json",
        cluster_graph = work_dir / f"{run_name}_cluster_graph.tsv",
        cluster_functions = work_dir / f"{run_name}_cluster_graph_functions.tsv",
        human_readable = work_dir / f"{run_name}_human_readable.txt",
        go_map = work_dir / f"{run_name}_GO_map.csv",
        protein_csv = work_dir / f"{run_name}_protein_annotations.csv"
    log:
        work_dir / "logs" / "PHILHARMONIC.log",
    output:
        zipfile = work_dir / f"{run_name}.zip"
    shell:
        "zip --junk-paths {output.zipfile} {input.network} {input.clusters} {input.go_map} {input.human_readable} {input.cluster_graph} {input.cluster_functions} {input.protein_csv} > {log} 2>&1"


rule download_required_files:
    output:
        go_database = work_dir / "go.obo",
        go_slim = work_dir / "goslim_generic.obo",
        pfam_database_zipped = work_dir / "Pfam-A.hmm.gz",
        pfam_gobp = work_dir / "pfam_gobp_most_specific.txt",
    log:
        work_dir / "logs" / "01_download_required_files.log",
    shell:
        """
        mkdir -p {work_dir} && \
        curl https://current.geneontology.org/ontology/go.obo -o {output.go_database} && \
        curl https://current.geneontology.org/ontology/subsets/goslim_generic.obo -o {output.go_slim} && \
        curl https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz -o {output.pfam_database_zipped} && \
        curl https://godm.loria.fr/data/pfam_gobp_most_specific.txt -o {output.pfam_gobp} \
        > {log} 2>&1
        """


rule prepare_hmmdb:
    input:
        pfam_database_zipped = work_dir / "Pfam-A.hmm.gz",
    output:
        pfam_database = work_dir / "Pfam-A.hmm",
        pfam_h3m = temp(work_dir / "Pfam-A.hmm.h3m"),
        pfam_h3i = temp(work_dir / "Pfam-A.hmm.h3i"),
        pfam_h3f = temp(work_dir / "Pfam-A.hmm.h3f"),
        pfam_h3p = temp(work_dir / "Pfam-A.hmm.h3p"),
    log:
        work_dir / "logs" / "02_prepare_hmmdb.log",
    shell:
        """
        gunzip -c {input.pfam_database_zipped} > {output.pfam_database} 2> {log} && \
        hmmpress {output.pfam_database} >> {log} 2>&1
        """


rule annotate_seqs_pfam:
    input:
        sequences = config["sequence_path"],
        pfam_database = work_dir / "Pfam-A.hmm",
        pfam_h3m = work_dir / "Pfam-A.hmm.h3m",
        pfam_h3i = work_dir / "Pfam-A.hmm.h3i",
        pfam_h3f = work_dir / "Pfam-A.hmm.h3f",
        pfam_h3p = work_dir / "Pfam-A.hmm.h3p",
    output:
        pfam_map = work_dir / f"{run_name}_hmmscan.tblout",
        hmmscan_out = work_dir / f"{run_name}_hmmscan.out",
        domtblout = work_dir / f"{run_name}_hmmscan.domtblout",
    threads: config["hmmscan"]["threads"]
    params:
        hmmscan_path = config["hmmscan"]["path"]
    log:
        work_dir / "logs" / "03_annotate_seqs_pfam.log",
    shell:
        "{params.hmmscan_path} --cpu {threads} -o {output.hmmscan_out} --tblout {output.pfam_map} --domtblout {output.domtblout} --acc --noali --notextw --cut_ga {input.pfam_database} {input.sequences} > {log} 2>&1"

rule annotate_seqs_go:
    input:
        hhtblout = work_dir / f"{run_name}_hmmscan.tblout",
        pfam_gobp = work_dir / "pfam_gobp_most_specific.txt",
    output:
        go_map = work_dir / f"{run_name}_GO_map.csv",
    log:
        work_dir / "logs" / "04_annotate_seqs_go.log",
    shell:
        "philharmonic build-go-map -o {output.go_map} --hhtblout {input.hhtblout} --pfam_go_files {input.pfam_gobp} > {log} 2>&1"


rule generate_candidates:
    input:
        sequences = config["sequence_path"],
        go_database = work_dir / "go.obo",
        go_map = work_dir / f"{run_name}_GO_map.csv",
    output:
        kept_proteins = work_dir / f"{run_name}_kept_proteins.txt",
        candidates = work_dir / f"{run_name}_candidates.tsv",
    params:
        n_pairs = config["dscript"]["n_pairs"],
        seed = config["seed"],
        go_filter = f"--go_filter={config['go_filter_path']}" if config.get('go_filter_path') else ""
    log:
        work_dir / "logs" / "05_generate_candidates.log",
    shell:
        "philharmonic generate-candidates --paircount {params.n_pairs} -o {output.candidates} --seq_out {output.kept_proteins} --go_map {input.go_map} --go_database {input.go_database} {params.go_filter} --sequences {input.sequences} --seed {params.seed} > {log} 2>&1"

rule predict_network:
    input:
        sequences = config["sequence_path"],
        candidates = work_dir / f"{run_name}_candidates.tsv",
    output:
        network = work_dir / f"{run_name}_network.positive.tsv",
    params:
        dscript_path = config["dscript"]["path"],
        dscript_model = config['dscript']['model'],
        device = config["dscript"]["device"],
        thresh = config["dsd"]["t"],
        outfile_prefix = work_dir / f"{run_name}_network",
    resources:
        nvidia_gpu=1
    log:
        work_dir / "logs" / "06_predict_network.log",
    shell:
        "{params.dscript_path} predict --pairs {input.candidates} --seqs {input.sequences} --model {params.dscript_model} --outfile {params.outfile_prefix} --device {params.device} --thresh {params.thresh} > {log}"


rule compute_distances:
    input:
        network = work_dir / f"{run_name}_network.positive.tsv",
    output:
        distances = work_dir / f"{run_name}_distances.DSD1",
    params:
        dsd_path = config["dsd"]["path"],
        t = config["dsd"]["t"],
        confidence = "-c" if config["dsd"]["confidence"] else "",
        outfile_prefix = work_dir / f"{run_name}_distances",
    log:
        work_dir / "logs" / "07_compute_distances.log",
    shell:
        "{params.dsd_path} {params.confidence} --converge -t {params.t} --outfile {params.outfile_prefix} {input.network} > {log} 2>&1"


rule cluster_network:
    input:
        network = work_dir / f"{run_name}_network.positive.tsv",
        distances = work_dir / f"{run_name}_distances.DSD1",
    output:
        clusters = temp(work_dir / f"{run_name}_clusters.disconnected.json"),
    params:
        min_cluster_size = config["clustering"]["min_cluster_size"],
        cluster_divisor = config["clustering"]["cluster_divisor"],
        init_k = config["clustering"]["init_k"],
        sparsity = config["clustering"]["sparsity_thresh"],
        seed = config["seed"]
    log:
        work_dir / "logs" / "08_cluster_network.log",
    shell:
        "philharmonic cluster-network --network_file {input.network} --dsd_file {input.distances} --output {output.clusters} --min_cluster_size {params.min_cluster_size} --cluster_divisor {params.cluster_divisor} --init_k {params.init_k} --sparsity {params.sparsity} --random_seed {params.seed} > {log} 2>&1"


rule reconnect_recipe:
    input:
        clusters = work_dir / f"{run_name}_clusters.disconnected.json",
        network = work_dir / f"{run_name}_network.positive.tsv",
    output:
        clusters_connected = temp(work_dir / f"{run_name}_clusters.recipe.json"),
    params:
        lr = config["recipe"]["lr"],
        cthresh = config["recipe"]["cthresh"],
        max_proteins = config["recipe"]["max_proteins"],
        metric = config["recipe"]["metric"]
    log:
        work_dir / "logs" / "09_reconnect_recipe.log",
    shell:
        "recipe-cluster cook --network-filepath {input.network} --cluster-filepath {input.clusters} --lr {params.lr} -cthresh {params.cthresh} --max {params.max_proteins} --metric {params.metric} --outfile {output.clusters_connected} > {log}"

rule add_cluster_functions:
    input:
        clusters = work_dir / f"{run_name}_clusters.recipe.json",
        go_map = work_dir / f"{run_name}_GO_map.csv",
    output:
        clusters_functional = work_dir / f"{run_name}_clusters.functional.json",
    log:
        work_dir / "logs" / "10_add_cluster_functions.log",
    shell:
        "philharmonic add-cluster-functions -o {output.clusters_functional} -cfp {input.clusters} --go-map {input.go_map} > {log} 2>&1"

rule cluster_graph:
    input:
        clusters = work_dir / f"{run_name}_clusters.functional.json",
        network = work_dir / f"{run_name}_network.positive.tsv",
        go_map = work_dir / f"{run_name}_GO_map.csv",
        go_database = work_dir / "goslim_generic.obo",
    output:
        graph = work_dir / f"{run_name}_cluster_graph.tsv",
        coc_functions = work_dir / f"{run_name}_cluster_graph_functions.tsv",
    params:
        recipe_metric = config["recipe"]["metric"],
        recipe_cthresh = config["recipe"]["cthresh"],
    log:
        work_dir / "logs" / "11_cluster_graph.log",
    shell:
        "philharmonic build-cluster-graph -o {output.graph} -coc {output.coc_functions} -cfp {input.clusters} -nfp {input.network} --go_map {input.go_map} --go_db {input.go_database} --recipe_metric {params.recipe_metric} --recipe_cthresh {params.recipe_cthresh} > {log} 2>&1"

rule summarize_clusters:
    input:
        clusters = work_dir / f"{run_name}_clusters.functional.json",
        go_database = work_dir / "go.obo",
    output:
        human_readable = work_dir / f"{run_name}_human_readable.txt",
        readable_json = work_dir / f"{run_name}_clusters.json",
    params:
        api_key = f"--api-key {os.environ['OPENAI_API_KEY']}" if config["use_llm"] else "",
        llm_model = config["llm"]["model"],
        do_llm_naming = "--llm-name" if config["use_llm"] else ""
    log:
        work_dir / "logs" / "12_summarize_clusters.log",
    shell:
        "philharmonic summarize-clusters {params.do_llm_naming} --model {params.llm_model} {params.api_key} -o {output.human_readable} --json {output.readable_json} --go_db {input.go_database} -cfp {input.clusters} > {log}"

rule write_protein_table:
    input:
        sequences = config["sequence_path"],
        readable_json = work_dir / f"{run_name}_clusters.json",
        go_map = work_dir / f"{run_name}_GO_map.csv",
        gff=config["gff"]["path"]
    output:
        protein_csv = work_dir / f"{run_name}_protein_annotations.csv"
    params:
        species=config["taxonomy"]["species"],
        genus=config["taxonomy"]["genus"],
        family=config["taxonomy"]["family"],
        order=config["taxonomy"]["order"],
        class_name=config["taxonomy"]["class"],
        phylum=config["taxonomy"]["phylum"],
        domain=config["taxonomy"]["domain"],
        gff_feature=config["gff"].get("feature_type", "CDS"),
        gff_attr=config["gff"].get("attribute_key", "protein_id"),
        gff_id_preface=config["gff"].get("id_preface", "")
    log:
        work_dir / "logs" / "13_json_to_csv.log"
    shell:
        """
        philharmonic write-protein-table \
            {input.readable_json} \
            {input.go_map} \
            {input.sequences} \
            {input.gff} \
            {output.protein_csv} \
            --species {params.species:q} \
            --genus {params.genus:q} \
            --family {params.family:q} \
            --order {params.order:q} \
            --class {params.class_name:q} \
            --phylum {params.phylum:q} \
            --domain {params.domain:q} \
            --gff-feature-type {params.gff_feature:q} \
            --gff-attribute-key {params.gff_attr:q} \
            --gff-id-preface "{params.gff_id_preface:q}" \
            > {log} 2>&1
        """
