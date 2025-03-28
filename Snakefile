rule PHILHARMONIC:
    input:
        network = f"{config['work_dir']}/{config['run_name']}_network.positive.tsv",
        clusters = f"{config['work_dir']}/{config['run_name']}_clusters.json",
        cluster_graph = f"{config['work_dir']}/{config['run_name']}_cluster_graph.tsv",
        cluster_functions = f"{config['work_dir']}/{config['run_name']}_cluster_graph_functions.tsv",
        human_readable = f"{config['work_dir']}/{config['run_name']}_human_readable.txt",
        go_map = f"{config['work_dir']}/{config['run_name']}_GO_map.csv",
    log:
        f"{config['work_dir']}/logs/PHILHARMONIC.log",
    output:
        zipfile = f"{config['work_dir']}/{config['run_name']}.zip"
    params:
    shell:
        "zip --junk-paths {output.zipfile} {input.network} {input.clusters} {input.go_map} {input.human_readable} {input.cluster_graph} {input.cluster_functions} > {log} 2>&1"


rule download_required_files:
    output:
        go_database = f"{config['work_dir']}/go.obo",
        go_slim = f"{config['work_dir']}/goslim_generic.obo",
        pfam_database_zipped = f"{config['work_dir']}/Pfam-A.hmm.gz",
        pfam_gobp = f"{config['work_dir']}/pfam_gobp_most_specific.txt",
    log:
        f"{config['work_dir']}/logs/01_download_required_files.log",
    run:
        commands = [
            "mkdir -p {config[work_dir]} > {log} 2>&1",
            "curl https://current.geneontology.org/ontology/go.obo -o {config[work_dir]}/go.obo > {log} 2>&1",
            "curl https://current.geneontology.org/ontology/subsets/goslim_generic.obo -o {config[work_dir]}/goslim_generic.obo > {log} 2>&1",
            "curl ftp://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz -o {config[work_dir]}/Pfam-A.hmm.gz > {log} 2>&1",
            "curl https://godm.loria.fr/data/pfam_gobp_most_specific.txt -o {config[work_dir]}/pfam_gobp_most_specific.txt > {log} 2>&1",
        ]
        for c in commands:
            shell(c)


rule prepare_hmmdb:
    input:
        pfam_database_zipped = f"{config['work_dir']}/Pfam-A.hmm.gz",
    output:
        pfam_database = f"{config['work_dir']}/Pfam-A.hmm",
        pfam_h3m = temp(f"{config['work_dir']}/Pfam-A.hmm.h3m"),
        pfam_h3i = temp(f"{config['work_dir']}/Pfam-A.hmm.h3i"),
        pfam_h3f = temp(f"{config['work_dir']}/Pfam-A.hmm.h3f"),
        pfam_h3p = temp(f"{config['work_dir']}/Pfam-A.hmm.h3p"),
    log:
        f"{config['work_dir']}/logs/02_prepare_hmmdb.log",
    run:
        commands = [
            "gunzip -c {config[work_dir]}/Pfam-A.hmm.gz > {config[work_dir]}/Pfam-A.hmm 2> {log}",
            "hmmpress {config[work_dir]}/Pfam-A.hmm > {log} 2>&1",
        ]
        for c in commands:
            shell(c)


rule annotate_seqs_pfam:
    input:
        sequences = f"{config['sequence_path']}",
        pfam_database = f"{config['work_dir']}/Pfam-A.hmm",
        pfam_h3m = f"{config['work_dir']}/Pfam-A.hmm.h3m",
        pfam_h3i = f"{config['work_dir']}/Pfam-A.hmm.h3i",
        pfam_h3f = f"{config['work_dir']}/Pfam-A.hmm.h3f",
        pfam_h3p = f"{config['work_dir']}/Pfam-A.hmm.h3p",
    output:
        pfam_map = f"{config['work_dir']}/{config['run_name']}_hmmscan.tblout",
    threads: config["hmmscan"]["threads"]
    params:
        work_dir = config["work_dir"],
        run_name = config["run_name"],
        hmmscan_path = config["hmmscan"]["path"]
    log:
        f"{config['work_dir']}/logs/03_annotate_seqs_pfam.log",
    shell:
        "{params.hmmscan_path} --cpu {threads} -o {params.work_dir}/{params.run_name}_hmmscan.out --tblout {params.work_dir}/{params.run_name}_hmmscan.tblout --domtblout {params.work_dir}/{params.run_name}_hmmscan.domtblout --acc --noali --notextw --cut_ga {input.pfam_database} {input.sequences} > {log} 2>&1"

rule annotate_seqs_go:
    input:
        hhtblout = f"{config['work_dir']}/{config['run_name']}_hmmscan.tblout",
        pfam_gobp = f"{config['work_dir']}/pfam_gobp_most_specific.txt",
    output:
        go_map = f"{config['work_dir']}/{config['run_name']}_GO_map.csv",
    log:
        f"{config['work_dir']}/logs/04_annotate_seqs_go.log",
    shell:
        "philharmonic build-go-map -o {output.go_map} --hhtblout {input.hhtblout} --pfam_go_files {input.pfam_gobp} > {log} 2>&1"


rule generate_candidates:
    input:
        sequences = f"{config['sequence_path']}",
        go_database = f"{config['work_dir']}/go.obo",
        go_map = f"{config['work_dir']}/{config['run_name']}_GO_map.csv",
    output:
        kept_proteins = f"{config['work_dir']}/{config['run_name']}_kept_proteins.txt",
        candidates = f"{config['work_dir']}/{config['run_name']}_candidates.tsv",
    params:
        n_pairs = config["dscript"]["n_pairs"],
        seed = config["seed"],
        go_filter = lambda wildcards, input: f"--go_filter={config['go_filter_path']}" if 'go_filter_path' in config else ""
    log:
        f"{config['work_dir']}/logs/05_generate_candidates.log",
    shell:
        "philharmonic generate-candidates --paircount {params.n_pairs} -o {output.candidates} --seq_out {output.kept_proteins} --go_map {input.go_map} --go_database {input.go_database} {params.go_filter} --sequences {input.sequences} --seed {params.seed} > {log} 2>&1"

rule predict_network:
    input:
        sequences = f"{config['sequence_path']}",
        candidates = f"{config['work_dir']}/{config['run_name']}_candidates.tsv",
    output:
        network = f"{config['work_dir']}/{config['run_name']}_network.positive.tsv",
    params:
        dscript_path = config["dscript"]["path"],
        dscript_model = config['dscript']['model'],
        work_dir = config["work_dir"],
        run_name = config["run_name"],
        device = config["dscript"]["device"],
        t = config["dsd"]["t"],
    resources:
        nvidia_gpu=1
    log:
        f"{config['work_dir']}/logs/06_predict_network.log",
    shell:
        "{params.dscript_path} predict --pairs {input.candidates} --seqs {input.sequences} --model {params.dscript_model} --outfile {params.work_dir}/{params.run_name}_network --device {params.device} --thresh {params.t} > {log}"


rule compute_distances:
    input:
        network = f"{config['work_dir']}/{config['run_name']}_network.positive.tsv",
    output:
        distances = f"{config['work_dir']}/{config['run_name']}_distances.DSD1",
    params:
        dsd_path = config["dsd"]["path"],
        work_dir = config["work_dir"],
        run_name = config["run_name"],
        t = config["dsd"]["t"],
        confidence = "-c" if config["dsd"]["confidence"] else ""
    log:
        f"{config['work_dir']}/logs/07_compute_distances.log",
    shell:
        "{params.dsd_path} {params.confidence} --converge -t {params.t} --outfile {params.work_dir}/{params.run_name}_distances {input.network} > {log} 2>&1"


rule cluster_network:
    input:
        network = f"{config['work_dir']}/{config['run_name']}_network.positive.tsv",
        distances = f"{config['work_dir']}/{config['run_name']}_distances.DSD1",
    output:
        clusters = temp(f"{config['work_dir']}/{config['run_name']}_clusters.disconnected.json"),
    params:
        work_dir = config["work_dir"],
        run_name = config["run_name"],
        min_cluster_size = config["clustering"]["min_cluster_size"],
        cluster_divisor = config["clustering"]["cluster_divisor"],
        init_k = config["clustering"]["init_k"],
        sparsity = config["clustering"]["sparsity_thresh"],
        seed = config["seed"]
    log:
        f"{config['work_dir']}/logs/08_cluster_network.log",
    shell:
        "philharmonic cluster-network --network_file {input.network} --dsd_file {input.distances} --output {output.clusters} --min_cluster_size {params.min_cluster_size} --cluster_divisor {params.cluster_divisor} --init_k {params.init_k} --sparsity {params.sparsity} --random_seed {params.seed} > {log} 2>&1"


rule reconnect_recipe:
    input:
        clusters = f"{config['work_dir']}/{config['run_name']}_clusters.disconnected.json",
        network = f"{config['work_dir']}/{config['run_name']}_network.positive.tsv",
    output:
        clusters_connected = temp(f"{config['work_dir']}/{config['run_name']}_clusters.recipe.json"),
    params:
        lr = config["recipe"]["lr"],
        cthresh = config["recipe"]["cthresh"],
        max_proteins = config["recipe"]["max_proteins"],
        metric = config["recipe"]["metric"]
    log:
        f"{config['work_dir']}/logs/09_reconnect_recipe.log",
    shell:
        "recipe-cluster cook --network-filepath {input.network} --cluster-filepath {input.clusters} --lr {params.lr} -cthresh {params.cthresh} --max {params.max_proteins} --metric {params.metric} --outfile {output.clusters_connected} > {log}"

rule add_cluster_functions:
    input:
        clusters = f"{config['work_dir']}/{config['run_name']}_clusters.recipe.json",
        go_map = f"{config['work_dir']}/{config['run_name']}_GO_map.csv",
    output:
        clusters_functional = f"{config['work_dir']}/{config['run_name']}_clusters.functional.json",
    log:
        f"{config['work_dir']}/logs/10_add_cluster_functions.log",
    shell:
        "philharmonic add-cluster-functions -o {output.clusters_functional} -cfp {input.clusters} --go-map {input.go_map} > {log} 2>&1"

rule cluster_graph:
    input:
        clusters = f"{config['work_dir']}/{config['run_name']}_clusters.functional.json",
        network = f"{config['work_dir']}/{config['run_name']}_network.positive.tsv",
        go_map = f"{config['work_dir']}/{config['run_name']}_GO_map.csv",
        go_database = f"{config['work_dir']}/goslim_generic.obo",
    output:
        graph = f"{config['work_dir']}/{config['run_name']}_cluster_graph.tsv",
        coc_functions = f"{config['work_dir']}/{config['run_name']}_cluster_graph_functions.tsv",
    params:
        recipe_metric = config["recipe"]["metric"],
        recipe_cthresh = config["recipe"]["cthresh"],
    log:
        f"{config['work_dir']}/logs/11_cluster_graph.log",
    shell:
        "philharmonic build-cluster-graph -o {output.graph} -coc {output.coc_functions} -cfp {input.clusters} -nfp {input.network} --go_map {input.go_map} --go_db {input.go_database} --recipe_metric {params.recipe_metric} --recipe_cthresh {params.recipe_cthresh} > {log} 2>&1"

rule summarize_clusters:
    input:
        clusters = f"{config['work_dir']}/{config['run_name']}_clusters.functional.json",
        go_database = f"{config['work_dir']}/go.obo",
    output:
        human_readable = f"{config['work_dir']}/{config['run_name']}_human_readable.txt",
        readable_json = f"{config['work_dir']}/{config['run_name']}_clusters.json",
    params:
        api_key = f"--api-key {os.environ['OPENAI_API_KEY']}" if config["use_llm"] else "",
        llm_model = config["llm"]["model"],
        do_llm_naming = "--llm-name" if config["use_llm"] else ""
    log:
        f"{config['work_dir']}/logs/12_summarize_clusters.log",
    shell:
        "philharmonic summarize-clusters {params.do_llm_naming} --model {params.llm_model} {params.api_key} -o {output.human_readable} --json {output.readable_json} --go_db {input.go_database} -cfp {input.clusters} > {log}"
