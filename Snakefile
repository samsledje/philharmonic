rule PHILHARMONIC:
    input:
        network = f"{config['work_dir']}/{config['run_name']}_network.positive.tsv",
        clusters = f"{config['work_dir']}/{config['run_name']}_clusters.json",
        cluster_graph = f"{config['work_dir']}/{config['run_name']}_cluster_graph.json",
        cytoscape = f"{config['work_dir']}/{config['run_name']}_cytoscape_session.cys" if config["use_cytoscape"] else [],
        human_readable = f"{config['work_dir']}/{config['run_name']}_human_readable.txt",
    output:
        zipfile = f"{config['work_dir']}/{config['run_name']}.zip"
    params:
    shell:
        "zip {output.zipfile} {input.network} {input.clusters} {input.cluster_graph} {input.human_readable} {input.cytoscape}"


rule download_required_files:
    output:
        go_database = f"{config['work_dir']}/go.obo",
        go_slim = f"{config['work_dir']}/goslim_generic.obo",
        pfam_database_zipped = f"{config['work_dir']}/Pfam-A.hmm.gz",
        pfam_gobp = f"{config['work_dir']}/pfam_gobp_most_specific.txt",
    log:
        "logs/download_required_files.log",
    run:
        commands = [
            "mkdir -p {config[work_dir]}",
            "curl https://current.geneontology.org/ontology/go.obo -o {config[work_dir]}/go.obo",
            "curl https://current.geneontology.org/ontology/subsets/goslim_generic.obo -o {config[work_dir]}/goslim_generic.obo",
            "curl ftp://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz -o {config[work_dir]}/Pfam-A.hmm.gz",
            "curl https://godm.loria.fr/data/pfam_gobp_most_specific.txt -o {config[work_dir]}/pfam_gobp_most_specific.txt",
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
        "logs/prepare_hmmdb.log",
    run:
        commands = [
            "gunzip -c {config[work_dir]}/Pfam-A.hmm.gz > {config[work_dir]}/Pfam-A.hmm",
            "hmmpress {config[work_dir]}/Pfam-A.hmm",
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
        "logs/annotate_seqs_pfam.log",
    conda:
        "environment.yml",
    shell:  "{params.hmmscan_path} --cpu {threads} -o {params.work_dir}/{params.run_name}_hmmscan.out --tblout {params.work_dir}/{params.run_name}_hmmscan.tblout --domtblout {params.work_dir}/{params.run_name}_hmmscan.domtblout --acc --noali --notextw --cut_ga {input.pfam_database} {input.sequences}"

rule annotate_seqs_go:
    input:
        hhtblout = f"{config['work_dir']}/{config['run_name']}_hmmscan.tblout",
        pfam_gobp = f"{config['work_dir']}/pfam_gobp_most_specific.txt",
    output:
        go_map = f"{config['work_dir']}/{config['run_name']}_GO_map.csv",
    log:
        "logs/annotate_seqs_go.log",
    conda:
        "environment.yml",
    shell:
        "python src/build_go_map.py -o {output.go_map} --hhtblout {input.hhtblout} --pfam_go_files {input.pfam_gobp}"


rule generate_candidates:
    input:
        sequences = f"{config['sequence_path']}",
        go_database = f"{config['work_dir']}/go.obo",
        go_filter = "assets/go_filter.txt",
        go_map = f"{config['work_dir']}/{config['run_name']}_GO_map.csv",
    output:
        kept_proteins = f"{config['work_dir']}/{config['run_name']}_kept_proteins.txt",
        candidates = f"{config['work_dir']}/{config['run_name']}_candidates.tsv",
    params:
        n_pairs = config["dscript"]["n_pairs"],
        seed = config["seed"],
    log:
        "logs/generate_candidates.log",
    conda:
        "environment.yml",
    shell:  "python src/generate_candidates.py --paircount {params.n_pairs} -o {output.candidates} --seq-out {output.kept_proteins} --go_map {input.go_map} --go_database {input.go_database} --go_filter {input.go_filter} --sequences {input.sequences} --seed {params.seed}"

rule predict_network:
    input:
        sequences = f"{config['sequence_path']}",
        candidates = f"{config['work_dir']}/{config['run_name']}_candidates.tsv",
    output:
        network = f"{config['work_dir']}/{config['run_name']}_network.positive.tsv",
    params:
        dscript_path = config["dscript"]["path"],
        dscript_model = f"{config['work_dir']}/{config['dscript']['model']}",
        work_dir = config["work_dir"],
        run_name = config["run_name"],
        device = config["dscript"]["device"],
        t = config["dsd"]["t"],
    resources:
        nvidia_gpu=1
    log:
        "logs/predict_network.log",
    conda:
        "environment.yml",
    shell:  "{params.dscript_path} predict --pairs {input.candidates} --seqs {input.sequences} --model {params.dscript_model} --outfile {params.work_dir}/{params.run_name}_network --device {params.device} --thresh {params.t}"


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
        "logs/compute_distances.log",
    conda:
        "environment.yml",
    shell: "{params.dsd_path} {params.confidence} --converge -t {params.t} --outfile {params.work_dir}/{params.run_name}_distances {input.network}"


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
        "logs/cluster_network.log",
    conda:
        "environment.yml",
    shell:  "python src/cluster_network.py --network_file {input.network} --dsd_file {input.distances} --output {output.clusters} --min_cluster_size {params.min_cluster_size} --cluster_divisor {params.cluster_divisor} --init_k {params.init_k} --sparsity {params.sparsity} --random_seed {params.seed}"


rule reconnect_recipe:
    input:
        clusters = f"{config['work_dir']}/{config['run_name']}_clusters.disconnected.json",
        network = f"{config['work_dir']}/{config['run_name']}_network.positive.tsv",
        dsd = f"{config['work_dir']}/{config['run_name']}_distances.DSD1",
        go_map = f"{config['work_dir']}/{config['run_name']}_GO_map.csv",
        go_db = f"{config['work_dir']}/go.obo"
    output:
        clusters_connected = temp(f"{config['work_dir']}/{config['run_name']}_clusters.recipe.json"),
    params:
        lr = config["recipe"]["lr"],
        cthresh = config["recipe"]["cthresh"],
        max_proteins = config["recipe"]["max_proteins"],
        metric = config["recipe"]["metric"]
    log:
        "logs/reconnect_recipe.log",
    conda:
        "environment.yml",
    shell: "recipe-cluster cook --network-filepath {input.network}  --cluster-filepath {input.clusters} --dsd-file {input.dsd} --protein-go-filepaths {input.go_map} --go-db-file {input.go_db} --lr {params.lr} -cthresh {params.cthresh} --max {params.max_proteins} --metric {params.metric} --cluster_outfile {output.clusters_connected}"
#    shell: "cp {input.clusters} {output.clusters_connected}"

rule add_cluster_functions:
    input:
        clusters = f"{config['work_dir']}/{config['run_name']}_clusters.recipe.json",
        go_map = f"{config['work_dir']}/{config['run_name']}_GO_map.csv",
    output:
        clusters_functional = f"{config['work_dir']}/{config['run_name']}_clusters.functional.json",
    log:
        "logs/add_cluster_functions.log",
    conda:
        "environment.yml",
    shell: "python src/add_cluster_functions.py -o {output.clusters_functional} -cfp {input.clusters} --go_map {input.go_map}"

rule cluster_graph:
    input:
        clusters = f"{config['work_dir']}/{config['run_name']}_clusters.functional.json",
        network = f"{config['work_dir']}/{config['run_name']}_network.positive.tsv",
        go_map = f"{config['work_dir']}/{config['run_name']}_GO_map.csv",
        go_database = f"{config['work_dir']}/go.obo",
    output:
        graph = f"{config['work_dir']}/{config['run_name']}_cluster_graph.json",
    log:
        "logs/cluster_graph.log",
    conda:
        "environment.yml",
    shell:  "python src/build_cluster_graph.py -o {output.graph} -cfp {input.clusters} -nfp {input.network} --go_map {input.go_map} --go_db {input.go_database}"

rule summarize_clusters:
    input:
        clusters = f"{config['work_dir']}/{config['run_name']}_clusters.functional.json",
        go_database = f"{config['work_dir']}/go.obo",
    output:
        human_readable = f"{config['work_dir']}/{config['run_name']}_human_readable.txt",
        readable_json = f"{config['work_dir']}/{config['run_name']}_clusters.json",
    params:
        api_key = f"--api_key {os.environ['OPENAI_API_KEY']}" if config["use_langchain"] else "",
        langchain_model = config["langchain"]["model"],
        do_llm_naming = "--llm_name" if config["use_langchain"] else ""
    log:
        "logs/summarize_clusters.log",
    conda:
        "environment.yml",
    shell:  "python src/summarize_clusters.py {params.do_llm_naming} --model {params.langchain_model} {params.api_key} -o {output.human_readable} --json {output.readable_json} --go_db {input.go_database} -cfp {input.clusters}"

rule vizualize_network:
    input:
        clusters = f"{config['work_dir']}/{config['run_name']}_clusters.functional.json",
        network = f"{config['work_dir']}/{config['run_name']}_network.positive.tsv",
        style = "assets/philharmonic_styles.xml"
    output:
        cytoscape = f"{config['work_dir']}/{config['run_name']}_cytoscape_session.cys",
    params:
        run_name = config["run_name"]
    log:
        "logs/create_cytoscape_session.log",
    conda:
        "environment.yml",
    shell:  "python src/build_cytoscape.py -s {input.style} -o {output.cytoscape} -cfp {input.clusters} -nfp {input.network} --name {params.run_name}"
