configfile:
    "config.yml"

envvars:
    "OPENAI_API_KEY"

rule PHILHARMONIC:
    input:  
        clusters = f"{config['work_dir']}/{config['run_name']}_clusters.json",
        cytoscape = f"{config['work_dir']}/{config['run_name']}_cytoscape_session.cys" if config["build_cytoscape"] else [],
        human_readable = f"{config['work_dir']}/{config['run_name']}_human_readable.txt",


rule download_required_files:
    output:
        go_database = f"{config['work_dir']}/go.obo",
        pfam_database_zipped = f"{config['work_dir']}/Pfam-A.hmm.gz",
        pfam_gomf = f"{config['work_dir']}/pfam_gomf_most_specific.txt",
        pfam_gobp = f"{config['work_dir']}/pfam_gobp_most_specific.txt",
        pfam_gocc = f"{config['work_dir']}/pfam_gocc_most_specific.txt",
        dscript_model = f"{config['work_dir']}/{config['dscript']['model']}",
    log:
        "logs/download_required_files.log",
    run:
        commands = [
            "mkdir -p {config[work_dir]}",
            "curl https://current.geneontology.org/ontology/go.obo -o {config[work_dir]}/go.obo",
            "curl ftp://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz -o {config[work_dir]}/Pfam-A.hmm.gz",
            "curl https://godm.loria.fr/data/pfam_gomf_most_specific.txt -o {config[work_dir]}/pfam_gomf_most_specific.txt",
            "curl https://godm.loria.fr/data/pfam_gobp_most_specific.txt -o {config[work_dir]}/pfam_gobp_most_specific.txt",
            "curl https://godm.loria.fr/data/pfam_gocc_most_specific.txt -o {config[work_dir]}/pfam_gocc_most_specific.txt",
            "curl http://cb.csail.mit.edu/cb/dscript/data/models/{config[dscript][model]} -o {config[work_dir]}/{config[dscript][model]}",
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
        pfam_map = temp(f"{config['work_dir']}/{config['run_name']}_hmmscan.tblout"),
    threads: config["hmmscan"]["threads"]
    params:
        work_dir = config["work_dir"],
        run_name = config["run_name"],
        hmmscan_path = config["hmmscan"]["path"]
    log:
        "logs/annotate_seqs_pfam.log",
    conda:
        "env.yml",
    shell:  "{params.hmmscan_path} --cpu {threads} -o {params.work_dir}/{params.run_name}_hmmscan.out --tblout {params.work_dir}/{params.run_name}_hmmscan.tblout --domtblout {params.work_dir}/{params.run_name}_hmmscan.domtblout --acc --noali --notextw --cut_ga {input.pfam_database} {input.sequences}"

rule annotate_seqs_go:
    input:
        hhtblout = f"{config['work_dir']}/{config['run_name']}_hmmscan.tblout",
        pfam_gomf = f"{config['work_dir']}/pfam_gomf_most_specific.txt",
        pfam_gobp = f"{config['work_dir']}/pfam_gobp_most_specific.txt",
        pfam_gocc = f"{config['work_dir']}/pfam_gocc_most_specific.txt",
    output:
        go_map = f"{config['work_dir']}/{config['run_name']}_GO_map.csv",
    log:
        "logs/annotate_seqs_go.log",
    shell:
        "python src/build_go_map.py -o {output.go_map} --hhtblout {input.hhtblout} --pfam_go_files {input.pfam_gomf} {input.pfam_gobp} {input.pfam_gocc}"


rule generate_candidates:
    input:
        sequences = f"{config['sequence_path']}",
        protein_shortlist = f"{config['protein_shortlist']}",
        go_database = f"{config['work_dir']}/go.obo",
        go_shortlist = f"{config['go_shortlist']}",
        go_filter = "assets/go_level2_marked-up.csv",
        go_map = f"{config['work_dir']}/{config['run_name']}_GO_map.csv",
    output:
        candidates = f"{config['work_dir']}/{config['run_name']}_candidates.tsv",
    params:
        n_pairs = config["n_pairs"],
        manual_annot_wt = config["manual_annot_wt"],
    log:
        "logs/generate_candidates.log",
    conda:
        "env.yml",
    shell:  "python src/generate_candidates.py --manual_annot_wt {params.manual_annot_wt} --paircount {params.n_pairs} -o {output.candidates} --go_map {input.go_map} --go_database {input.go_database} --go_filter {input.go_filter} --go_list {input.go_shortlist} --protein_list {input.protein_shortlist} --sequences {input.sequences}"

rule predict_network:
    input:
        sequences = f"{config['sequence_path']}",
        candidates = f"{config['work_dir']}/{config['run_name']}_candidates.tsv",
        dscript_model = f"{config['work_dir']}/{config['dscript']['model']}",
    output:
        network = f"{config['work_dir']}/{config['run_name']}_network.positive.tsv",
    params:
        work_dir = config["work_dir"],
        run_name = config["run_name"],
        device = config["dscript"]["device"],
        t = config["dsd"]["t"]
    resources:
        nvidia_gpu=1
    log:
        "logs/predict_network.log",
    conda:
        "env.yml",
    shell:  "dscript predict --pairs {input.candidates} --seqs {input.sequences} --model {input.dscript_model} --outfile {params.work_dir}/{params.run_name}_network --device {params.device} --thresh {params.t}"


rule compute_distances:
    input: 
        network = f"{config['work_dir']}/{config['run_name']}_network.positive.tsv",
    output:
        distances = temp(f"{config['work_dir']}/{config['run_name']}_distances.DSD1"),
    params:
        work_dir = config["work_dir"],
        run_name = config["run_name"],
        t = config["dsd"]["t"]
    log:
        "logs/compute_distances.log",
    conda:
        "env.yml",
    shell: "fastdsd -c --converge -t {params.t} --outfile {params.work_dir}/{params.run_name}_distances {input.network}"


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
        init_k = config["clustering"]["init_k"]
    log:
        "logs/cluster_network.log",
    conda:
        "env.yml",
    shell:  "python src/cluster_network.py --network_file {input.network} --dsd_file {input.distances} --output {output.clusters} --min_cluster_size {params.min_cluster_size} --cluster_divisor {params.cluster_divisor} --init_k {params.init_k}"


rule reconnect_recipe:
    input:
        clusters = f"{config['work_dir']}/{config['run_name']}_clusters.disconnected.json",
        network = f"{config['work_dir']}/{config['run_name']}_network.positive.tsv",
    output:
        clusters_connected = temp(f"{config['work_dir']}/{config['run_name']}_clusters.recipe.json"),
    # shell: "recipe --outfile {output.clusters_connected} -cfp {input.clusters}  -nfp {input.network}"
    log:
        "logs/reconnect_recipe.log",
    conda:
        "env.yml",
    shell: "cp {input.clusters} {output.clusters_connected}"

rule add_cluster_functions:
    input:
        clusters = f"{config['work_dir']}/{config['run_name']}_clusters.recipe.json",
        go_map = f"{config['work_dir']}/{config['run_name']}_GO_map.csv",
    output:
        clusters_functional = f"{config['work_dir']}/{config['run_name']}_clusters.functional.json",
    log:
        "logs/add_cluster_functions.log",
    conda:
        "env.yml",
    shell: "python src/add_cluster_functions.py -o {output.clusters_functional} -cfp {input.clusters} --go_map {input.go_map}"

rule summarize_clusters:
    input:
        clusters = f"{config['work_dir']}/{config['run_name']}_clusters.functional.json",
        go_database = f"{config['work_dir']}/go.obo",
    output:
        human_readable = f"{config['work_dir']}/{config['run_name']}_human_readable.txt",
        readable_json = f"{config['work_dir']}/{config['run_name']}_clusters.json",
    params:
        api_key = os.environ["OPENAI_API_KEY"],
        langchain_model = config["langchain"]["model"],
        do_llm_naming = "--llm_name" if config["langchain"]["do_llm"] else ""
    log:
        "logs/summarize_clusters.log",
    conda:
        "env.yml",
    shell:  "python src/summarize_clusters.py {params.do_llm_naming} --model {params.langchain_model} --api_key {params.api_key} -o {output.human_readable} --json {output.readable_json} --go_db {input.go_database} -cfp {input.clusters}"

rule create_cytoscape_session:
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
        "env.yml",
    shell:  "python src/build_cytoscape.py -s {input.style} -o {output.cytoscape} -cfp {input.clusters} -nfp {input.network} --name {params.run_name}"

