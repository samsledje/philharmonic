configfile:
    "config.yml"

envvars:
    "OPENAI_API_KEY"

rule PHILHARMONIC:
    input:  
        clusters = f"{config['work_dir']}/{config['run_name']}_clusters.readable.json",
        cytoscape = f"{config['work_dir']}/{config['run_name']}_cytoscape_session.cys",
        human_readable = f"{config['work_dir']}/{config['run_name']}_human_readable.txt",


rule download_required_files:
    output:
        go_database = f"{config['work_dir']}/go.obo",
        pfam_database_zipped = f"{config['work_dir']}/Pfam-A.hmm.gz",
        pfam_gomf = f"{config['work_dir']}/pfam_gomf_most_specific.txt",
        pfam_gobp = f"{config['work_dir']}/pfam_gobp_most_specific.txt",
        pfam_gocc = f"{config['work_dir']}/pfam_gocc_most_specific.txt",
        dscript_model = f"{config['work_dir']}/{config['dscript']['model']}",
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
        pfam_h3m = f"{config['work_dir']}/Pfam-A.hmm.h3m",
        pfam_h3i = f"{config['work_dir']}/Pfam-A.hmm.h3i",
        pfam_h3f = f"{config['work_dir']}/Pfam-A.hmm.h3f",
        pfam_h3p = f"{config['work_dir']}/Pfam-A.hmm.h3p",
    run:
        commands = [
            "gunzip {config[work_dir]}/Pfam-A.hmm.gz",
            "hmmpress {config[work_dir]}/Pfam-A.hmm",
        ]
        for c in commands:
            shell(c)


rule annotate_seqs_pfam:
    input:
        sequences = f"{config['sequence_path']}",
        pfam_database = f"{config['work_dir']}/Pfam-A.hmm",
    output:
        pfam_map = f"{config['work_dir']}/{config['run_name']}_hmmscan.tblout",
    # shell: "{config[hmmscan][hmmscan_path]} {output.pfam_map} {input.sequences}"
    shell: "{config[hmmscan][hmmscan_path]} -o {config[work_dir]}/{config[run_name]}_hmmscan.out --tblout {config[work_dir]}/{config[run_name]}_hmmscan.tblout --domtblout {config[work_dir]}/{config[run_name]}_hmmscan.domtblout --acc --noali --notextw --cut_ga {input.pfam_database} {input.sequences}"

rule annotate_seqs_go:
    input:
        hhtblout = f"{config['work_dir']}/{config['run_name']}_hmmscan.tblout",
        pfam_gomf = f"{config['work_dir']}/pfam_gomf_most_specific.txt",
        pfam_gobp = f"{config['work_dir']}/pfam_gobp_most_specific.txt",
        pfam_gocc = f"{config['work_dir']}/pfam_gocc_most_specific.txt",
    output:
        go_map = f"{config['work_dir']}/{config['run_name']}_GO_map.csv",
    shell: "cp data_old/Pdam_GO_map.csv {output.go_map}"
    # shell: "python src/build_go_map.py -o {output.go_map} --hhtblout {input.pfam_map} --pfam_go_files {input.pfam_gomf} {input.pfam_gobp} {input.pfam_gocc}"

rule generate_candidates:
    input:
        sequences = f"{config['sequence_path']}",
        protein_shortlist = f"{config['protein_shortlist']}",
        go_shortlist = f"{config['go_shortlist']}",
        # go_map = expand("{work_dir}/{run_name}_go.tsv",
            # work_dir=config["work_dir"], run_name=config["run_name"]
            # )
    output:
        candidates = f"{config['work_dir']}/{config['run_name']}_candidates.tsv",
    shell: "python src/generate_candidates.py -o {output.candidates} --go_list {input.go_shortlist} --protein_list {input.protein_shortlist} {input.sequences}"

rule predict_network:
    input:
        sequences = f"{config['sequence_path']}",
        candidates = f"{config['work_dir']}/{config['run_name']}_candidates.tsv",
        dscript_model = f"{config['work_dir']}/{config['dscript']['model']}",
    output:
        network = f"{config['work_dir']}/{config['run_name']}_network.positive.tsv",
    shell: "dscript predict --pairs {input.candidates} --seqs {input.sequences} --model {input.dscript_model} --outfile {config[work_dir]}/{config[run_name]}_network --device {config[device]} --thresh {config[dsd][t]}"


rule compute_distances:
    input: 
        network = f"{config['work_dir']}/{config['run_name']}_network.positive.tsv",
    output:
        distances = f"{config['work_dir']}/{config['run_name']}_distances.DSD1",
    shell: "fastdsd -c --converge -t {config[dsd][t]} --outfile {config[work_dir]}/{config[run_name]}_distances {input.network}"


rule cluster_network:
    input:
        network = f"{config['work_dir']}/{config['run_name']}_network.positive.tsv",
        distances = f"{config['work_dir']}/{config['run_name']}_distances.DSD1",
    output:
        clusters = f"{config['work_dir']}/{config['run_name']}_clusters_disconnected.json",
    shell: "python src/cluster_network.py --network_file {input.network} --dsd_file {input.distances} --output {output.clusters} --min_cluster_size {config[clustering][min_cluster_size]} --cluster_divisor {config[clustering][cluster_divisor]} --init_k {config[clustering][init_k]}"


rule reconnect_recipe:
    input:
        clusters = f"{config['work_dir']}/{config['run_name']}_clusters_disconnected.json",
        network = f"{config['work_dir']}/{config['run_name']}_network.positive.tsv",
    output:
        clusters_connected = f"{config['work_dir']}/{config['run_name']}_clusters.json",
    # shell: "recipe --outfile {output.clusters_connected} -cfp {input.clusters}  -nfp {input.network}"
    shell: "cp {input.clusters} {output.clusters_connected}"

rule add_cluster_functions:
    input:
        clusters = f"{config['work_dir']}/{config['run_name']}_clusters.json",
        go_map = f"{config['work_dir']}/{config['run_name']}_GO_map.csv",
    output:
        clusters_functional = f"{config['work_dir']}/{config['run_name']}_clusters.functional.json",
    shell: "python src/add_cluster_functions.py -o {output.clusters_functional} -cfp {input.clusters} --go_map {input.go_map}"

rule summarize_clusters:
    input:
        clusters = f"{config['work_dir']}/{config['run_name']}_clusters.functional.json",
        go_database = f"{config['work_dir']}/go.obo",
    output:
        human_readable = f"{config['work_dir']}/{config['run_name']}_human_readable.txt",
        readable_json = f"{config['work_dir']}/{config['run_name']}_clusters.readable.json",
    params:
        api_key=os.environ["OPENAI_API_KEY"]
    shell: "python src/summarize_clusters.py --api_key {params.api_key} -o {output.human_readable} --json {output.readable_json} --go_db {input.go_database} -cfp {input.clusters}"

rule create_cytoscape_session:
    input:
        clusters = f"{config['work_dir']}/{config['run_name']}_clusters.functional.json",
        network = f"{config['work_dir']}/{config['run_name']}_network.positive.tsv",
        style = "assets/philharmonic_styles.xml"
    output:
        cytoscape = f"{config['work_dir']}/{config['run_name']}_cytoscape_session.cys",
    shell: "python src/build_cytoscape.py -s {input.style} -o {output.cytoscape} -cfp {input.clusters} -nfp {input.network} --name {config['run_name']}"

