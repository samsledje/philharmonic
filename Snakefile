configfile:
    "config.yml"

envvars:
    "OPENAI_API_KEY"

rule PHILHARMONIC:
    input:  
        clusters = expand(
            "{work_dir}/{run_name}_clusters.functional.json",
            work_dir=config["work_dir"], run_name=config["run_name"]
            ),
        cytoscape = expand(
            "{work_dir}/{run_name}_cytoscape_session.cys",
            work_dir=config["work_dir"], run_name=config["run_name"]
            ),
        langchain = expand("{work_dir}/{run_name}_human_readable.txt",
            work_dir=config["work_dir"], run_name=config["run_name"]
            )


rule download_required_files:
    output:
        go_database = expand("{work_dir}/go.obo", work_dir=config["work_dir"]),
        pfam_database_zipped = expand("{work_dir}/Pfam-A.hmm.gz", work_dir=config["work_dir"]),
        pfam_gomf = expand("{work_dir}/pfam_gomf_most_specific.txt", work_dir=config["work_dir"]),
        pfam_gobp = expand("{work_dir}/pfam_gobp_most_specific.txt", work_dir=config["work_dir"]),
        pfam_gocc = expand("{work_dir}/pfam_gocc_most_specific.txt", work_dir=config["work_dir"]),
        dscript_model = expand("{work_dir}/{dscript_model}", work_dir=config["work_dir"], dscript_model=config["dscript"]["model"]),
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
        pfam_database_zipped = expand("{work_dir}/Pfam-A.hmm.gz", work_dir=config["work_dir"]),
    output:
        pfam_database = expand("{work_dir}/Pfam-A.hmm", work_dir=config["work_dir"]),
        pfam_h3m = expand("{work_dir}/Pfam-A.hmm.h3m", work_dir=config["work_dir"]),
        pfam_h3i = expand("{work_dir}/Pfam-A.hmm.h3i", work_dir=config["work_dir"]),
        pfam_h3f = expand("{work_dir}/Pfam-A.hmm.h3f", work_dir=config["work_dir"]),
        pfam_h3p = expand("{work_dir}/Pfam-A.hmm.h3p", work_dir=config["work_dir"]),
    run:
        commands = [
            "gunzip {config[work_dir]}/Pfam-A.hmm.gz",
            "hmmpress {config[work_dir]}/Pfam-A.hmm",
        ]
        for c in commands:
            shell(c)


rule annotate_seqs_pfam:
    input:
        sequences = expand("{sequence_path}", sequence_path=config["sequence_path"]),
        pfam_database = expand("{work_dir}/Pfam-A.hmm", work_dir=config["work_dir"]),
    output:
        pfam_map = expand("{work_dir}/{run_name}_hmmscan.tblout",
            work_dir=config["work_dir"], run_name=config["run_name"]
            )
    # shell: "{config[hmmscan][hmmscan_path]} {output.pfam_map} {input.sequences}"
    shell: "{config[hmmscan][hmmscan_path]} -o {config[work_dir]}/{config[run_name]}_hmmscan.out --tblout {config[work_dir]}/{config[run_name]}_hmmscan.tblout --domtblout {config[work_dir]}/{config[run_name]}_hmmscan.domtblout --acc --noali --notextw --cut_ga {input.pfam_database} {input.sequences}"

rule annotate_seqs_go:
    input:
        pfam_map = expand("{work_dir}/{run_name}_hmmscan.tblout",
            work_dir=config["work_dir"], run_name=config["run_name"]
            ),
        go_database = expand("{work_dir}/go.obo", work_dir=config["work_dir"]),
        pfam_gomf = expand("{work_dir}/pfam_gomf_most_specific.txt", work_dir=config["work_dir"]),
        pfam_gobp = expand("{work_dir}/pfam_gobp_most_specific.txt", work_dir=config["work_dir"]),
        pfam_gocc = expand("{work_dir}/pfam_gocc_most_specific.txt", work_dir=config["work_dir"]),
    output:
        go_map = expand("{work_dir}/{run_name}_GO_map.csv",
            work_dir=config["work_dir"], run_name=config["run_name"]
            )
    shell: "python src/philharmonic_go.py -o {output.go_map} {input.pfam_map} {input.pfam_gomf} {input.pfam_gobp} {input.pfam_gocc}"

rule generate_candidates:
    input:
        sequences = expand("{sequence_path}", sequence_path=config["sequence_path"]),
        protein_shortlist = expand("{protein_shortlist}", protein_shortlist=config["protein_shortlist"]),
        go_shortlist = expand("{go_shortlist}", go_shortlist=config["go_shortlist"]),
        # go_map = expand("{work_dir}/{run_name}_go.tsv",
            # work_dir=config["work_dir"], run_name=config["run_name"]
            # )
    output:
        candidates = expand("{work_dir}/{run_name}_candidates.tsv",
            work_dir=config["work_dir"], run_name=config["run_name"]
            )
    shell: "python src/generate_candidates.py -o {output.candidates} --go_list {input.go_shortlist} --protein_list {input.protein_shortlist} {input.sequences}"

rule predict_network:
    input:
        sequences = expand("{sequence_path}", sequence_path=config["sequence_path"]),
        candidates = expand("{work_dir}/{run_name}_candidates.tsv",
            work_dir=config["work_dir"], run_name=config["run_name"]
            ),
        dscript_model = expand("{work_dir}/{dscript_model}", work_dir=config["work_dir"], dscript_model=config["dscript"]["model"]),
    output:
        network = expand("{work_dir}/{run_name}_network.positive.tsv",
            work_dir=config["work_dir"], run_name=config["run_name"]
            )
    shell: "dscript predict --pairs {input.candidates} --seqs {input.sequences} --model {input.dscript_model} --outfile {config[work_dir]}/{config[run_name]}_network --device {config[device]} --thresh {config[dsd][t]}"


rule compute_distances:
    input: 
        network = expand("{work_dir}/{run_name}_network.positive.tsv",
            work_dir=config["work_dir"], run_name=config["run_name"]
            )
    output:
        distances = expand("{work_dir}/{run_name}_distances.DSD1",
            work_dir=config["work_dir"], run_name=config["run_name"]
            )
    shell: "fastdsd -c --converge -t {config[dsd][t]} --outfile {config[work_dir]}/{config[run_name]}_distances {input.network}"


rule cluster_network:
    input:
        network = expand("{work_dir}/{run_name}_network.positive.tsv",
            work_dir=config["work_dir"], run_name=config["run_name"]
            ),
        distances = expand("{work_dir}/{run_name}_distances.DSD1",
            work_dir=config["work_dir"], run_name=config["run_name"]
            )
    output:
        clusters =expand("{work_dir}/{run_name}_clusters_disconnected.json",
            work_dir=config["work_dir"], run_name=config["run_name"]
            )
    shell: "python src/cluster_network.py --network_file {input.network} --dsd_file {input.distances} --output {output.clusters} --min_cluster_size {config[clustering][min_cluster_size]} --cluster_divisor {config[clustering][cluster_divisor]} --init_k {config[clustering][init_k]}"


rule reconnect_recipe:
    input:
        clusters = expand("{work_dir}/{run_name}_clusters_disconnected.json",
            work_dir=config["work_dir"], run_name=config["run_name"]
            ),
        network = expand("{work_dir}/{run_name}_network.positive.tsv",
            work_dir=config["work_dir"], run_name=config["run_name"]
            )
    output:
        clusters_connected = expand("{work_dir}/{run_name}_clusters.json",
            work_dir=config["work_dir"], run_name=config["run_name"]
            )
    # shell: "recipe --outfile {output.clusters_connected} -cfp {input.clusters}  -nfp {input.network}"
    shell: "cp {input.clusters} {output.clusters_connected}"

rule add_cluster_functions:
    input:
        clusters = expand("{work_dir}/{run_name}_clusters.json",
            work_dir=config["work_dir"], run_name=config["run_name"]
            ),
        go_map = expand("{work_dir}/{run_name}_GO_map.csv",
            work_dir=config["work_dir"], run_name=config["run_name"]
            ),
    output:
        clusters_functional = expand("{work_dir}/{run_name}_clusters.functional.json",
            work_dir=config["work_dir"], run_name=config["run_name"]
            )
    shell: "python src/add_cluster_functions.py -o {output.clusters_functional} -cfp {input.clusters} --go_map {input.go_map}"

rule name_clusters_langchain:
    input:
        clusters = expand("{work_dir}/{run_name}_clusters.functional.json",
            work_dir=config["work_dir"], run_name=config["run_name"]
            ),
        go_database = expand("{work_dir}/go.obo", work_dir=config["work_dir"]),
    output:
        human_readable = expand("{work_dir}/{run_name}_human_readable.txt",
        work_dir=config["work_dir"], run_name=config["run_name"])
    params:
        api_key=os.environ["OPENAI_API_KEY"]
    shell: "python src/name_clusters.py --api_key {params.api_key} -o {output.human_readable} --go_db {input.go_database} -cfp {input.clusters}"

rule create_cytoscape_session:
    input:
        clusters = expand("{work_dir}/{run_name}_clusters.functional.json",
            work_dir=config["work_dir"], run_name=config["run_name"]
            ),
        style = "assets/philharmonic_styles.xml"
    output:
        cytoscape = expand("{work_dir}/{run_name}_cytoscape_session.cys",
        work_dir=config["work_dir"], run_name=config["run_name"])
    shell: "python src/build_cytoscape.py -s {input.style} -o {output.cytoscape} {input.clusters}"

