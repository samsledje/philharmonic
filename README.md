<p align="center">
    <img src="img/philharmonic_logo.png" width="400"/>
</p>

# Decoding the Functional Networks of Non-Model Organisms

[![PHILHARMONIC](https://img.shields.io/github/v/release/samsledje/philharmonic?include_prereleases)](https://github.com/samsledje/philharmonic/releases)
[![License](https://img.shields.io/github/license/samsledje/philharmonic)](https://github.com/samsledje/philharmonic/blob/main/LICENSE)
[![Ruff](https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/astral-sh/ruff/main/assets/badge/v2.json)](https://github.com/astral-sh/ruff)
<!-- [![DOI](https://zenodo.org/badge/308463847.svg)](https://zenodo.org/badge/latestdoi/308463847) -->

Protein interaction networks are a fundamental tool for modeling cellular and molecular function, and a large and  sophisticated toolbox has been developed to leverage the structure and topological organization of these networks to predict the functional roles of many under-studied genes, proteins and pathways. However, the overwhelming majority of experimental PPIs from which such networks are constructed come from  humans plus a small number of well-studied model organisms.

We introduce PHILHARMONIC: Protein Human-Transferred Interactome Learns Homology And Recapitulates Model Organism Network Interaction Clusters: a novel computational pipeline for de novo network inference and functional annotation in non-model organisms. PHILHARMONIC uses the D-SCRIPT deep learning method, trained on human PPIs, to learn how to predict PPIs directly from amino acid sequence alone to predict interactions genome-wide, then employs DSD coupled with Spectral Clustering  followed by a new method, Recipe to reconnect clusters. While the predicted PPIs will not individually be completely accurate, the clustering step allows us to aggregate the weaker pairwise signal into confident higher-level organization. We show that these clusters have substantial functional coherence, and we apply our method to predict functionally meaningful modules of proteins in the Coral Holobiont, finding interesting clusters in both the coral animal and symbiont.

## Table of Contents

1. [Installation](#installation)
2. [Usage](#usage)
3. [Interpreting Results](#interpreting-the-results)
4. [Workflow Overview](#workflow-overview)
5. [Detailed Configuration](#detailed-configuration)
6. [Citation](#citation)
7. [FAQ/Known Issues](#issues)
8. [Contributing](#issues)

## Installation

```bash
git clone https://github.com/samsledje/philharmonic.git
cd philharmonic
mamba create -f environment.yml
mamba activate philharmonic
pip install -e .
```

You may also want to install [Cytoscape](https://cytoscape.org/) for visualizing the networks.

## Usage

### Required data

The only data that PHILHARMONIC requires is a set of protein sequences in `.fasta` format. We provide a set of high-level GO terms on which to filter proteins prior to candidate generation and network prediction. You may optionally provide your own set of GO terms, as the `go_filter_path` argument in the configuration file.

### Setting up the config

The `config.yml` file is where you will specify the parameters for the pipeline. We provide a [sample config](config.yml) in this repository
with recommended parameters. You will need to specify the paths to your protein sequences. You can find an explanation for all parameters [below](#detailed-configuration). If you've installed Cytoscape, make sure it is open and running before you start the pipeline. Otherwise, set `build_cytoscape=false` in the configuration. If you use a different configuration file name or location, you can specify it with the `--configfile` flag when running Snakemake.

```yaml
# User Specified
run_name: [identifier for this run]
sequence_path: [path to protein sequences in .fasta format]
work_dir: [path to working directory]
...
```

### Running the pipeline

Once your configuration file is set up, you can invoke the pipeline with

```bash
snakemake -c {number of cores} --configfile {config file}
```

### Pipeline Outputs

We provide a zip of the most relevant output files in `[run].zip`, which contains the following files

```bash
run.zip
|
|-- run_human_readable.txt # Easily readable/scannable list of clusters
|-- run_network.positive.tsv # All edges predicted by D-SCRIPT
|-- run_clusters.json # Main result file, contains all clusters, edges, and functions
|-- run_cluster_graph.tsv # Graph of clusters, where edges are weighted by the number of connections between clusters
|-- run_cluster_graph_functions.tsv # Table of high-level cluster functions from GO Slim
|-- run_GO_map.tsv # Mapping between proteins and GO function labels
```

## Interpreting the Results

### 1. Result Summary

<a target="_blank" href="https://colab.research.google.com/github/samsledje/philharmonic/blob/main/nb/01_result_summary.ipynb">
  <img src="https://colab.research.google.com/assets/colab-badge.svg" alt="Open In Colab"/>
</a>

Using the `clusters.json` file, the `network.positive.tsv` file, the `GO map.tsv` file, and a [GO Slim](https://current.geneontology.org/ontology/subsets/goslim_generic.obo) database, you can view the overall network, a summary of the clustering, and explore individual clusters.

<table border="0" class="dataframe" align="center">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>Network</th>
    </tr>
    <tr>
      <th></th>
      <th></th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>Nodes</th>
      <td>8960</td>
    </tr>
    <tr>
      <th>Edges</th>
      <td>455490</td>
    </tr>
    <tr>
      <th>Degree (Med)</th>
      <td>37.0</td>
    </tr>
    <tr>
      <th>Degree (Avg)</th>
      <td>101.671875</td>
    </tr>
    <tr>
      <th>Sparsity</th>
      <td>0.005674</td>
    </tr>
  </tbody>
</table>
 
<p align="center">
 <img src="img/readme_sample_cluster.jpg" width="400"/>
</p>

```bash
Cluster of 21 proteins [pdam_00022258-RA, pdam_00005419-RA, pdam_00017455-RA, ...] (hash 549662403768153899)
7 proteins re-added by ReCIPE (degree, 0.75)
Edges: 0
Triangles: 0
Max Degree: 0
Top Terms:
		GO:0071380 - <cellular response to prostaglandin E stimulus> (14)
		GO:0022900 - <electron transport chain> (14)
		GO:0019233 - <sensory perception of pain> (14)
		GO:0008542 - <visual learning> (13)
		GO:0010759 - <positive regulation of macrophage chemotaxis> (13)
```

### 2. Functional Permutation Testing

<a target="_blank" href="https://colab.research.google.com/github/samsledje/philharmonic/blob/main/nb/02_functional_permutation_test.ipynb">
  <img src="https://colab.research.google.com/assets/colab-badge.svg" alt="Open In Colab"/>
</a>

Using the same files, you can run a statistical test of cluster function by permuting cluster labels, and computing the [Jaccard similarity](https://en.wikipedia.org/wiki/Jaccard_index) between terms in the same cluster.

![function enrichment](img/readme_function_enrichment.png)

### 3. Full network in Cytoscape

1. Load `network.positive.tsv` using `File -> Import -> Network from File`

### 4. Cluster graph in Cytoscape

1. Load `cluster_graph.tsv` using `File -> Import -> Network from File`
2. Load `cluster_graph_functions.tsv` using `File -> Import -> Table from File`
3. Add a `Column filter` on the `Edge: weight` attribute, selecting edges greater than ~50-100 weight
4. `Select -> Nodes -> Nodes Connected by Selected Edges` to subset the nodes
5. Create the subgraph with `File -> New Network -> From Selected Nodes, Selected Edges`
6. Layout the network with your layout of choice, we recommend `Layout -> Prefuse Force Directed Layout -> weight`
7. Add node colors using the [PHILHARMONIC style](assets/philharmonic_styles.xml), imported with `File -> Import -> Styles from File`

![cluster graph](img/readme_cluster_graph.svg)

## Workflow Overview

A detailed overview of PHILHARMNONIC can be found in the [manuscript](#citation). We briefly outline the pipeline below.

Each of these steps can be invoked independently by running `snakemake -c {number of cores} {target}`. The `{target}` is shown in parentheses following each step below.

![snakemake pipeline](img/pipeline.png)

1. Download necessary files (`download_required_files`)
2. Run [hmmscan](http://hmmer.org/) on protein sequences to annotate pfam domains (`annotate_seqs_pfam`)
3. Use pfam-go associations to add [GO terms](https://geneontology.org/) to sequences (`annotate_seqs_go`)
4. Generate candidate pairs (`generate_candidates`)
5. Use [D-SCRIPT](https://dscript.csail.mit.edu/) to predict network (`predict_network`)
6. Compute node distances with [FastDSD](https://github.com/samsledje/fastDSD) (`compute_distances`)
7. Cluster the network with [spectral clustering](https://scikit-learn.org/stable/modules/generated/sklearn.cluster.SpectralClustering.html) (`cluster_network`)
8. Use [ReCIPE](https://pypi.org/project/recipe-cluster/) to reconnect clusters (`reconnect_recipe`)
9. Annotate clusters with functions (`add_cluster_functions`)
10. Compute cluster graph (`cluster_graph`)
11. Name and describe clusters with [Langchain](https://www.langchain.com/) (`summarize_clusters`)
<!-- 12. Visualize in [Cytoscape](https://cytoscape.org/) (`vizualize_network`) -->

## Detailed Configuration

The `config.yml` file contains various parameters that control the behavior of the PHILHARMONIC pipeline. Below is a detailed explanation of each parameter, including default values:

### User Specified

- `run_name`: Identifier for this run [required]
- `sequence_path`: Path to protein sequences in .fasta format [required]
- `work_dir`: Path to the working directory where results will be stored (default: "results")
- `use_cytoscape`: Boolean flag to enable/disable Cytoscape visualization (default: false)
- `use_langchain`: Boolean flag to enable/disable Langchain for cluster summarization (default: false)

### General Parameters

- `seed`: Random seed for reproducibility (default: 6191998)

### hmmscan Parameters

- `hmmscan.path`: Path to the hmmscan executable (default: "hmmscan")
- `hmmscan.threads`: Number of threads to use for hmmscan (default: 16)

### D-SCRIPT Parameters

- `dscript.path`: Path to the D-SCRIPT executable (default: "dscript")
- `dscript.n_pairs`: Number of protein pairs to predict (-1 for all pairs) (default: -1)
- `dscript.model`: Pre-trained D-SCRIPT model to use. (default: "samsl/human_v1")
- `dscript.device`: GPU device to use (-1 for CPU) (default: 0)

### DSD Parameters

- `dsd.path`: Path to the FastDSD executable (default: "fastdsd")
- `dsd.t`: Edge existence threshold for DSD algorithm (default: 0.5)
- `dsd.confidence`: Boolean flag to use confidence scores (default: true)

### Clustering Parameters

- `clustering.init_k`: Initial number of clusters for spectral clustering (default: 500)
- `clustering.min_cluster_size`: Minimum size of a cluster (default: 3)
- `clustering.cluster_divisor`: Divisor used to determine the final number of clusters (default: 20)
- `clustering.sparsity_thresh`: Sparsity threshold for filtering edges (default: 1e-5)

### ReCIPE Parameters

- `recipe.lr`: Linear ratio for ReCIPE algorithm (default: 0.1)
- `recipe.cthresh`: Connectivity threshold to add proteins until for ReCIPE (default: 0.75)
- `recipe.max_proteins`: Maximum number of proteins to add to a cluster in ReCIPE (default: 20)
- `recipe.metric`: Metric to use for ReCIPE (default: "degree")

### Langchain Parameters

- `langchain.model`: Language model to use for cluster summarization (default: "gpt-4o")

## Citation

```bibtex
TBD
```

## Issues

- On Linux, the package `plac` may not install properly with the included `environment.yml`. If you are seeing the error `No module names 'asyncore'`, try running `mamba update plac`

## Contributing

```bash
git clone https://github.com/samsledje/philharmonic.git
cd philharmonic
mamba create -f environment.yml
mamba activate philharmonic
pip install -e .
mamba install -c conda-forge pre-commit
pre-commit install
git checkout -b [feature branch]
```
