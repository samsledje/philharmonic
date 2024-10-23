# python src/name_clusters.py --api_key {params.api_key} -o {output.human_readable} --go_db {input.go_database} -cfp {input.clusters}
import json
import os
import shlex
import subprocess as sp

import regex as re
import typer
from loguru import logger
from tqdm import tqdm

from .utils import load_cluster_json, parse_GO_database, print_cluster

app = typer.Typer()

task_instruction = """
You are an expert biologist with a deep understanding of the Gene Ontology. Your job is to give short, intuitive, high-level names to clusters of proteins, given a set of GO terms associated with the proteins and their frequencies.
"""

confidence_score = """
In addition to the name, you should indicate how confident you are that your label is correct, and that it is representative of the function of the cluster. This score should be None, Low, Medium, or High. If you cannot find a connection between the functions in the cluster, give a name of Unknown with a confidence of None.
"""

format_instruction = """
Your response should start with only the name of the cluster on the first line. Then, provide a short paragraph including your explanation for why you gave the cluster that name. Finally, on a new line, provide your confidence score. You should not put blank lines between these sections.
"""

analytical_approach = """
You should try to be as specific as possible in your naming, to avoid overlapping names with other similar clusters. However, the names should still be short and human readable, ideally fewer than 10 words. You should consider the most common GO terms in the cluster, and try to find a common theme or function that ties them together. If you cannot find a common theme, you should give the cluster a name of Unknown.
"""

one_shot_example = """
For example, given the cluster description below:

Cluster of 14 [pdam_00002129-RA,pdam_00001718-RA,...] (hash 1332063120138743063)
Triangles: 27.0
Max Degree: 8
Top Terms:
    GO:0071502 - <cellular response to temperature stimulus> (11)
    GO:0019233 - <sensory perception of pain> (11)
    GO:0042493 - <response to drug> (10)
    GO:0007603 - <phototransduction, visible light> (10)
    GO:0004876 - <complement component C3a receptor activity> (9)

We would name this cluster Temperature, Pain, and Drug Response because there is a high representation for GO terms related to temperature, drug, and pain response.
"""

request = """
Please name the following cluster:
"""

LLM_SYSTEM_TEMPLATE = (
    task_instruction
    + confidence_score
    + format_instruction
    + analytical_approach
    + one_shot_example
    + request
)


def llm_name_cluster(description, model="4o-mini", api_key=None):
    os.environ["OPENAI_API_KEY"] = api_key
    cmd = f"llm --system '{LLM_SYSTEM_TEMPLATE}' -m {model} '{description}' "

    proc = sp.Popen(shlex.split(cmd), stdout=sp.PIPE, stderr=sp.PIPE)
    out, _ = proc.communicate()

    reg = re.compile("(.*$)\s+(.*$)\s+(.*$)", re.MULTILINE)
    regsearch = reg.search(out.decode("utf-8"))
    name = regsearch.group(1)
    explanation = regsearch.group(2)
    confidence = regsearch.group(3)

    logger.info(f"Name: {name}")
    logger.info(f"Explanation: {explanation}")
    logger.info(f"Confidence: {confidence}")

    return name, explanation, confidence


@app.command()
def main(
    output: str = typer.Option(..., "-o", "--output", help="Output file"),
    json_output: str = typer.Option(None, "-j", "--json", help="Output in JSON format"),
    cluster_file_path: str = typer.Option(
        ..., "-cfp", "--cluster-file-path", help="Cluster file"
    ),
    go_db: str = typer.Option(..., "--go_db", help="GO database"),
    llm_name: bool = typer.Option(
        False, help="Use a large language model to name clusters"
    ),
    model: str = typer.Option("Meta-Llama-3-8B-Instruct", help="Language model to use"),
    api_key: str = typer.Option(None, help="OpenAI API key"),
):
    """Summarize clusters"""
    clusters = load_cluster_json(cluster_file_path)
    go_database = parse_GO_database(go_db)

    if llm_name:
        for _, clust in tqdm(clusters.items()):
            if not hasattr(clust, "llm_name"):
                hr = print_cluster(clust, go_database, return_str=True)
                try:
                    name, explanation, confidence = llm_name_cluster(
                        hr, model=model, api_key=api_key
                    )
                    clust["llm_name"] = name
                    clust["llm_explanation"] = explanation
                    clust["llm_confidence"] = confidence
                except Exception as e:
                    logger.error(f"Error: {e}")
                    clust["llm_name"] = "Unknown"
                    clust["llm_explanation"] = "Unknown"
                    clust["llm_confidence"] = "Unknown"

                clust["human_readable"] = print_cluster(
                    clust, go_database, return_str=True
                )
    else:
        for _, clust in tqdm(clusters.items()):
            clust["human_readable"] = print_cluster(clust, go_database, return_str=True)

    if json_output:
        with open(json_output, "w") as f:
            json.dump(clusters, f, indent=4)

    with open(output, "w") as f:
        for _, clust in clusters.items():
            description_string = print_cluster(clust, go_database, return_str=True)
            f.write(str(description_string))
            f.write("\n\n")


if __name__ == "__main__":
    app()
