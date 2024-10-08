# python src/name_clusters.py --api_key {params.api_key} -o {output.human_readable} --go_db {input.go_database} -cfp {input.clusters}
import json
import os

import typer
from langchain.prompts import ChatPromptTemplate
from langchain.schema import StrOutputParser
from langchain_openai import ChatOpenAI
from loguru import logger
from tqdm import tqdm

from .utils import load_cluster_json, parse_GO_database, print_cluster

app = typer.Typer()

task_instruction = """
You are an expert biologist with a deep understanding of the Gene Ontology. Your job is to give short, intuitive, high-level names to clusters of proteins, given a set of GO terms associated with the proteins and their frequencies.

"""

confidence_score = """
In addition to the name, you should indicate how confident you are that your label is correct, and that it is representative of the function of the cluster. This score should be between 0 and 1. If you cannot find a connection between the functions in the cluster, give a name of "Unknown" with a confidence of 0.
"""

format_instruction = """
Your response should be a single JSON formatted response, with no other text other than the JSON string. The response should have the following fields:

"hash": [cluster hash],
"short_name": [your intuitive name],
"explanation": [a one paragraph justification for the intuitive name]
"confidence_score": [a score between 0.0 and 1.0 indicating how confident you are that this name accurately reflects the function of the cluster]
"""

analytical_approach = """
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

We would name this cluster "Sensory Response" because there is a high representation for GO terms related to temperature, drug, and pain response.
"""

request = """
Please name the following cluster:
"""

llm_system_template = (
    task_instruction
    + confidence_score
    + format_instruction
    + analytical_approach
    + one_shot_example
    + request
)


def llm_summarize(cluster, go_database, model, sys_template, api_key):
    parser = StrOutputParser()
    prompt_template = ChatPromptTemplate.from_messages(
        [("system", sys_template), ("user", "{text}")]
    )
    chain = prompt_template | model | parser
    r = chain.invoke({"text": print_cluster(cluster, go_database, return_str=True)})

    # return r

    return json.loads(r)


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
    model: str = typer.Option(None, help="Language model to use"),
    api_key: str = typer.Option(None, help="OpenAI API key"),
):
    """Summarize clusters"""
    clusters = load_cluster_json(cluster_file_path)
    go_database = parse_GO_database(go_db)

    if llm_name:
        # with open("OPENAI_API_KEY.txt","r") as f:
        # passwd = f.read().strip()
        os.environ["OPENAI_API_KEY"] = api_key

        model = ChatOpenAI(
            model="gpt-4o",
            temperature=0,
            # max_tokens=200
        )

        for _, clust in tqdm(clusters.items()):
            if not hasattr(clust, "llm_name"):
                try:
                    llm_summary = llm_summarize(
                        clust, go_database, model, llm_system_template, api_key
                    )
                    clust["llm_name"] = llm_summary["short_name"]
                    clust["llm_explanation"] = llm_summary["explanation"]
                    clust["llm_confidence"] = llm_summary["confidence_score"]
                except Exception as e:
                    logger.error(f"Error: {e}")
                    clust["llm_name"] = "Unknown"
                    clust["llm_explanation"] = "Unknown"
                    clust["llm_confidence"] = "Unknown"

    if json_output:
        for _, clust in clusters.items():
            clust["human_readable"] = print_cluster(clust, go_database, return_str=True)
        with open(json_output, "w") as f:
            json.dump(clusters, f, indent=4)

    with open(output, "w") as f:
        for _, clust in clusters.items():
            description_string = print_cluster(clust, go_database, return_str=True)
            f.write(description_string)
            f.write("\n\n")


if __name__ == "__main__":
    app()
