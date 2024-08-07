# python src/name_clusters.py --api_key {params.api_key} -o {output.human_readable} --go_db {input.go_database} -cfp {input.clusters}
import os
import re
import json
import argparse
from tqdm import tqdm
from pathlib import Path
from langchain_openai import ChatOpenAI
from langchain_core.output_parsers import StrOutputParser
from langchain_core.prompts import ChatPromptTemplate

from utils import load_cluster_json, parse_GO_database, Cluster, log

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

llm_system_template = task_instruction + confidence_score + format_instruction + analytical_approach + one_shot_example + request

def llm_summarize(cluster, go_database, model, sys_template, api_key):
    
    parser = StrOutputParser()
    prompt_template = ChatPromptTemplate.from_messages(
        [("system", sys_template), ("user", "{text}")]
    )
    chain = prompt_template | model | parser
    r = chain.invoke({
            "text": cluster_format(cluster, go_database)
        })
    
    # return r
    
    return json.loads(r)

def cluster_format(clust, go_database, n_terms = 5):

    clust = Cluster(clust)
    description_string = ""

    if hasattr(clust, 'llm_name'):
        description_string += f"Cluster Name: {clust.llm_name}\n"

    members = clust._member_list()
    description_string += f"Cluster of {len(clust.members)} proteins [{members}] (hash {hash(clust)})\n"

    if hasattr(clust, 'G'):
        description_string += f"Edges: {len(clust.G.edges())}\n"
        description_string += f"Triangles: {clust.triangles()}\n"
        description_string += f"Max Degree: {0 if not len(clust.G.edges()) else max(clust.G.degree(), key=lambda x: x[1])[1]}\n"

    if hasattr(clust, 'GO_terms'):
        top_terms = clust.get_top_terms(n_terms)
        description_string += "Top Terms:\n"
        for gid, freq in top_terms:
            try:
                go_name = go_database[gid]
            except KeyError:
                go_name = "Unknown"
            description_string += f"\t\t{gid} - <{go_name}> ({freq})\n"

    if hasattr(clust, 'llm_explanation'):
        description_string += f"LLM Explanation: {clust.llm_explanation}\n"
        description_string += f"LLM Confidence: {clust.llm_confidence}\n"

    return description_string


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Create a human readable summary for results of a cluster analysis')
    parser.add_argument('-o', '--output', required=True, type=str, help='Output file')
    parser.add_argument('--json', type=str, default=None, help='Output in JSON format')
    parser.add_argument('-cfp', '--cluster_file_path', required=True, type=str, help='Cluster file')
    parser.add_argument('--go_db', required=True, type=str, help='GO database')
    parser.add_argument('--llm_name', action='store_true', help='Use a large language model to name clusters')
    parser.add_argument('--model', type=str, help='Language model to use')
    parser.add_argument('--api_key', type=str, help='OpenAI API key')

    args = parser.parse_args()

    clusters = load_cluster_json(args.cluster_file_path)
    go_database = parse_GO_database(args.go_db)

    if args.llm_name:
        # with open("OPENAI_API_KEY.txt","r") as f:
            # passwd = f.read().strip()
        os.environ["OPENAI_API_KEY"] = args.api_key

        model = ChatOpenAI(model="gpt-4o")

        for k, clust in tqdm(clusters.items()):
            if not hasattr(clust, 'llm_name'):
                try:
                    llm_summary = llm_summarize(clust, go_database, model, llm_system_template, args.api_key)
                    clust["llm_name"] = llm_summary["short_name"]
                    clust["llm_explanation"] = llm_summary["explanation"]
                    clust["llm_confidence"] = llm_summary["confidence_score"]
                except Exception as e:
                    log(f"Error: {e}")
                    clust["llm_name"] = "Unknown"
                    clust["llm_explanation"] = "Unknown"
                    clust["llm_confidence"] = "Unknown"
    
    if args.json:
        for k, clust in clusters.items():
            clust["human_readable"] = cluster_format(clust, go_database)
        with open(args.json, 'w') as f:
            json.dump(clusters, f, indent=4)
    
    with open(args.output, 'w') as f:
        for k, clust in clusters.items():
            description_string = cluster_format(clust, go_database)
            f.write(description_string)
            f.write("\n\n")