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

from utils import load_cluster_json, parse_GO_database, Cluster

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

def llm_chain_invoker(model, cluster):
    

    parser = StrOutputParser()
    prompt_template = ChatPromptTemplate.from_messages(
        [("system", system_template), ("user", "{text}")]
    )
    chain = prompt_template | model | parser
    r = chain.invoke({
            "text": str(c)
        })
    
    return json.loads(r)

def llm_format(clust, go_database, n_terms = 5):

    clust = Cluster(clust)

    members = clust._member_list()
    description_string = f"Cluster of {len(members)} proteins\n"

    if hasattr(clust, 'G'):
        description_string += f"Triangles: {clust.triangles()}\n"
        description_string += f"Max Degree: {max(clust.G.degree(), key=lambda x: x[1])[1]}\n"

    if hasattr(clust, 'GO_terms'):
        top_terms = clust.get_top_terms(n_terms)
        description_string += "Top Terms:\n"
        for gid, freq in top_terms:
            go_name = go_database[gid]["name"]
            description_string += f"\t{gid} - <{go_name}> ({freq})\n"

    return description_string


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Create a human readable summary for results of a cluster analysis')
    parser.add_argument('-o', '--output', required=True, type=str, help='Output file')
    parser.add_argument('-cfp', '--cluster_file_path', required=True, type=str, help='Cluster file')
    parser.add_argument('--go_db', required=True, type=str, help='GO database')
    parser.add_argument('--api_key', required=True, type=str, help='OpenAI API key')

    args = parser.parse_args()

    clusters = load_cluster_json(args.cluster_file_path)
    go_database = parse_GO_database(args.go_db)

    # with open("OPENAI_API_KEY.txt","r") as f:
    #     passwd = f.read().strip()
    # os.environ["OPENAI_API_KEY"] = args.api_key

    # model = ChatOpenAI(model="gpt-4")

    # results = []
    # REGEX = r"Cluster of \d+ \[.*?\] \(hash \d+\)\nTriangles: \d+.\d+\nMax Degree: \d+\nTop Terms:\n(?:\tGO.*?\n)+"
    # clusters = re.findall(REGEX, cluster_description)

    # def divide_chunks(l, n):

    #     # looping till length l
    #     for i in range(0, len(l), n):
    #         yield l[i:i + n]

    # for c in tqdm(divide_chunks(clusters, 50), total=len(clusters) // 50):
        
    #     results.append(r)

    # with open(f"{cluster_file}.llm_functions.json","w+") as f:
    #     f.write("[")
    #     for r in results:
    #         f.write(r)
    #         f.write("\n")
    #     f.write("]")

    # clusters_human_readable = {}
    # for k, clust in tqdm(clusters.items()):
    #     c = Cluster(clust)
    #     llm_response = llm_chain_invoker(str(c), model)
    #     clust.update(llm_response)
    #     clusters_human_readable[k] = clust
        # f.write(json.dumps({"hash": k, "short_name": short_name, "explanation": explanation, "confidence_score": confidence_score}))
        # f.write("\n")

    # with open(args.output, 'w') as f:
    #     json.dump(clusters, f)