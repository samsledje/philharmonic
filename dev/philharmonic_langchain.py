import getpass
import sys
import os
import re
from tqdm import tqdm

with open("OPENAI_API_KEY.txt","r") as f:
    passwd = f.read().strip()
os.environ["OPENAI_API_KEY"] = passwd

from langchain_openai import ChatOpenAI
from langchain_core.messages import HumanMessage, SystemMessage
from langchain_core.output_parsers import StrOutputParser
from langchain_core.prompts import ChatPromptTemplate

model = ChatOpenAI(model="gpt-4")

task_instruction = """
You are an expert biologist with a deep understanding of the Gene Ontology. Your job is to give short, intuitive, high-level names to clusters of proteins, given a set of GO terms associated with the proteins and their frequencies.

"""

confidence_score = """
In addition to the name, you should indicate how confident you are that your label is correct, and that it is representative of the function of the cluster. This score should be between 0 and 1. If you cannot find a connection between the functions in the cluster, give a name of "Unknown" with a confidence of 0.
"""

format_instruction = """
Your response should be a single JSON formatted file for all clusters, with each cluster having the following fields:

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

system_template = task_instruction + confidence_score + format_instruction + analytical_approach + one_shot_example + request

parser = StrOutputParser()
prompt_template = ChatPromptTemplate.from_messages(
    [("system", system_template), ("user", "{text}")]
)
chain = prompt_template | model | parser

cluster_file = sys.argv[1]
with open(cluster_file, "r") as f:
    cluster_description = f.read()

results = []
REGEX = r"Cluster of \d+ \[.*?\] \(hash \d+\)\nTriangles: \d+.\d+\nMax Degree: \d+\nTop Terms:\n(?:\tGO.*?\n)+"
clusters = re.findall(REGEX, cluster_description)

def divide_chunks(l, n):

    # looping till length l
    for i in range(0, len(l), n):
        yield l[i:i + n]

for c in tqdm(divide_chunks(clusters, 50), total=len(clusters) // 50):
    r = chain.invoke({
        "text": c
    })
    results.append(r)

with open(f"{cluster_file}.llm_functions.json","w+") as f:
    f.write("[")
    for r in results:
        f.write(r)
        f.write("\n")
    f.write("]")
