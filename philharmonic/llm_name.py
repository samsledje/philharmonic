import json
import typer
import shlex
import subprocess as sp
from tqdm import tqdm
from loguru import logger

system = "You are an expert biologist. Give a short, descriptive, human-readable name to the following cluster, based on your expertise in gene function and the gene ontology. Respond with only the name of the cluster in the first line, and nothing else"

app = typer.Typer()

@app.command()
def main(
        file_path: str = typer.Option(..., "-f", help="Path to JSON"),
        model_name: str = typer.Option("Meta-Llama-3-8B-Instruct","-m","--model",help="Model to use"),
        first_only: bool = typer.Option(..., help="Only run on the first cluster")
        ):

        with open(file_path,"r") as f:
            cdict = json.load(f)

        llm_names = {}
        for k, clust in tqdm(cdict.items(), total=len(cdict)):
            hr = clust["human_readable"]
            cmd = f"llm --system '{system}' -m {model_name} "
            proc = sp.Popen(shlex.split(cmd) + [ hr ], stdout=sp.PIPE, stderr=sp.PIPE)
            out,err = proc.communicate()
            name = out.decode('utf-8').split('<')[0]
            llm_names[k] = name
            logger.info(f"{k} - {name}")
            clust["llm_name"] = name

        with open(f"{file_path}.llm","w+") as f:
            json.dump(cdict, f)

        with open(f"{file_path}.llm_alone","w+") as f:
            json.dump(llm_names, f)


if __name__ == "__main__":
    app()
