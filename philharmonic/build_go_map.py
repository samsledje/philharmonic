from collections import defaultdict
from typing import List

import pandas as pd
import typer
from Bio import SearchIO
from loguru import logger

app = typer.Typer()


@app.command()
def main(
    output: str = typer.Option(..., "-o", "--output", help="Output file"),
    hhtblout: str = typer.Option(..., "--hhtblout", help="hhsearch tblout file"),
    pfam_go_files: List[str] = typer.Option(
        ..., "--pfam_go_files", help="Pfam GO files"
    ),
):
    """Build a GO map from a hhsearch tblout file and a list of Pfam GO files"""
    logger.info("Building GO map...")
    logger.info(f"Output file: {output}")
    logger.info(f"hhsearch tblout file: {hhtblout}")
    logger.info(f"Pfam GO files: {pfam_go_files}")

    hits = defaultdict(list)
    with open(hhtblout) as handle:
        for queryresult in SearchIO.parse(handle, "hmmer3-tab"):
            for hit in queryresult.hits:
                prot_id = hit.query_id
                pfam_id = hit.accession.split(".")[0]
                hits[prot_id].append(pfam_id)

    pfam_go_map = pd.concat(
        [pd.read_csv(f, sep=";", header=0) for f in pfam_go_files],
        ignore_index=True,
    )

    protein_go_map = []
    for prot_id, pfam_hits in hits.items():
        go_terms = pfam_go_map[pfam_go_map["PFAM"].isin(pfam_hits)]["GO"].tolist()
        protein_go_map.append((prot_id, "", ";".join(pfam_hits), ";".join(go_terms)))
    protein_go_dataframe = pd.DataFrame(
        protein_go_map, columns=["prot_id", "manual annot", "pfam_list", "GO_list"]
    )

    protein_go_dataframe.to_csv(output, index=False)


if __name__ == "__main__":
    app()
