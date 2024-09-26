import typer
from importlib.metadata import version as get_version
from . import (
    build_go_map,
    build_cytoscape,
    summarize_clusters,
    add_cluster_functions,
    cluster_network,
    build_cluster_graph,
    generate_candidates,
)

app = typer.Typer()

# ... existing imports ...


def version_callback(value: bool):
    if value:
        typer.echo(f"Philharmonic version: {get_version('philharmonic')}")
        raise typer.Exit()


@app.callback()
def main(
    version: bool = typer.Option(
        None,
        "--version",
        "-v",
        callback=version_callback,
        is_eager=True,
        help="Show the version and exit.",
    ),
):
    """
    Philharmonic: Decoding functional organization in non-model organisms
    """


# Register subcommands
app.command("build-go-map")(build_go_map.main)
app.command("build-cytoscape")(build_cytoscape.main)
app.command("summarize-clusters")(summarize_clusters.main)
app.command("add-cluster-functions")(add_cluster_functions.main)
app.command("cluster-network")(cluster_network.main)
app.command("build-cluster-graph")(build_cluster_graph.main)
app.command("generate-candidates")(generate_candidates.main)

if __name__ == "__main__":
    app()
