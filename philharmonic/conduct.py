import shlex
import subprocess as sp
import sys
from importlib.metadata import version as get_version
from pathlib import Path

import typer
from loguru import logger

from .utils import download_file_safe

app = typer.Typer()


def get_snakefile_remote_path(commit: str = "current"):
    """
    Get the path to the snakefile
    """

    if commit == "current":
        installed_version = f"v{get_version('philharmonic')}"
        snakefile_remote_path = f"https://raw.githubusercontent.com/samsledje/philharmonic/refs/tags/{installed_version}/Snakefile"
    else:
        snakefile_remote_path = f"https://raw.githubusercontent.com/samsledje/philharmonic/{commit}/Snakefile"
    return snakefile_remote_path


def download_snakefile(
    download_loc: Path = Path("Snakefile"), commit: str = "current"
) -> Path:
    """
    Download the snakefile from the repo
    """
    snakefile_remote_path = get_snakefile_remote_path(commit)
    local_path = download_loc.resolve()

    if local_path.exists():
        logger.info(f"Snakefile already exists at {local_path}")
        return local_path
    elif download_file_safe(snakefile_remote_path, download_loc):
        return local_path
    else:
        raise FileNotFoundError(
            f"Failed to download snakefile from {snakefile_remote_path}"
        )


def build_snakemake_command(snakefile: Path, config: Path, cores: int, **kwargs):
    """
    Build the snakemake command
    """
    return shlex.split(
        f"snakemake -s {snakefile} --configfile {config} --cores {cores} {' '.join(kwargs)}"
    )


@app.command(
    context_settings={"allow_extra_args": True, "ignore_unknown_options": True}
)
def main(
    ctx: typer.Context,
    config_file: str = typer.Option(..., "-cf", "--config-file", help="Config file"),
    cores: int = typer.Option(8, "-c", "--cores", help="Number of cores to use"),
    log_file: str = typer.Option(
        "philharmonic.log", "-l", "--log-file", help="Log file"
    ),
):
    """
    Main entrypoint into philharmonic
    """
    logger.info(f"Running snakemake with config file: {config_file}")
    logger.info(f"Running snakemake with n_cores: {cores}")

    snakefile = download_snakefile(Path(__file__).parent / "philharmonic_snakefile")
    config_path = Path(config_file).resolve()
    cmd = build_snakemake_command(snakefile, config_path, cores, *ctx.args)

    logger.info(f"Running command: {' '.join(cmd)}")

    with (
        sp.Popen(cmd, stdout=sp.PIPE, stderr=sp.STDOUT) as p,
        open(log_file, "ab") as file,
    ):
        assert p.stdout is not None
        for line in p.stdout:  # b'\n'-separated lines
            sys.stdout.buffer.write(line)  # pass bytes as is
            file.write(line)
            file.flush()


if __name__ == "__main__":
    app()
