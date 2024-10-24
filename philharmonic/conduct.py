import shlex
import subprocess as sp
import sys

import typer
from loguru import logger

app = typer.Typer()
PHILHARMONIC_SNAKEFILE = "Snakefile"


def build_snakemake_command(snakefile: str, config: str, cores: int, **kwargs):
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
    cmd = build_snakemake_command(PHILHARMONIC_SNAKEFILE, config_file, cores, *ctx.args)

    logger.info(f"Running command: {' '.join(cmd)}")
    with sp.Popen(cmd, stdout=sp.PIPE, stderr=sp.STDOUT) as p, open(
        log_file, "ab"
    ) as file:
        assert p.stdout is not None
        for line in p.stdout:  # b'\n'-separated lines
            sys.stdout.buffer.write(line)  # pass bytes as is
            file.write(line)
            file.flush()


if __name__ == "__main__":
    app()
