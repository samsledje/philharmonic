pipeline_figure: config.yml
	touch sample_sequences.fasta;
	snakemake --configfile config.yml --filegraph | dot -Tpng > img/pipeline.png;
	snakemake --configfile config.yml --filegraph | dot -Tsvg > img/pipeline.svg;
	rm sample_sequences.fasta;

format:
	ruff check --select I --fix
	ruff format