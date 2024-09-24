pipeline_figure: config.yml
	touch sample_sequences.fasta;
	snakemake --configfile config.yml --filegraph | dot -Tpng > pipeline.png;
	snakemake --configfile config.yml --filegraph | dot -Tsvg > pipeline.svg;
	rm sample_sequences.fasta;