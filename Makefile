install: pyproject.toml poetry.lock
	poetry install

publish: pyproject.toml poetry.lock
	poetry publish --build

pipeline_figure: config.yml
	touch sample_sequences.fasta;
	snakemake --configfile config.yml --filegraph | dot -Tpng > img/pipeline.png;
	snakemake --configfile config.yml --filegraph | dot -Tsvg > img/pipeline.svg;
	rm sample_sequences.fasta;

format:
	poetry run ruff check --fix .
	poetry run ruff format .

test:
	mypy philharmonic --ignore-missing-imports
	pytest --cov=philharmonic --cov-report=term-missing tests
