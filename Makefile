CURRENT_TAG := $(shell git describe --tags --abbrev=0)
NEXT_MAJOR := $(shell poetry version major --dry-run -s)
NEXT_MINOR := $(shell poetry version minor --dry-run -s)
NEXT_PATCH := $(shell poetry version patch --dry-run -s)
NEXT_DYNAMIC := $(shell dunamai from git --style pep440 --bump)

.PHONY : help major_version minor_version patch_version dynamic_version install pipeline_figure format test

version:
	@echo "Current version: $(CURRENT_TAG)"
	@echo "Dynamic version: $(NEXT_DYNAMIC) (make install)"
	@echo "Next major: $(NEXT_MAJOR) (make bump-major)"
	@echo "Next minor: $(NEXT_MINOR) (make bump-minor)"
	@echo "Next patch: $(NEXT_PATCH) (make bump-patch)"

bump-major:
	@echo "Current version will be bumped to: $(NEXT_MAJOR)"
	poetry version $(NEXT_MAJOR)
	git add pyproject.toml
	git commit -m "build: bump version to $(NEXT_MAJOR)"
	git tag -a "v$(NEXT_MAJOR)" -m "Release version $(NEXT_MAJOR)"

bump-minor:
	@echo "Current version will be bumped to: $(NEXT_MINOR)"
	poetry version $(NEXT_MINOR)
	git add pyproject.toml
	git commit -m "build: bump version to $(NEXT_MINOR)"
	git tag -a "v$(NEXT_MINOR)" -m "Release version $(NEXT_MINOR)"

bump-patch:
	@echo "Current version will be bumped to: $(NEXT_PATCH)"
	poetry version $(NEXT_PATCH)
	git add pyproject.toml
	git commit -m "build: bump version to $(NEXT_PATCH)"
	git tag -a "v$(NEXT_PATCH)" -m "Release version $(NEXT_PATCH)"

bump-dynamic:
	@echo "Current dynamic version is: $(NEXT_DYNAMIC)"
	poetry version $(NEXT_DYNAMIC)

install: bump-dynamic
	poetry install

pipeline_figure: config.yml
	touch sample_sequences.fasta;
	snakemake --configfile config.yml --filegraph | dot -Tpng > img/pipeline.png;
	snakemake --configfile config.yml --filegraph | dot -Tsvg > img/pipeline.svg;
	rm sample_sequences.fasta;

format:
	@echo "Formatting project with ruff"
	poetry run ruff check --fix .
	poetry run ruff format .

test:
	@echo "Running tests and type checks"
	mypy philharmonic --ignore-missing-imports
	pytest --cov=philharmonic --cov-report=term-missing tests
