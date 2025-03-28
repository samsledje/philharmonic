#################################################################################
# GLOBALS                                                                       #
#################################################################################

PROJECT_NAME = philharmonic
PYTHON_VERSION = 3.11
PYTHON_INTERPRETER = python

CURRENT_TAG := $(shell git describe --tags --abbrev=0)
NEXT_MAJOR := $(shell poetry version major --dry-run -s)
NEXT_MINOR := $(shell poetry version minor --dry-run -s)
NEXT_PATCH := $(shell poetry version patch --dry-run -s)
NEXT_DYNAMIC := $(shell dunamai from git --style pep440 --bump)

## Show the current version and the next versions
.PHONY: version
version:
	@echo "Current version: $(CURRENT_TAG)"
	@echo "Dynamic version: $(NEXT_DYNAMIC) (make install)"
	@echo "Next major: $(NEXT_MAJOR) (make bump-major)"
	@echo "Next minor: $(NEXT_MINOR) (make bump-minor)"
	@echo "Next patch: $(NEXT_PATCH) (make bump-patch)"

## Bump major version
.PHONY: bump-major
bump-major:
	@echo "Current version will be bumped to: $(NEXT_MAJOR)"
	poetry version $(NEXT_MAJOR)
	git add pyproject.toml
	git commit -m "build: bump version to $(NEXT_MAJOR)"
	git tag -a "v$(NEXT_MAJOR)" -m "Release version $(NEXT_MAJOR)"

## Bump minor version
.PHONY: bump-minor
bump-minor:
	@echo "Current version will be bumped to: $(NEXT_MINOR)"
	poetry version $(NEXT_MINOR)
	git add pyproject.toml
	git commit -m "build: bump version to $(NEXT_MINOR)"
	git tag -a "v$(NEXT_MINOR)" -m "Release version $(NEXT_MINOR)"

## Bump patch version
.PHONY: bump-patch
bump-patch:
	@echo "Current version will be bumped to: $(NEXT_PATCH)"
	poetry version $(NEXT_PATCH)
	git add pyproject.toml
	git commit -m "build: bump version to $(NEXT_PATCH)"
	git tag -a "v$(NEXT_PATCH)" -m "Release version $(NEXT_PATCH)"

## Bump the version using dunamai
.PHONY: bump-dynamic
bump-dynamic:
	@echo "Current dynamic version is: $(NEXT_DYNAMIC)"
	poetry version $(NEXT_DYNAMIC)

## Lock poetry dependencies
.PHONY: lock
lock:
	@echo "Locking poetry dependencies"
	export PYTHON_KEYRING_BACKEND=keyring.backends.fail.Keyring
	poetry lock -v

## Install the project using poetry
.PHONY: install
install: pyproject.toml poetry.lock
	export PYTHON_KEYRING_BACKEND=keyring.backends.fail.Keyring
	poetry install

## Publish the project to PyPI
.PHONY: publish
publish: pyproject.toml poetry.lock
	poetry publish --build

## Create Snakemake DAG figure for README
.PHONY: pipeline_figure
pipeline_figure: config.yml
	touch sample_sequences.fasta;
	snakemake --configfile config.yml --filegraph | dot -Tpng > img/pipeline.png;
	snakemake --configfile config.yml --filegraph | dot -Tsvg > img/pipeline.svg;
	rm sample_sequences.fasta;

## Run Ruff formatting and linting
.PHONY: format
format:
	@echo "Formatting project with ruff"
	poetry run ruff check --fix .
	poetry run ruff format .

## Run type checking, testing, and coverage reports
.PHONY: test
test:
	@echo "Running tests and type checks"
	mypy philharmonic --ignore-missing-imports
	pytest --cov=philharmonic --cov-report=term-missing tests

.DEFAULT_GOAL := help

define PRINT_HELP_PYSCRIPT
import re, sys; \
lines = '\n'.join([line for line in sys.stdin]); \
matches = re.findall(r'\n## (.*)\n[\s\S]+?\n([a-zA-Z_-]+):', lines); \
print(matches); \
print('Available rules:\n'); \
print('\n'.join(['{:25}{}'.format(*reversed(match)) for match in matches]))
endef
export PRINT_HELP_PYSCRIPT

## Show this help message
help:
	@$(PYTHON_INTERPRETER) -c "${PRINT_HELP_PYSCRIPT}" < $(MAKEFILE_LIST)
