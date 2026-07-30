.PHONY: help install dev lint format test bump-version build deploy clean

# Default shell
SHELL := /bin/bash

# Package metadata
PACKAGE_NAME := wormcat3
INIT_FILE := wormcat3/__init__.py
PYPROJECT_FILE := pyproject.toml

.DEFAULT_GOAL := help

help: ## Display available commands
	@echo "Usage: make [target]"
	@echo ""
	@echo "Targets:"
	@awk 'BEGIN {FS = ":.*?## "} /^[a-zA-Z_-]+:.*?## / {printf "  \033[36m%-15s\033[0m %s\n", $$1, $$2}' $(MAKEFILE_LIST)

install: ## Bootstrap environment and install all dependencies via uv
	@echo "==> Bootstrapping virtual environment and dependencies using uv..."
	uv sync --all-groups
	uv pip install -e .

dev: ## Start local development environment / Jupyter notebook
	@echo "==> Starting local development environment..."
	uv run jupyter notebook notebooks/start_here.ipynb

lint: ## Execute static analysis and formatting checks (ruff, mypy)
	@echo "==> Running ruff check..."
	uv run ruff check $(PACKAGE_NAME) tests
	@echo "==> Running ruff format check..."
	uv run ruff format --check $(PACKAGE_NAME) tests
	@echo "==> Running mypy type checking..."
	uv run mypy $(PACKAGE_NAME)

format: ## Format code and fix auto-fixable lint issues
	@echo "==> Formatting code with ruff..."
	uv run ruff format $(PACKAGE_NAME) tests
	uv run ruff check --fix $(PACKAGE_NAME) tests

test: ## Execute test suite via pytest
	@echo "==> Running test suite..."
	uv run pytest

bump-version: ## Bump patch version (or specify explicit VERSION=X.Y.Z)
	@if [ -n "$(VERSION)" ]; then \
		NEW_VERSION="$(VERSION)"; \
		echo "==> Bumping version to explicit: $$NEW_VERSION"; \
	else \
		CURRENT_VERSION=$$(grep "__version__ =" $(INIT_FILE) | cut -d'"' -f2); \
		IFS='.' read -r -a PARTS <<< "$$CURRENT_VERSION"; \
		MAJOR=$${PARTS[0]}; \
		MINOR=$${PARTS[1]}; \
		PATCH=$${PARTS[2]}; \
		NEW_VERSION="$$MAJOR.$$MINOR.$$((PATCH + 1))"; \
		echo "==> Auto-bumping patch version: $$CURRENT_VERSION → $$NEW_VERSION"; \
	fi; \
	sed -i '' "s/__version__ = .*/__version__ = \"$$NEW_VERSION\"/" $(INIT_FILE); \
	sed -i '' "s/^version = .*/version = \"$$NEW_VERSION\"/" $(PYPROJECT_FILE)

build: clean ## Build distribution packages (wheel and sdist)
	@echo "==> Building distribution packages..."
	uv run python -m build

deploy: bump-version build ## Bump version, build dist, validate, upload to PyPI & push git tags
	@echo "==> Validating distribution artifacts with twine..."
	uv run twine check dist/*
	@CURRENT_VERSION=$$(grep "__version__ =" $(INIT_FILE) | cut -d'"' -f2); \
	echo "==> Uploading $(PACKAGE_NAME) v$$CURRENT_VERSION to PyPI..."; \
	uv run twine upload --repository pypi dist/*; \
	git add $(INIT_FILE) $(PYPROJECT_FILE); \
	git commit -m "Bump version to v$$CURRENT_VERSION"; \
	git tag "v$$CURRENT_VERSION"; \
	git push && git push --tags

clean: ## Remove build, distribution, and cache artifacts
	@echo "==> Cleaning build and cache artifacts..."
	rm -rf build/ dist/ *.egg-info .pytest_cache .ruff_cache .mypy_cache
	find . -type d -name "__pycache__" -exec rm -rf {} +
	find . -type f -name "*.pyc" -delete
