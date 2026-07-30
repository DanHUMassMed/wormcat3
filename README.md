# WormCat 3 (`wormcat3`)

[![PyPI version](https://img.shields.io/pypi/v/wormcat3.svg?color=blue)](https://pypi.org/project/wormcat3/)
[![Python Version](https://img.shields.io/badge/python-3.13%2B-blue.svg)](https://www.python.org/downloads/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Code Style: Ruff](https://img.shields.io/badge/code%20style-ruff-000000.svg)](https://github.com/astral-sh/ruff)
[![Type Checked: Mypy](https://img.shields.io/badge/mypy-checked-blue.svg)](https://mypy-lang.org/)

> **An advanced Python tool for annotating, analyzing, and visualizing gene set enrichment data from *C. elegans* microarray, RNA-seq, or RNAi screen experiments.**

---

## Table of Contents

- [Overview & Key Features](#overview--key-features)
- [Original Publication & Abstract](#original-publication--abstract)
- [Online Web Tool & Workflow Diagrams](#online-web-tool--workflow-diagrams)
- [Installation](#installation)
- [Quick Start & Usage Examples](#quick-start--usage-examples)
  - [1. Gene Set Enrichment Analysis & Visualization](#1-gene-set-enrichment-analysis--visualization)
  - [2. Pre-ranked GSEA Analysis](#2-pre-ranked-gsea-analysis)
  - [3. Batch Execution (Excel / Directory of CSVs)](#3-batch-execution-excel--directory-of-csvs)
- [Development](#development)
- [License & Citation](#license--citation)

---

## Overview & Key Features

`wormcat3` is the modern Python implementation of **WormCat**, a computational framework for *C. elegans* functional gene annotation and enrichment analysis. Unlike standard Gene Ontology (GO) tools where up to 30% of *C. elegans* genes lack representation, WormCat provides near-complete annotation using a **3-level nested category strategy** (Cat1, Cat2, Cat3).

### Key Features
- 🧬 **Near-Complete *C. elegans* Annotation**: Annotated across 3 nested hierarchical levels for broad (Cat1) to fine-grained (Cat3) functional insights.
- 📊 **Enrichment & Statistical Testing**: Supports Fisher's Exact Test / Hypergeometric testing with multiple testing adjustments (Bonferroni, Benjamini-Hochberg FDR).
- 📈 **Pre-ranked GSEA**: Built-in support for GSEA analysis on DESeq2 or custom pre-ranked differential expression inputs.
- 🎨 **Rich Visualizations**: Automatically generates interactive HTML Sunburst charts and high-resolution Bubble charts.
- 📁 **Batch Processing & Excel Reports**: Single-command batch execution for multi-sheet Excel files or CSV folders with formatted Excel export.
- ⚡ **High Performance & Modern Python Stack**: Built for Python 3.13+, managed seamlessly with `uv` and `pandas` 2.x / `plotnine`.

---

## Original Publication & Abstract

### Analysis of genome-scale data with WormCat identifies novel enriched gene categories in studies from metabolic, tissue-specific, and lifespan-drug data

#### Authors: Amy Holdorf, Daniel Higgins, Anne Hart, Peter Boag, Gregory Pazour, Marian Walhout,and Amy Walker

[GENETICS February 1, 2020 vol. 214 no. 2 279-294;](https://academic.oup.com/genetics/article/214/2/279/5930455)

### Abstract
The emergence of large sets of gene regulation data has revealed the need for improved tools to 1) identify enriched functional gene categories and 2) visualize enrichment patterns across comparative datasets.  Gene ontogeny enrichment (GO) has several limitations for C. elegans analysis. First, around 30% of C. elegans genes are not represented in commonly used search engines. Second, it is difficult to compare multiple GO analyses. To allow visualization and categorization of C. elegans gene sets, we have developed a web-based tool, WormCat.  This tool uses a near complete annotation of C. elegans genes to determine category enrichment and define potential co-regulated or co-functioning gene sets. Then WormCat provides a scaled heat map for visualization along with enrichment statistics and annotation of each input gene. We have developed an annotation strategy based on a nested category approach where each gene is annotated at three levels.  Enrichment scores are generated at each level, allowing both broad (Cat1) and more detailed analysis (Cat2, Cat3).  Using WormCat on published RNA seq datasets from metabolic, tissue-specific or after treatment with lifespan-increasing drugs, we show that WormCat finds major categories appearing in GO searches and also identifies additional enriched categories that are informative for interpreting phenotypes or predicting biological function.  Thus, WormCat is a powerful tool that will allow a sophisticated analysis of gene enrichment in different types of C. elegans datasets.

---

## Online Web Tool & Workflow Diagrams

## Overview Wormcat
Wormcat is also available as an online tool at [www.wormcat.com](http://www.wormcat.com); the online version greatly simplifies the use of Wormcat and is maintained by the [Walker Lab at UMASS Medical School](https://amywalkerlab.com/).

###### The diagram below shows the flow of the Wormcat process:
<img src="http://wormcat.com/static/images/WormCat-Flow.png" alt="Flow" width="700"/>

###### The diagrams below shows sample output from a Wormcat.com run:
<img src="http://wormcat.com/static/images/results_screen.png" alt="Results" width="700"/>

##### Starburst view of categorical data
<img src="http://wormcat.com/static/images/sunburst.png" alt="starburst" width="400"/>

---

## Installation

### Using `uv` (Recommended)

```bash
uv add wormcat3
# Or in a virtual environment:
uv pip install wormcat3
```

### Using `pip`

```bash
pip install wormcat3
```

**Requirements**: Python `>= 3.13`.

---

## Quick Start & Usage Examples

### 1. Gene Set Enrichment Analysis & Visualization

Run enrichment analysis on a list of gene identifiers (WormBase IDs or Sequence Names) and automatically generate plots:

```python
from wormcat3 import Wormcat, PAdjustMethod

# Initialize Wormcat runner
wc = Wormcat(title="my_experiment", email="user@example.com")

# Execute enrichment analysis and create visual charts
wc.analyze_and_visualize_enrichment(
    gene_set_input="path/to/gene_list.csv",  # or a python list of gene IDs
    p_adjust_method=PAdjustMethod.BONFERRONI,
    p_adjust_threshold=0.05,
)
```

### 2. Pre-ranked GSEA Analysis

Perform pre-ranked Gene Set Enrichment Analysis (GSEA) on DESeq2 differential expression results:

```python
from wormcat3 import Wormcat

wc = Wormcat(title="gsea_experiment")
wc.perform_gsea_analysis(deseq2_input="path/to/deseq2_results.csv")
```

### 3. Batch Execution (Excel / Directory of CSVs)

Process multi-tab Excel files or entire directories of CSV gene lists at once, auto-generating combined Excel reports:

```python
from wormcat3 import Wormcat

wc = Wormcat(title="batch_run")
wc.wormcat_batch(input_data="path/to/multi_sheet_data.xlsx")
```

---

## Development

The project uses [`uv`](https://github.com/astral-sh/uv) and `make` for dependency management and developer workflows.

### Developer Setup

```bash
# Clone the repository
git clone https://github.com/DanHUMassMed/wormcat3.git
cd wormcat3

# Install dependencies and bootstrap virtual environment
make install
```

### Common Commands

- **Run Tests**: `make test`
- **Lint & Type Check**: `make lint`
- **Format Code**: `make format`
- **Launch Dev Notebook**: `make dev`
- **Build Distribution Packages**: `make build`

---

## License & Citation

### License
This project is licensed under the [MIT License](LICENSE).

### Citation
If you use `wormcat3` or WormCat in your research, please cite the original publication:

> Holdorf AD, Higgins DP, Hart AC, Boag PR, Pazour GJ, Walhout AJM, Walker AK. **Analysis of genome-scale data with WormCat identifies novel enriched gene categories in studies from metabolic, tissue-specific, and lifespan-drug data.** *Genetics*. 2020 Feb;214(2):279-294. doi: [10.1534/genetics.119.302919](https://academic.oup.com/genetics/article/214/2/279/5930455).
