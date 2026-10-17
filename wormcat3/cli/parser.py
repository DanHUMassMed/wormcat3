"""
wormcat3.cli.parser
~~~~~~~~~~~~~~~~~~~
Argparse CLI definition supporting modern subcommands and global flags.
"""

import argparse

import wormcat3.constants as cs
from wormcat3 import __version__


def build_parser() -> argparse.ArgumentParser:
    """Build and return the top-level CLI argument parser."""
    parser = argparse.ArgumentParser(
        prog="wormcat",
        description="WormCat 3: C. elegans Gene Set Enrichment Analysis & Visualization",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )

    # Global options
    parser.add_argument(
        "-V",
        "--version",
        action="version",
        version=f"%(prog)s {__version__}",
    )
    parser.add_argument(
        "-v",
        "--verbose",
        action="store_true",
        help="Enable verbose debug logging",
    )
    parser.add_argument(
        "-q",
        "--quiet",
        action="store_true",
        help="Silence standard logging output",
    )
    parser.add_argument(
        "--log-file",
        type=str,
        default=None,
        help="File path to record log output",
    )

    subparsers = parser.add_subparsers(
        dest="subcommand",
        title="Commands",
        metavar="<command>",
    )

    # --- enrich subcommand ---
    p_enrich = subparsers.add_parser(
        "enrich",
        help="Run enrichment analysis and visualization on a gene set",
        description="Run Fisher's exact enrichment test on a gene set and generate visualizations.",
    )
    p_enrich.add_argument(
        "-i",
        "--input",
        required=True,
        help="Input gene list file (CSV or plain text)",
    )
    p_enrich.add_argument(
        "-b",
        "--background",
        default=None,
        help="Custom background gene list file",
    )
    p_enrich.add_argument(
        "-a",
        "--annotation",
        default=cs.DEFAULT_ANNOTATION_FILE_NAME,
        help="Annotation dataset name or path",
    )
    p_enrich.add_argument(
        "-m",
        "--method",
        choices=["bonferroni", "fdr_bh"],
        default="bonferroni",
        help="P-value adjustment method",
    )
    p_enrich.add_argument(
        "-t",
        "--threshold",
        type=float,
        default=cs.DEFAULT_P_ADJUST_THRESHOLD,
        help="Significance threshold",
    )
    p_enrich.add_argument(
        "-g",
        "--gene-type",
        choices=[cs.GENE_TYPE_WORMBASE_ID, cs.GENE_TYPE_SEQUENCE_ID],
        default=None,
        help="Gene ID format (auto-detected if not specified)",
    )
    p_enrich.add_argument(
        "-o",
        "--output-dir",
        default=cs.DEFAULT_WORKING_DIR_PATH,
        help="Output destination folder",
    )
    p_enrich.add_argument(
        "--title",
        default=cs.DEFAULT_RUN_PREFIX,
        help="Run title identifier",
    )
    p_enrich.add_argument(
        "--email",
        default=None,
        help="User email for metadata tracking",
    )
    p_enrich.add_argument(
        "--no-viz",
        action="store_true",
        help="Skip bubble and sunburst chart creation",
    )

    # --- batch subcommand ---
    p_batch = subparsers.add_parser(
        "batch",
        help="Batch process multi-sheet Excel or folder of CSVs",
        description="Execute batch enrichment runs and compile a summary Excel report.",
    )
    p_batch.add_argument(
        "-i",
        "--input",
        required=True,
        help="Path to multi-sheet Excel file (.xlsx, .xls) or directory of CSVs",
    )
    p_batch.add_argument(
        "-b",
        "--background",
        default=None,
        help="Custom background gene list file",
    )
    p_batch.add_argument(
        "-a",
        "--annotation",
        default=cs.DEFAULT_ANNOTATION_FILE_NAME,
        help="Annotation dataset name or path",
    )
    p_batch.add_argument(
        "-m",
        "--method",
        choices=["bonferroni", "fdr_bh"],
        default="bonferroni",
        help="P-value adjustment method",
    )
    p_batch.add_argument(
        "-t",
        "--threshold",
        type=float,
        default=cs.DEFAULT_P_ADJUST_THRESHOLD,
        help="Significance threshold",
    )
    p_batch.add_argument(
        "-g",
        "--gene-type",
        choices=[cs.GENE_TYPE_WORMBASE_ID, cs.GENE_TYPE_SEQUENCE_ID],
        default=None,
        help="Gene ID format (auto-detected if not specified)",
    )
    p_batch.add_argument(
        "-o",
        "--output-dir",
        default=cs.DEFAULT_WORKING_DIR_PATH,
        help="Output destination folder",
    )
    p_batch.add_argument(
        "--title",
        default="batch_run",
        help="Batch run identifier",
    )
    p_batch.add_argument(
        "--email",
        default=None,
        help="User email for metadata tracking",
    )

    # --- gsea subcommand ---
    p_gsea = subparsers.add_parser(
        "gsea",
        help="Pre-ranked GSEA analysis from DESeq2 input",
        description="Run pre-ranked Gene Set Enrichment Analysis on differential expression results.",
    )
    p_gsea.add_argument(
        "-i",
        "--input",
        required=True,
        help="DESeq2 results CSV file (ID, log2FoldChange, pvalue)",
    )
    p_gsea.add_argument(
        "-a",
        "--annotation",
        default=cs.DEFAULT_ANNOTATION_FILE_NAME,
        help="Annotation dataset name or path",
    )
    p_gsea.add_argument(
        "-o",
        "--output-dir",
        default=cs.DEFAULT_WORKING_DIR_PATH,
        help="Output destination directory",
    )
    p_gsea.add_argument(
        "--title",
        default="gsea_run",
        help="Run title identifier",
    )
    p_gsea.add_argument(
        "--email",
        default=None,
        help="User email for metadata tracking",
    )
    p_gsea.add_argument(
        "--min-size",
        type=int,
        default=15,
        help="Minimum gene set size",
    )
    p_gsea.add_argument(
        "--max-size",
        type=int,
        default=500,
        help="Maximum gene set size",
    )
    p_gsea.add_argument(
        "--permutations",
        type=int,
        default=1000,
        help="Number of permutations",
    )
    p_gsea.add_argument(
        "--weight",
        type=float,
        default=1.0,
        help="Enrichment weighting score",
    )
    p_gsea.add_argument(
        "--seed",
        type=int,
        default=123,
        help="Random seed for reproducibility",
    )
    p_gsea.add_argument(
        "--threads",
        type=int,
        default=4,
        help="Number of parallel worker threads",
    )

    # --- gmt subcommand ---
    p_gmt = subparsers.add_parser(
        "gmt",
        help="Export WormCat annotations to Broad GMT format",
        description="Generate Category 1, 2, and 3 Broad GMT files from WormCat annotations.",
    )
    p_gmt.add_argument(
        "-o",
        "--output-dir",
        required=True,
        help="Directory to save generated .gmt files",
    )
    p_gmt.add_argument(
        "-a",
        "--annotation",
        default=cs.DEFAULT_ANNOTATION_FILE_NAME,
        help="Annotation dataset name or path",
    )
    p_gmt.add_argument(
        "-g",
        "--gene-type",
        choices=[cs.GENE_TYPE_WORMBASE_ID, cs.GENE_TYPE_SEQUENCE_ID],
        default=cs.GENE_TYPE_WORMBASE_ID,
        help="Gene ID column to export",
    )
    p_gmt.add_argument(
        "-p",
        "--prefix",
        default="wormcat",
        help="Filename prefix for generated GMT files",
    )

    # --- annotations subcommand ---
    p_annot = subparsers.add_parser(
        "annotations",
        help="Inspect available annotations",
        description="Inspect and list available WormCat annotation datasets.",
    )
    annot_sub = p_annot.add_subparsers(
        dest="annot_action",
        metavar="<action>",
    )
    annot_sub.add_parser(
        "list",
        help="List all available annotation CSV datasets",
    )

    # --- extract subcommand ---
    p_extract = subparsers.add_parser(
        "extract",
        help="Extract worksheets from an Excel workbook into CSV files",
        description="Extract all sheets from a multi-tab Excel workbook into separate CSV files.",
    )
    p_extract.add_argument(
        "-i",
        "--input",
        required=True,
        help="Path to input Excel file (.xlsx, .xls, .xlsm)",
    )
    p_extract.add_argument(
        "-o",
        "--output-dir",
        required=True,
        help="Directory to store extracted CSV files",
    )

    return parser
