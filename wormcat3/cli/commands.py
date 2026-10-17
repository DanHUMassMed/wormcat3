"""
wormcat3.cli.commands
~~~~~~~~~~~~~~~~~~~~~
Application service handlers executing CLI operations.
"""

from wormcat3.annotations_manger import AnnotationsManager
from wormcat3.cli.models import (
    BatchConfig,
    EnrichmentConfig,
    ExtractSheetsConfig,
    GMTExportConfig,
    GSEAConfig,
)
from wormcat3.logger import get_logger
from wormcat3.wormcat import Wormcat
from wormcat3.wormcat_excel import WormcatExcel

logger = get_logger(__name__)


def handle_enrich(config: EnrichmentConfig) -> int:
    """Execute single gene set enrichment analysis."""
    wc = Wormcat(
        working_dir_path=str(config.output_dir),
        title=config.title,
        annotation_file_name=config.annotation_file,
        email=config.email,
    )
    if config.no_viz:
        wc.perform_enrichment_analysis(
            gene_set_input=str(config.input_file),
            background_input=str(config.background_file) if config.background_file else None,
            p_adjust_method=config.method,
            p_adjust_threshold=config.threshold,
            gene_type=config.gene_type,
        )
    else:
        wc.analyze_and_visualize_enrichment(
            gene_set_input=str(config.input_file),
            background_input=str(config.background_file) if config.background_file else None,
            p_adjust_method=config.method,
            p_adjust_threshold=config.threshold,
            gene_type=config.gene_type,
        )
    return 0


def handle_batch(config: BatchConfig) -> int:
    """Execute batch gene set enrichment analysis."""
    wc = Wormcat(
        working_dir_path=str(config.output_dir),
        title=config.title,
        annotation_file_name=config.annotation_file,
        email=config.email,
    )
    wc.wormcat_batch(
        input_data=str(config.input_path),
        background_input=str(config.background_file) if config.background_file else None,
        p_adjust_method=config.method,
        p_adjust_threshold=config.threshold,
        gene_type=config.gene_type,
    )
    return 0


def handle_gsea(config: GSEAConfig) -> int:
    """Execute pre-ranked GSEA analysis on differential expression input."""
    wc = Wormcat(
        working_dir_path=str(config.output_dir),
        title=config.title,
        annotation_file_name=config.annotation_file,
        email=config.email,
    )
    wc.perform_gsea_analysis(deseq2_input=str(config.input_file))
    return 0


def handle_gmt(config: GMTExportConfig) -> int:
    """Export annotation datasets to Broad GMT format."""
    manager = AnnotationsManager(config.annotation_file)
    manager.create_gmt_for_annotations(
        output_dir_path=str(config.output_dir),
        id_col_nm=config.gene_type,
        output_file_nm_prefix=config.prefix,
    )
    return 0


def handle_annotations_list() -> int:
    """List all available annotation files."""
    files = AnnotationsManager.available_annotation_files()
    print("Available WormCat Annotation Files:")
    for f in files:
        print(f"  - {f}")
    return 0


def handle_extract(config: ExtractSheetsConfig) -> int:
    """Extract sheets from an Excel workbook into individual CSV files."""
    WormcatExcel.extract_csv_files(str(config.input_file), str(config.output_dir))
    print(f"Extracted CSV sheets to {config.output_dir}")
    return 0
