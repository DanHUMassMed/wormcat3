"""
tests/test_cli.py
~~~~~~~~~~~~~~~~~
Unit and integration tests for wormcat3 CLI.
"""

from pathlib import Path

import pytest
from pydantic import ValidationError

from wormcat3 import __version__
from wormcat3.cli.main import main
from wormcat3.cli.models import (
    BatchConfig,
    EnrichmentConfig,
    ExtractSheetsConfig,
    GSEAConfig,
)
from wormcat3.cli.parser import build_parser
from wormcat3.constants import PAdjustMethod


def test_cli_version_flag(capsys: pytest.CaptureFixture[str]) -> None:
    """Verify that -V and --version exit cleanly and display the current version."""
    parser = build_parser()
    with pytest.raises(SystemExit) as exc_info:
        parser.parse_args(["--version"])
    assert exc_info.value.code == 0
    captured = capsys.readouterr()
    assert f"wormcat {__version__}" in captured.out or f"wormcat {__version__}" in captured.err


def test_cli_help_flag(capsys: pytest.CaptureFixture[str]) -> None:
    """Verify that --help exits cleanly and lists subcommands."""
    parser = build_parser()
    with pytest.raises(SystemExit) as exc_info:
        parser.parse_args(["--help"])
    assert exc_info.value.code == 0
    captured = capsys.readouterr()
    assert "enrich" in captured.out
    assert "batch" in captured.out
    assert "gsea" in captured.out
    assert "gmt" in captured.out
    assert "annotations" in captured.out
    assert "extract" in captured.out


def test_cli_no_args_prints_help(capsys: pytest.CaptureFixture[str]) -> None:
    """Verify running wormcat without subcommands prints help to stderr and returns 1."""
    code = main([])
    assert code == 1
    captured = capsys.readouterr()
    assert "usage: wormcat" in captured.err


def test_cli_annotations_list(capsys: pytest.CaptureFixture[str]) -> None:
    """Verify 'wormcat annotations list' lists installed annotation datasets."""
    code = main(["annotations", "list"])
    assert code == 0
    captured = capsys.readouterr()
    assert "Available WormCat Annotation Files:" in captured.out
    assert "whole_genome_v2_nov-11-2021.csv" in captured.out


def test_validation_enrich_missing_input_file() -> None:
    """Verify ValidationError when input file does not exist."""
    with pytest.raises(ValidationError) as exc_info:
        EnrichmentConfig(input_file=Path("nonexistent_file_xyz.csv"))
    assert "Input gene list file not found" in str(exc_info.value)


def test_validation_enrich_threshold_out_of_bounds(tmp_path: Path) -> None:
    """Verify ValidationError when p_adjust threshold is outside (0, 1]."""
    test_file = tmp_path / "genes.csv"
    test_file.write_text("Gene\nWBGene00000001\nWBGene00000002\n")

    with pytest.raises(ValidationError):
        EnrichmentConfig(input_file=test_file, threshold=0.0)

    with pytest.raises(ValidationError):
        EnrichmentConfig(input_file=test_file, threshold=1.5)


def test_validation_enrich_method_coercion(tmp_path: Path) -> None:
    """Verify PAdjustMethod coercion from string."""
    test_file = tmp_path / "genes.csv"
    test_file.write_text("Gene\nWBGene00000001\nWBGene00000002\n")

    cfg1 = EnrichmentConfig(input_file=test_file, method="bonferroni")  # type: ignore[arg-type]
    assert cfg1.method == PAdjustMethod.BONFERRONI

    cfg2 = EnrichmentConfig(input_file=test_file, method="fdr_bh")  # type: ignore[arg-type]
    assert cfg2.method == PAdjustMethod.FDR


def test_validation_batch_invalid_path() -> None:
    """Verify ValidationError when batch input path does not exist."""
    with pytest.raises(ValidationError) as exc_info:
        BatchConfig(input_path=Path("nonexistent_path_xyz"))
    assert "Batch input path does not exist" in str(exc_info.value)


def test_validation_batch_invalid_file_extension(tmp_path: Path) -> None:
    """Verify ValidationError when a file passed to batch is not an Excel workbook."""
    txt_file = tmp_path / "invalid.txt"
    txt_file.write_text("some content")

    with pytest.raises(ValidationError) as exc_info:
        BatchConfig(input_path=txt_file)
    assert "Batch input file must be an Excel workbook" in str(exc_info.value)


def test_validation_gsea_missing_file() -> None:
    """Verify ValidationError when GSEA input file does not exist."""
    with pytest.raises(ValidationError) as exc_info:
        GSEAConfig(input_file=Path("nonexistent_deseq2.csv"))
    assert "DESeq2 input file not found" in str(exc_info.value)


def test_validation_extract_invalid_file(tmp_path: Path) -> None:
    """Verify ValidationError when extract input is not an Excel file."""
    txt_file = tmp_path / "invalid.txt"
    txt_file.write_text("some content")

    with pytest.raises(ValidationError) as exc_info:
        ExtractSheetsConfig(input_file=txt_file, output_dir=tmp_path / "out")
    assert "Input file must be an Excel workbook" in str(exc_info.value)


def test_cli_extract_sheets_integration(tmp_path: Path) -> None:
    """Verify 'wormcat extract -i <excel> -o <dir>' extracts sheets."""
    excel_path = Path("example_data/Murphy_TS.xlsx")
    if not excel_path.is_file():
        pytest.skip("example_data/Murphy_TS.xlsx not found")

    out_dir = tmp_path / "extracted_csvs"
    code = main(["extract", "-i", str(excel_path), "-o", str(out_dir)])
    assert code == 0
    extracted_files = list(out_dir.glob("*.csv"))
    assert len(extracted_files) > 0


def test_cli_gmt_integration(tmp_path: Path) -> None:
    """Verify 'wormcat gmt -o <dir>' produces GMT files."""
    out_dir = tmp_path / "gmt_out"
    code = main(["gmt", "-o", str(out_dir), "-p", "test_wormcat"])
    assert code == 0
    assert (out_dir / "test_wormcat_cat_1.gmt").is_file()
    assert (out_dir / "test_wormcat_cat_2.gmt").is_file()
    assert (out_dir / "test_wormcat_cat_3.gmt").is_file()


def test_cli_enrich_no_viz_integration(tmp_path: Path) -> None:
    """Verify 'wormcat enrich -i <csv> -o <dir> --no-viz' executes enrichment without plots."""
    sample_file = Path("example_data/DESeq2_Output/DESeq2-EC-ES/EC-ES-Up.csv")
    if not sample_file.is_file():
        pytest.skip("Sample file not found")

    out_dir = tmp_path / "enrich_out"
    code = main(
        [
            "enrich",
            "-i",
            str(sample_file),
            "-o",
            str(out_dir),
            "--title",
            "test_enrich",
            "--no-viz",
            "-m",
            "bonferroni",
            "-t",
            "0.05",
        ]
    )
    assert code == 0
    # Find the run directory created inside out_dir
    run_dirs = [d for d in out_dir.iterdir() if d.is_dir()]
    assert len(run_dirs) == 1
    run_dir = run_dirs[0]
    # Check that CSV outputs were generated
    csv_outputs = list(run_dir.glob("*.csv"))
    assert len(csv_outputs) > 0
    # Check that run params were recorded
    run_params = list(run_dir.glob("run_params_*.txt"))
    assert len(run_params) == 1


def test_cli_main_catches_validation_error() -> None:
    """Verify main returns exit code 1 when a ValidationError occurs."""
    code = main(["enrich", "-i", "nonexistent_file_12345.csv"])
    assert code == 1


def test_cli_main_catches_unknown_command(capsys: pytest.CaptureFixture[str]) -> None:
    """Verify main exits with error code 2 on invalid command syntax."""
    with pytest.raises(SystemExit) as exc_info:
        main(["unknown_subcommand"])
    assert exc_info.value.code == 2
