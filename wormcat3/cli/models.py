"""
wormcat3.cli.models
~~~~~~~~~~~~~~~~~~~
Type-safe Pydantic v2 schemas for CLI input validation.
"""

from pathlib import Path
from typing import Any, Literal, Optional

from pydantic import BaseModel, ConfigDict, Field, field_validator

from wormcat3.constants import (
    DEFAULT_ANNOTATION_FILE_NAME,
    DEFAULT_P_ADJUST_THRESHOLD,
    DEFAULT_RUN_PREFIX,
    DEFAULT_WORKING_DIR_PATH,
    PAdjustMethod,
)


class BaseRunConfig(BaseModel):
    """Base configuration shared across execution commands."""

    model_config = ConfigDict(arbitrary_types_allowed=True)

    output_dir: Path = Field(default_factory=lambda: Path(DEFAULT_WORKING_DIR_PATH))
    title: str = Field(default=DEFAULT_RUN_PREFIX)
    annotation_file: str = Field(default=DEFAULT_ANNOTATION_FILE_NAME)
    email: Optional[str] = None
    verbose: bool = False
    quiet: bool = False
    log_file: Optional[Path] = None


class EnrichmentConfig(BaseRunConfig):
    """Configuration schema for single gene set enrichment analysis."""

    input_file: Path
    background_file: Optional[Path] = None
    method: PAdjustMethod = Field(default=PAdjustMethod.BONFERRONI)
    threshold: float = Field(default=DEFAULT_P_ADJUST_THRESHOLD, gt=0.0, le=1.0)
    gene_type: Optional[Literal["Wormbase.ID", "Sequence.ID"]] = None
    no_viz: bool = False

    @field_validator("method", mode="before")
    @classmethod
    def coerce_method(cls, v: Any) -> PAdjustMethod:
        if isinstance(v, PAdjustMethod):
            return v
        if isinstance(v, str):
            return PAdjustMethod.from_str(v)
        raise ValueError(f"Invalid PAdjustMethod: {v}")

    @field_validator("input_file")
    @classmethod
    def validate_input_file(cls, v: Path) -> Path:
        if not v.is_file():
            raise ValueError(f"Input gene list file not found: {v}")
        return v

    @field_validator("background_file")
    @classmethod
    def validate_background_file(cls, v: Optional[Path]) -> Optional[Path]:
        if v is not None and not v.is_file():
            raise ValueError(f"Background gene list file not found: {v}")
        return v


class BatchConfig(BaseRunConfig):
    """Configuration schema for batch gene set enrichment analysis."""

    input_path: Path
    background_file: Optional[Path] = None
    method: PAdjustMethod = Field(default=PAdjustMethod.BONFERRONI)
    threshold: float = Field(default=DEFAULT_P_ADJUST_THRESHOLD, gt=0.0, le=1.0)
    gene_type: Optional[Literal["Wormbase.ID", "Sequence.ID"]] = None

    @field_validator("method", mode="before")
    @classmethod
    def coerce_method(cls, v: Any) -> PAdjustMethod:
        if isinstance(v, PAdjustMethod):
            return v
        if isinstance(v, str):
            return PAdjustMethod.from_str(v)
        raise ValueError(f"Invalid PAdjustMethod: {v}")

    @field_validator("input_path")
    @classmethod
    def validate_input_path(cls, v: Path) -> Path:
        if not v.exists():
            raise ValueError(f"Batch input path does not exist: {v}")
        if v.is_file() and v.suffix.lower() not in [".xlsx", ".xls", ".xlsm"]:
            raise ValueError(f"Batch input file must be an Excel workbook (.xlsx, .xls, .xlsm): {v}")
        return v

    @field_validator("background_file")
    @classmethod
    def validate_background_file(cls, v: Optional[Path]) -> Optional[Path]:
        if v is not None and not v.is_file():
            raise ValueError(f"Background gene list file not found: {v}")
        return v


class GSEAConfig(BaseRunConfig):
    """Configuration schema for pre-ranked GSEA analysis."""

    input_file: Path
    min_size: int = Field(default=15, ge=1)
    max_size: int = Field(default=500, ge=1)
    permutation_num: int = Field(default=1000, ge=10)
    weight: float = Field(default=1.0, ge=0.0)
    seed: int = Field(default=123)
    threads: int = Field(default=4, ge=1)

    @field_validator("input_file")
    @classmethod
    def validate_input_file(cls, v: Path) -> Path:
        if not v.is_file():
            raise ValueError(f"DESeq2 input file not found: {v}")
        return v


class GMTExportConfig(BaseModel):
    """Configuration schema for exporting annotations to GMT format."""

    output_dir: Path
    annotation_file: str = Field(default=DEFAULT_ANNOTATION_FILE_NAME)
    gene_type: Literal["Wormbase.ID", "Sequence.ID"] = "Wormbase.ID"
    prefix: str = "wormcat"


class ExtractSheetsConfig(BaseModel):
    """Configuration schema for extracting worksheets from an Excel workbook."""

    input_file: Path
    output_dir: Path

    @field_validator("input_file")
    @classmethod
    def validate_input_file(cls, v: Path) -> Path:
        if not v.is_file():
            raise ValueError(f"Input Excel file not found: {v}")
        if v.suffix.lower() not in [".xlsx", ".xls", ".xlsm"]:
            raise ValueError(f"Input file must be an Excel workbook (.xlsx, .xls, .xlsm): {v}")
        return v
