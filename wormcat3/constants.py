from enum import Enum

# Enum for p-value adjustment methods
class PAdjustMethod(Enum):
    BONFERRONI = 'bonferroni'
    FDR = 'fdr_bh'
    
    @classmethod
    def from_str(cls, value: str) -> "PAdjustMethod":
        for method in cls:
            if method.value == value:
                return method
        raise ValueError(f"Invalid PAdjustMethod: {value}")

# Wormcat Configuration
DEFAULT_WORKING_DIR_PATH = "./wormcat_out"
DEFAULT_RUN_PREFIX = "run"

# Annotations Management Configuration
DEFAULT_P_ADJUST_THRESHOLD = 0.05
DEFAULT_ANNOTATION_FILE_NAME = "whole_genome_v2_nov-11-2021.csv"

# Gene Set Enrichment Analysis
DEFAULT_GSEA_RESULTS_DIR = "./gsea_results"

GENE_TYPE_WORMBASE_ID = "Wormbase.ID"
GENE_TYPE_SEQUENCE_ID = "Sequence.ID"

# Bubble Chart Configuration
DEFAULT_TITLE = "RGS"
DEFAULT_WIDTH = 6
DEFAULT_HEIGHT = 5.5
