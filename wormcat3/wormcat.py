from . import file_util
from .annotations_manger import AnnotationsManager
from .statistical_analysis import EnrichmentAnalyzer, PAdjustMethod
from pathlib import Path
from typing import Union

import os

    
class Wormcat:
    """
    Main class that coordinates file handling, annotation management,
    and statistical analysis for gene enrichment.
    """
    
    def __init__(self, working_dir_path="./wormcat_out", annotation_file_name="whole_genome_v2_nov-11-2021.csv"):
        """Initialize Wormcat with working directory and annotation file."""
        
        ### Create the working directory 
        self.run_number = file_util.generate_5_digit_hash(prefix="run_")
        working_dir_path = Path(working_dir_path) / self.run_number
        self.working_dir_path = file_util.validate_directory_path(working_dir_path)
        
        # Setup annotation manager
        self.annotation_manager = AnnotationsManager(annotation_file_name)
        

    
    def enrichment_test(self, gene_set_input: Union[str, list], background_list=None, *, p_adjust_method=PAdjustMethod.BONFERRONI, p_adjust_threshold=0.01):
        """Perform enrichment test on the gene set."""
        
        assert 0 < p_adjust_threshold <= 1, "p_adjust_threshold must be between 0 and 1 (exclusive lower, inclusive upper)."

        if not isinstance(p_adjust_method, PAdjustMethod):
            raise ValueError(f"Invalid p_adjust_method: {p_adjust_method}. Must be a valid PAdjustMethod.")


        if isinstance(gene_set_input, str):
            gene_set_list = file_util.read_gene_set_file(gene_set_input)
        else:
            gene_set_list = gene_set_input
        
        # Preprocess gene set list
        gene_set_list = self.annotation_manager.dedup_list(gene_set_list)
        gene_type = self.annotation_manager.get_gene_id_type(gene_set_list)
        
        # Add annotations
        gene_set_and_categories_df = self.annotation_manager.add_annotations(gene_set_list, gene_type)
        
        # Save intermediate results
        rgs_and_categories_path = Path(self.working_dir_path) / f"input_annotated_{self.run_number}.csv"
        gene_set_and_categories_df.to_csv(rgs_and_categories_path, index=False)


        # Preprocess background list
        if background_list is not None:
            background_list = self.annotation_manager.dedup_list(background_list)
            background_type = self.annotation_manager.get_gene_id_type(background_list)
            if background_type != gene_type:
                raise ValueError("Gene Set Type and Background Type MUST be the same. {gene_type}!={background_type}")
            background_df, background_not_matched_df = self.annotation_manager.split_background_on_annotation_match(background_list, background_type)

            if not background_not_matched_df.empty:
                background_not_matched_path = Path(self.working_dir_path) / "background_not_matched.csv"
                background_not_matched_df.to_csv(background_not_matched_path, index=False)
        else:
            background_df = self.annotation_manager.annotations_df
        
        
        # Setup statistical analyzer
        self.analyzer = EnrichmentAnalyzer(
            background_df, 
            self.working_dir_path,
            self.run_number
        )
        
        # Run enrichment analysis
        return self.analyzer.perform_enrichment_test(
            gene_set_and_categories_df,
            p_adjust_method=p_adjust_method,
            p_adjust_threshold=p_adjust_threshold
        )


    