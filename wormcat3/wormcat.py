import os
import pandas as pd
from pathlib import Path
from typing import Union, List, Dict
from wormcat3 import file_util
from wormcat3.annotations_manger import AnnotationsManager
from wormcat3.statistical_analysis import EnrichmentAnalyzer
from wormcat3.gsea_analyzer import GSEAAnalyzer
from wormcat3.constants import PAdjustMethod
from wormcat3.bubble_chart import create_bubble_chart
from wormcat3.sunburst import create_sunburst
from wormcat3.wormcat_excel import WormcatExcel
import wormcat3.constants as cs

class Wormcat:
    """
    Main class that coordinates file handling, annotation management,
    and statistical analysis for gene enrichment.
    """
    
    def __init__(self, 
                 working_dir_path = cs.DEFAULT_WORKING_DIR_PATH, 
                 run_prefix = cs.DEFAULT_RUN_PREFIX, 
                 annotation_file_name = cs.DEFAULT_ANNOTATION_FILE_NAME):
        """Initialize Wormcat with working directory and annotation file."""
        
        ### Create the working directory 
        self.run_number = file_util.generate_5_digit_hash(prefix=run_prefix + "_")
        working_dir_path = Path(working_dir_path) / self.run_number
        self.working_dir_path = file_util.validate_directory_path(working_dir_path)
        
        # Setup annotation manager
        self.annotation_manager = AnnotationsManager(annotation_file_name)


    def perform_gsea_analysis(self, deseq2_input: Union[str, pd.DataFrame]):
        
        if isinstance(deseq2_input, str):
            deseq2_df = file_util.read_deseq2_file(deseq2_input)
        else:
            deseq2_df = deseq2_input

        gsea_analyzer = GSEAAnalyzer(self.working_dir_path)
        ranked_list_df = gsea_analyzer.create_ranked_list(deseq2_df)

        for category in [1,2,3]:
            gmt_format = self.annotation_manager.category_to_gmt_format(category)
            results_name = f"gsea_category_{category}_{self.run_number}"
            results_df = gsea_analyzer.run_preranked_gsea(ranked_list_df , gmt_format, results_name)
            # Save the results_df
            gsea_category_path = Path(self.working_dir_path) / f"{results_name}.csv"
            results_df.to_csv(gsea_category_path, index=False)

        
    def perform_enrichment_analysis(
            self, 
            gene_set_input: Union[str, list], 
            background_input: Union[str, list] = None, 
            *, 
            p_adjust_method = PAdjustMethod.BONFERRONI, 
            p_adjust_threshold = cs.DEFAULT_P_ADJUST_THRESHOLD
            
        )-> List[Dict[str, pd.DataFrame]]:
        """Perform enrichment test on the gene set."""
        
        if isinstance(gene_set_input, str):
            gene_set_list = file_util.read_gene_set_file(gene_set_input)
        else:
            gene_set_list = gene_set_input
        
        if isinstance(background_input, str):
            background_list = file_util.read_gene_set_file(background_input)
        else:
            background_list = background_input

        if not isinstance(p_adjust_method, PAdjustMethod):
            raise ValueError(f"Invalid p_adjust_method: {p_adjust_method}. Must be a valid PAdjustMethod.")

        assert 0 < p_adjust_threshold <= 1, "p_adjust_threshold must be between 0 and 1 (exclusive lower, inclusive upper)."

        
        # Preprocess gene set list
        gene_set_list = self.annotation_manager.dedup_list(gene_set_list)
        gene_type = self.annotation_manager.get_gene_id_type(gene_set_list)
        
        # Add annotations
        gene_set_and_categories_df = self.annotation_manager.add_annotations(gene_set_list, gene_type)
        
        # Save the annotated input gene set
        rgs_and_categories_path = Path(self.working_dir_path) / f"input_annotated_{self.run_number}.csv"
        gene_set_and_categories_df.to_csv(rgs_and_categories_path, index=False)


        # Preprocess background list
        if background_list is not None:  
            background_list = self.annotation_manager.dedup_list(background_list)
            background_type = self.annotation_manager.get_gene_id_type(background_list)
            if background_type != gene_type:
                raise ValueError("Gene Set Type and Background Type MUST be the same. {gene_type}!={background_type}")
            background_df, background_not_matched_df = self.annotation_manager.split_background_on_annotation_match(background_list, background_type)

            # Save the annotated background input
            background_annotated_path = Path(self.working_dir_path) / f"background_annotated_{self.run_number}.csv"
            background_df.to_csv(background_annotated_path, index=False)

            if not background_not_matched_df.empty:
                background_not_matched_path = Path(self.working_dir_path) / "background_not_matched.csv"
                background_not_matched_df.to_csv(background_not_matched_path, index=False)
        else:
            # If no background is provided we use the whole genome  
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

    def analyze_and_visualize_enrichment(self,
            gene_set_input: Union[str, list], 
            background_input: Union[str, list] = None, 
            *, 
            p_adjust_method = PAdjustMethod.BONFERRONI, 
            p_adjust_threshold = cs.DEFAULT_P_ADJUST_THRESHOLD):
        
        test_results = self.perform_enrichment_analysis(gene_set_input, background_input, p_adjust_method=p_adjust_method, p_adjust_threshold=p_adjust_threshold)
        for test_result in test_results:
            result_file_path, result_df = next(iter(test_result.items()))
            data_file_nm = os.path.basename(result_file_path)
            base_dir_path = os.path.dirname(result_file_path)
            plot_title = data_file_nm[:-10]
            create_bubble_chart(base_dir_path, data_file_nm, plot_title = plot_title)
            
        run_number = os.path.basename(base_dir_path)
        create_sunburst(base_dir_path,run_number)
        
    def wormcat_batch(self,
            input_data: str, 
            background_input: Union[str, list] = None, 
            *, 
            p_adjust_method = PAdjustMethod.BONFERRONI, 
            p_adjust_threshold = cs.DEFAULT_P_ADJUST_THRESHOLD):
        
        input_path = Path(input_data)
        
        # Check if path exists
        if not input_path.exists():
            raise FileNotFoundError(f"Path not found: {input_data}")
        
        if input_path.is_file():
                # Check if it's an Excel file
                if input_path.suffix.lower() in ['.xlsx', '.xls', '.xlsm']:
                    try:
                        csv_file_path = Path(self.working_dir_path) /f"{input_path.stem}_CSVs"
                        WormcatExcel.extract_csv_files(input_data, csv_file_path)
                    except Exception as e:
                        print(f"Invalid Excel file: {input_path}. Error: {str(e)}")
                        return
                else:
                    print(f"File is not an Excel file: {input_path}")
                    return
            
        # Check if it's a directory
        elif input_path.is_dir():
            csv_file_path = input_path
        else:
            print(f"input_data is neither a valid Excel file nor a directory with CSV files: {input_data}")
            return
                    
        # Look for CSV files
        csv_files = list(csv_file_path.glob('*.csv'))  
        if csv_files:
            for file in csv_files:
                wormcat = Wormcat(working_dir_path=self.working_dir_path,run_prefix=file.stem)
                wormcat.analyze_and_visualize_enrichment(str(file), background_input, p_adjust_method = p_adjust_method, p_adjust_threshold = p_adjust_threshold)
        else:
            print(f"Directory doesn't contain any CSV files: {input_path}")
            return 
        

        annotation_file_path = self.annotation_manager.annotation_file_path
        wormcat_excel = WormcatExcel()
        working_dir_path = Path(self.working_dir_path)
        wormcat_excel.create_summary_spreadsheet(self.working_dir_path, annotation_file_path, f"{working_dir_path}/{working_dir_path.stem}.xlsx")