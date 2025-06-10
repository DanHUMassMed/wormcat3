import os
import pandas as pd
from datetime import datetime
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
from wormcat3.wormcat_error import WormcatError, ErrorCode

class Wormcat:
    """
    Main class that coordinates file handling, annotation management,
    and statistical analysis for gene enrichment.
    """
    
    def __init__(self, 
                 working_dir_path = cs.DEFAULT_WORKING_DIR_PATH, 
                 title = cs.DEFAULT_RUN_PREFIX, 
                 annotation_file_name = cs.DEFAULT_ANNOTATION_FILE_NAME,
                 email = None):
        """Initialize Wormcat with working directory and annotation file."""

        self.email = email
        self.title = title
        ### Create the working directory 
        ### It is possible that when we generate the hash we have a conflict with a current output directory
        ### If there is a conflict try again
        is_validate_directory = False
        while(not is_validate_directory):
            self.run_number = file_util.generate_5_digit_hash(prefix=title + "_")
            working_dir_path = Path(working_dir_path) / self.run_number
            is_validate_directory, self.working_dir_path = file_util.validate_directory_path(working_dir_path, validation_indicator=True)
        
        # Setup annotation manager
        self.annotation_manager = AnnotationsManager(annotation_file_name)

    def _run_params(self, params: dict):
        # Define the output file name
        output_path = Path(self.working_dir_path) / f"run_params_{self.run_number}.txt"
        
        if "Timestamp" not in params:
            params["Timestamp"] = datetime.now().strftime("%Y-%m-%d %H:%M:%S")

        # Write the dictionary manually as "key: value"
        with open(output_path, 'w') as file:
            for key, value in params.items():
                file.write(f"{key}: {value}\n")
        
    def perform_gsea_analysis(self, deseq2_input: Union[str, pd.DataFrame]):
        
        if isinstance(deseq2_input, str):
            deseq2_df = file_util.read_deseq2_file(deseq2_input)
        else:
            deseq2_df = deseq2_input

        
        gsea_analyzer = GSEAAnalyzer(self.working_dir_path)
        removed_rows_df, deseq2_df = gsea_analyzer.clean_input_data(deseq2_df)
        
        # Save the removed rows
        if not removed_rows_df.empty:            
            removed_file_name = f"genes_removed_from_analysis_{self.run_number}.csv"
            removed_path = Path(self.working_dir_path) / removed_file_name
            removed_rows_df.to_csv(removed_path, index=False)
            
        ranked_list_df = gsea_analyzer.create_ranked_list(deseq2_df)

        first_3_ids = ranked_list_df['Gene'].head(3).tolist()
        gene_type = self.annotation_manager.get_gene_id_type(first_3_ids)

        for category in [1,2,3]:
            gmt_format = self.annotation_manager.category_to_gmt_format(category, id_col_nm=gene_type)
            results_name = f"gsea_category_{category}_{self.run_number}"
            results_df = gsea_analyzer.run_preranked_gsea(ranked_list_df , gmt_format, results_name)
            # Save the results_df
            gsea_category_path = Path(self.working_dir_path) / f"{results_name}.csv"
            results_df.to_csv(gsea_category_path, index=False)
        
        run_params = {
            "Email":self.email,
            "Title":self.title,
            "Annotation File":self.annotation_manager.annotation_name(),
            "Gene Type":gene_type
            }
        self._run_params(run_params)
        
        print(f"Analysis complete. Output can be found at {self.working_dir_path}")

        
    def perform_enrichment_analysis(
            self, 
            gene_set_input: Union[str, list], 
            background_input: Union[str, list] = None, 
            *, 
            p_adjust_method = PAdjustMethod.BONFERRONI, 
            p_adjust_threshold = cs.DEFAULT_P_ADJUST_THRESHOLD,
            gene_type = None
            
        )-> List[Dict[str, pd.DataFrame]]:
        """Perform enrichment test on the gene set."""
        
        if isinstance(gene_set_input, (str, Path)):
            gene_set_list = file_util.read_gene_set_file(Path(gene_set_input))
        elif isinstance(gene_set_input, list):
            gene_set_list = gene_set_input
        else:
            raise WormcatError(
                "Invalid type: gene_set_input must be a file name or a list",
                ErrorCode.INVALID_TYPE.to_dict()
                )

        if background_input is not None:
            if isinstance(background_input, (str, Path)):
                background_list = file_util.read_gene_set_file(Path(background_input))
            elif isinstance(background_input, list):
                background_list = background_input
            else:
                raise WormcatError(
                    "Invalid type: background_list must be a file name or a list",
                    ErrorCode.INVALID_TYPE.to_dict()
                    )
        else:
            background_list = None
        
        
        # Preprocess gene set list
        gene_set_list = self.annotation_manager.dedup_list(gene_set_list)
        
        if gene_type is None:
            gene_type = self.annotation_manager.get_gene_id_type(gene_set_list)
            
        if gene_type is not cs.GENE_TYPE_WORMBASE_ID or cs.GENE_TYPE_SEQUENCE_ID:
            raise WormcatError(
                "Invalid Gene Type: gene_type must be {cs.GENE_TYPE_WORMBASE_ID} or {cs.GENE_TYPE_SEQUENCE_ID}.",
                ErrorCode.INVALID_TYPE.to_dict()
                )
            
        # Add annotations
        gene_set_and_categories_df, genes_not_matched_df = self.annotation_manager.segment_genes_by_annotation_match(gene_set_list, gene_type)
        
        # Save the annotated input gene set
        rgs_and_categories_path = Path(self.working_dir_path) / f"input_annotated_{self.run_number}.csv"
        gene_set_and_categories_df.to_csv(rgs_and_categories_path, index=False)
        
        if not genes_not_matched_df.empty:
                genes_not_annotated_path = Path(self.working_dir_path) / f"genes_not_annotated_{self.run_number}.csv"
                genes_not_matched_df.to_csv(genes_not_annotated_path, index=False)


        # Preprocess background list
        if background_list is not None:  
            background_list = self.annotation_manager.dedup_list(background_list)
            
            background_type = self.annotation_manager.get_gene_id_type(background_list)
            if background_type != gene_type:
                raise WormcatError("Gene Set Type and Background Type MUST be the same. {gene_type}!={background_type}", 
                                    ErrorCode.CONSTRAINT_VIOLATION.to_dict())
                
            background_df, background_not_annotated_df = self.annotation_manager.segment_genes_by_annotation_match(background_list, background_type)

            # Save the annotated background input
            background_annotated_path = Path(self.working_dir_path) / f"background_annotated_{self.run_number}.csv"
            background_df.to_csv(background_annotated_path, index=False)

            if not background_not_annotated_df.empty:
                background_not_annotated_path = Path(self.working_dir_path) / f"background_not_annotated_{self.run_number}.csv"
                background_not_annotated_df.to_csv(background_not_annotated_path, index=False)
        else:
            # If no background is provided we use the whole genome  
            background_df = self.annotation_manager.annotations_df
        
        
        # Setup statistical analyzer
        self.analyzer = EnrichmentAnalyzer(
            background_df, 
            self.working_dir_path,
            self.run_number
        )
        
        run_params = {
            "Email":self.email,
            "Title":self.title,
            "Annotation File":self.annotation_manager.annotation_name(),
            "Significance Method":p_adjust_method.value,
            "Significance Threshold":p_adjust_threshold,
            "Domain Scope":"All Genes" if background_list is None else "Custom Background",
            "Domain Gene Set Count":f"{len(background_df):,}",
            "Input Gene Set Count":f"{len(gene_set_and_categories_df):,}",
            "Gene Type":gene_type
            }
        self._run_params(run_params)
        
        # Run enrichment analysis
        return self.analyzer.calculate_category_enrichment_scores(
            gene_set_and_categories_df,
            p_adjust_method=p_adjust_method,
            p_adjust_threshold=p_adjust_threshold
        )

    def analyze_and_visualize_enrichment(self,
            gene_set_input: Union[str, list], 
            background_input: Union[str, list] = None, 
            *, 
            p_adjust_method = PAdjustMethod.BONFERRONI, 
            p_adjust_threshold = cs.DEFAULT_P_ADJUST_THRESHOLD,
            gene_type = None):
        
        test_results = self.perform_enrichment_analysis(gene_set_input, background_input, 
                                                        p_adjust_method=p_adjust_method, p_adjust_threshold=p_adjust_threshold, gene_type=gene_type)
        for test_result in test_results:
            result_file_path, result_df = next(iter(test_result.items()))
            data_file_nm = os.path.basename(result_file_path)
            base_dir_path = os.path.dirname(result_file_path)
            plot_title = data_file_nm[:-10]
            plot_title = plot_title.replace("_", " ")
            plot_title = plot_title[:1].upper() + plot_title[1:]
            create_bubble_chart(base_dir_path, data_file_nm, plot_title = plot_title)
            
        run_number = os.path.basename(base_dir_path)
        create_sunburst(base_dir_path,run_number)
        print(f"Analysis complete. Output can be found at {self.working_dir_path}")
        
    def wormcat_batch(self,
            input_data: str, 
            background_input: Union[str, list] = None, 
            *, 
            p_adjust_method = PAdjustMethod.BONFERRONI, 
            p_adjust_threshold = cs.DEFAULT_P_ADJUST_THRESHOLD,
            gene_type = None):
        
        input_path = Path(input_data)
        
        # Check if path exists
        if not input_path.exists():
            raise WormcatError(f"Path not found: {input_data}", ErrorCode.FILE_NOT_FOUND.to_dict())
        
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
                wormcat = Wormcat(working_dir_path = self.working_dir_path, 
                                  annotation_file_name = self.annotation_manager.annotation_file_path, 
                                  title = file.stem)
                wormcat.analyze_and_visualize_enrichment(str(file), background_input, 
                                                         p_adjust_method = p_adjust_method, 
                                                         p_adjust_threshold = p_adjust_threshold,
                                                         gene_type = gene_type)
        else:
            print(f"Directory doesn't contain any CSV files: {input_path}")
            return 
        

        annotation_file_path = self.annotation_manager.annotation_file_path
        wormcat_excel = WormcatExcel()
        working_dir_path = Path(self.working_dir_path)
        wormcat_excel.create_summary_spreadsheet(self.working_dir_path, annotation_file_path, f"{working_dir_path}/{working_dir_path.stem}.xlsx")
        print(f"Analysis complete. Output can be found at {self.working_dir_path}")
