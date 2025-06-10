import pandas as pd
import os
from pathlib import Path
from wormcat3 import file_util
import wormcat3.constants as cs
from wormcat3.wormcat_error import WormcatError, ErrorCode
    

class AnnotationsManager:
    """ Manages gene annotations and preprocessing. """
    
    def __init__(self, annotation_file=cs.DEFAULT_ANNOTATION_FILE_NAME):
        """Initialize with the path to the annotation file."""
        if file_util.is_file_path(annotation_file):
            self.annotation_file_path = annotation_file
        else:
            self.annotation_file_path = file_util.find_file_path(annotation_file)
            if not self.annotation_file_path:
                raise WormcatError(f"Annotation file not found: {annotation_file}", ErrorCode.FILE_NOT_FOUND.to_dict())
                        
        self.annotations_df = self._load_annotations()
            
     
        
    def _load_annotations(self):
        """ Load annotations from file. """
        try:
            df = pd.read_csv(self.annotation_file_path)
            df.columns = df.columns.str.replace(' ', '.')
            if df.empty:
                raise WormcatError(f"Annotation file '{self.annotation_file_path}' is empty.", ErrorCode.FILE_LOAD_FAILED.to_dict())
            return df
        except Exception as e:
            raise WormcatError(f"Failed to load annotation file: {e}", 
                                ErrorCode.FILE_LOAD_FAILED.to_dict(), 
                                origin="AnnotationsManager._load_annotations")
            
    @staticmethod
    def available_annotation_files():
        """
        Return a list of .csv file names located in:
        - the 'extdata' directory (next to this script), and
        - the directory defined by the WORMCAT_DATA_PATH environment variable (if set).
        """
        file_path = Path(__file__).resolve()
        extdata_path = file_path.parent / "extdata"

        search_dirs = [extdata_path]

        env_path = os.environ.get("WORMCAT_DATA_PATH")
        if env_path:
            env_path = Path(env_path)
            if env_path.exists() and env_path.is_dir():
                search_dirs.append(env_path)

        found_files = set()

        for directory in search_dirs:
            if directory.exists():
                for f in directory.glob("*.csv"):
                    if f.is_file():
                        found_files.add(f.name)

        return sorted(found_files)

    
    def annotation_name(self):
        # Get the base name of the file (with extension, no directory)
        if isinstance(self.annotation_file_path, (str, os.PathLike)) and self.annotation_file_path:
            return os.path.basename(self.annotation_file_path)
        else:
            return "[NO_ANNOTATION_NAME]"

    
    def get_gene_id_type(self, gene_set):
        """ Determine the gene ID type from the gene set. """
        
        if len(gene_set) < 3:
            raise WormcatError("At least two genes are required for comparison.", ErrorCode.INVALID_VALUE.to_dict())
        
        # Check if the first two genes start with "WBGene"
        if gene_set[1].startswith("WBGene") and gene_set[2].startswith("WBGene"):
            return cs.GENE_TYPE_WORMBASE_ID
        elif not gene_set[1].startswith("WBGene") and not gene_set[2].startswith("WBGene"):
            return cs.GENE_TYPE_SEQUENCE_ID
        else:
            raise WormcatError("Invalid gene data: One gene starts with 'WBGene', but the other does not.", ErrorCode.INVALID_VALUE.to_dict())
    
    @staticmethod
    def dedup_list(input_list):
        """ Deduplicate a list while preserving order. """
        
        seen = set()
        deduped_list = []
        for item in input_list:
            if item not in seen:
                deduped_list.append(item)
                seen.add(item)
        return deduped_list
    
    def add_annotations(self, gene_set_list, gene_type):
        """ Add annotations to the gene set. """
        
        gene_set_df = pd.DataFrame(gene_set_list, columns=[gene_type])
        
        # Verify if 'gene_type' is a column in the DataFrame
        if gene_type not in self.annotations_df.columns:
            raise WormcatError(f"Column '{gene_type}' not found in the DataFrame.", ErrorCode.INVALID_VALUE.to_dict())
        
        return pd.merge(gene_set_df, self.annotations_df, on=gene_type, how='left')


    def segment_genes_by_annotation_match(self, gene_set_list, gene_type):
        """ Split genes into those with and without annotations. """
        
        gene_set_df = pd.DataFrame(gene_set_list, columns=[gene_type])
        
        # Check if gene_type is in both dataframes
        if gene_type not in gene_set_df.columns:
            raise WormcatError(f"'{gene_type}' MUST be a column in the Gene Set DataFrame.", ErrorCode.MISSING_FIELD.to_dict())
        if gene_type not in self.annotations_df.columns:
            raise WormcatError(f"Incorrect '{gene_type}' name. Column not found in {self.annotation_name()} file.", ErrorCode.MISSING_FIELD.to_dict())
        
        # Perform the left merge
        merged_df = pd.merge(gene_set_df, self.annotations_df, on=gene_type, how='left')
        
        # Split based on presence of annotation (assuming at least one non-key column in annotations_df)
        annotation_columns = [col for col in self.annotations_df.columns if col != gene_type]
        
        genes_matched_df = merged_df.dropna(subset=annotation_columns)
        genes_not_matched_df = merged_df[merged_df[annotation_columns].isnull().all(axis=1)]
        genes_not_matched_df = genes_not_matched_df[[gene_type]]

        return genes_matched_df, genes_not_matched_df

    def create_gmt_for_annotations(self, output_dir_path, id_col_nm="Wormbase.ID", output_file_nm_prefix="wormcat"):
        """ Create GMT formatted files for all categories. """
        for category in [1,2,3]:
            gmt_format = self.category_to_gmt_format(category, id_col_nm)
            output_file_path = f"{output_dir_path}/{output_file_nm_prefix}_cat_{category}.gmt"
            self._save_gmt_to_file(gmt_format, output_file_path)    

    def category_to_gmt_format(self, category, id_col_nm="Wormbase.ID"):
        """ Convert an annotation dataframe category to GMT format. """
        category_col = f"Category.{category}"
        id_col = "Function.ID"
                
        category_df = self.annotations_df[[id_col_nm, category_col]]
        category_df = category_df.rename(columns={category_col: id_col})
        
        # Validate that required columns exist
        required_cols = [id_col, id_col_nm]        
        missing_cols = [col for col in required_cols if col not in category_df.columns]
        if missing_cols:
            raise WormcatError(f"Incorrect column names: {', '.join(missing_cols)} {self.annotation_name()} file does not have these columns.", 
                                ErrorCode.MISSING_FIELD.to_dict())
        
        if category_df[id_col].isna().any():
            raise WormcatError(f"{self.annotation_name()} dataframe Column '{id_col}' contains NaN values", 
                                ErrorCode.INVALID_VALUE.to_dict())
        
        # Use ID column as description if desc_col is None
        grouped = category_df.groupby([id_col])[id_col_nm].apply(list).reset_index()
        
        if len(grouped) == 0:
            raise WormcatError("No gene sets found after grouping", ErrorCode.INVALID_VALUE.to_dict())
        
        # Create GMT formatted dictionary
        gmt_format = {}
        for _, row in grouped.iterrows():            
            # Filter out any None or NaN values from gene list
            gene_list = [str(gene) for gene in row[id_col_nm] if pd.notna(gene)]
            
            # Only include if there are genes in the set
            if gene_list:
                gmt_format[row[id_col]] = gene_list
                
        return gmt_format

    def _save_gmt_to_file(self, gmt_format, output_file_path='wormcat.gmt'):
        """ Write GMT formatted dictionary to disk. """
        # Ensure the output directory exists
        output_dir_path = os.path.dirname(output_file_path)
        file_util.validate_directory_path(output_dir_path, not_empty_check = False) 
        
        # Write to GMT file
        with open(output_file_path, 'w') as file:
            for gene_id, gene_list in gmt_format.items():
                description = gene_id
                line = f"{gene_id}\t{description}\t" + '\t'.join(gene_list)
                file.write(line + '\n')
        
        if not os.path.exists(output_file_path):
            print(f"[Warning] Output file {output_file_path} was not created.")
        elif os.path.getsize(output_file_path) == 0:
            print(f"[Warning] Output file {output_file_path} is empty.")
        else:
            print(f"Successfully created GMT file: {output_file_path}")
        
        return output_file_path
