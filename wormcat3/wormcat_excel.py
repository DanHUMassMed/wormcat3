"""
Create the summary Excel file based on the individual Wormcat runs
"""
import os
import pandas as pd
from pathlib import Path
from typing import List, Dict, Union, Any, Optional, Tuple, Set
import warnings
from wormcat3 import file_util

warnings.filterwarnings('ignore', category=UserWarning, module='openpyxl')
warnings.simplefilter(action='ignore', category=FutureWarning)


class ExcelConversionError(Exception):
    """Custom exception for Excel to CSV conversion errors"""
    pass


class CategoryProcessingError(Exception):
    """Custom exception for errors during category processing"""
    pass


class WormcatExcel:
    """Class to handle Wormcat Excel processing and summary creation"""
    
    # Class constants
    DEFAULT_CATEGORIES = [1, 2, 3]
    DEFAULT_P_VALUE_THRESHOLDS = [0.0000000001, 0.00000001, 0.000001, 0.0001, 0.05]
    DEFAULT_COLOR_FORMATS = [
        {'bg_color': 'white', 'font_color': 'black', 'num_format': '0'},
        {'bg_color': '#244162', 'font_color': 'white', 'num_format': '0.000E+00'},
        {'bg_color': '#355f91', 'font_color': 'white', 'num_format': '0.000E+00'},
        {'bg_color': '#95b3d7', 'font_color': 'black', 'num_format': '0.000E+00'},
        {'bg_color': '#b8cce4', 'font_color': 'black', 'num_format': '0.000E+00'},
        {'bg_color': '#f4f2fe', 'font_color': 'black', 'num_format': '0.000E+00'}
    ]
    
    def __init__(self, 
                 significance_threshold: float = 0.05,
                 p_value_thresholds: Optional[List[float]] = None):
        """
        Initialize the WormcatExcel processor with customizable parameters.
        
        Args:
            significance_threshold: P-value threshold for statistical significance (default: 0.05)
            p_value_thresholds: Thresholds for p-value color formatting (default: see DEFAULT_P_VALUE_THRESHOLDS)
        """
        self.significance_threshold = significance_threshold
        self.p_value_thresholds = p_value_thresholds or self.DEFAULT_P_VALUE_THRESHOLDS
        self.categories = self.DEFAULT_CATEGORIES
        self.color_formats = self.DEFAULT_COLOR_FORMATS
        self.excel_formats = None
        
    # File operations methods
    @staticmethod
    def extract_csv_files(excel_path: str, csv_file_path: str) -> Dict[str, str]:
        """
        Create CSV files from the given Excel spreadsheet.
        
        Args:
            excel_path: Path to the Excel file
            csv_file_path: Directory to save CSV files
            
        Returns:
            Dictionary mapping sheet names to output CSV file paths
            
        Raises:
            FileNotFoundError: If the Excel file doesn't exist
            ExcelConversionError: If there's an error processing the Excel file
        """
        # Validate inputs
        if not os.path.exists(excel_path):
            raise FileNotFoundError(f"Excel file not found at path: {excel_path}")
        
        csv_file_path = file_util.validate_directory_path(csv_file_path)
        
        # Process the Excel file
        try:
            input_excel = pd.ExcelFile(excel_path)
        except ValueError:
            raise ExcelConversionError(f"File [{excel_path}] is not a valid Excel file.")
        except Exception as e:
            raise ExcelConversionError(f"Failed to open Excel file: {str(e)}")

        # Process each sheet
        results = {}
        try:
            for sheet in input_excel.sheet_names:
                sheet_df = input_excel.parse(sheet)
                output_file = os.path.join(csv_file_path, f"{sheet}.csv")
                sheet_df.to_csv(output_file, index=False)
                results[sheet] = output_file
                
            return results
        except Exception as e:
            raise ExcelConversionError(f"Failed during conversion process: {str(e)}")
        
    def create_summary_spreadsheet(self, wormcat_out_path: str, 
                                  annotation_file: str, 
                                  out_xsl_file_nm: str) -> None:
        """
        After all the wormcat runs have been executed, create a summary Excel spreadsheet.
        This function collects data from all category output files and creates a summary.
        
        Args:
            wormcat_out_path: Path to the directory containing wormcat output folders
            annotation_file: Path to the annotation CSV file
            out_xsl_file_nm: Path for the output Excel file
            
        Raises:
            FileNotFoundError: If input directories don't exist
            Exception: For other processing errors
        """
        wormcat_path = Path(wormcat_out_path)
        if not wormcat_path.exists():
            raise FileNotFoundError(f"Wormcat output path not found: {wormcat_out_path}")
        
        process_lst = self._collect_category_files(wormcat_path)
            
        if not process_lst:
            print("No valid category files found to process")
            return
            
        df_process = pd.DataFrame(process_lst, columns=['sheet', 'category', 'file', 'label'])
        self._process_category_files(df_process, annotation_file, out_xsl_file_nm)



    def _create_category_summary(self, data: pd.DataFrame, category_name: str) -> pd.DataFrame:
        """
        Create a summary DataFrame for a specific category.
        
        Args:
            data: The input annotation DataFrame
            category_name: Name of the category column to summarize
            
        Returns:
            DataFrame with category values and their counts, sorted by category name
        """
        if category_name not in data.columns:
            raise ValueError(f"Category '{category_name}' not found in data columns")
            
        category = data[category_name].value_counts()
        category = pd.DataFrame({
            category_name: category.index, 
            'Count': category.values
        })
        category = category.sort_values(by=[category_name])
        return category

    # Data processing methods
    def _significant(self, value: float) -> Union[float, str]:
        """
        Convert p-values to significance indicators.
        
        Args:
            value: P-value to evaluate
        
        Returns:
            Original value if significant, 'NS' if not significant, 'NV' if not a value
           
        Notes:
            NV = Not a Value
            NS = Not Significant
        """
        if pd.isna(value):
            return 'NV'
        return value if value < self.significance_threshold else 'NS'

    def _process_category_file_row(self, row: pd.Series, sheet: pd.DataFrame) -> pd.DataFrame:
        """
        Process a single category file and merge with the existing sheet DataFrame.
        
        Args:
            row: Series containing file information (file path, category, label)
            sheet: Existing DataFrame to merge results into
            
        Returns:
            Updated DataFrame with merged results
            
        Raises:
            FileNotFoundError: If the category file doesn't exist
            CategoryProcessingError: If there's an error processing the category file
        """
        file_name = row['file']
        if not Path(file_name).exists():
            raise FileNotFoundError(f"Category file not found: {file_name}")
            
        try:
            label_category = f"Category {row['category']}"
            label_pvalue = f"{row['label']}_PValue"
            label_rgs = f"{row['label']}_RGS"

            cat_results = pd.read_csv(file_name)
            cat_results.rename(columns={
                'Category': label_category,
                'RGS': label_rgs, 
                'PValue': label_pvalue
            }, inplace=True)
            
            # Clean up unnecessary columns
            columns_to_drop = []
            for col in ['Unnamed: 0', 'AC']:
                if col in cat_results.columns:
                    columns_to_drop.append(col)
                    
            if columns_to_drop:
                cat_results.drop(columns_to_drop, axis=1, inplace=True)

            # Merge with existing sheet
            sheet = pd.merge(sheet, cat_results, on=label_category, how='outer')
            
            # Apply significance threshold and RGS formatting
            sheet[label_pvalue] = sheet[label_pvalue].apply(self._significant)
            sheet[label_rgs] = sheet[label_rgs].apply(lambda x: x if pd.notna(x) and x > 0 else 0)
            
            return sheet
        except Exception as e:
            raise CategoryProcessingError(f"Error processing file {file_name}: {str(e)}")

    # Excel formatting methods
    def _get_excel_formats(self, writer: pd.ExcelWriter) -> List[Any]:
        """
        Generate Excel format configurations for conditional formatting.
        
        Args:
            writer: Excel writer object
            
        Returns:
            List of Excel format objects
        """
        return [writer.book.add_format(fmt) for fmt in self.color_formats]

    def _create_legend(self, writer: pd.ExcelWriter) -> None:
        """
        Creates a simple sheet page as a Legend in the Excel file.
        
        Args:
            writer: Excel writer object
        """
        data = {'Color Code <=': self.p_value_thresholds[::-1]}
        legend_sheet = pd.DataFrame(data)

        legend_sheet.to_excel(writer, sheet_name='Legend', index=False)
        worksheet = writer.sheets['Legend']
        num_rows, num_columns = legend_sheet.shape
        sheet_range = f"A1:{chr(num_columns + 64)}{num_rows+1}"

        # Apply conditional formatting
        worksheet.conditional_format(sheet_range, {
            'type': 'cell', 
            'criteria': '=', 
            'value': 0, 
            'format': self.excel_formats[0]
        })
        
        for index, value in enumerate(self.p_value_thresholds):
            worksheet.conditional_format(sheet_range, {
                'type': 'cell', 
                'criteria': '<=', 
                'value': value, 
                'format': self.excel_formats[index+1]
            })

        worksheet.autofit()

    def _apply_conditional_formatting(self, worksheet, sheet_range: str) -> None:
        """
        Apply conditional formatting to a worksheet range.
        
        Args:
            worksheet: Excel worksheet object
            sheet_range: Cell range to format (e.g. "B1:Z100")
        """
        print(f"{worksheet=} {sheet_range=}")
        # Format for zero values
        worksheet.conditional_format(sheet_range, {
            'type': 'cell', 
            'criteria': '=', 
            'value': 0, 
            'format': self.excel_formats[0]
        })
        
        # Format for p-value thresholds
        for index, value in enumerate(self.p_value_thresholds):
            worksheet.conditional_format(sheet_range, {
                'type': 'cell', 
                'criteria': '<=', 
                'value': value, 
                'format': self.excel_formats[index+1]
            })

    # Core processing methods
    def _process_category_files(self, files_to_process: pd.DataFrame, 
                              annotation_file: str, 
                              out_data_xlsx: str) -> None:
        """
        Processes each category file and creates the corresponding Excel summary sheets.
        
        Args:
            files_to_process: DataFrame with file information (sheet, category, file, label)
            annotation_file: Path to the annotation CSV file
            out_data_xlsx: Path for the output Excel file
            
        Raises:
            FileNotFoundError: If input files don't exist
            Exception: For other processing errors
        """
        annotation_path = Path(annotation_file)
        if not annotation_path.exists():
            raise FileNotFoundError(f"Annotation file not found: {annotation_file}")
        
        try:
            # Load annotation data
            data = pd.read_csv(annotation_file)
            
            # Create Excel writer
            with pd.ExcelWriter(out_data_xlsx, engine='xlsxwriter') as writer:
                sheets = files_to_process['sheet'].unique()
                self.excel_formats = self._get_excel_formats(writer)

                # Create legend sheet
                self._create_legend(writer)

                # Process each sheet
                for sheet_label in sheets:
                    self._process_sheet(files_to_process, data, writer, sheet_label)
            
        except Exception as e:
            print(f"Error processing category files: {str(e)}")
            raise


    def _get_excel_col_name(self, n: int) -> str:
        """Convert 1-based column number to Excel column letter."""
        result = ""
        while n > 0:
            n, remainder = divmod(n - 1, 26)
            result = chr(65 + remainder) + result
        return result

    def _process_sheet(self, files_to_process: pd.DataFrame, 
                      data: pd.DataFrame, 
                      writer: pd.ExcelWriter, 
                      sheet_label: str) -> None:
        """
        Process a single sheet in the output Excel file.
        
        Args:
            files_to_process: DataFrame with file information
            data: Annotation data
            writer: Excel writer object
            sheet_label: Name of the sheet to process
        """
        cat_files = files_to_process[files_to_process['sheet'] == sheet_label]
        if cat_files.empty:
            print(f"No files to process for sheet {sheet_label}")
            return
                
        label_category = f"Category {cat_files['category'].iloc[0]}"
        
        try:
            # Create the initial summary sheet
            category_sheet = self._create_category_summary(data, label_category)
            
            # Process each file for this category
            for _, row in cat_files.iterrows():
                try:
                    file_path = Path(row['file'])
                    if not file_path.exists():
                        print(f"File not found: {row['file']}")
                        continue
                    
                    category_sheet = self._process_category_file_row(row, category_sheet)
                except (FileNotFoundError, CategoryProcessingError) as e:
                    print(str(e))
                    continue

            # Write the sheet to Excel
            category_sheet.to_excel(writer, sheet_name=sheet_label, index=False)
            worksheet = writer.sheets[sheet_label]
            
            # Apply formatting
            num_rows, num_columns = category_sheet.shape
            #sheet_range = f"B1:{chr(num_columns + 64)}{num_rows+1}"
            col_name = self._get_excel_col_name(num_columns)
            sheet_range = f"B1:{col_name}{num_rows+1}"
            self._apply_conditional_formatting(worksheet, sheet_range)
            
            worksheet.autofit()
            
        except Exception as e:
            print(f"Error processing sheet {sheet_label}: {str(e)}")
            raise

    def _collect_category_files(self, wormcat_path: Path) -> List[Dict[str, Any]]:
        """
        Collect all category files from the wormcat output directory.
        
        Args:
            wormcat_path: Path to the wormcat output directory
            
        Returns:
            List of dictionaries with file information
        """
        process_lst = []
        
        # Process each directory in the wormcat output path
        for dir_item in wormcat_path.iterdir():
            if not dir_item.is_dir():
                continue
                
            dir_nm = dir_item.name
            
            # Process each category
            for cat_num in self.categories:
                rgs_fisher_path = dir_item / f"category_{cat_num}_fisher_{dir_nm}.csv"
                
                # Only add files that exist
                if rgs_fisher_path.exists():
                    cat_nm = f"Cat{cat_num}"
                    row = {
                        'sheet': cat_nm, 
                        'category': cat_num,
                        'file': str(rgs_fisher_path), 
                        'label': dir_nm
                    }
                    process_lst.append(row)
                else:
                    print(f"File not found: {rgs_fisher_path}")
                    
        return process_lst


if __name__ == "__main__":
    # Example usage
    
    # Create a WormcatExcel instance with custom settings
    wormcat_excel = WormcatExcel(
        significance_threshold=0.01,  # Custom threshold
    )
    
    # Create a summary spreadsheet
    # wormcat_excel.create_summary_spreadsheet(
    #     "/path/to/wormcat_output", 
    #     "/path/to/annotation.csv", 
    #     "/path/to/output.xlsx"
    # )