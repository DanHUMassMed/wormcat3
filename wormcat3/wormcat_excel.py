"""
Create the summary Excel file based on the individual Wormcat runs
"""
import os
import pandas as pd
from pathlib import Path
from typing import List, Dict, Union, Any


def create_category_summary(data: pd.DataFrame, category_name: str) -> pd.DataFrame:
    """
    Create a summary DataFrame for a specific category.
    
    Args:
        data: The input annotation DataFrame
        category_name: Name of the category column to summarize
        
    Returns:
        DataFrame with category values and their counts, sorted by category name
    """
    category = data[category_name].value_counts()
    category = pd.DataFrame(
        {category_name: category.index, 'Count': category.values})
    category = category.sort_values(by=[category_name])
    return category


def process_category_file_row(row: pd.Series, sheet: pd.DataFrame) -> pd.DataFrame:
    """
    Process a single category file and merge with the existing sheet DataFrame.
    
    Args:
        row: Series containing file information (file path, category, label)
        sheet: Existing DataFrame to merge results into
        
    Returns:
        Updated DataFrame with merged results
    """
    try:
        file_name = row['file']
        label_category = f"Category {row['category']}"
        label_pvalue = f"{row['label']}_PValue"
        label_rgs = f"{row['label']}_RGS"

        cat_results = pd.read_csv(file_name)
        cat_results.rename(columns={'Category': label_category,
                           'RGS': label_rgs, 'PValue': label_pvalue}, inplace=True)
        
        # Only drop columns if they exist
        if 'Unnamed: 0' in cat_results.columns:
            cat_results.drop('Unnamed: 0', axis=1, inplace=True)
        if 'AC' in cat_results.columns:
            cat_results.drop('AC', axis=1, inplace=True)

        sheet = pd.merge(sheet, cat_results, on=label_category, how='outer')
        sheet[label_pvalue] = sheet[label_pvalue].apply(significant)
        sheet[label_rgs] = sheet[label_rgs].apply(lambda x: x if x > 0 else 0)
        return sheet
    except Exception as e:
        print(f"Error processing file {row['file']}: {str(e)}")
        return sheet


def significant(value: float) -> Union[float, str]:
    """
    Convert p-values to significance indicators.
    
    Args:
        value: P-value to evaluate
        
    Returns:
        Original value if p < 0.05, 'NS' if not significant, 'NV' if not a value
        
    Notes:
        NV = Not a Value
        NS = Not Significant
    """
    if pd.isna(value):
        return 'NV'
    return value if value < 0.05 else 'NS'


def create_legend(writer: pd.ExcelWriter, values: List[float], excel_formats: List[Any]) -> None:
    """
    Creates a simple sheet page as a Legend in the Excel file.
    
    Args:
        writer: Excel writer object
        values: List of threshold values for conditional formatting
        excel_formats: List of Excel format objects for different thresholds
    """
    data = {'Color Code <=': values[::-1]}
    legend_sheet = pd.DataFrame(data)

    legend_sheet.to_excel(writer, sheet_name='Legend', index=False)
    worksheet = writer.sheets['Legend']
    num_rows, num_columns = legend_sheet.shape
    sheet_range = f"A1:{chr(num_columns + 64)}{num_rows+1}"

    worksheet.conditional_format(sheet_range, {
                                 'type': 'cell', 'criteria': '=', 'value': 0, 'format': excel_formats[0]})
    for index, value in enumerate(values):
        worksheet.conditional_format(sheet_range, {
                                     'type': 'cell', 'criteria': '<=', 'value': value, 'format': excel_formats[index+1]})

    worksheet.autofit()


def get_excel_formats(writer):
    """
    Generate Excel format configurations for conditional formatting.
    
    Args:
        writer: Excel writer object
        
    Returns:
        List of Excel format objects
    """
    formats = [
        {'bg_color': 'white', 'font_color': 'black', 'num_format': '0'},
        {'bg_color': '#244162', 'font_color': 'white', 'num_format': '0.000E+00'},
        {'bg_color': '#355f91', 'font_color': 'white', 'num_format': '0.000E+00'},
        {'bg_color': '#95b3d7', 'font_color': 'black', 'num_format': '0.000E+00'},
        {'bg_color': '#b8cce4', 'font_color': 'black', 'num_format': '0.000E+00'},
        {'bg_color': '#f4f2fe', 'font_color': 'black', 'num_format': '0.000E+00'}
    ]
    return [writer.book.add_format(f) for f in formats]


def process_category_files(files_to_process: pd.DataFrame, annotation_file: str, out_data_xlsx: str) -> None:
    """
    Processes each category file and creates the corresponding Excel summary sheets.
    
    Args:
        files_to_process: DataFrame with file information (sheet, category, file, label)
        annotation_file: Path to the annotation CSV file
        out_data_xlsx: Path for the output Excel file
    """
    try:
        # Validate input file exists
        if not Path(annotation_file).exists():
            raise FileNotFoundError(f"Annotation file not found: {annotation_file}")
            
        data = pd.read_csv(annotation_file)
        writer = pd.ExcelWriter(out_data_xlsx, engine='xlsxwriter')

        sheets = files_to_process['sheet'].unique()
        values = [0.0000000001, 0.00000001, 0.000001, 0.0001, 0.05]
        excel_formats = get_excel_formats(writer)

        create_legend(writer, values, excel_formats)

        for sheet_label in sheets:
            cat_files = files_to_process[files_to_process['sheet'] == sheet_label]
            if cat_files.empty:
                continue
                
            label_category = f"Category {cat_files['category'].iloc[0]}"
            category_sheet = create_category_summary(data, label_category)
            
            for index, row in cat_files.iterrows():
                # Validate file exists
                if not Path(row['file']).exists():
                    print(f"Warning: File not found: {row['file']}")
                    continue
                category_sheet = process_category_file_row(row, category_sheet)

            category_sheet.to_excel(writer, sheet_name=sheet_label, index=False)
            worksheet = writer.sheets[sheet_label]
            num_rows, num_columns = category_sheet.shape
            sheet_range = f"B1:{chr(num_columns + 64)}{num_rows+1}"

            worksheet.conditional_format(sheet_range, {
                                         'type': 'cell', 'criteria': '=', 'value': 0, 'format': excel_formats[0]})
            for index, value in enumerate(values):
                worksheet.conditional_format(sheet_range, {
                                             'type': 'cell', 'criteria': '<=', 'value': value, 'format': excel_formats[index+1]})

            worksheet.autofit()

        writer.close()
        print(f"Successfully created Excel summary: {out_data_xlsx}")
        
    except Exception as e:
        print(f"Error processing category files: {str(e)}")


def create_summary_spreadsheet(wormcat_out_path: str, annotation_file: str, out_xsl_file_nm: str) -> None:
    """
    After all the wormcat runs have been executed, create a summary Excel spreadsheet.
    This function collects data from all category output files and creates a summary.
    
    Args:
        wormcat_out_path: Path to the directory containing wormcat output folders
        annotation_file: Path to the annotation CSV file
        out_xsl_file_nm: Path for the output Excel file
    """
    try:
        wormcat_path = Path(wormcat_out_path)
        if not wormcat_path.exists():
            raise FileNotFoundError(f"Wormcat output path not found: {wormcat_out_path}")
            
        process_lst = []
        
        # Process each directory in the wormcat output path
        for dir_item in wormcat_path.iterdir():
            if dir_item.is_dir():
                dir_nm = dir_item.name
                
                # Process each category
                for cat_num in [1, 2, 3]:
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
                        print(f"Warning: File not found: {rgs_fisher_path}")

        if not process_lst:
            print("No valid category files found to process")
            return
            
        df_process = pd.DataFrame(process_lst, columns=['sheet', 'category', 'file', 'label'])
        process_category_files(df_process, annotation_file, out_xsl_file_nm)
        
    except Exception as e:
        print(f"Error creating summary spreadsheet: {str(e)}")


if __name__ == "__main__":
    # Example usage
    # create_summary_spreadsheet("/path/to/wormcat_output", "/path/to/annotation.csv", "/path/to/output.xlsx")
    pass