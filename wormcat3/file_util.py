import os
import uuid
from pathlib import Path
import pandas as pd
import re
import zipfile
from wormcat3.wormcat_error import WormcatError, ErrorCode

def validate_directory_path(directory_path, validation_indicator=False, not_empty_check = True):
    """
    Ensure the directory exists and is empty or can be created.
    Returns the validated path.
    """
    path = Path(directory_path)

    if path.exists():
        if not path.is_dir():
            return (False, str(path)) if validation_indicator else str(path)
        if not_empty_check and any(path.iterdir()):
            return (False, str(path)) if validation_indicator else str(path)
    else:
        try:
            path.mkdir(parents=True, exist_ok=True)
        except OSError as e:
            print(f"[Error] Failed to create directory '{path}': {e}")
            return (False, str(path)) if validation_indicator else str(path)

    return (True, str(path)) if validation_indicator else str(path)


def find_file_path(file_name, additional_search_paths=[]):
    """Search for a file in the given list of directories."""
    
    file_path = Path(__file__).resolve()
    extdata_path = file_path.parent / "extdata"
    
    search_paths = [
        os.getcwd(),  # Current working directory
        extdata_path  # Defaults location
    ]
    
    # Check if the environment variable is set, and add it to the front of the path if it exists
    wormcat_data_path = os.environ.get("WORMCAT_DATA_PATH")
    if wormcat_data_path:
        search_paths.insert(0, Path(wormcat_data_path))  # Add it as a Path object
    
    search_paths += additional_search_paths
    for directory in search_paths:
        file_path = Path(directory) / file_name
        if file_path.exists():
            return str(file_path)  # Return the first found file path
    return None  # File not found in any directory

def get_name_from_file_path(file_path):
    # Get the base name of the file (without directory and extension)
    if isinstance(file_path, (str, os.PathLike)) and file_path:
        return os.path.splitext(os.path.basename(file_path))[0]
    else:
        return "[NO_NAME_FOUND]"
        
def read_deseq2_file(file_path):
    """Read deseq2 file"""
    
    if not Path(file_path).exists():
        raise WormcatError(f"Attempting to read a non-existent file: {file_path}", ErrorCode.FILE_NOT_FOUND.to_dict())
    
    # Read the CSV file into a DataFrame
    deseq2_df = pd.read_csv(file_path)

    # Ensure required columns are present
    required_columns = {'ID', 'log2FoldChange', 'pvalue'}
    missing_columns = required_columns - set(deseq2_df.columns)
    if missing_columns:
        raise WormcatError(f"{get_name_from_file_path(file_path)} file is missing required columns: {missing_columns}", ErrorCode.MISSING_FIELD.to_dict())
    
    return deseq2_df


def read_gene_set_file(file_path):
    """Read the first column of a CSV file as a list."""
    
    path = Path(file_path)
    if not path.is_file():
        raise WormcatError(
            message=f"Attempting to read a non-existent or invalid file: {file_path}",
            code=ErrorCode.FILE_NOT_FOUND.to_dict(),
            origin="read_gene_set_file"
        )

    try:
        df = pd.read_csv(file_path)
        return df.iloc[:, 0].tolist()
    except Exception as e:
        raise WormcatError(
            message=f"Unexpected error while reading file: {file_path}",
            code=ErrorCode.INTERNAL_ERROR.to_dict(),
            origin="read_gene_set_file",
            detail={"error": str(e)}
        )

def is_file_path(input_string: str) -> bool:
    """
    Check if the given string is a valid file path or just a file name.
    """
    input_path = Path(input_string)
    
    # Check if it contains directories (if it contains any path separator, it's considered a file path)
    if input_path.is_absolute() or os.path.sep in input_string:
        return True
    return False
    
def sanitize(text: str) -> str:
    # Replace spaces with underscores
    text = text.replace(" ", "_")
    # Remove invalid directory characters: / \ : * ? " < > | 
    text = re.sub(r'[\/:*?"<>|]', '', text)
    return text

def generate_5_digit_hash(prefix: str = "", suffix: str = "") -> str:
    prefix = sanitize(prefix)
    suffix = sanitize(suffix)
    core = str(uuid.uuid4().int % 100000).zfill(5)
    return f"{prefix}{core}{suffix}"

def extract_run_number(data_file_nm):
    # Define the regex pattern to capture the run number (5 digits after "run_")
    pattern = r"run_(\d{5})\.csv$"
    
    # Use re.search to find the match
    match = re.search(pattern, data_file_nm)
    
    # If a match is found, extract the run number
    if match:
        run_number = match.group(1)  # Extract the 5-digit run number
        return run_number
    else:
        raise WormcatError(f"Invalid file name: {data_file_nm}. It must end with 'run_?????.csv' where '?????' are any 5 digits.", 
                            ErrorCode.INVALID_NAME.to_dict())


def zip_dir(path, output_zip_path=None):
    """
    Recursively zips the contents of a directory.
    """
    
    path = os.path.abspath(path)

    if output_zip_path is None:
        output_zip_path = f"{path}.zip"
    else:
        output_zip_path = os.path.abspath(output_zip_path)
        output_dir = os.path.dirname(output_zip_path)
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)

    with zipfile.ZipFile(output_zip_path, 'w', zipfile.ZIP_DEFLATED) as zipf:
        for root, _, files in os.walk(path):
            for file in files:
                abs_filepath = os.path.join(root, file)
                rel_filepath = os.path.relpath(abs_filepath, start=path)
                zipf.write(abs_filepath, arcname=rel_filepath)

    return output_zip_path