import os
import uuid
from pathlib import Path
import re

def validate_directory_path(directory_path, not_empty_check = True):
    """
    Ensure the directory exists and is empty or can be created.
    Returns the validated path.
    """
    path = Path(directory_path)
    if path.exists():
        if path.is_dir():
            if not_empty_check and any(path.iterdir()):  # Directory is not empty
                raise ValueError(f"The directory '{path}' exists and is not empty.")
        else:
            raise ValueError(f"The path '{path}' exists and is not a directory.")
    else:
        path.mkdir(parents=True)
    return str(path)

def find_file_path(file_name, additional_search_paths=[]):
    """Search for a file in the given list of directories."""
    
    file_path = Path(__file__).resolve()
    extdata_path = file_path.parent / "extdata"
    
    search_paths = [
        os.getcwd(),  # Current working directory
        extdata_path  # Defaults location
    ]
    
    # Check if the environment variable is set, and add it first if it exists
    wormcat_data_path = os.environ.get("WORMCAT_DATA_PATH")
    if wormcat_data_path:
        search_paths.insert(0, Path(wormcat_data_path))  # Add it as a Path object
    
    search_paths += additional_search_paths
    for directory in search_paths:
        file_path = Path(directory) / file_name
        if file_path.exists():
            return str(file_path)  # Return the first found file path
    return None  # File not found in any directory

def read_gene_set_file(file_path):
    """Read the first column of a CSV file as a list."""
    import pandas as pd
    
    if not Path(file_path).exists():
        raise FileNotFoundError(f"The file {file_path} does not exist.")
    
    # Read the CSV file into a DataFrame
    df = pd.read_csv(file_path)
    
    # Extract the first column as a list
    return df.iloc[:, 0].tolist()

def is_file_path(input_string: str) -> bool:
    """
    Check if the given string is a valid file path or just a file name.
    """
    input_path = Path(input_string)
    
    # Check if it contains directories (if it contains any path separator, it's considered a file path)
    if input_path.is_absolute() or os.path.sep in input_string:
        return True
    return False
    
def generate_5_digit_hash(prefix: str = "", suffix: str = "") -> str:
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
        raise ValueError(f"Invalid file name: {data_file_nm}. It must end with 'run_00000.csv' where '00000' are any 5 digits.")


