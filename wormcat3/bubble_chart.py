import pandas as pd
import numpy as np
from plotnine import *
from pathlib import Path
import wormcat3.constants as cs
import warnings

# Suppress UserWarnings and DeprecationWarnings from plotnine
warnings.filterwarnings('ignore', category=UserWarning, module='plotnine')
warnings.filterwarnings('ignore', category=DeprecationWarning, module='plotnine')
    
def preprocess_bubble_data(data_file_path, add_calibration=False):
    """
    Preprocess CSV data for bubble plot visualization.
    
    Args:
        csv_file_name (str): Path to the CSV file containing the data.
        
    Returns:
        pd.DataFrame: Processed DataFrame ready for plotting.
    """
    try:
        # Read the CSV file
        bubbles_data = pd.read_csv(data_file_path)
        column_names = bubbles_data.columns.tolist()
        correction_method = 'Bonferroni' if 'Bonferroni' in column_names else 'FDR'
        
        # Remove rows with missing values
        bubbles_data = bubbles_data.dropna()
        
        # # Add calibration points
        if add_calibration:
            calibration_points = pd.DataFrame({
                "Category": ["calibration high", "calibration low"],
                "RGS": [250, 1],
                "AC": [0, 0],
                "PValue": [1.00E-50, 1],
                correction_method: [1.00E-50, 1]
            })
            bubbles_data = pd.concat([bubbles_data, calibration_points], ignore_index=True)
        
        # Compute normalized value as a placeholder
        scale_for_visualization=0.001
        bubbles_data['bubbles_z'] = round(scale_for_visualization * (bubbles_data['PValue'] - bubbles_data['PValue'].mean()) / bubbles_data['PValue'].std(), 2)
        
        # Create color coding based on Bonferroni-corrected p-values
        conditions = [     
            (bubbles_data[correction_method] < 1e-40),
            (bubbles_data[correction_method] < 1e-20),
            (bubbles_data[correction_method] < 1e-10),
            (bubbles_data[correction_method] < 1e-5),
            (bubbles_data[correction_method] < 1e-2),
            (bubbles_data[correction_method] < 5e-2),
            (bubbles_data[correction_method] >= 5e-2)
        ]
        
        values = ["Col1", "Col2", "Col3", "Col4", "Col5", "Col6", "NS"]
        bubbles_data['p_value_type'] = np.select(conditions, values, default='')

        
        # Create size coding based on RGS values
        size_conditions = [
            (bubbles_data['RGS'] > 150),
            (bubbles_data['RGS'] > 100),
            (bubbles_data['RGS'] > 75),
            (bubbles_data['RGS'] > 50),
            (bubbles_data['RGS'] > 25),
            (bubbles_data['RGS'] > 10),
            (bubbles_data['RGS'] > 5),
            (bubbles_data['RGS'] > 2),
            (bubbles_data['RGS'] <= 2)
        ]
        
        size_values = ["Size1", "Size2", "Size3", "Size4", "Size5", "Size6", "Size7", "Size8", "Size9"]
        bubbles_data['RGS_size'] = np.select(size_conditions, size_values, default='')
        
        bubbles_data = order_categories_by_column(
            bubbles_data, 
            category_column = 'Category', 
            order_by_column = correction_method, 
            ascending=False
        )
        
        return bubbles_data
        
    except FileNotFoundError:
        raise FileNotFoundError(f"CSV file not found: {data_file_path}")
    except pd.errors.EmptyDataError:
        raise ValueError(f"The CSV file is empty: {data_file_path}")
    except pd.errors.ParserError:
        raise ValueError(f"Error parsing CSV file: {data_file_path}")
    except KeyError as e:
        raise KeyError(f"Required column missing in CSV file: {e}")


def order_categories_by_column(df, category_column='Category', order_by_column='PValue', ascending=True):
    """
    Convert a column to a categorical type ordered by another column.
    
    Args:
        df (pd.DataFrame): DataFrame with category and ordering columns.
        category_column (str): Name of the column to be converted to categorical.
        order_by_column (str): Name of the column to order by.
        ascending (bool): Sort order - True for ascending, False for descending.
        
    Returns:
        pd.DataFrame: DataFrame with category column converted to categorical type.
    """
    if df is None or df.empty:
        raise ValueError("DataFrame is empty or None")
        
    if category_column not in df.columns:
        raise KeyError(f"Category column '{category_column}' not found in DataFrame")
        
    if order_by_column not in df.columns:
        raise KeyError(f"Order column '{order_by_column}' not found in DataFrame")
    
    # Sort the dataframe by the specified column
    sorted_df = df.sort_values(by=order_by_column, ascending=ascending)
    
    # Get the categories in the order determined by the sorting column
    category_order = sorted_df[category_column].unique()
    
    # Convert category column to categorical type with the specified order
    df[category_column] = pd.Categorical(
        df[category_column], 
        categories=category_order,
        ordered=True
    )
    
    return df
    
def generate_bubble_plot(df, svg_file_path, plot_title=cs.DEFAULT_TITLE, width=cs.DEFAULT_WIDTH, height=cs.DEFAULT_HEIGHT):
    """
    Generate and save a bubble plot as an SVG using plotnine.

    """
    if df is None or df.empty:
        raise ValueError("DataFrame is empty or None")
        
    required_columns = ['Category', 'bubbles_z', 'RGS_size', 'p_value_type']
    missing_columns = [col for col in required_columns if col not in df.columns]
    if missing_columns:
        raise KeyError(f"Required columns missing in DataFrame: {missing_columns}")

    # Define color and size mapping dictionaries
    color_mapping = {
        "Col1": "darkred",
        "Col2": "orangered",
        "Col3": "tomato",
        "Col4": "orange",
        "Col5": "gold",
        "Col6": "yellow",
        "NS": "dimgray"
    }
    
    size_mapping = {
        "Size1": 10,
        "Size2": 9,
        "Size3": 8,
        "Size4": 7,
        "Size5": 6,
        "Size6": 5,
        "Size7": 2.5,
        "Size8": 1,
        "Size9": 0.1
    }

    # Manual color and size scales
    color_scale = scale_color_manual(
        name="P value",
        labels=["10-40", "10-20", "10-10", "10-5", "0.001", "0.05", "NS"],
        values=color_mapping,
        limits=list(color_mapping.keys())
    )

    size_scale = scale_size_manual(
        name="Count",
        labels=["150", "100", "75", "50", "25", "10", "5", "2", "Zero"],
        values=size_mapping,
        limits=list(size_mapping.keys())
    )

    # Build the plot
    p = (
        ggplot(df, aes(x='Category', y='bubbles_z', size='RGS_size', color='p_value_type')) +
        ggtitle(title=plot_title) +
        geom_point(stat='identity') +
        color_scale +
        size_scale +
        labs(y=' ') +
        ylim(-1.0, 1.0) +
        coord_flip() +
        theme(
            panel_grid=element_blank(),
            panel_background=element_blank(),
            legend_key=element_rect(fill="white", color="white"),
            text=element_text(family="Arial", size=14),
            axis_text_x=element_blank(),   # Remove x-axis tick labels
            axis_ticks_x=element_blank(),  # Remove x-axis tick marks
            plot_title=element_text(margin={"b": 3, "units": "lines"})
        )
    )

    try:
        p.save(filename=svg_file_path, format = "svg", width = width, height = height)
        print(f"Saved plot to: {svg_file_path}")
    except Exception as e:
        raise IOError(f"Failed to save plot: {e}")


def create_bubble_chart(dir_path: str, data_file_nm: str, plot_title="RGS", add_calibration=False) -> None:
    try:
        data_file_path = Path(dir_path) / data_file_nm
        svg_file_nm = data_file_nm[:-3] + 'svg'
        svg_file_path = Path(dir_path) / svg_file_nm
        
        bubbles_data = preprocess_bubble_data(data_file_path, add_calibration = add_calibration)
        
        # Scale the height of the bubble chart based on the number of items
        bubbles_data_len = len(bubbles_data)
        height_fun = lambda x: min(0.1 * x + 6.3, 25)
        height = height_fun(bubbles_data_len)
        
        generate_bubble_plot(bubbles_data, svg_file_path, plot_title = plot_title,height=height, width=9)
    except Exception as e:
        print(f"Error: {e}")


if __name__ == "__main__":
    dir_path="/Users/dan/Code/Python/wormcat3/notebooks/wormcat_out/EC-ES-UP_73980"
    data_file_nm= "category_2_padj_bon_EC-ES-UP_73980.csv"

    create_bubble_chart(dir_path, data_file_nm, plot_title="Category 2 Up")
    
