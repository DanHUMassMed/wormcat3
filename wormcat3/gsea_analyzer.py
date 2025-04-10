import gseapy as gp
import pandas as pd
import numpy  as np
import os
from typing import Union, Optional, Dict, List, Tuple


class GSEAAnalyzer:
    """
    A class to perform and manage Gene Set Enrichment Analysis (GSEA) using gseapy.
    """
    
    def __init__(self, output_dir: str = 'gsea_results'):
        """
        Initialize the GSEAAnalyzer.
        
        Parameters:
        -----------
        output_dir : str, optional
            Directory where GSEA results will be saved, default is 'gsea_results'
        """
        self.output_dir = output_dir
        self._ensure_output_directory()
        self.results = None
    
    def _ensure_output_directory(self) -> None:
        """Create output directory if it doesn't exist."""
        if not os.path.exists(self.output_dir):
            os.makedirs(self.output_dir)
    
    def run_preranked_gsea(self, 
                           ranked_genes: Union[str, pd.DataFrame], 
                           gene_sets: Union[str, Dict], 
                           min_size: int = 15, 
                           max_size: int = 500,
                           permutation_num: int = 1000, 
                           weighted_score_type: int = 1,
                           seed: int = 123, 
                           processes: int = 4,
                           verbose: bool = True) -> pd.DataFrame:
        """
        Perform pre-ranked GSEA analysis and return results as a DataFrame.
        
        Parameters:
        -----------
        ranked_genes : str or pd.DataFrame
            Ranked gene list. Can be a file path or a pandas DataFrame with 'Gene' and 'Rank' columns.
        gene_sets : str or dict
            Gene sets to analyze. Can be a GMT file path or a dictionary.
        min_size : int, optional
            Minimum size of gene sets to analyze (default: 15).
        max_size : int, optional
            Maximum size of gene sets to analyze (default: 500).
        permutation_num : int, optional
            Number of permutations (default: 1000).
        weighted_score_type : int, optional
            Weight type for the score (0 or 1, default: 1).
        seed : int, optional
            Random seed for reproducibility (default: 123).
        processes : int, optional
            Number of processes to use (default: 4).
        verbose : bool, optional
            Whether to display detailed output (default: True).
        
        Returns:
        --------
        pd.DataFrame
            DataFrame containing the GSEA results sorted by FDR.
        
        Raises:
        -------
        FileNotFoundError
            If the gene_sets file doesn't exist.
        ValueError
            If ranked_genes DataFrame doesn't have required columns.
        RuntimeError
            If GSEA analysis fails.
        """
        # Validate inputs
        if isinstance(gene_sets, str) and not os.path.exists(gene_sets):
            raise FileNotFoundError(f"Gene sets file not found: {gene_sets}")
        
        if isinstance(ranked_genes, pd.DataFrame):
            required_columns = {'Gene', 'Rank'}
            if not required_columns.issubset(ranked_genes.columns):
                raise ValueError(f"ranked_genes DataFrame must contain columns: {required_columns}")
        
        try:
            # Run pre-ranked GSEA
            prerank_results = gp.prerank(
                rnk=ranked_genes,
                gene_sets=gene_sets,
                outdir=self.output_dir,
                min_size=min_size,
                max_size=max_size,
                permutation_num=permutation_num,
                weighted_score_type=weighted_score_type,
                seed=seed,
                processes=processes,
                verbose=verbose
            )
            
            # Store the full results object
            self.results = prerank_results
            
            # Extract relevant results into a DataFrame
            results_list = []
            for term in list(prerank_results.results):
                term_results = prerank_results.results[term]
                results_list.append([
                    term,
                    term_results['fdr'],
                    term_results['es'],
                    term_results['nes'],
                    term_results['pval'],
                    term_results['tag %']
                ])
            
            # Create and sort the results DataFrame
            results_df = pd.DataFrame(
                results_list, 
                columns=['Term', 'FDR', 'ES', 'NES', 'P-value', 'Tag %']
            ).sort_values('FDR').reset_index(drop=True)
            
            return results_df
            
        except Exception as e:
            raise RuntimeError(f"GSEA analysis failed: {str(e)}")
    
    def get_enriched_terms(self, fdr_threshold: float = 0.25) -> pd.DataFrame:
        """
        Extract significantly enriched terms based on FDR threshold.
        
        Parameters:
        -----------
        fdr_threshold : float, optional
            FDR threshold for significance (default: 0.25).
        
        Returns:
        --------
        pd.DataFrame
            DataFrame containing only significant terms.
        
        Raises:
        -------
        ValueError
            If no analysis has been run yet.
        """
        if self.results is None:
            raise ValueError("No GSEA analysis has been run yet. Call run_preranked_gsea first.")
        
        results_df = self.run_preranked_gsea(None, None)  # This will reuse stored results
        return results_df[results_df['FDR'] <= fdr_threshold]
    
    def get_leading_edge_genes(self, term: str) -> List[str]:
        """
        Extract leading edge genes for a specific term.
        
        Parameters:
        -----------
        term : str
            The pathway or gene set term.
        
        Returns:
        --------
        List[str]
            List of leading edge genes.
        
        Raises:
        -------
        ValueError
            If no analysis has been run or term doesn't exist.
        """
        if self.results is None:
            raise ValueError("No GSEA analysis has been run yet. Call run_preranked_gsea first.")
        
        if term not in self.results.results:
            raise ValueError(f"Term '{term}' not found in GSEA results.")
        
        return self.results.results[term]['lead_genes']

    @staticmethod
    def create_ranked_list(deseq2_output_df):
        """
        Generates a ranked gene list from DESeq2 differential expression results.
        
        The ranking score is calculated as: sign(log2FoldChange) * -log10(pvalue)
        This produces a score that considers both the direction and magnitude of change,
        as well as the statistical significance.
        
        Parameters:
        -----------
        deseq2_output_df : pd.DataFrame
            DataFrame containing DESeq2 results with required columns:
            - 'ID': Gene identifiers
            - 'log2FoldChange': Log2 fold change values
            - 'pvalue': P-values indicating statistical significance
        
        Returns:
        --------
        pd.DataFrame
            DataFrame with 'Gene' and 'Rank' columns, sorted by 'Rank' 
            in descending order (most upregulated and significant genes at the top).
        
        Raises:
        -------
        ValueError
            If the input DataFrame doesn't contain the required columns.
        AssertionError
            If input data doesn't meet expected format or contains invalid values.
        """
        # Input validation
        assert isinstance(deseq2_output_df, pd.DataFrame), "Input must be a pandas DataFrame"
        assert not deseq2_output_df.empty, "Input DataFrame cannot be empty"
        
        # Create a copy to avoid modifying the original DataFrame
        deseq2_copy = deseq2_output_df.copy()
        
        # Ensure required columns are present
        required_columns = {'ID', 'log2FoldChange', 'pvalue'}
        missing_columns = required_columns - set(deseq2_copy.columns)
        if missing_columns:
            raise ValueError(f"Input DataFrame is missing required columns: {missing_columns}")
        
        # Validate identifier column
        assert not deseq2_copy['ID'].isna().any(), "Column 'ID' contains NaN values, which are not allowed for identifiers"
        assert deseq2_copy['ID'].duplicated().sum() == 0, "Column 'ID' contains duplicate values, which are not allowed"
        
        # Validate p-values are in the correct range
        valid_pvalues = deseq2_copy['pvalue'].dropna()
        if not valid_pvalues.empty:
            assert (valid_pvalues >= 0).all() and (valid_pvalues <= 1).all(), "P-values must be between 0 and 1"
        
        # Handle zero or NaN p-values to avoid issues with log transformation
        min_nonzero_pvalue = deseq2_copy.loc[deseq2_copy['pvalue'] > 0, 'pvalue'].min()
        if pd.isna(min_nonzero_pvalue):
            min_nonzero_pvalue = 1e-300  # Fallback if no non-zero p-values exist
        
        # Replace zeros and NaNs in p-values
        deseq2_copy['pvalue'] = deseq2_copy['pvalue'].replace(0, min_nonzero_pvalue)
        deseq2_copy['pvalue'] = deseq2_copy['pvalue'].fillna(1.0)  # Missing p-values get a non-significant value
        
        # Handle NaNs in log2FoldChange
        deseq2_copy['log2FoldChange'] = deseq2_copy['log2FoldChange'].fillna(0)
        
        # Rename 'ID' to 'Gene' for output consistency
        deseq2_copy = deseq2_copy.rename(columns={'ID': 'Gene'})
        
        # Calculate the ranking score
        deseq2_copy['Rank'] = np.sign(deseq2_copy['log2FoldChange']) * -np.log10(deseq2_copy['pvalue'])
        
        # Verify no NaN or infinite values in rank
        assert not deseq2_copy['Rank'].isna().any(), "Rank calculation produced NaN values"
        assert not np.isinf(deseq2_copy['Rank']).any(), "Rank calculation produced infinite values"
        
        # Sort the DataFrame by ranking score in descending order
        ranked_list = deseq2_copy[['Gene', 'Rank']].sort_values(by='Rank', ascending=False)
        
        # Final validation of output
        assert ranked_list.shape[0] == deseq2_copy.shape[0], "Output has different number of rows than input"
        
        return ranked_list
    
def convert_dataframe_to_gmt(input_df, output_file='output.gmt', id_col='Function.ID', 
                            desc_col=None, gene_col='Wormbase.ID'):
    """
    Convert a pandas DataFrame to GMT file format.
    
    Parameters:
    -----------
    input_df : pandas.DataFrame
        Input DataFrame containing gene set information
    output_file : str
        Path to output GMT file
    id_col : str
        Column name for the gene set ID (also used as description if desc_col is None)
    desc_col : str or None
        Column name for the gene set description (if None, id_col will be used instead)
    gene_col : str
        Column name for the gene identifiers
    
    Returns:
    --------
    str
        Path to the created GMT file
    """
    # Assert input is a DataFrame
    assert isinstance(input_df, pd.DataFrame), "Input must be a pandas DataFrame"
    assert not input_df.empty, "Input DataFrame cannot be empty"
    
    # Ensure the output directory exists
    output_dir = os.path.dirname(output_file)
    if output_dir and not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    # Validate that required columns exist
    required_cols = [id_col, gene_col]
    if desc_col is not None:
        required_cols.append(desc_col)
    
    missing_cols = [col for col in required_cols if col not in input_df.columns]
    if missing_cols:
        raise ValueError(f"Missing required columns: {', '.join(missing_cols)}")
    
    # Assert ID column has no NaN values
    assert not input_df[id_col].isna().any(), f"Column '{id_col}' contains NaN values"
    
    # Group by required columns
    if desc_col is not None:
        # Use description column if provided
        grouped = input_df.groupby([id_col, desc_col])[gene_col].apply(list).reset_index()
        # Count unique gene sets
        num_gene_sets = len(input_df[[id_col, desc_col]].drop_duplicates())
    else:
        # Use ID column as description if desc_col is None
        grouped = input_df.groupby([id_col])[gene_col].apply(list).reset_index()
        # Add ID column as description column
        grouped[id_col + '_desc'] = grouped[id_col]
        desc_col = id_col + '_desc'
        # Count unique gene sets
        num_gene_sets = len(input_df[id_col].drop_duplicates())
    
    # Assert there's at least one gene set
    assert len(grouped) > 0, "No gene sets found after grouping"
    
    # Count total genes before filtering
    total_genes = sum(len(genes) for genes in grouped[gene_col])
    
    # Write to GMT file
    gene_sets_written = 0
    genes_written = 0
    
    with open(output_file, 'w') as file:
        for _, row in grouped.iterrows():
            # Handle potential NaN values in the description
            description = row[desc_col] if pd.notna(row[desc_col]) else row[id_col]
            
            # Filter out any None or NaN values from gene list
            gene_list = [str(gene) for gene in row[gene_col] if pd.notna(gene)]
            
            # Only write if there are genes in the set
            if gene_list:
                line = f"{row[id_col]}\t{description}\t" + '\t'.join(gene_list) + '\n'
                file.write(line)
                gene_sets_written += 1
                genes_written += len(gene_list)
    
    # Final assertions to ensure data was written
    assert gene_sets_written > 0, "No gene sets were written to the output file"
    assert os.path.exists(output_file), f"Output file {output_file} was not created"
    assert os.path.getsize(output_file) > 0, f"Output file {output_file} is empty"
    
    print(f"Successfully created GMT file: {output_file}")
    print(f"Processed {num_gene_sets} gene sets, wrote {gene_sets_written} sets with at least one gene")
    print(f"Total genes in input: {total_genes}, total genes written: {genes_written}")
    
    return output_file