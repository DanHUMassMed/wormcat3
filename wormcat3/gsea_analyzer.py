import gseapy as gp
import pandas as pd
import numpy  as np
import os
from pathlib import Path
from typing import Union, Dict, List
import wormcat3.constants as cs
from wormcat3 import file_util
from wormcat3.wormcat_error import WormcatError, ErrorCode

class GSEAAnalyzer:
    """
    A class to perform and manage Gene Set Enrichment Analysis (GSEA) using gseapy.
    """
    
    def __init__(self, output_dir: str = cs.DEFAULT_GSEA_RESULTS_DIR):
        """
        Initialize the GSEAAnalyzer.
        """
        self.required_columns = ['ID', 'log2FoldChange', 'pvalue']

        self.output_dir = output_dir
        self._ensure_output_directory()
        self.results = None
        #last_dir = Path(working_dir_path).name  # "output"
        
    
    def _ensure_output_directory(self) -> None:
        """Create output directory if it doesn't exist."""
        if not os.path.exists(self.output_dir):
            os.makedirs(self.output_dir)
    
    def run_preranked_gsea(self, 
                           ranked_genes: Union[str, pd.DataFrame], 
                           gene_sets: Union[str, Dict],
                           output_dir: str,
                           *,
                           min_size: int = 15, 
                           max_size: int = 500,
                           permutation_num: int = 1000, 
                           weight: float = 1.0,
                           seed: int = 123, 
                           threads: int = 4,
                           verbose: bool = False) -> pd.DataFrame:
        """
        Perform pre-ranked GSEA analysis and return results as a DataFrame.
        """
        # Validate inputs
        if isinstance(gene_sets, str) and not os.path.exists(gene_sets):
            raise WormcatError(f"Gene sets file not found: {gene_sets}", ErrorCode.FILE_NOT_FOUND.to_dict())
        
        if isinstance(ranked_genes, pd.DataFrame):
            required_columns = {'Gene', 'Rank'}
            if not required_columns.issubset(ranked_genes.columns):
                raise WormcatError(f"ranked_genes DataFrame must contain columns: {required_columns}", ErrorCode.MISSING_FIELD.to_dict())
        
        outdir = file_util.validate_directory_path(Path(self.output_dir)/output_dir)

        try:
            # Run pre-ranked GSEA
            prerank_results = gp.prerank(
                rnk = ranked_genes,
                gene_sets = gene_sets,
                outdir = outdir,
                min_size = min_size,
                max_size = max_size,
                permutation_num = permutation_num,
                weight = weight,
                seed = seed,
                threads = threads,
                verbose = verbose
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
            raise WormcatError(f"GSEA analysis failed: {str(e)}", ErrorCode.INTERNAL_ERROR.to_dict(), origin="GSEAAnalyzer.run_preranked_gsea")
    
    def get_enriched_terms(self, fdr_threshold: float = 0.25) -> pd.DataFrame:
        """
        Extract significantly enriched terms based on FDR threshold.
        """
        if self.results is None:
            raise WormcatError("No GSEA analysis has been run yet. Call run_preranked_gsea first.", ErrorCode.CONSTRAINT_VIOLATION.to_dict())
        
        results_df = self.run_preranked_gsea(None, None)  # This will reuse stored results
        return results_df[results_df['FDR'] <= fdr_threshold]
    
    def get_leading_edge_genes(self, term: str) -> List[str]:
        """
        Extract leading edge genes for a specific term.
        """
        if self.results is None:
            raise WormcatError("No GSEA analysis has been run yet. Call run_preranked_gsea first.", ErrorCode.CONSTRAINT_VIOLATION.to_dict())
        
        if term not in self.results.results:
            raise WormcatError(f"Term '{term}' not found in GSEA results.", ErrorCode.INVALID_VALUE.to_dict())
        
        return self.results.results[term]['lead_genes']


    def create_ranked_list(self, deseq2_output_df):
        """
        Generates a ranked gene list from DESeq2 differential expression results.
        
        The ranking score is calculated as: sign(log2FoldChange) * -log10(pvalue)
        This produces a score that considers both the direction and magnitude of change,
        as well as the statistical significance.

        """
        # Input validation
        assert isinstance(deseq2_output_df, pd.DataFrame), "Input must be a pandas DataFrame"
        assert not deseq2_output_df.empty, "Input DataFrame cannot be empty"
        
        # Create a copy to avoid modifying the original DataFrame
        deseq2_copy = deseq2_output_df.copy()
        
        # Ensure required columns are present
        missing_columns = set(self.required_columns) - set(deseq2_copy.columns)
        if missing_columns:
            raise WormcatError(f"Input DataFrame is missing required columns: {missing_columns}", ErrorCode.MISSING_FIELD.to_dict())
        
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
        
        ranked_list = GSEAAnalyzer._make_ranks_unique(ranked_list)
        
        # Check for duplicates in the 'Gene' column
        duplicate_genes = ranked_list[ranked_list['Gene'].duplicated(keep=False)]

        # Print any duplicates found
        if not duplicate_genes.empty:
            print(f"Found {len(duplicate_genes['Gene'].unique())} genes with duplicates:")
            duplicate_percent = (len(duplicate_genes) / len(ranked_list)) * 100
            print(f"Duplicated genes represent {duplicate_percent:.2f}% of the dataset")

        # Remove duplicates, keeping the entry with the highest rank for each gene
        # Since the list is already sorted by Rank (descending), we keep the first occurrence
        ranked_list_no_duplicates = ranked_list.drop_duplicates(subset='Gene', keep='first')

        # Verify duplicates were removed
        assert ranked_list_no_duplicates['Gene'].duplicated().sum() == 0, "Duplicates still exist!"

        # Final validation of output
        assert ranked_list.shape[0] == deseq2_copy.shape[0], "Output has different number of rows than input"
        
        return ranked_list

    @staticmethod
    def _make_ranks_unique(ranked_list):
        """
        Check for duplicates in the 'Rank' column and make them unique by adding small values
        that won't change the overall sorting order.
        """
        # Check for duplicates in the 'Rank' column
        duplicate_ranks = ranked_list['Rank'].duplicated(keep=False)
        num_duplicates = duplicate_ranks.sum()
        
        if num_duplicates > 0:
            
            # Store original order to ensure we don't change it
            ranked_list['original_order'] = range(len(ranked_list))
            
            # Group by Rank value to find duplicates
            for rank, group in ranked_list[duplicate_ranks].groupby('Rank'):
                
                # Calculate a small value that won't change the sorting order
                # Get the difference to the next smaller rank value
                idx = ranked_list['Rank'].unique().tolist().index(rank)
                if idx < len(ranked_list['Rank'].unique()) - 1:
                    next_smaller_rank = sorted(ranked_list['Rank'].unique(), reverse=True)[idx + 1]
                    epsilon = abs(rank - next_smaller_rank) / (len(group) * 2)
                else:
                    # If it's the smallest rank, use a small fraction of its absolute value
                    epsilon = abs(rank) / 1000000 if rank != 0 else 0.0000001
                    
                # Add incrementally smaller values to maintain the original order
                for i, (idx, row) in enumerate(group.iterrows()):
                    # Smaller adjustment for higher ranks (preserves descending order)
                    adjustment = epsilon * (i + 1) / (len(group) + 1)
                    ranked_list.at[idx, 'Rank'] = rank - adjustment
                    
            # Resort to ensure order is maintained
            ranked_list = ranked_list.sort_values(by='Rank', ascending=False)
            
            # Remove the temporary column
            ranked_list = ranked_list.drop('original_order', axis=1)
            
            # Verify no duplicates remain
            assert ranked_list['Rank'].duplicated().sum() == 0, "Failed to remove all duplicate rank values"
            print("Successfully made all Rank values unique while preserving order")
        else:
            print("No duplicate Rank values found")
            
        return ranked_list
    

    def clean_input_data(self, deseq2_df: pd.DataFrame):
        """
        Clean a DESeq2-like DataFrame by removing rows with NaNs in key columns
        and duplicate IDs (keeping the one with the lowest pvalue).

        Args:
            deseq2_df (pd.DataFrame): The input DataFrame with columns 'ID', 'log2FoldChange', 'pvalue'.

        Returns:
            Tuple[pd.DataFrame, pd.DataFrame]: (removed_rows_df, cleaned_deseq2_df)
        """
        # Ensure required columns are present
        missing_columns = set(self.required_columns) - set(deseq2_df.columns)
        if missing_columns:
            raise WormcatError(
                f"Input DataFrame is missing required columns: {missing_columns}",
                ErrorCode.MISSING_FIELD.to_dict()
            )

        # Step 1: Remove rows with NaNs in key columns
        na_mask = deseq2_df[self.required_columns].isna().any(axis=1)
        na_removed_df = deseq2_df[na_mask].copy()
        na_removed_df['reason'] = 'has_na'
        no_na_df = deseq2_df[~na_mask]

        # Step 2: Remove duplicates by keeping row with lowest pvalue
        sorted_df = no_na_df.sort_values('pvalue', ascending=True)
        dedup_df = sorted_df.drop_duplicates(subset='ID', keep='first')

        # Identify duplicates that were removed (not already marked as NaN)
        removed_duplicates_df = sorted_df[~sorted_df.index.isin(dedup_df.index)].copy()
        removed_duplicates_df['reason'] = 'duplicate'

        # Combine all removed rows with reason
        removed_rows_df = pd.concat([na_removed_df, removed_duplicates_df], ignore_index=True)

        return removed_rows_df, dedup_df