from pathlib import Path
import pandas as pd
from scipy.stats import fisher_exact
from statsmodels.stats.multitest import multipletests
import wormcat3.constants as cs
from wormcat3.constants import PAdjustMethod
from wormcat3.wormcat_error import WormcatError, ErrorCode
from typing import Callable
    
class EnrichmentAnalyzer:
    """Performs statistical enrichment analysis."""
    
    def __init__(self, annotations_df, output_dir, run_number=""):
        """Initialize with annotation dataframe and output directory."""
        self.annotations_df = annotations_df
        self.output_dir = output_dir
        self.run_number = run_number
        self.contingency_func = self._default_create_contingency
        self.categories = [1, 2, 3]  # Wormcat Categories
    
    def set_contingency_func(self, contingency_func: Callable):
        self.contingency_func = contingency_func
        
    def calculate_category_enrichment_scores(self, gene_set_and_categories_df, p_adjust_method=PAdjustMethod.BONFERRONI, p_adjust_threshold=0.01):
        """Run enrichment test for all categories."""
        
        if not isinstance(p_adjust_method, PAdjustMethod):
            raise WormcatError(f"Invalid p_adjust_method: {p_adjust_method}. Must be a valid PAdjustMethod.", ErrorCode.INVALID_VALUE.to_dict())

        try:
            p_adjust_threshold = float(p_adjust_threshold)
        except (TypeError, ValueError):
            raise WormcatError(
                f"Invalid p_adjust_threshold: {p_adjust_threshold}. Must be convertible to a float.",
                ErrorCode.INVALID_VALUE.to_dict()
            )

        if not (0 < p_adjust_threshold <= 1):
            raise WormcatError(
                f"Invalid p_adjust_threshold: {p_adjust_threshold}. Must be > 0 and ≤ 1.",
                ErrorCode.INVALID_VALUE.to_dict()
            )
                
        enrichment_scores_list = []
        
        for category in self.categories:
            fisher_cat_df = self._run_fisher_test(gene_set_and_categories_df, category)
            
            fisher_cat_adjusted_df = self._adjust_pvalues(
                fisher_cat_df, 
                category, 
                method=p_adjust_method.value, 
                threshold=p_adjust_threshold
            )
            enrichment_scores_list.append(fisher_cat_adjusted_df)
            
        return enrichment_scores_list
    
    def _run_fisher_test(self, gene_set_and_categories_df, category):
        """Run Fisher's exact test for a specific category."""
        category_column = f"Category.{category}"
        total_annotations_count = len(self.annotations_df)
        total_gene_set_count = len(gene_set_and_categories_df)
        
        gene_set_category_df = gene_set_and_categories_df[category_column].value_counts().reset_index()
        annotated_category_df = self.annotations_df[category_column].value_counts().reset_index()
        
        merged_categories_df = pd.merge(
            gene_set_category_df, 
            annotated_category_df, 
            how="left", 
            left_on=category_column, 
            right_on=category_column
        )
        merged_categories_df = merged_categories_df.rename(
            columns={category_column: "Category", "count_x": "RGS", "count_y": "AC"}
        )
        
        # Initialize results dataframe
        fisher_cat_df = pd.DataFrame(columns=["Category", "RGS", "AC", "PValue"])
        
        # Iterate over each row of the merged data
        for _, row in merged_categories_df.iterrows():
            rgs_value = row["RGS"]
            ac_value = row["AC"]
            
            if pd.isna(rgs_value) or pd.isna(ac_value):
                pvalue = None
            else:
                contingency_table = self.contingency_func(
                    rgs_value, 
                    total_gene_set_count, 
                    ac_value, 
                    total_annotations_count
                )
                _, pvalue = fisher_exact(contingency_table, alternative="greater")
            
            df_row = {"Category": row["Category"], "RGS": rgs_value, "AC": ac_value, "PValue": pvalue}

            fisher_cat_df.loc[len(fisher_cat_df)] = df_row
        
        # Sort and save
        fisher_cat_df = fisher_cat_df.sort_values(by=["PValue", "Category"],key=lambda col: col.str.lower() if col.name == "Category" else col)
        fisher_cat_file_path = Path(self.output_dir) / f"category_{category}_fisher_{self.run_number}.csv"
        fisher_cat_df.to_csv(fisher_cat_file_path, index=False)
        
        return fisher_cat_df
    
    def _adjust_pvalues(self, fisher_cat_df, category, *, method='bonferroni', threshold=0.01):
        """Adjust p-values using the specified method."""
        
        if method not in {'bonferroni', 'fdr_bh'}:
            raise WormcatError("Invalid method. Choose either 'bonferroni' or 'fdr_bh'.", ErrorCode.INVALID_VALUE.to_dict())
        
        padj_col = 'Bonferroni' if method == 'bonferroni' else 'FDR'
        
        fisher_cat_adjusted_df = fisher_cat_df.copy(deep=True)
        
        # Drop rows with any missing values
        fisher_cat_adjusted_df.dropna(inplace=True)
        
        if fisher_cat_adjusted_df.empty:
            return fisher_cat_adjusted_df
        
        # Sort by p-value
        fisher_cat_adjusted_df.sort_values(by=["PValue", "Category"], key=lambda col: col.str.lower() if col.name == "Category" else col, inplace=True)
        
        # Apply correction
        _, corrected_pvals, _, _ = multipletests(fisher_cat_adjusted_df['PValue'], method=method)
        fisher_cat_adjusted_df[padj_col] = corrected_pvals
        
        # Filter by threshold
        fisher_cat_adjusted_df = fisher_cat_adjusted_df[fisher_cat_adjusted_df[padj_col] < threshold]
        
        # Save results
        output_file_path = Path(self.output_dir) / f"category_{category}_padj_{self.run_number}.csv"
        fisher_cat_adjusted_df.to_csv(output_file_path, index=False)
        
        return {output_file_path: fisher_cat_adjusted_df}
    
    
    @staticmethod
    def _default_create_contingency(genes_in_both, gene_set_size, category_size, background_size):
        """
        Create a proper 2x2 contingency table for Fisher's exact test.
        
        Parameters:
        - genes_in_both: Number of genes in both the gene set and category
        - gene_set_size: Total number of genes in the gene set
        - category_size: Total number of genes in the category
        - background_size: Total number of genes in the background
        
        Returns:
        - A 2x2 contingency table as a list of lists
        """
        
        a = genes_in_both  # In both gene set and category
        b = gene_set_size - a  # In gene set but not in category
        c = category_size - a  # In category but not in gene set
        d = background_size - a - b - c  # In neither
                
        errors = []

        if not all(isinstance(x, int) and x >= 0 for x in [genes_in_both, gene_set_size, category_size, background_size]):
            errors.append("All input values must be non-negative integers.")

        if genes_in_both > category_size:
            errors.append("genes_in_both cannot exceed category_size.")

        if genes_in_both > gene_set_size:
            errors.append("genes_in_both cannot exceed gene_set_size.")

        if category_size > background_size:
            errors.append("category_size cannot exceed background_size.")

        if gene_set_size > background_size:
            errors.append("gene_set_size cannot exceed background_size.")

        if d < 0:
            errors.append("Calculated 'In neither' is negative. Check your inputs.")

        total_check = a + b + c + d # Final sanity check
        if total_check != background_size:
            errors.append(f"Sanity check failed: table total ({total_check}) != background_size ({background_size})")

        if errors:
            raise WormcatError(" ".join(errors), ErrorCode.INVALID_VALUE.to_dict())
        
        return [[a, b], [c, d]]
