import os
import gzip
import pandas as pd
import matplotlib.pyplot as plt
import scipy.stats as stats
import seaborn as sns
import matplotlib.pyplot as plt

# file path
GZ_FPATH = os.path.join('data', 'GSE134269_xeno_tpm.txt.gz')

# Step 1: Load data from .gz file
def load_data():
    # Print the file path for debugging
    print(f"Looking for file: {os.path.abspath(GZ_FPATH)}")
    
    # Check if the file exists
    if not os.path.exists(GZ_FPATH):
        print(f"ERROR: File not found: {GZ_FPATH}")
        raise FileNotFoundError(f"The file '{GZ_FPATH}' does not exist.")
    
    # Load the file if it exists
    with gzip.open(GZ_FPATH, 'rt') as gz:
        df = pd.read_csv(gz, delim_whitespace=True, index_col=0)
    
    print(f"Loaded dataset. Shape: {df.shape}")
    print(df.head())
    return df

# Step 2: Analyze coexpression data for specified genes
def analyze_expression(df, genes, output_path=None):
    print(f"\nAnalyzing coexpression for genes: {genes}")
    
    # Check for missing genes
    missing_genes = [gene for gene in genes if gene not in df.index]
    if missing_genes:
        print(f"WARNING: {missing_genes} are missing from the dataset")
    
    # Filter df for only available genes
    available_genes = [gene for gene in genes if gene in df.index]
    if available_genes:
        gene_data = df.loc[available_genes, :]
        print("\nExpression data for available genes:")
        print(gene_data)
        
        # Export the filtered data to a CSV file if `output_path` is provided
        if output_path:
            gene_data.to_csv(output_path)
            print(f"\nFiltered gene data exported to: {output_path}")
    else:
        print("\nNo genes of interest found in the dataset.")
        
      
def analyze_apoe_prdm16_correlation(df, output_path="apoe_prdm16_plot.png"):
    # Extract APOE and PRDM16 expression levels
    apoe_expression = df.loc["APOE"]
    prdm16_expression = df.loc["PRDM16"]

    # Create filtered dataset where both genes are non-zero
    apoe_nonzero = apoe_expression != 0
    prdm16_nonzero = prdm16_expression != 0
    filtered_indices = apoe_nonzero & prdm16_nonzero  # Boolean mask for rows where both genes != 0

    # Filter the data for the second subplot
    apoe_filtered = apoe_expression[filtered_indices]
    prdm16_filtered = prdm16_expression[filtered_indices]

    # Counts for legend
    apoe_zero_count = (apoe_expression == 0).sum()
    prdm16_zero_count = (prdm16_expression == 0).sum()
    apoe_0_prdm16_nonzero = ((apoe_expression == 0) & (prdm16_expression != 0)).sum()
    prdm16_0_apoe_nonzero = ((prdm16_expression == 0) & (apoe_expression != 0)).sum()
    correlation, p_value = stats.pearsonr(apoe_expression, prdm16_expression)

    # Create subplots
    fig, axs = plt.subplots(1, 2, figsize=(16, 6))

    ### All Data Scatter Plot ###
    axs[0].scatter(apoe_expression, prdm16_expression, color='blue', alpha=0.7)
    axs[0].set_title("APOE vs PRDM16 (All Data)")
    axs[0].set_xlabel("APOE Expression (TPM)")
    axs[0].set_ylabel("PRDM16 Expression (TPM)")
    axs[0].grid(True)

    ### Filtered Data Scatter Plot (Non-Zero Expressions) ###
    axs[1].scatter(apoe_filtered, prdm16_filtered, color='green', alpha=0.7)
    axs[1].set_title("APOE vs PRDM16 (Filtered: Both Non-Zero)")
    axs[1].set_xlabel("APOE Expression (TPM)")
    axs[1].set_ylabel("PRDM16 Expression (TPM)")
    axs[1].grid(True)

    # Add legends with additional details
    legend_text = f"""
    Dataset Summary:
    - Total zeros in APOE: {apoe_zero_count}
    - Total zeros in PRDM16: {prdm16_zero_count}
    - Cases APOE = 0, PRDM16 ≠ 0: {apoe_0_prdm16_nonzero}
    - Cases PRDM16 = 0, APOE ≠ 0: {prdm16_0_apoe_nonzero}
    - Pearson Correlation (All Data): {correlation:.2f}, p-value: {p_value:.2e}
    """
    fig.text(0.5, -0.1, legend_text, ha='center', fontsize=10)

    # Save the figure to file
    plt.tight_layout()
    plt.savefig(output_path, bbox_inches="tight", dpi=300)
    print(f"Figure saved to {output_path}")

    # Show the figure
    plt.show()
    

def plot_pairwise_scatter(df, genes, output_path="pairwise_relationships.png"):
    # Filter genes of interest and transpose dataset
    gene_data = df.loc[genes]
    gene_data_transposed = gene_data.T

    # Filter rows where both APOE and PRDM16 are non-zero
    filtered_data = gene_data_transposed[(gene_data_transposed['APOE'] != 0) & (gene_data_transposed['PRDM16'] != 0)]

    ### Pairplot for All Data ###
    pairplot_all = sns.pairplot(gene_data_transposed)
    pairplot_all.fig.suptitle("Pairwise Relationships (All Data)", y=1.03)  # Add title to the pairplot
    
    # Export the first pairplot
    pairplot_all.savefig("pairwise_all_data.png", dpi=300)
    print("Pairwise plot for ALL data saved as 'pairwise_all_data.png'")

    ### Pairplot for Filtered Non-Zero Data ###
    pairplot_filtered = sns.pairplot(filtered_data)
    pairplot_filtered.fig.suptitle("Pairwise Relationships (Filtered Data: Both Non-Zero)", y=1.03)  # Add title to the pairplot

    # Export the second pairplot
    pairplot_filtered.savefig("pairwise_filtered_data.png", dpi=300)
    print("Pairwise plot for FILTERED (non-zero) data saved as 'pairwise_filtered_data.png'")

    # Option to show both plots
    plt.show()
    
    


def main():
    print("Program started...")
  
    # Step 1: Load the dataset directly from the `.gz` file
    df = load_data()
    print("Data loaded...")
    
    # Step 2: Analyze specific genes of interest
    primary_genes = ["PRDM16", "APOE", "GFRA2"]
    synaptic_genes = ["DLG4", "SYN1", "NRXN1", "SLC17A7", "VGLUT1", "PSD95"]
    
    genes_of_interest = primary_genes
      # would like to concatenate with synaptic genes later
      # PSD95 missing from the dataset
    
    # Specify output file path for export
    OUTPUT_PATH = os.path.join('data', 'filtered_gene_data.csv')
    analyze_expression(df, genes_of_interest)

    plot_pairwise_scatter(df, genes_of_interest, output_path="pairwise_relationships.png")
    
    

if __name__ == "__main__":
    main()
