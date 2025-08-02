# Gene Expression Connectivity Correlation tool

This script analyzes the relationship between the genes, SOD1 and C9orf72 values, structural brain connectivity in 5 regions of the brain, and regional brain atrophy rates of 5 regions of the brain.
This script is intended for use in neurodegenerative diseases like ALS and others. 

# Parts of the script:
-  Loads 3 csv files
    -    Gene expression data between values of SOD1 and C9orf72
    -    A connectome matrix between 5 regions of the brain
    -    Regional brain atrophy rates
-  Extracts, cleans and processes the data
-  Computes correlations:
    -    Gene expression and connectivity
    -    Gene expression and atrophy
-  Visualizes a scatter plot and regression line for the relationships between the correlations
-  Saves the computed data and correlation values to a CSV file

# Input files:
-  Gene Expression CSV
   - Required columns: Region, SOD1, C9orf72
- Connectome CSV
    - Required columns: 5 regions(Region A-E)
    - Required Rows: 5 regions(Region A-E)
- Atrophy Rates CSV
    - Required columns: Region, Atrophy Rate

# Output:
-  Csv filed called scores.csv contains region data, expression values, connectivity strength, and atrophy rates
-  Console:
    - Correlation values between expression/connectivity and expression/atrophy
- Four scatter plots:
    - SOD1 vs. Connectivity Strength
    - C9orf72 vs. Connectivity Strength
    - SOD1 vs. Atrophy Rate
    - C9orf72 vs. Atrophy Rate

# Dependencies:
  pandas numpy matplotlib scipy

# Notes:
- The script assumes clean and complete data
- NaN or missing values are skipped
- Outputs will overwrite existing scores.csv file
