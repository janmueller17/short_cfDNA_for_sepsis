#!/usr/bin/env python3

"""
analyze_correlations.py

This script analyzes the correlations calculated by a previous script. It performs multiple testing adjustments,
plots distributions of p-values and correlation values, and filters results based on a specified p-value threshold.

Usage:
    python analyze_correlations.py -a <alpha> -v <corr_vals.csv> -p <corr_p.csv>
"""

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import argparse
from statsmodels.stats.multitest import multipletests
import os

def read_csv(file_path, index_col=0):
    """Read a CSV file into a DataFrame."""
    if not os.path.exists(file_path):
        raise FileNotFoundError(f"File not found: {file_path}")
    return pd.read_csv(file_path, index_col=index_col)

def adjust_p_values(df):
    """Adjust p-values for multiple testing using the Benjamini-Hochberg method."""
    adj_df = pd.DataFrame()
    for column in df.columns:
        pvalues = df[column].dropna().values
        adjusted_values = multipletests(pvalues, method='fdr_bh')[1]
        adjusted_column = pd.Series(np.nan, index=df[column].index)
        adjusted_column[~df[column].isna()] = adjusted_values
        adj_df[column] = adjusted_column
    return adj_df

def plot_histogram(df, output_path, log=False, range=None, bins=20):
    """Plot histogram of the DataFrame values and save to a file."""
    fig, ax = plt.subplots(figsize=(10, 7))
    df.hist(range=range, bins=bins, ax=ax, grid=False, log=log, bottom=1 if log else 0)
    plt.subplots_adjust(top=1, hspace=1, wspace=0.75)
    plt.savefig(output_path, bbox_inches='tight')
    plt.close(fig)

def filter_by_p_value(df, p_val_df, threshold):
    """Filter rows of df based on p-value threshold in p_val_df."""
    passing_filter = p_val_df[(p_val_df <= threshold).any(axis=1)].index
    return df.loc[passing_filter]

def main(args):
    # Base variables
    p_thresh = args.alpha
    corr_vals_file = args.corr_vals
    corr_p_file = args.corr_p

    # Extract base names for output files
    base_dir = os.path.dirname(corr_vals_file)
    output_basename_corr_vals = os.path.splitext(os.path.basename(corr_vals_file))[0]
    output_basename_p_vals = os.path.splitext(os.path.basename(corr_p_file))[0]

    # Read correlation values and p-values
    corr_val_in_df = read_csv(corr_vals_file)
    corr_p_in_df = read_csv(corr_p_file)
    corr_p_in_df = 10**(-corr_p_in_df)  # Convert -log10(p) to p-values

    # Multiple testing adjustment
    corr_p_adj_df = adjust_p_values(corr_p_in_df)

    # Plot distributions
    plot_histogram(corr_p_in_df, f"{base_dir}/{output_basename_p_vals}_distros.pdf", range=(0, 1))
    plot_histogram(corr_p_adj_df, f"{base_dir}/{output_basename_p_vals}_adj_distros.pdf", log=True, range=(0, 1))
    plot_histogram(corr_val_in_df, f"{base_dir}/{output_basename_corr_vals}_distros.pdf", log=True, range=(-1, 1), bins=40)

    # Filtering for adjusted p-value
    corr_val_filt_df = filter_by_p_value(corr_val_in_df, corr_p_adj_df, p_thresh)
    corr_p_filt_df = filter_by_p_value(corr_p_in_df, corr_p_adj_df, p_thresh)
    corr_p_adj_filt_df = filter_by_p_value(corr_p_adj_df, corr_p_adj_df, p_thresh)

    # Output filtered results
    corr_p_adj_df.to_csv(f"{base_dir}/{output_basename_p_vals}_adj.csv", index=True)
    corr_val_filt_df.to_csv(f"{base_dir}/{output_basename_corr_vals}_alpha_{p_thresh}.csv", index=True)
    corr_p_filt_df.to_csv(f"{base_dir}/{output_basename_p_vals}_alpha_{p_thresh}.csv", index=True)
    corr_p_adj_filt_df.to_csv(f"{base_dir}/{output_basename_p_vals}_adj_alpha_{p_thresh}.csv", index=True)

if __name__ == "__main__":
    # Define command line arguments
    parser = argparse.ArgumentParser(description="Analyze correlations from counts matrix data and clinical metrics.")
    parser.add_argument("-a", "--alpha", help="Alpha, p-value threshold", type=float, required=True)
    parser.add_argument("-v", "--corr_vals", help="Path to the correlation values CSV file", type=str, required=True)
    parser.add_argument("-p", "--corr_p", help="Path to the correlation p-values CSV file", type=str, required=True)
    args = parser.parse_args()

    main(args)