#!/usr/bin/env python3

"""
This script processes clinical and read count data to create linear regression models
predicting clinical metrics based on read count signals. The results are plotted and
saved as a PDF, showing the performance of the models for both identification and test
samples.

Usage:
python linear_model_clinical_predictor.py --clinical_metrics METRIC1 METRIC2 \
                    --base_path /path/to/output \
                    --test_samples_csv /path/to/test_samples.csv \
                    --ident_samples_csv /path/to/ident_samples.csv \
                    --ident_avg_correlation_csv /path/to/ident_avg_correlation.csv \
                    --ident_adj_p_csv /path/to/ident_adj_p.csv \
                    --ident_comb_corr_csv /path/to/ident_comb_corr.csv \
                    --clinical_metrics_csv /path/to/clinical_metrics.csv \
                    --alpha 0.05 \
                    --test_samples_count_mat_csv /path/to/test_samples_count_matrix.csv \
                    --ident_samples_count_mat_csv /path/to/ident_samples_count_matrix.csv
"""

import os
import pandas as pd
import scipy.stats as stats
import numpy as np
from sklearn.linear_model import LinearRegression
from sklearn.metrics import mean_absolute_error
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import argparse

def parse_arguments():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser()
    parser.add_argument("--clinical_metrics", help="Clinical metric(s) to test", type=str, nargs='+')
    parser.add_argument("--base_path", help="Path where the output should be directed", type=str)
    parser.add_argument("--test_samples_csv", help="CSV file with test samples metadata", type=str)
    parser.add_argument("--ident_samples_csv", help="CSV file with ident samples metadata", type=str)
    parser.add_argument("--ident_avg_correlation_csv", help="CSV file with average correlation from patients", type=str)
    parser.add_argument("--ident_adj_p_csv", help="CSV file with adjusted p-values of combined correlation", type=str)
    parser.add_argument("--ident_comb_corr_csv", help="CSV file with correlation values of combined correlation", type=str)
    parser.add_argument("--clinical_metrics_csv", help="CSV file with clinical metrics values", type=str)
    parser.add_argument("--alpha", help="Significance threshold", type=float)
    parser.add_argument("--test_samples_count_mat_csv", help="CSV file with read count matrix of test samples", type=str)
    parser.add_argument("--ident_samples_count_mat_csv", help="CSV file with read count matrix of ident samples", type=str)
    return parser.parse_args()

def load_data(args):
    """Load data from CSV files."""
    ident_avg_correlation = pd.read_csv(args.ident_avg_correlation_csv)
    ident_adj_p = pd.read_csv(args.ident_adj_p_csv)
    ident_comb_corr = pd.read_csv(args.ident_comb_corr_csv)
    clinical_metrics = pd.read_csv(args.clinical_metrics_csv)
    test_samples_metadata = pd.read_csv(args.test_samples_csv)
    ident_samples_metadata = pd.read_csv(args.ident_samples_csv)
    test_samples_count_signal_df = pd.read_csv(args.test_samples_count_mat_csv).transpose()
    ident_samples_count_signal_df = pd.read_csv(args.ident_samples_count_mat_csv).transpose()

    # Transpose and set correct headers
    test_samples_count_signal_df.columns = test_samples_count_signal_df.iloc[0]
    test_samples_count_signal_df = test_samples_count_signal_df[1:]
    test_samples_count_signal_df.sort_index(axis=1, inplace=True)

    ident_samples_count_signal_df.columns = ident_samples_count_signal_df.iloc[0]
    ident_samples_count_signal_df = ident_samples_count_signal_df[1:]
    ident_samples_count_signal_df.sort_index(axis=1, inplace=True)

    return (ident_avg_correlation, ident_adj_p, ident_comb_corr, clinical_metrics, 
            test_samples_metadata, ident_samples_metadata, 
            test_samples_count_signal_df, ident_samples_count_signal_df)

def merge_clinical_metrics(clinical_metrics, samples_metadata):
    """Merge clinical metrics with sample metadata."""
    merged_df = clinical_metrics.merge(samples_metadata, on=['Seq_ID', 'Day'], how='right')
    merged_df.set_index('Seq_ID', inplace=True)
    return merged_df

def filter_test_samples(test_samples_clinical_metrics, ident_samples_clinical_metrics, clinical_metric):
    """Filter test samples with clinical values within the data range of ident samples."""
    clinical_min_ident = ident_samples_clinical_metrics[clinical_metric].min()
    clinical_max_ident = ident_samples_clinical_metrics[clinical_metric].max()
    filtered_test_samples = test_samples_clinical_metrics[
        (test_samples_clinical_metrics[clinical_metric] >= clinical_min_ident) & 
        (test_samples_clinical_metrics[clinical_metric] <= clinical_max_ident)
    ]
    return filtered_test_samples

def get_selected_loci(df_comb_filt, metric_of_interest):
    """Select loci based on significance threshold and sort them."""
    if df_comb_filt.shape[0] >= 50:
        selected_loci = df_comb_filt.sort_values(by=metric_of_interest, na_position='last', ascending=False).head(50).Locus.tolist()
        print(f"More than 50 significant loci found for {metric_of_interest}.")
        
    elif df_comb_filt.shape[0] == 0:
        print(f"No significant loci found for {metric_of_interest}.")
        return []
    else:
        selected_loci = df_comb_filt.sort_values(by=metric_of_interest, na_position='last', ascending=False).Locus.tolist()
        print(f"Less than 50 significant loci found for {metric_of_interest}.")
    return selected_loci

def create_linear_model_and_plots(ident_samples_clinical_metrics, ident_samples_count_signal_df, 
                                  test_samples_clinical_metrics, test_samples_count_signal_df, 
                                  selected_loci, clinical_metric_of_interest, count_signal_name, base_path):
    """Create linear models, predictions, and plots."""
    output_pdf_path = os.path.join(base_path, f"{count_signal_name}_{clinical_metric_of_interest}_LinMod_predict.pdf")
    pdf_pages = PdfPages(output_pdf_path)

    for count_signal_of_interest in selected_loci:
        ident_true_df = pd.DataFrame(ident_samples_clinical_metrics[clinical_metric_of_interest]).merge(
            pd.DataFrame(ident_samples_count_signal_df[count_signal_of_interest]), how='right', left_index=True, right_index=True)
        ident_true_df = ident_true_df.dropna()

        test_true_df = pd.DataFrame(test_samples_clinical_metrics[clinical_metric_of_interest]).merge(
            pd.DataFrame(test_samples_count_signal_df[count_signal_of_interest]), how='right', left_index=True, right_index=True)
        test_true_df = test_true_df.dropna()

        ident_true_clin = np.array([ident_true_df[clinical_metric_of_interest].sort_index()]).reshape((-1, 1))
        ident_true_count = np.array([ident_true_df[count_signal_of_interest].sort_index()]).reshape((-1, 1))
        test_true_clin = np.array([test_true_df[clinical_metric_of_interest].sort_index()]).reshape((-1, 1))
        test_true_count = np.array([test_true_df[count_signal_of_interest].sort_index()]).reshape((-1, 1))

        model = LinearRegression().fit(X=ident_true_count, y=ident_true_clin)
        slope = model.coef_[0][0]
        intercept = model.intercept_[0]

        ident_pred_clin = model.predict(X=ident_true_count)
        ident_r = stats.pearsonr(y=ident_true_clin.ravel(), x=ident_true_count.ravel())[0]
        ident_mae = mean_absolute_error(y_true=ident_true_clin, y_pred=ident_pred_clin)

        test_pred_clin = model.predict(X=test_true_count)
        test_r = stats.pearsonr(y=test_true_clin.ravel(), x=test_true_count.ravel())[0]
        test_mae = mean_absolute_error(y_true=test_true_clin, y_pred=test_pred_clin)

        plot_results(ident_true_clin, ident_true_count, test_true_clin, test_true_count, 
                     test_pred_clin, slope, intercept, ident_r, ident_mae, test_r, test_mae, 
                     clinical_metric_of_interest, count_signal_name, count_signal_of_interest, pdf_pages)

    pdf_pages.close()

def plot_results(ident_true_clin, ident_true_count, test_true_clin, test_true_count, 
                 test_pred_clin, slope, intercept, ident_r, ident_mae, test_r, test_mae, 
                 clinical_metric_of_interest, count_signal_name, count_signal_of_interest, pdf_pages):
    """Plot results and save to PDF."""
    clinical_max = float(np.max(np.concatenate([ident_true_clin, test_true_clin, test_pred_clin]), axis=0)[0])
    clinical_min = float(np.min(np.concatenate([ident_true_clin, test_true_clin, test_pred_clin]), axis=0)[0])
    count_max = float(np.max(np.concatenate([ident_true_count, test_true_count]), axis=0)[0])
    count_min = float(np.min(np.concatenate([ident_true_count, test_true_count]), axis=0)[0])
    buffer = 0.1
    count_buffer_size = buffer * (count_max - count_min)
    clinical_buffer_size = buffer * (clinical_max - clinical_min)

    fig, axes = plt.subplots(1, 2, figsize=(14 / 2.54, 7 / 2.54))
    fig.suptitle(f'{count_signal_name}\n{count_signal_of_interest}')

    axes[0].scatter(y=ident_true_clin.tolist(), x=ident_true_count.tolist(), color='black', s=10)
    axes[0].plot(ident_true_count.tolist(), slope * np.array(ident_true_count.tolist()) + intercept, color='darkred', label='Linear Model')
    axes[0].set_ylabel(clinical_metric_of_interest)
    axes[0].set_xlabel("Normalized readcount")
    axes[0].set_xlim(count_min - count_buffer_size, count_max + count_buffer_size)
    axes[0].set_ylim(clinical_min - clinical_buffer_size, clinical_max + clinical_buffer_size)
    axes[0].xaxis.label.set_fontweight('bold')
    axes[0].yaxis.label.set_fontweight('bold')
    axes[0].set_title(f'Ident\nr = {ident_r:.3f}, MAE = {ident_mae:.1f}')

    axes[1].scatter(y=test_true_clin.tolist(), x=test_true_count.tolist(), color='darkblue', s=10, label='Actual')
    axes[1].plot(test_true_count.tolist(), slope * np.array(test_true_count.tolist()) + intercept, color='darkred', label='Predicted')
    axes[1].set_ylabel(clinical_metric_of_interest)
    axes[1].set_xlabel("Normalized readcount")
    axes[1].set_xlim(count_min - count_buffer_size, count_max + count_buffer_size)
    axes[1].set_ylim(clinical_min - clinical_buffer_size, clinical_max + clinical_buffer_size)
    axes[1].xaxis.label.set_fontweight('bold')
    axes[1].yaxis.label.set_fontweight('bold')
    axes[1].set_title(f'Test\nr = {test_r:.3f}, MAE = {test_mae:.1f}')

    plt.tight_layout()
    pdf_pages.savefig(fig, bbox_inches='tight')
    plt.close(fig)

def main():
    """Main function to execute the script."""
    args = parse_arguments()
    (ident_avg_correlation, ident_adj_p, ident_comb_corr, clinical_metrics, 
     test_samples_metadata, ident_samples_metadata, 
     test_samples_count_signal_df, ident_samples_count_signal_df) = load_data(args)

    count_signal_name = os.path.basename(args.test_samples_count_mat_csv).split('.')[0]

    # Combine dataframes for computation
    df_comb = ident_avg_correlation + ident_comb_corr
    df_comb['Locus'] = ident_comb_corr.Locus

    test_samples_clinical_metrics = merge_clinical_metrics(clinical_metrics, test_samples_metadata)
    ident_samples_clinical_metrics = merge_clinical_metrics(clinical_metrics, ident_samples_metadata)

    # Split the clinical metrics string into individual metrics
    if args.clinical_metrics and len(args.clinical_metrics) == 1:
        args.clinical_metrics = args.clinical_metrics[0].split()

    for clinical_metric_of_interest in args.clinical_metrics:
        filtered_test_samples_clinical_metrics = filter_test_samples(test_samples_clinical_metrics, ident_samples_clinical_metrics, clinical_metric_of_interest)

        df_comb_filt = df_comb[ident_adj_p[clinical_metric_of_interest] <= args.alpha]
        df_comb_filt.loc[:, clinical_metric_of_interest] = df_comb_filt[clinical_metric_of_interest].abs()
        
        selected_loci = get_selected_loci(df_comb_filt, clinical_metric_of_interest)
        if selected_loci:
            create_linear_model_and_plots(ident_samples_clinical_metrics, ident_samples_count_signal_df, 
                                          filtered_test_samples_clinical_metrics, test_samples_count_signal_df, 
                                          selected_loci, clinical_metric_of_interest, count_signal_name, args.base_path)

if __name__ == "__main__":
    main()