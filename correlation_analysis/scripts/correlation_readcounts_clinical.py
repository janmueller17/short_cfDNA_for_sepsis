#!/usr/bin/env python3

"""
This script calculates correlations between count matrix data and clinical metrics for timecourse samples of patients (individually, averaged, and combined).

Usage: python matrix_clinical_correlation.py -t <threads> -m <metadata.csv> -p <patients.csv> -c pearson -i <counts_matrix.csv> -v Leukocytes Erythrocytes Creatinine Bilirubin Albumine Urea ALT AST ALP Hemoglobine Thrombocytes Quick aPTT CRP PCT -o <output_directory>
"""

import pandas as pd
import scipy.stats as stats
import numpy as np
import multiprocessing as mp
import argparse
import os

def parse_args():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(description="Calculate correlations between counts matrix data and clinical metrics for patients and average the results.")
    parser.add_argument("-t", "--threads", help="Number of threads for parallel correlation calculation", type=int, required=True)
    parser.add_argument("-m", "--metadata", help="Path to the metadata CSV file", type=str, required=True)
    parser.add_argument("-p", "--patients", help="Path to the CSV file containing patient information", type=str, required=True)
    parser.add_argument("-c", "--CorrMethod", help="Correlation method: pearson or spearman", type=str, choices=['pearson', 'spearman'], required=True)
    parser.add_argument("-i", "--input", help="Path to the counts matrix CSV file", type=str, required=True)
    parser.add_argument("-v", "--variables", type=str, help="List of clinical metrics", required=True)
    parser.add_argument("-o", "--output", help="Base path for output files", type=str, required=True)
    return parser.parse_args()

def read_data(args):
    """Read metadata and patient information CSVs."""
    metadata = pd.read_csv(args.metadata)
    patients_df = pd.read_csv(args.patients)
    patients_to_analyze = patients_df['Individual'].unique().tolist()

    print("Metadata:")
    print(metadata.head())
    print("Patients DataFrame:")
    print(patients_df.head())

    count_data = pd.read_csv(args.input).transpose()
    count_data.columns = count_data.iloc[0]
    count_data = count_data[1:]

    print("Transposed Count Data (first 5 rows):")
    print(count_data.head())

    count_data_dict = {}
    for patient in patients_to_analyze:
        patient_seq_ids = patients_df[patients_df['Individual'] == patient]['Seq_ID']
        patient_days = patients_df[patients_df['Individual'] == patient]['Day'].values
        patient_count_data = count_data[count_data.index.isin(patient_seq_ids)]
        patient_count_data.index = patient_days
        count_data_dict[patient] = patient_count_data

        print(f"Processed count data for patient {patient}:")
        print(patient_count_data)

    return metadata, patients_df, patients_to_analyze, count_data_dict

def subset_clinical_data(metadata, patients_to_analyze, variables):
    """Subset clinical data and format for correlation calculation."""
    temp_clinical_metrics = metadata[metadata.Pat_ID.isin(patients_to_analyze)]
    temp_clinical_metrics = temp_clinical_metrics.set_index("Day")
    return temp_clinical_metrics[["Pat_ID"] + variables]

def prepare_data_combined(patients_to_analyze, count_data_dict, temp_clinical_metrics):
    """Prepare combined data for all patients."""
    clinical_df = pd.DataFrame(columns=temp_clinical_metrics.columns)
    count_df = pd.DataFrame(columns=count_data_dict[patients_to_analyze[0]].columns)

    for patient in patients_to_analyze:
        temp_clinical_df = temp_clinical_metrics[temp_clinical_metrics.Pat_ID == patient]
        missing_days = set(count_data_dict[patient].index) - set(temp_clinical_df.index)
        if missing_days:
            print(f"Missing days for patient {patient}: {missing_days}")
            count_data_dict[patient] = count_data_dict[patient].drop(missing_days, errors='ignore')
        
        temp_clinical_df = temp_clinical_df.loc[count_data_dict[patient].index].sort_index()
        count_data_dict[patient] = count_data_dict[patient].sort_index()
        
        if not temp_clinical_df.empty:
            clinical_df = pd.concat([clinical_df, temp_clinical_df], ignore_index=True)
        if not count_data_dict[patient].empty:
            count_df = pd.concat([count_df, count_data_dict[patient]], ignore_index=True)

    return clinical_df.drop("Pat_ID", axis=1), count_df

def prepare_data_individual(patient, count_data_dict, temp_clinical_metrics):
    """Prepare data for an individual patient."""
    temp_clinical_df = temp_clinical_metrics[temp_clinical_metrics.Pat_ID == patient]
    missing_days = set(count_data_dict[patient].index) - set(temp_clinical_df.index)
    if missing_days:
        print(f"Missing days for patient {patient}: {missing_days}")
        count_data_dict[patient] = count_data_dict[patient].drop(missing_days, errors='ignore')
    
    temp_clinical_df = temp_clinical_df.loc[count_data_dict[patient].index].sort_index()
    count_data_dict[patient] = count_data_dict[patient].sort_index()
    
    return temp_clinical_df.drop("Pat_ID", axis=1), count_data_dict[patient]

def calculate(method, count_df, temp_clinical_df, total_it_j, i):
    """Calculate correlations for a single column."""
    temp_count_df = count_df.iloc[:, i]
    corr_val_list = []
    corr_p_list = []
    
    for j in range(total_it_j):
        temp_clinical_df_j = temp_clinical_df.iloc[:, j].dropna()
        temp_count_df_j = temp_count_df.loc[temp_clinical_df_j.index]
        if method == "pearson":
            temp_corr, temp_p = stats.pearsonr(temp_clinical_df_j, temp_count_df_j)
        elif method == "spearman":
            temp_corr, temp_p = stats.spearmanr(temp_clinical_df_j, temp_count_df_j)
        corr_val_list.append(temp_corr)
        corr_p_list.append(temp_p)
    
    return corr_val_list, corr_p_list

def calculate_correlations(method, count_df, temp_clinical_df, total_it_i, total_it_j, threads):
    """Calculate correlations using the specified method."""
    with mp.Pool(threads) as pool:
        results = pool.starmap(calculate, [(method, count_df, temp_clinical_df, total_it_j, i) for i in range(total_it_i)])
    return results

def main():
    args = parse_args()
    output_basename = os.path.splitext(os.path.basename(args.input))[0]
    variables = args.variables.split()
    print("Clinical Metrics:", variables)
    metadata, patients_df, patients_to_analyze, count_data_dict = read_data(args)
    temp_clinical_metrics = subset_clinical_data(metadata, patients_to_analyze, variables)

    # Combined analysis for all patients
    temp_clinical_df_combined, count_df_combined = prepare_data_combined(patients_to_analyze, count_data_dict, temp_clinical_metrics)
    total_it_i_combined = count_df_combined.shape[1]
    total_it_j_combined = temp_clinical_df_combined.shape[1]

    results_combined = calculate_correlations(args.CorrMethod, count_df_combined, temp_clinical_df_combined, total_it_i_combined, total_it_j_combined, args.threads)

    corr_val_arr_combined = np.array([result[0] for result in results_combined])
    corr_p_arr_combined = np.array([result[1] for result in results_combined])

    corr_val_df_combined = pd.DataFrame(data=corr_val_arr_combined, index=count_df_combined.columns, columns=temp_clinical_df_combined.columns).rename_axis("Locus").round(5)
    corr_p_df_combined = pd.DataFrame(data=corr_p_arr_combined, index=count_df_combined.columns, columns=temp_clinical_df_combined.columns).rename_axis("Locus")
    corr_log_p_df_combined = -np.log10(corr_p_df_combined).round(5)

    corr_val_df_combined.to_csv(f"{args.output}/{output_basename}_combined_{args.CorrMethod}_corr_vals.csv", index=True)
    corr_log_p_df_combined.to_csv(f"{args.output}/{output_basename}_combined_{args.CorrMethod}_corr_p.csv", index=True)
    print("Outputting combined results")

    # Individual analysis for each patient
    individ_corr_val_dict = {}
    for patient in patients_to_analyze:
        temp_clinical_df_individual, count_df_individual = prepare_data_individual(patient, count_data_dict, temp_clinical_metrics)

        total_it_i_individual = count_df_individual.shape[1]
        total_it_j_individual = temp_clinical_df_individual.shape[1]

        results_individual = calculate_correlations(args.CorrMethod, count_df_individual, temp_clinical_df_individual, total_it_i_individual, total_it_j_individual, args.threads)

        corr_val_arr_individual = np.array([result[0] for result in results_individual])
        corr_p_arr_individual = np.array([result[1] for result in results_individual])

        corr_val_df_individual = pd.DataFrame(data=corr_val_arr_individual, index=count_df_individual.columns, columns=temp_clinical_df_individual.columns).rename_axis("Locus").round(5)
        corr_p_df_individual = pd.DataFrame(data=corr_p_arr_individual, index=count_df_individual.columns, columns=temp_clinical_df_individual.columns).rename_axis("Locus")
        corr_log_p_df_individual = -np.log10(corr_p_df_individual).round(5)

        corr_val_df_individual.to_csv(f"{args.output}/{output_basename}_{patient}_{args.CorrMethod}_corr_vals.csv", index=True)
        corr_log_p_df_individual.to_csv(f"{args.output}/{output_basename}_{patient}_{args.CorrMethod}_corr_p.csv", index=True)
        print(f"Outputting results for patient {patient}")

        # Store individual results in a dictionary
        individ_corr_val_dict[patient] = corr_val_df_individual

    # Averaging individual results
    averaged_corr_val_df = sum(df for df in individ_corr_val_dict.values()) / len(individ_corr_val_dict)
    averaged_corr_val_df.to_csv(f"{args.output}/{output_basename}_{args.CorrMethod}_individuals_avg_corr_vals.csv", index=True)
    print("Outputting averaged individual results")

if __name__ == "__main__":
    main()