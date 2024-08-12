import pandas as pd
import matplotlib.pyplot as plt
import argparse
import matplotlib.cm as cm
import os
import numpy as np
from scipy.stats import pearsonr

def parse_args():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser()

    # Arguments from the first script
    parser.add_argument("--cfoi", help="Clinical feature of interest to be plotted", type=str, nargs='+')
    parser.add_argument("--locus_interest", help="One specific locus to be plotted", type=str)
    parser.add_argument("--count_data_csv", help="Path to the Liquid footprinting signal data CSV file", type=str, required=True)
    parser.add_argument("--samples_metadata", help="Path to the samples metadata CSV file with Seq_IDs", type=str, required=True)
    parser.add_argument("--corr_val_df", help="Path to the combined correlation values CSV file", type=str, required=True)
    parser.add_argument("--corr_p_adj_df", help="Path to the adjusted p-values CSV file", type=str, required=True)
    parser.add_argument("--clinical_metrics", help="Path to the samples metadata CSV file with clinical feature values", type=str, required=True)
    parser.add_argument("--output_base_path", help="Base path for the output files", type=str, required=True)

    # Additional argument for specific patient plotting
    parser.add_argument("-p", "--patient", help="One specific patient ID to be analyzed", type=str)

    return parser.parse_args()

def read_data(args):
    """Read and process data."""
    samples_metadata = pd.read_csv(args.samples_metadata)
    patients_to_analyze = samples_metadata['Individual'].unique().tolist()

    corr_val_df = pd.read_csv(args.corr_val_df, index_col=0)
    corr_p_adj_df = pd.read_csv(args.corr_p_adj_df, index_col=0)

    count_data = pd.read_csv(args.count_data_csv).transpose()
    count_data.columns = count_data.iloc[0]
    count_data = count_data[1:]
    count_data.index.name = 'Seq_ID'
    count_data = count_data.reset_index()

    clinical_metrics = pd.read_csv(args.clinical_metrics)

    return corr_val_df, corr_p_adj_df, count_data, clinical_metrics, samples_metadata, patients_to_analyze

def plot_patient_specific(args, corr_val_df, count_data, clinical_metrics, samples_metadata, highest):
    """Plot data for a specific patient."""
    patient_to_analyze = args.patient
    features_of_interest = args.cfoi

    # Filter data for the specific patient
    temp_clinical_metrics = clinical_metrics[clinical_metrics.Pat_ID == patient_to_analyze]
    temp_clinical_metrics = temp_clinical_metrics.set_index("Day")
    temp_clinical_metrics = temp_clinical_metrics.loc[:, temp_clinical_metrics.columns.isin(corr_val_df.columns)]

    # Map Seq_ID to Day for the patient
    sample_days = samples_metadata[samples_metadata["Individual"] == patient_to_analyze].set_index("Seq_ID")["Day"]
    count_data = count_data[count_data["Seq_ID"].isin(sample_days.index)].set_index("Seq_ID")
    count_data.index = count_data.index.map(sample_days)
    count_data = count_data.sort_index()

    temp_clinical_metrics = temp_clinical_metrics.loc[count_data.index]
    temp_clinical_metrics = temp_clinical_metrics.sort_index()

    if (count_data.index == temp_clinical_metrics.index).all():
        font = {'fontname': 'Arial'}
        for feature in features_of_interest:
            selected_corr = count_data.loc[:, count_data.columns.isin(highest.index)]

            fig, axes = plt.subplots(1, 3, figsize=(7, 3))

            axes[0].plot(temp_clinical_metrics[feature].index, temp_clinical_metrics[feature], color="grey")
            axes[0].scatter(temp_clinical_metrics[feature].index, temp_clinical_metrics[feature], color="blue")
            axes[0].set_ylim(bottom=0)
            axes[0].set_xlabel("Days")
            axes[0].set_ylabel(feature)
            axes[0].xaxis.label.set_fontweight('bold')
            axes[0].yaxis.label.set_fontweight('bold')

            axes[1].plot(selected_corr.index, selected_corr.iloc[:, 0], color="grey")
            axes[1].scatter(selected_corr.index, selected_corr.iloc[:, 0], color="red")
            axes[1].set_ylim(bottom=0)
            axes[1].set_xlabel("Days")
            axes[1].set_ylabel("Normalized readcount")
            axes[1].xaxis.label.set_fontweight('bold')
            axes[1].yaxis.label.set_fontweight('bold')
            
            # Recalculate Pearson correlation
            x = temp_clinical_metrics[feature]
            y = selected_corr.iloc[:, 0]
            corr_val, _ = pearsonr(x, y)

            axes[2].scatter(temp_clinical_metrics[feature], selected_corr.iloc[:, 0], color="purple")
            axes[2].set_xlabel(feature)
            axes[2].set_ylabel("Normalized readcount")
            axes[2].xaxis.label.set_fontweight('bold')
            axes[2].yaxis.label.set_fontweight('bold')
            axes[2].set_title("r=" + str(round(corr_val, ndigits=5)))

            fig.suptitle(patient_to_analyze, fontweight='bold')
            plt.tight_layout()
            plt.rcParams['font.family'] = font['fontname']
            
            # Generate the output file name and path
            count_data_basename = os.path.basename(args.count_data_csv).replace('.csv', '')
            output_file_name = f"{count_data_basename}_{feature}_{args.locus_interest}_corr_{patient_to_analyze}_patient.pdf"
            output_file_path = os.path.join(args.output_base_path, output_file_name)
            output_file_path_clean = output_file_path.replace(':', '.')
            
            plt.savefig(output_file_path_clean, format='pdf', transparent=True, bbox_inches='tight')
            plt.close()
    else:
        print("Value pairs are not matched.")

def main():
    args = parse_args()
    corr_val_df, corr_p_adj_df, count_data, clinical_metrics, samples_metadata, patients_to_analyze = read_data(args)

    features_of_interest = args.cfoi  # Mover esta línea fuera del bucle

    for feature in features_of_interest:
        locus_of_interest = args.locus_interest
    
        if locus_of_interest is None:
            highest = corr_val_df[feature].abs().nlargest(1)
        else:
            highest = corr_val_df.loc[[locus_of_interest], feature]
        
        if args.patient:  # Plot for a specific patient if provided
            plot_patient_specific(args, corr_val_df, count_data, clinical_metrics, samples_metadata, highest)
        else:
            temp_clinical_metrics = clinical_metrics[clinical_metrics.Pat_ID.isin(patients_to_analyze)]
            temp_clinical_metrics = temp_clinical_metrics.set_index("Seq_ID")
            temp_clinical_metrics = temp_clinical_metrics[["Pat_ID", feature]]

            common_index = count_data.set_index("Seq_ID").index.intersection(temp_clinical_metrics.index)
            temp_clinical_df = temp_clinical_metrics.loc[common_index].sort_index()
            count_data = count_data.set_index("Seq_ID").loc[common_index].sort_index()

            if locus_of_interest is None:
                highest = corr_val_df[feature].abs().nlargest(1)
            else:
                highest = corr_val_df.loc[[locus_of_interest], feature]

            color_map = plt.get_cmap('Accent')
            font = {'fontname': 'Arial'}

            # Create a colormap for patient IDs
            unique_patients = temp_clinical_df['Pat_ID'].unique()
            colors = color_map(np.linspace(0, 1, len(unique_patients)))
            patient_color_map = {patient: colors[i] for i, patient in enumerate(unique_patients)}

            fig, axes = plt.subplots(figsize=(4, 3))

            selected_corr = count_data.loc[:, count_data.columns.isin(highest.index)]

            if temp_clinical_df[feature].shape[0] != selected_corr.iloc[:, 0].shape[0]:
                print("Error: The sizes of x and y do not match.")
                continue

            # Scatter plot with colors based on patient IDs
            for patient_id in unique_patients:
                patient_data = temp_clinical_df[temp_clinical_df['Pat_ID'] == patient_id]
                axes.scatter(patient_data[feature], selected_corr.loc[patient_data.index, selected_corr.columns[0]],
                             color=patient_color_map[patient_id], label=patient_id)

            axes.set_xlabel(feature, fontweight='bold')
            axes.set_ylabel("Normalized readcount", fontweight='bold')
            axes.set_title(f"r={corr_val_df.loc[selected_corr.columns[0], feature]} p_adj={round(corr_p_adj_df.loc[selected_corr.columns[0], feature], 5)}")

            plt.legend(title='Patient ID')
            plt.tight_layout()
            plt.rcParams['font.family'] = font['fontname']

            # Generate the output file name
            count_data_basename = os.path.basename(args.count_data_csv).replace('.csv', '')
            output_file_name = f"{count_data_basename}_{feature}_{args.locus_interest}_combined_corr_col_by_patient.pdf"
            output_file_path = os.path.join(args.output_base_path, output_file_name)
            output_file_path_clean = output_file_path.replace(':', '.')

            plt.savefig(output_file_path_clean, format='pdf', transparent=True, bbox_inches='tight')
            plt.close()

if __name__ == "__main__":
    main()