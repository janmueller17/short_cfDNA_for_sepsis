import matplotlib.pyplot as plt
import pandas as pd
import os
import sys

def validate_arguments():
    if len(sys.argv) != 2:
        print("Usage: python plot_piecharts.py <path_to_annotation_proportions_file>")
        sys.exit(1)

def validate_file(in_file):
    if not os.path.isfile(in_file):
        print(f"Error: File {in_file} does not exist.")
        sys.exit(1)
    try:
        annotation_proportions = pd.read_csv(in_file, sep='\t', header=None)
    except Exception as e:
        print(f"Error reading the file: {e}")
        sys.exit(1)
    return annotation_proportions

def plot_piecharts(ax, row, col, labels, sizes, explode, colors, distances):
    ax[row, col].pie(sizes, labels=labels, autopct='%1.1f%%', startangle=90, explode=explode, colors=colors, pctdistance=distances[0], labeldistance=distances[1])

def main():
    validate_arguments()
    in_file = sys.argv[1]
    in_name = os.path.splitext(os.path.basename(in_file))[0]
    in_path = os.path.dirname(in_file)
    
    annotation_proportions = validate_file(in_file)
    
    # Print the structure of the DataFrame for debugging
    print(annotation_proportions)

    fig, ax = plt.subplots(nrows=3, ncols=2, figsize=(10, 15))
    myexplode = [0, 0.1]
    mycolors = ["#999999", "#333333"]
    distances = [1.5, 1, 1]
    
    pairs = [
        ("Annotated", "Non-annotated"),
        ("TFBS", "Non-TFBS"),
        ("Genes", "Non-genes"),
        ("CRE", "Non-CRE"),
        ("CpG", "Non-CpG")
    ]

    for index, (row, col) in enumerate([(0, 0), (1, 0), (1, 1), (2, 0), (2, 1)]):
        label1, label2 = pairs[index]
        labels = [label1, label2]
        sizes = annotation_proportions[annotation_proportions[0].isin(labels)][1].values
        plot_piecharts(ax, row, col, labels, sizes, myexplode, mycolors, distances)
    
    ax[0, 1].axis('off')
    
    fig.savefig(os.path.join(in_path, f"{in_name}_piechart.pdf"), bbox_inches="tight")

if __name__ == "__main__":
    main()