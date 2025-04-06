import numpy as np
import matplotlib.pyplot as plt
import os
import sys
import json
import seaborn as sns
import pandas as pd
from functools import reduce
import utils_matanya as um
import subprocess
import argparse


def parse_args():
    parser = argparse.ArgumentParser(description="create heatmaps and clustermaps from a CSV file for overall TF binding in Loci.")
    parser.add_argument("--csv_file", type=str, required=True, help="Path to the CSV file.")
    parser.add_argument("--output_dir", type=str, required=True, help="Directory to save the output files.")
    parser.add_argument("--file_name", type=str, required=True, help="name of the file to save.")
    parser.add_argument("--title", type=str, required=True, help="title of the plot.")

    parser.add_argument("--heatmap", type=bool, default=False, help="Should we create a heatmap?")
    parser.add_argument("--clustermap", type=bool, default=False, help="Should we create a clustermap?")
    return parser.parse_args()

def check_args(args):
    if not os.path.exists(args.csv_file):
        raise FileNotFoundError(f"{args.mpra_file} does not exist.")
    if not args.heatmap and not args.clustermap:
        raise ValueError("You must choose at least one of heatmap or clustermap.")


def create_overall_binding_heatmap(csv_file, title, output_dir):
    sns.set(style="whitegrid")
    df = pd.read_csv(csv_file, index_col=0)

    df_clipped = np.clip(df, -10, 10)

    # --- Plot & Save Derived Heatmap ---
    plt.figure(figsize=(16, 10))
    sns.heatmap(
        df_clipped,
        cmap="coolwarm",
        center=0,
        xticklabels=False,
        # linewidths=0.5,
        # linecolor='gray'
    )
    plt.title(title, fontsize=14)
    plt.xlabel("Locus")
    plt.ylabel("Transcription Factor")
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, f"{title.replace(' ' , '_')}.png"))
    plt.close()

def create_overall_binding_clustermap(csv_file, title, output_dir, file_name):
    sns.set(style="white")
    df = pd.read_csv(csv_file, index_col=0)

    # Create clustered heatmap
    clustermap = sns.clustermap(
        df,              # TFs = rows, loci = columns
        cmap="coolwarm",
        center=0,
        figsize=(16, 10),
        xticklabels=False,        # Hide x labels if too many loci
        yticklabels=True,         # Show TF names
        cbar_kws={'label': 'Binding Strength'},
        method='average',         # Clustering method (can be 'ward', 'single', etc.)
        metric='euclidean'        # Distance metric
    )
    

    # clustermap.fig.suptitle("Clustermap of Max TF Binding (All Loci)", fontsize=16)
    clustermap.fig.suptitle(title, fontsize=16)

    # Save the figure,,
    clustermap.savefig(os.path.join(output_dir, f"{file_name}.png"))
    plt.close()


def main():
    args = parse_args()
    check_args(args)

    # Read the CSV file
    # df = pd.read_csv(args.csv_file, index_col=0)

    # Create the output directory if it doesn't exist
    os.makedirs(args.output_dir, exist_ok=True)

    # Create heatmap and clustermap if specified
    if args.heatmap:
        create_overall_binding_heatmap(args.csv_file, args.title, args.output_dir)
    if args.clustermap:
        create_overall_binding_clustermap(args.csv_file, args.title, args.output_dir, args.file_name)
    print(f"Plots saved in {args.output_dir}")
    print("Done!")
if __name__ == "__main__":
    main()