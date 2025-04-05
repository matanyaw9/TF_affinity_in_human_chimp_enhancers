import os
import json
import argparse
import numpy as np
import pandas as pd
import utils_matanya as um

# I made this file to be able to run it as jobs for the wexac :) 
# It is essentially the same as the notebook file

# todo: These will be used for titles and for file names
TF_AFFINITY_INSIDE_ANCESTRAL_SEQ_UNFILTERED = "Unfiltered TF Affinity inside Ancestral Sequence"
TF_AFFINITY_INSIDE_DERIVED_SEQ_UNFILTERED = "Unfiltered TF Affinity inside Derived Sequence"
TF_AFFINITY_INSIDE_ANCESTRAL_SEQ_FILTERED = "Filtered TF Affinity inside Ancestral Sequence"
TF_AFFINITY_INSIDE_DERIVED_SEQ_FILTERED = "Filtered TF Affinity inside Derived Sequence"
TF_AFFINITY_DIFFERENCE_DERIVED_ANCESTRAL_UNFILTERED = "Unfiltered TF Affinity Difference Derived - Ancestral"
TF_AFFINITY_DIFFERENCE_DERIVED_ANCESTRAL_FILTERED = "Filtered TF Affinity Difference Derived - Ancestral"


def parse_args():
    parser = argparse.ArgumentParser(description="Create TF to loci files.")
    parser.add_argument("--mpra_file", type=str, required=True, help="Path to the MPRA file.")
    parser.add_argument("--pbm_file", type=str, required=True, help="Path to the PBM file.")
    parser.add_argument("--output_dir", type=str, required=True, help="Directory to save the output files.")
    parser.add_argument("--window_size", type=int, default=8, help="Size of seq TF attaches to.")
    parser.add_argument("--start_index", type=int, default=0, help="Index in the MPRA file to start from.")
    parser.add_argument("--end_index", type=int, default=100, help="Index in the MPRA file to end at. For the whole file set to -1.")
    parser.add_argument("--escore_threshold", type=float, required=False, help="E-score threshold to filter on. Default is 0.35.")
    return parser.parse_args()

def check_args(args):
    if not os.path.exists(args.mpra_file):
        raise FileNotFoundError(f"{args.mpra_file} does not exist.")
    if not os.path.exists(args.pbm_file):
        raise FileNotFoundError(f"{args.pbm_file} does not exist.")
    
    if args.window_size < 1:
        raise ValueError("Window size must be greater than 0.")
    if args.start_index < 0:
        raise ValueError("Start index must be greater than or equal to 0.")
    

def save_with_metadata(df, path, metadata: dict):
    with open(path, 'w') as f:
        for key, value in metadata.items():
            f.write(f"# {key}: {value}\n")
        df.to_csv(f)



def create_tf_to_loci_files(MPRA_FILE:str, PBM_FILE:str, OUTPUT_DIR:str, WINDOW_SIZE:int, START_INDEX:int, END_INDEX:int, escore_threshold:float=0.35):
    """
    This function creates files for each locus in the MPRA file. The files contain the z-scores of the TFs for the 8-mers
    in the locus. The files are saved in the OUTPUT_DIR.
    I prefered to show the DFs Transposed to have the TFs as rows and the location on the enhancer as the x-axis.
    :param MPRA_FILE: Path to the MPRA file.
    :param PBM_FILE: Path to the PBM file.
    :param OUTPUT_DIR: Directory to save the output files.
    :param WINDOW_SIZE: Size of the window to consider.
    :param START_INDEX: Index in the MPRA file to start from.
    :param END_INDEX: Index in the MPRA file to end at. For the whole file set to -1.
    :param escore_threshold: E-score threshold to filter on.
    :return None
    """
    os.makedirs(OUTPUT_DIR, exist_ok=True)

    all_8mer_pbm_df = pd.read_csv(PBM_FILE, index_col=0, header=[0,1])
    escore_df, zscore_df = all_8mer_pbm_df['E-score'], all_8mer_pbm_df['Z-score']

    # Read the MPRA file
    columns_needed = ['oligo', 'sequence_ancestral', 'sequence_derived']
    mpra_df = pd.read_csv(MPRA_FILE, usecols=columns_needed)    # We will have a look only at the differencially expressing oligos

    for i, row in mpra_df.iloc[START_INDEX:END_INDEX, :].iterrows():
        clean_locus = um.clean_locus_name(row['oligo'])
        real_locus = row['oligo']

        metadata = {"original_locus": real_locus,
                    "ancestral_sequence": row['sequence_ancestral'],
                    "derived_sequence": row['sequence_derived'],
                    "window_size": WINDOW_SIZE,
                    "escore_threshold": escore_threshold,}
        
        # Save the TF z-scores for the two versions of the locus
        locus_dir = os.path.join(OUTPUT_DIR, clean_locus)
        os.makedirs(locus_dir, exist_ok=True)

        meta_path = os.path.join(locus_dir, "metadata.json")
        with open(meta_path, "w") as f:
            json.dump(metadata, f, indent=4)

        # Get the ancestral and derived DNA sequences
        anc_seq = row['sequence_ancestral'].upper()
        der_seq = row['sequence_derived'].upper()
        
        # Break the sequences into 8-mer windows
        anc_seq_windows = um.get_running_window_numpy(seq=anc_seq, window_size=WINDOW_SIZE)
        der_seq_windows = um.get_running_window_numpy(seq=der_seq, window_size=WINDOW_SIZE)
        
        # Choose for each 8-mer either itself or its reverse complement
        anc_min_windows = um.get_min_seq(anc_seq_windows)
        der_min_windows = um.get_min_seq(der_seq_windows)

        # Get the z-scores for the 8-mers and TFs
        anc_zscore_df_unfiltered = zscore_df.loc[anc_min_windows]
        der_zscore_df_unfiltered = zscore_df.loc[der_min_windows]
        anc_zscore_df_unfiltered.index = anc_seq_windows
        der_zscore_df_unfiltered.index = der_seq_windows


        anc_zscore_df_unfiltered_T = anc_zscore_df_unfiltered.T
        der_zscore_df_unfiltered_T = der_zscore_df_unfiltered.T
        anc_zscore_df_unfiltered_T.to_csv(os.path.join(locus_dir, "Ancestral_PBM_unfiltered.csv" ))
        der_zscore_df_unfiltered_T.to_csv(os.path.join(locus_dir, "Derived_PBM_unfiltered.csv" ))

        anc_zscore_df_unfiltered_T.columns = range(anc_zscore_df_unfiltered_T.shape[1])
        der_zscore_df_unfiltered_T.columns = range(der_zscore_df_unfiltered_T.shape[1])

        diff_pbm_df_unfiltered_T = der_zscore_df_unfiltered_T - anc_zscore_df_unfiltered_T
        diff_pbm_df_unfiltered_T.to_csv(os.path.join(locus_dir, "Difference_PBM_unfiltered.csv" ))

        um.create_heatmap(anc_zscore_df_unfiltered_T, title="TF Affinity in Ancestral Sequence (unfiltered)", locus=clean_locus, path=locus_dir)
        um.create_heatmap(der_zscore_df_unfiltered_T, title="TF Affinity in Derived Sequence (unfiltered)", locus=clean_locus, path=locus_dir)
        um.create_heatmap(diff_pbm_df_unfiltered_T, title="TF Affinity Difference: Derived - Ancestral (unfiltered)", locus=clean_locus, path=locus_dir)

        # Filter on E-score. E-score must be above 0.35 for at least on of the sequences to be considered
        # Filtering starts here:
        anc_escore_df = escore_df.loc[anc_min_windows] > escore_threshold
        der_escore_df = escore_df.loc[der_min_windows] > escore_threshold
        anc_escore_df.index = anc_seq_windows
        der_escore_df.index = der_seq_windows

        anc_zscore_df_filtered_T = anc_zscore_df_unfiltered.where(anc_escore_df, 0).T
        der_zscore_df_filtered_T = der_zscore_df_unfiltered.where(der_escore_df, 0).T
        anc_zscore_df_filtered_T.to_csv(os.path.join(locus_dir, "Ancestral_PBM_filtered.csv" ))
        der_zscore_df_filtered_T.to_csv(os.path.join(locus_dir, "Derived_PBM_filtered.csv" ))

        # ? Might need to first change the colums from 8mers to numbers

        anc_escore_df_T = anc_escore_df.T.set_axis(range(anc_escore_df.shape[0]), axis=1)
        der_escore_df_T = der_escore_df.T.set_axis(range(der_escore_df.shape[0]), axis=1)
        der_anc_diff_mask_T = np.logical_or(anc_escore_df_T, der_escore_df_T)
        diff_pbm_df_filtered_T = diff_pbm_df_unfiltered_T.where(der_anc_diff_mask_T, 0)
        diff_pbm_df_filtered_T.to_csv(os.path.join(locus_dir, "Difference_PBM_filtered.csv" ))

        um.create_heatmap(anc_zscore_df_filtered_T, title="TF Affinity in Ancestral Sequence (Filtered)", locus=clean_locus, path=locus_dir)
        um.create_heatmap(der_zscore_df_filtered_T, title="TF Affinity in Derived Sequence (Filtered)", locus=clean_locus, path=locus_dir)
        um.create_heatmap(diff_pbm_df_filtered_T, title="TF Affinity Difference: Derived - Ancestral (Filtered)", locus=clean_locus, path=locus_dir)

        print(f"Done with {i+1}/{mpra_df.shape[0]}")


def main():
    args = parse_args()
    check_args(args)
    MPRA_FILE = args.mpra_file
    PBM_FILE = args.pbm_file
    OUTPUT_DIR = args.output_dir
    WINDOW_SIZE = args.window_size
    START_INDEX = args.start_index
    END_INDEX = args.end_index if args.end_index != -1 else None
    escore_threshold = args.escore_threshold if args.escore_threshold else 0.35
    create_tf_to_loci_files(MPRA_FILE, PBM_FILE, OUTPUT_DIR, WINDOW_SIZE, START_INDEX, END_INDEX, escore_threshold)


if __name__ == "__main__":
    main()
    # print("Hello World!")