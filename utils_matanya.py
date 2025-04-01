import os
import re
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.colors import CenteredNorm
import seaborn as sns
from ortals_code.TF_review_in_seq import *
import ortals_code.TF_review_in_seq as ortals
from functools import reduce


def parse_8mer_tf_file(file_name):
    """
    This function takes a file with the PBM results on all 8mer sequences and a specific tf. It returns a dataframe
    with the 8mer, the reverse complement, the E-score, the median and the Z-score.
    :param file_name: str
    :return: pd.DataFrame if it's an 8-mer file, and the name of the TF. None otherwise.
    """
    # print("This got modified")
    df = pd.read_csv(file_name, sep="\t", header=None)
    if "AAA" not in df.iloc[0, 0]:
        df = df.iloc[1:].reset_index(drop=True)
    if df.shape[0] != 32896:  
        return None, None
    df = df.iloc[:, :5]
    TF_name = file_name.split(SLASH)[7]
    df.columns = ["8-mer", "reverse_complement", f"E-score_{TF_name}", f"Median_{TF_name}", f"Z-score_{TF_name}"]
    
    return df, TF_name

def get_reverse_complement(seq):
    """
    This function takes a sequence and returns its reverse complement.
    :param seq: str
    :return: str
    """
    trans_table = str.maketrans("ATGC", "TACG")
    if isinstance(seq, pd.Series):
        return seq.str.translate(trans_table).str[::-1]

    elif isinstance(seq, str):
        return seq.translate(trans_table)[::-1]
    else: 
        raise TypeError("Input must be a string or a Pandas Series of strings")
    

def get_running_window_numpy(seq, window_size=8):
    if len(seq) < window_size:
        return pd.Series([])
    dna_array = np.array(list(seq))
    from numpy.lib.stride_tricks import sliding_window_view
    windows = sliding_window_view(dna_array, window_size)
    return pd.Series(["".join(window) for window in windows])


def get_min_seq(seq):
    """
    This function takes a sequence or a Pandas Series of sequences and returns 
    the lexicographically smaller sequence between itself and its reverse complement.
    
    Used to get the correct sequence that appears in the DataFrame index.
    
    :param seq: str or pd.Series of str
    :return: str or pd.Series of str
    """
    reverse_complement_seq = get_reverse_complement(seq)

    if isinstance(seq, pd.Series):
        return pd.concat([seq, reverse_complement_seq], axis=1).min(axis=1)
    
    elif isinstance(seq, str):
        return min(seq, reverse_complement_seq)

    else:
        raise TypeError("Input must be a string or a Pandas Series of strings")


def convert_df_into_multiindex(df, id_cols=["8-mer", "reverse_complement"]):
    """
    Convert the 8-mer dataframe into a multiindex dataframe.
    """
    multi_columns = []
    for col in df.columns:
        if col in id_cols:
            multi_columns.append(("ID", col))
        else:
            parts = col.split("_", 1)
            if len(parts) == 2:
                category, tf_name = parts
            else:
                category, tf_name = col, "unkown"
            multi_columns.append((category, tf_name))
    new_df = df.copy()
    new_df.columns = pd.MultiIndex.from_tuples(multi_columns)
    return new_df


def get_running_window_numpy(seq, window_size=8):
    dna_array = np.array(list(seq))
    from numpy.lib.stride_tricks import sliding_window_view
    windows = sliding_window_view(dna_array, window_size)
    return pd.Series(["".join(window) for window in windows])

def locus_to_file_name(locus):
    pass


def create_heatmap(tf_2_seq,
                   title="Transcription Factor Affinity Heatmap", 
                   x_label="Starting Index", 
                   y_label="Transcription Factors",
                   figsize=(15, 15), 
                   path=None, 
                   locus=''):
    
    fig, ax = plt.subplots(figsize=figsize)
    norm = CenteredNorm(vcenter=0)  

    sns.heatmap(tf_2_seq, cmap='RdBu_r', annot=False, ax=ax, norm=norm)

    # Adjust layout to prevent overlapping
    fig.subplots_adjust(top=0.85)  # Moves everything down slightly
    
    # If locus is provided, add it as a secondary title above the main title
    if locus:
        fig.suptitle(f"Locus: {locus}", fontsize=14, fontweight='bold', y=0.98)  # Moves it to the very top
    
    # Main title below the secondary title
    ax.set_title(title, fontsize=16, pad=20)

    # Axis labels
    ax.set_xlabel(x_label, fontsize=12)
    ax.set_ylabel(y_label, fontsize=12)

    # Ensure path exists if provided
    if path:
        if not os.path.exists(path):
            os.makedirs(path)
        save_path = os.path.join(path, re.sub(r'[^a-zA-Z0-9_-]', '_', title) + ".png")
        plt.savefig(save_path)
        # print(f"Heatmap saved to: {save_path}")

    # plt.show()  # Display the heatmap
    plt.close()




# Computer Science Nerd Alert! 
# The following two functions are used to convert a sequence into it's index. 
# This is done by converting the sequence into a base-4 number, and then converting that number into a base-10 number.

# UPDATE: This code is not helpful as the sequences are not exactly in order because of the palindromic nature of the sequences.
# I will keep it here for now, but I will not use it.
# ! cool idea, redundent code.

def into_base_4(seq):
    """
    This function takes a sequence and returns its base 4 representation.
    :param seq: str
    :return: str
    """
    return int("".join([{"A": "0", "C": "1", "G": "2", "T": "3"}[base] for base in seq]))

def base4_to_base10(base_4_int):
    """
    This function takes a base 4 integer and returns its base 10 representation.
    :param base_4_int: int
    :return: int
    """
    return int(str(base_4_int), 4)

def seq_to_int(seq):
    """
    This function takes a sequence and returns its index.
    :param seq: str
    :return: int
    """
    return base4_to_base10(into_base_4(seq))

def int_to_seq(index, seq_length=8):
    """
    This function takes an index and returns the corresponding sequence.
    :param index: int
    :return: str
    """
    seq =  "".join([{"0": "A", "1": "C", "2": "G", "3": "T"}[base] for base in np.base_repr(index, 4)])
    # the rjust is to ensure we get a sequence of the right length, in base-4 we treat A as a 0. 
    return seq.rjust(seq_length, "A")


    