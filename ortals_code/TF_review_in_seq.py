#########################################################################
# This script go over the sequence, for each TF in the library          #
# For each location - save the PBM parameters of the suc-seq binding    #
# Export table of the binding in a csv file                             #
#########################################################################

# Imports:
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
import matplotlib.image as mpimg
from Bio import SeqIO
import pandas as pd
import os
from PIL import Image


# CONSTS:
DATA_TF_PATH = "/home/labs/davidgo/ortalh/TF_analysis_data/data_TF_for_analysis"
PATH_SAVE_RESULTS = "/home/labs/davidgo/ortalh/TF_analysis_results/"
INPUTS_FOLDER = "/home/labs/davidgo/ortalh/TF_analysis_inputs/"
INPUT_CSV = "TF_analysis_inputs.csv"
TF_BINDING_THRESHOLD = 0.35
TF_BINDING_MINIMUM = 0.2
SLASH = "/"
COLORSMAPS = [plt.get_cmap('tab20b_r'), plt.get_cmap('gist_rainbow'), plt.get_cmap('tab20c_r'), plt.get_cmap('tab20_r')]
COLORS = []
K_NARROW = 20
HUMAN_GENOME_FILE = "hg19.fa"
OPP_NUC = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}


def open_txt(file_name):
    with open(file_name, "r") as f:
        txt_data = f.read()
    return txt_data


def get_all_txt_in_folder(folder_path):
    folder_files = os.listdir(folder_path)
    folder_files_txt = []
    for root, dirs, files in os.walk(folder_path):
        for file in files:
            if file.endswith('.txt'):
                full_path = os.path.join(root, file)
                folder_files_txt.append(full_path)
    return folder_files_txt


def filter_TF_txt_files(list_txt_files):

    
    # filtering out all files that are not in the analysis
    # scope - checking the full mer files (not 1..1 or etc. files)
    list_filtered_text_files = []
    uniqe_TF = []
    all_unique_TF = []
    allowed_columns_len = [5, 7, 9]
    for path in list_txt_files:
        path_splited = path.split(SLASH)
        if path_splited[7] not in all_unique_TF:
            all_unique_TF.append(path_splited[7])
        if "mer" in path and ".1" not in path and not "top_enrichment" in path and not "pwm" in path:
            if path_splited[7] not in uniqe_TF:
                if len(open_txt(path).split("\n")[0].split("\t")) in allowed_columns_len:
                    uniqe_TF.append(path_splited[7])
                    list_filtered_text_files.append(path)
    for TF in all_unique_TF:
        if TF not in uniqe_TF:
            print(TF)
    return list_filtered_text_files


def creating_TF_list(TF_data):
    data_split = TF_data.split("\n")
    data_first_row = data_split[0]
    allowed_columns_len = [5,7,9]
    #if len(data_first_row.split("\t")) == 5:
    if len(data_first_row.split("\t")) in allowed_columns_len:
        if len(data_first_row.split("\t")) == 5 or len(data_first_row.split("\t")) == 7:
            if "mer" in data_split[0] and "E-score" in data_split[0] and "Median" in data_split[0] and "Z-score" in data_split[0]:
                data_split = data_split[1:]
            data_TF_org = {}
            for TF_mer in data_split:
                if TF_mer != "":
                    curr_TF_mer = TF_mer.split("\t")
                    data_TF_org[curr_TF_mer[0]] = {"mer-opp": curr_TF_mer[1],
                                                   "E-score": curr_TF_mer[2],
                                                   "Median": curr_TF_mer[3],
                                                   "Z-score": curr_TF_mer[4],
                                                   "mer-len": len(curr_TF_mer[0])}
        else:
            if "mer" in data_split[0] and "E-score" in data_split[0] and "Median" in data_split[0] and "Z-score" in data_split[0]:
                data_split = data_split[1:]
            data_TF_org = {}
            for TF_mer in data_split:
                if TF_mer != "":
                    curr_TF_mer = TF_mer.split("\t")
                    data_TF_org[curr_TF_mer[0]] = {"mer-opp": curr_TF_mer[1],
                                                   "E-score": curr_TF_mer[3],
                                                   "Median": curr_TF_mer[2],
                                                   "Z-score": curr_TF_mer[4],
                                                   "mer-len": len(curr_TF_mer[0])}
        return data_TF_org
    return False


def creating_TF_list_new(TF_data, TF_name, TF_dict_all, TF_dict_all_opp, k_list):
    data_split = TF_data.split("\n")
    data_first_row = data_split[0]
    sum_intensity, count_rows = 0, 0
    allowed_columns_len = [5, 7, 9]
    # if len(data_first_row.split("\t")) == 5:
    if len(data_first_row.split("\t")) in allowed_columns_len:
        if "mer" in data_split[0] and "score" in data_split[0]:
            data_split = data_split[1:]
        if len(data_first_row.split("\t")) == 9:
            escore_col = 3
            median_col = 2
        else:
            escore_col = 2
            median_col = 3

        for TF_mer in data_split:
            if TF_mer != "":
                curr_TF_mer = TF_mer.split("\t")
                if curr_TF_mer[0] not in TF_dict_all:
                    TF_dict_all[curr_TF_mer[0]] = []
                TF_dict_all[curr_TF_mer[0]].append({TF_name: {"mer-opp": curr_TF_mer[1],
                                                              "E-score": curr_TF_mer[escore_col],
                                                              "Median": curr_TF_mer[median_col],
                                                              "Z-score": curr_TF_mer[4],
                                                              "mer-len": len(curr_TF_mer[0])}})
                if curr_TF_mer[1] not in TF_dict_all:
                    if curr_TF_mer[1] not in TF_dict_all_opp:
                        TF_dict_all_opp[curr_TF_mer[1]] = []
                    TF_dict_all_opp[curr_TF_mer[1]].append({TF_name: {"mer": curr_TF_mer[0],
                                                                  "E-score": curr_TF_mer[escore_col],
                                                                  "Median": curr_TF_mer[median_col],
                                                                  "Z-score": curr_TF_mer[4],
                                                                  "mer-len": len(curr_TF_mer[0])}})
                sum_intensity += float(curr_TF_mer[median_col])
                count_rows += 1
                if len(curr_TF_mer[0]) not in k_list:
                    k_list.append(len(curr_TF_mer[0]))
        mean_intensity = sum_intensity / count_rows
        return TF_dict_all, TF_dict_all_opp, k_list, mean_intensity
    print(TF_name, str(len(data_first_row.split("\t"))))
    return False


def check_TF_in_seq_get_all_data_new(seq, all_TF, all_TF_opp, start_seq_index, end_seq_index, mut_location_range, k_list):
    if int(end_seq_index) - int(start_seq_index) + 1 != len(seq):
        print("seq is not correlative to start and end indexes")
    else:
        seq_TF_corralation = {}
        for k_len in k_list:
            for index_zero in range(len(seq)):
                curr_seq_index = index_zero + int(start_seq_index)
                if curr_seq_index not in seq_TF_corralation.keys():
                    seq_TF_corralation[curr_seq_index] = []
                seq_contain_mut = False
                if set(mut_location_range).issubset(set(range(curr_seq_index, curr_seq_index + k_len))):
                    seq_contain_mut = True
                curr_motif = seq[index_zero: index_zero + k_len]
                if curr_motif in all_TF:
                    for TF_dict in all_TF[curr_motif]:
                        TF_name = list(TF_dict.keys())[0]
                        seq_TF_corralation[int(curr_seq_index)].append({'TF Name': TF_name,
                                                                       'Nuc': seq[index_zero],
                                                                       'E-score': float(TF_dict[TF_name]["E-score"]),
                                                                       'Median': float(TF_dict[TF_name]["Median"]),
                                                                       'Z-score': float(TF_dict[TF_name]["Z-score"]),
                                                                       'mer': k_len,
                                                                       'motif': curr_motif,
                                                                       'motif-opp': TF_dict[TF_name]["mer-opp"],
                                                                       'Contain Mut': seq_contain_mut,
                                                                       'Is opp': False})
                if curr_motif in all_TF_opp:
                    for TF_dict in all_TF_opp[curr_motif]:
                        TF_name = list(TF_dict.keys())[0]
                        seq_TF_corralation[int(curr_seq_index)].append({'TF Name': TF_name,
                                                                       'Nuc': seq[index_zero],
                                                                       'E-score': float(TF_dict[TF_name]["E-score"]),
                                                                       'Median': float(TF_dict[TF_name]["Median"]),
                                                                       'Z-score': float(TF_dict[TF_name]["Z-score"]),
                                                                       'mer': k_len,
                                                                       'motif': TF_dict[TF_name]["mer"],
                                                                       'motif-opp': curr_motif,
                                                                       'Contain Mut': seq_contain_mut,
                                                                       'Is opp': True})
        return seq_TF_corralation


def mutation_human(start_index_gene, reg_data_seq, mutation_start, mutation_end, before_mut_nuc, after_mut_nuc):
    start_mut = mutation_start - int(start_index_gene)
    end_mut = mutation_end - int(start_index_gene) + 1
    mut_data_seq = reg_data_seq[:start_mut] + after_mut_nuc + reg_data_seq[end_mut:]
    return mut_data_seq


def check_TF_in_seq_get_all_data(seq, all_TF, start_seq_index, end_seq_index, mut_location_range):
    if int(end_seq_index) - int(start_seq_index) + 1 != len(seq):
        print("seq is not correlative to start and end indexes")
    else:
        seq_TF_corralation = {}
        for index_zero in range(len(seq)):
            curr_seq_index = index_zero + int(start_seq_index)
            if curr_seq_index not in seq_TF_corralation.keys():
                seq_TF_corralation[curr_seq_index] = []
            ## was hereeee seq_contain_mut = False
            for TF_name in all_TF:
                curr_mer_len = int(len(list(all_TF[TF_name].keys())[0]))
                seq_contain_mut = False
                if set(mut_location_range).issubset(set(range(curr_seq_index, curr_seq_index + curr_mer_len))):
                    seq_contain_mut = True
                curr_motif = seq[index_zero: index_zero + curr_mer_len]
                if curr_motif in (all_TF[TF_name]).keys():
                    seq_TF_corralation[int(curr_seq_index)].append({'TF Name': TF_name,
                                                                    'Nuc': seq[index_zero],
                                                                    'E-score': float(all_TF[TF_name][curr_motif]["E-score"]),
                                                                    'Median': float(all_TF[TF_name][curr_motif]["Median"]),
                                                                    'Z-score': float(all_TF[TF_name][curr_motif]["Z-score"]),
                                                                    'mer': curr_mer_len,
                                                                    'motif': curr_motif,
                                                                    'motif-opp': all_TF[TF_name][curr_motif]["mer-opp"],
                                                                    'Contain Mut': seq_contain_mut,
                                                                    'Is opp': False})
                else:
                    for key, inner_dict in all_TF[TF_name].items():
                        if curr_motif == inner_dict["mer-opp"]:
                            seq_TF_corralation[int(curr_seq_index)].append({'TF Name': TF_name,
                                                                            'Nuc': seq[index_zero],
                                                                            'E-score': float(inner_dict["E-score"]),
                                                                            'Median': float(inner_dict["Median"]),
                                                                            'Z-score': float(inner_dict["Z-score"]),
                                                                            'mer': curr_mer_len,
                                                                            'motif': curr_motif,
                                                                            'motif-opp': key,
                                                                            'Contain Mut': seq_contain_mut,
                                                                            'Is opp': True})
        return seq_TF_corralation


def count_TF_contain_mut(TF_seq_corralation):
    # count how much TF we have, and how much is near to the mutation
    TF_list = []
    TF_list_with_mutation = []
    for location in TF_seq_corralation:
        if TF_seq_corralation[location] != []:
            for TF_bind_in_loc in TF_seq_corralation[location]:
                curr_TF_name = TF_bind_in_loc["TF Name"]
                if curr_TF_name not in TF_list:
                    TF_list.append(curr_TF_name)
                in_mutation = TF_bind_in_loc["Contain Mut"]
                escore_above_threshold = TF_bind_in_loc["E-score"]
                if in_mutation and curr_TF_name not in TF_list_with_mutation:
                    TF_list_with_mutation.append(curr_TF_name)
    return TF_list, TF_list_with_mutation


def create_organ_csv(TF_in_organ, species, derived_or_ancestral):
    str_to_save = ""
    enter = "\n"
    for location in TF_in_organ:
        for TF in TF_in_organ[location]:
            str_to_save += (str(location) + "," +
                            str((int(location) + len(TF["motif"]))) + "," +
                            TF["TF Name"] + "," +
                            TF["Nuc"] + "," +
                            TF["motif"] + "," +
                            TF["motif-opp"] + "," +
                            str(TF["Is opp"]) + "," +
                            str(TF["E-score"]) + "," +
                            str(TF["Median"]) + "," +
                            str(TF["Z-score"]) + "," +
                            str(TF["Contain Mut"]) + "," +
                            species + "," +
                            derived_or_ancestral)
            str_to_save += enter
    return str_to_save


def save_data_in_csv(TF_in_human, TF_in_chimp, save_path, derived_species, ancestral_species, opp_mpra):
    enter = "\n"
    str_to_save = "Location Start, Location end, TF Name, Nuc start, motif, motif opp, Is opp, E-score, Intensity, Z-score, location contain mut, Speices, derived/ancestral"
    str_to_save += enter
    if opp_mpra == "FALSE":
        str_to_save += create_organ_csv(TF_in_human, derived_species, "derived")
        str_to_save += create_organ_csv(TF_in_chimp, ancestral_species, "ancestral")
    else:
        str_to_save += create_organ_csv(TF_in_human, derived_species, "ancestral")
        str_to_save += create_organ_csv(TF_in_chimp, ancestral_species, "derived")

    with open("{}all_TF_info.csv".format(save_path), 'w') as file:
        file.write(str(str_to_save))


##########################################################################################################
########################################### Start: Get Genome seq ########################################
def get_genome_name(genome_name):
    genome_name = genome_name.replace('.fasta','')
    genome_name = genome_name.replace('.fa','')
    return genome_name

def get_sequence_records(genome_file,records=0,genome_loaded=0):
    if genome_loaded == 0:
        genomes_directory = '/home/labs/davidgo/Collaboration/Genomes/Genome_fastas/'
        genome_full_path = genomes_directory + genome_file
        if os.path.exists(genome_full_path):
            records = list(SeqIO.parse(genome_full_path, "fasta"))
            print(f'Imported genome {genome_full_path}')
        else:
            print(f'No Genome file {genome_full_path}')
    return records

def seq_retrieve_interval(chromosome, start, end, genome_file, records = 0, genome_loaded = 0):
    """
    This function extracts the desired sequence from sequence records by chromosome and position for ONE interval
    records: sequence records
    """
    records = get_sequence_records(genome_file,records, genome_loaded)
    curr_chr = [record for record in records if record.id == chromosome] # this is how to exctract the record that matches our chromosome
    sequence = curr_chr[0].seq[start - 1:end] # note - the coords here are similar to bed format: start is -1 (starting with 0), and end is not included
    return sequence

def seq_retrieve_df(df, chr_column, start_column, end_column, genome_file):
    # This function gets a df and column names for chr start end and a sequence records objects and returns the sequence for each df line as a new column
    genome_name = get_genome_name(genome_file)
    records = get_sequence_records(genome_file,0,0)
    df[genome_name] = df.apply(func=lambda x: seq_retrieve_interval(x[chr_column], x[start_column], x[end_column], genome_file, records, 1), axis=1).apply(''.join)
    return df
########################################### End Get Genome seq ###########################################
##########################################################################################################


def get_parameters_from_input(input_file):
    with open(input_file, 'r') as f:
        data = f.read()
    data = data.split("\n")
    gene_name = data[0]
    human_chr = (data[1]).split(" ")[0]
    human_seq_start = int((data[1]).split(" ")[1])
    human_seq_end = int((data[1]).split(" ")[2])
    mutation_loc_start = int((data[2]).split(" ")[0])
    mutation_loc_end = int((data[2]).split(" ")[1])
    nuc_human = (data[2]).split(" ")[2]
    nuc_ancestral = (data[2]).split(" ")[3]
    return gene_name, human_chr, human_seq_start, human_seq_end, range(mutation_loc_start, mutation_loc_end), nuc_human, nuc_ancestral


def get_parameters_from_input_csv(input_row):
    data = input_row.split(",")
    gene_name = data[0]
    human_chr = data[1]
    human_seq_start = int(data[2])
    human_seq_end = int(data[3])
    mutation_loc_start = int(data[4])
    mutation_loc_end = int(data[5])
    nuc_human = data[6]
    nuc_variant = data[7]
    derived_species = data[8]
    ancestral_species = data[9]
    opp_mpra = data[10]
    return gene_name, human_chr, human_seq_start, human_seq_end, range(mutation_loc_start, mutation_loc_end), nuc_human, nuc_variant, derived_species, ancestral_species, opp_mpra


def get_mean_intensity_TF_chip_PBM(TF_dict_all):
    # get mean intensity of a PBM chip for every TF in the analysis
    TFs_mean_intensity = {}
    for TF_name in TF_dict_all:
        sum_intensity = []
        for seq_mer in TF_dict_all[TF_name]:
            sum_intensity.append(float(TF_dict_all[TF_name][seq_mer]["Median"]))
        mean_intensity = np.mean(sum_intensity)
        TFs_mean_intensity[TF_name] = mean_intensity
    return TFs_mean_intensity


def get_mean_intensity_TF_chip_PBM_new(TF_dict_all, TF_dict_all_opp):
    # get mean intensity of a PBM chip for every TF in the analysis
    TFs_mean_intensity = {}
    for TF_name in TF_dict_all:
        sum_intensity = []
        for seq_mer in TF_dict_all[TF_name]:
            sum_intensity.append(float(TF_dict_all[TF_name][seq_mer]["Median"]))
        mean_intensity = np.mean(sum_intensity)
        TFs_mean_intensity[TF_name] = mean_intensity
    return TFs_mean_intensity


def save_mean_intensity_for_TF(TFs_mean_intensity, path_save):
    str_csv = "TF, mean_intensity\n"
    for TF in TFs_mean_intensity:
        str_csv += TF + ',' + str(TFs_mean_intensity[TF]) + '\n'
    with open("{}mean_intensity.csv".format(path_save), 'w') as f:
        f.write(str_csv)


def save_to_candidate_TF(TF_in_human, TF_in_variant, save_path, derived_species, ancestral_species, opp_mpra):
    enter = "\n"
    str_to_save = "TF Name, Location Start, Location End, mer, is_opp, diff escore (derived-ancestral), diff intensity (derived-ancestral), diff zcore (derived-ancestral), escore human, escore variant, intensity human, intensity variant, zscore human, zscore variant"
    str_to_save += enter
    opp_param = 1
    if opp_mpra == "TRUE":
        opp_param = -1
    for location in TF_in_human:
        for TF in TF_in_human[location]:
            if TF["Contain Mut"]:
                ancestral_dict = next((d for d in TF_in_variant[location] if d.get('TF Name') == TF["TF Name"]), None)
                if float(TF["E-score"]) > TF_BINDING_THRESHOLD or float(ancestral_dict["E-score"]) > TF_BINDING_THRESHOLD:
                    if TF["Is opp"]:
                        mer = TF["motif-opp"]
                    else:
                        mer = TF["motif"]
                    str_to_save += (TF["TF Name"] + "," +
                                    str(location) + "," +
                                    str(location + TF['mer']) + "," +
                                    str(mer) + "," +
                                    str(TF["Is opp"]) + "," +
                                    str(float(float(TF["E-score"]) - float(ancestral_dict["E-score"])) * opp_param) + "," +
                                    str(float(float(TF["Median"]) - float(ancestral_dict["Median"])) * opp_param) + "," +
                                    str(float(float(TF["Z-score"]) - float(ancestral_dict["Z-score"]))* opp_param) + "," +
                                    str(TF["E-score"]) + "," +
                                    str(ancestral_dict["E-score"]) + "," +
                                    str(TF["Median"]) + "," +
                                    str(ancestral_dict["Median"]) + "," +
                                    str(TF["Z-score"]) + "," +
                                    str(ancestral_dict["Z-score"]))
                    str_to_save += enter

    with open("{}TF_near_mut_data.csv".format(save_path), 'w') as f:
        f.write(str_to_save)

def import_gals_data(tsv_row, counter_gene):
    gene_name = 'Seq{}'.format(counter_gene)
    human_chr = tsv_row[1]
    human_seq_start = int(tsv_row[4])
    human_seq_end = int(tsv_row[5])
    mutation_loc_start = int(tsv_row[2]) + 1
    mutation_loc_end = int(tsv_row[3]) + 1
    nuc_human = tsv_row[6]
    nuc_ancestral = tsv_row[7]

    return gene_name, human_chr, human_seq_start, human_seq_end, range(mutation_loc_start,
                                                                       mutation_loc_end), nuc_human, nuc_ancestral


def save_all_TFs_in_analysis(list_TF_path_filtered, TF_list_filtered_out, path_save):
    csv_path_str = "TF name, PBM result used\n"
    for curr_TF_path in list_TF_path_filtered:
        if curr_TF_path not in TF_list_filtered_out:
            csv_path_str += curr_TF_path.split(SLASH)[7] + "," + curr_TF_path + "\n"
    with open("{}TFs_pbm_results_paths.csv".format(path_save), 'w') as f:
        f.write(csv_path_str)


def main():
    """
    print("### getting Gal's data ###")
    gals_data = np.loadtxt("{}gal_sample_variants_df.tsv".format(INPUTS_FOLDER), delimiter='\t', dtype=str)
    """

    """
    print("##### Getting all variants input paths #####")
    input_files = []
    for path in os.listdir(INPUTS_FOLDER):
        input_files.append(path)
    """

    print("##### Getting all variants input paths #####")
    input_file_path = "{}{}".format(INPUTS_FOLDER, INPUT_CSV)
    data_input_file = open_txt(input_file_path)
    input_files = data_input_file.split("\n")[1:]

    print("##### Extracting txt files from TF path #####")
    print("##### Reading TF files #####")
    list_TF_path_files = get_all_txt_in_folder(DATA_TF_PATH)
    list_TF_path_filtered = filter_TF_txt_files(list_TF_path_files)

    TF_dict_all = {}
    count_all_TF = 0
    TF_dict_all_opp = {}
    mean_intensity_for_TF = {}
    TF_list_filtered_out = []
    k_list = []
    for curr_TF_path in list_TF_path_filtered:
        curr_TF_data = open_txt(curr_TF_path)
        TF_name = curr_TF_path.split(SLASH)[7]
        result = creating_TF_list_new(curr_TF_data, TF_name, TF_dict_all, TF_dict_all_opp, k_list)
        if result:
            TF_dict_all, TF_dict_all_opp, k_list, mean_intensity = result[0], result[1], result[2], result[3]
            mean_intensity_for_TF[TF_name] = mean_intensity
            count_all_TF += 1
        else:
            print("####### filtered out! file not in the right format! : {}".format(curr_TF_path))
            TF_list_filtered_out.append(curr_TF_path)
    print("####### overall not in the right format: {}".format(len(TF_list_filtered_out)))
    print("####### overall in the right format: {}".format(count_all_TF))

    print()
    counter = 0
    #for curr_variant_input_file in input_files:
    #for curr_variant_input_file in gals_data[1:]:
    for curr_variant_input_row in input_files:
        if curr_variant_input_row != "":
            #input_file = "{}{}".format(INPUTS_FOLDER, curr_variant_input_file)
            #gene_name, human_chr, human_seq_start, human_seq_end, mutation_loc_range, nuc_human, nuc_ancestral = get_parameters_from_input(input_file)

            gene_name, human_chr, human_seq_start, human_seq_end, mutation_loc_range, nuc_human, nuc_ancestral, derived_species, ancestral_species, opp_mpra = get_parameters_from_input_csv(curr_variant_input_row)

            ## Gal's
            # counter += 1
            # gene_name, human_chr, human_seq_start, human_seq_end, mutation_loc_range, nuc_human, nuc_ancestral = import_gals_data(curr_variant_input_file, counter)
            # PATH_SAVE_RESULTS = "/home/labs/davidgo/ortalh/TF_analysis_results/Gal_sample_variants/"
            ## End Gal's
    
            print("##### Getting human seq from human genome #####")
            print("#### Seq: {}, {}: {}-{}, mutation: {}-{} {}->{} ####".format(HUMAN_GENOME_FILE.split(".")[0], human_chr,
                                                                             human_seq_start, human_seq_end, mutation_loc_range[0], mutation_loc_range[-1],
                                                                             nuc_human, nuc_ancestral))
            human_seq = (seq_retrieve_interval(human_chr, human_seq_start, human_seq_end, HUMAN_GENOME_FILE)).upper()
            correct_mutation = True
            if len(mutation_loc_range) == 1:
                if human_seq[mutation_loc_range[0] - human_seq_start] != nuc_human:
                    print("Please check your mutation location!")
                    correct_mutation = False
            else:
                if human_seq[mutation_loc_range[0] - human_seq_start: mutation_loc_range[-1] - human_seq_start + 1] != nuc_human:
                    print("Please check your mutation location!")
                    correct_mutation = False
            if correct_mutation:
                print("### Starting simulation ###")
    
                print("##### Checking Which TFs bind the seq")
                print("####### Creating Mutated seq")
                mut_data_seq = mutation_human(human_seq_start, human_seq, mutation_loc_range[0], mutation_loc_range[-1],
                                              nuc_human, nuc_ancestral)
    
                print("### Checking mutation creation ###")
                print("## human seq len: {}, mut seq len: {}".format(len(human_seq), len(mut_data_seq)))
                print("## human seq near mut: {}".format(str(human_seq[mutation_loc_range[0] - human_seq_start - 6: mutation_loc_range[0] - human_seq_start + 6])))
                print("## mut   seq near mut: {}".format(str(mut_data_seq[mutation_loc_range[0] - human_seq_start - 6: mutation_loc_range[0] - human_seq_start + 6])))
    
                print("####### Getting data for the human seq")
                # human_all_motif_bind = check_TF_in_seq_get_all_data(human_seq, TF_dict_all, human_seq_start, human_seq_end, mutation_loc_range)
                human_all_motif_bind = check_TF_in_seq_get_all_data_new(human_seq, TF_dict_all, TF_dict_all_opp, human_seq_start, human_seq_end, mutation_loc_range, k_list)
                human_TF_list, human_TF_list_in_mutation = count_TF_contain_mut(human_all_motif_bind)
                print(
                    "##### Overall TFs binding to human seq: {} out of {} TFs in data. {} of them has a binding site in the mutation area.".format(
                        len(human_TF_list), count_all_TF, len(human_TF_list_in_mutation)))
    
                print("####### Getting data for the variant seq")
                # chimp_all_motif_bind = check_TF_in_seq_get_all_data(mut_data_seq, TF_dict_all, human_seq_start, human_seq_end, mutation_loc_range)
                variant_all_motif_bind = check_TF_in_seq_get_all_data_new(mut_data_seq, TF_dict_all, TF_dict_all_opp,
                                                                        human_seq_start, human_seq_end, mutation_loc_range,
                                                                        k_list)
                chimp_TF_list, chimp_TF_list_in_mutation = count_TF_contain_mut(variant_all_motif_bind)
                print(
                    "##### Overall TFs binding to variant seq: {} out of {} TFs in data. {} of them has a binding site in the mutation area.".format(
                        len(chimp_TF_list), count_all_TF, len(chimp_TF_list_in_mutation)))
    
                seq_folder = "{}{}_{}_{}_{}-{}_mutation_{}-{}_{}_{}".format(PATH_SAVE_RESULTS, HUMAN_GENOME_FILE.split(".")[0], gene_name,
                                                                      human_chr, human_seq_start, human_seq_end, mutation_loc_range[0], mutation_loc_range[-1],
                                                                      nuc_human, nuc_ancestral)
    
                try:
                    os.mkdir(seq_folder)
                    os.mkdir(seq_folder + SLASH + "summarized_plots")
                    os.mkdir(seq_folder + SLASH + "TFs_binding_csv")
                    os.mkdir(seq_folder + SLASH + "TFs_binding_plots")
                    print("###### creating folder {}".format(seq_folder))
                except FileExistsError:
                    print("###### folder already exists! {}".format(seq_folder))
    
                print("##### Saving Binding Data for each TF")
                csv_save_path = "{}{}TFs_binding_csv".format(seq_folder, SLASH)
                try:
                    os.mkdir(csv_save_path)
                    print(f"## Saving TFs binding csv here: '{csv_save_path}' (directory created)")
                except FileExistsError:
                    print(f"## Saving TFs binding csv here: '{csv_save_path}' (directory already exists)")
                save_data_in_csv(human_all_motif_bind, variant_all_motif_bind, csv_save_path + SLASH, derived_species, ancestral_species, opp_mpra)
                save_to_candidate_TF(human_all_motif_bind, variant_all_motif_bind, csv_save_path + SLASH, derived_species, ancestral_species, opp_mpra)
                save_all_TFs_in_analysis(list_TF_path_filtered, TF_list_filtered_out, csv_save_path + SLASH)
                print("## CSV saved! {}".format(csv_save_path))
    
                print("####### getting mean for each PBM chip")
                # mean_intensity_for_TF = get_mean_intensity_TF_chip_PBM(TF_dict_all)
                save_mean_intensity_for_TF(mean_intensity_for_TF, csv_save_path + SLASH)
                print("## CSV saved! {}".format(csv_save_path))
                print()


if __name__ == "__main__":
    main()
