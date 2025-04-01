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
from scipy import stats
import scipy.stats as stats
from statsmodels.stats.multitest import multipletests
import matplotlib.patches as patches


# CONSTS:
DATA_TF_PATH = "/home/labs/davidgo/ortalh/TF_analysis_data/data_TF_for_analysis"
PATH_SAVE_RESULTS = "/home/labs/davidgo/ortalh/TF_analysis_results/"
INPUTS_FOLDER = "/home/labs/davidgo/ortalh/TF_analysis_inputs/"
COLOR_CONFIG_FILE = "/home/labs/davidgo/Collaboration/backup/EvoRegions/General/colors_config.txt"
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


def get_colors():
    for cmap in COLORSMAPS:
        num_colors_needed = 200 // len(COLORSMAPS)
        for i in range(num_colors_needed - 5):
            color = cmap((i + 5) / num_colors_needed)
            COLORS.append(color)


def select_color_from_gradient(value, gradient):
    """Select a color from a gradient array based on a normalized value between 0 and 1.
    The number of bins matches the number of RGBs in the gradient."""
    if pd.isna(value):
        return (0, 0, 0)  # Return black for missing data
    index = int(value * (len(gradient) - 1))
    index = max(0, min(index, len(gradient) - 1))  # Ensure index is within the gradient's range
    return gradient[index]


def get_diff_mut_value(locations_human, locations_mut, human, mut):
    # store only the unique values of the mutations, if the value appears in human seq not containing it
    diff_mut_locations = []
    diff_mut_value = []
    for loc in locations_mut:
        index_mut = locations_mut.index(loc)
        if loc in locations_human:
            index_human = locations_human.index(loc)
            if human[index_human] != mut[index_mut]:
                diff_mut_locations.append(loc)
                diff_mut_value.append(mut[index_mut])
        else:
            print("unique location in mutated seq")
            diff_mut_locations.append(loc)
            diff_mut_value.append(mut[index_mut])
    return diff_mut_locations, diff_mut_value


def plot_a_graph_for_TF(sub_rows, sub_cols, sub_curr, plot_id, feature_name, TF_name, feature_list, keys_list, color_plot, color_loc_mut, mut_location_range, seq_id_txt, y_label):
    plt.subplot(sub_rows, sub_cols, sub_curr)
    plt.title("{} of TF: {} - {}".format(feature_name, TF_name, plot_id))
    plt.plot(keys_list, feature_list, 'o', color=color_plot)
    for mut_location in mut_location_range:
        plt.axvline(x=mut_location, linestyle='--', color=color_loc_mut)
    plt.xlabel(seq_id_txt)
    plt.ylabel(y_label)


def plot_graph_for_TF_diff_mut(sub_rows, sub_cols, sub_curr, keys_mut_diff, feature_reg, feature_list_diff, color_plot, color_loc_mut, color_diff, mut_location_range, keys_reg, seq_id_txt, nuc_list_mut, minus_seq_binding, keys_int_mut, TF_significant_locations, significant_color):
    plt.subplot(sub_rows, sub_cols, sub_curr)
    plt.plot(keys_mut_diff, feature_list_diff, 'o', color=color_plot)
    for mut_location in mut_location_range:
        plt.axvline(x=mut_location, linestyle='--', color=color_loc_mut)
    for i in range(len(keys_mut_diff)):
        if keys_mut_diff[i] in keys_reg:
            human_index = keys_reg.index(keys_mut_diff[i])
            plt.plot([keys_mut_diff[i], keys_mut_diff[i]],
                     [feature_list_diff[i], feature_reg[human_index]],
                     color=color_diff)
    plt.xticks([])
    plt.subplot(sub_rows, sub_cols, sub_curr + sub_cols)
    plt.plot(keys_mut_diff, feature_list_diff, 'o', color=color_plot)
    for mut_location in mut_location_range:
        plt.axvline(x=mut_location, linestyle='--', color=color_loc_mut)
    for i in range(len(keys_mut_diff)):
        if keys_mut_diff[i] in keys_reg:
            human_index = keys_reg.index(keys_mut_diff[i])
            if keys_mut_diff[i] in TF_significant_locations:
                plt.plot([keys_mut_diff[i], keys_mut_diff[i]],
                         [feature_list_diff[i], feature_reg[human_index]],
                         color=significant_color)
            else:
                plt.plot([keys_mut_diff[i], keys_mut_diff[i]],
                         [feature_list_diff[i], feature_reg[human_index]],
                         color=color_diff)

    plt.xticks(keys_int_mut, nuc_list_mut)
    labels = plt.gca().get_xticklabels()
    for i, label in enumerate(labels):
        if i in minus_seq_binding:
            label.set_color('purple')
    plt.xlim(mut_location_range[0] - 20, mut_location_range[0] + 20)


def plot_TF_escore_and_median_in_seq_nice(human_TF_list, mut_TF_list, mut_location, mean_intensity_for_TF, seq_id_txt, save_name, TF_significant_dict, derived_scpecies, ancestral_species):
    mut_color = 'palevioletred'
    reg_color = 'steelblue'
    color_diff = 'paleturquoise'
    color_diff_significant = 'sandybrown'
    mut_loc_color = 'silver'

    unique_TFs_human = set(TF['TF Name'] for location in human_TF_list.values() for TF in location)
    unique_TFs_mut = set(TF['TF Name'] for location in mut_TF_list.values() for TF in location)
    shared_mut = []

    for index, TF in enumerate(unique_TFs_human):
        filtered_data = {location: TF_dict for location, list_TF_in_location in human_TF_list.items() for TF_dict in list_TF_in_location if TF_dict['TF Name'] == TF}
        keys = list(filtered_data.keys())
        keys_int_human = [int(k) for k in keys]
        escores_human = [filtered_data[k]['E-score'] for k in keys]
        intensities_human = [(np.log(filtered_data[k]['Median']) - np.log(mean_intensity_for_TF[TF])) for k in keys]
        zscores_human = [filtered_data[k]['Z-score'] for k in keys]
        plt.figure(figsize=(20, 10))
        plt.suptitle("Escore, Median Intensity and Z-score results for {}".format(TF), fontsize=14)
        # plot Escore - all
        plot_a_graph_for_TF(2, 3, 1, "all seq", "E-score", TF, escores_human, keys_int_human,
                            reg_color, mut_loc_color, mut_location, seq_id_txt, "E-score")
        plt.xticks([])
        # plot Intensity - all
        plot_a_graph_for_TF(2, 3, 2, "all seq", "Median Intensity", TF, intensities_human, keys_int_human,
                            reg_color, mut_loc_color, mut_location, seq_id_txt, "ln[Intensity] - ln[chip_average]")
        plt.xticks([])
        # plot Z-score - all
        plot_a_graph_for_TF(2, 3, 3, "all seq", "Z-score", TF, zscores_human, keys_int_human,
                            reg_color, mut_loc_color, mut_location, seq_id_txt, "Z-score")
        plt.xticks([])
        # plot Escore - +-30
        nuc_list = [filtered_data[k]['Nuc'] for k in keys]
        plot_a_graph_for_TF(2, 3, 4, "mut location \xB120", "E-score", TF, escores_human, keys_int_human,
                            reg_color, mut_loc_color, mut_location, seq_id_txt, "E-score")
        plt.xticks(keys_int_human, nuc_list)
        plt.xlim(mut_location[0] - 20, mut_location[0] + 20)
        # plot Intensity - +- 30
        plot_a_graph_for_TF(2, 3, 5, "mut location \xB120", "Median Intensity", TF, intensities_human, keys_int_human,
                            reg_color, mut_loc_color, mut_location, seq_id_txt, "ln[Intensity] - ln[chip_average]")
        plt.xticks(keys_int_human, nuc_list)
        plt.xlim(mut_location[0] - 20, mut_location[0] + 20)
        # plot Z-score - +- 30
        plot_a_graph_for_TF(2, 3, 6, "mut location \xB120", "Z-score", TF, zscores_human, keys_int_human,
                            reg_color, mut_loc_color, mut_location, seq_id_txt, "Z-score")
        plt.xticks(keys_int_human, nuc_list)
        plt.xlim(mut_location[0] - 20, mut_location[0] + 20)

        # plot mut diff
        if TF in unique_TFs_mut:
            shared_mut.append(TF)
            filtered_data_mut = {location: TF_dict for location, list_TF_in_location in mut_TF_list.items() for TF_dict in
                             list_TF_in_location if TF_dict['TF Name'] == TF}
            keys_mut = list(filtered_data_mut.keys())
            escores_mut = [filtered_data_mut[k]['E-score'] for k in keys]
            intensities_mut = [(np.log(filtered_data_mut[k]['Median']) - np.log(mean_intensity_for_TF[TF])) for k in keys]
            zscores_mut = [filtered_data_mut[k]['Z-score'] for k in keys]
            keys_int_mut = [int(k) for k in keys_mut]
            keys_mut_diff_escore, escores_mut_diff = get_diff_mut_value(keys_int_human, keys_int_mut, escores_human, escores_mut)
            keys_mut_diff_intensity, intensity_mut_diff = get_diff_mut_value(keys_int_human, keys_int_mut, intensities_human, intensities_mut)
            keys_mut_diff_zscore, zscores_mut_diff = get_diff_mut_value(keys_int_human, keys_int_mut,
                                                                             zscores_human, zscores_mut)
            nuc_list_mut = [filtered_data_mut[k]['Nuc'] for k in keys]
            minus_seq_binding = [i for i, k in enumerate(keys) if filtered_data_mut[k]['Is opp']]
            if TF not in TF_significant_dict:
                significant_list = []
            else:
                significant_list = TF_significant_dict[TF]
            if keys_mut_diff_escore != []:
                plot_graph_for_TF_diff_mut(2, 3, 1, keys_mut_diff_escore, escores_human,
                                           escores_mut_diff, mut_color, mut_loc_color, color_diff, mut_location, keys_int_human,
                                           seq_id_txt, nuc_list_mut, minus_seq_binding, keys_int_mut, significant_list, color_diff_significant)
            if keys_mut_diff_intensity != []:
                plot_graph_for_TF_diff_mut(2, 3, 2, keys_mut_diff_intensity, intensities_human,
                                           intensity_mut_diff, mut_color, mut_loc_color, color_diff, mut_location, keys_int_human,
                                           seq_id_txt, nuc_list_mut, minus_seq_binding, keys_int_mut, significant_list, color_diff_significant)
            if keys_mut_diff_zscore != []:
                plot_graph_for_TF_diff_mut(2, 3, 3, keys_mut_diff_zscore, zscores_human,
                                           zscores_mut_diff, mut_color, mut_loc_color, color_diff, mut_location, keys_int_human,
                                           seq_id_txt, nuc_list_mut, minus_seq_binding, keys_int_mut, significant_list, color_diff_significant)

            reg_legend = mlines.Line2D([], [], marker='o', color=reg_color, markersize=10, label=derived_scpecies)
            mut_legend = mlines.Line2D([], [], marker='o', color=mut_color, markersize=10, label=ancestral_species)
            mut_loc_legend = mlines.Line2D([], [], linestyle='--', color=mut_loc_color, markersize=10, label='Mut Loc')

            #plt.legend(handles=[reg_legend, mut_legend, mut_loc_legend], bbox_to_anchor=(-0.1, -0.5), fontsize=10, loc='center', ncol=3)
            #plt.legend(handles=[reg_legend, mut_legend, mut_loc_legend], bbox_to_anchor=(0.1, -0.2), fontsize=10, loc='upper center', ncol=3)
            plt.legend(handles=[reg_legend, mut_legend, mut_loc_legend], bbox_to_anchor=(1, 1), fontsize=10,
                       loc='upper left')
            plt.tight_layout()
            save_path = os.path.join("{}{}.jpg".format(save_name, TF))
            plt.savefig(save_path, bbox_inches='tight')
            plt.close()


def concatenate_TF_plot_pwn(pwm_png_dict, TF_folder_path_jpg):
    pwm_with_plot_path = "{}plots_with_pwm{}".format(TF_folder_path_jpg, SLASH)
    try:
        os.mkdir(pwm_with_plot_path)
        print("##### plots with PWM in the folder: {} (folder created)".format(pwm_with_plot_path))
    except FileExistsError:
        print("##### plots with PWM in the folder: {} (folder exists)".format(pwm_with_plot_path))
    TF_jpg = {}
    for item in os.listdir(TF_folder_path_jpg):
        full_path = os.path.join(TF_folder_path_jpg, item)
        if item.endswith('.jpg'):
            TF_name = (full_path.split('/')[8]).split('.')[0]
            if TF_name not in TF_jpg:
                TF_jpg[TF_name] = Image.open(full_path)
    for TF_name in TF_jpg.keys():
        TF_jpg_img = TF_jpg[TF_name]
        if TF_name in pwm_png_dict:
            TF_pwm_png_img = pwm_png_dict[TF_name]
            total_width = max(TF_jpg_img.width, TF_pwm_png_img.width)
            total_height = TF_jpg_img.height + TF_pwm_png_img.height
            new_image = Image.new('RGBA', (total_width, total_height),
                                  (255, 255, 255, 0))  # RGBA for transparency, or RGB for no transparency

            new_image.paste(TF_jpg_img, (0, 0))
            x_offset = (total_width - TF_pwm_png_img.width) // 2
            new_image.paste(TF_pwm_png_img, (x_offset, TF_jpg_img.height))
            rgb_image = new_image.convert('RGB')
            rgb_image.save('{}{}.jpg'.format(pwm_with_plot_path, TF_name), 'JPEG')
    print("##### concatenate figures finished!")

def create_organ_csv(TF_in_organ, species):
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
                           str(TF["Contain Mut"]) + "," +
                           species)
            str_to_save += enter
    return str_to_save


def summarized_figure_for_seq(human_all_motif_bind, chimp_all_motif_bind, mut_location, save_path, seq_info):
    dict_TF_escore = {}

    for location in human_all_motif_bind:
        for TF in human_all_motif_bind[location]:
            curr_Escore = float(TF["E-score"])
            tf_name = TF["TF Name"]
            if tf_name not in dict_TF_escore.keys():
                dict_TF_escore[tf_name] = TF | {"location": location}
            else:
                max_dict = dict_TF_escore[tf_name]
                if curr_Escore > float(max_dict["E-score"]):
                    dict_TF_escore[tf_name] = TF | {"location": location}

    plt.figure(figsize=(15, 10))

    sorted_dict_TF_escore = dict(sorted(dict_TF_escore.items(), key=lambda x: x[1]['E-score'], reverse=True))

    TF_description = [f"{tf_data['TF Name']}, {tf_data['motif-opp']}, {tf_data['location']}" for tf_data in sorted_dict_TF_escore.values()]
    TF_escore = [tf_data['E-score'] for tf_data in sorted_dict_TF_escore.values()]

    plt.bar(TF_description, TF_escore, width=0.5, color='cornflowerblue', edgecolor='blue')
    plt.axhline(y=0.4, color='r', linestyle='--')
    plt.xticks(rotation=45, ha='right')
    plt.xlabel('Transcription Factors (TFs)')
    plt.ylabel('Max Escore')
    plt.title('Max Escore for Each TF, seq: {}'.format(seq_info))

    plt.tight_layout()
    plt.savefig("{}summarized.jpg".format(save_path), bbox_inches='tight')
    plt.close()
    for TF_max in dict_TF_escore:
        if dict_TF_escore[TF_max]["Contain Mut"]:
            print(TF_max, "Contain Mut!")

def get_escore_info_per_TF(organ_dict):
    tf_names = []
    location_per_TF = {}
    Escore_per_TF = {}
    for loc, TF_list_in_loc in organ_dict.items():
        for TF_dict_in_loc in TF_list_in_loc:
            if TF_dict_in_loc["E-score"] >= TF_BINDING_THRESHOLD:
                TF_name = TF_dict_in_loc["TF Name"]
                if TF_name not in tf_names:
                    tf_names.append(TF_name)
                if TF_name not in location_per_TF:
                    location_per_TF[TF_name] = []
                if TF_name not in Escore_per_TF:
                    Escore_per_TF[TF_name] = []
                location_per_TF[TF_name].append(loc)
                Escore_per_TF[TF_name].append(TF_dict_in_loc["E-score"])
    return location_per_TF, Escore_per_TF

def summarized_figure_for_seq_all_binding_dots(human_all_motif_bind, chimp_all_motif_bind, mutation_location_range, save_path, seq_info):
    location_per_TF_human, Escore_per_TF_human = get_escore_info_per_TF(human_all_motif_bind)
    location_per_TF_chimp, Escore_per_TF_chimp = get_escore_info_per_TF(chimp_all_motif_bind)

    plt.figure(figsize=(20, 15))
    plt.suptitle("Escore of TFs binding to seq: {}".format(seq_info))
    plt.subplot(2, 1, 1)
    color_counter = 0
    for TF_name in location_per_TF_human:
        plt.plot(location_per_TF_human[TF_name], Escore_per_TF_human[TF_name], 'o', label=TF_name, color=COLORS[color_counter])
        color_counter += 1
    plt.axhline(y=TF_BINDING_THRESHOLD, linestyle='--', color=COLORS[color_counter + 5], label='Escore binding Threshold')
    for mut_location in mutation_location_range:
        plt.axvline(x=mut_location, linestyle='--', color=COLORS[color_counter + 10], label='mutation location')
    plt.xlabel("Sequence")
    plt.ylabel("E-score")
    plt.title("TF binding to Derived seq".format(seq_info))

    plt.subplot(2, 1, 2)
    color_counter = 0
    for TF_name in location_per_TF_chimp:
        plt.plot(location_per_TF_chimp[TF_name], Escore_per_TF_chimp[TF_name], 'o', label=TF_name, color=COLORS[color_counter])
        color_counter += 1
    plt.axhline(y=TF_BINDING_THRESHOLD, linestyle='--', color=COLORS[color_counter + 5], label='Escore binding Threshold')
    for mut_location in mutation_location_range:
        plt.axvline(x=mut_location, linestyle='--', color=COLORS[color_counter + 10], label='mutation location')
    plt.xlabel("Sequence")
    plt.ylabel("E-score")
    plt.title("TF binding to Ancestral seq".format(seq_info))
    plt.legend(loc='upper center', bbox_to_anchor=(0.5, -0.2), fontsize=12, frameon=False, ncol=7)
    plt.tight_layout()
    plt.savefig("{}summarized_dots.jpg".format(save_path), bbox_inches='tight')
    plt.close()


def get_pwm_paths(folder_path):
    folder_files_png = []
    for root, dirs, files in os.walk(folder_path):
        for file in files:
            if file.endswith('.png'):
                full_path = os.path.join(root, file)
                if not 'RC' in full_path and not 'secondary' in full_path:
                    folder_files_png.append(full_path)
    return folder_files_png


def get_pwm_dict(folder_path):
    pwm_paths = get_pwm_paths(folder_path)
    pwm_dict = {}
    for pwm in pwm_paths:
        TF_name = pwm.split('/')[7]
        pwm_dict[TF_name] = Image.open(pwm)
    return pwm_dict


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
    return gene_name, human_chr, human_seq_start, human_seq_end, range(mutation_loc_start,
                                                                       mutation_loc_end), nuc_human, nuc_ancestral


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
    return gene_name, human_chr, human_seq_start, human_seq_end, range(mutation_loc_start, mutation_loc_end), nuc_human, nuc_variant, derived_species, ancestral_species


def create_max_diff_per_TF(path_csv):
    TF_above_threshold = {}
    with open("{}TF_near_mut_data.csv".format(path_csv), 'r') as f:
        TF_above_threshold_data = f.read()

    TF_above_threshold_data = TF_above_threshold_data.split("\n")[1:]

    for TF_row in TF_above_threshold_data:
        TF_row = TF_row.split(",")
        if TF_row != ['']:
            TF_name = TF_row[0]
            TF_location_start = TF_row[1]
            TF_location_end = TF_row[2]
            TF_diff_zscore = float(TF_row[7])

            if TF_name not in TF_above_threshold:
                TF_above_threshold[TF_name] = {"TF_diff_zscore": 0.0}

            if TF_above_threshold[TF_name]["TF_diff_zscore"] < abs(TF_diff_zscore):
                TF_above_threshold[TF_name]["TF_location_start"] = TF_location_start
                TF_above_threshold[TF_name]["TF_location_end"] = TF_location_end
                TF_above_threshold[TF_name]["TF_diff_zscore"] = TF_diff_zscore

    return TF_above_threshold


def find_all_TF_above_tresh(TF_in_seq_dict, TF_ancestral_seq_dict):
    TF_to_escore_above_thresh = {}
    for location in TF_in_seq_dict.keys():
        for TF in TF_in_seq_dict[location]:
            if TF["E-score"] > TF_BINDING_THRESHOLD:
                if TF["TF Name"] not in TF_to_escore_above_thresh:
                    TF_to_escore_above_thresh[TF["TF Name"]] = []
                TF_to_escore_above_thresh[TF["TF Name"]].append((location, len(TF["motif"])))
    for location in TF_ancestral_seq_dict.keys():
        for TF in TF_in_seq_dict[location]:
            if TF["E-score"] > TF_BINDING_THRESHOLD:
                if TF["TF Name"] not in TF_to_escore_above_thresh:
                    TF_to_escore_above_thresh[TF["TF Name"]] = []
                if (location, len(TF["motif"])) not in TF_to_escore_above_thresh[TF["TF Name"]]:
                    TF_to_escore_above_thresh[TF["TF Name"]].append((location, len(TF["motif"])))
    return TF_to_escore_above_thresh


def plot_zscore_diff(max_diff_TF, save_path):
    up_diff = "powderblue"
    down_diff = "thistle"

    TF_names = list(max_diff_TF.keys())
    values = [max_diff_TF[TF]["TF_diff_zscore"] for TF in TF_names]

    combined = list(zip(TF_names, values))
    sorted_combined = sorted(combined, key=lambda x: abs(x[1]), reverse=True)
    sorted_TF_names, sorted_values = zip(*sorted_combined)
    sorted_TF_names = list(sorted_TF_names)
    sorted_values = list(sorted_values)

    colors = [up_diff if val > 0 else down_diff for val in sorted_values]
    sorted_TF_names = [name.replace("_mus_musculus", "") for name in sorted_TF_names]

    plt.figure(figsize=(20, 15))
    plt.bar(sorted_TF_names, sorted_values, color=colors, width=0.3, align='center')
    plt.axhline(y=0, color="black")
    plt.xlabel("TFs with E_score above {} binding near variant location".format(TF_BINDING_THRESHOLD))
    plt.ylabel("Z-score diff values")
    plt.title("Z score diff of TFs near mutation")

    if len(TF_names) == 1:
        plt.xlim(-0.5, 0.5)

    up_legend = mlines.Line2D([], [], marker='o', color=up_diff, markersize=10, label='Positive difference')
    down_legend = mlines.Line2D([], [], marker='o', color=down_diff, markersize=10, label='Negative difference')

    plt.legend(handles=[up_legend, down_legend], bbox_to_anchor=(1, 0.5), loc='center left', fontsize=10, ncol=1)
    plt.xticks(rotation=45, ha='right', fontsize=8)

    plt.savefig("{}only_near_mut_zscore_change.jpg".format(save_path), bbox_inches='tight')
    plt.close()


def summarized_all_TF_binds_to_seq(TF_in_seq_dict, TF_ancestral_seq_dict, mut_location_range, seq_info, save_path, significant_TF, color_gradient, csv_path):
    TF_to_escore_above_thresh = find_all_TF_above_tresh(TF_in_seq_dict, TF_ancestral_seq_dict)

    TF_above_threshold_max = create_max_diff_per_TF(csv_path)
    max_diff_color = {}
    for TF in TF_above_threshold_max:
        max_diff_color[TF] = select_color_from_gradient(TF_above_threshold_max[TF]["TF_diff_zscore"], color_gradient)

    plt.figure(figsize=(15, 10))
    counter = 1
    color_counter = 0
    color_background = ['white', 'lightgrey']
    all_locations = list(TF_in_seq_dict.keys())
    for TF in TF_to_escore_above_thresh.keys():
        TF_locations = [loc_mer[0] for loc_mer in TF_to_escore_above_thresh[TF]]
        TF_mer = [loc_mer[1] for loc_mer in TF_to_escore_above_thresh[TF]]
        plt.plot(all_locations, [counter] * len(all_locations), linestyle='None', marker='s', color=color_background[counter % 2])
        for index in range(len(TF_locations)):
            plt.plot(range(TF_locations[index], TF_locations[index] + TF_mer[index]), [counter]*TF_mer[index], linestyle='None', marker='s', color=COLORS[color_counter])
            if TF in significant_TF:
                if TF_locations[index] in significant_TF[TF]:
                    print("location significant", TF, TF_locations[index])
                    plt.plot([TF_locations[index], TF_locations[index]], [counter - 0.5, counter + 0.5], color="black", zorder=3)
                    plt.plot([TF_locations[index] + TF_mer[index], TF_locations[index] + TF_mer[index]], [counter - 0.5, counter + 0.5], color="black", zorder=3)
                    plt.plot([TF_locations[index], TF_locations[index] + TF_mer[index]], [counter - 0.5, counter - 0.5], color="black", zorder=3)
                    plt.plot([TF_locations[index], TF_locations[index] + TF_mer[index]], [counter + 0.5, counter + 0.5], color="black", zorder=3)
        counter += 1
        color_counter += 1
    for mut_location in mut_location_range:
        plt.axvline(x=mut_location, linestyle='--', color=COLORS[color_counter + 10], label='mutation location')
    plt.xlabel(seq_info)
    plt.ylabel("TFs binding to Seq")
    plt.xticks([])
    TF_list_for_plot = list(TF_to_escore_above_thresh.keys())
    TF_list_for_plot = [name.replace("_mus_musculus", "") for name in TF_list_for_plot]
    plt.yticks(range(1, counter), TF_list_for_plot, fontsize=8)
    plt.title("TFs binding to seq (E-score above {})".format(TF_BINDING_THRESHOLD))
    plt.savefig("{}summarized_all_TFs.jpg".format(save_path), bbox_inches='tight')
    plt.close()

    # plotting Z score diff
    plot_zscore_diff(TF_above_threshold_max, save_path)


def read_csv_TF_bind_to_seq(csv_path):
    TF_bind_dict_human = {}
    TF_bind_dict_chimp = {}
    with open(csv_path, 'r') as f:
        all_data = f.read()
    all_data = all_data.split('\n')[1:]
    for row in all_data:
        row = row.split(',')
        if row[-1] == 'derived':
            if int(row[0]) not in TF_bind_dict_human:
                TF_bind_dict_human[int(row[0])] = []
            TF_bind_dict_human[int(row[0])].append({'TF Name': row[2],
                                                    'Nuc': row[3],
                                                    'E-score': float(row[7]),
                                                    'Median': float(row[8]),
                                                    'Z-score': float(row[9]),
                                                    'motif': row[4],
                                                    'motif-opp': row[5],
                                                    'Contain Mut': bool(row[10] == 'True'),
                                                    'Is opp': bool(row[6]) == 'True'})
        if row[-1] == 'ancestral':
            if int(row[0]) not in TF_bind_dict_chimp:
                TF_bind_dict_chimp[int(row[0])] = []
            TF_bind_dict_chimp[int(row[0])].append({'TF Name': row[2],
                                                    'Nuc': row[3],
                                                    'E-score': float(row[7]),
                                                    'Median': float(row[8]),
                                                    'Z-score': float(row[9]),
                                                    'motif': row[4],
                                                    'motif-opp': row[5],
                                                    'Contain Mut': row[10] == 'True',
                                                    'Is opp': row[6] == 'True'})
    return TF_bind_dict_human, TF_bind_dict_chimp


def read_csv_of_mean_intensity(csv_path):
    mean_intensity = {}
    with open(csv_path, 'r') as f:
        all_data = f.read()
    all_data = all_data.split('\n')[1:]
    for row in all_data:
        if row != '':
            row = row.split(',')
            mean_intensity[row[0]] = float(row[1])
    return mean_intensity


def save_statistical_values(human_all_motif_bind, chimp_all_motif_bind, csv_save_path):
    TF_list_pval_escore_check = {}
    summary_csv_stat = "TF,location,escore_derived,escore_ansectral,diff_zscore,p_val,p_val_cdf,corrected_pval\n"
    for location in human_all_motif_bind:
        for TF in human_all_motif_bind[location]:
            if TF["Contain Mut"]:
                ancestral_dict = next((d for d in chimp_all_motif_bind[location] if d['TF Name'] == TF["TF Name"]), None)
                if TF["TF Name"] not in TF_list_pval_escore_check:
                    TF_list_pval_escore_check[TF["TF Name"]] = {}
                    TF_list_pval_escore_check[TF["TF Name"]]["p_values"] = []
                    TF_list_pval_escore_check[TF["TF Name"]]["p_values_cdf"] = []
                    TF_list_pval_escore_check[TF["TF Name"]][("escore_derived")] = []
                    TF_list_pval_escore_check[TF["TF Name"]][("escore_ancestral")] = []
                    TF_list_pval_escore_check[TF["TF Name"]][("zscore_diff")] = []
                    TF_list_pval_escore_check[TF["TF Name"]]["locations"] = []
                zscore_diff = float(TF["Z-score"]) - float(ancestral_dict["Z-score"])
                zscores_pvalue = stats.norm.sf(np.abs(zscore_diff)) * 2
                Zscore_cdf = 2 * (1 - stats.norm.cdf(abs(zscore_diff)))
                TF_list_pval_escore_check[TF["TF Name"]]["p_values"].append(float(zscores_pvalue))
                TF_list_pval_escore_check[TF["TF Name"]]["p_values_cdf"].append(float(Zscore_cdf))
                TF_list_pval_escore_check[TF["TF Name"]][("escore_derived")].append(TF["E-score"])
                TF_list_pval_escore_check[TF["TF Name"]][("escore_ancestral")].append(ancestral_dict["E-score"])
                TF_list_pval_escore_check[TF["TF Name"]][("zscore_diff")].append(zscore_diff)
                TF_list_pval_escore_check[TF["TF Name"]]["locations"].append(location)

    significant_TF = {}
    for TF in TF_list_pval_escore_check:
        _, zscores_p_values_corrected, _, _ = multipletests(TF_list_pval_escore_check[TF]["p_values"], method='fdr_bh')
        TF_list_pval_escore_check[TF]["p_values_corrected"] = zscores_p_values_corrected
        for index in range(len(zscores_p_values_corrected)):
            if zscores_p_values_corrected[index] < 0.05 and (float(TF_list_pval_escore_check[TF]["escore_derived"][index]) >= TF_BINDING_THRESHOLD or float(TF_list_pval_escore_check[TF]["escore_ancestral"][index]) >= TF_BINDING_THRESHOLD):
                if TF not in significant_TF:
                    significant_TF[TF] = []
                significant_TF[TF].append(TF_list_pval_escore_check[TF]["locations"][index])
                summary_csv_stat += (TF + "," +
                                     str(TF_list_pval_escore_check[TF]["locations"][index]) + "," +
                                     str(TF_list_pval_escore_check[TF]["escore_derived"][index]) + "," +
                                     str(TF_list_pval_escore_check[TF]["escore_ancestral"][index]) + "," +
                                     str(TF_list_pval_escore_check[TF][("zscore_diff")][index]) + "," +
                                     str(TF_list_pval_escore_check[TF]["p_values"][index]) + "," +
                                     str(zscores_p_values_corrected[index]) + "," +
                                     str(TF_list_pval_escore_check[TF]["p_values_cdf"][index]) + "\n")

    with open("{}TF_pass_statistics.csv".format(csv_save_path), 'w') as f:
        f.write(summary_csv_stat)

    return significant_TF


def write_pwm_paths_to_csv(csv_path, pwm_dict):
    with open(csv_path, 'r') as f_read:
        f_read.read()
    f_read_split = f_read.split('\n')
    f_write = f_read_split[:1] + "pwm_path" + "\n"
    for row in f_read_split[1:]:
        row_split = row.split(',')
        TF_name, pbm_path = row_split[0], row_split[1]
        if TF_name in pwm_dict:
            f_write += TF_name + "," + pbm_path + "," + pwm_dict[TF_name] + "\n"
        else:
            f_write += TF_name + "," + pbm_path + "\n"

    with open(csv_path, 'w') as f:
        f.write(f_write)


def main():
    get_colors()

    print("### Getting color scale ###")
    gradients = {}
    with open(COLOR_CONFIG_FILE, "r") as file:
        config_content = file.read()
        exec(config_content, gradients)

    gradient_1 = gradients['gradient_1']
    """
    print("##### Getting human seq from human genome #####")

    print("##### Getting all variants input paths #####")
    input_files = []
    for path in os.listdir(INPUTS_FOLDER):
        input_files.append(path)
    """

    print("##### Getting all variants input paths #####")
    input_file_path = "{}{}".format(INPUTS_FOLDER, INPUT_CSV)
    data_input_file = open_txt(input_file_path)
    input_files = data_input_file.split("\n")[1:]

    #curr_variant_input_file = "chr6_396418_A_G.txt"
    #curr_variant_input_file = "chr6_396321_C_T.txt"

    print("##### Extracting pwm files from TF folders #####")
    dict_pwm = get_pwm_dict(DATA_TF_PATH)
    print()

    for curr_variant_input_row in input_files:
        if curr_variant_input_row != "":
            """
            input_file = "{}{}".format(INPUTS_FOLDER, curr_variant_input_file)
            gene_name, human_chr, human_seq_start, human_seq_end, mutation_loc_range, nuc_human, nuc_ancestral = get_parameters_from_input(input_file)
            """
            gene_name, human_chr, human_seq_start, human_seq_end, mutation_loc_range, nuc_human, nuc_ancestral, derived_species, ancestral_species = get_parameters_from_input_csv(
                curr_variant_input_row)

            seq_info = "Seq: {}, {}: {}-{}, mutation: {}-{} {}->{}".format(HUMAN_GENOME_FILE.split(".")[0], human_chr,
                                                                             human_seq_start, human_seq_end, mutation_loc_range[0], mutation_loc_range[-1],
                                                                             nuc_human, nuc_ancestral)

            seq_folder = "{}{}_{}_{}_{}-{}_mutation_{}-{}_{}_{}".format(PATH_SAVE_RESULTS, HUMAN_GENOME_FILE.split(".")[0],
                                                                        gene_name,
                                                                        human_chr, human_seq_start, human_seq_end,
                                                                        mutation_loc_range[0], mutation_loc_range[-1],
                                                                        nuc_human, nuc_ancestral)

            csv_path = seq_folder + SLASH + "TFs_binding_csv" + SLASH + "all_TF_info.csv"

            human_all_motif_bind, chimp_all_motif_bind = read_csv_TF_bind_to_seq(csv_path)
            print("#### Seq: {}, {}: {}-{}, mutation: {}-{} {}->{} ####".format(HUMAN_GENOME_FILE.split(".")[0], human_chr,
                                                                                human_seq_start, human_seq_end,
                                                                                mutation_loc_range[0],
                                                                                mutation_loc_range[-1],
                                                                                nuc_human, nuc_ancestral))

            csv_mean_path = seq_folder + SLASH + "TFs_binding_csv" + SLASH + "mean_intensity.csv"
            mean_intensity_for_TF = read_csv_of_mean_intensity(csv_mean_path)

            # seq_folder = "{}{}_{}_{}-{}_mutation_{}_{}_{}".format(PATH_SAVE_RESULTS, HUMAN_GENOME_FILE.split(".")[0], human_chr, human_seq_start, human_seq_end, mutation_loc, nuc_human, nuc_ancestral)

            significant_TF = save_statistical_values(human_all_motif_bind, chimp_all_motif_bind,
                                          seq_folder + SLASH + "TFs_binding_csv" + SLASH)

            print("##### Plotting figure for each TF")
            plots_save_path = "{}{}TFs_binding_plots".format(seq_folder, SLASH)

            plot_TF_escore_and_median_in_seq_nice(human_all_motif_bind, chimp_all_motif_bind, mutation_loc_range, mean_intensity_for_TF,
                                             seq_info, plots_save_path + SLASH, significant_TF, derived_species, ancestral_species)
            concatenate_TF_plot_pwn(dict_pwm, plots_save_path + SLASH)

            print("##### Saving summarized figure")
            summarized_save_path = "{}{}summarized_plots".format(seq_folder, SLASH)

            # summarized_figure_for_seq(human_all_motif_bind, chimp_all_motif_bind, mutation_loc, summarized_save_path + SLASH, seq_info)
            summarized_figure_for_seq_all_binding_dots(human_all_motif_bind, chimp_all_motif_bind, mutation_loc_range, summarized_save_path + SLASH,
                                                       seq_info)
            summarized_all_TF_binds_to_seq(human_all_motif_bind, chimp_all_motif_bind, mutation_loc_range, seq_info, summarized_save_path + SLASH, significant_TF, gradient_1, seq_folder + SLASH + "TFs_binding_csv" + SLASH)

            print()

if __name__ == "__main__":
    main()
