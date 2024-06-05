# import modules
import datetime
import sys
import time
import os
import re
import pysam
from variant_extractor import VariantExtractor
from functions import compare, dict_generator, chr_converter, get_windows  # local scripts
import matplotlib.pyplot as plt
import numpy as np
from statistics import mode, StatisticsError
from enum import Enum


class SomaticVariationType(Enum):
    UNCLASSIFIED = 0
    NORMAL_SINGLE_READ_VARIANT = 1
    TUMORAL_SINGLE_READ_VARIANT = 2
    NORMAL_ONLY_VARIANT = 3
    TUMORAL_ONLY_VARIANT = 4
    TUMORAL_NORMAL_VARIANT = 5


class CalledGenomicVariant:
    TYPE_SNV = 'SNV'
    TYPE_INDEL = 'INDEL'

    def __init__(self, seq_name, pos, end, var_type, length, allele):
        self.seq_name: str = seq_name
        self.pos: int = pos
        self.end: int = end
        self.var_type: str = var_type
        self.length = length
        self.allele: str = allele
        self.somatic_variation_type = SomaticVariationType.UNCLASSIFIED
        self.supporting_reads: int = 1

    def add_supporting_read(self):
        self.supporting_reads += 1

    def __eq__(self, var2):
        if self.seq_name != var2.seq_name:
            return False
        if self.var_type != var2.var_type:
            return False
        if self.pos != var2.pos:
            return False
        if self.end != var2.end:
            return False
        if self.length != var2.length:
            return False
        if self.allele != var2.allele:
            return False
        return True


# Record the start time
start_time = time.time()

# Define input files
# CLI Inputs #
in_N_bam = sys.argv[1]
in_T_bam = sys.argv[2]
in_ref = sys.argv[3]
in_vcf = sys.argv[4]

# Define output files
try:
    out_dir = sys.argv[5]
    os.makedirs(out_dir)
except FileExistsError:
    print('Output files have been generated at an existing directory:')
    print(out_dir)

# Open samfiles with pysam and vcf with VariantExtractor
samfile_N = pysam.AlignmentFile(in_N_bam, 'rb')  # read a bam file
samfile_T = pysam.AlignmentFile(in_T_bam, 'rb')
extractor = VariantExtractor(in_vcf)
ref_genome = pysam.Fastafile(in_ref)

# Create a dict to assign an integer value to each chr name

print("Begin dict.")
chr_dict = dict_generator(ref_genome)
print("End dict.")


def anno_cigar2(indel_dict, snv_dict, cigar_str, md_positions_list, current_ref_pos, aln_chr, datasetIdx, ref_genome,
                aln):
    """This function takes a cigar string and an md_positions_list and annotates the potential SNVs and Indels in
    two different dictionaries.
    datasetIdx: 0 if tumoral, 1 if normal"""

    regexp = r"(?<=[a-zA-Z=])(?=[0-9])|(?<=[0-9])(?=[a-zA-Z=])"  # regexp to split cigar string
    ref_consuming = ['M', 'D', 'N', '=', 'X']  # stores reference consuming cigar operations
    read_consuming_only = ['S', 'H', 'I']  # stores read consuming cigar operations
    cigar_indels = ["I", "D"]  # cigar operations to report as INDEL
    dict_cigar = {"X": "MM", "I": "INS", "D": "DEL"}
    dict_tuples = {0: "M", 1: "I", 2: "D", 3: "N", 4: "S", 5: "H", 6: "P", 7: "=", 8: "X", 9: "B"}
    read_sequence = aln.query_sequence
    ref_sequence = aln.get_reference_sequence()
    read_consumed_bases = 0
    cigar_list = re.split(regexp, cigar_str)  # we split the string into a list
    # we parse the cigar list to get the INDELs
    cigar_list_idx = 0  # to parse cigar_list
    current_ref_cigar_len = 0  # to track the position in the genome
    # We parse the md tag to get the mismatches
    ref_mismatch_positions = []  # list to store the mismatch positions
    md_pos = 0  # to parse md_list
    md_length = 0  # to track the position in the genome
    for symbol in md_positions_list:
        if symbol == "0":  # md separator
            pass  # ignore
        elif symbol[0] == "^":  # deletions
            md_length += len(symbol) - 1
        elif re.match(r'^\d', symbol):  # matches
            md_length += int(symbol)
        else:  # mismatches
            md_length += 1
            ref_mismatch_positions.append(md_length)
        md_pos += 1
    mm_pos_idx = 0  # index to parse ref_mismatch_positions list
    for cigar_list_idx, symbol in enumerate(cigar_list):
        if symbol.isdigit():
            cigar_op = cigar_list[cigar_list_idx + 1]
            if cigar_op in cigar_indels:
                true_pos = current_ref_pos + current_ref_cigar_len
                pos_init = true_pos
                pos_end = int(true_pos) + int(symbol)
                var_type = cigar_op
                length = symbol
                called_indel = CalledGenomicVariant(variant_chr, pos_init, pos_end, CalledGenomicVariant.TYPE_INDEL,
                                                    int(symbol), "")
                # we store the variation in the dict
                if called_indel.pos not in indel_dict:
                    indel_dict[called_indel.pos] = []
                indel_pos_list = indel_dict[called_indel.pos]
                indel_exists = False
                for var_indel in indel_pos_list:  # we check whether the indel exists
                    if called_indel == var_indel:
                        called_indel = var_indel
                        indel_exists = True
                        break
                if not indel_exists:
                    if datasetIdx == 0:
                        called_indel.somatic_variation_type = SomaticVariationType.TUMORAL_SINGLE_READ_VARIANT
                    if datasetIdx == 1:
                        called_indel.somatic_variation_type = SomaticVariationType.NORMAL_SINGLE_READ_VARIANT
                    indel_dict[called_indel.pos].append(called_indel)
                else:  # if the variation is already annotated in the dict
                    var_code = called_indel.somatic_variation_type
                    called_indel.add_supporting_read()  # we add 1 supporting read
                    if datasetIdx == 0:
                        if var_code == SomaticVariationType.NORMAL_SINGLE_READ_VARIANT or var_code == SomaticVariationType.NORMAL_ONLY_VARIANT:
                            called_indel.somatic_variation_type = SomaticVariationType.TUMORAL_NORMAL_VARIANT
                        if var_code == SomaticVariationType.TUMORAL_SINGLE_READ_VARIANT:
                            called_indel.somatic_variation_type = SomaticVariationType.TUMORAL_ONLY_VARIANT
                    if datasetIdx == 1:
                        if var_code == SomaticVariationType.TUMORAL_SINGLE_READ_VARIANT or var_code == SomaticVariationType.TUMORAL_ONLY_VARIANT:
                            called_indel.somatic_variation_type = SomaticVariationType.TUMORAL_NORMAL_VARIANT
                        if var_code == SomaticVariationType.NORMAL_SINGLE_READ_VARIANT:
                            called_indel.somatic_variation_type = SomaticVariationType.NORMAL_ONLY_VARIANT
            if cigar_op in ref_consuming:
                current_ref_cigar_len += int(symbol)
            if cigar_op == 'M' and len(
                    ref_mismatch_positions) > 0:  # works only if mismatches and matches are coded as M
                if len(ref_mismatch_positions) == 0 or mm_pos_idx >= len(ref_mismatch_positions):
                    continue
                mm_ref_pos = ref_mismatch_positions[mm_pos_idx]
                while mm_ref_pos < current_ref_cigar_len and mm_pos_idx < len(ref_mismatch_positions):
                    pos_in_read = mm_ref_pos + read_consumed_bases
                    substitution = read_sequence[pos_in_read - 1]
                    pos_snv = current_ref_pos + mm_ref_pos
                    # we store the variation in the dict
                    if substitution != 'N':  # N = unknown base
                        called_snv = CalledGenomicVariant(variant_chr, pos_snv, pos_snv, CalledGenomicVariant.TYPE_SNV, 1,
                                                          substitution)
                        if called_snv.pos not in snv_dict:
                            snv_dict[called_snv.pos] = []
                        snv_pos_list = snv_dict[called_snv.pos]
                        snv_exists = False
                        for var_snv in snv_pos_list:  # we check whether the snv exists
                            if called_snv == var_snv:
                                called_snv = var_snv
                                snv_exists = True
                                break
                        if not snv_exists:
                            if datasetIdx == 0:
                                called_snv.somatic_variation_type = SomaticVariationType.TUMORAL_SINGLE_READ_VARIANT
                            if datasetIdx == 1:
                                called_snv.somatic_variation_type = SomaticVariationType.NORMAL_SINGLE_READ_VARIANT
                            snv_pos_list.append(called_snv)
                        else:
                            var_code = called_snv.somatic_variation_type
                            called_snv.add_supporting_read()  # we add 1 supporting read
                            if datasetIdx == 0:
                                if var_code == SomaticVariationType.NORMAL_SINGLE_READ_VARIANT or var_code == SomaticVariationType.NORMAL_ONLY_VARIANT:
                                    called_snv.somatic_variation_type = SomaticVariationType.TUMORAL_NORMAL_VARIANT
                                if var_code == SomaticVariationType.TUMORAL_SINGLE_READ_VARIANT:
                                    called_snv.somatic_variation_type = SomaticVariationType.TUMORAL_ONLY_VARIANT
                            if datasetIdx == 1:
                                if var_code == SomaticVariationType.TUMORAL_SINGLE_READ_VARIANT or var_code == SomaticVariationType.TUMORAL_ONLY_VARIANT:
                                    called_snv.somatic_variation_type = SomaticVariationType.TUMORAL_NORMAL_VARIANT
                                if var_code == SomaticVariationType.NORMAL_SINGLE_READ_VARIANT:
                                    called_snv.somatic_variation_type = SomaticVariationType.NORMAL_ONLY_VARIANT
                    mm_pos_idx += 1
                    if mm_pos_idx < len(ref_mismatch_positions):
                        mm_ref_pos = ref_mismatch_positions[mm_pos_idx]
            if cigar_op in read_consuming_only:
                read_consumed_bases += int(symbol)
            if cigar_op == 'D':  # D = DEL
                read_consumed_bases -= int(symbol)


# Create a 2000 kb window for each variant
windows = get_windows(extractor, ref_genome, window_size=2000)

total_windows = 0
for chr in windows:
    total_windows += len(windows[chr])
# print(total_windows)

INDEL_counts = ([], [], [], [], [], [])  # tuple to store the INDELs
SNV_counts = ([], [], [], [], [], [])  # tuple to store the SNVs
INDEL_supp_reads_count = ([], [], [], [], [], [])
SNV_supp_reads_count = ([], [], [], [], [], [])

# Parse vcf and sam files in parallel
dict_indel_count = {}  # initialize dict to store INDELs count
dict_snv_count = {}  # init dict to store SNVs count

def pick_variant(var_list):
    """This function is meant to give priority to variants based on its variation somatic type
    in cases when there are multiple variant types in the same position."""
    return max(var_list, key=lambda var: var.somatic_variation_type.value)

for seq_name, seq_windows in windows.items():
    # Print the chromosome name
    print('Chr:{}'.format(seq_name))

    for window in seq_windows:

        #print(window)
        start_pos, end_pos, variant = window
        variant_chr = variant.contig
        variant_pos = variant.pos
        variant_end = variant.end
        tumor_aln_iter = samfile_T.fetch(seq_name, start_pos, end_pos)
        normal_aln_iter = samfile_N.fetch(seq_name, start_pos, end_pos)
        for iterIdx, current_iter in enumerate((tumor_aln_iter, normal_aln_iter)):
            for current_aln in current_iter:
                if current_aln.cigarstring is None:
                    continue  # go back to for loop beginning
                cigar = str(current_aln.cigarstring)  # we get the cigar string
                md_tag = current_aln.get_tag("MD", with_value_type=False)
                pattern_md = r'0|\^[A-Z]+|[A-Z]|[0-9]+'
                md_list = re.findall(pattern_md, md_tag)
                ref_pos = current_aln.reference_start  # get the reference position where the cigar begins
                aln_chr = current_aln.reference_name  # we get the chr
                anno_cigar2(dict_indel_count, dict_snv_count, cigar, md_list, ref_pos, aln_chr, iterIdx, ref_genome,
                            current_aln)
            if iterIdx == 1:
                single_read_variant_counts = 0
                tumoral_only_counts = 0
                tumoral_normal_variant_counts = 0
                normal_only_variants_counts = 0
                for k, v in dict_indel_count.items():
                    var_to_count: CalledGenomicVariant = pick_variant(v)
                    if (var_to_count.somatic_variation_type == SomaticVariationType.NORMAL_SINGLE_READ_VARIANT or
                            var_to_count.somatic_variation_type == SomaticVariationType.TUMORAL_SINGLE_READ_VARIANT):
                        single_read_variant_counts += 1
                    if var_to_count.somatic_variation_type == SomaticVariationType.TUMORAL_ONLY_VARIANT:
                        tumoral_only_counts += 1
                    if var_to_count.somatic_variation_type == SomaticVariationType.TUMORAL_NORMAL_VARIANT:
                        tumoral_normal_variant_counts += 1
                    if var_to_count.somatic_variation_type == SomaticVariationType.NORMAL_ONLY_VARIANT:
                        normal_only_variants_counts += 1
                    INDEL_supp_reads_count[var_to_count.somatic_variation_type.value].append(var_to_count.supporting_reads)
                INDEL_counts[SomaticVariationType.TUMORAL_SINGLE_READ_VARIANT.value].append(single_read_variant_counts)
                INDEL_counts[SomaticVariationType.TUMORAL_ONLY_VARIANT.value].append(tumoral_only_counts)
                INDEL_counts[SomaticVariationType.TUMORAL_NORMAL_VARIANT.value].append(tumoral_normal_variant_counts)
                INDEL_counts[SomaticVariationType.NORMAL_ONLY_VARIANT.value].append(normal_only_variants_counts)
                dict_indel_count = {}  # re-init dict
                single_read_variant_counts = 0
                tumoral_only_counts = 0
                tumoral_normal_variant_counts = 0
                normal_only_variants_counts = 0
                for k, v in dict_snv_count.items():
                    var_to_count: CalledGenomicVariant = pick_variant(v)
                    if (var_to_count.somatic_variation_type == SomaticVariationType.NORMAL_SINGLE_READ_VARIANT or
                            var_to_count.somatic_variation_type == SomaticVariationType.TUMORAL_SINGLE_READ_VARIANT):
                        single_read_variant_counts += 1
                    if var_to_count.somatic_variation_type == SomaticVariationType.TUMORAL_ONLY_VARIANT:
                        tumoral_only_counts += 1
                    if var_to_count.somatic_variation_type == SomaticVariationType.TUMORAL_NORMAL_VARIANT:
                        tumoral_normal_variant_counts += 1
                    if var_to_count.somatic_variation_type == SomaticVariationType.NORMAL_ONLY_VARIANT:
                        normal_only_variants_counts += 1
                    SNV_supp_reads_count[var_to_count.somatic_variation_type.value].append(var_to_count.supporting_reads)
                SNV_counts[SomaticVariationType.TUMORAL_SINGLE_READ_VARIANT.value].append(single_read_variant_counts)
                SNV_counts[SomaticVariationType.TUMORAL_ONLY_VARIANT.value].append(tumoral_only_counts)
                SNV_counts[SomaticVariationType.TUMORAL_NORMAL_VARIANT.value].append(tumoral_normal_variant_counts)
                SNV_counts[SomaticVariationType.NORMAL_ONLY_VARIANT.value].append(normal_only_variants_counts)
                dict_snv_count = {}  # re-init dict

# Define the data to plot
data1 = INDEL_counts[SomaticVariationType.TUMORAL_NORMAL_VARIANT.value]
data2 = SNV_counts[SomaticVariationType.TUMORAL_NORMAL_VARIANT.value]
data3 = INDEL_counts[
    SomaticVariationType.TUMORAL_SINGLE_READ_VARIANT.value]  # It's both tumoral and normal single reads
data4 = SNV_counts[
    SomaticVariationType.TUMORAL_SINGLE_READ_VARIANT.value]  # It's both tumoral and normal single reads
data5 = INDEL_counts[SomaticVariationType.TUMORAL_ONLY_VARIANT.value]
data6 = SNV_counts[SomaticVariationType.TUMORAL_ONLY_VARIANT.value]
data7 = INDEL_counts[SomaticVariationType.NORMAL_ONLY_VARIANT.value]
data8 = SNV_counts[SomaticVariationType.NORMAL_ONLY_VARIANT.value]

# Process data to plot (remove outliers): std method

raw_data = [data1, data2, data3, data4, data5, data6, data7, data8]
processed_data = []
multiplier = 1.5

for i in raw_data:
    # Calculate mean and standard deviation
    mean = np.mean(i)
    std_dev = np.std(i)

    # Define thresholds
    lower_threshold = mean - multiplier * std_dev
    upper_threshold = mean + multiplier * std_dev

    # Filter samples based on the standard deviation method
    filtered_samples = [x for x in i if x >= lower_threshold and x <= upper_threshold]

    # Add filtered samples to processed data
    processed_data.append(filtered_samples)

# Process data to remove 0s

data_non0 = []
for i in processed_data:
    # Remove zeros from the dataset using list comprehension
    non_zero_data = [x for x in i if x != 0]

    # Add non-zero data to the data_non0 list
    data_non0.append(non_zero_data)

# Normalize the data with log10 transformation

normalized_data = []
for data in processed_data:
    normalized_data.append([np.log10(x) if x > 0 else 0 for x in data])

# Create a single figure with subplots
fig, axes = plt.subplots(4, 2, figsize=(24, 24))

# Plot for potential germinal INDELs
#maxim1 = np.max(processed_data[0])
moda1 = mode(processed_data[0])
axes[0][0].hist(processed_data[0], bins=range(int(min(processed_data[0])), int(max(processed_data[0])) + 2), color='blue', alpha=0.7)
axes[0][0].set_xlim(left=0, right=15)
axes[0][0].scatter(moda1, 0, color="red", marker='o', label="Mode: {}".format(moda1))
#axes[0][0].scatter(maxim1, 0, color="red", marker='o', label="Max: {}".format(maxim1))
axes[0][0].set_xlabel('Number of INDELs/window', fontsize=16)
axes[0][0].set_ylabel('Count', fontsize=16)
axes[0][0].set_title('Potential germline INDELs', fontsize=20)
axes[0][0].legend()

# Plot for potential germinal SNVs
#maxim2 = np.max(processed_data[1])
moda2 = mode(processed_data[1])
axes[0][1].hist(processed_data[1], bins=range(int(min(processed_data[1])), int(max(processed_data[1])) + 2), color='blue', alpha=0.7)
axes[0][1].set_xlim(left=0, right=500)
axes[0][1].scatter(moda2, 0, color="red", marker='o', label="Mode: {}".format(moda2))
#axes[0][1].scatter(maxim2, 0, color="red", marker='o', label="Max: {}".format(maxim2))
axes[0][1].set_xlabel('Number of SNVs/window', fontsize=16)
axes[0][1].set_ylabel('Count', fontsize=16)
axes[0][1].set_title('Potential germline SNVs', fontsize=20)
axes[0][1].legend()

# Plot for single-read INDELs
#maxim3 = np.max(processed_data[2])
moda3 = mode(processed_data[2])
axes[1][0].hist(processed_data[2], bins=range(int(min(processed_data[2])), int(max(processed_data[2])) + 2), color='blue', alpha=0.7)
axes[1][0].set_xlim(left=0, right=50)
axes[1][0].scatter(moda3, 0, color="red", marker='o', label="Mode: {}".format(moda3))
#axes[1][0].scatter(maxim3, 0, color="red", marker='o', label="Max: {}".format(maxim3))
axes[1][0].set_xlabel('Number of INDELs/window', fontsize=16)
axes[1][0].set_ylabel('Count', fontsize=16)
axes[1][0].set_title('Single-read INDELs', fontsize=20)
axes[1][0].legend()

# Plot for single-read SNVs
#maxim4 = np.max(processed_data[3])
moda4 = mode(processed_data[3])
axes[1][1].hist(processed_data[3], bins=range(int(min(processed_data[3])), int(max(processed_data[3])) + 2), color='blue', alpha=0.7)
axes[1][1].set_xlim(left=0, right=500)
axes[1][1].scatter(moda4, 0, color="red", marker='o', label="Mode: {}".format(moda4))
#axes[1][1].scatter(maxim4, 0, color="red", marker='o', label="Max: {}".format(maxim4))
axes[1][1].set_xlabel('Number of SNVs/window', fontsize=16)
axes[1][1].set_ylabel('Count', fontsize=16)
axes[1][1].set_title('Single-read SNVs', fontsize=20)
axes[1][1].legend()

# Plot for tumoral-only INDELs
#maxim5 = np.max(processed_data[4])
moda5 = mode(processed_data[4])
axes[2][0].hist(processed_data[4], bins=range(int(min(processed_data[4])), int(max(processed_data[4])) + 2), color='blue', alpha=0.7)
axes[2][0].set_xlim(left=0, right=20)
axes[2][0].scatter(moda5, 0, color="red", marker='o', label="Mode: {}".format(moda5))
#axes[2][0].scatter(maxim5, 0, color="red", marker='o', label="Max: {}".format(maxim5))
axes[2][0].set_xlabel('Number of INDELs/window', fontsize=16)
axes[2][0].set_ylabel('Count', fontsize=16)
axes[2][0].set_title('Tumoral-only INDELs', fontsize=20)
axes[2][0].legend()

# Plot for tumoral-only SNVs
#maxim6 = np.max(processed_data[5])
moda6 = mode(processed_data[5])
axes[2][1].hist(processed_data[5], bins=range(int(min(processed_data[5])), int(max(processed_data[5])) + 2), color='blue', alpha=0.7)
axes[2][1].set_xlim(left=0, right=250)
axes[2][1].scatter(moda6, 0, color="red", marker='o', label="Mode: {}".format(moda6))
#axes[2][1].scatter(maxim6, 0, color="red", marker='o', label="Max: {}".format(maxim6))
axes[2][1].set_xlabel('Number of SNVs/window', fontsize=16)
axes[2][1].set_ylabel('Count', fontsize=16)
axes[2][1].set_title('Tumoral-only SNVs', fontsize=20)
axes[2][1].legend()

# Plot for normal-only INDELs
#maxim7 = np.max(processed_data[6])
moda7 = mode(processed_data[6])
axes[3][0].hist(processed_data[6], bins=range(int(min(processed_data[6])), int(max(processed_data[6])) + 2), color='blue', alpha=0.7)
axes[3][0].set_xlim(left=0, right=10)
axes[3][0].scatter(moda7, 0, color="red", marker='o', label="Mode: {}".format(moda7))
#axes[3][0].scatter(maxim7, 0, color="red", marker='o', label="Max: {}".format(maxim7))
axes[3][0].set_xlabel('Number of INDELs/window', fontsize=20)
axes[3][0].set_ylabel('Count', fontsize=16)
axes[3][0].set_title('Normal-only INDELs', fontsize=20)
axes[3][0].legend()

# Plot for normal-only SNVs
#maxim8 = np.max(processed_data[7])
moda8 = mode(processed_data[7])
axes[3][1].hist(processed_data[7], bins=range(int(min(processed_data[7])), int(max(processed_data[7])) + 2), color='blue', alpha=0.7)
axes[3][1].set_xlim(left=0, right=250)
axes[3][1].scatter(moda8, 0, color="red", marker='o', label="Mode: {}".format(moda8))
#axes[3][1].scatter(maxim8, 0, color="red", marker='o', label="Max: {}".format(maxim8))
axes[3][1].set_xlabel('Number of SNVs/window', fontsize=16)
axes[3][1].set_ylabel('Count', fontsize=16)
axes[3][1].set_title('Normal-only SNVs', fontsize=20)
axes[3][1].legend()

# Adjust spacing between subplots
plt.subplots_adjust(hspace=0.5)

# Save the combined plot
plt.savefig('{}/hist_combined.png'.format(out_dir))

# Create a single figure with subplots (raw data)
fig, axes = plt.subplots(4, 2, figsize=(24, 24))

# Plot for potential germinal INDELs
maxim1 = np.max(raw_data[0])
moda1 = mode(raw_data[0])
axes[0][0].hist(raw_data[0], color='blue', alpha=0.7)
# axes[0][0].set_xlim(left=0, right=25)
axes[0][0].scatter(moda1, 0, color="green", marker='o', label="Mode: {}".format(moda1))
axes[0][0].scatter(maxim1, 0, color="red", marker='o', label="Max: {}".format(maxim1))
axes[0][0].set_xlabel('Number of INDELs/window')
axes[0][0].set_ylabel('Count')
axes[0][0].set_title('Potential germline INDELs')
axes[0][0].legend()

# Plot for potential germinal SNVs
maxim2 = np.max(raw_data[1])
moda2 = mode(raw_data[1])
axes[0][1].hist(raw_data[1], color='blue', alpha=0.7)
# axes[0][1].set_xlim(left=0, right=1000)
axes[0][1].scatter(moda2, 0, color="green", marker='o', label="Mode: {}".format(moda2))
axes[0][1].scatter(maxim2, 0, color="red", marker='o', label="Max: {}".format(maxim2))
axes[0][1].set_xlabel('Number of SNVs/window')
axes[0][1].set_ylabel('Count')
axes[0][1].set_title('Potential germline SNVs')
axes[0][1].legend()

# Plot for single-read INDELs
maxim3 = np.max(raw_data[2])
moda3 = mode(raw_data[2])
axes[1][0].hist(raw_data[2], color='blue', alpha=0.7)
# axes[1][0].set_xlim(left=0, right=100)
axes[1][0].scatter(moda3, 0, color="green", marker='o', label="Mode: {}".format(moda3))
axes[1][0].scatter(maxim3, 0, color="red", marker='o', label="Max: {}".format(maxim3))
axes[1][0].set_xlabel('Number of INDELs/window')
axes[1][0].set_ylabel('Count')
axes[1][0].set_title('Single-read INDELs')
axes[1][0].legend()

# Plot for single-read SNVs
maxim4 = np.max(raw_data[3])
moda4 = mode(raw_data[3])
axes[1][1].hist(raw_data[3], color='blue', alpha=0.7)
# axes[1][1].set_xlim(left=0, right=1500)
axes[1][1].scatter(moda4, 0, color="green", marker='o', label="Mode: {}".format(moda4))
axes[1][1].scatter(maxim4, 0, color="red", marker='o', label="Max: {}".format(maxim4))
axes[1][1].set_xlabel('Number of SNVs/window')
axes[1][1].set_ylabel('Count')
axes[1][1].set_title('Single-read SNVs')
axes[1][1].legend()

# Plot for tumoral-only INDELs
maxim5 = np.max(raw_data[4])
moda5 = mode(raw_data[4])
axes[2][0].hist(raw_data[4], color='blue', alpha=0.7)
# axes[2][0].set_xlim(left=0, right=100)
axes[2][0].scatter(moda5, 0, color="green", marker='o', label="Mode: {}".format(moda5))
axes[2][0].scatter(maxim5, 0, color="red", marker='o', label="Max: {}".format(maxim5))
axes[2][0].set_xlabel('Number of INDELs/window')
axes[2][0].set_ylabel('Count')
axes[2][0].set_title('Tumoral-only INDELs')
axes[2][0].legend()

# Plot for tumoral-only SNVs
maxim6 = np.max(raw_data[5])
moda6 = mode(raw_data[5])
axes[2][1].hist(raw_data[5], color='blue', alpha=0.7)
# axes[2][1].set_xlim(left=0, right=1500)
axes[2][1].scatter(moda6, 0, color="green", marker='o', label="Mode: {}".format(moda6))
axes[2][1].scatter(maxim6, 0, color="red", marker='o', label="Max: {}".format(maxim6))
axes[2][1].set_xlabel('Number of SNVs/window')
axes[2][1].set_ylabel('Count')
axes[2][1].set_title('Tumoral-only SNVs')
axes[2][1].legend()

# Plot for normal-only INDELs
maxim7 = np.max(raw_data[6])
moda7 = mode(raw_data[6])
axes[3][0].hist(raw_data[6], color='blue', alpha=0.7)
# axes[3][0].set_xlim(left=0, right=100)
axes[3][0].scatter(moda7, 0, color="green", marker='o', label="Mode: {}".format(moda7))
axes[3][0].scatter(maxim7, 0, color="red", marker='o', label="Max: {}".format(maxim7))
axes[3][0].set_xlabel('Number of INDELs/window')
axes[3][0].set_ylabel('Count')
axes[3][0].set_title('Normal-only INDELs')
axes[3][0].legend()

# Plot for normal-only SNVs
maxim8 = np.max(raw_data[7])
moda8 = mode(raw_data[7])
axes[3][1].hist(raw_data[7], color='blue', alpha=0.7)
# axes[3][1].set_xlim(left=0, right=1500)
axes[3][1].scatter(moda8, 0, color="green", marker='o', label="Mode: {}".format(moda8))
axes[3][1].scatter(maxim8, 0, color="red", marker='o', label="Max: {}".format(maxim8))
axes[3][1].set_xlabel('Number of SNVs/window')
axes[3][1].set_ylabel('Count')
axes[3][1].set_title('Normal-only SNVs')
axes[3][1].legend()

# Adjust spacing between subplots
plt.subplots_adjust(hspace=0.5)

# Save the combined plot
plt.savefig('{}/hist_raw_data.png'.format(out_dir))

# Create a single figure with subplots (plot non0 data)
fig, axes = plt.subplots(4, 2, figsize=(24, 24))

# Plot for potential germinal INDELs
maxim1 = np.max(data_non0[0])
moda1 = mode(data_non0[0])
axes[0][0].hist(data_non0[0], bins=range(int(min(data_non0[0])), int(max(data_non0[0])) + 2), color='blue', alpha=0.7)
# axes[0][0].set_xlim(left=0, right=25)
axes[0][0].scatter(moda1, 0, color="green", marker='o', label="Mode: {}".format(moda1))
axes[0][0].scatter(maxim1, 0, color="red", marker='o', label="Max: {}".format(maxim1))
axes[0][0].set_xlabel('Number of INDELs/window')
axes[0][0].set_ylabel('Count')
axes[0][0].set_title('Potential germline INDEL counts for 2kb windows')
axes[0][0].legend()

# Plot for potential germinal SNVs
maxim2 = np.max(data_non0[1])
moda2 = mode(data_non0[1])
axes[0][1].hist(data_non0[1], bins=range(int(min(data_non0[1])), int(max(data_non0[1])) + 2), color='blue', alpha=0.7)
# axes[0][1].set_xlim(left=0, right=1000)
axes[0][1].scatter(moda2, 0, color="green", marker='o', label="Mode: {}".format(moda2))
axes[0][1].scatter(maxim2, 0, color="red", marker='o', label="Max: {}".format(maxim2))
axes[0][1].set_xlabel('Number of SNVs/window')
axes[0][1].set_ylabel('Count')
axes[0][1].set_title('Potential germline SNV counts for 2kb windows')
axes[0][1].legend()

# Plot for single-read INDELs
maxim3 = np.max(data_non0[2])
moda3 = mode(data_non0[2])
axes[1][0].hist(data_non0[2], bins=range(int(min(data_non0[2])), int(max(data_non0[2])) + 2), color='blue', alpha=0.7)
# axes[1][0].set_xlim(left=0, right=100)
axes[1][0].scatter(moda3, 0, color="green", marker='o', label="Mode: {}".format(moda3))
axes[1][0].scatter(maxim3, 0, color="red", marker='o', label="Max: {}".format(maxim3))
axes[1][0].set_xlabel('Number of INDELs/window')
axes[1][0].set_ylabel('Count')
axes[1][0].set_title('Single-read INDEL counts for 2kb windows')
axes[1][0].legend()

# Plot for single-read SNVs
maxim4 = np.max(data_non0[3])
moda4 = mode(data_non0[3])
axes[1][1].hist(data_non0[3], bins=range(int(min(data_non0[3])), int(max(data_non0[3])) + 2), color='blue', alpha=0.7)
# axes[1][1].set_xlim(left=0, right=1500)
axes[1][1].scatter(moda4, 0, color="green", marker='o', label="Mode: {}".format(moda4))
axes[1][1].scatter(maxim4, 0, color="red", marker='o', label="Max: {}".format(maxim4))
axes[1][1].set_xlabel('Number of SNVs/window')
axes[1][1].set_ylabel('Count')
axes[1][1].set_title('Single-read SNV counts for 2kb windows')
axes[1][1].legend()

# Plot for tumoral-only INDELs
maxim5 = np.max(data_non0[4])
moda5 = mode(data_non0[4])
axes[2][0].hist(data_non0[4], bins=range(int(min(data_non0[4])), int(max(data_non0[4])) + 2), color='blue', alpha=0.7)
# axes[2][0].set_xlim(left=0, right=100)
axes[2][0].scatter(moda5, 0, color="green", marker='o', label="Mode: {}".format(moda5))
axes[2][0].scatter(maxim5, 0, color="red", marker='o', label="Max: {}".format(maxim5))
axes[2][0].set_xlabel('Number of INDELs/window')
axes[2][0].set_ylabel('Count')
axes[2][0].set_title('Tumoral-only INDEL counts for 2kb windows')
axes[2][0].legend()

# Plot for tumoral-only SNVs
maxim6 = np.max(data_non0[5])
moda6 = mode(data_non0[5])
axes[2][1].hist(data_non0[5], bins=range(int(min(data_non0[5])), int(max(data_non0[5])) + 2), color='blue', alpha=0.7)
# axes[2][1].set_xlim(left=0, right=1500)
axes[2][1].scatter(moda6, 0, color="green", marker='o', label="Mode: {}".format(moda6))
axes[2][1].scatter(maxim6, 0, color="red", marker='o', label="Max: {}".format(maxim6))
axes[2][1].set_xlabel('Number of SNVs/window')
axes[2][1].set_ylabel('Count')
axes[2][1].set_title('Tumoral-only SNV counts for 2kb windows')
axes[2][1].legend()

# Plot for normal-only INDELs
maxim7 = np.max(data_non0[6])
moda7 = mode(data_non0[6])
axes[3][0].hist(data_non0[6], bins=range(int(min(data_non0[6])), int(max(data_non0[6])) + 2), color='blue', alpha=0.7)
# axes[3][0].set_xlim(left=0, right=100)
axes[3][0].scatter(moda7, 0, color="green", marker='o', label="Mode: {}".format(moda7))
axes[3][0].scatter(maxim7, 0, color="red", marker='o', label="Max: {}".format(maxim7))
axes[3][0].set_xlabel('Number of INDELs/window')
axes[3][0].set_ylabel('Count')
axes[3][0].set_title('Normal-only INDEL counts for 2kb windows')
axes[3][0].legend()

# Plot for normal-only SNVs
maxim8 = np.max(data_non0[7])
# moda8 = mode(data_non0[7])
axes[3][1].hist(data_non0[7], bins=range(int(min(data_non0[7])), int(max(data_non0[7])) + 2), color='blue', alpha=0.7)
# axes[3][1].set_xlim(left=0, right=1500)
# axes[3][1].scatter(moda8, 0, color="green", marker='o', label="Mode: {}".format(moda8))
axes[3][1].scatter(maxim8, 0, color="red", marker='o', label="Max: {}".format(maxim8))
axes[3][1].set_xlabel('Number of SNVs/window')
axes[3][1].set_ylabel('Count')
axes[3][1].set_title('Normal-only SNV counts for 2kb windows')
axes[3][1].legend()

# Adjust layout to prevent overlap
plt.tight_layout()

# Save the combined plot
plt.savefig('{}/hist_combined_data_non0.pdf'.format(out_dir))

# Create a single figure with subplots (plot log10 norm data)
fig, axes = plt.subplots(4, 2, figsize=(24, 24))

# Plot for potential germinal INDELs
maxim1 = np.max(normalized_data[0])
moda1 = mode(normalized_data[0])
axes[0][0].hist(normalized_data[0], bins=range(int(min(normalized_data[0])), int(max(normalized_data[0])) + 2), color='blue', alpha=0.7)
# axes[0][0].set_xlim(left=0, right=25)
axes[0][0].scatter(moda1, 0, color="green", marker='o', label="Mode: {}".format(moda1))
axes[0][0].scatter(maxim1, 0, color="red", marker='o', label="Max: {}".format(maxim1))
axes[0][0].set_xlabel('Number of INDELs/window')
axes[0][0].set_ylabel('Count')
axes[0][0].set_title('Potential germline INDEL counts for 2kb windows')
axes[0][0].legend()

# Plot for potential germinal SNVs
maxim2 = np.max(normalized_data[1])
moda2 = mode(normalized_data[1])
axes[0][1].hist(normalized_data[1], bins=range(int(min(normalized_data[1])), int(max(normalized_data[1])) + 2), color='blue', alpha=0.7)
# axes[0][1].set_xlim(left=0, right=1000)
axes[0][1].scatter(moda2, 0, color="green", marker='o', label="Mode: {}".format(moda2))
axes[0][1].scatter(maxim2, 0, color="red", marker='o', label="Max: {}".format(maxim2))
axes[0][1].set_xlabel('Number of SNVs/window')
axes[0][1].set_ylabel('Count')
axes[0][1].set_title('Potential germline SNV counts for 2kb windows')
axes[0][1].legend()

# Plot for single-read INDELs
maxim3 = np.max(normalized_data[2])
moda3 = mode(normalized_data[2])
axes[1][0].hist(normalized_data[2], bins=range(int(min(normalized_data[2])), int(max(normalized_data[2])) + 2), color='blue', alpha=0.7)
# axes[1][0].set_xlim(left=0, right=100)
axes[1][0].scatter(moda3, 0, color="green", marker='o', label="Mode: {}".format(moda3))
axes[1][0].scatter(maxim3, 0, color="red", marker='o', label="Max: {}".format(maxim3))
axes[1][0].set_xlabel('Number of INDELs/window')
axes[1][0].set_ylabel('Count')
axes[1][0].set_title('Single-read INDEL counts for 2kb windows')
axes[1][0].legend()

# Plot for single-read SNVs
maxim4 = np.max(normalized_data[3])
moda4 = mode(normalized_data[3])
axes[1][1].hist(normalized_data[3], bins=range(int(min(normalized_data[3])), int(max(normalized_data[3])) + 2), color='blue', alpha=0.7)
# axes[1][1].set_xlim(left=0, right=1500)
axes[1][1].scatter(moda4, 0, color="green", marker='o', label="Mode: {}".format(moda4))
axes[1][1].scatter(maxim4, 0, color="red", marker='o', label="Max: {}".format(maxim4))
axes[1][1].set_xlabel('Number of SNVs/window')
axes[1][1].set_ylabel('Count')
axes[1][1].set_title('Single-read SNV counts for 2kb windows')
axes[1][1].legend()

# Plot for tumoral-only INDELs
maxim5 = np.max(normalized_data[4])
moda5 = mode(normalized_data[4])
axes[2][0].hist(normalized_data[4], bins=range(int(min(normalized_data[4])), int(max(normalized_data[4])) + 2), color='blue', alpha=0.7)
# axes[2][0].set_xlim(left=0, right=100)
axes[2][0].scatter(moda5, 0, color="green", marker='o', label="Mode: {}".format(moda5))
axes[2][0].scatter(maxim5, 0, color="red", marker='o', label="Max: {}".format(maxim5))
axes[2][0].set_xlabel('Number of INDELs/window')
axes[2][0].set_ylabel('Count')
axes[2][0].set_title('Tumoral-only INDEL counts for 2kb windows')
axes[2][0].legend()

# Plot for tumoral-only SNVs
maxim6 = np.max(normalized_data[5])
moda6 = mode(normalized_data[5])
axes[2][1].hist(normalized_data[5], bins=range(int(min(normalized_data[5])), int(max(normalized_data[5])) + 2), color='blue', alpha=0.7)
# axes[2][1].set_xlim(left=0, right=1500)
axes[2][1].scatter(moda6, 0, color="green", marker='o', label="Mode: {}".format(moda6))
axes[2][1].scatter(maxim6, 0, color="red", marker='o', label="Max: {}".format(maxim6))
axes[2][1].set_xlabel('Number of SNVs/window')
axes[2][1].set_ylabel('Count')
axes[2][1].set_title('Tumoral-only SNV counts for 2kb windows')
axes[2][1].legend()

# Plot for normal-only INDELs
maxim7 = np.max(normalized_data[6])
moda7 = mode(normalized_data[6])
axes[3][0].hist(normalized_data[6], bins=range(int(min(normalized_data[6])), int(max(normalized_data[6])) + 2), color='blue', alpha=0.7)
# axes[3][0].set_xlim(left=0, right=100)
axes[3][0].scatter(moda7, 0, color="green", marker='o', label="Mode: {}".format(moda7))
axes[3][0].scatter(maxim7, 0, color="red", marker='o', label="Max: {}".format(maxim7))
axes[3][0].set_xlabel('Number of INDELs/window')
axes[3][0].set_ylabel('Count')
axes[3][0].set_title('Normal-only INDEL counts for 2kb windows')
axes[3][0].legend()

# Plot for normal-only SNVs
maxim8 = np.max(normalized_data[7])
# moda8 = mode(normalized_data[7]])
axes[3][1].hist(normalized_data[7], bins=range(int(min(normalized_data[7])), int(max(normalized_data[7])) + 2), color='blue', alpha=0.7)
# axes[3][1].set_xlim(left=0, right=1500)
# axes[3][1].scatter(moda8, 0, color="green", marker='o', label="Mode: {}".format(moda8))
axes[3][1].scatter(maxim8, 0, color="red", marker='o', label="Max: {}".format(maxim8))
axes[3][1].set_xlabel('Number of SNVs/window')
axes[3][1].set_ylabel('Count')
axes[3][1].set_title('Normal-only SNV counts for 2kb windows')
axes[3][1].legend()

# Adjust layout to prevent overlap
plt.tight_layout()

# Save the combined plot
plt.savefig('{}/hist_combined_data_log10norm.pdf'.format(out_dir))

# Make some stats
dataset = [(data1, "INDELs tumor-normal"),
           (data2, "SNVs tumor-normal"),
           (data3, "INDELs single reads"),
           (data4, "SNVs single reads"),
           (data5, "INDELs tumor-only"),
           (data6, "SNVs tumor-only"),
           (data7, "INDELs normal-only"),
           (data8, "SNVs normal-only")]

total_windows = 0
for chr in windows:
    total_windows += len(windows[chr])
#print(total_windows)

stats_file = "{}stats.txt".format(out_dir)
with open(stats_file, "w") as file:
    file.write('Number of windows: {}\n'.format(total_windows))
    file.write("\n")
    for data in dataset:
        counts = len(data[0])
        minim = np.min(data[0])
        maxim = np.max(data[0])
        mean = np.mean(data[0])
        median = np.median(data[0])
        try:
            moda = mode(data[0])
        except StatisticsError:
            moda = "Undefined"
        std = np.std(data[0])
        var = np.var(data[0])
        file.write('-------------{}-------------\n'.format(data[1]))
        file.write('Counts: {}\n'.format(counts))
        file.write('Min: {}\n'.format(minim))
        file.write('Max: {}\n'.format(maxim))
        file.write('Mean: {}\n'.format(mean))
        file.write('Median: {}\n'.format(median))
        file.write('Mode: {}\n'.format(moda))
        file.write('Standard deviation: {}\n'.format(std))
        file.write('Variance: {}\n'.format(var))

# Supporting reads stats

dataset_supporting = [(INDEL_supp_reads_count, "INDELs"),
                      (SNV_supp_reads_count, "SNVs")]

var_type = ["UNCLASSIFIED", "NORMAL_SINGLE_READ_VARIANT", "TUMORAL_SINGLE_READ_VARIANT", "NORMAL_ONLY_VARIANT",
            "TUMORAL_ONLY_VARIANT", "TUMORAL_NORMAL_VARIANT"]

supporting_stats_file = "{}supporting_reads_stats.txt".format(out_dir)
with open(supporting_stats_file, "w") as file:
    file.write('Supporting reads summary\n')
    file.write("\n")
    for data in dataset_supporting:
        file.write('--------------------{}--------------------\n'.format(data[1]))
        file.write("\n")
        for idx, i in enumerate(data[0]):
            counts = len(i)
            if len(i) == 0:  # Check if the list is empty
                minim = None
                maxim = None
                mean = None
                median = None
                moda = None
                std = None
                var = None
            else:
                minim = np.min(i)
                maxim = np.max(i)
                mean = np.mean(i)
                median = np.median(i)
                try:
                    moda = mode(i)
                except StatisticsError:
                    moda = "Undefined"
                std = np.std(i)
                var = np.var(i)
            # write the results
            file.write('-------------{}-------------\n'.format(var_type[idx]))
            file.write('Counts: {}\n'.format(counts))
            file.write('Min: {}\n'.format(minim))
            file.write('Max: {}\n'.format(maxim))
            file.write('Mean: {}\n'.format(mean))
            file.write('Median: {}\n'.format(median))
            file.write('Mode: {}\n'.format(moda))
            file.write('Standard deviation: {}\n'.format(std))
            file.write('Variance: {}\n'.format(var))
            file.write("\n")

# Window length analysis

window_length_qual = {"2,000":0, "2,000-100,000":0, "100,000-1M":0, ">1M":0}

for seq_name, seq_windows in windows.items():
    for window in seq_windows:
        length = window[1] - window[0]
        if int(length) <= 2000:
            window_length_qual["2,000"] += 1
        elif 2000 < int(length) <= 100000:
            window_length_qual["2,000-100,000"] += 1
        elif 100000 < int(length) <= 1000000:
            window_length_qual["100,000-1M"] += 1
        elif int(length) > 1000000:
            window_length_qual[">1M"] += 1
        #else:
            #print(length)

window_length_file = "{}window_length_qual.txt".format(out_dir)

with open(window_length_file, "w") as file:
    file.write('Windows length:\tCount:\n')
    for i in window_length_qual:
        line = str(i) + "\t" + str(window_length_qual[i]) + "\n"
        file.write(line)

# Record the end time
end_time = time.time()

# Calculate and print the total execution time
total_time_s = end_time - start_time
total_time_min = round(total_time_s / 60, 2)
total_time_h = round(total_time_s / 3600, 2)
print(f"Total execution time: {total_time_min} min.")
print(f"Total execution time: {total_time_h} h.")
