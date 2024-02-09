'''This script is meant to generate a text file with all the variants in SAM file alignment reads -opened with pysam-
found in a window of 2kb from a variant described in a VCF file -opened with VariantExtractor-'''

# import modules
import datetime
import time
import os
import re
import pysam
from variant_extractor import VariantExtractor
from functions import compare, anno_cigar, anno_cigar2, dict_generator, chr_converter  # local scripts
import matplotlib.pyplot as plt
import numpy as np

# Record the start time
start_time = time.time()

# Define input files
# in_N_bam = '/mnt/mn/mosaic_genome_PCAWG_1_N.bam'
# in_N_bai = '/mnt/mn/mosaic_genome_PCAWG_1_N.bam.bai'
# in_ref = '/mnt/mn/genome.fa'
in_ref = '/home/bscuser/Desktop/Obfuscation/SAM_files/mosaic_genome_PCAWG_1/genomic_data/genome.fa'
# in_T_bam = '/mnt/mn/mosaic_genome_PCAWG_1_T.bam'
in_T_bam = '/home/bscuser/Desktop/Obfuscation/SAM_files/mosaic_genome_PCAWG_1/genomic_data/mosaic_genome_PCAWG_1_T.bam'
in_T_bai = '/mnt/mn/mosaic_genome_PCAWG_1_T.bam.bai'
# in_vcf = '/mnt/mn/mosaic_genome_PCAWG_1.vcf'
in_vcf = '/home/bscuser/Desktop/Obfuscation/SAM_files/mosaic_genome_PCAWG_1/genomic_data/mosaic_genome_PCAWG_1.vcf'
in_100_vcf = '/mnt/mn/mosaic_genome_PCAWG_1_100.vcf'
in_25_vcf = '/mnt/mn/mosaic_genome_PCAWG_1_25.vcf'

# Define output files
try:
    out_dir = "/home/bscuser/Desktop/Obfuscation/SAM_files/mosaic_genome_PCAWG_1/output/"
    os.makedirs(out_dir)
except FileExistsError:
    print('Output files have been generated at an existing directory:')
    print(out_dir)

outN = "/home/bscuser/Desktop/Obfuscation/SAM_files/mosaic_genome_PCAWG_1/output/variation_N.txt"
outT = "/home/bscuser/Desktop/Obfuscation/SAM_files/mosaic_genome_PCAWG_1/output/variation_T.txt"
outNT = "/home/bscuser/Desktop/Obfuscation/SAM_files/mosaic_genome_PCAWG_1/output/variation_N-T.txt"

# Generate output files and define the result structure
out_files = [outN, outT, outNT]
for i in out_files:
    with open(i, 'w') as file:
        file.write("##date={}\n".format(datetime.datetime.now()))
        file.write("##script={}\n".format(os.path.basename(__file__)))
        file.write("##input_vcf={}\n".format(in_vcf))
        file.write("##input_sam_N={}\n".format(in_T_bam))
        file.write("##input_sam_T={}\n\n".format(in_T_bam))
        file.write("#CHR\tPOS\tREF\tALT\tTYPE\tVCF_TYPE\tVCF_POS\n")  # we define the header

# Open samfiles with pysam and vcf with VariantExtractor
# samfile_N = pysam.AlignmentFile(in_N_bam, 'rb')  # read a bam file
samfile_T = pysam.AlignmentFile(in_T_bam, 'rb')
extractor = VariantExtractor(in_vcf)

# Define some things related to cigar management
ref_consuming = ['M', 'D', 'N', '=', 'X']  # stores reference consuming cigar operations
regexp = r"(?<=[a-zA-Z=])(?=[0-9])|(?<=[0-9])(?=[a-zA-Z=])"  # regexp to split cigar string
symbol_add = ["X", "I", "D"]  # cigar operations to report
dict_cigar = {"X": "MM", "I": "INS", "D": "DEL"}
dict_tuples = {0: "M", 1: "I", 2: "D", 3: "N", 4: "S", 5: "H", 6: "P", 7: "=", 8: "X", 9: "B"}

# Create a dict to assign an integer value to each chr name

print("Begin dict.")
chr_dict = dict_generator(in_ref)
print("End dict.")

# Create a 2000 kb window for each variant

window_list = []  # create empty list

for variant in extractor:
    chr_name = chr_converter(chr_dict, variant.contig)  # we assign the integer value
    window_list.append((chr_name, variant.pos - 1000, variant.pos + 1000, variant.variant_type.name))

window_index = 0  # index to parse the window list
# window_vars_list = []  # list lo store the variations in each window
# window_vars = 0  # counter for the number of variations in each window
window_SNV_list = []  # list lo store the SNVs in each window
window_SNV = 0  # counter for the number of SNVs in each window
window_INDEL_list = []  # list lo store the INDELs in each window
window_INDEL = 0  # counter for the number of INDELs in each window

# Parse vcf and sam files in parallel

aln_iter = samfile_T.__iter__()
aln = next(aln_iter, None)  # if there is no next alignment, then aln == None

while aln is not None and window_index < len(window_list):

    if aln.cigarstring is None:
        aln = next(aln_iter, None)
        continue  # go back to while loop beginning

    # print('cigar:', aln.cigarstring)
    # print(aln.reference_name, aln.reference_start, aln.reference_end, window_list[window_index][0], window_list[window_index][1], window_list[window_index][2])
    chr_aln = chr_converter(chr_dict, aln.reference_name)  # we assign the integer value
    cmp = compare(chr_aln, aln.reference_start, aln.reference_end,
                  window_list[window_index][0], window_list[window_index][1], window_list[window_index][2])
    # print('cmp:', cmp)

    # if cmp is None:
    # continue  # go back to while loop beginning

    if -1 <= cmp <= 1:  # if intersection

        cigar = str(aln.cigarstring)  # we get the cigar string
        cigar_tuple = aln.cigartuples  # we get the cigar_tuple
        # ToDo: try with cigartuples?
        # print(cigar)
        md_tag = aln.get_tag("MD", with_value_type=False)  # we get the md tag
        # print(md_tag, cigar)

        ref_pos = aln.reference_start  # get the reference position where the cigar begins
        aln_chr = aln.reference_name  # we get the chr
        window = window_list[window_index]
        total_vars = anno_cigar2(outT, cigar, md_tag, regexp, ref_pos, aln_chr, symbol_add, dict_cigar, ref_consuming,
                                 window)
        # we annotate the variation and get the amount of variants in the aln
        # window_vars += total_vars
        window_SNV += total_vars[0]
        window_INDEL += total_vars[1]
        # print("total_vars:", total_vars)
        # print("window_SNVs:", window_SNV)
        # print("window_INDELs:", window_INDEL)
        # print("window_SNV_list:", window_SNV_list)
        # print("window_INDEL_list:", window_INDEL_list)
        aln = next(aln_iter, None)  # change to next alignment

    if cmp < -1:  # cmp == -2 or -3
        aln = next(aln_iter, None)  # change to next alignment

    if cmp > 1:  # cmp == 2 or 3
        # if window_index < len(window_list)-1:  # if not the last window
        window_index += 1  # change to next window maintaining the alignment
        # window_vars_list.append(window_vars)  # store the variation amount in the list
        # window_vars = 0  # restore the counter
        window_SNV_list.append(window_SNV)  # store the variation amount in the list
        window_SNV = 0  # restore the counter
        window_INDEL_list.append(window_INDEL)  # store the variation amount in the list
        window_INDEL = 0  # restore the counter
        # else:  # if there are no more windows and no intersection
        #   break

# Trial prints
# print(window_vars_list)
# print("Number of vcf variants checked:", len(window_vars_list))

# Plot the hist
# plt.hist(window_vars_list, bins=range(min(window_vars_list), max(window_vars_list) + 1), edgecolor='black')
# Add labels and title
# plt.xlabel('Number of variations/window')
# plt.ylabel('Frequency')
# plt.title('Distribution of the variation among the 2kb windows')
# Save the histogram as a PDF file
# plt.savefig('{}histogram_plot.pdf'.format("out_dir"))
# Display the histogram
# plt.show()

# Create a figure with two subplots
fig, (ax1, ax2) = plt.subplots(nrows=2, ncols=1, sharex=False, figsize=(8, 6))

# Create the histogram for the first subplot
ax1.hist(window_SNV_list, bins=range(min(window_SNV_list), max(window_SNV_list) + 1),
         alpha=0.5, color='blue', label='SNVs')
ax1.set_ylabel('Frequency')
ax1.set_title('Distribution of the SNVs among the 2kb windows')

# Create the histograms for the second subplot
ax2.hist(window_INDEL_list, bins=range(min(window_INDEL_list), max(window_INDEL_list) + 1),
         alpha=0.5, color='orange', label='INDELs')
ax2.set_xlabel('Number of variations/window')
ax2.set_ylabel('Frequency')
ax2.set_title('Distribution of the INDELs among the 2kb windows')

# Add a common title for the entire figure
# fig.suptitle('Two Histograms')

# Add legend
ax1.legend()
ax2.legend()

# Adjust layout to prevent overlapping titles and labels
plt.tight_layout()

# Save the histogram as a PDF file
plt.savefig('{}hist/hist_var_types.pdf'.format(out_dir))

# Show the plot
# plt.show()

# Statistics

SNV_counts = len(window_SNV_list)  # should be same as len(window_list) ?
SNV_min = np.min(window_SNV_list)
SNV_max = np.max(window_SNV_list)
SNV_mean = np.mean(window_SNV_list)
SNV_median = np.median(window_SNV_list)
SNV_std = np.std(window_SNV_list)
SNV_var = np.var(window_SNV_list)

INDEL_counts = len(window_INDEL_list)  # should be same as len(window_list) ?
INDEL_min = np.min(window_INDEL_list)
INDEL_max = np.max(window_INDEL_list)
INDEL_mean = np.mean(window_INDEL_list)
INDEL_median = np.median(window_INDEL_list)
INDEL_std = np.std(window_INDEL_list)
INDEL_var = np.var(window_INDEL_list)

stats_file = '/home/bscuser/Desktop/Obfuscation/SAM_files/mosaic_genome_PCAWG_1/output/stats.txt'
with open(stats_file, "w") as file:
    file.write('Number of windows: {}\n'.format(len(window_list)))
    file.write("\n")
    file.write('-------------SNVs-------------\n')
    file.write('Counts: {}\n'.format(SNV_counts))
    file.write('Min: {}\n'.format(SNV_min))
    file.write('Max: {}\n'.format(SNV_max))
    file.write('Mean: {}\n'.format(SNV_mean))
    file.write('Median: {}\n'.format(SNV_median))
    file.write('Standard deviation: {}\n'.format(SNV_std))
    file.write('Variance: {}\n'.format(SNV_var))
    file.write('-------------INDELs-------------\n')
    file.write('Counts: {}\n'.format(INDEL_counts))
    file.write('Min: {}\n'.format(INDEL_min))
    file.write('Max: {}\n'.format(INDEL_max))
    file.write('Mean: {}\n'.format(INDEL_mean))
    file.write('Median: {}\n'.format(INDEL_median))
    file.write('Standard deviation: {}\n'.format(INDEL_std))
    file.write('Variance: {}\n'.format(INDEL_var))

# Record the end time
end_time = time.time()

# Calculate and print the total execution time
total_time_s = end_time - start_time
total_time_h = round(total_time_s / 3600, 2)
# print(f"Total execution time: {total_time_s} seconds.")
print(f"Total execution time: {total_time_h} h.")