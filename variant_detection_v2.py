'''This script is meant to generate a text file with all the variants in SAM file alignment reads -opened with pysam-
found in a window of 2kb from a variant described in a VCF file -opened with VariantExtractor-'''

# import modules
import datetime
import sys
import time
import os
import re
import pysam
from variant_extractor import VariantExtractor
from functions import compare, anno_cigar2, dict_generator, chr_converter  # local scripts
import matplotlib.pyplot as plt
import numpy as np
from statistics import mode

# Record the start time
start_time = time.time()

# Define input files
# CLI Inputs #
in_N_bam = sys.argv[1]
in_T_bam = sys.argv[2]
in_ref = sys.argv[3]
in_vcf = sys.argv[4]
"""
in_N_bam = '/mnt/mn/mosaic_genome_PCAWG_1_N.bam'
in_N_bai = '/mnt/mn/mosaic_genome_PCAWG_1_N.bam.bai'
# in_ref = '/mnt/mn/genome.fa'
in_ref = '/home/bscuser/Desktop/Obfuscation/SAM_files/mosaic_genome_PCAWG_1/genomic_data/genome.fa'
# in_T_bam = '/mnt/mn/mosaic_genome_PCAWG_1_T.bam'
in_T_bam = '/home/bscuser/Desktop/Obfuscation/SAM_files/mosaic_genome_PCAWG_1/genomic_data/mosaic_genome_PCAWG_1_T.bam'
in_T_bai = '/mnt/mn/mosaic_genome_PCAWG_1_T.bam.bai'
# in_vcf = '/mnt/mn/mosaic_genome_PCAWG_1.vcf'
in_vcf = '/home/bscuser/Desktop/Obfuscation/SAM_files/mosaic_genome_PCAWG_1/genomic_data/mosaic_genome_PCAWG_1.vcf'
in_100_vcf = '/mnt/mn/mosaic_genome_PCAWG_1_100.vcf'
in_25_vcf = '/mnt/mn/mosaic_genome_PCAWG_1_25.vcf'
"""

# Define output files
try:
    # out_dir = "/home/bscuser/Desktop/Obfuscation/SAM_files/mosaic_genome_PCAWG_1/output/"
    # out_dir = "/run/media/nicolasggg/TOSHIBA EXT/BSC_CG_group/Obfuscation/testing"
    out_dir = sys.argv[5]
    os.makedirs(out_dir)
except FileExistsError:
    print('Output files have been generated at an existing directory:')
    print(out_dir)
    print('{}/hist_combined.pdf'.format(out_dir))

#outN = f'{out_dir}variation_N.txt'
#outT = f'{out_dir}variation_T.txt'
#outNT = f'{out_dir}variation_N-T.txt'
# outN = "/home/bscuser/Desktop/Obfuscation/SAM_files/mosaic_genome_PCAWG_1/output/variation_N.txt"
# outT = "/home/bscuser/Desktop/Obfuscation/SAM_files/mosaic_genome_PCAWG_1/output/variation_T.txt"
# outNT = "/home/bscuser/Desktop/Obfuscation/SAM_files/mosaic_genome_PCAWG_1/output/variation_N-T.txt"

# Generate output files and define the result structure
# out_files = [outN, outT, outNT]
#for i in out_files:
 #   with open(i, 'w') as file:
  #      file.write("##date={}\n".format(datetime.datetime.now()))
   #     file.write("##script={}\n".format(os.path.basename(__file__)))
    #    file.write("##input_vcf={}\n".format(in_vcf))
     #   file.write("##input_sam_N={}\n".format(in_T_bam))
      #  file.write("##input_sam_T={}\n\n".format(in_T_bam))
       # file.write("#CHR\tPOS\tREF\tALT\tTYPE\tVCF_TYPE\tVCF_POS\n")  # we define the header

# Open samfiles with pysam and vcf with VariantExtractor
samfile_N = pysam.AlignmentFile(in_N_bam, 'rb')  # read a bam file
samfile_T = pysam.AlignmentFile(in_T_bam, 'rb')
extractor = VariantExtractor(in_vcf)
ref_genome = pysam.Fastafile(in_ref)


# Codes representing the status of variations according to the read and tumor normal sample context
TUMORAL_SINGLE_READ_VARIANT = 0
NORMAL_SINGLE_READ_VARIANT = 1
TUMORAL_ONLY_VARIANT = 2
NORMAL_ONLY_VARIANT = 3
TUMORAL_NORMAL_VARIANT = 4

# Create a dict to assign an integer value to each chr name

print("Begin dict.")
chr_dict = dict_generator(ref_genome)
print("End dict.")

# Create a 2000 kb window for each variant

window_list = []  # create empty list

for variant in extractor:
    chr_name = chr_converter(chr_dict, variant.contig)  # we assign the integer value
    window_list.append((chr_name, variant.pos - 1000, variant.pos + 1000, variant.variant_type.name))

# window_index = 0  # index to parse the window list
INDEL_counts = ([], [], [], [], [])  # we init a list to store the INDELs for each window
SNV_counts = ([], [], [], [], [])  # we init a list to store the SNVs for each window
# window_vars_list = []  # list lo store the variations in each window
# window_vars = 0  # counter for the number of variations in each window
# window_SNV_list = []  # list lo store the SNVs in each window
# window_SNV = 0  # counter for the number of SNVs in each window
# window_INDEL_list = []  # list lo store the INDELs in each window
# window_INDEL = 0  # counter for the number of INDELs in each window

# Parse vcf and sam files in parallel

tumor_aln_iter = samfile_T.__iter__()
normal_aln_iter = samfile_N.__iter__()
# tumor_aln = next(tumor_aln_iter, None)  # if there is no next alignment, then tumor_aln == None
# normal_aln = next(normal_aln_iter, None)
dict_indel_count = {}  # initialize dict to store INDELs count
dict_snv_count = {}  # init dict to store SNVs count
hold_iter = [False, False]  # to change the window while maintaining the aln
holded_alns = [0, 0]  # init the list

for window in window_list:
    for iterIdx, current_iter in enumerate((tumor_aln_iter, normal_aln_iter)):
        if hold_iter[iterIdx]:
            current_aln = holded_alns[iterIdx]
            hold_iter[iterIdx] = False
        else:
            current_aln = next(current_iter, None)
        while current_aln is not None:
            if current_aln.cigarstring is None:
                current_aln = next(current_iter, None)
                continue  # go back to while loop beginning
            # print('cigar:', tumor_aln.cigarstring)
            # print(current_aln.reference_name, current_aln.reference_start, current_aln.reference_end,
                  # window[0], window[1], window[2])
            chr_aln = chr_converter(chr_dict, current_aln.reference_name)  # we assign the integer value
            cmp = compare(chr_aln, current_aln.reference_start, current_aln.reference_end,
                          window[0], window[1], window[2])
            # print('cmp:', cmp)
            # if cmp is None:
            # continue  # go back to while loop beginning
            if -1 <= cmp <= 1:  # if intersection
                cigar = str(current_aln.cigarstring)  # we get the cigar string
                # cigar_tuple = tumor_aln.cigartuples  # we get the cigar_tuple
                # print(cigar)
                # we get the md tag and process it
                md_tag = current_aln.get_tag("MD", with_value_type=False)
                pattern_md = r'0|\^[A-Z]+|[A-Z]|[0-9]+'
                md_list = re.findall(pattern_md, md_tag)
                # print(md_list)
                ref_pos = current_aln.reference_start  # get the reference position where the cigar begins
                aln_chr = current_aln.reference_name  # we get the chr
                # window = window_list[window_index]
                # total_vars = anno_cigar2(dict_indel_count, outT, cigar, md_list, regexp, ref_pos, aln_chr, symbol_add, dict_cigar, ref_consuming, window)
                # anno_cigar2(dict_indel_count, dict_snv_count, outT, cigar, md_list, regexp, ref_pos, aln_chr,
                  #          symbol_add, dict_cigar,
                   #         ref_consuming, window)
                anno_cigar2(dict_indel_count, dict_snv_count, cigar, md_list, ref_pos, aln_chr, iterIdx, ref_genome, current_aln)
                current_aln = next(current_iter, None)  # change to next alignment
                # we annotate the variation and get the amount of variants in the tumor_aln
                # window_vars += total_vars
                # window_SNV += total_vars[0]
                # window_INDEL += total_vars[1]
                # print("total_vars:", total_vars)
                # print("window_SNVs:", window_SNV)
                # print("window_INDELs:", window_INDEL)
                # print("window_SNV_list:", window_SNV_list)
                # print("window_INDEL_list:", window_INDEL_list)
            if cmp < -1:  # cmp == -2 or -3
                current_aln = next(current_iter, None)  # change to next alignment
            if cmp > 1:  # cmp == 2 or 3
                # if window_index < len(window_list)-1:  # if not the last window
                if iterIdx == 1:
                    single_read_variant_counts = 0
                    tumoral_only_counts = 0
                    tumoral_normal_variant_counts = 0
                    normal_only_variants_counts = 0
                    for k, v in dict_indel_count.items():
                        if v == NORMAL_SINGLE_READ_VARIANT or v == TUMORAL_SINGLE_READ_VARIANT:
                            single_read_variant_counts += 1
                        if v == TUMORAL_ONLY_VARIANT:
                            tumoral_only_counts += 1
                        if v == TUMORAL_NORMAL_VARIANT:
                            tumoral_normal_variant_counts += 1
                        if v == NORMAL_ONLY_VARIANT:
                            normal_only_variants_counts += 1
                    INDEL_counts[TUMORAL_SINGLE_READ_VARIANT].append(single_read_variant_counts)
                    INDEL_counts[TUMORAL_ONLY_VARIANT].append(tumoral_only_counts)
                    INDEL_counts[TUMORAL_NORMAL_VARIANT].append(tumoral_normal_variant_counts)  # we store the number of variations in the window
                    INDEL_counts[NORMAL_ONLY_VARIANT].append(normal_only_variants_counts)
                    dict_indel_count = {}  # re-init dict
                    single_read_variant_counts = 0
                    tumoral_only_counts = 0
                    tumoral_normal_variant_counts = 0
                    normal_only_variants_counts = 0
                    for k, v in dict_snv_count.items():
                        if v == NORMAL_SINGLE_READ_VARIANT or v == TUMORAL_SINGLE_READ_VARIANT:
                            single_read_variant_counts += 1
                        if v == TUMORAL_ONLY_VARIANT:
                            tumoral_only_counts += 1
                        if v == TUMORAL_NORMAL_VARIANT:
                            tumoral_normal_variant_counts += 1
                        if v == NORMAL_ONLY_VARIANT:
                            normal_only_variants_counts += 1
                    SNV_counts[TUMORAL_SINGLE_READ_VARIANT].append(single_read_variant_counts)
                    SNV_counts[TUMORAL_ONLY_VARIANT].append(tumoral_only_counts)
                    SNV_counts[TUMORAL_NORMAL_VARIANT].append(tumoral_normal_variant_counts)
                    SNV_counts[NORMAL_ONLY_VARIANT].append(normal_only_variants_counts)
                    # SNV_counts.append(len(dict_snv_count))  # we store the number of variations in the window
                    dict_snv_count = {}  # re-init dict
                hold_iter[iterIdx] = True
                holded_alns[iterIdx] = current_aln
                break
                # window_index += 1  # change to next window maintaining the alignment
                # print(SNV_counts)
                # window_vars_list.append(window_vars)  # store the variation amount in the list
                # window_vars = 0  # restore the counter
                # window_SNV_list.append(window_SNV)  # store the variation amount in the list
                # window_SNV = 0  # restore the counter
                # window_INDEL_list.append(window_INDEL)  # store the variation amount in the list
                # window_INDEL = 0  # restore the counter
                # else:  # if there are no more windows and no intersection
                #   break
    # if hold_iter[0] and hold_iter[1]:
      #  continue

""" 
    while tumor_aln is not None and window_index < len(window_list):

    if tumor_aln.cigarstring is None:
        tumor_aln = next(tumor_aln_iter, None)
        continue  # go back to while loop beginning

    # print('cigar:', tumor_aln.cigarstring)
    # print(tumor_aln.reference_name, tumor_aln.reference_start, tumor_aln.reference_end, window_list[window_index][0], window_list[window_index][1], window_list[window_index][2])
    chr_aln = chr_converter(chr_dict, tumor_aln.reference_name)  # we assign the integer value
    cmp = compare(chr_aln, tumor_aln.reference_start, tumor_aln.reference_end,
                  window_list[window_index][0], window_list[window_index][1], window_list[window_index][2])
    # print('cmp:', cmp)

    # if cmp is None:
    # continue  # go back to while loop beginning

    if -1 <= cmp <= 1:  # if intersection

        cigar = str(tumor_aln.cigarstring)  # we get the cigar string
        #cigar_tuple = tumor_aln.cigartuples  # we get the cigar_tuple
        # print(cigar)
        # we get the md tag and process it
        md_tag = tumor_aln.get_tag("MD", with_value_type=False)
        pattern_md = r'0|\^[A-Z]+|[A-Z]|[0-9]+'
        md_list = re.findall(pattern_md, md_tag)
        # print(md_list)

        ref_pos = tumor_aln.reference_start  # get the reference position where the cigar begins
        aln_chr = tumor_aln.reference_name  # we get the chr
        window = window_list[window_index]
        #total_vars = anno_cigar2(dict_indel_count, outT, cigar, md_list, regexp, ref_pos, aln_chr, symbol_add, dict_cigar, ref_consuming, window)
        anno_cigar2(dict_indel_count, dict_snv_count, outT, cigar, md_list, regexp, ref_pos, aln_chr, symbol_add, dict_cigar,
                    ref_consuming, window)
        # we annotate the variation and get the amount of variants in the tumor_aln
        # window_vars += total_vars
        #window_SNV += total_vars[0]
        #window_INDEL += total_vars[1]
        # print("total_vars:", total_vars)
        # print("window_SNVs:", window_SNV)
        # print("window_INDELs:", window_INDEL)
        # print("window_SNV_list:", window_SNV_list)
        # print("window_INDEL_list:", window_INDEL_list)
        tumor_aln = next(tumor_aln_iter, None)  # change to next alignment

    if cmp < -1:  # cmp == -2 or -3
        tumor_aln = next(tumor_aln_iter, None)  # change to next alignment

    if cmp > 1:  # cmp == 2 or 3
        # if window_index < len(window_list)-1:  # if not the last window
        INDEL_counts.append(len(dict_indel_count))  # we store the number of variations in the window
        dict_indel_count = {}  # re-init dict
        SNV_counts.append(len(dict_snv_count))  # we store the number of variations in the window
        dict_snv_count = {}  # re-init dict
        window_index += 1  # change to next window maintaining the alignment
        #print(SNV_counts)
        # window_vars_list.append(window_vars)  # store the variation amount in the list
        # window_vars = 0  # restore the counter
        #window_SNV_list.append(window_SNV)  # store the variation amount in the list
        #window_SNV = 0  # restore the counter
        #window_INDEL_list.append(window_INDEL)  # store the variation amount in the list
        #window_INDEL = 0  # restore the counter
        # else:  # if there are no more windows and no intersection
        #   break
"""
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

# ToDO: añadir en cada plot la moda y el max

# Create a single figure with two subplots
fig, axes = plt.subplots(4, 2, figsize=(24, 8))
data1 = INDEL_counts[TUMORAL_NORMAL_VARIANT]
data2 = SNV_counts[TUMORAL_NORMAL_VARIANT]
data3 = INDEL_counts[TUMORAL_SINGLE_READ_VARIANT]  # For now, it both tumoral and normal single reads
data4 = SNV_counts[TUMORAL_SINGLE_READ_VARIANT]  # For now, it both tumoral and normal single reads
data5 = INDEL_counts[TUMORAL_ONLY_VARIANT]
data6 = SNV_counts[TUMORAL_ONLY_VARIANT]
data7 = INDEL_counts[NORMAL_ONLY_VARIANT]
data8 = SNV_counts[NORMAL_ONLY_VARIANT]

# Plot for potential germinal INDELs
axes[0][0].hist(data1, bins=range(int(min(data1)), int(max(data1)) + 2), color='blue', alpha=0.7)
#axes[0][0].set_xlim(left=0, right=25)
axes[0][0].set_xlabel('Number of INDELs/window')
axes[0][0].set_ylabel('Count')
axes[0][0].set_title('Potential germline INDEL counts for 2kb windows')

# Plot for potential germinal SNVs
axes[0][1].hist(data2, bins=range(int(min(data2)), int(max(data2)) + 2), color='blue', alpha=0.7)
#axes[0][1].set_xlim(left=0, right=1000)
axes[0][1].set_xlabel('Number of SNVs/window')
axes[0][1].set_ylabel('Count')
axes[0][1].set_title('Potential germline SNV counts for 2kb windows')

# Plot for single-read INDELs
axes[1][0].hist(data3, bins=range(int(min(data3)), int(max(data3)) + 2), color='blue', alpha=0.7)
#axes[1][0].set_xlim(left=0, right=100)
axes[1][0].set_xlabel('Number of INDELs/window')
axes[1][0].set_ylabel('Count')
axes[1][0].set_title('Single-read INDEL counts for 2kb windows')

# Plot for single-read SNVs
axes[1][1].hist(data4, bins=range(int(min(data4)), int(max(data4)) + 2), color='blue', alpha=0.7)
#axes[1][1].set_xlim(left=0, right=1500)
axes[1][1].set_xlabel('Number of SNVs/window')
axes[1][1].set_ylabel('Count')
axes[1][1].set_title('Single-read SNV counts for 2kb windows')

# Plot for tumoral-only INDELs
axes[2][0].hist(data5, bins=range(int(min(data5)), int(max(data5)) + 2), color='blue', alpha=0.7)
#axes[2][0].set_xlim(left=0, right=100)
axes[2][0].set_xlabel('Number of INDELs/window')
axes[2][0].set_ylabel('Count')
axes[2][0].set_title('Tumoral-only INDEL counts for 2kb windows')

# Plot for tumoral-only SNVs
axes[2][1].hist(data6, bins=range(int(min(data6)), int(max(data6)) + 2), color='blue', alpha=0.7)
#axes[2][1].set_xlim(left=0, right=1500)
axes[2][1].set_xlabel('Number of SNVs/window')
axes[2][1].set_ylabel('Count')
axes[2][1].set_title('Tumoral-only SNV counts for 2kb windows')

# Plot for normal-only INDELs
axes[3][0].hist(data7, bins=range(int(min(data7)), int(max(data7)) + 2), color='blue', alpha=0.7)
#axes[3][0].set_xlim(left=0, right=100)
axes[3][0].set_xlabel('Number of INDELs/window')
axes[3][0].set_ylabel('Count')
axes[3][0].set_title('Normal-only INDEL counts for 2kb windows')

# Plot for normal-only SNVs
axes[3][1].hist(data8, bins=range(int(min(data8)), int(max(data8)) + 2), color='blue', alpha=0.7)
#axes[3][1].set_xlim(left=0, right=1500)
axes[3][1].set_xlabel('Number of SNVs/window')
axes[3][1].set_ylabel('Count')
axes[3][1].set_title('Normal-only SNV counts for 2kb windows')

# Adjust layout to prevent overlap
plt.tight_layout()

# Save the combined plot
plt.savefig('{}/hist_combined.pdf'.format(out_dir))

"""# Create a figure with two subplots
fig, (ax1, ax2) = plt.subplots(nrows=2, ncols=1, sharex=False, figsize=(8, 6))

# Create the histogram for the first subplot
ax1.hist(window_SNV_list, bins=range(min(window_SNV_list), max(window_SNV_list) + 1),
         alpha=0.5, color='blue', label='SNVs')
ax1.set_ylabel('Frequency')
ax1.set_title('Distribution of the SNVs among the 2kb windows')

# Set x-axis limit for the first subplot
ax1.set_xlim(right=2500)
# ax1.set_ylim()

# Create the histograms for the second subplot
ax2.hist(window_INDEL_list, bins=range(min(window_INDEL_list), max(window_INDEL_list) + 1),
         alpha=0.5, color='orange', label='INDELs')
ax2.set_xlabel('Number of variations/window')
ax2.set_ylabel('Frequency')
ax2.set_title('Distribution of the INDELs among the 2kb windows')

# Set x-axis limit for the second subplot
ax2.set_xlim(right=250)

# Add a common title for the entire figure
# fig.suptitle('Two Histograms')

# Add legend
ax1.legend()
ax2.legend()

# Adjust layout to prevent overlapping titles and labels
plt.tight_layout()

# Save the histogram as a PDF file
plt.savefig('{}hist/hist_SNVs.pdf'.format(out_dir))

# Show the plot
# plt.show()"""

# Make some stats
dataset = [(data1, "INDELs tumor-normal"),
           (data2, "SNVs tumor-normal"),
           (data3, "INDELs single reads"),
           (data4, "SNVs single reads"),
           (data5, "INDELs tumor-only"),
           (data6, "SNVs tumor-only"),
           (data7, "INDELs normal-only"),
           (data8, "SNVs normal-only")]

stats_file = "{}stats.txt".format(out_dir)
with open(stats_file, "w") as file:
    for data in dataset:
        counts = len(data[0])
        minim = np.min(data[0])
        maxim = np.max(data[0])
        mean = np.mean(data[0])
        median = np.median(data[0])
        moda = mode(data[0])
        std = np.std(data[0])
        var = np.var(data[0])

    file.write('Number of windows: {}\n'.format(len(window_list)))
    file.write("\n")
    file.write('-------------{}-------------\n'.format(data[1]))
    file.write('Counts: {}\n'.format(counts))
    file.write('Min: {}\n'.format(minim))
    file.write('Max: {}\n'.format(maxim))
    file.write('Mean: {}\n'.format(mean))
    file.write('Median: {}\n'.format(median))
    file.write('Mode: {}\n'.format(moda))
    file.write('Standard deviation: {}\n'.format(std))
    file.write('Variance: {}\n'.format(var))


"""# Statistics

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
    file.write('Variance: {}\n'.format(INDEL_var))"""

# Record the end time
end_time = time.time()

# Calculate and print the total execution time
total_time_s = end_time - start_time
total_time_h = round(total_time_s / 3600, 2)
# print(f"Total execution time: {total_time_s} seconds.")
print(f"Total execution time: {total_time_h} h.")
