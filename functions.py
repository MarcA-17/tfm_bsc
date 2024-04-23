"""This script contains all the functions created for the variant_detection script"""

import re
import sys
import pysam
from variant_extractor import VariantExtractor

in_N_bam = sys.argv[1]
in_T_bam = sys.argv[2]
in_ref = sys.argv[3]
in_vcf = sys.argv[4]

samfile_N = pysam.AlignmentFile(in_N_bam, 'rb')  # read a bam file
samfile_T = pysam.AlignmentFile(in_T_bam, 'rb')
extractor = VariantExtractor(in_vcf)
ref_genome = pysam.Fastafile(in_ref)

TUMORAL_SINGLE_READ_VARIANT = 0
NORMAL_SINGLE_READ_VARIANT = 1
TUMORAL_ONLY_VARIANT = 2
NORMAL_ONLY_VARIANT = 3
TUMORAL_NORMAL_VARIANT = 4


def intersect(alignment, variant):
    """This function is meant to assess whether an alignment from a samfile -opened with pysam- intersects with a
      variant record from a vcf file -opened with VariantExtractor-"""

    # Check if they are on the same chromosome
    if alignment.reference_name == variant.contig:
        # Check if all relevant attributes are integers
        if all(isinstance(attr, int) for attr in [alignment.reference_start, variant.pos, alignment.reference_end]):
            # Check for intersection
            if alignment.reference_start <= variant.pos <= alignment.reference_end:
                return True

    return False


def intersect_window(alignment, window):
    """This function is meant to assess whether an alignment from a samfile -opened with pysam- intersects with an
    alignment window"""

    # Check if they are on the same chromosome
    if alignment.reference_name == window[0]:
        # Check if all relevant attributes are integers
        if all(isinstance(attr, int) for attr in [alignment.reference_start, alignment.reference_end,
                                                  window[1], window[2]]):
            # Check for intersection
            if alignment.reference_start >= window[1] or alignment.reference_end <= window[2]:
                return True

    return False


def compare(chr_al, init_al, end_al, chr_wd, init_wd, end_wd):
    """This function compares an alignment and a window and assesses whether the alignment falls before,
    within or after the window"""

    # Check if there is any relevant attribute that is not an integer
    # if any(not isinstance(attr, int) for attr in [init_al, init_wd, end_al, end_wd]):
    # print("At least one attribute is not an integer.")
    # print(init_al, init_wd, end_al, end_wd)  # tmp
    # return None

    # Define the overlap condition
    overlap = init_wd <= end_al and end_wd >= init_al

    if chr_al < chr_wd:
        return -3
    elif chr_al > chr_wd:
        return 3

    # For the remaining cases chr_al == chr_wd
    if end_al < end_wd:
        return -1 if overlap else -2
    elif end_wd < end_al:
        return 1 if overlap else 2

    # For the remaining cases end_al == end_wd
    if init_al < init_wd:
        return -1
    elif init_wd < init_al:
        return 1

    # Only remaining case: window and alignment share start and end
    return 0


def anno_cigar(out_file, cigar, regexp, ref_pos, aln_chr, symbol_add, dict_cigar, ref_consuming, window):
    """This function takes a cigar string and annotates all the changes that it contains in an output file, while
    it returns the amount of changes annotated. It returns the number of total variations within the window.
    regexp: splits the cigar string
    symbol_add: cigar operations to report
    dict_cigar: dictionary to translate cigar letters
    ref_consuming: stores reference consuming cigar operations"""

    with open(out_file, "a") as fileT:
        cigar_list = re.split(regexp, cigar)  # we split the string into a list
        # we parse the cigar list
        cigar_pos = 0  # to parse cigar_list
        cigar_length = 0  # to track the position in the genome
        total_vars = 0  # counter to store the variations found in this tumor_aln
        for i in cigar_list:
            if i.isdigit():
                if cigar_list[cigar_pos + 1] in symbol_add:
                    true_pos = ref_pos + cigar_length
                    # we report the variation
                    fileT.write("{}\t".format(aln_chr))  # we write the chr
                    fileT.write("{}\t".format(true_pos))  # we write the pos
                    fileT.write("-\t")  # we write the ref
                    fileT.write("-\t")  # we write the alt
                    fileT.write("{}\t".format(dict_cigar[cigar_list[cigar_pos + 1]]))  # we write the type
                    fileT.write("{}\t".format(window[3]))  # we write the VCF_TYPE
                    fileT.write("{}\n".format(int(window[1]) + 1000))  # we write the VCF_POS
                    total_vars += 1  # store the variation

                if cigar_list[cigar_pos + 1] in ref_consuming:  # if cigar operation consumes reference
                    cigar_length += int(i)
                cigar_pos += 1
            else:
                cigar_pos += 1

    return total_vars


def encode_var(pos_init, pos_end, var_type, length):
    "This function returns a unique string for each indel to be its key in a dict"

    unique_string = "{};{};{};{}".format(pos_init, pos_end, var_type, length)
    return unique_string


def dict_generator(fasta_file):
    """This function parses a fasta file (ideally the human reference genome) and returns a dictionary
        whose keys are the chr names and whose values are an integer representing their ordinal number"""

    chr_dict = {}  # create an empty dict
    # ordinal = 1  # initialize the ordinal sequence
    ref_sequences = fasta_file.references
    chr_dict = {k: v for (k, v) in zip(ref_sequences, range(len(ref_sequences)))}
    # with open(fasta, "r") as file:
    #    for line in file:
    #        if line.startswith(">"):
    #            key = line[1:].split()[0]  # we get the chr name
    #            chr_dict[key] = ordinal
    #            ordinal += 1
    return chr_dict


def chr_converter(chr_dict, chr_name):
    """This function translates the chr name from an alignment file to its integer value using the dictionary
    that has been created wit the dict_generator function."""

    chr_integer = chr_dict[chr_name]
    return chr_integer


chr_dict = dict_generator(ref_genome)


def get_windows(variants, ref_genome, window_size=2000):
    half_window = int(window_size / 2)
    windows = dict()
    for seq in ref_genome.references:
        # seq_converted = chr_converter(chr_dict, seq)
        windows[seq] = []
    for variant in variants:
        windows_in_seq = windows[variant.contig]
        if variant.variant_type == "INV":  # VariantType.INV
            if variant.pos + half_window > variant.end - half_window:
                windows_in_seq.append((variant.pos - half_window, variant.end + half_window, variant))
            else:
                windows_in_seq.append((variant.pos - half_window, variant.pos + half_window, variant))
                windows_in_seq.append((variant.end - half_window, variant.end + half_window, variant))
        elif variant.variant_type == "TRA":  # VariantType.TRA
            windows_in_seq.append((variant.pos - half_window, variant.pos + half_window, variant))
            windows_in_seq.append((variant.end - half_window, variant.end + half_window, variant))
        elif variant.variant_type == "SNV":  # VariantType.SNV
            windows_in_seq.append((variant.pos - half_window, variant.pos + half_window, variant))
        else:
            windows_in_seq.append((variant.pos - half_window, variant.end + half_window, variant))
    for seq_name, seq_windows in windows.items():
        seq_windows.sort(key=lambda x: (x[0], x[1]))
    return windows
