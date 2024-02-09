"""This script contains all the functions created for the variant_detection script"""


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


import re


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
        total_vars = 0  # counter to store the variations found in this aln
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


def anno_cigar2(out_file, cigar, md_tag, regexp, ref_pos, aln_chr, symbol_add, dict_cigar, ref_consuming, window):
    """This function takes a cigar string and annotates all the changes that it contains in an output file, while
    it returns the amount of changes annotated. It returns the number of potential SNVs and INDELs within the window.
    regexp: splits the cigar string
    symbol_add: cigar operations to report
    dict_cigar: dictionary to translate cigar letters
    ref_consuming: stores reference consuming cigar operations"""

    with open(out_file, "a") as fileT:
        cigar_list = re.split(regexp, cigar)  # we split the string into a list
        # we parse the cigar list
        cigar_pos = 0  # to parse cigar_list
        cigar_length = 0  # to track the position in the genome
        total_SNV = 0  # counter to store the SNVs found in this aln
        total_INDEL = 0  # counter to store the SNVs found in this aln
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

                    #if cigar_list[cigar_pos + 1] == "X":
                        #total_SNV += 1
                    if cigar_list[cigar_pos + 1] == "I" or "D":
                        total_INDEL += 1

                if cigar_list[cigar_pos + 1] in ref_consuming:  # if cigar operation consumes reference
                    cigar_length += int(i)
                cigar_pos += 1
            else:
                cigar_pos += 1

        # We parse the md tag to get the SNVs
        md_regexp = r'(?<![A-Z^])[A-Z]'  # this pattern matches any letter not preceded either by ^ or [A-Z]
        md_SNV_list = re.findall(md_regexp, md_tag)
        total_SNV = len(md_SNV_list)

    return total_SNV, total_INDEL

def var_anno(out_file, cigar_tuple, md_list, regexp, ref_pos, aln_chr, symbol_add, dict_cigar, ref_consuming, window):
    """This function takes a cigar string and a md_tag and annotates all the changes that it contains in an output file,
    while it returns the amount of changes annotated. It returns the number of potential SNVs and INDELs within the
    window.
    regexp: splits the cigar string
    symbol_add: cigar operations to report
    dict_cigar: dictionary to translate cigar letters
    ref_consuming: stores reference consuming cigar operations"""

    with open(out_file, "a") as fileT:
        # we parse the cigar tuple
        cigar_pos = 0  # to parse cigar_list
        cigar_length = 0  # to track the position in the genome
        total_SNV = 0  # counter to store the SNVs found in this aln
        total_INDEL = 0  # counter to store the SNVs found in this aln
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

                    #if cigar_list[cigar_pos + 1] == "X":
                        #total_SNV += 1
                    if cigar_list[cigar_pos + 1] == "I" or "D":
                        total_INDEL += 1

                if cigar_list[cigar_pos + 1] in ref_consuming:  # if cigar operation consumes reference
                    cigar_length += int(i)
                cigar_pos += 1
            else:
                cigar_pos += 1

        # We parse the md tag to get the SNVs
        md_regexp = r'(?<![A-Z^])[A-Z]'  # this pattern matches any letter not preceded either by ^ or [A-Z]
        md_SNV_list = re.findall(md_regexp, md_tag)
        total_SNV = len(md_SNV_list)

    return total_SNV, total_INDEL


def dict_generator(fasta):
    """This function parses a fasta file (ideally the human reference genome) and returns a dictionary
        whose keys are the chr names and whose values are an integer representing their ordinal number"""

    chr_dict = {}  # create an empty dict
    ordinal = 1  # initialize the ordinal sequence

    with open(fasta, "r") as file:
        for line in file:
            if line.startswith(">"):
                key = line[1:].split()[0]  # we get the chr name
                chr_dict[key] = ordinal
                ordinal += 1
    return chr_dict


def chr_converter(chr_dict, chr_name):
    """This function translates the chr name from an alignment file to its integer value using the dictionary
    that has been created wit the dict_generator function."""

    chr_integer = chr_dict[chr_name]

    return chr_integer
