"""This script contains all the functions created for the variant_detection script"""

import re
import pysam


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


def anno_cigar2(dict_indel_count, dict_snv_count, cigar, md_list, ref_pos, aln_chr, datasetIdx, ref_genome, aln):
    """This function takes a cigar string and annotates all the changes that it contains in an output file, while
    it returns the amount of changes annotated. It returns the number of potential SNVs and INDELs within the window.
    regexp: splits the cigar string
    cigar_indels: cigar operations to report
    dict_cigar: dictionary to translate cigar letters
    datasetIdx: 0 if tumoral, 1 if normal
    ref_consuming: stores reference consuming cigar operations"""

    #with open(out_file, "a") as fileT:
    regexp = r"(?<=[a-zA-Z=])(?=[0-9])|(?<=[0-9])(?=[a-zA-Z=])"  # regexp to split cigar string
    ref_consuming = ['M', 'D', 'N', '=', 'X']  # stores reference consuming cigar operations
    # read_consuming_only = ['S', 'H', 'M', 'I', 'N', '=', 'X'] # stores read consuming cigar operations
    read_consuming_only = ['S', 'H', 'I']  # stores read consuming cigar operations
    cigar_indels = ["I", "D"]  # cigar operations to report
    dict_cigar = {"X": "MM", "I": "INS", "D": "DEL"}
    dict_tuples = {0: "M", 1: "I", 2: "D", 3: "N", 4: "S", 5: "H", 6: "P", 7: "=", 8: "X", 9: "B"}
    # read_sequence = aln.get_forward_sequence()
    read_sequence = aln.query_sequence
    ref_sequence = aln.get_reference_sequence()
    read_consumed_bases = 0
    cigar_list = re.split(regexp, cigar)  # we split the string into a list
    # we parse the cigar list to get the INDELs
    cigar_list_idx = 0  # to parse cigar_list
    current_ref_cigar_len = 0  # to track the position in the genome
    # We parse the md tag to get the mismatches
    ref_mismatch_positions = []
    md_pos = 0  # to parse md_list
    md_length = 0  # to track the position in the genome
    for symbol in md_list:
        if symbol == "0":  # md separator
            pass  # ignore
        elif symbol[0] == "^":  # deletions
            md_length += len(symbol) - 1
        elif re.match(r'^\d', symbol):  # matches
            md_length += int(symbol)
        else:  # mismatches
            md_length += 1
            # pos_snv = ref_pos + md_length
            # substitution = symbol
            ref_mismatch_positions.append(md_length)
            """unique_string_snv = "{};{}".format(pos_snv, substitution)
            if unique_string_snv not in dict_snv_count:
                if datasetIdx == 0:
                    dict_snv_count[unique_string_snv] = TUMORAL_SINGLE_READ_VARIANT
                    # we report the variation
                    # fileT.write("{}\t".format(aln_chr))  # we write the chr
                    # fileT.write("{}\t".format(true_pos))  # we write the pos
                    # fileT.write("-\t")  # we write the ref
                    # fileT.write("-\t")  # we write the alt
                    # fileT.write("{}\t".format(dict_cigar[cigar_list[cigar_list_idx + 1]]))  # we write the type
                    # fileT.write("{}\t".format(window[3]))  # we write the VCF_TYPE
                    # fileT.write("{}\n".format(int(window[1]) + 1000))  # we write the VCF_POS
                if datasetIdx == 1:
                    dict_snv_count[unique_string_snv] = NORMAL_SINGLE_READ_VARIANT
            else:  # if the variation is already annotated in the dict
                var_code = dict_snv_count[unique_string_snv]
                print(f'Idx: {datasetIdx} SNV: {unique_string_snv} code: {var_code}')
                if datasetIdx == 0:
                    if var_code == NORMAL_SINGLE_READ_VARIANT or var_code == NORMAL_ONLY_VARIANT:
                        dict_snv_count[unique_string_snv] = TUMORAL_NORMAL_VARIANT
                    if var_code == TUMORAL_SINGLE_READ_VARIANT:
                        dict_snv_count[unique_string_snv] = TUMORAL_ONLY_VARIANT
                if datasetIdx == 1:
                    if var_code == TUMORAL_SINGLE_READ_VARIANT or var_code == TUMORAL_ONLY_VARIANT:
                        dict_snv_count[unique_string_snv] = TUMORAL_NORMAL_VARIANT
                    if var_code == NORMAL_SINGLE_READ_VARIANT:
                        dict_snv_count[unique_string_snv] = NORMAL_ONLY_VARIANT
                print(f'Idx: {datasetIdx} SNV: {unique_string_snv} code: {dict_snv_count[unique_string_snv]}')
                # dict_snv_count[unique_string_snv] += 1"""
        md_pos += 1
    mm_pos_idx = 0
    # mm_ref_pos = ref_mismatch_positions[mm_pos_idx]
    for cigar_list_idx, symbol in enumerate(cigar_list):
        if symbol.isdigit():
            cigar_op = cigar_list[cigar_list_idx + 1]
            if cigar_op in cigar_indels:
                true_pos = ref_pos + current_ref_cigar_len
                # we create the unique key for the dict_indel_count
                pos_init = true_pos
                pos_end = int(true_pos) + int(symbol)
                var_type = cigar_op
                length = symbol
                unique_string_indel = encode_var(pos_init, pos_end, var_type, length)
                # we store the variation in the dict
                if unique_string_indel not in dict_indel_count:
                    if datasetIdx == 0:
                        dict_indel_count[unique_string_indel] = TUMORAL_SINGLE_READ_VARIANT
                    if datasetIdx == 1:
                        dict_indel_count[unique_string_indel] = NORMAL_SINGLE_READ_VARIANT
                else:  # if the variation is already annotated in the dict
                    var_code = dict_indel_count[unique_string_indel]
                    if datasetIdx == 0:
                        if var_code == NORMAL_SINGLE_READ_VARIANT or var_code == NORMAL_ONLY_VARIANT:
                            dict_indel_count[unique_string_indel] = TUMORAL_NORMAL_VARIANT
                        if var_code == TUMORAL_SINGLE_READ_VARIANT:
                            dict_indel_count[unique_string_indel] = TUMORAL_ONLY_VARIANT
                    if datasetIdx == 1:
                        if var_code == TUMORAL_SINGLE_READ_VARIANT or var_code == TUMORAL_ONLY_VARIANT:
                            dict_indel_count[unique_string_indel] = TUMORAL_NORMAL_VARIANT
                        if var_code == NORMAL_SINGLE_READ_VARIANT:
                            dict_indel_count[unique_string_indel] = NORMAL_ONLY_VARIANT
            if cigar_op in ref_consuming:
                current_ref_cigar_len += int(symbol)
            if cigar_op == 'M' and len(ref_mismatch_positions) > 0: # works only if mismatches and matches are coded as M
                if len(ref_mismatch_positions) == 0 or mm_pos_idx >= len(ref_mismatch_positions):
                    continue
                mm_ref_pos = ref_mismatch_positions[mm_pos_idx]
                while mm_ref_pos < current_ref_cigar_len and mm_pos_idx < len(ref_mismatch_positions):
                    pos_in_read = mm_ref_pos + read_consumed_bases
                    # print(f'{aln.query_name} : {pos_in_read} : {len(read_sequence)}')
                    # if aln.query_name == '634d6d382503821a082c4ca8e34ffcf2':
                        # print(f'{mm_pos_idx} : {mm_ref_pos} : {cigar}')
                    substitution = read_sequence[pos_in_read - 1]
                    pos_snv = ref_pos + mm_ref_pos
                    # print(f'{pos_snv} : {substitution}')
                    unique_string_snv = "{};{}".format(pos_snv, substitution)
                    if unique_string_snv not in dict_snv_count:
                        if datasetIdx == 0:
                            dict_snv_count[unique_string_snv] = TUMORAL_SINGLE_READ_VARIANT
                        if datasetIdx == 1:
                            dict_snv_count[unique_string_snv] = NORMAL_SINGLE_READ_VARIANT
                    else:
                        var_code = dict_snv_count[unique_string_snv]
                        # print(f'Idx: {datasetIdx} SNV: {unique_string_snv} code: {var_code}')
                        if datasetIdx == 0:
                            if var_code == NORMAL_SINGLE_READ_VARIANT or var_code == NORMAL_ONLY_VARIANT:
                                dict_snv_count[unique_string_snv] = TUMORAL_NORMAL_VARIANT
                            if var_code == TUMORAL_SINGLE_READ_VARIANT:
                                dict_snv_count[unique_string_snv] = TUMORAL_ONLY_VARIANT
                        if datasetIdx == 1:
                            if var_code == TUMORAL_SINGLE_READ_VARIANT or var_code == TUMORAL_ONLY_VARIANT:
                                dict_snv_count[unique_string_snv] = TUMORAL_NORMAL_VARIANT
                            if var_code == NORMAL_SINGLE_READ_VARIANT:
                                dict_snv_count[unique_string_snv] = NORMAL_ONLY_VARIANT
                        # print(
                          #   f'Idx: {datasetIdx} SNV: {unique_string_snv} code: {dict_snv_count[unique_string_snv]}')
                    mm_pos_idx += 1
                    if mm_pos_idx < len(ref_mismatch_positions):
                        mm_ref_pos = ref_mismatch_positions[mm_pos_idx]
            if cigar_op in read_consuming_only:
                read_consumed_bases += int(symbol)
            if cigar_op == 'D':
                read_consumed_bases -= int(symbol)



def dict_generator(fasta_file):
    """This function parses a fasta file (ideally the human reference genome) and returns a dictionary
        whose keys are the chr names and whose values are an integer representing their ordinal number"""

    chr_dict = {}  # create an empty dict
    # ordinal = 1  # initialize the ordinal sequence
    ref_sequences = fasta_file.references
    chr_dict = {k: v for (k, v) in zip(ref_sequences, range(len(ref_sequences)))}
    #with open(fasta, "r") as file:
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
