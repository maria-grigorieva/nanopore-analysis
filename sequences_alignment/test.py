from collections import Counter
from Bio import SeqIO
import glob
import re
import numpy as np
import pprint
import pandas as pd
from Bio import pairwise2
from difflib import SequenceMatcher

def merge_sequencies():
    """
    Merge sequencies from a specified directory directory
    :return:
    """
    return [str(rec.seq) for file_path in glob.glob('../data_samples/barcodes/FAP38830_pass_barcode12_bcc3428d_0.fastq') for rec in
            SeqIO.parse(file_path, "fastq")]

sequences = merge_sequencies()

def get_nmers(sequence, n):
    nmers = []
    for i in range(len(sequence) - n + 1):
        nmer = sequence[i:i+n]
        nmers.append(nmer)
    return nmers

def generate_versions(initial_sequence):
    versions = set()
    required_similarity = 0.9
    sequence_length = len(initial_sequence)

    for i in range(sequence_length):
        # Deletion
        version = initial_sequence[:i] + initial_sequence[i+1:]
        similarity = SequenceMatcher(None, initial_sequence, version).ratio()
        if similarity >= required_similarity:
            versions.add(version)

        for c in ['A', 'C', 'G', 'T']:
            # Modification
            version = initial_sequence[:i] + c + initial_sequence[i+1:]
            similarity = SequenceMatcher(None, initial_sequence, version).ratio()
            if similarity >= required_similarity:
                versions.add(version)

            # Insertion at position i
            version = initial_sequence[:i] + c + initial_sequence[i:]
            similarity = SequenceMatcher(None, initial_sequence, version).ratio()
            if similarity >= required_similarity:
                versions.add(version)

    return versions

initial_sequence = "GGCTTCTGGACTACCTATGC"
versions = generate_versions(initial_sequence)
print("Number of versions:", len(versions))
print(versions)

# # Subsequence and its known position
# test_sequence = "ATTATTACTTCGTTCAGTACATTGCTAAGGTTAATCCATTCCCTCCGATAGATGAAGCCAGCACCTCTCCTCTGACTGTAACCACGGCACATACAGTGGCATACACGTTTGACATGGGCATAGGTAGTCCAGAAGCCGGCTTCTGGACTACCTATGCATTTACCACGTTCGTGCAACCAATATACATCGTGGTTACAGTCAGAAATTGCTGGTTTCTCGATGGATGTACTGGCAA"
# # subsequence = "GGCTTCTGGACTAACTATGC" # edited
# subsequence = "GGCTTCTGGACTACCTATGC" # initial
#
# for a in pairwise2.align.globalms(subsequence, test_sequence, 2, -1, -2, -1):
#     print(pairwise2.format_alignment(*a))
#
# n = 20
# result = get_nmers(test_sequence, n)
#
# # for s in result:
# #     for a in pairwise2.align.globalms(subsequence, s, 1, -1, -0.5, -0.1):
# #         print(pairwise2.format_alignment(*a))
# # alignment = pairwise2.align.localxx(subsequence, test_sequence)
# # # print(alignment[0])
# # print(pairwise2.format_alignment(*alignment[0]))
# position = 31
#
# # Align sequences based on the known subsequence and position
# alignments = []
#
# # Define a custom scoring system for matches, mismatches, gaps, and opening/extension penalties
# match_score = 2
# mismatch_score = -1
# gap_score = -2
# open_penalty = -1
# # extend_penalty = -0.5
#
# alignments = []
# for sequence in sequences:
#     aligns =  pairwise2.align.globalms(subsequence,
#                                       sequence,
#                                       match_score,
#                                       mismatch_score,
#                                       gap_score,
#                                       open_penalty)
#     print(pairwise2.format_alignment(*aligns[0]))
#     alignments.append(aligns[0])
#     substrings = []
#     current_substring = ""
#
#     for char in aligns[0].seqA:
#         if char == "-":
#             current_substring += char
#         else:
#             if len(current_substring) <= 3:
#                 current_substring += char
#             elif current_substring:
#                 substrings.append(current_substring)
#                 current_substring = char
#             else:
#                 current_substring = char
#
#     # Add the last substring if it exists
#     if current_substring:
#         substrings.append(current_substring)
#
#     # Remove '-' characters between substrings with letters
#     final_substrings = []
#     temp_substring = ""
#     for substring in substrings:
#         if "-" in substring:
#             temp_substring += substring.replace("-", "")
#         else:
#             if temp_substring:
#                 final_substrings.append(temp_substring)
#                 temp_substring = ""
#             final_substrings.append(substring)
#
#     # Append the last temp_substring if it exists
#     if temp_substring:
#         final_substrings.append(temp_substring)
#
#     print(final_substrings)
#
# # # Print the alignments
# # for alignment in alignments:
# #     print(pairwise2.format_alignment(*alignment))
#
#
#
#
