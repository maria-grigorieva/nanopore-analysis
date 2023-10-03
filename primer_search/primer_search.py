from collections import Counter
from Bio import SeqIO
import glob
import re
import numpy as np
import pprint
import pandas as pd
from Bio import pairwise2
import math
import HTSeq



def merge_sequencies():
    """
    Merge sequencies from a specified directory directory
    :return:
    """
    return [str(rec.seq) for file_path in glob.glob('../input_data/*.fastq') for rec in
                SeqIO.parse(file_path, "fastq")]


def find_primer_sequences():
    # Join all sequences into one string for easy searching
    seq_list = merge_sequencies()
    all_sequences = ''.join(seq_list)
    all_20mers = []

    # Find all 20-mers and count their frequencies
    for sequence in seq_list:
        all_20mers.extend(re.findall(r'(?=(.{20}))', sequence))

    frequency_counter = Counter(all_20mers)

    # Select 20-mers with frequency >= 70% of the total number of sequencies
    min_frequency = 0.6 * len(seq_list)
    selected_20mers = [seq for seq, freq in frequency_counter.items() if freq >= min_frequency]

    # Test the selected 20-mers for compliance with the conditions
    primers = []
    for seq in selected_20mers:
        seq_freq = 0
        distances = []
        for sequence in seq_list:
            # Find the occurrences of the 20-mer in each initial sequence
            occurrences = [m.start() for m in re.finditer(seq, sequence)]
            for i in range(len(occurrences)-1):
                dist = abs(occurrences[i] - occurrences[i+1])
                if dist <= 71:
                    seq_freq += 1
                    distances.append(dist)
        if seq_freq > 0:
            primers.append({'seq': seq,
                'freq': frequency_counter.get(seq),
                'distance_mean': np.mean(distances),
                'distance_min': np.min(distances),
                'distance_max': np.max(distances),
                'seq_freq': seq_freq}
               )

    return primers


# primers = find_primer_sequences()
# pprint.pprint(primers)
# df = pd.DataFrame(primers)
# print(df)

def convert_to_fasta(sequences, filenames):
    with open(filenames, 'w') as file:
        for i, sequence in enumerate(sequences):
            file.write(f'>Sequence {i+1}\n')
            file.write(f'{sequence}\n')


def search_neighbouring_sequences(sequences, selected_seq):
    neighbours = []
    for n_one in selected_seq:
        for n_two in selected_seq:
            freq = 0
            for seq in sequences:
                occurrences_one = [m.start() for m in re.finditer(n_one, seq)]
                occurrences_two = [m.start() for m in re.finditer(n_two, seq)]
                for one in occurrences_one:
                    for two in occurrences_two:
                        if abs(one - two) == 20:
                            freq += 1
            neighbours.append({'pair': (n_one, n_two),
                               'freq': freq}
                              )
    df = pd.DataFrame(neighbours)
    df = df.sort_values(by=['freq'], ascending=False)
    return df


def process_fastq_file(file_path):
    with HTSeq.FastqReader(file_path) as f:
        for read in f:
            print(read.seq)
            print(read.qualstr)
            print(read.qual)
            line = ''
            for i in range(len(read.seq)):
                if read.qual[i] >= 10:
                    line += str(read.seq)[i]
                else:
                    line += '*'
            print(line)
# Example usage:
fastq_file_path = '../input_data/FAP38830_pass_barcode04_bcc3428d_0.fastq'
process_fastq_file(fastq_file_path)

#
# seq_list = merge_sequencies()
# #
# # selected_seq = [seq for seq in seq_list if 'GATAGATGAAACCAGCACCT' in seq]
#
# all_sequences = ''.join(seq_list)
# all_20mers = []
#
# # Find all 20-mers and count their frequencies
# for sequence in seq_list:
#     all_20mers.extend(re.findall(r'(?=(.{20}))', sequence))
#
# frequency_counter = Counter(all_20mers)
#
# # Select 20-mers with frequency >= 70% of the total number of sequencies
# min_frequency = 0.6 * len(seq_list)
# selected_20mers = [seq for seq, freq in frequency_counter.items() if freq >= min_frequency]
#
# # convert_to_fasta(selected_seq, 'selected.fasta')

