from collections import Counter
from Bio import SeqIO
import glob
import re
import numpy as np
import pprint
import pandas as pd

def merge_sequencies():
    """
    Merge sequencies from a specified directory directory
    :return:
    """
    return [str(rec.seq) for file_path in glob.glob('../input_data/*.fastq') for rec in
                SeqIO.parse(file_path, "fastq")]

def find_consecutive_20mer_pairs(long_sequences, kmer_list):
    consecutive_pairs = []

    for sequence in long_sequences:
        for i in range(len(sequence) - 19):
            current_kmer = sequence[i:i+20]
            next_kmer = sequence[i+20:i+40]

            if current_kmer in kmer_list and next_kmer in kmer_list:
                consecutive_pairs.append((current_kmer, next_kmer))

    return consecutive_pairs



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
                'distance': np.median(distances),
                'seq_freq': seq_freq}
               )

    # Find consecutive 20-mers
    consecutive_20mers = find_consecutive_20mer_pairs(seq_list, selected_20mers)

    return primers, consecutive_20mers


primers, consecutive_20mers = find_primer_sequences()
pprint.pprint(primers)
df = pd.DataFrame(primers)
pprint.pprint(df)

pprint.pprint(consecutive_20mers)