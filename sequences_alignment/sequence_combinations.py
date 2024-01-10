import random
import itertools
import pprint

import random

def generate_deletions(sequence):
    modified_sequences = []
    for i in range(len(sequence)):
        for nucleotide in 'ACGT':
            modified_sequence = sequence[:i] + nucleotide + sequence[i+1:]
            modified_sequences.append(modified_sequence)
    return modified_sequences

def generate_insertions(sequence):
    modified_sequences = []
    for i in range(len(sequence) + 1):
        for nucleotide in 'ACGT':
            modified_sequence = sequence[:i] + nucleotide + sequence[i:]
            modified_sequence = modified_sequence[:len(sequence)]  # Cut to initial sequence length
            modified_sequences.append(modified_sequence)
    return modified_sequences

def generate_modifications(sequence):
    modified_sequences = []
    for i in range(len(sequence)):
        for nucleotide in 'ACGT':
            if sequence[i] != nucleotide:
                modified_sequence = sequence[:i] + nucleotide + sequence[i+1:]
                modified_sequences.append(modified_sequence)
    return modified_sequences

initial_sequence = "GGCTTCTGG"

deletion_sequences = generate_deletions(initial_sequence)
insertion_sequences = generate_insertions(initial_sequence)
modification_sequences = generate_modifications(initial_sequence)
# shift_sequences = generate_shifts(initial_sequence)

similar_sequences = [initial_sequence]
similar_sequences.extend(deletion_sequences)
similar_sequences.extend(insertion_sequences)
similar_sequences.extend(modification_sequences)
# similar_sequences.extend(shift_sequences)

similar_sequences = set(similar_sequences)

for i in similar_sequences:
    print(i)

print(len(similar_sequences))