import numpy as np
import pandas as pd
from extract_from_fastq import get_list_of_sequences
import matplotlib.pyplot as plt
import re
import os


ref = 'TAGGGAAACACGATAGAATCCGAACAGCACC'
ref_20 = 'TTACGTATTGCTAAGGTTAA'
ref_31 = 'TAGGGAAACACGATAGAATCCAAACACCCCC'

array_of_candidates= ['GCTAAGGTTAATAGGGAAACACGATAGAATC',
                        'TGCTAAGGTTAATAGGGAAACACGATAGAAT',
                        'TTGCTAAGGTTAATAGGGAAACACGATAGAA',
                        'ATTGCTAAGGTTAATAGGGAAACACGATAGA',
                        'TATTGCTAAGGTTAATAGGGAAACACGATAG',
                        'GTATTGCTAAGGTTAATAGGGAAACACGATA',
                        'AGGTTAATAGGGAAACACGATAGAATCCGAA',
                        'AAGGTTAATAGGGAAACACGATAGAATCCGA',
                        'TTAATAGGGAAACACGATAGAATCCGAACAG',
                        'TAAGGTTAATAGGGAAACACGATAGAATCCG',
                        'GTTAATAGGGAAACACGATAGAATCCGAACA',
                        'TAATAGGGAAACACGATAGAATCCGAACAGC',
                        'AATAGGGAAACACGATAGAATCCGAACAGCA',
                        'ATAGGGAAACACGATAGAATCCGAACAGCAC',
                        'TAGGGAAACACGATAGAATCCGAACAGCACC',
                        'AGGGAAACACGATAGAATCCGAACAGCACCT',
                        'GGTTAATAGGGAAACACGATAGAATCCGAAC',
                        'CTAAGGTTAATAGGGAAACACGATAGAATCC']

sequences = get_list_of_sequences('new_data/FAP38830_pass_barcode04_bcc3428d_0.fastq')


def extract_70_subsequences(sequences, ref):
    # aptamer is in the middle
    subsequences = []
    flank_length = 20
    aptamer_length = 31
    fragment_length = aptamer_length + flank_length*2

    for sequence in sequences:
        min_position = flank_length
        max_position = len(sequence) - (aptamer_length + flank_length + 1)  # Max possible position to ensure a 70-nucleotide subsequence
        start_position = sequence.find(ref)
        if start_position != -1:
            if start_position < max_position or start_position > min_position:
                subsequence = sequence[start_position - flank_length:start_position - flank_length + fragment_length]
                if len(subsequence) == fragment_length:
                    subsequences.append(subsequence)
    split_sequencies = [list(s) for s in subsequences]
    matrix = np.array(split_sequencies)
    return matrix


def extract_70_subsequences_primers(sequences, ref, offset=0):
    # start with primer to understand that it is the primer (sequence with the highest frequency)
    subsequences = []
    flank_length = 20
    aptamer_length = 31
    fragment_length = aptamer_length + flank_length*2

    for sequence in sequences:
        min_position = offset
        max_position = len(sequence) - (fragment_length + offset + 1)  # Max possible position to ensure a 70-nucleotide subsequence
        start_position = sequence.find(ref)
        if start_position != -1:
            if start_position < max_position or start_position >= min_position:
                subsequence = sequence[start_position - offset:start_position - offset + fragment_length]
                if len(subsequence) == fragment_length:
                    subsequences.append(subsequence)
    split_sequencies = [list(s) for s in subsequences]
    matrix = np.array(split_sequencies)
    return matrix


def extract_51_subsequences(sequences, ref, offset=0):
    # start with primer and add 31 letters to the right (for the left primer)
    # allows us to find the wildcard for our aptamer
    subsequences = []
    flank_length = 20
    aptamer_length = 31
    fragment_length = aptamer_length + flank_length

    for sequence in sequences:
        min_position = offset
        max_position = len(sequence) - (fragment_length + offset + 1)  # Max possible position to ensure a 70-nucleotide subsequence
        start_position = sequence.find(ref)
        if start_position != -1:
            if start_position < max_position or start_position >= min_position:
                subsequence = sequence[start_position - offset:start_position - offset + fragment_length]
                if len(subsequence) == fragment_length:
                    subsequences.append(subsequence)
    split_sequencies = [list(s) for s in subsequences]
    matrix = np.array(split_sequencies)
    return matrix


def extract_aptamer_substrings(strings, wildcard):
    substring_length = len(wildcard)

    extracted_substrings = []
    for string in strings:
        string_length = len(string)
        if substring_length > string_length:
            continue

        for i in range(string_length - substring_length + 1):
            substring = string[i: i + substring_length]
            match = True
            for j in range(substring_length):
                if wildcard[j] != '*' and wildcard[j] != substring[j]:
                    match = False
                    break
            if match:
                extracted_substrings.append(substring)

    split_sequencies = [list(s) for s in extracted_substrings]
    matrix = np.array(split_sequencies)
    print(matrix)
    return matrix


def calculate_probabilities(matrix):
    unique_letters = ['A','C','G','T']
    frequencies = {}
    for letter in unique_letters:
        column_counts = np.sum(matrix == letter, axis=0)
        frequencies[letter] = column_counts / matrix.shape[0]
    df = pd.DataFrame(frequencies).T
    print(df)
    return df


def highest_probability_sequence(df):
    # Infer the sequence of letters for each position
    sequence = ''
    for column in df.columns:
        max_letter = df[column].idxmax()
        sequence += max_letter

    print(f"The sequence of letters with the highest probability in each position: {sequence}")

    # Infer the sequence of letters with 100% probability, assuming it is continuous
    sequence = ''
    sequence_positions = []
    for column in df.columns:
        probabilities = df[column]
        if max(prob == 1.0 for prob in probabilities):
            max_letter = probabilities.idxmax()
            sequence += max_letter
            sequence_positions.append(1)
        else:
            sequence_positions.append(0)

    print(f"The sequence of letters with 100% continuous probability: {sequence}")
    if len(sequence) == 31 and sequence_positions.index(1) == 20:
        print("Aptamer has been found!!!")


def generate_aptamer_wildcard(df):
    # Search for continuous sequences starting from columns with index larger than 20
    # only for left primer
    columns = df.columns[df.columns.astype(int) >= 20]
    sequence = ''
    current_seq = ''

    for col in columns:
        probabilities = df[col]
        if max(prob > 0.5 for prob in probabilities):
            max_letter = probabilities.idxmax()
            current_seq += max_letter
        else:
            current_seq += '*'

    # Include the last sequence
    if current_seq:
        sequence += current_seq

    print(sequence)
    return sequence


def plot_aptamer(df, ref):
    # Set the colors for each letter
    colors = {'A': 'red', 'C': 'green', 'G': 'blue', 'T': 'orange'}

    # Plot the probabilities for each letter
    for letter in df.index:
        plt.plot(df.loc[letter], label=letter, color=colors[letter])

    # Increase the width of the plot to 20
    fig = plt.gcf()
    fig.set_size_inches(20, fig.get_figheight())

    # Set the legend, x-label, and y-label
    plt.legend()
    plt.xlabel('Position')
    plt.ylabel('Probability')

    # Add labels with the titles of the letters with the highest probability
    for column in df.columns:
        max_letter = df[column].idxmax()
        max_prob = df[column].max()
        plt.text(int(column), max_prob, max_letter, ha='center', va='bottom', fontsize=10)

    # Display the plot
    # plt.show()
    plt.savefig(f'new_data/output/plots/{ref}.png')
    plt.clf()
    plt.cla()
    plt.close()


# for c in array_of_candidates:
#     result = extract_70_subsequences(sequences, c)
#     freq = calculate_probabilities(result)
#     # plot_aptamer(freq, c)
#     highest_probability_sequence(freq)


# aptamers = extract_aptamer_substrings(sequences, 'TAGGGAAACAC*ATA*A**C**AA***C*C*')
# freq = calculate_probabilities(aptamers)
# plot_aptamer(freq, 'TAGGGAAACAC*ATA*A**C**AA***C*C*')
# highest_probability_sequence(freq)
# result = extract_51_subsequences(sequences, ref_20)
# print(result)
# freq = calculate_probabilities(result)
# print(freq)
# plot_aptamer(freq, ref_20)
# highest_probability_sequence(freq)
# generate_aptamer_wildcard(freq)
