import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from Bio import SeqIO
from collections import Counter
import glob
import re
import os

DATA_DIR = 'new_data'
OUTPUT_DIR = 'new_data/output'
INPUT_DATA = 'FAP38830_pass_barcode04_bcc3428d_0'
LETTERS = ['A', 'C', 'G', 'T']
file_pattern = '*.fastq'

APTAMER_LENGTH = 31
PRIMER_LENGTH = 20
TOTAL_LENGTH = APTAMER_LENGTH + PRIMER_LENGTH * 2
APTAMER_PRIMER = APTAMER_LENGTH + PRIMER_LENGTH
NUMBER_OF_STEPS = 20


right_primer = 'GCATAGGTAGTCCAGAAGCC'
left_primer = 'CTCCTCTGACTGTAACCACG'
ref = 'GGCTTCTGG'


def merge_sequencies():
    """
    Merge sequencies from a specified directory directory
    :return:
    """
    return [str(rec.seq) for file_path in glob.glob(f'{DATA_DIR}/{file_pattern}') for rec in
                SeqIO.parse(file_path, "fastq")]


def select_ref_sequencies(seq_list, reference):
    """
    Choose the sequences that include the specified reference sequence
    at the given position from all available sequences.
    :param seq_list: list of al sequencies
    :param reference: dict with info about reference, i.e.
    {'seq': 'GGCTTCTGG','start_pos': 31}
    :return: numpy array, i.e.
    [['G' 'C' 'G' ... 'T' 'G' 'C']
     ['G' 'C' 'A' ... 'G' 'C' 'A']
     ['C' 'C' 'G' ... 'T' 'G' 'C']
     ...
    """
    result = []
    seq = reference['seq']
    right = APTAMER_PRIMER - reference['start_pos'] - len(seq)
    left = reference['start_pos']
    pattern = rf'.{{{left}}}{re.escape(seq)}.{{{right}}}'

    for seq in seq_list:
        matches = re.findall(pattern, seq)
        for match in matches:
            result.append(match)

    return np.array([list(s) for s in result])

# def select_ref_sequencies(seq_list, reference):
#
#     result = []
#
#     reference_seq = reference['seq']
#     ref_start_pos = reference['start_pos']
#
#     for seq in seq_list:
#         index = seq.find(reference_seq)
#         if index >= ref_start_pos and index - ref_start_pos + TOTAL_LENGTH <= len(seq):
#             result.append(seq[index - ref_start_pos:index - ref_start_pos + TOTAL_LENGTH])
#     # print(f'Number of extracted sequencies with {reference_seq} at position {ref_start_pos} is: {len(result)}')
#
#     return np.array([list(s) for s in result])


# get k-mers frequences -- primers
def count_kmers(sequences, kmer_length=PRIMER_LENGTH, offset=1):
    kmer_frequencies = Counter(
        sequence[i:i + kmer_length] for sequence in sequences
        for i in range(0, len(sequence) - kmer_length + 1, offset))
    return pd.DataFrame(kmer_frequencies.items(),
                        columns=['kmer', 'frequency']).\
        sort_values('frequency', ascending=False)


# Save k-mer frequencies to Excel file
def save_kmers(df, input_data):
    kmer_len = len(df['kmer'].values[0])
    output_filename = f'{OUTPUT_DIR}/{input_data}_{kmer_len}_mer_freq.xlsx'
    df.to_excel(output_filename)


def calculate_probabilities(matrix):
    """
    Calculate probabilities of the appearence of each letter (A,C,G,T)
    at each position of a sequence
    :param matrix: the output of select_ref_sequencies()
    :return: pandas DataFrame with the following structure:
    Letter | 0    | 1    | 2   | ....| N
       A   | 0.12 | 0.98 | 1.0 | ... | ...
       C   | ...  | ...  | ... | ... | ...
       ...
    """
    return pd.DataFrame({letter: np.sum(matrix == letter, axis=0) / matrix.shape[0] for letter in LETTERS}).T


def update_reference(df, seq_list, reference, direction=-1):
    """
    Shift the current reference sequence to the left or right
    and generate a new matrix (numpy array) by extracting
    sequences from the array of all sequences.
    :param df: output of calculate_probabilities() - DataFrame with probabilities
    :param seq_list: list of all sequencies
    :param reference: reference dictionary
    :param direction: directions = {'left': -1, 'right': 1}
    :return: dictionary with the following structure:
       'seq' : reference sequence,
       'start_pos' : start position of the reference in aptamer,
       'sequencies' : list of extracted sequencies,
       'n_seq' : number of extracted sequencies
    """
    new_start = reference['start_pos'] + direction
    letters_probability = df[new_start].to_dict()
    sorted_dict = dict(sorted(letters_probability.items(), key=lambda item: item[1], reverse=True))
    refs = []
    for k, v in sorted_dict.items():
        new_ref = {'seq': k + reference['seq'][:-1] if direction == -1 else reference['seq'][1:] + k,
                   'start_pos': new_start}
        sequencies = select_ref_sequencies(seq_list, new_ref)
        new_ref['n_seqs'] = len(sequencies)
        new_ref['sequencies'] = sequencies
        refs.append(new_ref)
    return max(refs, key=lambda x: x['n_seqs'])


def highest_probability_sequence(df):
    return ''.join(df[column].idxmax() for column in df.columns)


def plot_probabilities(df, reference):
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
        max_letter, max_prob = df[column].idxmax(), df[column].max()
        plt.text(int(column), max_prob, max_letter, ha='center', va='bottom', fontsize=10)

    ref_seq = reference['seq']
    ref_pos = reference['start_pos']
    plt.savefig(f'new_data/output/plots/references/{ref_pos}:{ref_seq}.png')
    plt.clf()
    plt.cla()
    plt.close()


def generate_wildcard(df):
    columns = df.columns[df.columns.astype(int) >= 0]
    return ''.join([df[col].idxmax() if df[col].max() >= 0.5 else '*' for col in columns])


def wildcard_sequencies(sequencies, wildcard):
    length = len(wildcard)
    extracted_substrings = [seq[i: i + length] for seq in sequencies for i in
                            range(len(seq) - length + 1) if
                            not any(wildcard[j] != '*' and wildcard[j] != seq[i + j] for j in range(length))
                            if length <= len(seq)]
    return np.array([list(s) for s in extracted_substrings])



def aptamer_search():

    seq_list = merge_sequencies()
    # seq_list = [str(rec.seq) for rec in SeqIO.parse('new_data/FAP38830_pass_barcode11_bcc3428d_0.fastq', "fastq")]

    print(f'Number of sequences in {INPUT_DATA}: {len(seq_list)}')

    writer = pd.ExcelWriter('new_data/steps.xlsx')

    reference = {'seq': 'GGCTTCTGG',
                 'start_pos': 31
                 }

    candidates = []
    probabilities = []

    refs = select_ref_sequencies(seq_list, reference)

    print(refs)

    freq = calculate_probabilities(refs)

    freq.to_excel(writer, sheet_name=reference['seq'], index=True, header=True)
    pd.DataFrame(refs).to_excel(writer, sheet_name=reference['seq'], startrow=freq.shape[0] + 2,
                                                   index=True, header=True)

    plot_probabilities(freq, reference)
    candidates.append(highest_probability_sequence(freq))

    # move left
    while reference['start_pos'] > 0:
        reference = update_reference(freq, seq_list, reference, direction=-1)
        #refs = select_ref_sequencies(seq_list, reference)
        print(reference['seq'])
        print(reference['n_seqs'])
        print(reference['start_pos'])
        freq = calculate_probabilities(reference['sequencies'])
        #plot_probabilities(freq, reference)

        freq.to_excel(writer, sheet_name=reference['seq'], index=True, header=True)
        pd.DataFrame(reference['sequencies']).to_excel(writer, sheet_name=reference['seq'], startrow=freq.shape[0] + 2,
                                                       index=True, header=True)
        candidates.append(highest_probability_sequence(freq))
        probabilities.append(freq)

    # Reset the reference value to its initial state
    reference = {'seq': 'GGCTTCTGG',
                 'start_pos': 31
                 }

    # move right
    while reference['start_pos'] < APTAMER_PRIMER - len(reference['seq']):
        reference = update_reference(freq, seq_list, reference, direction=1)
        print(reference['seq'])
        print(reference['n_seqs'])
        print(reference['start_pos'])
        #refs = select_ref_sequencies(seq_list, reference)
        freq = calculate_probabilities(reference['sequencies'])

        freq.to_excel(writer, sheet_name=reference['seq'], index=True, header=True)
        pd.DataFrame(reference['sequencies']).to_excel(writer, sheet_name=reference['seq'], startrow=freq.shape[0] + 2,
                                                       index=True, header=True)

        #plot_probabilities(freq, reference)
        candidates.append(highest_probability_sequence(freq))
        probabilities.append(freq)

    print(candidates)
    writer.save()

    # Create an empty dictionary to store the merged data
    merged_data = {}

    # Merge the DataFrames
    for df in probabilities:
        # Iterate over each row
        for row in ['A', 'C', 'G', 'T']:
            # Check if the row exists in the merged data
            if row not in merged_data:
                merged_data[row] = df.loc[row]  # Create the row in the merged data
            else:
                merged_data[row] = pd.concat([merged_data[row], df.loc[row]],
                                             axis=1)  # Concatenate the data for the row
    writer = pd.ExcelWriter('output.xlsx')
    for row in ['A', 'C', 'G', 'T']:
        # Transpose the merged data for the row
        transposed_data = merged_data[row].transpose()
        # transposed_data = transposed_data.set_index('index', inplace=True)
        transposed_data.reset_index(drop=True, inplace=True)
        transposed_data.to_excel(writer, sheet_name=f"{row}", index=False)
        transposed_data.describe().to_excel(writer, sheet_name=f"{row}-stats", index=False)

    writer.save()
        #
        # transposed_data.to_csv(f'new_data/output/{row}.csv')


    # print(merged_data)

    for df in probabilities:
        # Iterate over each row
        for row in ['A', 'C', 'G', 'T']:
            # Check if the row exists in the merged data
            if row not in merged_data:
                merged_data[row] = df.loc[row]  # Create the row in the merged_data dictionary
            else:
                # Compare the probabilities element-wise and keep the maximum values
                merged_data[row] = pd.concat([merged_data[row], df.loc[row]], axis=1).max(axis=1)


    # Convert the merged data dictionary back to a DataFrame
    merged_df = pd.DataFrame(merged_data).T



    # Print the merged DataFrame
    print(merged_df)
    plot_probabilities(merged_df, {'seq':'X', 'start_pos':0})
    print(highest_probability_sequence(merged_df))


#
aptamer_search()
