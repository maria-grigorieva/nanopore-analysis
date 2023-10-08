import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from Bio import SeqIO
from Bio.Seq import Seq
import glob
import re
import argparse
import os
from collections import Counter

# ref = 'GGCTTCTGG'
# pos = 51

# Create an ArgumentParser object
parser = argparse.ArgumentParser(description='Search an aptamer among low-quality sequencies with known length and primers')

# Add arguments to the parser
parser.add_argument('-al', '--alen', type=int, help='Length of an aptamer', default=31)
parser.add_argument('-i', '--input', type=str, help='Directory with input data', default='input_data')
parser.add_argument('-o', '--output', type=str, help='Directory with output data', default='output')
parser.add_argument('-pl', '--left_primer', type=str, help='Left Primer')
parser.add_argument('-pr', '--right_primer', type=str, help='Right Primer')
parser.add_argument('-r', '--ref', type=str, help='Initial reference sequence', default='auto')
parser.add_argument('-p', '--pos', type=int, help='Start position of the reference sequence', default=-1)
parser.add_argument('-s', '--save', type=bool, help='Save to excel (True/False)', default=False)

# Parse the arguments
args = parser.parse_args()

# Access the argument values
APTAMER_LENGTH = args.alen
PL = args.left_primer
PR = args.right_primer
PRIMER_LENGTH = len(PL)
TOTAL_LENGTH = APTAMER_LENGTH + PRIMER_LENGTH * 2
APTAMER_PRIMER = APTAMER_LENGTH + PRIMER_LENGTH
DATA_DIR = args.input
OUTPUT_DIR = f'{DATA_DIR}/{args.output}'
PLOTS_DIR = f'{OUTPUT_DIR}/plots/references'
os.makedirs(OUTPUT_DIR, exist_ok=True)
os.makedirs(PLOTS_DIR, exist_ok=True)
REFERENCE_LENGTH = len(args.ref)
START_POS = args.pos
PRIMER_TYPE = 'left' if START_POS <= PRIMER_LENGTH - REFERENCE_LENGTH else 'right'
SAVE = args.save

# Set all variations of primers
R_PL = PL[::-1]
R_PR = PR[::-1]
C_PL = str(Seq(PL).complement())
C_PR = str(Seq(PR).complement())
RC_PL = C_PL[::-1]
RC_PR = C_PR[::-1]

def main():

    seq_list = merge_sequencies()

    if args.ref == 'auto' and args.pos == -1:
        ref, pos = auto_reference(seq_list, offset = 1)
        INIT_REFERENCE = {'seq': ref, 'start_pos': pos}
    else:
        INIT_REFERENCE = {'seq': args.ref, 'start_pos': args.pos}

    print(f'Number of sequences in {DATA_DIR}: {len(seq_list)}')

    # Print out the values of all the arguments
    print(f"Length of an aptamer: {APTAMER_LENGTH}")
    print(f"Length of a primer: {PRIMER_LENGTH}")
    print(f"Left Primer: {PL}")
    print(f"Right Primer: {PR}")
    print(f"Directory with input data: {DATA_DIR}")
    print(f"Directory for the output: {OUTPUT_DIR}")
    print(f"Initial reference length: {REFERENCE_LENGTH}")
    print(f"Initial reference sequence: {INIT_REFERENCE['seq']}")
    print(f"Start position of the initial reference sequence: {INIT_REFERENCE['start_pos']}")

    candidates, probabilities = [], []

    print(f'''Choosing the sequences that include the initial reference {INIT_REFERENCE['seq']} 
            at {INIT_REFERENCE['start_pos']} position...''')

    INIT_REFERENCE['sequencies'] = select_ref_sequencies(seq_list, INIT_REFERENCE)

    print(f'''{len(INIT_REFERENCE['sequencies'])} have been selected with 
                {INIT_REFERENCE['seq']} at {INIT_REFERENCE['start_pos']}''')

    print(f'''Calculating probabilities of the occurence of letters at each position...''')
    freq = calculate_probabilities(INIT_REFERENCE, True)

    if SAVE:
        fname = f'''steps_{INIT_REFERENCE['seq']}'''
        steps_writer = pd.ExcelWriter(f'{OUTPUT_DIR}/{fname}.xlsx')
        write_steps_excel(freq, INIT_REFERENCE, steps_writer)

    print('Save plots with frequencies')
    plot_probabilities(freq, INIT_REFERENCE)

    reference = INIT_REFERENCE

    candidates.append(highest_probability_sequence(freq))

    print('Moving slicing window...')
    move_slicing_window(seq_list,
                        reference,
                        INIT_REFERENCE,
                        freq,
                        candidates,
                        probabilities,
                        steps_writer if SAVE else None)

    if SAVE:
        steps_writer.close()

    # Create an empty dictionary to store the merged data
    merged_data = {}

    print('Merging data from all letters probabilities cases...')
    # Merge the probabilities of occurrences for each letter
    for df in probabilities:
        for row in ['A', 'C', 'G', 'T']:
            merged_data[row] = df.loc[row] if row not in merged_data else \
                pd.concat([merged_data[row], df.loc[row]], axis=1)

    if SAVE:
        fname = f'''output_{INIT_REFERENCE['seq']}'''
        output_writer = pd.ExcelWriter(f'{OUTPUT_DIR}/{fname}.xlsx')

    print(f'Calculating statistics...')
    final_composition = []
    for row in ['A', 'C', 'G', 'T']:
        # Transpose the merged data for the row
        transposed_data = merged_data[row].transpose()
        transposed_data.reset_index(drop=True, inplace=True)
        if SAVE:
            transposed_data.to_excel(output_writer,
                                     sheet_name=f"{row}",
                                     index=True,
                                     header=True)
            stats = transposed_data.describe()
            stats.to_excel(output_writer,
                            sheet_name=f"{row}-stats",
                            index=True,
                            header=True)
            final_composition.append(stats.iloc[5].values)
    final_composition_df = pd.DataFrame(final_composition, index=['A', 'C', 'G', 'T'])
    final_composition_df.to_excel(output_writer,
                            sheet_name="final-composition",
                            index=True,
                            header=True)
    if SAVE:
        output_writer.close()

    print('Infer the sequence having the highest probability:')
    for df in probabilities:
        # Iterate over each row
        for row in ['A', 'C', 'G', 'T']:
            # Check if the row exists in the merged data
            if row not in merged_data:
                merged_data[row] = df.loc[row]  # Create the row in the merged_data dictionary
            else:
                # Compare the probabilities element-wise and keep the maximum values
                merged_data[row] = pd.concat([merged_data[row], df.loc[row]], axis=1).quantile(q=0.75, axis=1)

    # Convert the merged data dictionary back to a DataFrame
    merged_df = pd.DataFrame(merged_data).T

    result = highest_probability_sequence(merged_df.copy())
    plot_probabilities(merged_df, {'seq':result, 'start_pos':0})
    if PRIMER_TYPE == 'right':
        print(f'''Found candidate with reference {INIT_REFERENCE['seq']} at position {INIT_REFERENCE['start_pos']}:
                {result[:APTAMER_LENGTH]}---
                {result[APTAMER_LENGTH:APTAMER_LENGTH+PRIMER_LENGTH]}''')
    elif PRIMER_TYPE == 'left':
        print(f'''Found candidate with reference {INIT_REFERENCE['seq']} at position {INIT_REFERENCE['start_pos']}:
                {result[:PRIMER_LENGTH]}---
                {result[PRIMER_LENGTH:PRIMER_LENGTH+APTAMER_LENGTH]}''')


def merge_sequencies():
    """
    Merge sequencies from a specified directory directory
    :return:
    """
    return [str(rec.seq) for file_path in glob.glob(f'{DATA_DIR}/*.fastq') for rec in
                SeqIO.parse(file_path, "fastq")]


def auto_reference(seq_list, offset = 1):
    kmer_frequencies = Counter(
        sequence[i:i + REFERENCE_LENGTH] for sequence in seq_list
        for i in range(0, len(sequence) - REFERENCE_LENGTH + 1, offset))
    df = pd.DataFrame(kmer_frequencies.items(),
                        columns=['kmer', 'frequency']). \
        sort_values('frequency', ascending=False)
    ref = df['kmer'].values[0]

    kmers = list(set(seq[i:i + TOTAL_LENGTH] for seq in seq_list
                     for i in range(0, len(seq) - TOTAL_LENGTH + 1, offset)))
    result = [{'seq': s, 'pos': s.find(ref)} for s in kmers if s.find(ref) != -1]
    pos_counts = Counter(item['pos'] for item in result)

    max_count = max(pos_counts.values())  # Get the maximum count
    max_count_pos = [pos for pos, count in pos_counts.items() if count == max_count]
    print(f'The automatically selected reference with length {REFERENCE_LENGTH} is {ref} at position {max_count_pos[0]}')
    return ref, max_count_pos[0]


def select_ref_sequencies(seq_list, reference):
    """
    Choose the sequences that include the specified reference sequence
    at the given position from all available sequences.
    :param seq_list: list of all sequencies
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

    left = reference['start_pos'] if PRIMER_TYPE == 'left' else reference['start_pos'] - PRIMER_LENGTH
    right = TOTAL_LENGTH - PRIMER_LENGTH - reference['start_pos'] - len(seq) \
            if PRIMER_TYPE == 'left' else TOTAL_LENGTH - reference['start_pos'] - len(seq)
    print(f'Positions: {left}-{seq}-{right}')
    # right = TOTAL_LENGTH - reference['start_pos'] - len(seq)
    # left = reference['start_pos']
    pattern = rf'.{{{left}}}{re.escape(seq)}.{{{right}}}'

    for seq in seq_list:
        matches = re.findall(pattern, seq)
        for match in matches:
            if not detect_glued_primers(match):
                result.append(match)

    return np.array([list(s) for s in result]) if len(result) > 0 else None


def detect_glued_primers(ref, length=6):
    return True if C_PR[-length:] + PR[:length] in ref else False

def calculate_probabilities(reference, weights = False):
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
    matrix = reference['sequencies']
    df = pd.DataFrame({letter: np.sum(matrix == letter, axis=0) / matrix.shape[0]
                         for letter in ['A', 'C', 'G', 'T']}).T
    if weights:
        start_index = reference['start_pos'] - PRIMER_LENGTH if PRIMER_TYPE == 'right' \
                        else reference['start_pos']
        end_index = reference['start_pos'] - PRIMER_LENGTH + REFERENCE_LENGTH if PRIMER_TYPE == 'right' \
                        else reference['start_pos'] + REFERENCE_LENGTH
        weights = add_weights(start_index, end_index)
        # Multiply numbers in each row with the weights
        df.iloc[:, 1:] = df.iloc[:, 1:].mul(weights, axis=1)

    return df


def add_weights(start_index, end_index):
    left_step, right_step = 1 / (start_index - 1), 1 / (PRIMER_LENGTH + APTAMER_LENGTH - end_index)
    right = [0] * (PRIMER_LENGTH + APTAMER_LENGTH - end_index)
    ref = [1] * (end_index - start_index)
    left = [0] * (start_index - 1)
    right = [right_step * (i + 1) for i in range(len(right))]
    left = [left_step * (i + 1) for i in range(len(left))]
    weights = [*left, *ref, *reversed(right)]
    return weights


def update_reference(df, seq_list, reference, direction = -1):
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
    # Choose the probabilities of letters at position 'new_start'
    # new_start - PRIMER_LENGTH + REFERENCE_LENGTH - 1
    index = {
        (-1, 'left'): new_start,
        (-1, 'right'): new_start - PRIMER_LENGTH,
        (1, 'left'): new_start + REFERENCE_LENGTH - 1,
        (1, 'right'): new_start + REFERENCE_LENGTH - PRIMER_LENGTH - 1
    }.get((direction, PRIMER_TYPE), 0)
    letters_probability = df[index].to_dict()

    # And sort them in descending order
    sorted_dict = dict(sorted(letters_probability.items(), key=lambda item: item[1], reverse=True))
    # Initialize an empty array for storing variations of references with different letters at position 'new_start'
    refs = []
    for k, v in sorted_dict.items():
        new_ref = {'seq': k + reference['seq'][:-1] if direction == -1 else reference['seq'][1:] + k,
                   'start_pos': new_start}
        try:
            sequencies = select_ref_sequencies(seq_list, new_ref)
            new_ref['n_seqs'], new_ref['sequencies'] = len(sequencies), sequencies
            refs.append(new_ref)
            new_ref = max(refs, key=lambda x: x['n_seqs'])
            print(
                f'''{len(new_ref['sequencies'])} have been selected with {new_ref['seq']} at {new_ref['start_pos']}''')
            return new_ref
        except Exception as e:
            continue



def write_steps_excel(freq, reference, writer=None):

    freq.to_excel(writer,
                  sheet_name=reference['seq'],
                  index=True,
                  header=True)
    pd.DataFrame(reference['sequencies']).to_excel(writer,
                                                   sheet_name=reference['seq'],
                                                   startrow=freq.shape[0] + 2,
                                                   index=True,
                                                   header=True)


def move_slicing_window(seq_list,
                        reference,
                        init_ref,
                        freq,
                        candidates,
                        probabilities,
                        writer):
    # move left
    print('Moving left...')
    left_limit = PRIMER_LENGTH if PRIMER_TYPE == 'right' else 0
    right_limit = TOTAL_LENGTH - REFERENCE_LENGTH - 1 if PRIMER_TYPE == 'right' \
        else TOTAL_LENGTH - PRIMER_LENGTH - REFERENCE_LENGTH - 1
    while reference['start_pos'] > left_limit:
        try:
            reference = update_reference(freq, seq_list, reference, direction=-1)
            freq = calculate_probabilities(reference)
            if SAVE:
                write_steps_excel(freq, reference, writer)
            candidates.append(highest_probability_sequence(freq))
            if reference['n_seqs'] >= 10:
                probabilities.append(freq)
            print('Save plots with frequencies')
            plot_probabilities(freq, reference)
        except Exception as e:
            continue

    # Reset the reference value to its initial state
    reference = init_ref
    freq = calculate_probabilities(reference)
    print('The reference sequence has been reset to the initial value')

    print('Moving right...')
    # move right
    while reference['start_pos'] < right_limit:
        try:
            reference = update_reference(freq, seq_list, reference, direction=1)
            freq = calculate_probabilities(reference)
            if SAVE:
                write_steps_excel(freq, reference, writer)
            candidates.append(highest_probability_sequence(freq))
            if reference['n_seqs'] >= 10:
                probabilities.append(freq)
            print('Save plots with frequencies')
            plot_probabilities(freq, reference)
        except Exception as e:
            continue


def highest_probability_sequence(df):
    """
    Choose sequences that have the greatest likelihood of letters appearing in all positions
    :param df:
    :return:
    """
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
    plt.savefig(f'{PLOTS_DIR}/{ref_pos}:{ref_seq}.png')
    plt.clf()
    plt.cla()
    plt.close()


if __name__ == "__main__":
    main()

