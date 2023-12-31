import pandas as pd
import logging
import matplotlib.pyplot as plt
import numpy as np
from Bio import SeqIO
from Bio.Seq import Seq
import glob
import re
import argparse
import os
from fuzzywuzzy import fuzz
import itertools


# Define the argument parser
parser = argparse.ArgumentParser(description='Search an aptamer among low-quality sequences with known length and primers')

# Define the arguments
args_info = [
    ['-al', '--alen', int, 'Length of an aptamer', 31],
    ['-i', '--input', str, 'Directory with input data', 'input_data'],
    ['-o', '--output', str, 'Directory with output data', 'output'],
    ['-pl', '--left_primer', str, 'Left Primer', None],
    ['-pr', '--right_primer', str, 'Right Primer', None],
    ['-r', '--ref', str, 'Initial reference sequence', 'auto'],
    ['-p', '--pos', int, 'Start position of the reference sequence', -1],
    ['-f', '--fuzzy', bool, 'Add fuzzy search', False],
    ['-s', '--save', bool, 'Save to excel (True/False)', False],
    ['-c', '--complement', bool, 'Add complementary primer', False],
    ['-w', '--weights', bool, 'Take into account nucleotide phred scores', False]
]

# Add arguments to the parser
for arg_info in args_info:
    parser.add_argument(arg_info[0], arg_info[1], type=arg_info[2], help=arg_info[3], default=arg_info[4])

# Parse the arguments
args = parser.parse_args()

# Access the argument values
APTAMER_LENGTH, PL, PR = args.alen, args.left_primer, args.right_primer
PRIMER_LENGTH = len(PL)
TOTAL_LENGTH = APTAMER_LENGTH + PRIMER_LENGTH * 2

DATA_DIR = args.input
OUTPUT_DIR = f'{DATA_DIR}/{args.output}'
PLOTS_DIR = f'{OUTPUT_DIR}/plots/references'
os.makedirs(OUTPUT_DIR, exist_ok=True)
os.makedirs(PLOTS_DIR, exist_ok=True)

REFERENCE = args.ref
REFERENCE_LENGTH = len(REFERENCE)
START_POS = args.pos
PRIMER_TYPE = 'left' if START_POS <= PRIMER_LENGTH - REFERENCE_LENGTH else 'right'
FUZZY = args.fuzzy
COMPLEMENT = args.complement
SAVE_FILES = args.save
PHRED_SCORES = args.weights

# Set all variations of primers
R_PL, R_PR = PL[::-1], PR[::-1]
C_PL, C_PR = str(Seq(PL).complement()), str(Seq(PR).complement())
RC_PL, RC_PR = C_PL[::-1], C_PR[::-1]


def main():

    seq_list = merge_sequences()

    INIT_REFERENCE = {'seq': args.ref,
                      'complement': str(Seq(PR).complement()[::-1][0:REFERENCE_LENGTH]),
                      'start_pos': args.pos}

    INIT_REFERENCE['sequences'],\
    INIT_REFERENCE['scores'] = select_ref_sequences(seq_list, INIT_REFERENCE)

    print_info(seq_list, INIT_REFERENCE)

    probabilities = []
    weights_list = []

    print(f'''Choosing the sequences that include the initial reference {INIT_REFERENCE['seq']} 
            at {INIT_REFERENCE['start_pos']} position...''')

    print(f'''{len(INIT_REFERENCE['sequences'])} have been selected with 
                {INIT_REFERENCE['seq']} at {INIT_REFERENCE['start_pos']}''')

    print(f'''Calculating probabilities of the occurence of letters at each position...''')
    freq = calculate_probabilities(INIT_REFERENCE['sequences'])
    probabilities.append(freq)

    weights = calculate_weights(INIT_REFERENCE['scores'])
    weights_list.append(weights)

    fname = f'''steps_{INIT_REFERENCE['seq']}'''

    if SAVE_FILES:
        steps_writer = pd.ExcelWriter(f'{OUTPUT_DIR}/{fname}.xlsx')
        write_steps_excel(freq, INIT_REFERENCE, steps_writer)
    else:
        steps_writer = None

    plot_probabilities(freq, INIT_REFERENCE)

    reference = INIT_REFERENCE

    print('Moving slicing window...')
    move_slicing_window(seq_list,
                        reference,
                        INIT_REFERENCE,
                        freq,
                        probabilities,
                        weights_list,
                        steps_writer)

    if SAVE_FILES:
        steps_writer.close()

    # Create an empty dictionary to store the merged data
    merged_data = {}

    print('Merging data from all letters probabilities cases...')
    for df in probabilities:
        for row in ['A', 'C', 'G', 'T']:
            merged_data[row] = df.loc[row] if row not in merged_data else \
                pd.concat([merged_data[row], df.loc[row]], axis=1)

    # Create an empty dict to store the merged weights data
    weights_df = pd.concat(weights_list)
    weights_df = scale_dataframe(weights_df)
    weights_df.reset_index(drop=True, inplace=True)

    fname = f'''output_{INIT_REFERENCE['seq']}'''
    if SAVE_FILES:
        output_writer = pd.ExcelWriter(f'{OUTPUT_DIR}/{fname}.xlsx')

    print(f'Calculating statistics...')
    final_composition_median, final_composition_low, final_composition_high = {}, {}, {}

    for row in ['A', 'C', 'G', 'T']:
        # Transpose the merged data for the row
        transposed_data = merged_data[row].transpose()
        transposed_data.reset_index(drop=True, inplace=True)

        if SAVE_FILES:
            transposed_data.to_excel(output_writer, sheet_name=f"{row}", index=True, header=True)

        stats = weighted_statistics(transposed_data, weights_df) if PHRED_SCORES else transposed_data.describe()

        if SAVE_FILES:
            stats.to_excel(output_writer, sheet_name=f"{row}-stats", index=True, header=True)

        final_composition_low[row] = stats.iloc[4].to_numpy()
        final_composition_median[row] = stats.iloc[5].to_numpy()
        final_composition_high[row] = stats.iloc[6].to_numpy()

    final_composition_low_df = pd.DataFrame(final_composition_low).T
    final_composition_median_df = pd.DataFrame(final_composition_median).T
    final_composition_high_df = pd.DataFrame(final_composition_high).T

    dfs = [final_composition_low_df, final_composition_median_df, final_composition_high_df]
    sheet_names = ["final-composition-low", "final-composition-median", "final-composition-high"]

    if SAVE_FILES:
        for df, sheet_name in zip(dfs, sheet_names):
            df.to_excel(output_writer, sheet_name=sheet_name, index=True, header=True)

        weights_df.to_excel(output_writer, sheet_name='Weights', index=True, header=True)

        output_writer.close()

    infer_sequence('0.25', final_composition_low_df, INIT_REFERENCE)
    infer_sequence('median', final_composition_median_df, INIT_REFERENCE)
    infer_sequence('0.75', final_composition_high_df, INIT_REFERENCE)

def weighted_statistics(df_A, df_B, threshold=0.4):
    statistics = {
        'count': [],
        'mean': [],
        'std': [],
        'min': [],
        '25%': [],
        '50%': [],
        '75%': [],
        'max': []
    }

    for column in df_A.columns:
        a_values = df_A.loc[df_B[column] >= threshold, column]
        statistics['count'].append(np.count_nonzero(~np.isnan(a_values)))
        statistics['mean'].append(np.mean(a_values))
        statistics['std'].append(np.std(a_values))
        statistics['min'].append(np.min(a_values))
        statistics['25%'].append(np.percentile(a_values, 25))
        statistics['50%'].append(np.percentile(a_values, 50))
        statistics['75%'].append(np.percentile(a_values, 75))
        statistics['max'].append(np.max(a_values))

    return pd.DataFrame(statistics).T


def print_info(seq_list, init_ref):
    # Print out the values of all the arguments
    print(f'Number of sequences in {DATA_DIR}: {len(seq_list)}')
    print(f"Length of an aptamer: {APTAMER_LENGTH}")
    print(f"Length of a primer: {PRIMER_LENGTH}")
    print(f"Left Primer: {PL}")
    print(f"Right Primer: {PR}")
    print(f"Directory with input data: {DATA_DIR}")
    print(f"Directory for the output: {OUTPUT_DIR}")
    print(f"Initial reference length: {REFERENCE_LENGTH}")
    print(f"Initial reference sequence: {init_ref['seq']}")
    print(f"Start position of the initial reference sequence: {init_ref['start_pos']}")

def infer_sequence(probability, df, init_ref):
    result = highest_probability_sequence(df.copy())
    plot_probabilities(df, {'seq': result, 'start_pos': 0}, probability)
    print_result(result, init_ref)

def print_result(result, ref):
    primer = result[:APTAMER_LENGTH] if PRIMER_TYPE == 'right' else result[:PRIMER_LENGTH]
    aptamer = result[APTAMER_LENGTH:APTAMER_LENGTH + PRIMER_LENGTH] \
        if PRIMER_TYPE == 'right' else result[PRIMER_LENGTH:PRIMER_LENGTH + APTAMER_LENGTH]

    print(f"Found candidate with reference {ref['seq']} at position {ref['start_pos']}:")
    print(f"{primer}---{aptamer}")

def merge_sequences():
    """
    Merge sequences from a specified directory directory
    """
    return [{'sequence': str(rec.seq),
            'score': rec.letter_annotations["phred_quality"]} for file_path in glob.glob(f'{DATA_DIR}/*.fastq') for rec in
                SeqIO.parse(file_path, "fastq")]

def extract_segment(sequence, score, pattern, matches, scores):
    interval = [{'start': m.start(0), 'end': m.end(0)} for m in re.finditer(pattern, sequence)]
    if len(interval) > 0:
        for i in interval:
            if not detect_glued_primers(sequence[i['start']:i['end']]):
                matches.append(sequence[i['start']:i['end']])
                scores.append(score[i['start']:i['end']])

def select_ref_sequences(seq_list, reference):
    """
    Choose the sequences that include the specified reference sequence
    at the given position from all available sequences.
    :param seq_list: list of all sequences
    :param reference: dict with info about reference, i.e.
    {'seq': 'GGCTTCTGG','start_pos': 31}
    :return: numpy array, i.e.
    [['G' 'C' 'G' ... 'T' 'G' 'C']
     ['G' 'C' 'A' ... 'G' 'C' 'A']
     ['C' 'C' 'G' ... 'T' 'G' 'C']
     ...
    """
    seq = reference['seq']
    left = reference['start_pos'] if PRIMER_TYPE == 'left' else reference['start_pos'] - PRIMER_LENGTH
    right = TOTAL_LENGTH - PRIMER_LENGTH - reference['start_pos'] - len(seq) \
            if PRIMER_TYPE == 'left' else TOTAL_LENGTH - reference['start_pos'] - len(seq)
    matches = []
    scores = []
    if not FUZZY:
        for item in seq_list:
            extract_segment(item['sequence'],
                            item['score'],
                            rf'.{{{left}}}{re.escape(seq)}.{{{right}}}',
                            matches,
                            scores)

            if COMPLEMENT:
                comp_seq = reference['complement']
                extract_segment(item['sequence'],
                                item['score'],
                                rf'.{{{left}}}{re.escape(comp_seq)}.{{{right}}}',
                                matches,
                                scores)

        matches = np.array([list(s) for s in matches]) if len(matches) > 0 else None
        scores = np.array([list(s) for s in scores]) if len(matches) > 0 else None
    else:
        fuzzy_matches = list(
            set([item['sequence'][i:i + REFERENCE_LENGTH] for item in seq_list for i in range(len(item['sequence']) - REFERENCE_LENGTH + 1)
                 if fuzz.ratio(item['sequence'][i:i + REFERENCE_LENGTH], seq) >= 90]))

        for item in seq_list:
            for fuz in fuzzy_matches:
                extract_segment(item['sequence'],
                                item['score'],
                                rf'.{{{left}}}{re.escape(fuz)}.{{{right}}}',
                                matches,
                                scores)
        matches = np.array([list(s) for s in matches]) if len(matches) > 0 else None
        scores = np.array([list(s) for s in scores]) if len(matches) > 0 else None
            # matches.extend(re.findall(rf'.{{{left}}}{re.escape(f)}.{{{right}}}', s) for f in filtered_sequences if
            #                len(filtered_sequences) > 0)
        # matches = list(itertools.chain.from_iterable(matches))

    return matches, scores

def detect_glued_primers(ref, length=6):
    return RC_PR[-length:] + PR[:length] in ref or RC_PL[:length] + PR[:length] in ref


def scale_dataframe(df):
    # Get the minimum and maximum values for each column
    min_vals = df.min()
    max_vals = df.max()

    # Calculate the range for each column
    ranges = max_vals - min_vals

    # Scale the values of the DataFrame from 0 to 1
    scaled_df = (df - min_vals) / ranges

    return scaled_df


def calculate_weights(scores):

    return pd.DataFrame(scores.mean(axis=0)).T


def calculate_probabilities(sequences):
    """
    Calculate probabilities of the appearence of each letter (A,C,G,T)
    at each position of a sequence
    :param matrix: the output of sequences()
    :return: pandas DataFrame with the following structure:
    Letter | 0    | 1    | 2   | ....| N
       A   | 0.12 | 0.98 | 1.0 | ... | ...
       C   | ...  | ...  | ... | ... | ...
       ...
    """
    return pd.DataFrame({letter: np.sum(sequences == letter, axis=0) / sequences.shape[0]
                         for letter in ['A', 'C', 'G', 'T']}).T


def update_reference(df, seq_list, reference, direction = -1):
    """
    Shift the current reference sequence to the left or right and generate a new matrix (numpy array)
    by extracting sequences from the array of all sequences.
    :param df: output of calculate_probabilities() - DataFrame with probabilities
    :param seq_list: list of all sequences
    :param reference: reference dictionary
    :param direction: directions = {'left': -1, 'right': 1}
    :return: new reference record:
       'seq' : reference sequence,
       'start_pos' : start position of the reference in aptamer,
       'sequences' : list of extracted sequences,
       'n_seq' : number of extracted sequences
    """
    new_start = reference['start_pos'] + direction
    index = {
        (-1, 'left'): new_start,
        (-1, 'right'): new_start - PRIMER_LENGTH,
        (1, 'left'): new_start + REFERENCE_LENGTH - 1,
        (1, 'right'): new_start + REFERENCE_LENGTH - PRIMER_LENGTH - 1
    }.get((direction, PRIMER_TYPE), 0)

    letters_probability = df[index].to_dict()

    #sorted_dict = dict(sorted(letters_probability.items(), key=lambda item: item[1], reverse=True))

    refs = []
    for k, v in letters_probability.items():
        try:
            new_ref = {'seq': k + reference['seq'][:-1] if direction == -1 else reference['seq'][1:] + k,
                       'start_pos': new_start}
            complement = str(Seq(PR).complement()[::-1][0:REFERENCE_LENGTH])
            new_ref['complement'] = k + complement[:-1] if direction == -1 else complement[1:] + k
            sequences, scores = select_ref_sequences(seq_list, new_ref)
            new_ref['n_seqs'], \
            new_ref['sequences'], \
            new_ref['scores'] = len(sequences), sequences, scores
            ref_name = new_ref['seq']
            n_seqs = new_ref['n_seqs']
            print(f'Check reference {ref_name}: number of sequences is {n_seqs}')
            refs.append(new_ref)
        except Exception as e:
            continue
    res = max(refs, key=lambda x: x['n_seqs'])
    print(
        f'''{len(res['sequences'])} have been selected with {res['seq']} at {res['start_pos']}''')
    return res


def write_steps_excel(freq, reference, writer=None):
    freq.to_excel(writer,
                  sheet_name=reference['seq'],
                  index=True,
                  header=True)
    pd.DataFrame(reference['sequences']).to_excel(writer,
                                                   sheet_name=reference['seq'],
                                                   startrow=freq.shape[0] + 2,
                                                   index=True,
                                                   header=True)


def move_slicing_window(seq_list, reference, init_ref, freq, probabilities, weights_list, writer):
    left_limit = PRIMER_LENGTH if PRIMER_TYPE == 'right' else 0
    right_limit = TOTAL_LENGTH - REFERENCE_LENGTH - 1 if PRIMER_TYPE == 'right' \
        else TOTAL_LENGTH - PRIMER_LENGTH - REFERENCE_LENGTH - 1

    def update_window(direction):
        nonlocal reference, freq
        reference = update_reference(freq,
                                     seq_list,
                                     reference,
                                     direction=direction)
        freq = calculate_probabilities(reference['sequences'])
        weights = calculate_weights(reference['scores'])

        if SAVE_FILES:
            write_steps_excel(freq, reference, writer)

        probabilities.append(freq)
        weights_list.append(weights)

        plot_probabilities(freq, reference)

    print('Moving left...')
    while reference['start_pos'] > left_limit:
        update_window(direction=-1)

    reference = init_ref
    freq = calculate_probabilities(reference['sequences'])

    print('The reference sequence has been reset to the initial value')

    print('Moving right...')
    while reference['start_pos'] < right_limit:
        update_window(direction=1)

    print('Moving slicing window has been finished!')


def highest_probability_sequence(df):
    """
    Choose sequences that have the greatest likelihood of letters appearing in all positions
    """
    return ''.join(df[column].idxmax() for column in df.columns)


def plot_probabilities(df, reference, title=None):
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
    figname = f'{PLOTS_DIR}/{ref_pos}:{ref_seq}.png' if title is None else f'{PLOTS_DIR}/{ref_seq}-{title}.png'
    plt.savefig(figname)
    plt.clf()
    plt.cla()
    plt.close()


if __name__ == "__main__":
    main()

