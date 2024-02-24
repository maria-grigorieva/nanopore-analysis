import pandas as pd
import logging
import matplotlib.pyplot as plt
import numpy as np
from Bio import SeqIO
from Bio.Seq import Seq
import re
import argparse
import os
from fuzzywuzzy import fuzz
import tqdm
from PIL import Image, ImageDraw
from itertools import combinations
import uuid
from itertools import product
from collections import defaultdict
from Bio import pairwise2
from Bio.pairwise2 import format_alignment


# Define the argument parser
parser = argparse.ArgumentParser(description='Search an aptamer among low-quality sequences with known length and primers')

# Define the arguments
args_info = [
    ['-al', '--alen', int, 'Length of an aptamer', 31],
    ['-i', '--input', str, 'Path to the input fastq file', 'input_data'],
    ['-o', '--output', str, 'Directory with output data', '../results'],
    ['-pl', '--left_primer', str, 'Left Primer', None],
    ['-pr', '--right_primer', str, 'Right Primer', None],
    ['-r', '--ref', str, 'Initial reference sequence', 'auto'],
    ['-p', '--pos', int, 'Start position of the reference sequence', -1],
    ['-f', '--fuzzy', bool, 'Add fuzzy search', False],
    ['-s', '--save', bool, 'Save to excel (True/False)', False],
    ['-c', '--complement', bool, 'Add complementary primer', False],
    ['-w', '--weights', bool, 'Take into account nucleotide phred scores', False],
    ['-ph', '--cutoff', int, 'Phred score cut off', 15],
    ['-cl', '--initial_cleaning', bool, 'Initial cleaning of sequences from glued primers', True]
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

FILE_PATH = args.input
OUTPUT_DIR = f'../results/{args.output}'
PLOTS_DIR = f'{OUTPUT_DIR}/plots/references'
LOG_FILE = f'{OUTPUT_DIR}/app.log'
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
PHRED_CUTOFF = args.cutoff
CLEAN = args.initial_cleaning

# Set all variations of primers
R_PL, R_PR = PL[::-1], PR[::-1]
C_PL, C_PR = str(Seq(PL).complement()), str(Seq(PR).complement())
RC_PL, RC_PR = C_PL[::-1], C_PR[::-1]

PRIMERS = [{"primer_type": "PR", "primer_value": PR, "number_of_occurences": 0},
            {"primer_type": "PL", "primer_value": PL, "number_of_occurences": 0},
            {"primer_type": "R_PL", "primer_value": R_PL, "number_of_occurences": 0},
            {"primer_type": "R_PR", "primer_value": R_PR, "number_of_occurences": 0},
            {"primer_type": "C_PL", "primer_value": C_PL, "number_of_occurences": 0},
            {"primer_type": "C_PR", "primer_value": C_PR, "number_of_occurences": 0},
            {"primer_type": "RC_PR", "primer_value": RC_PR, "number_of_occurences": 0},
            {"primer_type": "RC_PL", "primer_value": RC_PL, "number_of_occurences": 0}]

GLUED_PRIMERS = []

logging.basicConfig(level=logging.INFO, filename=LOG_FILE, filemode="w",
                    format="%(asctime)s %(levelname)s %(message)s")

APTAMERS = []
CYCLES = 0

def main():

    get_glued_primers_combinations(length=7)

    print_info()

    seq_list = get_sequences(clean=CLEAN)
    logging.info(f'{len(seq_list)} sequences have been read!')
    # REFERENCE, REFERENCE_LENGTH, START_POS = search_best_primer(seq_list)

    # todo: remove whole sequences with glued primers?
    # seq_list = [s for s in seq_list if not glued_primers(s['sequence'], 9)]
    # print(len(seq_list))

    # result = ""
    #
    # while result is not None:
    #     INIT_REFERENCE = initialize(seq_list)
    #     result = searching(INIT_REFERENCE, seq_list, APTAMERS)
    #     if result is not None:
    #         APTAMERS.append(result)

    INIT_REFERENCE = initialize(seq_list)
    result = searching(INIT_REFERENCE, seq_list)
    APTAMERS.append(result)

def initialize(seq_list):
    ref = {'seq': args.ref,
          'complement': str(Seq(PR).complement()[::-1][0:REFERENCE_LENGTH]),
          'start_pos': args.pos}

    ref['sequences'], \
        ref['scores'] = select_ref_sequences(seq_list, ref)

    return ref

def merging_probabilities(probabilities):
    # Create an empty dictionary to store the merged data
    merged_data = {}
    logging.info('Merging data from all letters probabilities cases...')
    for df in probabilities:
        for row in ['A', 'C', 'G', 'T']:
            merged_data[row] = df.loc[row] if row not in merged_data else \
                pd.concat([merged_data[row], df.loc[row]], axis=1)
    return merged_data

def searching(ref, seq_list):
    probabilities = []
    weights_list = []

    logging.info(f'''Choosing the sequences that include the initial reference {ref['seq']} 
            at {ref['start_pos']} position...''')

    logging.info(f'''{len(ref['sequences'])} have been selected with 
                {ref['seq']} at {ref['start_pos']}''')

    logging.info(f'''Calculating probabilities of the occurence of letters at each position...''')
    freq = calculate_probabilities(ref['sequences'])
    probabilities.append(freq)

    weights = calculate_weights(ref['scores'])
    weights_list.append(weights)

    fname = f'''steps_{ref['seq']}'''

    if SAVE_FILES:
        steps_writer = pd.ExcelWriter(f'{OUTPUT_DIR}/{fname}.xlsx')
        write_steps_excel(freq, ref, steps_writer)
    else:
        steps_writer = None

    plot_probabilities(freq, ref)

    reference = ref

    logging.info('Moving slicing window...')
    try:
        move_slicing_window(seq_list,
                            reference,
                            ref,
                            freq,
                            probabilities,
                            weights_list,
                            steps_writer)

        if SAVE_FILES:
            steps_writer.close()

        # Create an empty dictionary to store the merged data
        merged_data = merging_probabilities(probabilities)
        # Create an empty dict to store the merged weights data
        weights_df = pd.concat(weights_list)
        #weights_df = scale_dataframe(weights_df)
        weights_df.reset_index(drop=True, inplace=True)

        logging.info(f'Calculating statistics...')
        result = calculate_statistics(merged_data, ref, weights_df)

        return result
    except:
        logging.warning('No more options')
        return None

def calculate_statistics(merged_data, ref, weights_df):

    fname = f'''output_{ref['seq']}'''
    if SAVE_FILES:
        output_writer = pd.ExcelWriter(f'{OUTPUT_DIR}/{fname}.xlsx')

    final_composition_median, final_composition_low, final_composition_high = {}, {}, {}

    for row in ['A', 'C', 'G', 'T']:
        # Transpose the merged data for the row
        transposed_data = merged_data[row].transpose()
        transposed_data.reset_index(drop=True, inplace=True)

        if SAVE_FILES:
            transposed_data.to_excel(output_writer, sheet_name=f"{row}", index=True, header=True)

        # stats = weighted_statistics(transposed_data, weights_df) if PHRED_SCORES else transposed_data.describe()
        stats = transposed_data.describe()

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

    result_25 = infer_sequence('0.25', final_composition_low_df, ref)
    result_median = infer_sequence('median', final_composition_median_df, ref)
    result_75 = infer_sequence('0.75', final_composition_high_df, ref)
    return result_median

def weighted_statistics(df_A, df_B, threshold=0.4):
    statistics = {
        'count': [], 'mean': [], 'std': [], 'min': [],
        '25%': [], '50%': [], '75%': [], 'max': []
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


def print_info():
    # Print out the values of all the arguments
    logging.info("CONFIGURATION PARAMETERS:")
    logging.info("===============================================================================")
    logging.info(f"Length of an aptamer: {APTAMER_LENGTH}")
    logging.info(f"Length of a primer: {PRIMER_LENGTH}")
    logging.info(f"Left Primer: {PL}")
    logging.info(f"Right Primer: {PR}")
    logging.info(f"Reversed Left Primer: {R_PL}")
    logging.info(f"Reversed Right Primer: {R_PR}")
    logging.info(f"Complementaty Left Primer: {C_PL}")
    logging.info(f"Complementary Right Primer: {C_PR}")
    logging.info(f"Complementaty Reversed Left Primer: {RC_PL}")
    logging.info(f"Complementary Reversed Right Primer: {RC_PR}")
    logging.info(f"Input file: {FILE_PATH}")
    logging.info(f"Directory for the output: {OUTPUT_DIR}")
    logging.info(f"Initial reference length: {REFERENCE_LENGTH}")
    logging.info(f"Initial reference sequence: {REFERENCE}")
    logging.info(f"Start position of the initial reference sequence: {START_POS}")
    logging.info(f"Type of a primer for searching: {PRIMER_TYPE}")
    logging.info(f"Fuzzy search: {FUZZY}")
    logging.info(f"Take into account nucleotide phred scores: {PHRED_SCORES}")
    logging.info(f"Phred scores cutoff: {PHRED_CUTOFF}")
    logging.info("===============================================================================")


def infer_sequence(probability, df, init_ref):
    result = highest_probability_sequence(df.copy())
    plot_probabilities(df, {'seq': result, 'start_pos': 0}, probability)
    print_result(result, init_ref, probability)
    return result

def print_result(result, ref, probability):
    primer = result[:APTAMER_LENGTH] if PRIMER_TYPE == 'right' else result[:PRIMER_LENGTH]
    aptamer = result[APTAMER_LENGTH:APTAMER_LENGTH + PRIMER_LENGTH] \
        if PRIMER_TYPE == 'right' else result[PRIMER_LENGTH:PRIMER_LENGTH + APTAMER_LENGTH]

    logging.info(f"Found candidate with reference {ref['seq']} at position {ref['start_pos']}:")
    logging.info(f'Probability: {probability}')
    logging.info(f"{primer}---{aptamer}")


def generate_fastq_file(data, output_filename):
    records = []
    for i, item in enumerate(data):
        sequence = item['sequence']
        score = item['score']
        record_id = str(uuid.uuid4())
        record = SeqIO.SeqRecord(Seq(sequence), id=record_id, name='', description='')
        record.letter_annotations["phred_quality"] = score
        records.append(record)

    SeqIO.write(records, output_filename, "fastq")
    return output_filename


def get_sequences(clean=True, mode='complete'):
    """
    Merge sequences from a specified directory
    """
    results = []
    records = list(SeqIO.parse(FILE_PATH, "fastq"))
    logging.info(f'Initial number of records in {FILE_PATH}: {len(records)}')
    logging.info('Reading file...')
    for rec in tqdm.tqdm(records, total=len(list(records)), leave=False, desc=FILE_PATH):
        sequence = str(rec.seq)
        score = rec.letter_annotations["phred_quality"]
        glued_pr = glued_primers(sequence)
        if glued_pr is not None:
            for r in split_sequence(sequence, score, glued_pr):
                results.append(r)
        else:
            results.append({'sequence': sequence,
                        'score': score})
    logging.info(f'{len(results)} sequences have been read!')
    create_image_with_colored_sequence(results[0:1000], f'{OUTPUT_DIR}/internal.png')
    generate_fastq_file(results, f'{OUTPUT_DIR}/internal.fastq')
    logging.info(f'Internal fastq file has been written.')

    logging.info("DETECTED GLUED PRIMERS:")
    logging.info("===============================================================================")
    for primer in GLUED_PRIMERS:
        type = primer['glued_primers_type']
        value = primer['glued_primers_value']
        n = primer['n_occurences']
        logging.info(f'Glued primers type: {type}, value = {value}, occurences = {n}')
    logging.info("===============================================================================")

    return results

def extract_segment(sequence, score, pattern, matches, scores, remove_incorrect=True, threshold=75):
    interval = [{'start': m.start(0), 'end': m.end(0)} for m in re.finditer(pattern, sequence)]
    if len(interval) > 0:
        for i in interval:
            s = sequence[i['start']:i['end']]

            if remove_incorrect:
                # remove all sequences with wrong initial reference at its start position
                # remove all incorrect sequences (i.e. AAAAAAA, ACACACACAC, GGGGGGGGG, ...)
                position = 0 if PRIMER_TYPE == 'left' else APTAMER_LENGTH
                primer = PR if PRIMER_TYPE == 'right' else PL
                #if calculate_similarity(s[position:position+PRIMER_LENGTH],primer) >= threshold:
                #if pairwise_similarity(s[position:position + PRIMER_LENGTH], primer) >= threshold:
                if fuzz.ratio(s[position:position+PRIMER_LENGTH],primer) >= threshold:
                    matches.append(s)
                    scores.append(score[i['start']:i['end']])
                #
                # incorrect_substrings = find_repeatable_substrings(s)
                # if len(incorrect_substrings) > 0:
                #     for s in incorrect_substrings:
                #         result_string = s
                #         for substring in incorrect_substrings:
                #             result_string = result_string.replace(substring, "")

            else:
                matches.append(s)
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
        # TODO: return back after debugging!!!!
        if PHRED_SCORES:
            matches = bad_phred_score_remover(matches, scores)
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
    return matches, scores

def bad_sequences_remover(matches, scores, threshold=15):
    for i in range(scores.shape[0] - 1, -1, -1):
        if np.mean(scores[i]) < threshold:
            scores = np.delete(scores, i, axis=0)
            # Remove the corresponding element in 'matches' using the same index
            matches = np.delete(matches, i, axis=0)
    return matches, scores

def bad_phred_score_remover(matches, scores, threshold=15):
    for i, subarray in enumerate(scores):
        for j, nucleotide in enumerate(subarray):
            if nucleotide <= threshold:
                matches[i][j] = None
    return matches

def split_sequence(sequence, score, indices):
    splitted = []
    start_index = 0
    for idx in indices:
        splitted.append({'sequence': sequence[start_index:idx],
                        'score': score[start_index:idx]})
        start_index = idx
    splitted.append({'sequence': sequence[start_index:],
                        'score': score[start_index:]})
    filtered_data = [item for item in splitted if len(item['sequence']) > TOTAL_LENGTH]
    #logging.info(f'Sequence {sequence} is split into {len(filtered_data)} subsequences as it contains glued primers {substr}:')
    return filtered_data


def get_glued_primers_combinations(length=round(PRIMER_LENGTH/2)):
    primers = [p["primer_value"] for p in PRIMERS]
    for i, j in combinations(primers, 2):
        glued_primers_type = []
        val1 = i[-length:] + j[:length]
        for p in PRIMERS:
            if i[-length:] in p['primer_value']:
                glued_primers_type.append(p['primer_type'])
            if j[:length] in p['primer_value']:
                glued_primers_type.append(p['primer_type'])
            if len(glued_primers_type) == 2:
                break
        GLUED_PRIMERS.append({'glued_primers_type': '--'.join(glued_primers_type),
                              'glued_primers_value': val1,
                              'n_occurences': 0})

        val2 = j[-length:] + i[:length]
        GLUED_PRIMERS.append({'glued_primers_type': '--'.join(glued_primers_type[::-1]),
                              'glued_primers_value': val2,
                              'n_occurences': 0})

def glued_primers(s, fuzzy=False, length=7, threshold=90):
    result = []
    for g in GLUED_PRIMERS:
        if fuzzy:
            substrings = [s[i:i + length*2] for i in range(len(s) - length*2 + 1)]
            tmp = [s.find(i) for i in substrings if fuzz.ratio(i,g['glued_primers_value']) >= threshold]
        else:
            tmp = [val for val in [i.start() for i in re.finditer(g['glued_primers_value'], s)]]
        if len(tmp) > 0:
            g['n_occurences'] += len(tmp)
        for val in tmp:
            result.append(round(val + len(g['glued_primers_value']) / 2))

    return sorted(result) if len(result)>0 else None


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

def get_combinations(length):
    letters = ['A', 'C', 'G', 'T']
    combinations = []

    for repeat in range(1, length+1):
        for combo in product(letters, repeat=repeat):
            if len(combo)==length:
                combinations.append(''.join(combo))
    return combinations

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

    # letters_probability = df[index].to_dict()
    refs = []
    step = abs(direction)
    # get letters combinations
    # comb = list(combinations(['A','C','G','T'], step))
    for letter in get_combinations(step):
    # for k, v in letters_probability.items():
    # for letter in ['A','C','G','T']:
        new_ref = {'seq': letter + reference['seq'][:-step] if direction == -step else reference['seq'][step:] + letter,
                   'start_pos': new_start}
        ref_name = new_ref['seq']
        complement = str(Seq(PR).complement()[::-step][0:REFERENCE_LENGTH])
        new_ref['complement'] = letter + complement[:-step] if direction == -step else complement[step:] + letter

        try:
            sequences, scores = select_ref_sequences(seq_list, new_ref)
            new_ref['n_seqs'], \
            new_ref['sequences'], \
            new_ref['scores'] = len(sequences), sequences, scores
            n_seqs = new_ref['n_seqs']

            #freq = calculate_probabilities(sequences)
            #inference = highest_probability_sequence(freq)
            new_ref['primer_score'] = evaluate_sequences_correctness(ref_name, sequences, scores)
            # new_ref['correctness'] = pairwise_similarity(sequences, scores)
            #new_ref['correctness'] = fuzz.ratio(inference[START_POS-PRIMER_LENGTH:], PR)
            correctness = new_ref['primer_score']
            new_ref['hits'] = n_seqs * correctness
            refs.append(new_ref)
        except Exception as e:
            continue
    res = max(refs, key=lambda x: x['hits'])

    logging.info(
        f'''{len(res['sequences'])} have been selected with {res['seq']} at {res['start_pos']}''')
    logging.info('*******************************************************************************')

    return res

def evaluate_sequences_correctness(ref_name, sequences, scores):
    seq = []
    sc = []
    aligned_primers = []
    for idx, row in enumerate(sequences):
        s = ''.join(row)
        position = 0 if PRIMER_TYPE == 'left' else APTAMER_LENGTH
        primer = PR if PRIMER_TYPE == 'right' else PL
        #seq.append(calculate_similarity(s[position:position+PRIMER_LENGTH], primer)/100)
        seq.append(pairwise_similarity(s[position:position+PRIMER_LENGTH], primer, aligned_primers)/100)
        # seq.append(round(fuzz.ratio(s[position:position+PRIMER_LENGTH], primer)/100,4))
        sc.append(1 - np.mean([pow(10, -i/10) for i in scores[idx][position:position+PRIMER_LENGTH]]))
        res = [np.mean([a,b]) for a, b in zip(seq, sc)]
    mean_sim_score = round(np.mean([i.score for i in aligned_primers]))
    n_seqs = len(sequences)
    primer_score = np.mean(res)
    hits = n_seqs * primer_score
    logging.info(f'{ref_name}: Number of sequences is {n_seqs}, correctness = {primer_score}, hits = {hits}')
    for i in aligned_primers:
        if i.score == mean_sim_score:
            for s in format_alignment(*i).split('\n'):
                logging.info(s)
            break
    return primer_score

def calculate_similarity(string1, string2):
    count_same = sum(c1 == c2 for c1, c2 in zip(string1, string2))
    percentage_similarity = (count_same / len(string1))*100
    return percentage_similarity

def pairwise_similarity(string1, string2, aligned_primers):
    alignments = pairwise2.align.globalxx(string1, string2)
    mean_score = np.mean([i.score for i in alignments])
    # save alignments for logging
    for a in alignments:
        if a.score == mean_score:
            aligned_primers.append(a)
            break
    return mean_score * 100 / len(string1)

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


def move_slicing_window(seq_list, reference, init_ref, freq, probabilities, weights_list, writer, step=1):
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

    logging.info("===============================================================================")
    logging.info('MOVING LEFT...')
    logging.info("===============================================================================")
    pbar = tqdm.tqdm(total=(reference['start_pos'] - left_limit))
    while reference['start_pos'] > left_limit + step:
        update_window(direction=-step)
        pbar.update(1)  # Increment the progress bar by 1
    pbar.close()  # Close the progress bar once the loop is finished

    reference = init_ref
    freq = calculate_probabilities(reference['sequences'])
    logging.info("===============================================================================")
    logging.info('The reference sequence has been reset to the initial value')
    logging.info('\n')
    logging.info("===============================================================================")
    logging.info('MOVING RIGHT...')
    logging.info("===============================================================================")
    pbar = tqdm.tqdm(total=(left_limit-reference['start_pos']))
    while reference['start_pos'] < right_limit - step:
        update_window(direction=step)
        pbar.update(1)  # Increment the progress bar by 1
    pbar.close()  # Close the progress bar once the loop is finished
    logging.info("===============================================================================")

    logging.info('Moving slicing window has been finished!')


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

def create_image_with_colored_sequence(records, output_file, limit=200):
    width = 20  # Width of each character box
    height = 20  # Height of each character box
    padding = 5  # Padding between character boxes
    max_score = 40  # Maximum score

    total_width = limit * (width + padding)
    total_height = height * len(records)

    image = Image.new("RGB", (total_width, total_height), "white")
    draw = ImageDraw.Draw(image)

    y = 0

    for record in records:
        sequence = record['sequence']
        scores = record['score']

        x = 0

        for i, score in enumerate(scores):
            character = sequence[i]

            # Calculate the color based on the score
            normalized_score = score / max_score  # Normalize the score between 0 and 1
            red = int(255 * (1 - normalized_score))
            green = int(255 * normalized_score)
            color = (red, green, 0)

            # Draw the character box with the corresponding color
            draw.rectangle([x, y, x + width, y + height], fill=color)
            draw.text((x+8, y), character, fill="black")

            x += width + padding

        y += height

    image.save(output_file)

def find_repeatable_substrings(string):
    repeatable_substrings = set()
    # Find repeatable characters
    repeatable_characters = re.findall(r'((\w)\2{5,})', string)
    repeatable_substrings.update([substring[0] for substring in repeatable_characters])
    # Find repeatable bigrams
    repeatable_bigrams = re.findall(r'((\w{2})\2{4,})', string)
    repeatable_substrings.update([substring[0] for substring in repeatable_bigrams])
    # Find repeatable trigrams
    repeatable_trigrams = re.findall(r'((\w{3})\2{3,})', string)
    repeatable_substrings.update([substring[0] for substring in repeatable_trigrams])
    return list(repeatable_substrings)


# def search_best_primer(sequences):
#     ref = ''
#     ref_len = 0
#     start_pos = 0
#     for item in PRIMERS:
#         total_occurrences = 0
#         for s in sequences:
#             occurrences = s['sequence'].count(item['primer_value'])
#             total_occurrences += occurrences
#         item['number_of_occurences'] = total_occurrences
#
#     best_primer = max(PRIMERS, key=lambda x: x['number_of_occurences'])
#     if best_primer['primer_type'] in ['PL','R_PL','C_PL','RC_PL']:
#         ref_len = round(PRIMER_LENGTH * 0.45)
#         start_pos = PRIMER_LENGTH - ref_len
#         ref = best_primer['primer_value'][-ref_len:]
#     elif best_primer['primer_type'] in ['PR', 'R_PR', 'C_PR', 'RC_PR']:
#         ref_len = round(PRIMER_LENGTH * 0.45)
#         start_pos = PRIMER_LENGTH + APTAMER_LENGTH
#         ref = best_primer['primer_value'][:ref_len]
#     return ref, ref_len, start_pos



if __name__ == "__main__":
    main()

