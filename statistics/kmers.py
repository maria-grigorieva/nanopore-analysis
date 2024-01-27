import pprint
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
import matplotlib.pyplot as plt
import argparse
import plotly.graph_objects as go
from collections import Counter
from difflib import SequenceMatcher
import re
from fuzzywuzzy import fuzz


parser = argparse.ArgumentParser(description='Statistics of fastq files')
args_info = [
    ['-i', '--input', str, 'Path to the fastq input file', 'input_file'],
    ['-k', '--kmer_length', int, 'Length of kmers', 'kmer_length', 30],
    ['-m', '--mode', str, 'Mode', 'search'],
    ['-n', '--n_seq', int, 'Number of subsequences', 10],
    ['-s', '--save', bool, 'Save to file', False],
    ['-pl', '--left_primer', str, 'Left Primer', None],
    ['-pr', '--right_primer', str, 'Right Primer', None],
    ['-seq', '--sequence_to_check', str, 'Sequence for checking', None]
]

for arg_info in args_info:
    parser.add_argument(arg_info[0], arg_info[1], type=arg_info[2], help=arg_info[3], default=arg_info[4])

args = parser.parse_args()

filename = args.input
kmer_length = args.kmer_length
save = args.save
mode = args.mode
n_seq = args.n_seq
PL = args.left_primer
PR = args.right_primer
SEQ = args.sequence_to_check

R_PL = R_PR = C_PL = C_PR = RC_PL = RC_PR = None

# Set all variations of primers
if PL is not None and PR is not None:
    R_PL, R_PR = PL[::-1], PR[::-1]
    C_PL, C_PR = str(Seq(PL).complement()), str(Seq(PR).complement())
    RC_PL, RC_PR = C_PL[::-1], C_PR[::-1]

def main():

    sequences = [str(rec.seq) for rec in SeqIO.parse(filename, "fastq")]
    print(len(sequences))
    # filter out sequences with glued primers
    # sequences = [s for s in sequences if not glued_primers(s, 9)]
    # print(len(sequences))
    search_candidates(sequences, mode)

def search_candidates(sequences, mode='search'):

    if mode == 'search':

        df = count_kmers(sequences, kmer_length)
        save_to_file(df) if save else pprint.pprint(df)

    elif mode == 'check' and SEQ is not None:
        df = count_sequence_frequency(sequences, SEQ)

    references = df['kmer'].iloc[0:n_seq]
    occurences = [{'ref': r, 'occurences': get_all_occurrences(r, sequences)} for r in references]

    peaks, diffs = {}, []
    for i in occurences:
        ref = i['ref']
        extremums = []
        get_extremums_by_occurrences(i['occurences'], extremums)
        peaks[ref] = extremums
        differences = [extremums[i + 1] - extremums[i] for i in range(len(extremums) - 1)]
        diffs.extend(differences)

    draw_histograms(occurences, 'plotly')
    pprint.pprint(peaks)

    counter = Counter(diffs)
    max_occurrence = max(counter, key=counter.get)
    print(f'Most frequent distance between occurences is: {max_occurrence}')

    # Filtering and sorting
    filtered_items = [item for item in peaks.items() if abs(item[1][0] - item[1][1]) <= max_occurrence]
    sorted_items = sorted(filtered_items, key=lambda x: x[1][0])

    # Creating a new dictionary with the filtered and sorted items
    filtered_dict = {item[0]: item[1] for item in sorted_items}
    print("Filtered and sorted dictionary:")
    print(filtered_dict)

    sorted_items = sorted(filtered_dict.items(), key=lambda x: x[1][0])
    sequences = [i[0] for i in sorted_items]
    print(merge_strings(sequences))

def merge_strings(strings):
    merged_string = strings[0]  # Start with the first string in the array

    for string in strings[1:]:
        match = SequenceMatcher(None, string, merged_string).find_longest_match()
        common = string[match.a:match.a + match.size]
        if merged_string[-len(common):] == common:
            residual = string[len(common):]
            merged_string += residual

    return merged_string

def glued_primers(ref, length=6, threshold=90):
    substrings = [ref[i:i + length*2] for i in range(len(ref) - length*2 + 1)]
    for s in substrings:
        if (fuzz.ratio(s, RC_PR[-length:] + PR[:length]) >= threshold or \
            fuzz.ratio(s, RC_PL[-length:] + PR[:length]) >= threshold or \
            fuzz.ratio(s, RC_PR[-length:] + PL[:length]) >= threshold):
            return True
    return False

def save_to_file(df):
    output_filename = f'{kmer_length}_mer_freq.xlsx'
    df.to_excel(output_filename)

def find_repeatable_substrings(string):
    repeatable_substrings = set()
    # Find repeatable characters
    repeatable_characters = re.findall(r'((\w)\2{6,})', string)
    repeatable_substrings.update([substring[0] for substring in repeatable_characters])
    # Find repeatable bigrams
    repeatable_bigrams = re.findall(r'((\w{2})\2{4,})', string)
    repeatable_substrings.update([substring[0] for substring in repeatable_bigrams])
    # Find repeatable trigrams
    repeatable_trigrams = re.findall(r'((\w{3})\2{3,})', string)
    repeatable_substrings.update([substring[0] for substring in repeatable_trigrams])

    return repeatable_substrings


def count_kmers(sequences, kmer_length=30, step=10, threshold=5):

    kmer_frequencies = {}
    for sequence in sequences:
        for i in range(0, len(sequence) - kmer_length + 1, step):
            kmer = sequence[i:i + kmer_length]
            # remove all incorrect sequences (i.e. AAAAAAA, ACACACACAC, GGGGGGGGG, ...)
            if any([i in kmer for i in find_repeatable_substrings(kmer)]):
                continue
                # remove all subsequences similar to primers
                # ratio_pl = SequenceMatcher(None, kmer, PL).ratio()
                # ratio_pr = SequenceMatcher(None, kmer, PR).ratio()
                # if ratio_pl <= 0.3:
                #     if ratio_pr <= 0.3:
            else:
                if kmer not in kmer_frequencies:
                    kmer_frequencies[kmer] = 0
            kmer_frequencies[kmer] += 1
    df = pd.DataFrame(kmer_frequencies.items(), columns=['kmer', 'frequency'])
    df = df.sort_values('frequency', ascending=False)
    return df[df['frequency'] > threshold]

def count_sequence_frequency(sequences, seq):
    count = 0
    for item in sequences:
        count += item.count(seq)
    index = pd.Index(range(len(seq)))  # Replace `seq` with the appropriate scalar value
    df = pd.DataFrame({'kmer': seq, 'frequency': count}, index=index)
    return df


def get_all_occurrences(reference, all_sequences):
    positions = []
    for s in all_sequences:
        start = 0
        while True:
            index = s.find(reference, start)
            if index == -1:
                break
            positions.append(index)
            start = index + len(reference)
    return sorted(positions)

def count_values(numbers):
    value_counts = {}
    for num in numbers:
        if num in value_counts:
            value_counts[num] += 1
        else:
            value_counts[num] = 1
    sorted_counts = sorted(value_counts.items(), key=lambda x: x[1], reverse=True)
    return sorted_counts

def get_extremums_by_occurrences(arr, extremums, eps=10, chunks=2):
    counts = Counter(arr)
    max_occurrence = max(counts.values())
    peaks = [num for num, count in counts.items() if count == max_occurrence]
    extremums.append(peaks[0])
    arr = [num for num in arr if num >= peaks[0]+eps]
    chunks = chunks - 1
    if chunks > 0:
        get_extremums_by_occurrences(arr, extremums, eps, chunks)


# def draw_histogram_sequence(data, reference, bins=50):
#     plt.hist(data, bins=bins, edgecolor='black')
#     plt.xlabel('Position')
#     plt.ylabel('Frequency')
#     plt.title(f'Distribution of start positions for {reference}')
#     plt.show()
#     if save:
#         output_filename = f'{reference}_mer_freq.png'
#         plt.savefig(output_filename)
#

def draw_histogram_sequence(data, reference, bins=50, save=False):
    fig = go.Figure()
    fig.add_trace(go.Histogram(x=data, nbinsx=bins))

    fig.update_layout(
        xaxis_title='Position',
        yaxis_title='Frequency',
        title=f'Distribution of start positions for {reference}'
    )

    fig.show()
    # if save:
    #     output_filename = f'{reference}_mer_freq.png'
    #     pio.write_image(fig, output_filename)

def draw_histograms(arrays, type='matplotlib'):

    if type == 'matplotlib':
        fig, ax = plt.subplots()

        for item in arrays:
            val = item['ref']
            ax.hist(item['occurences'], alpha=0.5, bins='auto', density=True, label=f'{val}')
        ax.set_xlabel('Position')
        ax.set_ylabel('Frequency')
        ax.set_title('Histograms of Subarrays')
        ax.legend()

        plt.show()
    elif type == 'plotly':
        fig = go.Figure()

        for item in arrays:
            val = item['ref']
            fig.add_trace(go.Histogram(x=item['occurences'], name=val, histnorm='probability'))

        fig.update_layout(
            xaxis_title='Position',
            yaxis_title='Frequency',
            title='Histograms of Subarrays',
            showlegend=True
        )

        fig.show()


if __name__ == "__main__":
    main()


