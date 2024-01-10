import pandas as pd
from difflib import SequenceMatcher
from itertools import combinations
import numpy as np
from Bio import SeqIO
import matplotlib.pyplot as plt
from Bio.Seq import Seq
from collections import Counter
import pprint

def count_kmers(sequences, kmer_length=70, step=1):

    kmer_frequencies = {}
    for sequence in sequences:
        for i in range(0, len(sequence) - kmer_length + 1, step):
            kmer = sequence[i:i + kmer_length]
            if kmer not in kmer_frequencies:
                kmer_frequencies[kmer] = 0
            kmer_frequencies[kmer] += 1
    df = pd.DataFrame(kmer_frequencies.items(), columns=['kmer', 'frequency'])
    df = df.sort_values('frequency', ascending=False)
    return df

def save_kmers(df, kmer_length=70):
    output_filename = f'{kmer_length}_mer_freq.xlsx'
    df.to_excel(output_filename)

def longest_common_sequence(mode='df', df=None, fpath=None, method='sequentive', length=20, threshold=5):
    if mode == 'file':
        df = pd.read_excel(fpath)
        df = df.sort_values(by='frequency', ascending=False)

    sequences = df[df['frequency']>threshold]['kmer'].tolist()
    aptamers = []
    if method == 'sequentive':
        start = sequences[0]
        for i in sequences[1:]:
            match = SequenceMatcher(None, start, i).find_longest_match()
            start = start[match.a:match.a + match.size]
            if match.size >= length:
                aptamers.append(start)
    elif method == 'pairwise':
        for pair in combinations(sequences, 2):
            string1, string2 = pair
            match = SequenceMatcher(None, string1, string2).find_longest_match(0, len(string1), 0, len(string2))
            if match.size >= length:
                aptamers.append(string1[match.a:match.a + match.size])

    result = list(set(aptamers))

    return result

def draw_histogram(data, left, right):
    # Plotting the histogram
    plt.hist(data, bins=50, edgecolor='black')
    # Adding labels and title
    plt.xlabel('Values')
    plt.ylabel('Frequency')
    plt.title(f'Histogram for {left} and {right}')
    # Displaying the histogram
    plt.show()


def draw_histogram_sequence(data, reference):
    # Plotting the histogram
    plt.hist(data, bins=50, edgecolor='black')
    # Adding labels and title
    plt.xlabel('Values')
    plt.ylabel('Frequency')
    plt.title(f'Histogram for {reference}')
    # Displaying the histogram
    plt.show()

def primers_distances(sequence, left, right):
    dist = []
    left_start = 0
    right_start = 0
    while True:
        left_index = sequence.find(left, left_start)
        if left_index == -1:
            break
        right_index = sequence.find(right, right_start)
        if right_index == -1:
            break
        if right_index > left_index:
            dist.append(right_index - left_index)
        left_start = left_index + len(left)
        right_start = right_index + len(right)
    return dist

def sequence_distances(sequence, reference):
    dist = []
    start = 0
    while True:
        index = sequence.find(reference, start)
        if index == -1:
            break
        dist.append(index)
        start = index + len(reference)
    return dist

def find_all_indices(string, substring):
    indices = []
    start = 0
    while True:
        index = string.find(substring, start)
        if index == -1:
            break
        indices.append(index)
        start = index + len(substring)
    return indices

def calculate_diffs(a, b):
    min_len = min(len(a), len(b))  # Find the minimum length between a and b
    diffs = []
    for i in range(min_len):
        diff = b[i] - a[i]
        diffs.append(diff)
    return diffs

def check_primers(all_sequences, mode='df', df=None, fpath=None, primers_length=20, threshold=5):
    if mode == 'file':
        df = pd.read_excel(fpath)
        df = df.sort_values(by='frequency', ascending=False)

    sequences = df[df['frequency']>threshold]['kmer'].tolist()

    for i in sequences:
        left_primer = i[0:primers_length][-10:]
        right_primer = i[-primers_length:][0:10]
        if left_primer == right_primer:
            break
        print(f'left: {left_primer}, right: {right_primer}')
        distances = []
        for s in all_sequences:
            # get distances between left and right primers in all sequences
            # left_occurences = find_all_indices(s, left_primer)
            # right_occurences = find_all_indices(s, right_primer)
            diffs = primers_distances(s, left_primer, right_primer)
            if len(diffs) > 0:
                #draw_histogram(diffs, left_primer, right_primer)
                distances.extend(diffs)
            #calculate_diffs(left_occurences, right_occurences, distances)
        print(distances)
        draw_histogram(distances, left_primer, right_primer)
        print(np.mean(distances))
        print(len(distances))


def count_values(numbers):
    value_counts = {}
    for num in numbers:
        if num in value_counts:
            value_counts[num] += 1
        else:
            value_counts[num] = 1
    sorted_counts = sorted(value_counts.items(), key=lambda x: x[1], reverse=True)
    return sorted_counts


def draw_histograms(arrays):
    colors = ['blue', 'green', 'red', 'purple', 'orange', 'yellow']
    fig, ax = plt.subplots()

    for i, array in enumerate(arrays):
        ax.hist(array, alpha=0.7, bins='auto', density=True, label=f'Subarray {i + 1}', color=colors[i % len(colors)])

    ax.set_xlabel('Value')
    ax.set_ylabel('Frequency')
    ax.set_title('Histograms of Subarrays')
    ax.legend()

    plt.show()


def get_distances(reference, all_sequences):
    distances = []
    positions = []
    for s in all_sequences:
        print(s)
        pos = sequence_distances(s, reference)
        print(pos)
        if len(pos) > 0:
            positions.extend(pos)
        if len(pos) > 1:
            d = [(pos[i+1] - pos[i]) for i in range(len(pos)-1)]
            distances.extend(d)
    print(positions)
    print(distances)
    # draw_histogram_sequence(distances, reference)
    draw_histogram_sequence(positions, reference)
    print(np.mean(positions))
    print(len(positions))
    print(count_values(positions))
    print(np.mean(distances))
    print(len(distances))


sequences = [str(rec.seq) for rec in SeqIO.parse('../input_data/merged.fastq', "fastq")]
# df = count_kmers(sequences, kmer_length=71)
# save_kmers(df, 71)
result = longest_common_sequence(mode='file', fpath='31_mer_freq.xlsx', method='sequentive', threshold=50)
pprint.pprint(result)
get_distances(result[0], sequences)
# check_primers(sequences, mode='file', fpath='71_mer_freq.xlsx', threshold=50)


