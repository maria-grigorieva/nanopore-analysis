import pprint
import pandas as pd
from Bio import SeqIO
import matplotlib.pyplot as plt
import argparse
import plotly.graph_objects as go
from collections import Counter


parser = argparse.ArgumentParser(description='Statistics of fastq files')
args_info = [
    ['-i', '--input', str, 'Path to the fastq input file', 'input_file'],
    ['-k', '--kmer_length', int, 'Length of kmers', 'kmer_length', 30],
    ['-m', '--mode', str, 'Mode', 'max'],
    ['-n', '--n_seq', int, 'Number of subsequences', 10],
    ['-s', '--save', bool, 'Save to file', False]
]

for arg_info in args_info:
    parser.add_argument(arg_info[0], arg_info[1], type=arg_info[2], help=arg_info[3], default=arg_info[4])

args = parser.parse_args()

filename = args.input
kmer_length = args.kmer_length
save = args.save
mode = args.mode
n_seq = args.n_seq

def main():
    sequences = [str(rec.seq) for rec in SeqIO.parse(filename, "fastq")]
    df = count_kmers(sequences, kmer_length)
    if save:
        save_to_file(df)
    else:
        pprint.pprint(df)

    if mode == 'max':
        # get the most frequent sequence
        reference = df['kmer'].iloc[0]
        occurences = get_all_occurrences(reference, sequences)
        value_counts = count_values(occurences)
        pprint.pprint(value_counts)
        draw_histogram_sequence(occurences, reference)
    elif mode == 'top':
        # get topN
        references = df['kmer'].iloc[0:n_seq]
        occurences = [{'ref': r, 'occurences': get_all_occurrences(r, sequences)} for r in references]
        peaks = {}
        diffs = []
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
        filtered_items = [item for item in peaks.items() if abs(item[1][0] - item[1][1]) <= 50]
        sorted_items = sorted(filtered_items, key=lambda x: x[1][0])
        # Creating a new dictionary with the filtered and sorted items
        filtered_dict = {item[0]: item[1] for item in sorted_items}
        print("Filtered and sorted dictionary:")
        print(filtered_dict)

        sorted_items = sorted(filtered_dict.items(), key=lambda x: x[1][0])
        aligned_string = sorted_items[0][0]  # Set the initial aligned string as the first key
        for i in range(1, len(sorted_items)):
            key, values = sorted_items[i]
            prev_key, prev_values = sorted_items[i - 1]
            shift = values[0] - prev_values[0]
            aligned_string += key[-shift:]

        print("Aligned string:")
        print(aligned_string)


def save_to_file(df):
    output_filename = f'{kmer_length}_mer_freq.xlsx'
    df.to_excel(output_filename)

def count_kmers(sequences, kmer_length=70, step=1, threshold=5):

    kmer_frequencies = {}
    for sequence in sequences:
        for i in range(0, len(sequence) - kmer_length + 1, step):
            kmer = sequence[i:i + kmer_length]
            if kmer not in kmer_frequencies:
                kmer_frequencies[kmer] = 0
            kmer_frequencies[kmer] += 1
    df = pd.DataFrame(kmer_frequencies.items(), columns=['kmer', 'frequency'])
    df = df.sort_values('frequency', ascending=False)
    return df[df['frequency'] > threshold]


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


def draw_histogram_sequence(data, reference, bins=50):
    plt.hist(data, bins=bins, edgecolor='black')
    plt.xlabel('Position')
    plt.ylabel('Frequency')
    plt.title(f'Distribution of start positions for {reference}')
    plt.show()
    if save:
        output_filename = f'{reference}_mer_freq.png'
        plt.savefig(output_filename)

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


