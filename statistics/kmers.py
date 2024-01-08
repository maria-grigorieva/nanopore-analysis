import pprint
import pandas as pd
from Bio import SeqIO
import matplotlib.pyplot as plt
import argparse

parser = argparse.ArgumentParser(description='Statistics of fastq files')
args_info = [
    ['-i', '--input', str, 'Path to the fastq input file', 'input_file'],
    ['-k', '--kmer_length', int, 'Length of kmers', 'kmer_length', 30],
    ['-m', '--mode', str, 'Mode', 'max'],
    ['-s', '--save', bool, 'Save to file', False]
]

for arg_info in args_info:
    parser.add_argument(arg_info[0], arg_info[1], type=arg_info[2], help=arg_info[3], default=arg_info[4])

args = parser.parse_args()

filename = args.input
kmer_length = args.kmer_length
save = args.save
mode = args.mode

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
        draw_histogram_sequence(occurences, reference)
    elif mode == 'top':
        # get top-10
        references = df['kmer'].iloc[0:10]
        occurences = [{'ref': r, 'occurences': get_all_occurrences(r, sequences)} for r in references]
        draw_histograms(occurences)

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
    return positions

def draw_histogram_sequence(data, reference, bins=50):
    plt.hist(data, bins=bins, edgecolor='black')
    plt.xlabel('Values')
    plt.ylabel('Frequency')
    plt.title(f'Distribution of start positions for {reference}')
    plt.show()
    if save:
        output_filename = f'{reference}_mer_freq.png'
        plt.savefig(output_filename)

def draw_histograms(arrays):
    fig, ax = plt.subplots()

    for item in arrays:
        val = item['ref']
        ax.hist(item['occurences'], alpha=0.5, bins='auto', density=True, label=f'{val}')
    ax.set_xlabel('Value')
    ax.set_ylabel('Frequency')
    ax.set_title('Histograms of Subarrays')
    ax.legend()

    plt.show()

if __name__ == "__main__":
    main()


