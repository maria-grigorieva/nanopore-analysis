import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from Bio import SeqIO
from collections import Counter
import glob

DATA_DIR = 'new_data'
OUTPUT_DIR = 'new_data/output'
FILETYPE = 'fastq'
INPUT_DATA = 'FAP38830_pass_barcode04_bcc3428d_0'
LETTERS = ['A', 'C', 'G', 'T']
file_pattern = '*.fastq'

APTAMER_LENGTH = 31
PRIMER_LENGTH = 20
    # 'FAP38830_pass_barcode04_bcc3428d_0'
    # 'FAP38830_pass_barcode05_bcc3428d_0'
    # 'FAP38830_pass_barcode11_bcc3428d_0'

def merge_sequencies():
    return [str(rec.seq) for file_path in glob.glob(f'{DATA_DIR}/{file_pattern}') for rec in
                SeqIO.parse(file_path, "fastq")]


# get k-mers frequences -- primers
def count_kmers(sequences, kmer_length=PRIMER_LENGTH, offset=1):
    kmer_frequencies = Counter(
        sequence[i:i + kmer_length] for sequence in sequences for i in range(0, len(sequence) - kmer_length + 1, offset))
    return pd.DataFrame(kmer_frequencies.items(), columns=['kmer', 'frequency']).sort_values('frequency', ascending=False)


# Save k-mer frequencies to Excel file
def save_kmers(df, input_data):
    kmer_len = len(df['kmer'].values[0])
    output_filename = f'{OUTPUT_DIR}/{input_data}_{kmer_len}_mer_freq.xlsx'
    df.to_excel(output_filename)


def calculate_probabilities(matrix):
    return pd.DataFrame({letter: np.sum(matrix == letter, axis=0) / matrix.shape[0] for letter in LETTERS}).T


def highest_probability_sequence(df):
    return ''.join(df[column].idxmax() for column in df.columns)


def primer_plus_aptamer(sequences, ref):
    # start with primer and add 31 letters to the right (for the left primer)
    # allows us to find the wildcard for our aptamer
    subsequences = []
    fragment_length = APTAMER_LENGTH + PRIMER_LENGTH
    for sequence in sequences:
        min_position, max_position = 0, len(sequence) - (fragment_length + 1)
        pos = sequence.find(ref)
        if pos != -1 and (pos < max_position or pos >= min_position):
            if len(sequence[pos:pos + fragment_length]) == 51:
                subsequences.append(sequence[pos:pos + fragment_length])
    return np.array([list(s) for s in subsequences])


def plot_probabilities(df, ref):
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

    plt.savefig(f'new_data/output/plots/{ref}.png')
    plt.clf()
    plt.cla()
    plt.close()


def generate_aptamer_wildcard(df):
    columns = df.columns[df.columns.astype(int) >= 20]
    return ''.join([df[col].idxmax() if df[col].max() >= 0.5 else '*' for col in columns])


def wildcard_sequencies(sequencies, wildcard):
    length = len(wildcard)
    extracted_substrings = [seq[i: i + length] for seq in sequencies for i in
                            range(len(seq) - length + 1) if
                            not any(wildcard[j] != '*' and wildcard[j] != seq[i + j] for j in range(length))
                            if length <= len(seq)]
    return np.array([list(s) for s in extracted_substrings])


def define_primer(df):
    # extract the first primer candidate with the highest frequency
    primer_tmp = df['kmer'].values[0]
    primer_freq = df[df['kmer'] == primer_tmp]['frequency'].values[0]
    print(f'Primer candidate: {primer_tmp} with frequency = {primer_freq}')
    return primer_tmp, primer_freq


def aptamer_search():
    #seq_list = [str(rec.seq) for rec in SeqIO.parse(f'{DATA_DIR}/{INPUT_DATA}.{FILETYPE}', "fastq")]
    seq_list = merge_sequencies()
    print(f'Number of sequences in {INPUT_DATA}: {len(seq_list)}')
    # get list of 20-mers (primers)
    df = count_kmers(seq_list, 20)
    primer_tmp, primer_freq = define_primer(df)

    # primer plus 31-mer sequence
    primer_aptamer = primer_plus_aptamer(seq_list, primer_tmp)
    freq = calculate_probabilities(primer_aptamer)
    plot_probabilities(freq, primer_tmp)

    aptamer_wildcard = generate_aptamer_wildcard(freq)
    print(f'Aptamer wildcard: {aptamer_wildcard}')

    aptamers = wildcard_sequencies(seq_list, aptamer_wildcard)
    print(f'Number of sequences matching the aptamer wildcard: {len(aptamers)}')
    freq = calculate_probabilities(aptamers)
    plot_probabilities(freq, aptamer_wildcard)

    aptamer = highest_probability_sequence(freq)
    print(f'Aptamer: {aptamer}')


aptamer_search()
