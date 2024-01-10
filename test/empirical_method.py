import pandas as pd
import matplotlib.pyplot as plt
import glob
from Bio import SeqIO
from Bio.Seq import Seq

import os
# from extract_from_fastq import get_list_of_sequences

sequences = [str(rec.seq) for rec in SeqIO.parse('../input_data/merged.fastq', "fastq")]
from collections import Counter
import re

def get_all_kmers(seq_list, kmer_length, step=1):
    kmers = [seq[i:i + kmer_length] for seq in seq_list for i in range(0, len(seq) - kmer_length + 1, step)]
    print(f'Number of unique {kmer_length}-fragments: {len(kmers)}')
    d = dict(Counter(kmers))
    return {key: value for key, value in d.items() if value >= 3000}

def find_equal_substrings(strings):

    candidates = get_all_kmers(strings, 8)
    print(candidates)
    found_substrings = []
    for c in candidates:
        print(c)
        for s in strings:
            matches = re.findall(rf'{re.escape(c)}.{{31}}{re.escape(c)}', s)
            found_substrings.extend(matches)
    substring_frequency = dict(Counter(found_substrings))
    return substring_frequency

# subs = find_equal_substrings(sequences)
# print(len(subs))
# print(subs)
# sorted_data = dict(sorted(subs.items(), key=lambda x: x[1], reverse=True))
# print(sorted_data)




DATA_DIR = 'new_data'
OUTPUT_DIR = 'new_data/output'
INPUT_DATA = 'FAP38830_pass_barcode04_bcc3428d_0'
FILETYPE = 'fastq'
OUTPUT_FILE = f'{INPUT_DATA}.xlsx'

# get k-mers frequences
def count_kmers(sequences, kmer_length=70, step=1):

    kmer_frequencies = {}
    # Process each sequence
    for sequence in sequences:
        # counted_kmers = set()  # Set to keep track of counted kmers in the current sequence
        # Iterate over each k-mer in the sequence
        for i in range(0, len(sequence) - kmer_length + 1, step):
            kmer = sequence[i:i + kmer_length]
            # if kmer not in counted_kmers:
            if kmer not in kmer_frequencies:
                kmer_frequencies[kmer] = 0
            kmer_frequencies[kmer] += 1
                # counted_kmers.add(kmer)
    df = pd.DataFrame(kmer_frequencies.items(), columns=['kmer', 'frequency'])
    df = df.sort_values('frequency', ascending=False)  # Sorting in descending order
    return df


# Save k-mer frequencies to Excel file
# def save_kmers(df, input_data):
#     kmer_len = len(df['kmer'].values[0])
#     output_filename = f'{OUTPUT_DIR}/{input_data}_{kmer_len}_mer_freq.xlsx'
#     df.to_excel(output_filename)


def plot_frequencies_dist(df, output_path):
    kmer_len = len(df['kmer'].values[0])
    output_file = f'{output_path}/freq_hist_{kmer_len}_mers.png'
    df['frequency'].plot(kind='hist', bins=10)
    plt.xlabel('Frequency')
    plt.ylabel('Count')
    plt.title(f'Distribution of Frequency for {kmer_len}-mers')
    plt.savefig(output_file)


def find_longest_substring(strings):
    from difflib import SequenceMatcher
    aptamers = []
    start = strings[0]
    for i in strings[1:]:
        match = SequenceMatcher(None, start, i).find_longest_match()
        print(match)  # -> Match(a=0, b=15, size=9)
        if match.size >= 31:
            start = start[match.a:match.a + match.size]
            if match.size == 31:
                aptamers.append(start)
            print(start[match.a:match.a + match.size])
        # print(i[match.b:match.b + match.size])
    print(set(aptamers))
    return list(set(aptamers))[0]

def get_primers(strings, aptamer):
    from collections import Counter
    left = []
    right = []

    for s in strings:
        if aptamer in s:
            start = s.find(aptamer)
            end = start + len(aptamer) - 1
            left_primer = s[start-20:start]
            right_primer = s[end:end+20]
            if len(left_primer) == 20:
                left.append(left_primer)
            if len(right_primer) == 20:
                right.append(right_primer)
    print(Counter(left))
    print(Counter(right))

def search_aptamers(seq_list):
                    # input_data, output_file):

    frequencies_70 = count_kmers(seq_list, 70)
    #save_kmers(frequencies_70, input_data)
    #
    # df = pd.read_excel('new_data/output/70_mer_frequencies.xlsx')
    df = frequencies_70.sort_values(by='frequency', ascending=False)
    seq_list = df[df['frequency']>=100]['kmer'].tolist()
    #df[df['frequency']>=10].to_excel('71-mers.xlsx')

    longest = find_longest_substring(df[df['frequency']>=5]['kmer'].tolist())
    print(longest)
    get_primers(df[df['frequency']>=5]['kmer'].tolist(),longest
                )
    # result_df = count_kmers(seq_list, 31)
    # print(result_df)

    # result_df.to_excel(f'{OUTPUT_DIR}/{output_file}')


search_aptamers(sequences)


# Example usage
string_list = ["abcdefg", "bcde", "cdef"]
longest_substring = find_longest_substring(string_list)
print("Longest common substring:", longest_substring)


def search_primers(seq_list, input_data, output_file):

    frequencies_70 = count_kmers(seq_list, 70)
    # save_kmers(frequencies_70, input_data)
    df = frequencies_70.sort_values(by='frequency', ascending=False)
    seq_list = df[df['frequency']>=5]['kmer'].tolist()

    result_df = count_kmers(seq_list, 20)
    print(result_df)

    result_df.to_excel(f'new_data/output/primers/{output_file}')


def fastq_iterator(mode = 'aptamer'):
    file_pattern = '*.fastq'

    # Iterate over the FASTQ files in the directory
    for file_path in glob.glob(f'{DATA_DIR}/{file_pattern}'):
        seq_list = [str(rec.seq) for rec in SeqIO.parse(file_path, "fastq")]
        # seq_list = get_list_of_sequences(file_path)
        input_data = os.path.basename(file_path).split('.')[0]
        output_file = f'{input_data}.xlsx'
        if mode == 'aptamer':
            search_aptamers(seq_list, input_data, output_file)
        elif mode == 'primer':
            search_primers(seq_list, input_data, output_file)


def complementary_sequence(sequence):
    seq = Seq(sequence)
    return seq.complement()


def get_stats():
    file_pattern = '*.fastq'

    initial_data = {
        'aptamer': Seq('TAGGGAAACACGATAGAATCCGAACAGCACC'),
        'aptamer_complementary': complementary_sequence('TAGGGAAACACGATAGAATCCGAACAGCACC'),
        'primer_left': Seq('CTCCTCTGACTGTAACCACG'),
        'primer_left_complementary': complementary_sequence('CTCCTCTGACTGTAACCACG'),
        'primer_right': Seq('GGCTTCTGGACTACCTATGC'),
        'primer_right_complementary': complementary_sequence('GGCTTCTGGACTACCTATGC'),
        'my_primer': Seq('TTACGTATTGCTAAGGTTAA'),
        'my_primer_complementary': complementary_sequence('TTACGTATTGCTAAGGTTAA')
    }

    print(initial_data)

    # Function to search for a sequence in a record
    def find_sequence_occurence(fastq_file, sequence, sequence_name):
        file_name = os.path.basename(file_path).split('.')[0]
        count = 0
        n_seq = 0
        for record in SeqIO.parse(fastq_file, "fastq"):
            n_seq += 1
            if sequence in record.seq:
                count += 1
        return {'file_name': file_name,
                'number_of_sequencies': n_seq,
                'sequence_name': sequence_name,
                'occurencies': count}

    result = []
    # Iterate over the FASTQ files in the directory
    for file_path in glob.glob(f'{DATA_DIR}/{file_pattern}'):
        for key, value in initial_data.items():
            result.append(find_sequence_occurence(file_path, value, key))

    df = pd.DataFrame(result)
    pivot_df = df.pivot(index=['file_name','number_of_sequencies'], columns='sequence_name', values='occurencies')
    print(pivot_df)
    pivot_df.to_excel('new_data/output/stats_v1.xlsx')

# get_stats()
# fastq_iterator('primer')

def excel_to_fasta_converter(input_file, output_file):
    df = pd.read_excel(input_file)
    sequences = df['kmer'].values.tolist()
    with open(output_file, 'w') as file:
        for i, s in enumerate(sequences):
            file.write(f'>sequence_{i+1}\n')
            file.write(s + '\n')

#excel_to_fasta_converter('new_data/output/FAP38830_pass_barcode04_bcc3428d_0_70_mer_freq.xlsx','new_data/output/kmer_31.fasta')


# seq_list = get_list_of_sequences(f'{DATA_DIR}/{INPUT_DATA}.{FILETYPE}')
# print(f'Number of sequences: {len(seq_list)}')

# def count_kmers(sequences, kmers):
#     num_kmers = len(kmers[0])
#
#     # Initialize the Count-min Sketch data structure
#     width = num_kmers  # Width of the sketch (depends on the number of kmers)
#     depth = 4  # Depth of the sketch (can be adjusted for desired accuracy)
#     sketch = np.zeros((depth, width), dtype=np.uint32)
#
#     # Generate hash functions for each row of the sketch
#     hash_functions = [mmh3.hash128] * depth
#
#     # Process each sequence
#     for sequence in sequences:
#         # Iterate over each k-mer in the sequence
#         for i in range(len(sequence) - num_kmers + 1):
#             kmer = sequence[i:i + num_kmers]
#
#             # Calculate the hash value for the k-mer
#             hashes = [hash_fn(kmer) for hash_fn in hash_functions]
#
#             # Increment the corresponding cells in the sketch
#             for j in range(depth):
#                 sketch[j, hashes[j] % width] += 1
#
#     # Calculate the estimated frequencies of each k-mer
#     frequencies = []
#     for kmer in kmers:
#         hashes = [hash_fn(kmer) for hash_fn in hash_functions]
#         counts = [sketch[j, hashes[j] % width] for j in range(depth)]
#         estimated_frequency = min(counts)  # Choose the minimum count across all cells
#         frequencies.append(estimated_frequency)
#
#     return frequencies
#
# freqs = count_kmers(seq_list, result_df['kmer'].values.tolist())
# result_df['total_frequencies'] = freqs
# print(freqs)





# kmers = get_all_kmers(seq_list, 70)
# frequencies = count_kmers(seq_list, kmers)
# print(frequencies)
# save_kmers(frequencies, 'new_data/output')
# plot_frequencies_dist(frequencies, 'new_data/output/plots')

# for k in range(5,55):
#     kmers = get_all_kmers(seq_list, kmer_length=k)
#     frequencies = count_kmers(seq_list, kmers)
#     save_kmers(frequencies, 'new_data/output')
#     plot_frequencies_dist(frequencies, 'new_data/output/plots')


#
# df = count_kmers(sequences, kmer_length=31, step=1)
# df.to_excel('31_mers.xlsx')