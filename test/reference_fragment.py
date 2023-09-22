from extract_from_fastq import get_list_of_sequences, extract_sequences_from_fastq
from collections import Counter
import pandas as pd
import numpy as np

file_name = 'data/FAP38830_pass_barcode02_bcc3428d_0.fastq'

# seq_list = extract_sequences_from_fastq('data')
seq_list = get_list_of_sequences(file_name)
print(f'Number of sequences: {len(seq_list)}')

def get_overlapped_seq_fragments(seq, fragment_size=9):
    return [seq[i:i+fragment_size] for i in range(len(seq) - fragment_size + 1)]

def get_all_kmers(seq_list, kmer_length=9):
    kmers = list(set(seq[i:i + kmer_length] for seq in seq_list for i in range(len(seq) - kmer_length + 1)))
    print(f'Number of unique {kmer_length}-fragments: {len(kmers)}')
    return kmers

def get_fragment_positions(fragment, seq_list):
    positions = []
    for seq in seq_list:
        if fragment in seq:
            fragment_positions = [pos for pos in range(len(seq)) if seq[pos:].startswith(fragment)]
            positions.extend([pos for pos in fragment_positions])
    return positions

unique_fragments = get_all_kmers(seq_list)

# Count the number of initial sequences where each unique fragment exists
fragment_count = Counter(fragment for fragment in unique_fragments
                         for seq in seq_list if fragment in seq)

# Sort the result in ascending order based on occurrences
sorted_fragments = sorted(fragment_count.items(), key=lambda x: x[1])

# Print the sorted result and fragment positions
for fragment, count in sorted_fragments:
    print(f"Fragment: {fragment}, Occurrences: {count}")
    # positions = get_fragment_positions(fragment, seq_list)
    # print(f"Fragment: {fragment}, Occurrences: {count}, Positions: {positions}")

# Save the results to Excel file
# Create a DataFrame from the data
df = pd.DataFrame([{'fragment': fragment,
                    'occurences': count} for fragment, count in sorted_fragments])
# df = pd.DataFrame([{'fragment': fragment,
#                     'occurences': count,
#                     'positions': get_fragment_positions(fragment, seq_list),
#                     'avg_position':np.mean(get_fragment_positions(fragment, seq_list))} for fragment, count in sorted_fragments])

# Save the DataFrame to an Excel file
df.to_excel('fragments.xlsx', index=False)






