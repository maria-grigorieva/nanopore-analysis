import pandas as pd
from extract_from_fastq import get_list_of_sequences
from collections import Counter
from itertools import combinations
from functools import reduce

sequences = get_list_of_sequences('data/FAP38830_pass_barcode02_bcc3428d_0.fastq')
fragments_df = pd.read_excel('fragments.xlsx')
fragments_df.sort_values(by='occurences', inplace=True, ascending=False)
fragments = fragments_df['fragment'].values[0:20]

def get_overlapped_seq_fragments(seq, fragment_size):
    return [seq[i:i+fragment_size] for i in range(len(seq) - fragment_size + 1)]

fragment_size = 9
max_num_fragments = 5  # Maximum number of fragments

# Generate overlapping fragments for different numbers of fragments
overlapped_fragments = []
for num_fragments in range(2, max_num_fragments + 1):
    fragments = get_overlapped_seq_fragments(sequences, fragment_size)
    combinations_fragments = list(combinations(fragments, num_fragments))
    overlapped_fragments.extend(combinations_fragments)

# Get unique fragments
unique_fragments = list(set(overlapped_fragments))
print(f'Number of unique fragments: {len(unique_fragments)}')

# Count the number of initial sequences where each unique fragment exists
fragment_counts = Counter(fragment for fragment in unique_fragments
                          for seq in sequences if fragment in seq)

# Sort the result in ascending order based on occurrences
sorted_fragments = sorted(fragment_counts.items(), key=lambda x: x[1])

# Print the sorted result and fragment positions
for fragment, count in sorted_fragments:
    print(f"Fragment: {fragment}, Occurrences: {count}")

