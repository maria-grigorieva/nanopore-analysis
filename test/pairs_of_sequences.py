from collections import Counter
import pandas as pd
from extract_from_fastq import get_list_of_sequences

sequences = get_list_of_sequences('data/FAP38830_pass_barcode02_bcc3428d_0.fastq')
fragments_df = pd.read_excel('fragments.xlsx')
fragments_df.sort_values(by='occurences', inplace=True, ascending=False)
fragments = fragments_df['fragment'].values[0:50]

def check_fragments(long_string, fragment1, fragment2):
    start1 = long_string.find(fragment1)
    start2 = long_string.find(fragment2)

    if start1 == -1 or start2 == -1:
        return False  # Either fragment is not present in the long string

    end1 = start1 + len(fragment1)
    end2 = start2 + len(fragment2)

    if end1 <= start2 or end2 <= start1:
        return True  # Fragments do not overlap
    else:
        # Check reversed order
        start1_reversed = long_string.find(fragment1[::-1])
        start2_reversed = long_string.find(fragment2[::-1])

        if start1_reversed != -1 and start2_reversed != -1:
            end1_reversed = start1_reversed + len(fragment1)
            end2_reversed = start2_reversed + len(fragment2)

            if end1_reversed <= start2_reversed or end2_reversed <= start1_reversed:
                return True  # Fragments in reversed order do not overlap

    return False  # Fragments overlap in both orders


def find_non_overlapping_fragments(long_sequences, fragments):
    result = []

    # Iterate through each long sequence
    for sequence in long_sequences:

        # Generate all possible pairs of fragments
        for i in range(len(fragments)-1):
            for j in range(i+1, len(fragments)):
                if fragments[i] == fragments[j]:
                    break
                fragment1 = fragments[i]
                fragment2 = fragments[j]

                # Check if the fragments do not overlap
                if check_fragments(sequence, fragment1, fragment2):
                    result.append((fragment1, fragment2))

    # Count occurrences of fragment pairs in all sequences
    counts = Counter(tuple(pair) for sequence_pairs in result for pair in zip(sequence_pairs[::2], sequence_pairs[1::2]))

    # Convert counts to pandas DataFrame
    df_counts = pd.DataFrame(list(counts.items()), columns=['Fragment Pairs', 'Count'])

    return df_counts

counts = find_non_overlapping_fragments(sequences, fragments)

print(counts)

counts.to_excel('FragmentPairs.xlsx')