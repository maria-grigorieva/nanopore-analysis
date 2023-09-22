import pandas as pd
from difflib import SequenceMatcher


df = pd.read_excel('new_data/output/70_mer_frequencies.xlsx')
df = df.sort_values(by='frequency', ascending=False)
seq_list = df[df['frequency']>5]['kmer'].tolist()

def longest_common_sequence(sequences):
    commons = []
    for i in range(0,len(sequences)-1):
        seqA = sequences[i]
        seqB = sequences[i+1]
        seqMatch = SequenceMatcher(None, seqA, seqB)
        match = seqMatch.find_longest_match(0, len(seqA), 0, len(seqB))
        if (match.size != 0):
            longest_sequence = seqA[match.a:match.a + match.size]
            if len(longest_sequence) > 31:
                commons.append(longest_sequence)
            else:
                break

    return commons

def deep_search(df):
    common = "A" * len(df['kmer'].values[0])
    frequencies = df['frequency'].unique()
    for i in frequencies:
        if len(common) > 31:
            seq_list = df[df['frequency'] >= i]['kmer'].tolist()
            common = longest_common_sequence(seq_list)
    return common

# Example usage
# result = deep_search(df)
result = longest_common_sequence(seq_list)
print("Longest common substring:", result)


