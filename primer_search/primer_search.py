from collections import Counter
from Bio import SeqIO
import glob
import re

def merge_sequencies():
    """
    Merge sequencies from a specified directory directory
    :return:
    """
    return [str(rec.seq) for file_path in glob.glob('../input_data/FAP38830_pass_barcode04_bcc3428d_0.fastq') for rec in
                SeqIO.parse(file_path, "fastq")]

def find_primer_sequences():
    # Join all sequences into one string for easy searching
    seq_list = merge_sequencies()
    all_sequences = ''.join(seq_list)

    # Find all 20-mers and count their frequencies
    all_20mers = re.findall(r'(?=(.{20}))', all_sequences)
    frequency_counter = Counter(all_20mers)

    # Select 20-mers with frequency >= 70% of the total number of sequencies
    min_frequency = 0.7 * len(seq_list)
    selected_20mers = [seq for seq, freq in frequency_counter.items() if freq >= min_frequency]

    # Test the selected 20-mers for compliance with the conditions
    primers = []
    for seq in selected_20mers:
        # Find the occurrences of the 20-mer in each initial sequence
        occurrences = [m.start() for m in re.finditer(seq, all_sequences)]

        for position in occurrences:
            # Check if there is at least one occurrence with a nearby 20-mer
            if any(abs(pos - position) <= 51 for pos in occurrences if pos != position):
                # Append the primer sequence to the results
                primers.append(all_sequences[position:position + 20])

    return primers


primers = find_primer_sequences()
print(primers)
