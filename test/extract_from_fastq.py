from Bio import SeqIO
from Bio.SeqIO import QualityIO
import numpy as np
from Bio.Seq import Seq
from Bio import AlignIO
from Bio import pairwise2
from Bio.pairwise2 import format_alignment
from Bio import Align
from Bio.SeqRecord import SeqRecord
import os


file_clusters = 'output/final_clusters.tsv'
fastq_file_path = 'data/FAP38830_pass_barcode02_bcc3428d_0.fastq'

def extract_cluster(filename, cluster_id):
    ids = []
    with open(filename, 'r') as file:
        for line in file:
            columns = line.strip().split('\t')
            if len(columns) >= 2 and columns[0] == str(cluster_id):
                string = columns[1]
                runid_index = string.find('_runid')
                if runid_index != -1:
                    extracted_string = string[:runid_index]  # Extracts the substring until "_runid"
                    ids.append(extracted_string)
    return ids


# Function to extract sequences from a fastq file
def extract_sequences(file_path, id_list):
    extracted_sequences = {}

    with open(file_path, 'r') as file:
        lines = file.readlines()

        for i in range(0, len(lines), 4):
            sequence_id = lines[i].strip()[1:].split(' ')[0]  # Remove the '@' symbol from the sequence id
            sequence = lines[i+1].strip()

            if sequence_id in id_list:
                extracted_sequences[sequence_id] = sequence

    # Print the extracted sequences
    for sequence_id, sequence in extracted_sequences.items():
        print(f"Sequence ID: {sequence_id}")
        print(f"Sequence: {sequence}")
        print()

    return extracted_sequences

def get_list_of_sequences(file_path):
    return [str(rec.seq) for rec in SeqIO.parse(file_path, "fastq")]

def convert_to_fasta(input_filename, output_filename):
    SeqIO.convert(input_filename, "fastq", output_filename, "fasta")


def extract_sequences_from_fastq(folder_path):
    merged_sequences = []

    # Iterate over all files in the specified folder
    for filename in os.listdir(folder_path):
        file_path = os.path.join(folder_path, filename)

        # Check if the file is a FASTQ file
        if file_path.endswith(".fastq") or file_path.endswith(".fq"):
            with open(file_path, "r") as file:
                # Parse the FASTQ file using SeqIO
                sequences = [record.seq for record in SeqIO.parse(file, "fastq")]

                # Add the sequences to the merged list
                merged_sequences.extend(sequences)

    return merged_sequences


def sort_fastq_by_length_phred(input_file, output_file):
    records = []

    # Read the FASTQ file and store the sequences in a list
    with open(input_file, "r") as f:
        for record in SeqIO.parse(f, "fastq"):
            records.append(record)

    # Sort the records based on sequence length (longest first)
    records.sort(key=lambda r: len(r.seq), reverse=True)

    # Sort the records based on average Phred score (highest first)
    records.sort(key=lambda r: sum(r.letter_annotations["phred_quality"]) / len(r.letter_annotations["phred_quality"]),
                 reverse=True)

    # Write the sorted records to the output file
    with open(output_file, "w") as f:
        SeqIO.write(records, f, "fastq")


print(len(extract_sequences_from_fastq('new_data')))
print(len(extract_sequences_from_fastq('data')))
# Usage example
# input_file = "data/FAP38830_pass_barcode02_bcc3428d_0.fastq"
# output_file = "output.fasta"
#
# convert_to_fasta(input_file, output_file)

#sort_fastq_by_length_phred(input_file, output_file)

#print(get_list_of_sequences('data/FAP38830_pass_barcode02_bcc3428d_0.fastq'))

# input_id_list = extract_cluster(file_clusters, 0)
# extracted_sequences = extract_sequences(fastq_file_path, input_id_list)
# convert_to_fasta(extracted_sequences, 'cluster_0_shortest.fasta')

# alignment = Align.MultipleSeqAlignment([rec for rec in SeqIO.parse('cluster_0.fasta', "fasta")])
# print(alignment)
#
# from Bio.Align.Applications import MuscleCommandline
# cline = MuscleCommandline(input="cluster_0.fasta", out="alignment_cluster_0.txt")
# print(cline)

