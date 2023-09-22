from Bio import SeqIO
import os
file_clusters = 'clustering_results/final_clusters.tsv'
fastq_file_path = 'new_data/FAP38830_pass_barcode05_bcc3428d_0.fastq'

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
