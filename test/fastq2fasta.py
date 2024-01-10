from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

# convert FASTQ to FASTA
sequences = []
for rec in SeqIO.parse('/Users/maria/PyCharmProjects/nanopore-analysis/input_data/merged.fastq', "fastq"):
    sequences.append(rec)
records = []
for rec in sequences:
    records.append(SeqRecord(rec.seq,rec.id))

SeqIO.write(records, "/Users/maria/PyCharmProjects/nanopore-analysis/input_data/merged.fasta", "fasta")
#
# sequences = [rec.seq for rec in SeqIO.parse('/Users/maria/PyCharmProjects/nanopore-analysis/input_data/merged.fastq', "fastq")]
# longest_length = max(len(s) for s in sequences)
# print(longest_length)
# padded_sequences = [str(s).ljust(longest_length, '-') for s in sequences]
# records = (SeqRecord(Seq(s)) for s in padded_sequences)
# SeqIO.write(records, "/Users/maria/PyCharmProjects/nanopore-analysis/input_data/merged.fasta", "fasta")