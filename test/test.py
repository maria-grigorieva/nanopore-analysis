from Bio import SeqIO
from Bio.SeqIO import QualityIO
import numpy as np
from Bio.Seq import Seq
from Bio import AlignIO
from Bio import pairwise2
from Bio.pairwise2 import format_alignment
from Bio import Align
from Bio.SeqRecord import SeqRecord
import pandas as pd
import plotly.graph_objects as go

df = pd.read_csv('../results/barcode05_GGCTTCTGG_51_TREE_80/tree.csv')

# df = df[(df['path'] == 'CCAGCACCT_42')]
# df['source'] = df['prev_seq'].astype(str) + df['prev_start_pos'].astype(str)
# df['target'] = df['seq'].astype(str) + df['start_pos'].astype(str)
df['source'] = df['prev_seq'].astype(str).str[0] + df['prev_start_pos'].astype(str)
df['target'] = df['letter'].astype(str) + df['start_pos'].astype(str)

labels = list(set(np.concatenate((df['source'].values, df['target'].values), axis=0)))
# Assign each element with its index
indexed_arr = [(index, value) for index, value in enumerate(labels)]
# def indexing(row, indexed_arr):
#     index = next((index for index, value in indexed_arr if value == row['source']), None)
#     return index
df['source_index'] = df.apply(lambda row: next((index for index, value in indexed_arr if value == row['source']), None), axis=1)
df['target_index'] = df.apply(lambda row: next((index for index, value in indexed_arr if value == row['target']), None), axis=1)

# Initialize sources, targets, and values
sources = df['source_index'].values[::-1]
targets = df['target_index'].values[::-1]
values = df['hits'].values[::-1]

# Create Sankey diagram
fig = go.Figure(data=[go.Sankey(
    node=dict(
        pad=15,
        thickness=20,
        line=dict(color="black", width=0.5),
        label=labels
    ),
    link=dict(
        source=sources,
        target=targets,
        value=values
    )
)])

fig.update_layout(title_text="Sankey Diagram")
fig.show()


#
#
#
# fig = go.Figure(data=[go.Sankey(
#     node = dict(
#       pad = 15,
#       thickness = 20,
#       line = dict(color = "black", width = 0.5),
#       label = df['seq'].values,
#       color = "blue"
#     ),
#     link = dict(
#       source = df['prev_start_pos'].values, # indices correspond to labels, eg A1, A2, A1, B1, ...
#       target = df['start_pos'].values,
#       value = df['hits'].values
#   ))])
#
# fig.update_layout(title_text="Basic Sankey Diagram", font_size=10)
# fig.show()


# count = 0
# for rec in SeqIO.parse(filename, "fastq"):
#     count += 1
# print("%i reads" % count)
#
# good_reads = (
#     rec
#     for rec in SeqIO.parse(filename, "fastq")
#     if max(rec.letter_annotations["phred_quality"]) >= 20
# )
# count = SeqIO.write(good_reads, "good_quality.fastq", "fastq")
# print("Saved %i reads" % count)

# print(errors)
# sizes = [len(rec.seq) for rec in SeqIO.parse(filename, "fastq")]
# len(sizes), min(sizes), max(sizes)

# import matplotlib as mpl
# import matplotlib.pyplot as plt
#
# plt.hist(sizes, bins=20)
# plt.show()
#
# plt.hist(phred_qualities, bins=20)
# plt.show()
#
# plt.hist(errors, bins=20)
# plt.show()

import pylab
#
# for subfigure in [1, 2]:
#     filename = "data/FAP38830_pass_barcode02_bcc3428d_0.fastq"
#     pylab.subplot(1, 2, subfigure)
#     for i, record in enumerate(SeqIO.parse(filename, "fastq")):
#         if i >= 50:
#             break  # trick!
#         pylab.plot(record.letter_annotations["phred_quality"])
#     pylab.ylim(0, 45)
#     pylab.ylabel("PHRED quality score")
#     pylab.xlabel("Position")
# pylab.savefig("FAP38830_pass_barcode02_bcc3428d_0.png")
# print("Done")

# for i, record in enumerate(SeqIO.parse(filename, "fastq")):
#     if i >= 50:
#         break  # trick!
#     pylab.plot(record.letter_annotations["phred_quality"])
#     # plt.plot(record.letter_annotations["phred_quality"])
#     pylab.ylim(0,45)
#     pylab.ylabel("PHRED quality score")
#     pylab.xlabel("Position")
# pylab.show()

# sequences = [rec.seq for rec in SeqIO.parse(filename, "fastq")]
#
# print(sequences)
#
# alignments = pairwise2.align.globalxx(sequences[1], sequences[2])
# #print(alignments)
#
# for alignment in alignments:
#     print(format_alignment(*alignment))
#
#
# aligner = Align.PairwiseAligner()
# print(aligner)


# convert FASTQ to FASTA
# sequences = []
# for rec in SeqIO.parse(filename, "fastq"):
#     sequences.append(rec)
#
# print(sequences)
#
# records = []
# for rec in sequences:
#     records.append(SeqRecord(rec.seq,rec.id))
#
# SeqIO.write(records, "example.fasta", "fasta")
#
# sequences = [rec.seq for rec in SeqIO.parse(filename, "fastq")]
# longest_length = max(len(s) for s in sequences)
# print(longest_length)
# padded_sequences = [str(s).ljust(longest_length, '-') for s in sequences]
# records = (SeqRecord(Seq(s)) for s in padded_sequences)
# SeqIO.write(records, "example.fasta", "fasta")
#
# alignment = AlignIO.read("example.fasta", "fasta")
# print(alignment)
#
# from Bio.Align.Applications import MuscleCommandline
# cline = MuscleCommandline(input="example.fasta", out="msa.txt")
# print(cline)