from Bio import SeqIO
from Bio.SeqIO import QualityIO
import numpy as np
from Bio.Seq import Seq
from Bio import AlignIO
from Bio import pairwise2
from Bio.pairwise2 import format_alignment
from Bio import Align
from Bio.SeqRecord import SeqRecord

filename = "data/FAP38830_pass_barcode02_bcc3428d_0.fastq"

count = 0
for rec in SeqIO.parse(filename, "fastq"):
    count += 1
print("%i reads" % count)

good_reads = (
    rec
    for rec in SeqIO.parse(filename, "fastq")
    if max(rec.letter_annotations["phred_quality"]) >= 20
)
count = SeqIO.write(good_reads, "good_quality.fastq", "fastq")
print("Saved %i reads" % count)

errors = []
phred_qualities = []
for rec in SeqIO.parse(filename, "fastq"):
    phred_quality = np.mean(rec.letter_annotations["phred_quality"])
    phred_qualities.append(phred_quality)
    error = 10 ** (- phred_quality / 10)
    errors.append(error)
    # print(rec.seq)

sizes = [len(rec.seq) for rec in SeqIO.parse(filename, "fastq")]
len(sizes), min(sizes), max(sizes)

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
sequences = []
for rec in SeqIO.parse(filename, "fastq"):
    sequences.append(rec)

print(sequences)

records = []
for rec in sequences:
    records.append(SeqRecord(rec.seq,rec.id))

SeqIO.write(records, "example.fasta", "fasta")

sequences = [rec.seq for rec in SeqIO.parse(filename, "fastq")]
longest_length = max(len(s) for s in sequences)
print(longest_length)
padded_sequences = [str(s).ljust(longest_length, '-') for s in sequences]
records = (SeqRecord(Seq(s)) for s in padded_sequences)
SeqIO.write(records, "example.fasta", "fasta")

alignment = AlignIO.read("example.fasta", "fasta")
print(alignment)

from Bio.Align.Applications import MuscleCommandline
cline = MuscleCommandline(input="example.fasta", out="msa.txt")
print(cline)