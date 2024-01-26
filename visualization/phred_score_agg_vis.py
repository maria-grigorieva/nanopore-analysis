import pylab as plt
import numpy as np
import pandas as pd
import matplotlib.patches as patches
from Bio import SeqIO
import sys
import argparse
from pathlib import Path

def main(args):
    parser = argparse.ArgumentParser(description='Visualization of phred score for fastq files')
    args_info = [
        ['-i', '--input', str, 'Path to the fastq input file', 'input_file'],
        ['-l', '--limit', int, 'Limit number of sequences', 'limit', 10000],
        ['-sl', '--seq_length', int, 'Sequence length', 100],
        ['-bp', '--boxplots', bool, 'Draw boxplots', False],
        ['-b', '--bins', int, 'Number of bins at X axis', 20],
        ['-s', '--save', bool, 'Save to file', False]
    ]

    for arg_info in args_info:
        parser.add_argument(arg_info[0], arg_info[1], type=arg_info[2], help=arg_info[3], default=arg_info[4])

    args = parser.parse_args()

    filename = args.input
    limit = args.limit
    seq_len = args.seq_length
    boxplots = args.boxplots
    bins = args.bins
    save = args.save

    plot_fastq_qualities(filename,
                         None,
                         limit,
                         seq_len,
                         boxplots,
                         bins,
                         save)


def plot_fastq_qualities(filename,
                         ax=None,
                         limit=20000,
                         seq_len=100,
                         boxplots=False,
                         bins=20,
                         save=False):

    fastq_parser = SeqIO.parse(filename, "fastq")

    # n_seq = len([str(rec.seq) for rec in fastq_parser])
    # print(n_seq)
    res=[]
    c=0
    finished = False

    while not finished:
        try:
            for record in fastq_parser:
                if boxplots:
                    if len(record.letter_annotations["phred_quality"]) >= seq_len:
                        score=record.letter_annotations["phred_quality"][:seq_len-1]
                        res.append(score)
                else:
                    score = record.letter_annotations["phred_quality"][:seq_len - 1]
                    res.append(score)
            c += 1
            if c > limit:
                break
            finished = True
        except ValueError as e:
            print(e)
    print(f'Selected sequences {len(res)}')
    df = pd.DataFrame(res)
    l = len(df.T)+1

    if ax==None:
        f,ax=plt.subplots(figsize=(12,5))
    rect = patches.Rectangle((0,0),l,20,linewidth=0,facecolor='r',alpha=.4)
    ax.add_patch(rect)
    rect = patches.Rectangle((0,20),l,8,linewidth=0,facecolor='yellow',alpha=.4)
    ax.add_patch(rect)
    rect = patches.Rectangle((0,28),l,12,linewidth=0,facecolor='g',alpha=.4)
    ax.add_patch(rect)
    df.mean().plot(ax=ax,c='black')
    if boxplots:
        boxprops = dict(linestyle='-', linewidth=1, color='black')
        df.plot(kind='box', ax=ax, grid=False, showfliers=False,
                color=dict(boxes='black',whiskers='black')  )
    N = bins  # Maximum number of values displayed on the x-axis

    # Calculate the step size based on the length of the axis and the maximum number of values
    step_size = max(int(l / N), 1)

    # Set the x-ticks and labels with the calculated step size
    ax.set_xticks(np.arange(0, l, step_size))
    ax.set_xticklabels(np.arange(0, l, step_size))

    ax.set_xlabel('position(bp)')
    ax.set_xlim((0,l))
    ax.set_ylim((0,40))
    ax.set_title('per base sequence quality')
    if not save:
        plt.show()
    else:
        plt.savefig(f'{Path(filename).stem}_phred_score.png')
    return

if __name__ == "__main__":
    main(sys.argv[1:])