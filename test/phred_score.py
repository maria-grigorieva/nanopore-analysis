import os
import math
import pylab as plt
import numpy as np
import pandas as pd
import matplotlib.patches as patches
from Bio import SeqIO
import re


def plot_fastq_qualities(filename, ax=None, limit=20000, seq_len=100, boxplots=False, bins=20):

    fastq_parser = SeqIO.parse(filename, "fastq")

    # n_seq = len([str(rec.seq) for rec in fastq_parser])
    # print(n_seq)
    res=[]
    c=0
    finished = False

    while not finished:
        try:
            for record in fastq_parser:
                #if len(record.letter_annotations["phred_quality"]) >= seq_len:
                score=record.letter_annotations["phred_quality"][:seq_len-1]
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
    #
    # ax.set_xticks(np.arange(0, l, 5))
    # ax.set_xticklabels(np.arange(0, l, 5))
    ax.set_xlabel('position(bp)')
    ax.set_xlim((0,l))
    ax.set_ylim((0,40))
    ax.set_title('per base sequence quality')
    plt.show()

    return


def fastq_to_dataframe(filename, size=1000):
    """Convert fastq to dataframe.
        size: limit to the first reads of total size
        Returns: dataframe with reads
    """

    ext = os.path.splitext(filename)[1]
    if ext=='.fastq' or ext=='.gz':
        fastq_parser = SeqIO.parse(filename, "fastq")
    else:
        fastq_parser = SeqIO.parse(open(filename, "r"), "fastq")
    i=0
    res=[]
    for fastq_rec in fastq_parser:
        #print (fastq_rec.seq)
        i+=1
        if i>size:
            break
        res.append([fastq_rec.id, str(fastq_rec.seq)])
    df = pd.DataFrame(res, columns=['id','seq'])
    df['length'] = df.seq.str.len()
    return df

def normpdf(x, mean, sd):
    """sample a normal distribution at given point"""

    var = float(sd)**2
    denom = (2*math.pi*var)**.5
    num = math.exp(-(float(x)-float(mean))**2/(2*var))
    return num/denom

def plot_fastq_gc_content(filename, ax=None, limit=50000):

  from Bio.SeqUtils import GC
  if ax==None:
      f,ax=plt.subplots(figsize=(12,5))
  df = fastq_to_dataframe(filename, size=limit)
  gc = df.seq.apply(lambda x: GC(x))
  gc.hist(ax=ax,bins=150,color='black',grid=False,histtype='step',lw=2)
  ax.set_xlim((0,100))
  x=np.arange(1,100,.1)
  f = [normpdf(i, gc.mean(), gc.std()) for i in x]
  ax2=ax.twinx()
  ax2.plot(x,f)
  ax2.set_ylim(0,max(f))
  ax.set_title('GC content',size=15)
  plt.show()
  return


def select_ref_sequencies():

    seq_list = [{'sequence': str(rec.seq),
                'score': rec.letter_annotations["phred_quality"]}
                for rec in SeqIO.parse('../input_data/merged.fastq', "fastq")]

    reference = {'seq': 'GGCTTCTGG', 'start_pos': 51}

    PRIMER_TYPE = 'right'
    PRIMER_LENGTH = 20
    TOTAL_LENGTH = 71

    seq = reference['seq']
    left = reference['start_pos'] if PRIMER_TYPE == 'left' else reference['start_pos'] - PRIMER_LENGTH
    right = TOTAL_LENGTH - PRIMER_LENGTH - reference['start_pos'] - len(seq) \
            if PRIMER_TYPE == 'left' else TOTAL_LENGTH - reference['start_pos'] - len(seq)

    pattern = rf'.{{{left}}}{re.escape(seq)}.{{{right}}}'
    print(f'Positions: {left}-{seq}-{right}')
    matches = []
    scores = []
    for item in seq_list:
        sequence = item['sequence']
        score = item['score']
        interval = [{'start': m.start(0), 'end': m.end(0)} for m in re.finditer(pattern, sequence)]
        if len(interval) > 0:
            matches.append(sequence[interval[0]['start']:interval[0]['end']])
            scores.append(score[interval[0]['start']:interval[0]['end']])

    result_matches = np.array([list(s) for s in matches]) if len(matches) > 0 else None
    result_phred = np.array([list(s) for s in scores]) if len(matches) > 0 else None
    mean_values = result_phred.mean(axis=0)
    mean_df = pd.DataFrame(mean_values, columns=['Mean'])
    print(mean_df)


# select_ref_sequencies()

plot_fastq_qualities('../bigfile/bigfile.fastq',
                     seq_len=1000)

# plot_fastq_gc_content('../input_data/FAP38830_pass_barcode04_bcc3428d_0.fastq',
#                      limit=100)