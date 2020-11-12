import warnings
warnings.simplefilter(action='ignore')
from Bio import SeqIO
import numpy as np
import os, random, re, sys, glob
import pandas as pd
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

# Read and extract dataframe of fasta file
def read_fasta(filename):
    df = pd.DataFrame()
    reads = []
    ids = []
    with open(filename) as f:
        seqs = SeqIO.parse(f, 'fasta')
        for seq in seqs:
            reads.append(str(seq.seq))
            ids.append(seq.id)
    df['ID'] = ids
    df['seqs'] = reads
    return df

def assign_labels(df, multi=False):
    viral = []
    for i in df.ID:
        if('gi' in i):
            viral.append(1)
        elif((multi==True) and ('chr' in i)):
            viral.append(2)
        else:
            viral.append(0)
    print(df.shape)
    df['viral'] = viral
    return df

# transform a list of reads into a list of one-hot encoded vectors:
# (A,C,G,N,T) --> (00001, 00010, 00100, 01000, 10000)
def seqs2onehot(seqs):
    def one_hot_encode(seq):
        mapping = dict(zip("ACGNT", range(5)))
        seq = [mapping[i] if i in ['A', 'T', 'C', 'G', 'N'] else mapping['N'] for i in seq]
        return np.eye(5)[seq]
    onehotvecs = []
    for i in seqs:
        if(len(i) < 150):
            continue
        onehotvecs.append(np.array(one_hot_encode(i), dtype='bool'))
    return onehotvecs