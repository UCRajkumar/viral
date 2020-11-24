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
# (A,C,G,T) --> (0001, 0010, 0100, 1000)
def seqs2onehot(seqs):
    def one_hot_encode(seq):
        mapping = dict(zip("ACGT", range(4)))
        seq = [mapping[i] if i in ['A', 'T', 'C', 'G'] else mapping['A'] for i in seq]
        return np.eye(4)[seq]
    onehotvecs = []
    for i in seqs:
        if(len(i) < 150):
            continue
        onehotvecs.append(np.array(one_hot_encode(i), dtype='bool'))
    return onehotvecs

# reference: https://www.geeksforgeeks.org/edit-distance-dp-5/
def editDistDP(str1, str2, m, n): 

    dp = [[0 for x in range(n + 1)] for x in range(m + 1)] 
  
    for i in range(m + 1): 
        for j in range(n + 1): 
            if i == 0: 
                dp[i][j] = j    # Min. operations = j 
            elif j == 0: 
                dp[i][j] = i    # Min. operations = i 
  
            elif str1[i-1] == str2[j-1]: 
                dp[i][j] = dp[i-1][j-1] 
            else: 
                dp[i][j] = 1 + min(dp[i][j-1],        # Insert 
                                   dp[i-1][j],        # Remove 
                                   dp[i-1][j-1])    # Replace 
  
    return dp[m][n] 

def onehotdecode(arr):
    str_builder = ''
    for i in arr:
        val = np.argmax(i)
        if(val == 0):
            str_builder += 'A'
        elif(val == 1):
            str_builder += 'C'
        if(val == 2):
            str_builder += 'G'
        if(val == 3):
            str_builder += 'T'
    return str_builder