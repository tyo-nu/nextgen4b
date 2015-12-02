import numpy as np
import pandas as pd
from matplotlib import pyplot as plt

import logging, sys

from Bio import AlignIO
from Bio import SeqIO

import alignfilter

#####################
# Dataframe creation and manipulation
#####################

def getPositionalMisinc(seqs, template, n, letterorder=['C', 'A', 'T', 'G']):
    mat = np.zeros([len(letterorder), len(letterorder), len(seqs)])
    
    for i in range(len(seqs)):
        if seqs[i][n] in letterorder and template[n] in letterorder:
            mat[letterorder.index(template[n])][letterorder.index(seqs[i][n])][i] += 1
            
    return mat

def getAllPositionMisincs(seqs, template, letterorder=['C', 'A', 'T', 'G']):
    mat = np.zeros([len(letterorder), len(letterorder), len(template)])
    
    for i in range(len(template)):
        posMat = getPositionalMisinc(seqs, template, i, letterorder=letterorder)
        mat[:,:,i] = np.sum(posMat, axis=2)
    
    return mat
    
def posMtoPandas(m, letterorder=['C', 'A', 'T', 'G']):
    # Generate column labels
    labels = []
    for i in range(len(letterorder)):
        for j in range(len(letterorder)):
            labels.append(letterorder[i]+'->'+letterorder[j])
    
    # Populate dataframe with Z-columns of matrix
    df = pd.DataFrame(data=np.zeros([m.shape[2], len(labels)]), columns=labels)
    for i in range(len(letterorder)):
        for j in range(len(letterorder)):
            df[letterorder[i]+'->'+letterorder[j]] = m[i,j,:]
    
    return df
    
def addSequenceColumn(df, template):
    tS = pd.Series(data=list(template), name='sequence', index=df.index)
    df['sequence'] = tS
    return df

############
# Main Routine Helper Functions
############
    
def processNeedle(needleObj, tempSeq):
    # filter based on alignment score.
    seqs = alignfilter.cullAlignments(needleObj)    
    return doAnalysis(seqs, tempSeq)
    
def doAnalysis(seqs, template):
    logging.info('Started Analysis')
    m = getAllPositionMisincs(seqs, template)
    df = addSequenceColumn(posMtoPandas(m), template)
    logging.info('Finished Analysis')
    return df

############
# Main Routine
############

if __name__ == '__main__':
    if len(sys.argv) > 1:
        ofilen = sys.argv[1]
    
    alnData = AlignIO.parse(open(ofilen), "emboss")
    tempData = list(SeqIO.parse(open('temptemplate.fa'), 'fasta'))[0]
    
    df = processNeedle(alnData, tempData)
    of = open('_misinc_data.csv', 'w')
    df.to_csv(of)
    of.close()
    
    