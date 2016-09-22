import logging
import sys

import numpy as np
import pandas as pd
import tqdm
import yaml
from Bio import SeqIO

from ..process.filter import cull_alignments


#####################
# Dataframe creation and manipulation
#####################

def get_positional_misinc(seqs, template, n, letterorder=['C', 'A', 'T', 'G']):
    mat = np.zeros([len(letterorder), len(letterorder), len(seqs)])
    
    for i in range(len(seqs)):
        if seqs[i][n] in letterorder and template[n] in letterorder:
            mat[letterorder.index(template[n])][letterorder.index(seqs[i][n])][i] += 1
            
    return mat

def get_all_position_misincs(seqs, template, letterorder=['C', 'A', 'T', 'G']):
    mat = np.zeros([len(letterorder), len(letterorder), len(template)])
    for i in range(len(template)):
        posMat = get_positional_misinc(seqs, template, i, letterorder=letterorder)
        mat[:,:,i] = np.sum(posMat, axis=2)
    return mat
    
def pos_mat_to_df(m, letterorder=['C', 'A', 'T', 'G']):
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
    
def add_sequence_column(df, template):
    tS = pd.Series(data=list(template), name='sequence', index=df.index)
    df['sequence'] = tS
    return df

############
# Main Routine Helper Functions
############
    
def do_analysis(seqs, template):
    m = get_all_position_misincs(seqs, template)
    df = add_sequence_column(pos_mat_to_df(m), template)
    return df

############
# Main Routines
############

def analyze_all_experiments(yf_name, data_dir='./'):
    """
    Given a folder of aligned fasta files from `filter`, output the old 
    misinc_data.csv files, along with a summary .csv of misincorporations
    at a given site.
    """
    with open(yf_name) as expt_f:
        expt_yaml = yaml.load(expt_f) # Should probably make this a class at some point...
    runs = expt_yaml['ngsruns']
    for run in tqdm.tqdm(runs.keys()):
        expts = runs[run]['experiments']
        for expt in expts:
            analyzed_data_fname = '%s_%s_misinc_data.csv' % (expt, run)
            template = expt_yaml['experiments'][expt]['template_seq']
            aln_seqs = list(SeqIO.parse('aln_seqs_%s_%s.fa' % (run, expt),
                                        'fasta'))
            data = do_analysis(aln_seqs, template)
            
            # Save dataframe
            with open(analyzed_data_fname, 'w') as of:
                data.to_csv(of)

