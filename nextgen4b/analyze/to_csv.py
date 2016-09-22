import argparse
import os
import sys

import numpy as np
import pandas as pd
from statsmodels.stats import proportion

##############
# CSV Loading Routines
##############
   
def get_exp_run_data(fname):
    tokens = fname.split('_')
    # (Experiment, Run)
    return (tokens[0], tokens[1])

def get_csv_file_names(directory='.', csv_suffix='misinc_data.csv'):
    csv_fnames = [f for f in os.listdir(directory) if os.path.isfile(f) and f.endswith(csv_suffix)]
    return csv_fnames
    
#####################
# Simple Stats 
#####################

def get_stats(df, cutoff=1):
    data = df[[c for c in df.columns if not c == 'sequence']]
    total_n = data.sum(axis=1) # All non-corresponding data should be zero
    correct_n = data[[n+'->'+n for n in ['A', 'C', 'T', 'G']]] # Get the columns that correspond to correct incorporations
    misinc_n = total_n - correct_n.sum(axis=1) # Same as above.
    
    # 
    rate = misinc_n / total_n
    lb, ub = proportion.proportion_confint(misinc_n, total_n, method='jeffrey')
    
    # Assemble output dataframe
    simp_df = pd.DataFrame()
    simp_df['rate'] = rate
    simp_df['lb'] = lb
    simp_df['ub'] = ub
    simp_df['n'] = total_n
    simp_df['sequence'] = df.sequence
    
    return simp_df
    
def write_all_simple_misinc(directory='.', letterorder=['C', 'A', 'T', 'G']):
    csv_fnames = get_csv_file_names(directory=directory)
    
    for fname in csv_fnames:
        df = pd.read_csv(fname, index_col=0)
        simple_df = get_stats(df)
        simple_df.to_csv('simple_'+fname)
        
#####################
# Single Position Stats
#####################

def get_pos_stats(df, nIdx, cutoff=1, expt=1, letterorder=['C', 'A', 'T', 'G']):
    # Get row of interest
    data = df[[c for c in df.columns if not c == 'sequence' and not c == 'index']].iloc[nIdx]
    nt = df['sequence'].iloc[nIdx]
    total_n = float(data.sum())
    
    # Set up dataframe
    ntCols = ['N->'+c for c in letterorder] + ['N->!N']
    outsCols = ['ct', '%', '%lb', '%ub']
    cols = [x+'_'+out for x in ntCols for out in outsCols] + ['total_n', 'sequence']
    out_df = pd.DataFrame(index=[expt], columns=cols)
    
    out_df['sequence'] = nt
    out_df['total_n'] = total_n   
    
    # Do individual nucleotide stats
    for n in letterorder:
        ct = data[nt+'->'+n]
        rate = ct / total_n
        lb, ub = proportion.proportion_confint(ct, total_n, method='jeffrey')
        
        out_df['N->'+n+'_ct'] = ct
        out_df['N->'+n+'_%'] = rate
        out_df['N->'+n+'_%lb'] = lb
        out_df['N->'+n+'_%ub'] = ub
    
    # Do aggregate misincorporation stats
    misinc_n = total_n - out_df['N->%c_ct' % nt]
    lb, ub = proportion.proportion_confint(misinc_n, total_n, method='jeffrey')
    out_df['N->!N_ct'] = misinc_n
    out_df['N->!N_%'] = misinc_n / total_n
    out_df['N->!N_%lb'] = lb
    out_df['N->!N_%ub'] = ub
    
    return out_df
    
def write_all_pos_stats(nIdx, directory='.', outfile='summary_all.csv',
                     letterorder=['C', 'A', 'T', 'G']):
    csv_fnames = get_csv_file_names(directory=directory)
    # expts = [int(getExpRunData(fname)[0][3:]) for fname in csv_fnames]
    
    # Set ids to be run first, then experiment
    f_ids = ['.'.join(reversed([str(x) for x in get_exp_run_data(fname)])) for fname in csv_fnames]
    
    # eeek this is a hack... wish there were a way to append using the correct columns
    # rather than just use .values and assume the columns are in the same order.
    ntCols = ['N->'+c for c in letterorder] + ['N->!N']
    outsCols = ['ct', '%', '%lb', '%ub']
    cols = [x+'_'+out for x in ntCols for out in outsCols] + ['total_n', 'sequence']
    
    out_df = pd.DataFrame(index=sorted(f_ids), columns=cols)
   
    for fname, f_id in zip(csv_fnames, f_ids):
        row = get_pos_stats(pd.read_csv(fname, index_col=0),
                            nIdx, expt=f_id, letterorder=letterorder)
        out_df.ix[f_id] = row.values

    # Sort by experiment number
    # NOTE: this is hard-coded for 
    out_df['sort_col'] = out_df.index.map(lambda x: int(x.split('.')[0][3:]))
    out_df.sort_values('sort_col', inplace=True)
    out_df.drop('sort_col', axis=1, inplace=True)
    
    out_df.to_csv(outfile)