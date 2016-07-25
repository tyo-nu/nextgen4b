############
# Imports
############

import numpy as np
import pandas as pd

from statsmodels.stats import proportion

import itertools
import re

import matplotlib.pyplot as plt

# Optional
import seaborn as sbn
sbn.set_style('white')
sbn.set_palette('muted', 10)


############
# Housekeeping
############

def define_words(alphabet, n_letters):
    return [''.join(x) for x in itertools.product(alphabet, repeat=n_letters)]

def replace_rare_base(word, rare_base='A'):
    return ''.join([c if (c is rare_base or c is '-') else 'X' for c in word])

def get_nonrare_words(word_list, rare_base='A'):
    # This only works for non-dash sequences right now.
    n = len(word_list[0])

    words = define_words([rare_base, 'X'], n)
    counts = {w: 0 for w in words}
    res = {w: re.compile(w) for w in words}

    for s in word_list:
        mod_s = replace_rare_base(s, rare_base=rare_base)
        for w in words:
            # try using regexes? might not be faster than simple comps though...
            if res[w].search(mod_s):
                counts[w] +=1

    return counts

def get_pos_cts(L, letter_order=['C','T','G','A']):
    n = len(L[0])
    cts = pd.DataFrame(np.zeros((n, len(letter_order))), columns=letter_order)

    for i in range(n):
        v = [w[i] for w in L]
        for c in letter_order:
            cts[c][i] = v.count(c)

    return cts


############
# Processing Functions
############

def process_sample(idx, data_path, prefix='words_', keep_del=False, keep_dash=True):
    """
    Read in a words file (specified by index and a datapath), outputs a list of
    the words, stripped of deletions and dashes if indicated.
    """
    with open('%s%s%i.txt' % (data_path, prefix, idx)) as f:
        data = [line.strip() for line in f.readlines()]
        if keep_dash:
            out_data = data
        elif keep_del:
            no_dash = [line for line in data[-1] if '-' not in line]
        else:
            no_dash = [line for line in data[-1] 
                       if '-' not in line and 'd' not in line]
    return out_data

def create_count_df(idxs, t_range, data_path, n, prefix='words_', rare_base='A'):
    """
    Take in a list of indicies and their timepoints, output a dataframe
    consisting of word counts at each timepoint. 
    """
    assert len(idxs) == len(t_range)

    counts = [get_nonrare_words(process_sample(idx, data_path, refix=prefix),
                                rare_base=rare_base)
              for idx in idxs]
    
    out_df = pd.DataFrame()
    out_df['Time'] = t_range
    for word in define_words(['A', 'X'], n):
        out_df[word] = [ct[word] for ct in counts_5]

    return out_df

def create_confint_df(count_df, ignored_cols=['Time']):
    """
    
    """
    words = [c for c in count_df.columns if c not in ignored_cols]
    ci_df = pd.DataFrame()

    for w in words:
        lb, ub = proportion.proportion_confint(counts_df[w], counts_df.sum(axis=1), method='jeffrey')
        ci_df['%s_lb' % w] = lb
        ci_df['%s_ub' % w] = ub

    return ci_df

def create_rate_df(count_df, ignored_cols=['Time']):
    """
    
    """
    select_col_df = count_df[[c for c in count_df.columns if c not in ignored_cols]]
    rate_df = count_df.div(select_col_df.sum(axis=1), axis=0)
    rate_df.Time = count_df.Time

    ci_df = create_confint_df(count_df, ignored_cols=ignored_cols)

    return rate_df, ci_df

def average_replicate_dfs(df1, df2):