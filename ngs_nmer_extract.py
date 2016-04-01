# Ted Cybulski
# 4/1/16

# Implemented to analyze Alex's 6mer dataset.
# Takes in a list of indices to define the sites that make up a motif
# Takes in another list of indices that define the sites we want to count at
# Finds unique motifs, then counts the incidence of nucleotides at each of the count sites

from Bio import SeqIO

import pandas as pd
import sys


def define_words(alphabet, n):
    return [''.join(x) for x in itertools.product(alphabet, repeat=n)]
    
def process_fa(f_name, mot_idxs, ct_idxs, rare_base='A', letter_order=['C','G','T','A']):
    seqs = SeqIO.parse(f_name, 'fasta')
    mot = []
    ct_nts = []
     
    for s in seqs:
        m = ''.join([s[idx] for idx in mot_idxs])
        if '-' not in m and rare_base not in m:
            mot.append(m)
            ct_nts.append([s[idx] for idx in ct_idxs])
    
    unique_mots = list(set(mot))
    mot_cts = [mot.count(m) for m in unique_mots]
    
    # Build df
    df = pd.DataFrame({'motif': unique_mots, 'total': mot_cts})
    for c_idx, i in zip(ct_idxs, range(len(ct_idxs))):
        site_letters = [l[i] for l in ct_nts]
        for letter in letter_order:
            site_mot_ct = [[site_letters[i] for i in range(len(mot)) if mot[i] == m].count(letter) for m in unique_mots]
            df['%s_%i_counts' % (letter, c_idx)] = site_mot_ct
    
    return df



if __name__ == '__main__':
    # 