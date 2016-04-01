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
    
def process_fa(f_name, mot_idxs, ct_idxs):
     seqs = SeqIO.parse(f_name)
     mot = []
     
     for s in seqs:
        mot.append(''.join([s[idx] for idx in mot_idxs]))
        ct_nts = [s[idx] for idx in ct_idxs]
        
        
    return 



if __name__ == '__main__':
    # 