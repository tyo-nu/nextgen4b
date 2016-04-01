# Ted Cybulski
# 4/1/16

"""
Implemented to analyze Alex's 6mer dataset. It takes in a list of indices 
to define the sites that make up a motif, and anotherlist of indices that 
define the sites we want to count at. Finds unique motifs, then counts the 
incidence of nucleotides at each of the count sites.
"""

from Bio import SeqIO
import pandas as pd

import sys, os, argparse
    
def nmer_extract_from_fa(f_name, mot_idxs, ct_idxs,
                         bad_chars=['A','-'], letter_order=['C','G','T','A']):
    """
    Open a .fasta file, extract a given set of bases as a motif 
    and a given set of bases as counted nts.
    
    Arguments:
    f_name          -- The name of the .fasta file to process
    mot_idxs        -- The 0-indexed indices of the letters that make up a 
                        motif, in order
    ct_idxs         -- The 0-indexed indices of the letters to be counted 
                        (for errors, etc.)
    
    Keywords:
    bad_chars       -- Discard motifs that include these characters
    letter_order    -- The letters we're looking for (in order)
    """
    
    seqs = SeqIO.parse(f_name, 'fasta')
    mot = []
    ct_nts = []
    max_idx = max(max(mot_idxs), max(ct_idxs))
     
    for s in seqs:
        # construct motif for each sequence
        if len(s) > max_idx:
            m = ''.join([s[idx] for idx in mot_idxs])
            
            # if it's acceptable, record the motif and count sites
            if all([c not in m for c in bad_chars]):
                mot.append(m)
                ct_nts.append([s[idx] for idx in ct_idxs])
    
    unique_mots = list(set(mot))
    
    # Get totals for each of the motifs
    mot_cts = [mot.count(m) for m in unique_mots]
    
    # Build df
    df = pd.DataFrame({'motif': unique_mots, 'total': mot_cts})
    for c_idx, i in zip(ct_idxs, range(len(ct_idxs))):
        site_letters = [l[i] for l in ct_nts] # Letters at each site
        for letter in letter_order:
            # Get desired letter count for each motif
            site_mot_ct = [[site_letters[i] for i in range(len(mot)) \
                            if mot[i] == m].count(letter) for m in unique_mots]
            df['%s_%i_counts' % (letter, c_idx)] = site_mot_ct
    
    return df


def get_all_fnames(directory='.', suffix='.fa'):
    """
    Find all files in the given directory with a given suffix.
    """
    fnames = [f for f in os.listdir(directory)
                if os.path.isfile(f) and f.endswith(suffix)]
    return fnames


if __name__ == '__main__':
    # Parse stuff
    parser = argparse.ArgumentParser(
        description='Get motifs and errors from all .fa files in directory.')
    parser.add_argument('-M', '--motifsites', nargs='+', metavar='M',
                        help='Sites to look for motif bases', required=True)
    parser.add_argument('-C', '--countsites', nargs='+', metavar='C',
                        help='Sites to count bases at', required=True)       
    args = parser.parse_args(sys.argv[1:])
    
    # Do work
    found_files = get_all_fnames(suffix='.fa')
    
    for f in found_files:
        df = nmer_extract_from_fa(f, args.motifsites, args.countsites)
        df.to_csv(''.join(f.split('.')[:-1])+'_motifs.csv')