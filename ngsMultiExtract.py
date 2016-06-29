from Bio import SeqIO
import yaml
import sys
import os

def replace_deletions(word, seq, idxs, del_letter='d'):
    """
    Replace any '-' in word with del_letter if the nucleotides next to it in
    seq are not '-'.


    """
    new_word = [c for c in word]

    for i, letter in enumerate(word):
        # assume we're not looking at the start or end of sequence
        idx = idxs[i]
        assert idx > 0 and idx < len(seq)

        if letter == '-':
            if seq[idx-1] != '-' and seq[idx+1] != '-':
                new_word[i] = del_letter

    return ''.join(new_word)

def get_positions(f_name, sites, keep_dashes=True, mark_deletions=True):
    """
    Reads in a fasta file of sequences (usually produced by nextgen_main.py)
    at location f_name, and pulls out the bases at the (0-start) indices in
    sites.

    Input:
        - f_name: str
        - sites: list (ints)
        - keep_dashes: bool
        - mark_deletions: bool
    Output:
        - words: list

    Options:
    keep_dashes: if this is false, get_positions will discard any words with a
        dash in them (generally denoting a deletion)
    mark_deletions: if this is true, deletions (dashes flanked by non-dashes on
        both sides) will be marked (with a 'd', but this should be
        customizable?)
    """

    words = []
    with open(f_name) as f:
        seq_iter = SeqIO.parse(f, 'fasta')
        for s in seq_iter:
            selected_letters = ''.join([str(s.seq[i]) for i in sites])
            if '-' not in selected_letters:
                words.append(selected_letters)
            elif keep_dashes:
                if mark_deletions:
                    words.append(replace_deletions(selected_letters, s, sites))
                else:
                    words.append(selected_letters)

    return words

if __name__ == '__main__':
    in_name = sys.argv[1]
    out_name = sys.argv[2]
    sites = [int(x) for x in sys.argv[3:]]

    words = get_positions(in_name, sites, keep_dashes=True)
    
    with open(out_name, 'w') as of:
        for word in words:
            of.write('%s\n' % word)