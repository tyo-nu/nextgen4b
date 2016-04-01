from Bio import SeqIO
import yaml
import sys
import os

def get_positions(f_name, sites, keep_dashes=True):
    words = []
    with open(f_name) as f:
        seq_iter = SeqIO.parse(f, 'fasta')
        for s in seq_iter:
            selected_letters = ''.join([str(s.seq[i]) for i in sites])
            if '-' not in selected_letters or keep_dashes:
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