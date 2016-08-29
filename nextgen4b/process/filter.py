"""
nextgen4b.process.filter
TC, 8/11/16

A collection of functions that read, filter, and output sequence data from
next-generation sequencing experiments.
"""
import gzip
import logging
import os
import re
import sys
import time
import uuid

import numpy as np
import yaml
from Bio import AlignIO, SeqIO
from Bio.Emboss.Applications import NeedleCommandline
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from tqdm import tqdm

__all__ = ['filter_sample', 'run_all_experiments']

#####################
# File Management
#####################

def load_ngs_file(fpath, ftype='fastq'):
    """
    Load a .fastq file to a SeqIO iterator, un-gzip if necessary.
    """
    if fpath.endswith('.gz'):
        seq_f = gzip.open(fpath, 'rb')
    elif fpath.endswith('.fastq'):
        seq_f = open(fpath, 'rb')
    else:
        raise ValueError('File does not end in .gz or .fastq; confirm file type.')
    f_iter = SeqIO.parse(seq_f, ftype)
    return f_iter

#####################
# Main Filter Code
#####################

def filter_sample(f_name, pe_name, bcs, templates, f_filt_seqs, r_filt_seqs):
    """
    Output filtered sequences as dictionary, indexed by barcode.
    Sequences will be aligned to the provided template.
    Parts of the template not represented will be '-'
    """
    logging.info('Started filtering routine for %s', f_name)

    # Compile regexes
    f_res = compile_res(f_filt_seqs)
    pe_res = compile_res(r_filt_seqs)

    # Load as generators, then filter
    logging.info('Loading Files')
    f_seqs = load_ngs_file(f_name)
    for regex in f_res:
        f_seqs = filter_seqs(f_seqs, regex)
    pe_seqs = load_ngs_file(pe_name)
    for regex in pe_res:
        pe_seqs = filter_seqs(pe_seqs, regex)

    # Barcode Filtering/Demux
    bc_seqs = barcodeDemux(f_seqs, bcs)

    # Sequence-based filtering
    for expt in bcs.keys():
        logging.info('Starting post-demux filtering for expt ID %s', expt)
        # Filter based on PE matches, only return the copied sequence
        # Assumes the first RE in f_res will terminate the copied sequence
        # copiedFuncGenerator's output should return all sequence before the adapter
        seqs = filter_pe_mismatch(bc_seqs[expt], pe_seqs,
                                  gen_copied_seq_function(f_res))

        # Quality filter
        if len(seqs) > 0:
            seqs = quality_filter(seqs) # Quality Filtering (needs to only have copied sequence)
        else:
            logging.info("""No sequences left, skipped quality score
                            filtering.\n*** Skipping remaining filtering for
                            expt ID %s***""", expt)
            bc_seqs[expt] = seqs
            continue

        seqs = strip_forward_barcodes(seqs, l_barcode=len(bcs[expt])) # Remove barcodes before align

        # Align filter
        if len(seqs) > 0:
            seqs = alignment_filter(seqs, templates[expt]) # Do alignment-based filtering
        else:
            logging.info("""No sequences left, skipped align filtering.\n***
                            Skipping remaining filtering for expt ID %s***""",
                         expt)
            bc_seqs[expt] = seqs
            continue

        # Length filtering
        if len(seqs) > 0:
            seqs = len_filter(seqs, l_barcode=len(bcs[expt])) # Length Filtering
        else:
            logging.info("""No sequences left, skipped length filtering.\n***
                            Skipping remaining filtering for expt ID %s***""",
                         expt)
            bc_seqs[expt] = seqs
            continue

        bc_seqs[expt] = seqs

    return bc_seqs

#####################
# F/R Regex Filtering
#####################

def filter_seqs(seqs, q_re):
    """
    Filter an iterator based on whether items match a regex object.
    """
	# sIO_iterator should be generated by Bio.SeqIO.parse
	# q_re should be a regex object generated by re.compile()
	# Outputs a list of Seq objects that have the desired sequence in them.
    logging.info('Started regex filter: %s', q_re.pattern)
    out_l = [s for s in seqs if q_re.search(str(s.seq))]
    logging.info('Finished regex filter. Kept %i sequences.', len(out_l))
    return out_l

def compile_res(seqs):
    """
    Compile regex for each string in a list, return list of regex objects.
    """
    # Takes a list of sequences you want to filter for
    # Outputs a list of regex objects that you can iterate over
    return [re.compile(s) for s in seqs]

#####################
# Barcode Filtering
#####################

# Need to do this before alignment.
def strip_forward_barcodes(f_seqs, l_barcode=5):
    return [s[l_barcode:] for s in f_seqs]

def barcodeDemux(seqs, bcs):
    """
    Takes lists of sequence objects, dict of barcodes (indexed by expt. ID)
    Demuxes based on the barcode the sequences start with
    Discards sequences that don't start with a barcode exact match
    Assumes forward read -> sequences start with a barcode
    """

    logging.info('Started barcode demuxing.')

    n_seqs = 0
    bc_filtered_data = {}
    for expt in bcs.keys():
        bc_filtered_data[expt] = [s for s in seqs if str(s.seq).startswith(bcs[expt])]
        n_seqs += len(bc_filtered_data[expt])

    logging.info('Finished barcode demuxing. Kept %i of %i sequences.', n_seqs, len(seqs))

    return bc_filtered_data

#####################
# Paired End Match Filtering
#####################

def get_coords(s):
    return ':'.join(s.description.split(' ')[0].split(':')[3:])

def get_sense(s):
    return s.description.split(' ')[1].split(':')[0]

def get_copied_seq(s, f_res):
    return s[f_res[0].search(str(s.seq)).end():list(f_res[1].finditer(str(s.seq)))[-1].start()]

def gen_copied_seq_function(f_res):
    return lambda s: get_copied_seq(s, f_res)

def filter_pe_mismatch(f_seqs, pe_seqs, copied_func):
    """
    Args:
        f_seqs - sequences from forward reads. Presumably filtered for the
                 required adatper(s).
        pe_seqs - the paired end sequences of f_seqs. Also presumably filtered
                  for the required adapter(s).
        copied_func - takes a sequence, should ouptut the DNA that we expect
                      to have been copied, i.e. that should be on the paired
                      end read.

    Outputs a list of forward sequences that pass two filters:
        * Have a coordinate match in the paired end reads
        * That coordinate match has the same sequence.

    Prunes the sequences down to what was actually copied (i.e. pre-adapter)
    """

    logging.info('Started Paired-End Filtering')

    # Some housekeeping stuff
    proc_ct = 0 # number of sequences processed
    co_ct = 0 # number of sequences with coordinate matches
    aln_ct = 0 # number of sequences that have paired end sequence matches
    matched_seq_list = []

    # Get coordinate list
    pe_coordL = [get_coords(s) for s in pe_seqs]
    for s in f_seqs:
        if pe_coordL.count(get_coords(s)): # Filter based on paired-end presence
            co_ct += 1
            copied = copied_func(s) # Get the part of the sequence that was actually copied
            if str(pe_seqs[0].reverse_complement().seq).find(str(copied.seq)): # Filter on PE match
                aln_ct += 1
                matched_seq_list.append(copied)

        proc_ct += 1
        if not (proc_ct % 5000):
            logging.info("Processed %i out of %i", proc_ct, len(f_seqs))

    logging.info("Finished Paired-End Filtering")
    logging.info("""Kept %i of %i forward sequences after coordinate
                    filtering""", co_ct, len(f_seqs))
    logging.info("""Kept %i of %i forward sequences after paired-end sequence
                    matching""", aln_ct, co_ct)

    return matched_seq_list

#####################
# Q-score Filtering
#####################

def quality_filter(seqs, q_cutoff=20):
    logging.info('Started Quality Score Filtering')
    out_l = [s for s in seqs
             if not any(s.letter_annotations['phred_quality']
                        < np.ones(len(s.letter_annotations['phred_quality']))*q_cutoff)]
    logging.info('Finished Quality Score Filtering. Kept %i of %i sequences.',
                 len(out_l), len(seqs))
    return out_l

#####################
# Length Filtering
#####################

def len_filter(seqs, l_cutoff=70, u_cutoff=200, l_barcode=0):
    """
    Return only sequence objects that have length between l_cutoff and
    u_cutoff
    """
    logging.info('Started Length Filtering')
    out_l = [s for s in seqs if (len(s.seq) >= (l_cutoff + l_barcode)) and
             (len(s.seq) <= (u_cutoff + l_barcode))]
    logging.info('Finished Length Filtering. Kept %i of %i sequences.',
                 len(out_l), len(seqs))
    return out_l

#####################
# Alignment Filtering
#####################

def alignment_filter(seqs, template, gapopen=10, gapextend=0.5, lo_cutoff=300,
                     hi_cutoff=1000, cleanup=True):
    logging.info('Started alignment-based filtering')
    start_n_seqs = len(seqs)

    # Save the template and sequences as temporary fasta files
    # Probably some hacking that can be done in the NeedleCommandline stuff
    template_f_name = 'temptemplate.fa'
    seqs_f_name = 'tempseq.fa'

    with open(seqs_f_name, 'w') as sh:
        SeqIO.write(seqs, sh, 'fastq')

    with open(template_f_name, 'w') as temp_seq_file:
        # Make temp sequence file for alignment
        temp_seq = SeqRecord(Seq(template), id='template', name='template')
        SeqIO.write(temp_seq, temp_seq_file, 'fasta')

    # Generate alignment command, run the alignment
    logging.info("""Began EMBOSS needle routine with settings:\ngapopen:
                    %i\ngapextend: %i\nlo_cutoff: %i\nhi_cutoff: %i""",
                 gapopen, gapextend, lo_cutoff, hi_cutoff)
    ofilen = 'temp_'+str(uuid.uuid4())+'.needle'
    needle_cline = NeedleCommandline(asequence=template_f_name,
                                     bsequence=seqs_f_name, gapopen=gapopen,
                                     gapextend=gapextend, outfile=ofilen)
    needle_cline()
    logging.info('Finished EMBOSS needle routine')

    aln_data = AlignIO.parse(open(ofilen), "emboss")
    new_seqs = cull_alignments(aln_data, lo_cutoff=lo_cutoff,
                               hi_cutoff=hi_cutoff)

    # Exit routine
    if cleanup:
        logging.info('Cleaning up temp files')
        os.remove(template_f_name)
        os.remove(seqs_f_name)
        os.remove(ofilen)
    logging.info("""Finished alignment-based filtering. Kept %i of %i
                    sequences.""", len(new_seqs), start_n_seqs)
    return new_seqs


def cull_alignments(aln_data, lo_cutoff=400, hi_cutoff=650):
    new_seqs = []

    for alignment in aln_data:
        if (alignment.annotations['score'] > lo_cutoff) \
            and (alignment.annotations['score'] < hi_cutoff):
            # Template should have no gaps, and should contain the whole
            # non-template sequence
            if not str(alignment[0].seq).count('-') > 0:
                new_seqs.append(alignment[1])
                new_seqs[-1].annotations['alnscore'] = alignment.annotations['score']
    return new_seqs

#####################
# Main Routines
#####################

def run_all_experiments(yf_name, save_intermediates=True):
    """
    Filters all sequences noted in the passed YAML file.
    """
    # Start Logging
    timestr = time.strftime("%Y%m%d-%H%M%S")
    logging.basicConfig(filename='ngs_'+timestr+'.log',
                        level=logging.DEBUG, format='%(asctime)s %(message)s')

    # Load YAML file
    with open(yf_name) as expt_f:
        expt_yaml = yaml.load(expt_f) # Should probably make this a class at some point...
    logging.info('Loaded YAML experiment file '+yf_name)

    runs = expt_yaml['ngsruns']
    logging.info('Found NGS Runs: '+', '.join(runs))

    for run in tqdm(runs.keys()):
        logging.info('Performing routine for NGS Run '+run)
        expts = runs[run]['experiments']
        logging.info('Found experiments '+', '.join(expts))

        # Get barcodes, templates for all experiments in the run
        bcs = {}
        templates = {}
        for expt in expts:
            bcs[expt] = expt_yaml['experiments'][expt]['barcode']
            templates[expt] = expt_yaml['experiments'][expt]['template_seq']

        # Do filtering
        logging.info('Starting filtering for run %s', run)
        aln_seqs = filter_sample(runs[run]['f_read_name'],
                                 runs[run]['pe_read_name'],
                                 bcs, templates,
                                 runs[run]['filter_seqs']['forward'],
                                 runs[run]['filter_seqs']['reverse'])
        if save_intermediates:
            for expt in aln_seqs.keys():
                with open('aln_seqs_%s_%s.fa' % (run, expt), 'w') as out_f:
                    SeqIO.write(aln_seqs[expt], out_f, 'fasta')
        logging.info('Finished filtering for run %s', run)

if __name__ == '__main__':
    if len(sys.argv) > 1:
        yaml_name = sys.argv[1]
    else:
        yaml_name = 'samples.yaml'

    run_all_experiments(yaml_name, save_intermediates=True)