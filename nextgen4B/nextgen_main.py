import nextgenfilter as ngf
import nextgenanalyze as nga
import ngsCSVAnalyze as csva

from Bio import SeqIO

import yaml
import cPickle as pickle
import sys, os

import pandas as pd

import logging
import time
from tqdm import tqdm

# Setup Logging
timestr = time.strftime("%Y%m%d-%H%M%S")
logging.basicConfig(filename='ngs_'+timestr+'.log',
    level=logging.DEBUG, format='%(asctime)s %(message)s')

def runAllExperiments(yfname, save_intermediates=True):
    # Load YAML file
    with open(yfname) as expt_f:
        expt_yaml = yaml.load(expt_f) # Should probably make this a class at some point...
    logging.info('Loaded YAML experiment file '+yfname)
    
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
        logging.info('Starting filtering for run %s' % run)
        aln_seqs = ngf.filterSample(runs[run]['f_read_name'], runs[run]['pe_read_name'],
                                    bcs, templates, 
                                    runs[run]['filter_seqs']['forward'],
                                    runs[run]['filter_seqs']['reverse'])
        if save_intermediates:
            for expt in aln_seqs.keys():
                with open('aln_seqs_%s_%s.fa' % (run, expt), 'w') as f:
                    SeqIO.write(aln_seqs[expt], f, 'fasta')
            # pickle.dump(aln_seqs, open('aln_seqs_'+run+'.pkl', 'w')) # Save aln_seqs
        logging.info('Finished filtering for run '+run)
        
        # For each experiment in each run, do analysis
        analyzed_data = {}
        for expt in expts:
            analyzed_data_fname = expt+'_'+run+'_misinc_data.csv'
            if not os.path.isfile(analyzed_data_fname):
                # Get positional misincorporations
                template = templates[expt]
                analyzed_data[expt] = nga.doAnalysis(aln_seqs[expt], template)
                
                # Save dataframe
                of = open(analyzed_data_fname, 'w')
                analyzed_data[expt].to_csv(of)
                of.close()
            else:
                infile = open(analyzed_data_fname)
                analyzed_data[expt] = pd.read_csv(infile)
                infile.close()
                logging.info('Found, loaded, pre-existing analysis file: '+ analyzed_data_fname)
            
            # Do analysis
            # ... not yet
            logging.info('Finished analysis for experiment '+expt)
            
    
if __name__ == '__main__':
    if len(sys.argv) > 1:
        yaml_name = sys.argv[1]
    else:
        yaml_name = 'samples.yaml'
    
    runAllExperiments(yaml_name, save_intermediates=True)