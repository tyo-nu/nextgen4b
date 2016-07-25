import os, sys
import yaml
from Bio import SeqIO

def getRun(fname):
    # Can modify this for whatever sequencing standard you use
    # This is for MiSeq maybe? It's from the UIC sequencing core

    tokens = fname.split('_')
    sample_number = int(tokens[1][1:])
    strand_number = int(tokens[3][1:])

    return (sample_number, strand_number)

def formatRun(f1, f2, fseqs, peseqs, expnames):
    runDict = {}
    runDict['f_read_name'] = f1
    runDict['pe_read_name'] = f2
    runDict['filter_seqs'] = {'forward': fseqs, 'reverse': peseqs}
    runDict['experiments'] = expnames
    return runDict

def formatExpt(name, template, bc='', exp_dict={}):
    exptDict = {}
    exptDict['name'] = name
    exptDict['template_seq'] = str(template.seq)
    exptDict['barcode'] = bc
    exptDict['exp_data'] = exp_dict
    return exptDict

def loadTemplates(fname, ftype='fasta'):
    return list(SeqIO.parse(open(fname), ftype))

def generateYAMLDict(fSeqs, templateFname, peSeqs=[''], oneExpPerSamp=True, oneTemplate=True, barcoded=False):
    # get only fastq/fastq.gz filenames in current directory
    fileFnames = [f for f in os.listdir('.') if os.path.isfile(f) and f.endswith(('fastq', 'gz'))]

    # Put forward, paired-end file names in a run-indexed dictionary
    runData = {}
    for f in fileFnames:
        (sample, strand) = getRun(f)
        if sample not in runData.keys():
            runData[sample] = {strand: f}
        else:
            runData[sample][strand] = f

    templates = loadTemplates(templateFname)

    # Get formatted dictionary for each run/experiments
    runs = {}
    exps = {}
    for sample in runData.keys():
        if oneExpPerSamp:
            # set up templates
            if oneTemplate:
                template = templates[0]
            else:
                template = templates[sample]

            # set up run, experiment dictionaries
            runName = 'run'+str(sample)
            expName = 'exp'+str(sample)
            runs[runName] = formatRun(runData[sample][1], runData[sample][2], fSeqs, peSeqs, [expName])
            exps[expName] = formatExpt(expName, template)

    yamlDict = {}
    yamlDict['ngsruns'] = runs
    yamlDict['experiments'] = exps

    return yamlDict


if __name__ == '__main__':
    # Haha... we should have documented this when we wrote it...
    # input 1 is the output filename, input 2 is a fasta file with your template sequence

    # Hardcoded for Alex's 4-base stuff. Should probably change this.
    # Now its the common sequences... should work for both of us?
    fSeqs = ['GGGCTAGTCGTCTGTATAGG','AGACCAAGTCTCTGCTACCGTA']
    peSeqs = [''] # why?

    yDict = generateYAMLDict(fSeqs, sys.argv[2], peSeqs=peSeqs)

    # Write YAML file
    with open(sys.argv[1], 'w') as of:
        of.write(yaml.dump(yDict))
    