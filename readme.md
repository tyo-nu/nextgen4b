# nextgen4b
### A library to analyze positional error rates in nextgen sequencing data

[![License](https://img.shields.io/badge/license-MIT-blue.svg?style=flat)](https://github.com/tcyb/nexgen4b/blob/master/license.txt)

A set of scripts to manipulate and process .fastq files from next-generation sequencing experiments for "rare-base" assays as conducted in the Tyo lab at Northwestern University. In general, these scripts will filter your data for bad sequences, align them to a given template, retrieve the bases at positions of interest, and analyze them in a number of ways.

If you have short, synthetic sequences and are interested in error rates at particular points, this might be the library for you.

## Installation

First, clone the repo:

    git clone https://github.com/tcyb/nextgen4b

Then install the required packages

    pip install -r requirements.txt

You can then install `nextgen4b` using `setup.py`:

    python setup.py install

If you plan on modifying or contributing code, you can install in-place via symlink:

    python setup.py develop

Lastly, you will need to install EMBOSS and a modified biopython

### Installing EMBOSS

This code relies on EMBOSS's optimized `needle` routine in order to perform sequence alignment. You can find instructions on installation [here](http://emboss.sourceforge.net/download/). You will need to add the EMBOSS `bin` directory to your path for Biopython to be able to access `needle`.

### Adding alignment metadata to Bio.AlignIO

For this code to run correctly, you need a modified AlignIO module that reads the alignment score into the Alignment object's 'annotations' field. This is addressed here: biopython/biopython#692

In order to install and modify biopython, you can do the following:

    git clone https://github.com/tcyb/biopython
    cd biopython
    git checkout d980d67d34329312174728427273f8ca063ca4aa

Then follow the instructions to install biopython.

	python setup.py build
	python setup.py test
	python setup.py install

Alternatively, you can directly modify the AlignIO code in your existing installation of biopython. Hopefully this will be resolved in a new release of biopython.

## Usage

### Data Assumptions

This code assumes paired-end reads in either FASTQ or gzipped FASTQ format.

### Get all your .fastq files into one folder

Typically, we start with sequencing data from the MiSeq at UIC's sequencing core. They give us back *.fastq.gz files in individual folders. Your mileage may vary, but most of the batch-processing routines here assume fastq.gz files all in one
common folder. How you get there doesn't make much difference, but for us, it can be a pain in the butt to do manually.

In case your data setup is the same as ours and you're working on Windows, you can get your data into the proper structure using the following command in the command line:

    for /r %i in (*.fastq.gz) do @move "%i" .

### Generate an experiment-describing YAML file

You will need a `*.yaml` file that describes the experiments and samples (or **runs**) in your dataset. This includes information on where to locate the sequence files, what sequence information to expect, and any other experimental notes. The file should look something like this:

    ngsruns:
        run1:
            experiments: [exp1]
            f_read_name: TC1_S1_L001_R1_001.fastq.gz
            filter_seqs:
                forward: &id002 [CATTGTCCCTAT, AGACCAAGTCTCTGCTACCGTA]
                reverse: &id003 ['']
            pe_read_name: TC1_S1_L001_R2_001.fastq.gz

        ...

        run20:
            experiments: [exp20]
            f_read_name: TC10_S10_L001_R1_001.fastq.gz
            filter_seqs:
                forward: *id002
                reverse: *id003
            pe_read_name: TC10_S10_L001_R2_001.fastq.gz

    experiments:
        exp1:
            barcode: ''
            exp_data: &id001 {}
            name: exp1
            template_seq: GGGCTAGTCGTCTGTATAGGTCTTGCTTCTATCTTTGGCTTCTGTATTTGTCGTCTTGCTTATTGTTCTTGTTCTTATGTTCTGTTCTGGTATTTCGGTT

        ...

        exp20:
            barcode: ''
            exp_data: *id001
            name: exp20
            template_seq: GGGCTAGTCGTCTGTATAGGTCTTGCTTCTATCTTTGGCTTCTGTATTTGTCGTCTTGCTTATTGTTCTTGTTCTTATGTTCTGTTCTGGTATTTCGGTT

#### Some documentation on nextgen4b YAML files

* `ngsruns` - Information about sequence files from either individual sequencing runs or data demultiplexed on the sequencer.
  * `runN` - Internal label for the run. Can be anything, must be unique within `ngsruns`
    * `experiments` - A list of what experiments are contained in the given sequencing run. References the keys in the `experiments` portion of the YAML file.
    * `f_read_name` - Where to locate the forward read data for the run. Should point to either a `fasta` or `fasta.gz` file.
    * `filter_seqs\forward` - A list of sequences to look for in forward reads. Absence of these sequences in a given read will cause the read to be filtered. The first sequence in this list will define where the copied DNA starts. The second sequence will define the end of copied DNA. 
    * `filter_seqs\reverse` - A list of sequences to look for in paired end reads. Syntax is identical to `filter_seqs\forward`
    * `pe_read_name` - Where to locate the paired-end read data for the run. Should point to either a `fasta` or `fasta.gz` file.
* `experiments` - Information about experiments represented in sequencing runs. Multiple experiments can exist in one run, and a given experiment can have multiple instances if it occurs in multiple runs (these are not combined).
  * `expN` - Internal label for the experiment. Can be anything, must be unique within `experiments`. These should be entries in various `ngsruns\runN\experiments` lists.
  * `barcode` - A sequence that would appear 5' of the extension primer for demultiplexing during early iterations of experiments.
  * `exp_data` - A dictionary to indicate experimental conditions for a given experiment. User-defined and not of consequence here.
  * `name` - The name of the given experiment. User-defined and can be any string.
  * `template_seq` - The expected sequence data. This will be what reads are aligned against.

This file should adhere to YAML specs. Tricks and pain points have included:
  * Indentation **must** use spaces. Hard tabs will cause an error when parsing.
  * Repeated values can be aliased using the `&id`/`*id` syntax
  * Fields that should not be used for filtering should be indicated using `''`

### Performing analysis

Usage has changed a bit, full documentation coming.

## Author

* Ted Cybulski

## License

MIT License &copy; 2016 Thaddeus Cybulski
