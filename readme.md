A set of scripts to manipulate and process .fastq files from next-generation sequencing experiments. The workflow is generally as follows.

# Get all your .fastq files into one folder
We typically use a MiSeq, which generates 
Getting files from BaseSpace folder format to all sequences in one directory:

	for /r %i in (*.fastq.gz) do @move "%i" .

# Generate an experiment-describing YAML file

Add the YAML file into the same directory.
Will document this someday. In short: Generally consists of **experiment** and **run** entries. **run** entries should correspond to individual .fastq (or .fastq.gz) files, and should contain information about forward and paired-end read file locations, experiments contained, and adapter sequences used. **experiment** entries contain information about experimental conditions, barcode sequences used if you're demuxing on those, template information, and a unique ID (name).

Then you can run the code.

# Performing analysis

Nextgen_main.py usage:

	python nextgen_main.py [samples.yaml]

Note: For this code to run correctly, you need a modified AlignIO module that reads the alignment score into the Alignment object's 'annotations' field.

This will generate a whole bunch of .csv's, one for each experiment in each run. Each .csv should have information about the observed reads at each position on the template supplied in the **experiment** entry in the YAML file.

You can then run

	python ngsCSVanalyze.py [n]

to get the summarized output for entry n for all of the .csv's in the directory.

There is also code in ngsCSVanalyze.py to get full summaries for each experiment, but that's been disabled for now.
