# IsoFinder
IsoFinder is a program that predicts which mRNA transcripts are transcript variants from the same gene within a transcriptome, and clusters them. For each cluster, it also predicts a gene model. This program is designed for identifying splice patterns of genes and gene models when an assembled genome is unavailable.


IsoFinder has been designed to work in conjunction with an unpublished program named **Virtual Genome Walker (VGW)**, which locally assemble genome reads based on mRNA transcripts. Please contact me to direct you to the developer of the VGW program.


## Dependencies 
Please install the following:
- BioPython
- mafft
- BLAST+.


## Install
Download source from github at:

	git clone https://github.com/BiGzGiB/IsoFinder.git


## Usage

```
$python all_at_once.py -help
usage: all_at_once.py [-h] [-p [FLOAT]] [-i [INT]] [-s] [-b]
                      filename.FASTA filename.FASTA

positional arguments:
  filename.FASTA  Genome Walked sequence file
  filename.FASTA  Original sequence file before running Genome Walker

optional arguments:
  -h, --help      show this help message and exit
  -p [FLOAT]      Minimum percentage identity for the grouping to occur based
                  on the BLAST result (default: 95.0)
  -i [INT]        Minimum number of intronic sequence match for transcripts to
                  be grouped together if transcripts only have one matching
                  exon (default: 30)
  -s, --separate  Creates a new folder that contains information about
                  individual gene models (default: False)
  -b, --blast     Does not delete the BLAST output (default: False)
```

First positional argument is the output file of VGW named “Final_sequences.fasta”.

Second positional argument is the mRNA sequence/transcriptome file that was inputted into VGW to process and create the “Final_sequences.fasta” file.


Example:
'''
\# After moving to the play_data1 directory
python isofinder.py 1iter_Final_sequences.fasta transcriptome1.txt
'''

