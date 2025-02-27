# PipelineProject
COMP 381 Python Pipeline

In order to run this pipeline, the following dependencies are required: biopython, NCBI's command line tools,
Bowtie2, spades, Blast+, and unzip. Run the following command to install these if you haven't already

pip install biopython install ncbi-datasets-cli

For Ubuntu, run the following command:

sudo apt-get install bowtie2 spades ncbi-blast+ unzip

on macOS, run the following command:

brew install bowtie2 spades blast unzip

To run this code, the following command line prompt format is used:

python Pipeline_wrapper.py --input (directory you downloaded files to) -f1 fastq1 -f2 fastq2 -f3 fastq3 -f4 fastq4

To run this with the given sample data, the command would use the following fastq files:

python Pipeline_wrapper.py --input (directory you downloaded files to) -f1 0030_1_sample.fastq -f2 0030_2_sample.fastq -f3 0033_1_sample.fastq -f4 0033_2_sample.fastq

To obtain the data included in this repo, I used the following commands to retrieve the given dataset from
NCBI and dump the data into paired-end fastq files

wget https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR5660030/SRR5660030

wget https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR5660033/SRR5660033

fasterq-dump SRR5660030

fasterq-dump SRR5660033

I retrieved the links to the SRR files from the following studies:

https://www.ncbi.nlm.nih.gov/sra/SRX2896360

https://www.ncbi.nlm.nih.gov/sra/SRX2896363

I went to the SRR runs for each of these and downloaded the full SRA archive data from the data access tab of NCBI.This gave me the entire dataset, which I ran the code on and stored the results in PipelineProject_full_set.log

To obtain the sample datasets provided, I ran the following command for each of the 4 paired-end files

head -n 40000 SRR5660030_1.fastq > 0030_1_sample.fastq

... and so on, for each of the files