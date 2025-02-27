#modules to import
import os
from Bio import Seq, SeqIO
from Bio.Seq import Seq
import sys
import argparse

#argument parser so different inputs can be used
def check_arg(args=None):
    parser = argparse.ArgumentParser(
    description="HW1 Problem 1 (shows on help -h)")
    parser.add_argument("-i", "--input",
    help="input directory",
    required=True)
    parser.add_argument("-f1", "--fastq1", help="First FASTQ file", required=True)
    parser.add_argument("-f2", "--fastq2", help="Second FASTQ file", required=True)
    parser.add_argument("-f3", "--fastq3", help="Third FASTQ file", required=True)
    parser.add_argument("-f4", "--fastq4", help="Fourth FASTQ file", required=True)
    return parser.parse_args(args)

arguments = check_arg(sys.argv[1:])
infile = arguments.input
#set the fastq inputs to the given files
fastq1 = arguments.fastq1
fastq2 = arguments.fastq2
fastq3 = arguments.fastq3
fastq4 = arguments.fastq4

#create the ouput directory (i put exist_ok in for when i was running and re-running my code)
dir_create = 'PipelineProject_Yusef_Golzar'
os.makedirs(dir_create, exist_ok=True)#change into the created directory
os.chdir(dir_create)
#set the accession number for the bowtie. i decided to include this so that this script could
#also be used for other genomes
accession = 'GCF_000845245.1'

#create the bowtie command line prompt and execute it with os.system
bowtie_dataset = f'datasets download genome accession {accession} --include gff3,rna,cds,protein,genome,seq-report'
os.system(bowtie_dataset)
#unzip the downloaded dataset
unzip_bowtie = f'unzip ncbi_dataset.zip'
os.system(unzip_bowtie)
#build the genome from the imported dataset
bowtie_build = f'bowtie2-build {infile}/PipelineProject_Yusef_Golzar/ncbi_dataset/data/{accession}/{accession}_ViralProj14559_genomic.fna HCMV'
os.system(bowtie_build)
#map the input fastq files to the imported genomes
bowtie_1_map = f'bowtie2 --quiet -x HCMV -1 {infile}/{fastq1} -2 {infile}/{fastq2} -S HCMV2dpimap_sample.sam -k 1 --no-unal'
bowtie_2_map = f'bowtie2 --quiet -x HCMV -1 {infile}/{fastq3} -2 {infile}/{fastq4} -S HCMV6dpimap_sample.sam -k 1 --no-unal'
os.system(bowtie_1_map)
os.system(bowtie_2_map) 

#print the number of reads before and after bowtie filtering to the log file
don_1_sample = f'echo "Donor 1 (2dpi) had $(($(wc -l < {infile}/{fastq1}) / 4)) read pairs before bowtie2 filtering and $(($(grep -v \'^@\' HCMV2dpimap_sample.sam | wc -l) / 2)) read pairs after" > PipelineProject.log'

don_2_sample = f'echo "Donor 1(6dpi) had $(($(wc -l < {infile}/{fastq3}) / 4)) read pairs before bowtie2 filtering and $(($(grep -v \'^@\' HCMV6dpimap_sample.sam | wc -l) / 2)) read pairs after" >> PipelineProject.log'

os.system(don_1_sample)
os.system(don_2_sample)

#create the spades command with all 4 paired-end files
spades_command = f'spades.py --only-assembler -k 99 -t 2 -1 {infile}/{fastq1} -2 {infile}/{fastq2} -1 {infile}/{fastq3} -2 {infile}/{fastq4} -o spades_assembly'
os.system(spades_command)
#write the spades command to the log
write_spades_command = f'echo {spades_command} >> PipelineProject.log'
os.system(write_spades_command)

#create a variable to find the contigs from spades
spades_contigs = 'spades_assembly/contigs.fasta'
#save the sequences from the contigs file 
records = list(SeqIO.parse(spades_contigs, 'fasta'))
sequences = [str(record.seq) for record in records]

#create variables to store the number of contigs with length greater than the given threshold
#and to store the total length of these contigs
seqs_over_thresh = 0
total_seq_length = 0
#loop through all the sequnces
for sequence in sequences:
    #store the length of the current sequence
    seq_length = len(sequence)
    #if the seq length is over the threshold, add 1 to the counter and add the length of the sequence to the total length
    if seq_length > 1000:
        seqs_over_thresh += 1
        total_seq_length += seq_length
#write the total number of sequences over the threshold and the total seq length to the log file
with open('PipelineProject.log', 'a') as outfile:
    outfile.write(f'There are {seqs_over_thresh} contigs > 1000 bp in the assembly \n')
    outfile.write(f'There are {total_seq_length} bp in the assembly \n')

#create an empty variable for the longest contig and set a vvariable for its length to 0
longest_contig = None
longest_contig_length = 0
#for every contig
for sequence in sequences:
    #save the length of the contig
    contig_length = len(sequence)
    #if the length of the contig is greater than the currently stored longest one,
    if contig_length > longest_contig_length:
        #the new longest contig is the current one,
        longest_contig = sequence
        #and the highest length is the length of that contig
        longest_contig_length = contig_length

#create a fasta file to store the contig so that it can be properly input into the blast query
with open('longest_contig.fasta', 'w') as outfile:
    outfile.write(longest_contig)
#set the query file to the file I just made
query_seqfile = 'longest_contig.fasta'
#download all of the Betaherpesvirinae datasets from ncbi to be used to create a local database.
download_betasets = 'datasets download virus genome taxon Betaherpesvirinae --include genome --filename Betaherpesvirinae_dataset.zip'
os.system(download_betasets)
#unzip the imported data
unzip_betasets = 'unzip Betaherpesvirinae_dataset.zip -d Betaherpesvirinae_data'
os.system(unzip_betasets)
#create the path to the file so that it can be changed if necessary
genomic_file_path = 'Betaherpesvirinae_data/ncbi_dataset/data/genomic.fna'
#create the local database using the data I imported
make_local_db = f'makeblastdb -in {genomic_file_path} -out Betaherpesvirinae -title Betaherpesvirinae -dbtype nucl'
os.system(make_local_db)

#perform a blastn query using the local database and the longest contig file, then retrieve a tab-delimited output with the specified parameters
blast_query = 'blastn -query longest_contig.fasta -db Betaherpesvirinae -out blast_results.txt -outfmt "6 sacc pident length qstart qend sstart send bitscore evalue stitle"'
os.system(blast_query)

#create a list to store the hits from the blast query
top_hits = []
#open the blast results file
with open('blast_results.txt', 'r') as infile:
    for line in infile:
        #strip each of the lines of their \t characters so they can be read properly
        fields = line.strip().split('\t')
        #append the current data into the top_hits list
        top_hits.append(fields)
        #once the top ten hits are stored in the list, stop looping (these are the top 10 that I need)
        if len(top_hits) >= 10:
            break
#open the log file
with open('PipelineProject.log', 'a') as outfile:
    #write the header containing each of the parameters I needed
    outfile.write("sacc\tpident\tlength\tqstart\tqend\tsstart\tsend\tbitscore\tevalue\tstitle\n")
    #write the top ten hits underneath the header
    for hit in top_hits:
        outfile.write('\t'.join(hit) + '\n')