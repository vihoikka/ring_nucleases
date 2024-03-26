'''
This script is meant to run after R has been run and a list of interesting loci has been extracted.
This list is a newline separated list of locus names that are for any reason deemed interesting.
The current script does the following:
1. Extracts all "unknown proteins" that have the locus name in the header
2. Subjects each such protein to a SAVED/CARF search
3. Extracts tidied file, then reopens it to run a hmmscan against local pfam
4. Outputs everything in one neat file
'''
import argparse
import pandas as pd
import subprocess
import sys
import os
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Blast.Applications import NcbiblastpCommandline

#parse command line arguments for the following: input genome list, path to unknown proteins fasta file, path to output file, project
parser = argparse.ArgumentParser(description='This script takes a list of loci and extracts the unknown proteins for each locus, then subjects them to a blastp search against the NCBI nr database.')
parser.add_argument('-i', '--input', type=str, metavar='input', required=True,
                    help='The input file containing the list of loci to extract unknown proteins for.')
parser.add_argument('-u', '--unknown_proteins', type=str, metavar='unknown_proteins', required=True,
                    help='The fasta file containing all unknown proteins.')
parser.add_argument('-o', '--output_basename', type=str, metavar='output', required=True,
                    help='The output file to write the results to.')
#parser.add_argument('-p', '--project', type=str, metavar='project', required=True)

args = parser.parse_args()

CARFSAVED_hmm_db_path = "/media/volume/st_andrews/databases/carfsaved/02_hmm_profiles/carfsaved.hmm"
pfams_hmm_db = "/media/volume/st_andrews/databases/pfam/Pfam-A.hmm"
project_root = "/media/volume/st_andrews/new_effectors"

#output files
all_against_pfam = args.output_basename + ".all.pfam.rawresults.txt"
carfsaved_raw_results = args.output_basename + ".CARFSAVED.rawresults.txt"
carfsaved_info = args.output_basename + ".CARFSAVED.seq.txt"
carfsaved_fasta = args.output_basename + ".CARFSAVED.fasta"
carfsaved_against_pfam = args.output_basename + ".CARFSAVED.pfam.rawresults.txt"
carfsaved_against_pfam_info = args.output_basename + ".CARFSAVED.pfam.info.txt"

E_value_pfam = "1e-05"
E_value_CARFSAVED = "1e-05"

#read in the list of loci
loci = []
with open(args.input, 'r') as f:
    for line in f:
        loci.append(line.strip())

#read in the unknown proteins fasta file
unknown_proteins = SeqIO.to_dict(SeqIO.parse(args.unknown_proteins, 'fasta'))

#extract any proteins from the fasta file that have the locus name in the header
unknown_proteins_to_blast = []
for locus in loci:
    for protein in unknown_proteins:
        if locus in protein:
            unknown_proteins_to_blast.append(unknown_proteins[protein])

#write the extracted proteins to a fasta file
SeqIO.write(unknown_proteins_to_blast, "unknown_proteins_to_blast.fasta", 'fasta')

#run blastp against the NCBI nr database
#blastp_cline = NcbiblastpCommandline(query="unknown_proteins_to_blast.fasta", db="nr", evalue=0.001, outfmt=6, out="blastp_results.txt")
#stdout, stderr = blastp_cline()

#run all unknown proteins against the pfam database
subprocess.run(['hmmscan', '--tblout', all_against_pfam, '--cpu', '40','-E',E_value_pfam, pfams_hmm_db, 'unknown_proteins_to_blast.fasta'])


#run hmmscan subprocess against our SAVED/CARF database
subprocess.run(['hmmscan', '--tblout', carfsaved_raw_results, '--cpu', '40','-E',E_value_CARFSAVED, CARFSAVED_hmm_db_path, 'unknown_proteins_to_blast.fasta'])

#open the hmmscan result file (hmmscan_results.txt) and extract the locus name, the protein name, the e-value and the score
#write these to a new file
with open(carfsaved_raw_results, 'r') as f:
    with open(carfsaved_info, 'w') as g:
        for line in f:
            if not line.startswith("#"):
                line = line.strip().split()
                locus = line[0]
                protein = line[2]
                evalue = line[4]
                score = line[5]
                g.write(locus + "\t" + protein + "\t" + evalue + "\t" + score + "\n")

##open the hmmscan result file (hmmscan_results.txt) and extract the locus name, the protein name, the e-value and the score and create a pandas dataframe
hmmscan_results = pd.read_csv(carfsaved_info, sep='\t', header=None)
print(hmmscan_results)
#define headers
hmmscan_results.columns = ["Target_name", "Query", "E-value", "Score"]

#create new header called Seq and set default value to nothing
hmmscan_results["Seq"] = ""
print(hmmscan_results)
carfsavedproteins = []

#open the fasta file and based on the "query", extract the sequence and insert it into the dataframe in a column called "Seq"
for index, row in hmmscan_results.iterrows():
    query = row["Query"]
    seq = unknown_proteins[query].seq
    hmmscan_results.at[index, "Seq"] = seq
    carfsavedproteins.append(SeqRecord(seq, id=query, description=""))

#write the dataframe to a file
hmmscan_results.to_csv(carfsaved_info, sep='\t', index=False)

#write the carfsaved proteins to a fasta file
SeqIO.write(carfsavedproteins, carfsaved_fasta, 'fasta')

#run hmmscan subprocess against our local pfam database using the above fasta file
subprocess.run(['hmmscan', '--tblout', carfsaved_against_pfam, '--cpu', '40','-E', E_value_pfam, pfams_hmm_db, carfsaved_fasta])

#tidy results 
with open(carfsaved_against_pfam, 'r') as f:
    with open(carfsaved_against_pfam_info, 'w') as g:
        for line in f:
            if not line.startswith("#"):
                line = line.strip().split()
                locus = line[0]
                protein = line[2]
                evalue = line[4]
                score = line[5]
                g.write(locus + "\t" + protein + "\t" + evalue + "\t" + score + "\n")

#load the results into a dataframe
carfsaved_against_pfam_results = pd.read_csv(carfsaved_against_pfam, sep='\t', header=None)