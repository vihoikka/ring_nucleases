'''
Finds all pairwise combinations of alignments and extracts only the sequences that are present in both alignments.
This version is locus-specific. The sample-specific version is otherwise compatible, except for
the fact that 16S headers do not contain locus-specific information. This is why when comparing 16S's,
we split the sample headers to contain only the sample part. This will create some pseudoreplication, 
so keep this in mind.
'''


import pandas as pd
import os
import argparse
from tqdm import tqdm
from concurrent.futures import ProcessPoolExecutor
import ast
import Bio

parser = argparse.ArgumentParser()
parser.add_argument("--cas10", required = True)
parser.add_argument("--cas7", required = True)
parser.add_argument("--cas5", required = True)
parser.add_argument("--cora", required = True)
parser.add_argument("--rna16s", required = True)
parser.add_argument("--outfolder", required = True)

args = parser.parse_args()
cas10 = args.cas10
cas7 = args.cas7
cas5 = args.cas5
cora = args.cora
rna16s= args.rna16s
outfolder = args.outfolder


import itertools

# Define a function to read a FASTA-format alignment file
def read_alignment(filename):
    sequences = {}
    with open(filename) as f:
        for line in f:
            if line.startswith(">"):
                name = line[1:].strip()
                sequences[name] = ""
            else:
                sequences[name] += line.strip()
    return sequences

# Define the list of alignment filenames
filenames = [cas10, cas7, cas5, cora, rna16s]

# Read each alignment file and store the sequences in a dictionary
alignments = {}
for filename in filenames:
    print("Reading " + str(filename))
    sequences = read_alignment(filename)
    alignments[filename] = sequences

# Find the pairwise combinations of alignments
pairs = list(itertools.combinations(filenames, 2))
print("Combinations: " + str(pairs))

# Loop over each pair of alignments
for pair in pairs:
    print("Comparing " + str(pair[0]) + " and " + str(pair[1]) + "...")
    # Get the sequences from each alignment
    seq1 = alignments[pair[0]]
    seq2 = alignments[pair[1]]

    new_dict_seq1 = {}
    new_dict_seq2 = {}

    #remove possible spaces from sequence header
    for key in seq1.keys():
        new_key = key.split()[0]
        #also, if the other set of sequences is 16S, split these headers to only contain sample part
        if pair[1] == rna16s:
            print("Since " + str(pair[1]) + " is 16S, splitting the header to only contain sample part.")
            #split by "_" and get the first two parts. Retain the "_" between them
            new_key = "_".join(new_key.split("_")[:2])
        new_dict_seq1[new_key] = seq1[key]


    for key in seq2.keys():
        new_key = key.split()[0]
        #also, if the other set of sequences is 16S, split these headers to only contain sample part.
        if pair[0] == rna16s:
            print("Since " + str(pair[0]) + " is 16S, splitting the header to only contain sample part.")
            #split by "_" and get the first two parts. Retain the "_" between them
            new_key = "_".join(new_key.split("_")[:2])
        new_dict_seq2[new_key] = seq2[key]

    seq1 = new_dict_seq1
    seq2 = new_dict_seq2

    # Get the set intersection of the entry names
    print("Finding intersection of entry names...")
    names1 = set(seq1.keys())
    names2 = set(seq2.keys())
    shared_names = names1.intersection(names2)

    print(shared_names)

    # Extract only the filename (without the directory path and the "alignment" part)
    pair_filenames = [os.path.basename(f)[:-4].split("_")[0] for f in pair]

    print("Pair filenames: " + str(pair_filenames))

    base_folder_name = f"{outfolder}/{pair_filenames[0]}_{pair_filenames[1]}"
    os.makedirs(base_folder_name)

    # Write new alignment files with only the shared entries
    for name in shared_names:
        print("Writing " + name + "...")
        with open(f"{base_folder_name}/{pair_filenames[0]}_{pair_filenames[1]}_{pair_filenames[0]}.afa", "a") as f1, \
                open(f"{base_folder_name}/{pair_filenames[0]}_{pair_filenames[1]}_{pair_filenames[1]}.afa", "a") as f2:
            f1.write(">" + name + "\n" + seq1[name] + "\n")
            f2.write(">" + name + "\n" + seq2[name] + "\n")

    with open(f"{base_folder_name}/{pair_filenames[0]}_{pair_filenames[1]}_common.csv", "a") as commonNamesFile:
        commonNamesFile.write("Genome\t" + pair_filenames[0] + "_" + pair_filenames[1] + "\n")
        for name in shared_names:
            commonNamesFile.write(name + ",True\n")

    with open(f"{outfolder}/{pair_filenames[0]}_{pair_filenames[1]}.done", "w") as donefile:
        donefile.write("done")