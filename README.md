# Ring nuclease search and characterisation pipeline
This repo builds on top of another repo ([https://github.com/vihoikka/hoikkala_etal_typeIII_effectors](url)).
This pipeline is modified towards finding ring nucleases (RNs) in type III CRISPR-Cas loci.
It also looks at phage genomes for RNs. The instructions below assume basic understanding of Snakemake, python and shell/bash environment.

## Requirements
* Unix based machine (Linux, MacOS)
* Conda
* Snakemake
* Python 3

*Other dependencies are installed by the pipeline when running with the --use-conda flag.*

### Databases
#### Databases included in this repository:
* A DIAMOND database of ring nucleases
* Cas10 HMM database
* Effector HMM database
* Ring nuclease HMM database
* A temperature database for prokaryotes ([TEMPURA](http://togodb.org/db/tempura))

#### Databases that require custom installation
* NCBI refseq genomes of bacteria and archaea (.fna, .faa and .gff files). I recommend downloading the database using ncbi-datasets (taxon 2 for bacteria; taxon 2157 for archaea):
```shell
datasets download genome taxon 2 --filename bacteria_genomes.zip --include gff3,genome,protein --dehydrated --annotated --assembly-level complete --assembly-source RefSeq
unzip bacteria_genomes.zip
datasets rehydrate --directory .
```
    * Repeat the download for archaea (taxon 2157). Next, *copy* the archaeal genome folders to the bacterial folder so that they are all one big mess (*leave original archaea folders as they are*). The pipeline will use the contents of the original archaea folder and the contents of bacteria (+archaea) folder to create a file that lists the domain of each genome entry. Not too classy, but the only way I could think of how to combine all prokaryotes in the same folder without losing track of domains.
* Millard phage database. Download the latest database from Millardlab.org. If using e.g. Sept. 2024 ([link](https://millardlab.org/bacteriophage-genomics/phage-genomes-sept-2024/)), you need the files 1Sep2024_genomes.fa, 1Sep2024_vConTACT2_proteins.faa and 1Sep2024_data.tsv.
* TMHMM academic license (get from [https://services.healthtech.dtu.dk/cgi-bin/sw_request?software=tmhmm&version=2.0c&packageversion=2.0c&platform=Linux]). Once registered, you will receive a model file called TMHMM2.0.model. Place the file in the folder data/TMHMM.

### Other modifications
You will need to modify multiple paths in the snakefile to make the pipeline work on your system. You will find these lines clearly marked in the beginning of snakefile rn.smk

**Example run**
```snakemake --snakefile rn.smk --use-conda --cores 40 --config protein_clustering="False" getGenomesBy="local" genome_mode="all" cas10_anchor="True" --rerun-triggers mtime```
Modify the core count to match your system.

> **Note on CRISPR subtype corrections**
> After the first run, you may need to perform manual corrections of CRISPR-Cas subtypes.
> The rule mastercombiner contains a python dictionary for correcting the subtypes.
> Note that the locus id given to each locus is arbitrary and stochastic – corrections therefore need to be applied every time the pipeline is run and the existing corrections in the existing codebase are very likely **not** applicable to your run. To make the manual corrections, I recommend producing a phylogenetic tree of Cas10s using the accompanying R script and manually checking if some loci seem out of place. Modify the dictionary in rule mastercombiner to include locus specific corrections and rerun starting from mastercombiner (rerunning means deleting the mastercombiner output file mastertable_v2.tsv and rerunning the pipeline. You may need to remove other downstream output files, such as the "done" file to trigger this rule again).

### The flow of the pipeline
1. Obtains bacterial and archaeal genomes from local files (there's also a remote option if needed). Also creates annotation file to hold information on which acc number corresponds to bacteria and which to archaea
2. Filters all genomes by presence of Cas10 using an HMM profile
3. Extracts Cas10 protein sequences and creates symlinks to the corresponding host genomes
4. Extracts host taxonomy infromation from the genome files (rule getTaxInfo)
5. Runs CRISPRCasTyper on Cas10-assocation genomes
6. Proceeds only with type III loci
7. Runs an in-house script cATyper on all loci, looking for pre-established effector proteins using HMM
8. Create HMM profiles based on pre-flagged HD-domain or GGDD-domain containing Cas10 (rules rule Cas10_HD_hmm_maker and rule Cas10_GGDD_hmm_maker). Note that different profiles are used for GGDD and GGDE seeded cyclase domains.
9. Create Cas10 tree using muscle -based alignment and fasttree
10. Similar to step 7, find and annotate ring nucleases based on HMM profiles
11. Also find ring nucleases on phage genomes (using Millard phage db).
12. An optional step is to characterise prophages using the rule dbscan_swa. To enable, uncomment the dbscan_swa line from the input of rule final.
12. Combine all results together into a mastertable
13. Generate website and other output files for downstream analysis (folder R contains data to be analysed using accompanying R scripts)