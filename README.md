# crispr_type_iii_effector_finder
A snakemake pipeline that characterises type III CRISPR-Cas loci and aims to discover new cOA dependent effectors. **Work in progress**.
Requires Snakemake and Conda. All other dependencies are installed when running with the --use-conda option.

Example run
```snakemake --snakefile new_effectors.smk --use-conda --cores 40 --config protein_clustering="False" getGenomesBy="local" genome_mode="all" cas10_anchor="True" --rerun-triggers mtime```

## The flow
1. Obtains bacterial and archaeal genomes from local files (there's also a remote option if needed). Also creates annotation file to hold information on which acc number corresponds to bacteria and which to archaea
2. Filters all genomes by presence of Cas10 using an HMM profile
3. Extracts Cas10 protein sequences and creates symlinks to the corresponding host genomes
4. Extracts host taxonomy infromation from the genome files (rule getTaxInfo)
5. Runs CRISPRCasTyper on Cas10-assocation genomes
6. Proceeds only with type III loci
7. Runs an in-house script cATyper on all loci, looking for pre-established effector proteins using HMM
8. Extract Cas10, Cas7, Cas5 and CorA fastas from type III loci. Also mark down if Cas10 contains HD or GGDD on a sequence comparison level (rule typeIII_characterizer)
9. Create HMM profiles based on pre-flagged HD-domain or GGDD-domain containing Cas10 (rules rule Cas10_HD_hmm_maker and rule Cas10_GGDD_hmm_maker). Note that different profiles are used for GGDD and GGDE seeded cyclase domains.
10. Looks for unexpected or poor E-value scoring proteins in CRISPR-Cas loci. These proteins are flagged as being potentially interesting (rule unknown_finder)
11. Various rules characterizing different proteins such as CorA and Cas7
12. Create Cas10 tree using muscle -based alignment and fasttree
13. Combine all results together into a mastertable
14. Extract loci that fulfill specific criteria (e.g. no HD, has GGDD, has potential unknown effectors in locus) and examine the unknown proteins in these loci by HMM searches against local SAVED/CARF and PFAM databases

A separate R project is then used to plot data and create models. R project currently not included here, but should be added in the future
