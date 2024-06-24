'''
This pipeline is focused in ring nucleases and transcription factors in type III CRISPR-Cas loci. 
Built on top of new_effectors.smk from a previous project
'''

project = "run1"

thread_hogger = 50 #number of threads dedicated to a single thread-hogging rule. Meant for rules that do not rely on wildcards.
thread_small = 5 #for wildcard tasks that benefit from multithreading
thread_ultrasmall = 1 #for single-thread wildcard based rules

ncbi_email = "ville.hoikkala@jyu.fi"

local_protein_blast_db = "/mnt/shared/apps/databases/ncbi/nr"

base_path = "/home/vhoikkal/scratch/private/runs/ring_nucleases" + "/" + project
program_root = "/home/vhoikkal/scratch/private/pipelines/ring_nucleases/crispr_type_iii_effector_finder"

webflags_db = base_path + "/home/vhoikkal/scratch/private/databases/webflags"

cas10_cluster_threshold = 0.9 #initial Cas10 clustering threshold
crispr_locus_interference_cutoff = 0 #cutoff for CRISPR loci interference completeness. Loci with less than this percentage of interference genes present are discarded

protein_clustering = str(config["protein_clustering"])
getGenomesBy = str(config["getGenomesBy"])
cas10_anchor = config["cas10_anchor"]
catyper_hmm_evalue = "1e-10"

hmm_msa_folder = "/home/vhoikkal/scratch/private/databases/custom_hmm/known_type_iii_effectors/profiles" #folder names within this folder are used for creating effector dictionary
hmm_database_folder = "/home/vhoikkal/scratch/private/databases/custom_hmm/known_type_iii_effectors/concatenated_profiles" #contains the concatenated and hmmpressed hmm profiles
hmm_database_file = "all.hmm" #filename for the concatenated hmmpressed hmm profiles

cas10_db = "/mnt/shared/scratch/vhoikkal/private/databases/custom_hmm/cas10/all_cas10s.hmm"
modified_cas10_hd = "/home/vhoikkal/scratch/private/databases/custom_hmm/cas10/HD_HMM.msa" #modified cas10 hmm profile for hmmsearch

carfsaved_db = "/home/vhoikkal/scratch/private/databases/carfsaved/carfsaved.hmm"

TM_path = "/mnt/shared/scratch/vhoikkal/private/databases/TM"
additional_cas_proteins = "/mnt/shared/scratch/vhoikkal/private/databases/cas_proteins"
validated_effectors_hmm_db = "/home/vhoikkal/scratch/private/databases/custom_hmm/validated_new_effectors/concatenated_profiles/all.hmm"
validated_effectors_folder = "/home/vhoikkal/scratch/private/databases/custom_hmm/validated_new_effectors"

ring_nuclease_db = "/home/vhoikkal/scratch/private/databases/custom_hmm/ring_nucleases/concatenated_profiles/all.hmm"
ring_nuclease_folder = "/home/vhoikkal/scratch/private/databases/custom_hmm/ring_nucleases"

#public databases (local)
cogs_db = "/home/vhoikkal/scratch/private/databases/cogs/COG_KOG"
pfam_db = "/home/vhoikkal/scratch/private/databases/pfam/Pfam-A.hmm"
pdb30_db = "/home/vhoikkal/scratch/private/databases/pdb30/pdb30"

temperature_data = "/home/vhoikkal/scratch/private/databases/tempura/200617_TEMPURA.csv"


genomes_folder = "/home/vhoikkal/scratch/private/databases/ncbi_genomes/bacteria/ncbi_dataset/data"
archaea_folder = "/home/vhoikkal/scratch/private/databases/ncbi_genomes/archaea/ncbi_dataset/data"
genome_count = 51920
subsampling_seed = 666

genomes_json = os.path.join(genomes_folder,"assembly_data_report.jsonl")

if cas10_anchor == False:
    prefiltering_wildcards = "05_host_genomes"
    prefiltering_host_genomes = "06_host_genomes"
elif cas10_anchor == True:
    prefiltering_wildcards = "02_host_wildcards"
    prefiltering_host_genomes = "03_host_genomes"

def aggregate_crisprcas(wildcards):
    if getGenomesBy == "local":
        checkpoint_output = checkpoints.expand_host_genomelist_to_wildcards.get(**wildcards).output[0]
        ivals = glob_wildcards(os.path.join(checkpoint_output,"{i}.txt")).i
    elif getGenomesBy == "remote":
        checkpoint_output = checkpoints.getHosts.get(**wildcards).output[0]
        ivals = glob_wildcards(os.path.join(checkpoint_output,"{sample}/{i}_genome.fna")).i
    return expand(base_path + "/07_cctyper/{i}/{i}_renaming.done", i=ivals)

def aggregate_host_genomes(wildcards):
    if getGenomesBy == "local":
        checkpoint_output = checkpoints.expand_host_genomelist_to_wildcards.get(**wildcards).output[0]
        ivals = glob_wildcards(os.path.join(checkpoint_output,"{i}.txt")).i
    elif getGenomesBy == "remote":
        checkpoint_output = checkpoints.getHosts.get(**wildcards).output[0]
        ivals = glob_wildcards(os.path.join(checkpoint_output,"{sample}/{i}_genome.fna")).i
    return expand(base_path + "/02_genome_wildcards/{i}.txt", i=ivals)

def aggregate_cas10_sequences(wildcards):
    checkpoint_output = checkpoints.type_iii_wildcarder.get(**wildcards).output[0]
    cvals = glob_wildcards(os.path.join(checkpoint_output,"{c}/cas_operons.tsv")).c
    return expand(base_path + "/09_crispr_iii_CorA/loci/{c}/{c}_Cas10.faa", c=cvals)


def aggregate_renamed_crisprs(wildcards):
    checkpoint_output = checkpoints.expand_host_genomelist_to_wildcards.get(**wildcards).output[0]
    ivals = glob_wildcards(os.path.join(checkpoint_output,"{j}.txt")).j
    return expand(base_path + "/07_cctyper/{j}/{j}_renaming.done", j=ivals)

def aggregate_download_genomes(wildcards):
    if getGenomesBy == "local":
        checkpoint_output = checkpoints.expand_host_genomelist_to_wildcards.get(**wildcards).output[0]
        ivals = glob_wildcards(os.path.join(checkpoint_output,"{i}.txt")).i
    elif getGenomesBy == "remote":
        checkpoint_output = checkpoints.getHosts.get(**wildcards).output[0]
        ivals = glob_wildcards(os.path.join(checkpoint_output,"{sample}/{i}_genome.fna")).i
    return expand(base_path + "/06_host_genomes/{i}/{i}_genome.fna", i=ivals)


    if getGenomesBy == "local":
        if cas10_anchor == True:
            checkpoint_output = checkpoints.expand_host_genomelist_to_wildcards_postcas10.get(**wildcards).output[0]
        else:
            checkpoint_output = checkpoints.expand_host_genomelist_to_wildcards.get(**wildcards).output[0]
        ivals = glob_wildcards(os.path.join(checkpoint_output,"{j}.txt")).j
    elif getGenomesBy == "remote":
        checkpoint_output = checkpoints.getHosts.get(**wildcards).output[0]
        ivals = glob_wildcards(os.path.join(checkpoint_output,"{sample}/{i}_genome.fna")).i
    print(expand(base_path + "/19_16S_seq/{j}_16S.fna", j=ivals))
    return expand(base_path + "/19_16S_seq/{j}_16S.fna", j=ivals)

def aggregate_cas10_booleans(wildcards):
    if getGenomesBy == "local":
        checkpoint_output = checkpoints.expand_host_genomelist_to_wildcards.get(**wildcards).output[0]
        ivals = glob_wildcards(os.path.join(checkpoint_output,"{i}.txt")).i
    elif getGenomesBy == "remote":
        checkpoint_output = checkpoints.getHosts.get(**wildcards).output[0]
        ivals = glob_wildcards(os.path.join(checkpoint_output,"{sample}/{i}_genome.fna")).i
    return expand(base_path + "/062_genomes_cas10/{i}/{i}_cas10_boolean.tsv", i=ivals)

def aggregate_cas10_seq_postboolean(wildcards):
    if getGenomesBy == "local":
        checkpoint_output = checkpoints.expand_host_genomelist_to_wildcards.get(**wildcards).output[0]
        ivals = glob_wildcards(os.path.join(checkpoint_output,"{i}.txt")).i
    elif getGenomesBy == "remote":
        checkpoint_output = checkpoints.getHosts.get(**wildcards).output[0]
        ivals = glob_wildcards(os.path.join(checkpoint_output,"{sample}/{i}_genome.fna")).i
    return expand(base_path + "/063_cas10_seq/{i}/{i}_cas10.faa", i=ivals)


def aggregate_cas10_sequences_prior_to_1st_clustering(wildcards):
    if getGenomesBy == "local":
        checkpoint_output = checkpoints.expand_host_genomelist_to_wildcards.get(**wildcards).output[0]
        ivals = glob_wildcards(os.path.join(checkpoint_output,"{i}.txt")).i
    elif getGenomesBy == "remote":
        checkpoint_output = checkpoints.getHosts.get(**wildcards).output[0]
        ivals = glob_wildcards(os.path.join(checkpoint_output,"{sample}/{i}_genome.fna")).i
    return expand(base_path + "/063_cas10_seq/{i}/{i}_cas10.faa", i=ivals)


def aggregate_cas10_genomes(wildcards):
    if getGenomesBy == "local":
        checkpoint_output = checkpoints.expand_host_genomelist_to_wildcards_postcas10.get(**wildcards).output[0]
        jvals = glob_wildcards(os.path.join(checkpoint_output,"{j}.txt")).j
    elif getGenomesBy == "remote":
        checkpoint_output = checkpoints.getHosts.get(**wildcards).output[0]
        jvals = glob_wildcards(os.path.join(checkpoint_output,"{sample}/{j}_genome.fna")).j
    return expand(base_path + "/03_postfiltering_genome_wildcards/{j}.txt", j=jvals)


    if getGenomesBy == "local":
        checkpoint_output = checkpoints.expand_host_genomelist_to_wildcards_postcas10.get(**wildcards).output[0]
        jvals = glob_wildcards(os.path.join(checkpoint_output,"{j}.txt")).j
    elif getGenomesBy == "remote":
        checkpoint_output = checkpoints.getHosts.get(**wildcards).output[0]
        jvals = glob_wildcards(os.path.join(checkpoint_output,"{sample}/{j}_genome.fna")).j
    return expand(base_path + "/30_clustertable/{j}/{j}.tsv", j=jvals)

def aggregate_typeIII_info(wildcards):
    checkpoint_output = checkpoints.type_iii_wildcarder.get(**wildcards).output[0]
    cvals = glob_wildcards(os.path.join(checkpoint_output,"{c}/cas_operons.tsv")).c
    return expand(base_path + "/09_crispr_iii_CorA/loci/{c}/{c}_crispr_iii_info.tsv", c=cvals)

def aggregate_taxInfo(wildcards):
    if getGenomesBy == "local":
        if cas10_anchor == True:
            checkpoint_output = checkpoints.expand_host_genomelist_to_wildcards_postcas10.get(**wildcards).output[0]
        else:
            checkpoint_output = checkpoints.expand_host_genomelist_to_wildcards.get(**wildcards).output[0]
        ivals = glob_wildcards(os.path.join(checkpoint_output,"{j}.txt")).j
    elif getGenomesBy == "remote":
        checkpoint_output = checkpoints.getHosts.get(**wildcards).output[0]
        ivals = glob_wildcards(os.path.join(checkpoint_output,"{sample}/{i}_genome.fna")).i
    return expand(base_path + "/06_host_genomes/{j}/{j}_taxon.txt", j=ivals)

def aggregate_renamed(wildcards):
    if getGenomesBy == "local":
        if cas10_anchor == True:
            checkpoint_output = checkpoints.expand_host_genomelist_to_wildcards_postcas10.get(**wildcards).output[0]
        else:
            checkpoint_output = checkpoints.expand_host_genomelist_to_wildcards.get(**wildcards).output[0]
        ivals = glob_wildcards(os.path.join(checkpoint_output,"{j}.txt")).j
    elif getGenomesBy == "remote":
        checkpoint_output = checkpoints.getHosts.get(**wildcards).output[0]
        ivals = glob_wildcards(os.path.join(checkpoint_output,"{sample}/{j}_genome.fna")).j
    return expand(base_path + "/07_cctyper/{j}/{j}_renaming.done", j=ivals)



    checkpoint_output = checkpoints.align_matcher.get(**wildcards).output[0]
    kvals = glob_wildcards(os.path.join(checkpoint_output,"{k}.done")).k
    return expand(base_path + "/42_matching_trees_distance_matrices/{k}/{k}_correlation.tsv", k=kvals)

def aggregate_cATyper_hmm(wildcards):
    checkpoint_output = checkpoints.type_iii_wildcarder.get(**wildcards).output[0]
    cvals = glob_wildcards(os.path.join(checkpoint_output,"{c}/cas_operons.tsv")).c
    return expand(base_path + "/10_cATyper_hmm/{c}/{c}_cATyper_hmm_rows.tsv", c=cvals)

def aggregate_ring_nuclease_fusions(wildcards):
    checkpoint_output = checkpoints.type_iii_wildcarder.get(**wildcards).output[0]
    cvals = glob_wildcards(os.path.join(checkpoint_output,"{c}/cas_operons.tsv")).c
    return expand(base_path + "/73_ring_nuclease_fusions/{c}/{c}_ring_nuclease_fusions.tsv", c=cvals)


def aggregate_validated_new_effectors_hmm(wildcards):
    checkpoint_output = checkpoints.type_iii_wildcarder.get(**wildcards).output[0]
    cvals = glob_wildcards(os.path.join(checkpoint_output,"{c}/cas_operons.tsv")).c
    return expand(base_path + "/60_validated_new_effectors/{c}/{c}_validated_new_effectors_hmm.tsv", c=cvals)

def aggregate_validated_new_effectors_analysis(wildcards):
    checkpoint_output = checkpoints.type_iii_wildcarder.get(**wildcards).output[0]
    cvals = glob_wildcards(os.path.join(checkpoint_output,"{c}/cas_operons.tsv")).c
    return expand(base_path + "/61_validated_new_effectors_analysis/{c}/{c}_cATyper_results.tsv", c=cvals)

def aggregate_ring_nucleases_hmm(wildcards):
    checkpoint_output = checkpoints.type_iii_wildcarder.get(**wildcards).output[0]
    cvals = glob_wildcards(os.path.join(checkpoint_output,"{c}/cas_operons.tsv")).c
    return expand(base_path + "/70_ring_nucleases/{c}/{c}_ring_nucleases_hmm.tsv", c=cvals)

def aggregate_ring_nucleases_analysis(wildcards):
    checkpoint_output = checkpoints.type_iii_wildcarder.get(**wildcards).output[0]
    cvals = glob_wildcards(os.path.join(checkpoint_output,"{c}/cas_operons.tsv")).c
    return expand(base_path + "/71_ring_nucleases_analysis/{c}/{c}_cATyper_results.tsv", c=cvals)

def aggregate_ring_nucleases_cas10_hmm(wildcards):
    checkpoint_output = checkpoints.type_iii_wildcarder.get(**wildcards).output[0]
    cvals = glob_wildcards(os.path.join(checkpoint_output,"{c}/cas_operons.tsv")).c
    return expand(base_path + "/101_ring_fusions_cas10/{c}/{c}_ring_nuclease_cas10_hmm.tsv", c=cvals)

def aggregate_cATyper_analysis(wildcards):
    checkpoint_output = checkpoints.type_iii_wildcarder.get(**wildcards).output[0]
    cvals = glob_wildcards(os.path.join(checkpoint_output,"{c}/cas_operons.tsv")).c
    return expand(base_path + "/11_cATyper_analysis/{c}/{c}_cATyper_results.tsv", c=cvals)

def aggregate_cATyper_etp(wildcards):
    checkpoint_output = checkpoints.type_iii_wildcarder.get(**wildcards).output[0]
    cvals = glob_wildcards(os.path.join(checkpoint_output,"{c}/cas_operons.tsv")).c
    return expand(base_path + "/11_cATyper_analysis/{c}/{c}_effector_to_protein.tsv", c=cvals)

def aggregate_cATyper_pte(wildcards):
    checkpoint_output = checkpoints.type_iii_wildcarder.get(**wildcards).output[0]
    cvals = glob_wildcards(os.path.join(checkpoint_output,"{c}/cas_operons.tsv")).c
    return expand(base_path + "/11_cATyper_analysis/{c}/{c}_protein_to_effector.tsv", c=cvals)

def aggregate_validated_new_effectors_etp(wildcards):
    checkpoint_output = checkpoints.type_iii_wildcarder.get(**wildcards).output[0]
    cvals = glob_wildcards(os.path.join(checkpoint_output,"{c}/cas_operons.tsv")).c
    return expand(base_path + "/61_validated_new_effectors_analysis/{c}/{c}_effector_to_protein.tsv", c=cvals)

def aggregate_validated_new_effectors_pte(wildcards):
    checkpoint_output = checkpoints.type_iii_wildcarder.get(**wildcards).output[0]
    cvals = glob_wildcards(os.path.join(checkpoint_output,"{c}/cas_operons.tsv")).c
    return expand(base_path + "/61_validated_new_effectors_analysis/{c}/{c}_protein_to_effector.tsv", c=cvals)

def aggregate_ring_nucleases_etp(wildcards):
    checkpoint_output = checkpoints.type_iii_wildcarder.get(**wildcards).output[0]
    cvals = glob_wildcards(os.path.join(checkpoint_output,"{c}/cas_operons.tsv")).c
    return expand(base_path + "/71_ring_nucleases_analysis/{c}/{c}_effector_to_protein.tsv", c=cvals)

def aggregate_ring_nucleases_pte(wildcards):
    checkpoint_output = checkpoints.type_iii_wildcarder.get(**wildcards).output[0]
    cvals = glob_wildcards(os.path.join(checkpoint_output,"{c}/cas_operons.tsv")).c
    return expand(base_path + "/71_ring_nucleases_analysis/{c}/{c}_protein_to_effector.tsv", c=cvals)


def aggregate_CorA_sequences_catyper(wildcards):
    checkpoint_output = checkpoints.type_iii_wildcarder.get(**wildcards).output[0]
    cvals = glob_wildcards(os.path.join(checkpoint_output,"{c}/cas_operons.tsv")).c
    #if a file does not exist, touch it
    list_of_coras = []
    for c in cvals:
        #check if this locus has a cas10.faa file. We don't want to study CorAs that are not associated with Cas10less loci or those
        #in which filtering removed Cas10 due to length
        #check if cas10 file exists
        #print("Size of cas10 file for locus " + c + " is " + str(os.path.getsize(base_path + "/09_crispr_iii_CorA/loci/" + c + "/" + c + "_Cas10.faa"))) 
        #if directory exists
        cas10_dir_path = base_path + "/09_crispr_iii_CorA"
        cora_dir_path = base_path + "/11_cATyper_analysis/" + c
        if (os.path.exists(cas10_dir_path)) & (os.path.exists(cora_dir_path)):
            if os.path.getsize(base_path + "/09_crispr_iii_CorA/loci/" + c + "/" + c + "_Cas10.faa") > 0:
                #print("Cas10 file is large enough for locus " + c + " so we will add it to the list of coras")
                if c == "GCF_003544875.1_0":
                    print("This is the weird one")
                list_of_coras.append(c)
                if not os.path.exists(base_path + "/11_cATyper_analysis/" + c + "/cora.faa"):
                    print("Touching an empty CorA")
                    open(base_path + "/11_cATyper_analysis/" + c + "/cora.faa", 'a').close()
            if c == "GCF_003544875.1_0":
                print("This is the weird one 2")
    print("Length of list of coras: " + str(len(list_of_coras)))
    #return the list of coras after adding the filepaths to them
    fullpath_coras = []
    for locus in list_of_coras:
        fullpath_coras.append(base_path + "/11_cATyper_analysis/" + locus + "/cora.faa")
    print(fullpath_coras)
    return fullpath_coras

def aggregate_CorA_sequences(wildcards):
    checkpoint_output = checkpoints.type_iii_wildcarder.get(**wildcards).output[0]
    cvals = glob_wildcards(os.path.join(checkpoint_output,"{c}/cas_operons.tsv")).c
    return expand(base_path + "/09_crispr_iii_CorA/loci/{c}/CorA.faa", c=cvals)


def aggregate_unknowns(wildcards):
    checkpoint_output = checkpoints.type_iii_wildcarder.get(**wildcards).output[0]
    cvals = glob_wildcards(os.path.join(checkpoint_output,"{c}/cas_operons.tsv")).c
    return expand(base_path + "/30_unknown_effectors/{c}/{c}_unknown_proteins.faa", c=cvals)

def aggregate_unknowns_locus_info(wildcards):
    checkpoint_output = checkpoints.type_iii_wildcarder.get(**wildcards).output[0]
    cvals = glob_wildcards(os.path.join(checkpoint_output,"{c}/cas_operons.tsv")).c
    return expand(base_path + "/30_unknown_effectors/{c}/{c}_locus_unknown_info.tsv", c=cvals)

def aggregate_csx19(wildcards):
    checkpoint_output = checkpoints.type_iii_wildcarder.get(**wildcards).output[0]
    cvals = glob_wildcards(os.path.join(checkpoint_output,"{c}/cas_operons.tsv")).c
    return expand(base_path + "/110_csx19/{c}/{c}_csx19.faa", c=cvals)

def aggregate_crispr_locus_proteins(wildcards):
    '''
    Aggregates the outputs of the crispr_locus_proteins rule.
    '''
    checkpoint_output = checkpoints.type_iii_wildcarder.get(**wildcards).output[0]
    cvals = glob_wildcards(os.path.join(checkpoint_output,"{c}/cas_operons.tsv")).c
    return expand(base_path + "/072_crispr_locus_proteins/{c}/{c}_crispr_locus_proteins.faa", c=cvals)

def aggregate_locus_viz(wildcards):
    '''
    Aggregates the outputs of the crispr_locus_proteins rule.
    '''
    checkpoint_output = checkpoints.type_iii_wildcarder.get(**wildcards).output[0]
    cvals = glob_wildcards(os.path.join(checkpoint_output,"{c}/cas_operons.tsv")).c
    return expand(base_path + "/90_locus_viz/{c}/{c}_viz.png", c=cvals)

def aggregate_group4_blasts(wildcards):
    '''
    Aggregates the outputs of the crispr_locus_proteins rule.
    '''
    checkpoint_output = checkpoints.type_iii_wildcarder.get(**wildcards).output[0]
    cvals = glob_wildcards(os.path.join(checkpoint_output,"{c}/cas_operons.tsv")).c
    return expand(base_path + "/group4_prober/{c}/{c}.blast", c=cvals)

def aggregate_known_effector_wildcards(wildcards):
    checkpoint_output = checkpoints.effector_wildcarder.get(**wildcards).output[0]
    effector_vals = glob_wildcards(os.path.join(checkpoint_output,"{effector}.eff")).effector
    return expand(base_path + "/45_effector_tree/{effector}_tree.txt", effector=effector_vals)

def aggregate_known_effector_wildcarder(wildcards):
    checkpoint_output = checkpoints.effector_wildcarder.get(**wildcards).output[0]
    effector_vals = glob_wildcards(os.path.join(checkpoint_output,"{effector}.eff")).effector
    return expand(base_path + "/40_known_effector_wildcards/{effector}.eff", effector=effector_vals)

def aggregate_known_effector_annotations(wildcards):
    checkpoint_output = checkpoints.effector_wildcarder.get(**wildcards).output[0]
    effector_vals = glob_wildcards(os.path.join(checkpoint_output,"{effector}.eff")).effector
    return expand(base_path + "/43_effector_hmmer_analysis/{effector}/{effector}_sorted_filtered_mapped.tsv", effector=effector_vals)

def aggregate_cATyper_hhsuite(wildcards):
    checkpoint_output = checkpoints.effector_wildcarder.get(**wildcards).output[0]
    effector_vals = glob_wildcards(os.path.join(checkpoint_output,"{effector}.eff")).effector
    return expand(base_path + "/42_effector_hhsuite/{effector}/{effector}_all.tsv", effector=effector_vals)

def aggregate_cATyper_hhsuite_parser(wildcards):
    checkpoint_output = checkpoints.effector_wildcarder.get(**wildcards).output[0]
    effector_vals = glob_wildcards(os.path.join(checkpoint_output,"{effector}.eff")).effector
    return expand(base_path + "/42_effector_hhsuite/{effector}/{effector}_hhsuite_parsed.tsv", effector=effector_vals)

def aggregate_cATyper_hhsuite_parser_cogs(wildcards):
    checkpoint_output = checkpoints.effector_wildcarder.get(**wildcards).output[0]
    effector_vals = glob_wildcards(os.path.join(checkpoint_output,"{effector}.eff")).effector
    return expand(base_path + "/42_effector_hhsuite_cogs/{effector}/{effector}_hhsuite_parsed_cogs.tsv", effector=effector_vals)

def aggregate_cora_neighbourhoods(wildcards):
    checkpoint_output = checkpoints.cora_neighbourhood_preparation.get(**wildcards).output[0]
    cora_loci = glob_wildcards(os.path.join(checkpoint_output,"{cora_locus}_crispr_locus_proteins.faa")).cora_locus
    return expand(base_path + "/52_cora_neighbour_analysis/{cora_locus}/neighbourhood_results.tsv", cora_locus=cora_loci)

def aggregate_webflags_cluster_small_unk_proteins_wildcarded(wildcards):
    checkpoint_output = checkpoints.small_unk_protein_clustered_wildcarder.get(**wildcards).output[0]
    proteinvals = glob_wildcards(os.path.join(checkpoint_output,"{protein}.txt")).protein
    return expand(base_path + "/96_small_unk_clustered_webflags_wildcards/{protein}/webflags.done", protein=proteinvals)


#Cas10 blast parameters
blast_headers = "qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle sacc"
blast_headers_group4 = "qseqid\tsseqid\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\tlocus"
cas10_e_value = "1e-20"
corA_hmm_value = "1e-20"

group4_pfam_headers = "target_name\taccession\tquery_name\taccession\tE-value\tscore\tbias\tE-value\tscore\tbias\texp\treg\tclu\tov\tenv\tdom\trep\tinc\tdescription_of_target"

print("Starting ring nuclease / transcription factor pipeline")

rule all: 
    input: base_path + "/done"
    #input: base_path + "/01_genomelist/annotated_genomelist.csv"
    #input: base_path + "/062_genomes_cas10/{i}/{i}_cas10.tsv"
    #input: base_path + "/064_cas10_clusters/cas10_all.faa"
     


rule write_down_genome_info:
    '''
    Outputs a text file genome names based on folders in genomes_folders
    '''
    output: base_path + "/01_genomelist/genomelist.txt"
    threads: thread_ultrasmall
    shell:
        """
        cd {genomes_folder}
        find * -maxdepth 0 -type d > {output}
        """

rule annotate_bacteria_and_archaea_domains:
    """
    Marks down whether a sample is an archaeon or a bacterium
    1. Create a .csv file in the **archaea folder**. This two-column file shows that each sample ID there is, is an archaeon
	2. Create another .csv file in the **bacteria folder**. Note that bacteria folder also contains the archaea (they have been copied there to simplify the fetching of genomes).
        This csv will claim that all samples, including the archaea, are bacteria.
	3. Use the archaea-csv file to **reannotate** the .csv file created in step two, making all archaea appear archaea
	4. The final output file will then be the one created in step three

    This table can then be used in subsequent rules by merging it to any output table using sample name as common identifier.
    *Note that this solution is not ideal and will break down easily if the input genomes are in a different format etc.
    So make sure that the bacteria folder contains **all** samples and that the archaea folder only contains archaeons*
    """
    output: base_path + "/01_genomelist/annotated_genomelist.csv"
    params: 
        archaea_csv = base_path + "/01_genomelist/archaea.csv",
        bacteria_csv = base_path + "/01_genomelist/bacteria.csv"
    threads: thread_ultrasmall
    shell:
        """
        cd {archaea_folder}
        find * -maxdepth 0 -type d > {params.archaea_csv}

        cd {genomes_folder}
        find * -maxdepth 0 -type d > {params.bacteria_csv}

        cd {program_root}
        python3 scripts/annotate_bacteria_and_archaea.py --archaea {params.archaea_csv} --bacteria {params.bacteria_csv} --out {output}
        """


if getGenomesBy == "local":
    checkpoint expand_host_genomelist_to_wildcards:
        '''
        Takes as input a list of genome accessions separated by newline. These genomes are then turned into wildcards.
        In random mode, uses subsampling to randomly sample a set of genomes.
        '''
        input:
            genome_info = rules.write_down_genome_info.output,
            domain_annotations = rules.annotate_bacteria_and_archaea_domains.output
        output: directory(base_path + "/" + prefiltering_wildcards)
        threads: thread_ultrasmall
        run:
            import pandas as pd
            if not os.path.exists(str(output)):
                os.makedirs(str(output))
            genomefiles = pd.read_csv(str(input.genome_info), sep = "\t", header = None)
            #print("HELLOO " + genomefiles)
            genomes = [] #final list of accessions

            #Random mode. Subsamples n genomes from all genomes randomly with seed.
            if config["genome_mode"] == "random":
                subsample_genomes = genomefiles.sample(n = genome_count, random_state=subsampling_seed)
                genomes = subsample_genomes[0]

            elif config["genome_mode"] == "all":
                genomes = genomefiles[0]

            #Taxid id mode. NCBI datasets json is probed to find representatives of a taxid.
            elif config["genome_mode"] == "taxid":
                json_df = pd.read_json(genomes_json, lines=True)
                json_df = json_df.join(pd.json_normalize(json_df['organism'])).drop('organism', axis='columns') #expand pandas dictionary column into new columns in the actual df
                taxidlistlist = []

                with open(config["taxidlistfile"]) as f: #read taxid list file
                    taxids_comma = f.read()
                    idlist = taxids_comma.split(",")

                taxidlist = pd.DataFrame(idlist) #convert list to df
                #taxidlist = pd.read_csv(config["taxidlistfile"], sep=",", header=None)
                taxidlist.columns = ["taxId"] #add column to the df
                json_df['taxId'] = json_df['taxId'].astype("string")
                taxidlist['taxId'] = taxidlist['taxId'].astype("string")
                chosen = pd.merge(json_df, taxidlist, on="taxId", how="inner")
                #chosen = json_df.loc[json_df["taxId"] == int(config["taxid"])]
                print("Shape of subsampled dataframe using TaxIDs from " + str(config["taxidlistfile"]) + ": " + str(chosen.shape))
                genomes = chosen["accession"].values.tolist()


            #add genome to list if .faa, .gff and .fna files are present
            for i in genomes:
                if (os.path.isfile(os.path.join(genomes_folder,i,"protein.faa")) and (os.path.isfile(os.path.join(genomes_folder,i,"genomic.gff")))): #make sure it's annotated
                    for fname in os.listdir(os.path.join(genomes_folder,i)): #checking fna is more complicated because the filename is not consistent
                        if fname.endswith('.fna'):
                            with open(str(output) + "/" + str(i) + '.txt', 'w') as f:
                                f.write(str(output))
        


    rule download_genomes:
        '''
        Based on the wildcards from expand_host_genomelist_to_wildcards,
        makes symlinks for the desired genomes from the genomes_folder
        into the working folder (prefiltering_host_genomes).
        '''
        input:
            ivalue = base_path + "/" + prefiltering_wildcards + "/{i}.txt",
        output: 
            fna = base_path + "/" + prefiltering_host_genomes + "/{i}/{i}_genome.fna",
            faa = base_path + "/" + prefiltering_host_genomes + "/{i}/{i}_proteins.faa",
            gff = base_path + "/" + prefiltering_host_genomes + "/{i}/{i}_features.gff"
        params: #these are mostly deprecated
            folder = base_path + "/" + prefiltering_host_genomes,
        conda: "envs/ncbidownload.yaml"
        threads: thread_ultrasmall
        log:
            out = base_path + "/logs/" + prefiltering_host_genomes + "/{i}.log",
            err = base_path + "/logs/" + prefiltering_host_genomes + "/{i}.err"
        shell:
            '''
            echo "Creating symlinks for chosen genomes"
            ln -s {genomes_folder}/{wildcards.i}/*.fna {output.fna}
            ln -s {genomes_folder}/{wildcards.i}/*.faa {output.faa}
            ln -s {genomes_folder}/{wildcards.i}/*.gff {output.gff}
            '''


elif getGenomesBy == "remote":
    species_space = str(config["species"])
    checkpoint getHosts:
        '''
        Downloads desired host genomes from NCBI using the wonderful python program ncbi-genome-download.
        Files are then gunzipped and named after their parent folder
        '''
        output: directory(base_path + "/" + prefiltering_host_genomes)
        conda: "envs/ncbidownload.yaml"
        params:
            parallels = thread_hogger,
            folder = base_path + "/06_host_genomes"
        threads: thread_ultrasmall
        shell:
            '''
            mkdir -p {params.folder}
            cd {params.folder}
            count=$(ncbi-genome-download --dry-run --genera "{species_space}" bacteria | wc -l)
            echo "Will download $count genomes from {species_space}. Note that the download will likely fail if the dataset is large (>1000 genomes). In such a case, just restart the script."
            printf 'Do you want to continue? (y/n)? '
            read -p "Do you want to continue? (y/n) " -n 1 -r
            if [[ $REPLY =~ ^[Yy]$ ]] ;then
                ncbi-genome-download --genera "{species_space}" bacteria --parallel {params.parallels} --formats fasta,protein-fasta,gff,assembly-report
                mv refseq/bacteria/* {params.folder}
                find {params.folder} -type f -name '*.gz' | tqdm | xargs gunzip
                echo renaming
                rename 's!(.*)/(.*)genomic\.fna!$1/$1_genome.fna!' */*genomic.fna
                rename 's!(.*)/(.*)protein\.faa!$1/$1_proteins.faa!' */*protein.faa
                rename 's!(.*)/(.*)genomic\.gff!$1/$1_features.gff!' */*genomic.gff
                rename 's!(.*)/(.*)assembly_report\.txt!$1/$1_report.txt!' */*assembly_report.txt
            else
                exit 1
            fi

            '''


if cas10_anchor == True: #if we filter analyzable genomes by the presence of Cas10. Uses Cas10 HMM profiles from CCtyper.
    print("Running Cas10 stuff")
    rule Cas10_genomes:
        '''
        First filter all sequences to get only >500 AA long proteins. This prevents the inclusion of
        truncated Cas10s or those that are cut artifically due to contig ending. Also reduces HMM searches.
        
        Then:
        Searches the filtered proteins of host against a user-specified preprepared hmm profile.
        Strips file of "#" -rows and checks if there is anything left. If there is, the gene has been found.
        Boolean output format:
        i,True/False,cas10_accession
        '''
        input:
            proteins = base_path + "/" + prefiltering_host_genomes + "/{i}/{i}_proteins.faa"
        output:
            hmm = base_path + "/062_genomes_cas10/{i}/{i}_cas10.tsv",
            boolean = base_path + "/062_genomes_cas10/{i}/{i}_cas10_boolean.tsv"
        params:
            proteins_filt = base_path + "/062_genomes_cas10/{i}/{i}_proteins_lengthfiltered.faa",
            out = base_path + "/062_genomes_cas10/{i}/{i}_temp.out",
            rows1 = base_path + "/062_genomes_cas10/{i}/{i}_temp_rows_1.out",
            cas10_db = cas10_db,
            rows = base_path + "/062_genomes_cas10/{i}/{i}_temp_rows.out",
            headers = base_path + "/062_genomes_cas10/{i}/{i}_headers.out",
            all_data = base_path + "/062_genomes_cas10/{i}/{i}_all_data.out",
        conda: "envs/hmmer.yaml"
        threads: thread_ultrasmall
        log:
            out = base_path + "/logs/062_genomes_cas10/{i}.log",
            err = base_path + "/logs/062_genomes_cas10/{i}.err"
        shell:
            '''
            echo "Running rule Cas10_genomes"
            echo -e 'target_name\taccession\ttlen\tquery_name\taccession\tqlen\tE-value_fullseq\tscore_fullseq\tbias_fullseq\t#_domain\tof_domain\tc-Evalue_domain\tscore_domain\tbias_domain\t_hmm_from\thmm_to\t_ali_from\tali_to\tenv_from\tenv_to\tacc\tdescription' > {output.hmm}
            if [ -s "{input.proteins}" ]; then
                cat {input.proteins} | seqkit seq -m 500 > {params.proteins_filt}
                if [ -s "{params.proteins_filt}" ]; then
                    hmmscan --domtblout {params.out} --cpu {threads} -E {corA_hmm_value} {params.cas10_db} {params.proteins_filt} 2> {log.err} 1> {log.out}
                    grep -v "#" {params.out} > {params.rows}||:
                    head -1 {params.rows} > {params.rows1}
                    echo "id,cas10_boolean,cas10_acc" > {output.boolean}
                    if [ -s {params.rows1} ]; then
                        cat {params.rows1} >> {output.hmm}
                        ACC=$(awk -F ' ' '{{print $4}}' {params.rows1})
                        echo "{wildcards.i},True","${{ACC}}" > {output.boolean}
                    else
                        echo "{wildcards.i},False,-" > {output.boolean}
                    fi
                    touch {output.hmm}
                else
                    echo "{wildcards.i},False,-" > {output.boolean}
                fi
            else
                echo "{wildcards.i},False,-" > {output.boolean}
            fi
            '''

        
    rule extract_cas10_sequence:
        input:
            boolean = rules.Cas10_genomes.output.boolean,
            proteins = base_path + "/" + prefiltering_host_genomes + "/{i}/{i}_proteins.faa"
        output: base_path + "/063_cas10_seq/{i}/{i}_cas10.faa"
        threads: thread_ultrasmall
        conda: "envs/hmmer.yaml"
        shell: 
            '''
            CAS10_ID=$(awk -F ',' '{{print $3}}' {input.boolean})
            if [ "$CAS10_ID" != "-" ]; then
                echo ">"{wildcards.i} > {output}
                seqkit grep -r -p ${{CAS10_ID}} {input.proteins} | seqkit seq -w 0 | tail -n +2 >> {output}
            else
                touch {output}
            fi
            '''

    rule concat_cluster_cas10_proteins:
        '''
        Takes all Cas10 sequences, concatenates them into one fasta and
        clusters this fasta at cas10_cluster_threshold (0-1) similarity.
        Greps only representative lines "marked by *"
        A python script extracts the accession from these lines into output.reps
        '''
        input: aggregate_cas10_seq_postboolean
        output:
            all_cas10 = base_path + "/064_cas10_clusters/cas10_all.faa",
            clusters = base_path + "/064_cas10_clusters/cas10_clust.faa.clstr",
            proteins = base_path + "/064_cas10_clusters/cas10_clust.faa",
            reps = base_path + "/064_cas10_clusters/cas10_unique_genomes.txt",
        params:
            clusterlines = base_path + "/064_cas10_clusters/clusterlines.txt"
        threads: thread_hogger
        conda: "envs/groupChar.yaml"
        shell:
            '''
            echo "concatenating cas10 sequences for clustering"
            find '{base_path}/063_cas10_seq' -maxdepth 2 -type f -wholename '*/*_cas10.faa' -print0 | xargs -0 cat > {output.all_cas10}
            echo clustering
            cd-hit -i {output.all_cas10} -o {output.proteins} -c {cas10_cluster_threshold} -n 5 -d 0 -M 16000 -T {threads}
            grep "*" {output.clusters} > {params.clusterlines}
            python scripts/getClusterRep.py --input {params.clusterlines} --output {output.reps}
            '''

    checkpoint expand_host_genomelist_to_wildcards_postcas10:
        '''
        '''
        input:
            genomes = rules.concat_cluster_cas10_proteins.output.reps
        output: directory(base_path + "/03_postfiltering_genome_wildcards")
        threads: thread_ultrasmall
        run:
            print("Expanding host list to wildcards...")

            list_of_genomes_to_remove = ["GCA_019977735.1"] #use this list to manually exclude any unwanted genomes

            if not os.path.exists(str(output)):
                os.makedirs(str(output))

            with open(str(input.genomes)) as filehandle:
                genomes = [line.rstrip('\n') for line in filehandle]

            for i in genomes:
                if i not in list_of_genomes_to_remove:
                    sample = str(str(i).split(",")[0])
                    with open(str(output) + "/" + sample + '.txt', 'w') as f:
                        f.write(sample)



    rule download_genomes_postCas10:
        '''
        After filtering by presence of Cas10, create new symlinks for such genomes.
        '''
        input: 
            genomes = base_path + "/03_postfiltering_genome_wildcards/{j}.txt",
        output: 
            fna = base_path + "/06_host_genomes/{j}/{j}_genome.fna",
            faa = base_path + "/06_host_genomes/{j}/{j}_proteins.faa",
            gff = base_path + "/06_host_genomes/{j}/{j}_features.gff"
        params: #these are mostly deprecated
            folder = base_path + "/06_host_genomes",
            zipped = base_path + "/06_host_genomes/{j}.zip",
            original_fna = base_path + "/06_host_genomes/{j}/*.fna", #for checking successful download
            original_prot = base_path + "/06_host_genomes/{j}/protein.faa",#for checking successful download
            original_gff = base_path + "/06_host_genomes/{j}/genomic.gff", #for checking successful download
        conda: "envs/ncbidownload.yaml"
        threads: thread_ultrasmall
        log:
            out = base_path + "/logs/06_host_genomes/{j}.log",
            err = base_path + "/logs/06_host_genomes/{j}.err"
        shell:
            '''
            ln -s {genomes_folder}/{wildcards.j}/*.fna {output.fna}
            ln -s {genomes_folder}/{wildcards.j}/*.faa {output.faa}
            ln -s {genomes_folder}/{wildcards.j}/*.gff {output.gff}
            '''

rule getTaxInfo:
    input: base_path + "/03_postfiltering_genome_wildcards/{j}.txt"
    output:
        taxon = base_path + "/06_host_genomes/{j}/{j}_taxon.txt"
    threads: thread_ultrasmall
    run:
        import pandas as pd
        json_df = pd.read_json(genomes_json, lines=True)
        json_df = json_df.join(pd.json_normalize(json_df['organism'])).drop('organism', axis='columns') #expand pandas dictionary column into new columns in the actual df
        json_df = json_df.join(pd.json_normalize(json_df['averageNucleotideIdentity'])).drop('averageNucleotideIdentity', axis='columns') #expand pandas dictionary column into new columns in the actual df 
        try:
            row = json_df.loc[json_df["accession"] == wildcards.j]
            species = row['submittedSpecies'].iloc[0].strip("][")
            genus = species.split(" ")[0]
            with open(output.taxon, "w") as f:
                f.write("Sample\tGenus\tSpecies\n")
                f.write(wildcards.j + "\t" + genus + "\t" + species + "\n")
        except:
            with open(output.taxon, "w") as f:
                            f.write(wildcards.j + "\tUnknown\tUnknown\n")
            pass


rule concat_taxInfo:
    input: aggregate_taxInfo
    output: base_path + "/06_host_genomes/taxInfo.txt"
    threads: thread_small
    shell:
        """
        echo "Sample\tgenus\tspecies\n" > {output}
        find '{base_path}/06_host_genomes/' -maxdepth 2 -type f -wholename '*_taxon.txt' -print0 | xargs -0 cat >> {output}
        """

rule CRISPRCasTyper:
    '''
    Types CRISPR-Cas loci in the genomes based on fasta.
    Takes output of the checkpoint. Automatically scales to the number of files
    produced by the checkpoint, so no need to call for the aggregator function here.

    '''
    input: base_path + "/06_host_genomes/{j}/{j}_genome.fna"
    output: base_path + "/07_cctyper/{j}/{j}.done"
    conda: "envs/CRISPRCasTyper.yaml"
    log: base_path + "/logs/07_cctyper/{j}/{j}_cctyper.log"
    params:
        outdir = base_path + "/07_cctyper/{j}"
    threads: thread_ultrasmall
    shell:
        '''
        rm -rf {params.outdir}
        cctyper '{input}' '{base_path}/07_cctyper/{wildcards.j}' --prodigal single 2>{log}
        touch '{output}'
        '''

rule CRISPRCasTyper_rename:
    '''
    UPDATED VERSION. This uses the CRISPR-Cas.tab file to include only fully functional CRISPR-Cas loci
    The sed adds assembly information to the results file, if it exists (which is originally based on just the nucleotide accessions)
    '''
    input: rules.CRISPRCasTyper.output
    output: base_path + "/07_cctyper/{j}/{j}_renaming.done"
    params:
        outdir = base_path + "/07_cctyper/{j}"
    threads: thread_ultrasmall
    shell:
        '''
        if test -e "{params.outdir}/CRISPR_Cas.tab";then
            mv {params.outdir}/CRISPR_Cas.tab {params.outdir}/CRISPR_Cas_old.tab
            sed '1s/$/\tassembly/; 2,$s/$/\t{wildcards.j}/' {params.outdir}/CRISPR_Cas_old.tab > {params.outdir}/CRISPR_Cas.tab
            rm {params.outdir}/CRISPR_Cas_old.tab
        fi
        touch {output}
        '''

rule concat_renamed_crisprs:
    '''
    This rule aggregates the outputs of the CRISPRCasTyper_rename rule.
    It is necessary to do this because the CRISPRCasTyper_rename rule is called by the checkpoint,
    and the checkpoint does not automatically aggregate the outputs of the rule.
    '''
    input: aggregate_renamed
    output: base_path + "/07_cctyper/renaming.done"
    threads: thread_ultrasmall
    shell:
        '''
        touch {output}
        ''' 

checkpoint type_iii_wildcarder:
    '''
    Extracts only type III loci from the CCTyper output
    This checkpoint examines outputs from CCTyper and creates new wildcards based on complete loci.
    This marks a point where the pipeline no longer looks at individual strains, but rather at the CRISPR-Cas loci individually.
    As outputs, this checkpoint extracts locus-specific information from CCTyper outputs (CRISPR-Cas.tab and cas_operons.tab)
    NOTE: now excludes all loci with >1 Cas10s
    '''
    input:
        cctyper_done = aggregate_renamed,
    output: directory(base_path + "/071_cctyper_loci")
    params:
        interference_cutoff = crispr_locus_interference_cutoff, #0-100. If the interference score is lower than this, the locus is discarded
        cctyper_folder = base_path + "/07_cctyper",
    threads: thread_hogger
    shell:
        '''
        python3 scripts/loci_wildcarder.py --input_folder {params.cctyper_folder} --output_folder {output} --interference_cutoff {params.interference_cutoff}
        '''


rule cATyper_hmm_search:
    '''
    Runs HMMscan for all proteins in each locus. As database we use the concatenated hmm database of all known effectors.
    The main output file is hmm_rows, which contains all hmm hits for all proteins in the locus.
    The column target_name is the hmm hit and query_name is the protein accession
    '''
    input:
        crispr_positive_samples = base_path + "/071_cctyper_loci/{c}/cas_operons.tsv",
    output:
        contig_proteins = base_path + "/10_cATyper_hmm/{c}/{c}_contig_proteins.faa",
        temp_rows = base_path + "/10_cATyper_hmm/{c}/{c}_temp_rows.out",
        hmm_rows = base_path + "/10_cATyper_hmm/{c}/{c}_cATyper_hmm_rows.tsv",
        cora = base_path + "/10_cATyper_hmm/{c}/CorA.faa",
    params:
        outdir = base_path + "/10_cATyper_hmm/{c}",
        hmm_msa_folder = hmm_msa_folder,
        host_genomes_folder = base_path + "/06_host_genomes",
        hmm_profile = hmm_database_folder + "/" + hmm_database_file,
        temp_hmm = base_path + "/10_cATyper_hmm/{c}/{c}_temp.out",
    conda: "envs/hmmer.yaml"
    log:
        out = base_path + "/10_cATyper_hmm/logs/{c}.out",
        err = base_path + "/10_cATyper_hmm/logs/{c}.err",
    threads: thread_ultrasmall
    shell:
        '''
        python scripts/catyper_prepper_10.py --locus {wildcards.c} --cas_operons_file {input.crispr_positive_samples} --output_folder {params.outdir} --host_genomes_folder {params.host_genomes_folder} --mode pre_hmm 2> {log.err} 1> {log.out}
        echo "Running hmmscan" >> {log.out}
        hmmscan --domtblout {params.temp_hmm} --cpu {threads} -E {catyper_hmm_evalue} {params.hmm_profile} {output.contig_proteins}  &> /dev/null
        echo "Removing commented rows" >> {log.out}
        grep -v "^#" {params.temp_hmm} > {output.temp_rows} ||:
        echo "Writing header" >> {log.out}
        echo -e 'target_name\taccession\ttlen\tquery_name\taccession\tqlen\tE-value_fullseq\tscore_fullseq\tbias_fullseq\t#_domain\tof_domain\tc-Evalue_domain\ti-Evalue_domain\tscore_domain\tbias_domain\t_hmm_from\thmm_to\t_ali_from\tali_to\tenv_from\tenv_to\tacc\tdescription' > {output.hmm_rows}
        echo "Checking if hits were found" >> {log.out}
        if [ -s {output.temp_rows} ]; then 
            echo "Hits found for {wildcards.c}" >> {log.out}
            cat {output.temp_rows} >> {output.hmm_rows}
        else
            echo "No hits found for {wildcards.c}" >> {log.out}
            touch {output.hmm_rows}
        fi
        touch {output.cora}
        '''

rule concatenate_cATyper_hmm:
    input: aggregate_cATyper_hmm
    output: base_path + "/10_cATyper_hmm/cATyper_all.tsv"
    threads: thread_ultrasmall
    shell:
        """
        awk '(NR == 1) || (FNR > 1)' {input} > {output}
        """

rule cATyper_analysis:
    '''
    This is the analysis step of cATyper. Note that we are still using the same python script as in the previous rule,
    but the mode is now set to "post_hmm".
    '''
    input:
        crispr_positive_samples = base_path + "/071_cctyper_loci/{c}/cas_operons.tsv",
        contig_proteins = base_path + "/10_cATyper_hmm/{c}/{c}_contig_proteins.faa",
        hmm_rows = base_path + "/10_cATyper_hmm/{c}/{c}_cATyper_hmm_rows.tsv",
    output:
        catyper = base_path + "/11_cATyper_analysis/{c}/{c}_cATyper_results.tsv",
        hmm_targets = base_path + "/11_cATyper_analysis/{c}/{c}_cATyper_hmm_targets.tsv",
        effector_to_protein = base_path + "/11_cATyper_analysis/{c}/{c}_effector_to_protein.tsv",
        protein_to_effector = base_path + "/11_cATyper_analysis/{c}/{c}_protein_to_effector.tsv",
        plottable_effector_positions = base_path + "/11_cATyper_analysis/{c}/{c}_plottable_effector_positions.tsv",
    params:
        outdir = base_path + "/11_cATyper_analysis/{c}",
        host_genomes_folder = base_path + "/06_host_genomes",
        hmm_msa_folder = hmm_msa_folder,
        tmhmm_model_path = TM_path + "/TMHMM2.0.model",
    conda: "envs/hmmer.yaml"
    log:
        out = base_path + "/11_cATyper_analysis/logs/{c}.out",
        err = base_path + "/11_cATyper_analysis/logs/{c}.err",
    threads: thread_ultrasmall
    shell:
        '''
        echo "Running cATyper analysis" >> {log.out}
        echo "Listing targets" >> {log.out}
        ls {params.hmm_msa_folder}/*/*.hmm > {output.hmm_targets}
        python scripts/catyper_prepper_10.py --locus {wildcards.c} --cas_operons_file {input.crispr_positive_samples} --output_folder {params.outdir} --host_genomes_folder {params.host_genomes_folder} --mode post_hmm --hmm_targets {output.hmm_targets} --hmm_rows {input.hmm_rows} --catyper_out {output.catyper} --catyper_type "known_effector" --effector_plot_data {output.plottable_effector_positions} --tmhmm_model_path {params.tmhmm_model_path} 2> {log.err} 1> {log.out}
        '''

rule concatenate_cATyper_analysis:
    input: aggregate_cATyper_analysis
    output: base_path + "/11_cATyper_analysis/cATyper_all.tsv"
    threads: thread_ultrasmall
    shell:
        """
        awk '(NR == 1) || (FNR > 1)' {base_path}/11_cATyper_analysis/*/*_cATyper_results.tsv > {output}
        touch {output}
        """

rule concatenate_cATyper_analysis_effector_scores:
    '''
    Concatenates info on the effector hits from cATyper analysis
    '''
    input:
        protein_to_effector = aggregate_cATyper_pte,
        effector_to_protein = aggregate_cATyper_etp
    output:
        protein_to_effector_concatenated = base_path + "/11_cATyper_analysis/cATyper_protein_to_effector.tsv",
        effector_to_protein_concatenated = base_path + "/11_cATyper_analysis/cATyper_effector_to_protein.tsv",
    threads: thread_small
    shell:
        '''
        cat {input.protein_to_effector} > {output.protein_to_effector_concatenated}
        cat {input.effector_to_protein} > {output.effector_to_protein_concatenated}
        '''

rule analyse_cATyper_effector_scores:
    '''
    The algorithm that analyses how well each hmm profile performs in cATyper.
    '''
    input:
        pte = rules.concatenate_cATyper_analysis_effector_scores.output.protein_to_effector_concatenated,
        etp = rules.concatenate_cATyper_analysis_effector_scores.output.effector_to_protein_concatenated,
    output:
        effector_scores = base_path + "/11_cATyper_analysis/cATyper_effector_scores.tsv",
        effector_scores_summary = base_path + "/11_cATyper_analysis/cATyper_effector_scores_summary.tsv",
        #effector_scores_plot1 = base_path + "/11_cATyper_analysis/cATyper_effector_scores_plot1.png",
    params:
        outdir = base_path + "/11_cATyper_analysis"
    conda: "envs/hmmer.yaml"
    threads: thread_hogger
    shell:
        '''
        python scripts/catyper_effector_scores.py --pte {input.pte} --etp {input.etp} --output {output.effector_scores} --output_summary {output.effector_scores_summary} --outdir {params.outdir}
        touch {output.effector_scores}
        touch {output.effector_scores_summary}
        '''

checkpoint effector_wildcarder:
    '''
    This rule creates wildcards for each previosly known effector.
    Wildcards are used downstream when creating phylogenetic trees etc for each effector.
    '''
    output:
        directory(base_path + "/40_known_effector_wildcards")
        #hmm_targets = base_path + "/40_known_effector_wildcards/{effector}/{effector}.txt"
    params:
        hmm_targets = base_path + "/40_known_effector_wildcards/alltargets.txt",
        hmm_msa_folder = hmm_msa_folder,
        outdir = base_path + "/40_known_effector_wildcards",
    threads: thread_small
    run:
        import glob
        #Equivalent of `mkdir`
        Path(params.outdir).mkdir(parents=True, exist_ok=True)
        
        #Equivalent of `ls {hmm_msa_folder}/*/*.hmm > {hmm_targets}`
        with open(params.hmm_targets, 'w') as f:
            files = glob.glob(f"{params.hmm_msa_folder}/*/*.hmm") #get all files in the hmm_msa_folder
            f.write('\n'.join(files))

        with open(params.hmm_targets, 'r') as f:
            lines = f.readlines()
            for line in lines:
                filename = os.path.basename(line) # get filename from path
                effector = re.split("\.", filename)[0] #remove extension
                if "#" not in line:
                    effector = re.split("_", effector)[1]
                elif "#" in line: #in case a single effector is characterised by multiple hmm profiles, we mark the filenames with # followed by the effector, e.g ca6_csm6-ca6_italicus#csm6-ca6
                    effector = re.split("#", effector)[1]
                print(effector)
                with open(f'{params.outdir}/{effector}.eff', 'w') as f:
                    f.write(effector)

rule effector_fetcher:
    '''
    Fetches protein sequences for an effector from catyper_analysis outputs
    '''
    input:
        aggregate_known_effector_wildcarder,
        catyper_done = rules.concatenate_cATyper_analysis.output, #this indicates we can safely probe the folders for effector sequences
    output:
        multifasta = base_path + "/41_known_effector_mf/{effector}.faa",
    threads: thread_ultrasmall
    shell:
        '''
        effector_name={wildcards.effector}
        echo "Fetching effector sequences for {wildcards.effector}"
        touch {output.multifasta}
        find {base_path}/11_cATyper_analysis/ -name \"*${{effector_name}}*.faa\" -exec cat {{}} \; >> {output.multifasta}
        '''

rule effector_hmmer:
    '''
    Uses hmmscan against pfam to annotate the discovered effectors' domains.
    Output used when plotting tree with R
    '''
    input:
        multifasta = rules.effector_fetcher.output.multifasta
    output:
        raw_hmm = base_path + "/42_effector_hmmer/{effector}.hmm",
        #filtered_hmm_tsv = base_path + "/42_effector_hmmer/{effector}.tsv",
    params:
        base_path = base_path + "/42_effector_hmmer",
        pfam_db = pfam_db
    threads: thread_small
    conda: "envs/hmmer.yaml"
    log:
        out = base_path + "/42_effector_hmmer/logs/{effector}.out",
        err = base_path + "/42_effector_hmmer/logs/{effector}.err",
    shell:
        '''
        python3 scripts/effector_hmmer.py --input {input.multifasta} --output_basepath {params.base_path} --effector {wildcards.effector} --pfam_path {params.pfam_db} 2> {log.err} 1> {log.out}
        '''

rule cATyper_hmm_search_hhsuite_pdb:
    '''
    HMM search against the PDB database using HHsuite.
    Using a custom made wrapper to first divide the multifasta into separate fastas for each protein as required by hhblits.
    The wrapper divides a multifasta (wildcard {c}) into its constituent proteins and then runs hhblits on each protein.
    The resulting files are .tsv files, each corresponding to a single protein within the effector wildcard (e.g. RelE protein X)
    These individual .tsv files are then concatenated into a single .tsv file for each effector and fed to the parse_hhsuite rule.
    '''
    input:
        multifasta = rules.effector_fetcher.output.multifasta
    output:
        hhsuite_concat = base_path + "/42_effector_hhsuite/{effector}/{effector}_all.tsv",
    params:
        outdir = base_path + "/42_effector_hhsuite/{effector}",
        pdb30 = pdb30_db + "/pdb30",
    conda: "envs/hhsuite.yaml"
    threads: thread_small
    log:
        out = base_path + "/42_effector_hhsuite/logs/{effector}.out",
        err = base_path + "/42_effector_hhsuite/logs/{effector}.err",
    shell:
        '''
        python3 scripts/hhblits_wrapper.py --input {input.multifasta} --output_basepath {params.outdir} --database {params.pdb30}
        cat {params.outdir}/hhblits/*.tsv > {output.hhsuite_concat}
        '''

rule parse_hhsuite:
    '''
    Parser output from HHSuite (PDB)
    '''
    input: rules.cATyper_hmm_search_hhsuite_pdb.output.hhsuite_concat
    output: base_path + "/42_effector_hhsuite/{effector}/{effector}_hhsuite_parsed.tsv"
    conda: "envs/hhsuite.yaml"
    threads: thread_small
    params:
        database = "PDB"
    shell:
        '''
        python3 scripts/hhsuite_parser.py --infile {input} --outfile {output}  --database {params.database}
        '''

rule concatenate_cATyper_hmm_hhsuite:
    input: aggregate_cATyper_hhsuite_parser
    output: base_path + "/42_effector_hhsuite/cATyper_all.tsv"
    threads: thread_small
    shell:
        """
        awk '(NR == 1) || (FNR > 1)' {input} > {output}
        """

rule cATyper_hmm_search_hhsuite_cog:
    '''
    HMM search against the COG database using HHsuite.
    Using a custom made wrapper to first divide the multifasta into separate fastas for each protein as required by hhblits.
    The wrapper divides a multifasta (wildcard {c}) into its constituent proteins and then runs hhblits on each protein.
    The resulting files are .tsv files, each corresponding to a single protein within the effector wildcard (e.g. RelE protein X)
    These individual .tsv files are then concatenated into a single .tsv file for each effector and fed to the parse_hhsuite rule.
    '''
    input:
        multifasta = rules.effector_fetcher.output.multifasta
    output:
        hhsuite_concat = base_path + "/42_effector_hhsuite_cogs/{effector}/{effector}_all.tsv",
    params:
        outdir = base_path + "/42_effector_hhsuite_cogs/{effector}",
        cogs = cogs_db,
    conda: "envs/hhsuite.yaml"
    threads: thread_small
    log:
        out = base_path + "/42_effector_hhsuite_cogs/logs/{effector}.out",
        err = base_path + "/42_effector_hhsuite_cogs/logs/{effector}.err",
    shell:
        '''
        python3 scripts/hhblits_wrapper.py --input {input.multifasta} --output_basepath {params.outdir} --database {params.cogs} 2> {log.err} 1> {log.out}
        #check if any .tsv files exist in {params.outdir}/hhblits
        if [ -s {params.outdir}/hhblits/*.tsv ]; then
            cat {params.outdir}/hhblits/*.tsv > {output.hhsuite_concat}
        else
            touch {output.hhsuite_concat}
        fi
        '''

rule parse_hhsuite_cogs:
    input: rules.cATyper_hmm_search_hhsuite_cog.output.hhsuite_concat
    output: base_path + "/42_effector_hhsuite_cogs/{effector}/{effector}_hhsuite_parsed_cogs.tsv"
    conda: "envs/hhsuite.yaml"
    params:
        database = "COGs",
        mapping = cogs_db + "/cog-20.def.tab"
    threads: thread_small
    log:
        out = base_path + "/42_effector_hhsuite_cogs/logs/{effector}.out",
        err = base_path + "/42_effector_hhsuite_cogs/logs/{effector}.err",
    shell:
        '''
        python3 scripts/hhsuite_parser.py --infile {input} --outfile {output} --database {params.database} --mapping {params.mapping} 2> {log.err} 1> {log.out}
        '''

rule concatenate_cATyper_hmm_hhsuite_cogs:
    input: aggregate_cATyper_hhsuite_parser_cogs
    output: base_path + "/42_effector_hhsuite_cogs/cATyper_all_hmm_cogs.tsv"
    threads: thread_small
    shell:
        """
        awk '(NR == 1) || (FNR > 1)' {input} > {output}
        """

rule effector_analyse_domains:
    '''
    Takes the HMM output from effector_hmmer and creates a simplified table of domains for each effector.
    This is only for the Pfam analysis.
    LEGACY: Links protein ID to each effector (the original effector fastas previously contained locus IDs, not protein IDs)
    '''
    input:
        raw_hmm = rules.effector_hmmer.output.raw_hmm,
        multifasta = rules.effector_fetcher.output.multifasta
    output:
        domains = base_path + "/43_effector_hmmer_analysis/{effector}/{effector}_sorted_filtered_mapped.tsv",
        filtered_hmm_tsv = base_path + "/43_effector_hmmer_analysis/{effector}/{effector}_sorted_filtered.tsv",
        sorted_hmm_tsv = base_path + "/43_effector_hmmer_analysis/{effector}/{effector}_sorted.tsv"
    params:
        base_path = base_path + "/43_effector_hmmer_analysis",
        effector_locus_map = base_path + "/43_effector_hmmer_analysis/{effector}/{effector}.locus_map.tsv",
    threads: thread_small
    log:
        out = base_path + "/43_effector_hmmer_analysis/logs/{effector}.out",
        err = base_path + "/43_effector_hmmer_analysis/logs/{effector}.err"
    shell:
        '''
        cat {base_path}/11_cATyper_analysis/*/*_protein_to_effector.tsv | grep "{wildcards.effector}" > {params.effector_locus_map}
        python3 scripts/effector_hmmer_analyse.py --input {input.raw_hmm} --multifasta {input.multifasta} --output_basepath {params.base_path} --effector {wildcards.effector} --effector_locus_map {params.effector_locus_map} 2> {log.err} 1> {log.out}
        '''

rule effector_align:
    input: rules.effector_fetcher.output.multifasta
    output: base_path + "/44_effector_alignments/{effector}.afa",
    conda: "envs/trees.yaml"
    log:
        out = base_path + "/44_effector_alignments/logs/{effector}.out",
        err = base_path + "/44_effector_alignments/logs/{effector}.err"
    threads: thread_small
    shell:
        '''
        muscle -super5 "{input}" -output "{output}" -threads {threads} 2> {log.err} 1> {log.out} 
        '''


rule effector_tree:
    input: rules.effector_align.output
    output: base_path + "/45_effector_tree/{effector}_tree.txt",
    conda: "envs/trees.yaml"
    threads: thread_small
    shell:
        '''
        FastTree -wag -gamma {input} > {output}
        '''


rule concatenate_effector_wildcards:
    input:
        trees = aggregate_known_effector_wildcards,
        annotations = aggregate_known_effector_annotations,
        checkpoint = aggregate_known_effector_wildcarder
    output: base_path + "/known_effectors_finished.txt"
    threads: thread_small
    shell:
        '''
        cat {input} > {output}
        '''

rule crispr_locus_proteins:
    '''
    Extracts all proteins from a CRISPR-Cas locus +- 4000 bp
    '''
    input:
        cas_operons_tsv = base_path + "/071_cctyper_loci/{c}/cas_operons.tsv",
    output:
        crispr_locus_proteins = base_path + "/072_crispr_locus_proteins/{c}/{c}_crispr_locus_proteins.faa",
    params:
        output_folder = base_path + "/072_crispr_locus_proteins/{c}",
        locus = "{c}",
        host_genomes_folder = base_path + "/06_host_genomes",
    threads: thread_ultrasmall
    run:
        import os
        import pandas as pd
        from Bio import SeqIO
        import gffutils
        import re
        import sys

        #when looking for effectors in the locus, the effector search range is the number of bases up or downstream of the cctyper defined cas operon boundaries
        effector_search_range = 4000

        #get sample name by splitting the locus name at the underscore and taking the first two parts
        sample = params.locus.split("_")[0] + "_" + params.locus.split("_")[1]

        #read in the cas operons file
        cas_operons_df = pd.read_csv(input.cas_operons_tsv, sep ="\t", header = 0)

        #find the gff file for the sample from the host_genomes_folder/samples folder
        gff_file = os.path.join(params.host_genomes_folder, sample, sample + "_features.gff")

        #Using gffutils, create a database of the gff file
        db_path = os.path.join(params.output_folder, params.locus)

        #get the path to the proteins fasta file for the current sample
        proteins_fasta = os.path.join(params.host_genomes_folder, sample, sample + "_proteins.faa")

        #create folder for the above path if it does not exist
        #if not os.path.exists(os.path.dirname(db_path)):
        #    print("Creating folder " + str(os.path.dirname(db_path)))
        #    os.makedirs(os.path.dirname(params.db_path))
        db = gffutils.create_db(gff_file, 
                                dbfn=db_path + ".db", 
                                force=True, 
                                keep_order=True, 
                                merge_strategy='merge', 
                                sort_attribute_values=True,
                                id_spec=["protein_id", "Name", "ID"])

        #read in the db
        db = gffutils.FeatureDB(db_path + ".db", keep_order=True)

        #get contig from row 0, column Contig
        contig = cas_operons_df.iloc[0]["Contig"]

        cas_operon_start = int(cas_operons_df['Start'][0]) - effector_search_range
        cas_operon_end = int(cas_operons_df['End'][0]) + effector_search_range

        #from the gff file, extract all the features that are on the contig. The featuretype must be "CDS"
        protein_ids = []
        proteins_on_contig = db.region(featuretype='CDS', seqid=contig)

        #then, extract all proteins whose coordinates are between the start and end of the cas operon +- the effector search range
        proteins_on_crispr_locus = db.region(featuretype='CDS', seqid=contig, start = cas_operon_start, end = cas_operon_end)

        #convert the returned generator to a list
        print("Extracting protein IDs from gff file")
        for i in proteins_on_crispr_locus:
            #check if the attributes of the feature has a key called 'protein_id'
            if "protein_id" in i.attributes:
                id_in_gff = str(i.attributes['protein_id'][0])
                #id = str(i.attributes['ID'][0]).split("-")[1]
                protein_ids.append(id_in_gff)

        #using biopython, extract the protein sequences using the list protein_ids from the proteins fasta file
        #the proteins fasta file is in the same folder as the gff file
        protein_seqs = []
        print("Reading protein fasta file from " + str(proteins_fasta))
        for record in SeqIO.parse(proteins_fasta, "fasta"):
            if record.id in protein_ids:
                protein_seqs.append(record)

        #write the protein sequences to a multifasta file
        print("Writing protein sequences to " + str(output.crispr_locus_proteins))
        SeqIO.write(protein_seqs, output.crispr_locus_proteins, "fasta")



rule aggregate_crispr_locus_proteins:
    '''
    Aggregates the outputs of the crispr_locus_proteins rule.
    '''
    input: aggregate_crispr_locus_proteins
    output: base_path + "/072_crispr_locus_proteins/crispr_locus_proteins_all.faa"
    threads: thread_small
    shell:
        '''
        cat {input} > {output}
        '''



rule typeIII_characterizer:
    '''
    Characterizes previously known genes from type III loci (currently Cas10, Cas7, Cas5, CorA).
    This version incorporates the HD domain search. HD domains are searched for in the Cas10 protein
    and within the first 5 to 30 residues. Also adds Cas10 length.
    Now also looks for GGDD domain in Cas10.
    '''
    input:
        cctyper = base_path + "/071_cctyper_loci/{c}/cas_operons.tsv", #the cas_operons file is trimmed down to only include this specific {c} locus
        #gff = base_path + "/06_host_genomes/{j}/{j}_features.gff",
        #proteins = base_path + "/06_host_genomes/{j}/{j}_proteins.faa",
    output:
        Cas10_fasta = base_path + "/09_crispr_iii_CorA/loci/{c}/{c}_Cas10.faa",
        Cas5_fasta = base_path + "/09_crispr_iii_CorA/loci/{c}/{c}_Cas5.faa",
        Cas7_fasta = base_path + "/09_crispr_iii_CorA/loci/{c}/{c}_Cas7.faa",
        #CorA_fasta = base_path + "/09_crispr_iii_CorA/loci/{c}/{c}_CorA.faa", deprecated as of 17.1.2024
        info = base_path + "/09_crispr_iii_CorA/loci/{c}/{c}_crispr_iii_info.tsv",
    params:
        this_folder = base_path + "/09_crispr_iii_CorA/loci/{c}",
        outputfolder = base_path + "/09_crispr_iii_CorA/loci/{c}",
        sample_folder = base_path + "/06_host_genomes/", #this is used to find the gff and proteins files that are strain-, not locus-specific
        cctyper_folder = base_path + "/07_cctyper"
    conda: "envs/gff_utils.yaml"
    log:
        out = base_path + "/09_crispr_iii_CorA/logs/{c}.out",
        err = base_path + "/09_crispr_iii_CorA/logs/{c}.err"
    threads: thread_ultrasmall
    shell:
        '''
        python3 scripts/type_iii_effector_finder_2.0_HD.py --locus_id {wildcards.c} --sample_folder {params.sample_folder} --this_folder {params.this_folder} --outputfolder {params.outputfolder} --cas_operons {input.cctyper} --info_out {output.info} --cctyper_path {params.cctyper_folder} 2> {log.err} 1> {log.out}
        touch {output.Cas10_fasta}
        touch {output.Cas5_fasta}
        touch {output.Cas7_fasta}
        '''

# rule Cas10_fusion_finder:
#     '''
#     Uses the same HMM library as cAtyper to look for additional domains in Cas10 proteins.
#     '''
#     input:
#         Cas10_fasta = base_path + "/09_crispr_iii_CorA/loci/{c}/{c}_Cas10.faa",
#     output:
#         catyper = base_path + "/10_cATyper_hmm/{c}/{c}_cATyper_results.tsv",
#         hmm_targets = base_path + "/10_cATyper_hmm/{c}/{c}_cATyper_hmm_targets.tsv",    
#         hmm_rows = base_path + "/10_cATyper_hmm/{c}/{c}_cATyper_hmm_rows.tsv",
#         temp_rows = base_path + "/10_cATyper_hmm/{c}/{c}_temp_rows.out",
#     params:
#         outdir = base_path + "/10_cATyper_hmm/{c}",
#         hmm_msa_folder = hmm_msa_folder,
#         host_genomes_folder = base_path + "/06_host_genomes",
#         hmm_profile = hmm_database_folder + "/" + hmm_database_file,
#         temp_hmm = base_path + "/10_cATyper_hmm/{c}/{c}_temp.out",
#         temp_rows = base_path + "/10_cATyper_hmm/{c}/{c}_temp_rows.out",
#     conda: "envs/hmmer.yaml"
#     log:
#         out = base_path + "/10_cATyper_hmm/logs/{c}.out",
#         err = base_path + "/10_cATyper_hmm/logs/{c}.err",
#         out2 = base_path + "/10_cATyper_hmm/logs/{c}.out2",
#         err2 = base_path + "/10_cATyper_hmm/logs/{c}.err2",
#     shell:
#         '''
#         python scripts/catyper_prepper_10.py --locus {wildcards.c} --cas_operons_file {input.crispr_positive_samples} --output_folder {params.outdir} --host_genomes_folder {params.host_genomes_folder} --mode pre_hmm 2> {log.err} 1> {log.out}
#         echo "Running hmmscan" >> {log.out}
#         hmmscan --domtblout {params.temp_hmm} --cpu 8 -E {catyper_hmm_evalue} {params.hmm_profile} {output.contig_proteins}  &> /dev/null
#         echo "Removing commented rows" >> {log.out}
#         grep -v "#" {params.temp_hmm} > {output.temp_rows} ||:
#         echo "Writing header" >> {log.out}
#         echo -e 'target_name\taccession\ttlen\tquery_name\taccession\tqlen\tE-value_fullseq\tscore_fullseq\tbias_fullseq\t#_domain\tof_domain\tc-Evalue_domain\ti-Evalue_domain\tscore_domain\tbias_domain\t_hmm_from\thmm_to\t_ali_from\tali_to\tenv_from\tenv_to\tacc\tdescription' > {output.hmm_rows}
#         echo "Checking if hits were found" >> {log.out}
#         if [ -s {output.temp_rows} ]; then 
#             echo "Hits found for {wildcards.c}" >> {log.out}
#             cat {output.temp_rows} >> {output.hmm_rows}
#         else
#             echo "No hits found for {wildcards.c}" >> {log.out}
#             touch {output.hmm_rows}
#         fi

#         echo "Listing targets" >> {log.out}
#         ls {params.hmm_msa_folder}/*/*.hmm > {output.hmm_targets}

#         echo "Running cATyper" >> {log.out}
#         python scripts/catyper_prepper_10.py --locus {wildcards.c} --cas_operons_file {input.crispr_positive_samples} --output_folder {params.outdir} --host_genomes_folder {params.host_genomes_folder} --mode post_hmm --hmm_targets {output.hmm_targets} --hmm_rows {output.hmm_rows} --catyper_out {output.catyper}  2> {log.err2} 1> {log.out2}
#         touch {output.cora}
#         '''

rule unknown_finder:
    '''
    Probes the CCTyper outputs for unknown genes and/or genes with low E-values in the HMM search.
    In the future, the info file could contain the following columns:
        - locus_id
        - protein id
        - protein length
        - protein sequence
        - presence of CARF/SAVED domain
        - presence of HEPN domain
        - transmembrane properties
    Also outputs a fasta file with the unknown proteins. The output files are completely empty if no unknown proteins are found.
    
    For transmembrane prediction, we use TMHMM. This is not a conda package, so is installed
    using pip. The pip freeze command + grep is used to check if the package is already installed.
    NOTE: TMHMM will fail out of the box. A file called /home/ubuntu/miniconda3/envs/environmentname/lib/python3.9/site-packages/tmhmm/api.py in the conda installation will need this 
    to bypass a non-crucial warning that causes the crash:
    import warnings
        warnings.filterwarnings('ignore', category=RuntimeWarning)
    '''
    input:
        cctyper = base_path + "/071_cctyper_loci/{c}/cas_operons.tsv",
        hmm = rules.cATyper_hmm_search.output.temp_rows
    output:
        unknown_proteins = base_path + "/30_unknown_effectors/{c}/{c}_unknown_proteins.faa",
        info = base_path + "/30_unknown_effectors/{c}/{c}_unknown_proteins_info.tsv", #each row is a protein
        locus_info = base_path + "/30_unknown_effectors/{c}/{c}_locus_unknown_info.tsv", #each row is a locus
    params:
        outputfolder = base_path + "/30_unknown_effectors/{c}",
        host_genomes_folder = base_path + "/06_host_genomes",
        cctyper_folder = base_path + "/07_cctyper",
        sample_folder = base_path + "/06_host_genomes/", #this is used to find the gff and proteins files that are strain-, not locus-specific
        extra_cas_db = additional_cas_proteins + "/new_effectors_additional_cas.dmnd",
    conda: "envs/gff_utils.yaml"
    log:
        out = base_path + "/30_unknown_effectors/logs/{c}.out",
        err = base_path + "/30_unknown_effectors/logs/{c}.err"
    threads: thread_ultrasmall
    shell:
        '''
        pip freeze | grep tmhmm.py || pip install tmhmm.py
        python3 scripts/unknown_finder.py --locus_id {wildcards.c} --cctyper {input.cctyper} --hmm {input.hmm} --additional_cas_db {params.extra_cas_db} --outputfolder {params.outputfolder} --info_out {output.info} --unknown_proteins_output {output.unknown_proteins} --unknown_proteins_output_150 {output.unknown_proteins_150}--host_genomes_folder {params.host_genomes_folder} --cctyper_path {params.cctyper_folder} --sample_folder {params.sample_folder} --output_locus_info {output.locus_info} 2> {log.err} 1> {log.out}
        touch {output.unknown_proteins}
        touch {output.info}
        '''

rule concatenate_unknowns:
    input: aggregate_unknowns
    output:
        proteins = base_path + "/30_unknown_effectors/unknowns.faa",
        info = base_path + "/30_unknown_effectors/unknowns_info.tsv"
    threads: thread_ultrasmall
    shell:
        '''
        find '{base_path}/30_unknown_effectors' -maxdepth 2 -type f -wholename '*/*_unknown_proteins.faa' -print0 | xargs -0 cat >> {output.proteins}
        echo "locus_id\tsample\tsequence\tlength\tevalue\tposition\tid\tcctyper" > {output.info}
        find '{base_path}/30_unknown_effectors' -maxdepth 2 -type f -wholename '*/*_unknown_proteins_info.tsv' -print0 | xargs -0 cat >> {output.info}
        '''

rule concatenate_unknowns_locus_info:
    '''
    Aggregates the locus specific regarding unknown effectors info files
    '''
    input: aggregate_unknowns_locus_info
    output:
        info = base_path + "/30_unknown_effectors/unknowns_info_loci.tsv"
    threads: thread_ultrasmall
    shell:
        '''
        echo "locus_id\tsample\tno_of_unknowns\tunknown_proteins" > {output.info}
        find '{base_path}/30_unknown_effectors' -maxdepth 2 -type f -wholename '*/*_locus_unknown_info.tsv' -print0 | xargs -0 cat >> {output.info}
        '''

rule cluster_unknowns:
    '''
    Makes clusters of the unknown proteins using CD-HIT.
    Choose of word size:
        -n 5 for thresholds 0.7 ~ 1.0
        -n 4 for thresholds 0.6 ~ 0.7
        -n 3 for thresholds 0.5 ~ 0.6
        -n 2 for thresholds 0.4 ~ 0.5
    '''
    input: rules.concatenate_unknowns.output.proteins
    output:
        proteins = base_path + "/31_unknowns_cluster/unknowns_cluster.faa",
        clusterinfo = base_path + "/31_unknowns_cluster/unknowns_cluster.faa.clstr",
    conda: "envs/trees.yaml"
    log:
        out = base_path + "/31_unknowns_cluster/logs/unknowns_cluster.out",
        err = base_path + "/31_unknowns_cluster/logs/unknowns_cluster.err"
    params:
        cluster_cutoff = 0.4
    threads: thread_hogger
    shell:
        '''
        cd-hit -i {input} -o {output.proteins} -c {params.cluster_cutoff} -n 2 -d 0 -M 16000 -T {threads}
        '''


rule concatenate_type_iii_info:
    '''
    1. Creates header in output
    2. Concatenates individual loci info files
    3. Removes the headers from the concatenated file using inverse grep
    4. Joins this modified file with the header file
    '''
    input: aggregate_typeIII_info
    output: base_path + "/09_crispr_iii_CorA/loci/type_iii_info.tsv"
    params:
        temp_out = base_path + "/09_crispr_iii_CorA/loci/temp.tsv" 
    threads: thread_ultrasmall
    shell:
        '''
        echo "Cas10\tCas5\tCas7\tLocus\tSample\tCas10_GGDD\tCas10_GGDD_coord\tCas10_GGDD_seq\tCas10_GGDE\tCas10_GGDE_coord\tCas10_GGDE_seq\tCas10_GGED\tCas10_GGED_coord\tCas10_GGED_seq\tCas10_HD\tCas10_HD_list\tCas10_DH\tCas10_HD_coord\tCas10_DH_coord\tCas10_coord\tCas10_length\tSubtype\tcyclase_literal" > {params.temp_out}
        find '{base_path}/09_crispr_iii_CorA/loci' -maxdepth 2 -type f -wholename '*/*_crispr_iii_info.tsv' -print0 | xargs -0 tail -n +2 >> {params.temp_out}
        grep -v "==>" {params.temp_out} > {output}
        '''

rule Cas10_concatenate:
    '''
    Concatenates Cas10 sequences output by rule concatenate_type_iii_info
    '''
    input: aggregate_cas10_sequences
    output: base_path + "/10_Cas10_cluster/cas10s.faa"
    threads: thread_ultrasmall
    shell:
        '''
        cat {input} > {output} 
        '''

rule Cas10_HD_hmm_maker:
    '''
    1. Extracts the first 10-35 AA of Cas10s that are marked as having an HD domain.
    2. Aligns the Cas10s using Muscle
    3. Creates a new HD-HMM profile from this alignment

    Also has the option to use a premade alignment. Use --use_existing_alignment to do this (does not need argument)
    '''
    input:
        info_table = rules.concatenate_type_iii_info.output,
        cas10_sequences = rules.Cas10_concatenate.output,
    output:
        hmm = base_path + "/Cas10_HD_hmm_profiles/HD_HMM.hmm",
        msa = base_path + "/Cas10_HD_hmm_profiles/HD_HMM.msa",
        faa = base_path + "/Cas10_HD_hmm_profiles/HD_HMM.faa",
    conda: "envs/hmmer.yaml"
    threads: thread_ultrasmall
    shell:
        '''
        python scripts/Cas10_HD_hmm_maker.py --input_table {input.info_table} --cas10_sequences {input.cas10_sequences} --output {output.hmm} --msa {output.msa} --faa {output.faa}  --existing_alignment_path {modified_cas10_hd}
        touch {output.msa}
        touch {output.faa}
        '''

rule Cas10_HD_hmmer:
    '''
    Using the HD-HMM profiles made in rule Cas10_HD_hmm_maker, searches for HD domains in all Cas10s
    '''
    input:
        hmm_db = rules.Cas10_HD_hmm_maker.output.hmm,
        cas10_sequences = rules.Cas10_concatenate.output
    output:
        raw_table = base_path + "/Cas10_HD_hmm_profiles/HD_HMM_hits.out"
    params:
        out = base_path + "/Cas10_HD_hmm_profiles/HD_HMM_hits_temp.out",
        rows = base_path + "/Cas10_HD_hmm_profiles/HD_HMM_hits_rows.out",
        rows1 = base_path + "/Cas10_HD_hmm_profiles/HD_HMM_hits_rows1.out",
        E = "1e-01"
    conda: "envs/hmmer.yaml"
    threads: thread_hogger
    shell:
        '''
        hmmscan --tblout {params.out} --cpu {threads} -E {params.E} {input.hmm_db} {input.cas10_sequences} &> /dev/null
        grep -v "#" {params.out} > {params.rows}||:
        echo "target_name accession query_name accession E-value score bias E-value_best_domain score bias exp reg clu over env dom rep inc description_of_target" > {output.raw_table}
        cat {params.rows} >> {output.raw_table}
        '''

rule merge_HD_hmm_with_info:
    '''
    Merges the HD-HMM hits with the type III info table
    '''
    input:
        info_table = rules.concatenate_type_iii_info.output,
        HD_hmm_hits = rules.Cas10_HD_hmmer.output.raw_table
    output:
        merged_table = base_path + "/Cas10_HD_hmm_profiles/HD_HMM_hits_merged.tsv"
    threads: thread_ultrasmall
    run:
        import pandas as pd
        info = pd.read_csv(str(input.info_table), sep = "\t")
        print(info)
        hmm_hits = pd.read_csv(str(input.HD_hmm_hits), sep = '\s+')
        print(hmm_hits)

        #reduce the hmm_hits to only query_name and E-value_best_domain
        hmm_hits = hmm_hits[["query_name", "E-value_best_domain"]]

        #rename the columns
        hmm_hits.columns = ["Locus", "HD_E-value"]
        #add column "HD_hmm_boolean" and set it to True
        hmm_hits["HD_hmm_boolean"] = True

        #merge the info table with the hmm_hits table
        merged = pd.merge(info, hmm_hits, on = "Locus", how = "left")
        #if some rows were left without a hit, set the HD_hmm_boolean to False
        merged["HD_hmm_boolean"] = merged["HD_hmm_boolean"].fillna(False)

        #from the final file, remove all columns but Locus, HD_E-value and HD_hmm_boolean
        merged = merged[["Locus", "HD_E-value", "HD_hmm_boolean"]]

        #save the merged table
        merged.to_csv(str(output.merged_table), sep = "\t", index = False)


rule R_HD:
    '''
    R script that produces a histogram of Cas10s and their HD domains
    '''
    input: rules.concatenate_type_iii_info.output
    output:
        HD_histogram = base_path + "/4_R_HD/cas10_HD_lengths.png"
    log:
        out = base_path + "/4_R_HD/logs/HD_hist.out",
        err = base_path + "/4_R_HD/logs/HD_hist.err"
    params:
        outputfolder = base_path + "/4_R_HD"
    conda:
        "envs/R.yaml"
    threads: thread_ultrasmall
    shell:
        '''
        Rscript R/HD_R.R --input {input} --output {params.outputfolder} 2> {log.out} 1> {log.err}
        '''


rule Cas10_GGDD_hmm_maker:
    '''
    1. Extracts the X AA of sequence around GGDDs in Cas10s.
    2. Aligns the Cas10s using Muscle
    3. Creates a new HD-HMM profile from this alignment

    As of 23.8.2023, running the script two times: once for GGDD and once for GGDE.
    '''
    input:
        info_table = rules.concatenate_type_iii_info.output,
        cas10_sequences = rules.Cas10_concatenate.output,
    output:
        hmm = base_path + "/Cas10_GGDD_hmm_profiles/GGDD_HMM.hmm",
        msa = base_path + "/Cas10_GGDD_hmm_profiles/GGDD_HMM.msa",
        faa = base_path + "/Cas10_GGDD_hmm_profiles/GGDD_HMM.faa",
        hmm_GGDE = base_path + "/Cas10_GGDD_hmm_profiles/GGDE_HMM.hmm",
        msa_GGDE = base_path + "/Cas10_GGDD_hmm_profiles/GGDE_HMM.msa",
        faa_GGDE = base_path + "/Cas10_GGDD_hmm_profiles/GGDE_HMM.faa",
    conda: "envs/hmmer.yaml"
    threads: thread_hogger
    shell:
        '''
        python scripts/Cas10_GGDD_hmm_maker.py --input_table {input.info_table} --cas10_sequences {input.cas10_sequences} --output {output.hmm} --msa {output.msa} --faa {output.faa} --motif "GGDD"
        python scripts/Cas10_GGDD_hmm_maker.py --input_table {input.info_table} --cas10_sequences {input.cas10_sequences} --output {output.hmm_GGDE} --msa {output.msa_GGDE} --faa {output.faa_GGDE} --motif "GGDE"
        touch {output.msa}
        touch {output.faa}
        '''

rule Cas10_GGDD_hmmer:
    '''
    Using the GGDD-HMM profiles made in rule Cas10_GGDD_hmm_maker, searches for GGDD domains in all Cas10s
    '''
    input:
        hmm_db_GGDD = rules.Cas10_GGDD_hmm_maker.output.hmm,
        hmm_db_GGDE = rules.Cas10_GGDD_hmm_maker.output.hmm_GGDE,
        cas10_sequences = rules.Cas10_concatenate.output
    output:
        raw_table_GGDD = base_path + "/Cas10_GGDD_hmm_profiles/GGDD_HMM_hits.out",
        raw_table_GGDE = base_path + "/Cas10_GGDD_hmm_profiles/GGDE_HMM_hits.out"
    params:
        out = base_path + "/Cas10_GGDD_hmm_profiles/GGDD_HMM_hits_temp.out",
        rows = base_path + "/Cas10_GGDD_hmm_profiles/GGDD_HMM_hits_rows.out",
        rows1 = base_path + "/Cas10_GGDD_hmm_profiles/GGDD_HMM_hits_rows1.out",
        E = "1e-03"
    conda: "envs/hmmer.yaml"
    threads: thread_hogger
    shell:
        '''
        hmmscan --tblout {params.out} --cpu {threads} -E {params.E} {input.hmm_db_GGDD} {input.cas10_sequences}
        grep -v "#" {params.out} > {params.rows}||:
        echo "target_name accession query_name accession E-value score bias E-value_best_domain score bias exp reg clu over env dom rep inc description_of_target" > {output.raw_table_GGDD}
        cat {params.rows} >> {output.raw_table_GGDD}

        hmmscan --tblout {params.out} --cpu {threads} -E {params.E} {input.hmm_db_GGDE} {input.cas10_sequences}
        grep -v "#" {params.out} > {params.rows}||:
        echo "target_name accession query_name accession E-value score bias E-value_best_domain score bias exp reg clu over env dom rep inc description_of_target" > {output.raw_table_GGDE}
        cat {params.rows} >> {output.raw_table_GGDE}
        '''

rule merge_GGDD_hmm_with_info:
    '''
    Merges the GGDD-HMM hits with the type III info table
    '''
    input:
        info_table = rules.concatenate_type_iii_info.output,
        GGDD_hmm_hits = rules.Cas10_GGDD_hmmer.output.raw_table_GGDD,
        GGDE_hmm_hits = rules.Cas10_GGDD_hmmer.output.raw_table_GGDE
    output:
        merged_table = base_path + "/Cas10_GGDD_hmm_profiles/GGDD_HMM_hits_merged.tsv"
    threads: thread_small
    run:
        import pandas as pd
        info = pd.read_csv(str(input.info_table), sep = "\t")
        hmm_hits_GGDD = pd.read_csv(str(input.GGDD_hmm_hits), sep = '\s+')
        hmm_hits_GGDE = pd.read_csv(str(input.GGDE_hmm_hits), sep = '\s+')

        #reduce the hmm_hits to only query_name and E-value_best_domain
        hmm_hits_GGDD = hmm_hits_GGDD[["query_name", "E-value_best_domain"]]
        hmm_hits_GGDE = hmm_hits_GGDE[["query_name", "E-value_best_domain"]]

        #rename the columns
        hmm_hits_GGDD.columns = ["Locus", "GGDD_E-value"]
        hmm_hits_GGDE.columns = ["Locus", "GGDE_E-value"]
        #add column "GGDD_hmm_boolean" and set it to True
        hmm_hits_GGDD["GGDD_hmm_boolean"] = True
        hmm_hits_GGDE["GGDE_hmm_boolean"] = True

        #merge the info table with the hmm_hits_GGDD table
        merged = pd.merge(info, hmm_hits_GGDD, on = "Locus", how = "left")
        merged = pd.merge(merged, hmm_hits_GGDE, on = "Locus", how = "left")

        #if some rows were left without a hit, set the GGDD_hmm_boolean to False
        merged["GGDD_hmm_boolean"] = merged["GGDD_hmm_boolean"].fillna(False)
        merged["GGDE_hmm_boolean"] = merged["GGDE_hmm_boolean"].fillna(False)

        #create new column cyclase which is true only if column GGDD_hmm_boolean or column GGDE_hmm_boolean is true and cyclase_literal is true
        merged["cyclase"] = merged["cyclase_literal"] & (merged["GGDD_hmm_boolean"] | merged["GGDE_hmm_boolean"])

        #from the final file, remove all columns but Locus, GGDD_E-value and GGDD_hmm_boolean
        merged = merged[["Locus", "GGDD_E-value", "GGDD_hmm_boolean", "GGDE_E-value", "GGDE_hmm_boolean", "cyclase", "cyclase_literal"]]

        #save the merged table
        merged.to_csv(str(output.merged_table), sep = "\t", index = False)


rule combine_GGDD_HMM_to_mastertable:
    '''
    Merges the HD and GGDD-HMM data with the master table. Also merges domains info.
    Note that this is not the final mastertable (see rule mastercombiner)
    '''
    input:
        info = rules.concatenate_type_iii_info.output,
        HD_hmm_hits = rules.merge_HD_hmm_with_info.output.merged_table,
        GGDD_hmm_hits = rules.merge_GGDD_hmm_with_info.output.merged_table,
        domains = rules.annotate_bacteria_and_archaea_domains.output,
    output:
        final_info_table = base_path + "/mastertable.tsv"
    threads: thread_ultrasmall
    run:
        import pandas as pd
        info = pd.read_csv(str(input.info), sep = "\t")
        domains = pd.read_csv(str(input.domains), sep = "\t")
        hmm_hits_HD = pd.read_csv(str(input.HD_hmm_hits), sep = "\t")
        hmm_hits_GGDD = pd.read_csv(str(input.GGDD_hmm_hits), sep = "\t")
        print(hmm_hits_GGDD)
        merged = pd.merge(info, hmm_hits_HD, on = "Locus", how = "left")
        merged = pd.merge(merged, hmm_hits_GGDD, on = "Locus", how = "left")
        merged = pd.merge(merged, domains, left_on = "Sample", right_on = "sample", how = "left")
        merged.to_csv(str(output.final_info_table), sep = "\t", index = False)


#a rule that uses the output file cora_type_iii_info from the rule above filter samples with CorA and further divide them into CRISPR-Cas subtypes
rule cora_plot_extractor:
    '''
    Gets cctyper generated plots for each sample with a CorA and a type III CRISPR-Cas system.
    '''
    input: 
        info = rules.concatenate_type_iii_info.output[0]
    output:
        done = base_path + "/xtra1_cora_iii_loci_plots/done.done",
    threads: thread_small
    run:
        import pandas as pd
        import shutil
        #using pandas, get the info from the cora_type_iii_info file and filter out samples that do not have a CorA
        cora_type_iii_info = pd.read_csv(input.info, sep = "\t")
        cora_type_iii_info = cora_type_iii_info[cora_type_iii_info["CorA"] == True]

        #create output folder if it does not exist
        if not os.path.exists(base_path + "/xtra1_cora_iii_loci_plots"):
            os.makedirs(base_path + "/xtra1_cora_iii_loci_plots")

        #For each sample in the pandas dataframe, extract the corresponding plot from the cctyper folder (path is base_path + "/07_cctyper/{j}/plot.png) and copy it to the output folder in a CRISPR-Cas subtype subfolder (e.g. III-A or III-B) depending on the Subtype column in the cora_type_iii_info file
        for index, row in cora_type_iii_info.iterrows():
            sample = row["Sample"]
            subtype = row["Subtype"]
            print(sample + "," + subtype)
            #check if subtype does not contain the substring "Hybrid" (this is because some samples have a hybrid subtype, e.g. III-A/B, and we don't want to create a III-A/B folder)
            if "Hybrid" not in subtype:
                #if the subtype folder does not exist, create it
                if not os.path.exists(base_path + "/xtra1_cora_iii_loci_plots/" + subtype):
                    os.makedirs(base_path + "/xtra1_cora_iii_loci_plots/" + subtype)
                shutil.copyfile(base_path + "/07_cctyper/" + sample + "/plot.png", base_path + "/xtra1_cora_iii_loci_plots/" + subtype + "/" + sample + "_plot.png")

        #create the done.done file to indicate that the rule has finished running
        open(output.done, "w").close()


rule CorA_concatenate:
    '''
    This version works on the cATyper outputs instead of CCTyper outputs
    '''
    input:
        coras = aggregate_CorA_sequences_catyper,
        type_iii_wildcarder_finished = rules.concatenate_type_iii_info.output
    output: base_path + "/09_crispr_iii_CorA/CorAs.faa"
    threads: thread_ultrasmall
    shell:
        '''
        cat {input.coras} > {output}
        '''

rule CorA_cluster:
    '''

    '''
    input: rules.CorA_concatenate.output
    output:
        proteins = base_path + "/15_CorA_cluster/CorA_cluster.faa",
        clusterinfo = base_path + "/15_CorA_cluster/CorA_cluster.faa.clstr",
    conda: "envs/trees.yaml"
    log:
        out = base_path + "/15_CorA_cluster/logs/CorA_align.out",
        err = base_path + "/15_CorA_cluster/logs/CorA_align.err"
    threads: thread_hogger
    shell:
        '''
        cd-hit -i {input} -o {output.proteins} -c 0.90 -n 5 -d 0 -M 16000 -T {threads}
        '''


rule CorA_align:
    input: rules.CorA_cluster.output.proteins
    output: base_path + "/16_CorA_align/CorA_alignment.afa",
    conda: "envs/trees.yaml"
    log:
        out = base_path + "/16_CorA_align/logs/CorA_align.out",
        err = base_path + "/16_CorA_align/logs/CorA_align.err"
    threads: thread_hogger
    shell:
        '''
        muscle -super5 "{input}" -output "{output}" -threads {threads} 2> {log.err} 1> {log.out}
        '''

rule CorA_align_unclustered:
    input: rules.CorA_concatenate.output
    output: base_path + "/16_CorA_align_unclustered/CorA_alignment_unclustered.afa",
    conda: "envs/trees.yaml"
    log:
        out = base_path + "/16_CorA_align/logs/CorA_align.out",
        err = base_path + "/16_CorA_align/logs/CorA_align.err"
    threads: thread_hogger
    shell:
        '''
        muscle -super5 "{input}" -output "{output}" -threads {threads} 2> {log.err} 1> {log.out}
        '''

rule CorA_tree:
    '''

    '''
    input: rules.CorA_align.output
    output: base_path + "/17_CorA_tree/CorA_tree.txt",
    conda: "envs/trees.yaml"
    threads: thread_hogger
    shell:
        '''
        FastTree -wag -gamma {input} > {output}
        '''

rule CorA_tree_unclustered:
    '''

    '''
    input: rules.CorA_align_unclustered.output
    output: base_path + "/17_CorA_tree_unclustered/CorA_tree_unclustered.txt",
    conda: "envs/trees.yaml"
    threads: thread_hogger
    shell:
        '''
        FastTree -wag -gamma {input} > {output}
        '''


rule Cas10_cluster:
    '''
    '''
    input: 
        cas10 = rules.Cas10_concatenate.output
    output:
        proteins = base_path + "/10_Cas10_cluster/cas10_cluster.faa",
        clusterinfo = base_path + "/10_Cas10_cluster/cas10_cluster.faa.clstr",
    conda: "envs/trees.yaml"
    log:
        out = base_path + "/10_Cas10_cluster/logs/cas10_align.out",
        err = base_path + "/10_Cas10_cluster/logs/cas10_align.err"
    threads: thread_hogger
    shell:
        '''
        echo {protein_clustering}
        if [ {protein_clustering} = "True" ]; then
            cd-hit -i {input.cas10} -o {output.proteins} -c 0.99 -n 5 -d 0 -M 16000 -T {threads}
        else
            cp {input.cas10} {output.proteins}
            touch {output.clusterinfo}
        fi
        '''


rule Cas10_align:
    '''
    
    '''
    input: rules.Cas10_cluster.output.proteins
    output: base_path + "/11_Cas10_align/cas10_alignment.afa",
    conda: "envs/trees.yaml"
    log:
        out = base_path + "/11_Cas10_align/logs/cas10_align.out",
        err = base_path + "/11_Cas10_align/logs/cas10_align.err"
    threads: thread_hogger
    shell:
        '''
        muscle -super5 "{input}" -output "{output}" -threads {threads} 2> {log.err} 1> {log.out}
        '''

rule Cas10_tree:
    '''

    '''
    input: rules.Cas10_align.output
    output: base_path + "/12_Cas10_tree/cas10_tree.txt",
    conda: "envs/trees.yaml"
    threads: thread_hogger
    shell:
        '''
        FastTree -wag -gamma {input} > {output}
        '''

rule validated_new_effectors:
    '''
    This rule is updated as potential new effectors from group4 emerge.
    Creates tables for each locus showing whether a validated effector is
    found in that locus. Does not interfere with running the group4 analysis
    and the validated effectors found here will also be listed in group4 proteins,
    as that is where they were originally found.
    '''
    input:
        proteins = rules.cATyper_hmm_search.output.contig_proteins
    output:
        #validated_new_effectors = base_path + "/60_validated_new_effectors/{c}/{c}_validated_new_effectors_orig.tsv",
        temp_rows = base_path + "/60_validated_new_effectors/{c}/{c}_validated_new_effectors_temp.tsv",
        hmm_rows = base_path + "/60_validated_new_effectors/{c}/{c}_validated_new_effectors_hmm.tsv"
    conda: "envs/hmmer.yaml"
    params:
        validated_effectors_hmm_db = validated_effectors_hmm_db,
        temp_hmm = base_path + "/60_validated_new_effectors/{c}/{c}_validated_new_effectors.temp",
        evalue = "1e-10"
    log:
        out = base_path + "/60_validated_new_effectors/logs/{c}/{c}_validated_new_effectors.out",
        err = base_path + "/60_validated_new_effectors/logs/{c}/{c}_validated_new_effectors.err"
    threads: thread_ultrasmall
    shell:
        '''
        hmmscan --domtblout {params.temp_hmm} --cpu {threads} -E {params.evalue} {params.validated_effectors_hmm_db} {input.proteins} &> /dev/null
        echo "Removing commented rows" >> {log.out}
        grep -v "^#" {params.temp_hmm} > {output.temp_rows} ||:
        echo "Writing header" >> {log.out}
        echo -e 'target_name\taccession\ttlen\tquery_name\taccession\tqlen\tE-value_fullseq\tscore_fullseq\tbias_fullseq\t#_domain\tof_domain\tc-Evalue_domain\ti-Evalue_domain\tscore_domain\tbias_domain\t_hmm_from\thmm_to\t_ali_from\tali_to\tenv_from\tenv_to\tacc\tdescription' > {output.hmm_rows}
        echo "Checking if hits were found" >> {log.out}
        if [ -s {output.temp_rows} ]; then 
            echo "Hits found for {wildcards.c}" >> {log.out}
            cat {output.temp_rows} >> {output.hmm_rows}
        else
            echo "No hits found for {wildcards.c}" >> {log.out}
            touch {output.hmm_rows}
        fi
        '''

rule concatenate_validate_new_effectors_hmm:
    input: aggregate_validated_new_effectors_hmm
    output: base_path + "/60_validated_new_effectors/validated_new_effectors_all_hmm.tsv"
    threads: thread_ultrasmall
    shell:
        """
        awk '(NR == 1) || (FNR > 1)' {input} > {output}
        """

rule validated_new_effectors_analysis:
    input:
        crispr_positive_samples = base_path + "/071_cctyper_loci/{c}/cas_operons.tsv",
        contig_proteins = base_path + "/10_cATyper_hmm/{c}/{c}_contig_proteins.faa",
        hmm_rows = base_path + "/60_validated_new_effectors/{c}/{c}_validated_new_effectors_hmm.tsv",
    output:
        catyper = base_path + "/61_validated_new_effectors_analysis/{c}/{c}_cATyper_results.tsv",
        hmm_targets = base_path + "/61_validated_new_effectors_analysis/{c}/{c}_cATyper_hmm_targets.tsv",
        effector_to_protein = base_path + "/61_validated_new_effectors_analysis/{c}/{c}_effector_to_protein.tsv",
        protein_to_effector = base_path + "/61_validated_new_effectors_analysis/{c}/{c}_protein_to_effector.tsv",
        plottable_effector_positions = base_path + "/61_validated_new_effectors_analysis/{c}/{c}_plottable_effector_positions.tsv",
    params:
        outdir = base_path + "/61_validated_new_effectors_analysis/{c}",
        host_genomes_folder = base_path + "/06_host_genomes",
        hmm_msa_folder = validated_effectors_folder + "/profiles",
        tmhmm_model_path = TM_path + "/TMHMM2.0.model",
    conda: "envs/groupChar.yaml"
    log:
        out = base_path + "/61_validated_new_effectors_analysis/logs/{c}.out",
        err = base_path + "/61_validated_new_effectors_analysis/logs/{c}.err",
    threads: thread_ultrasmall
    shell:
        '''
        echo "Running validated effector analysis" >> {log.out}
        echo "Listing targets" >> {log.out}
        ls {params.hmm_msa_folder}/*/*.hmm > {output.hmm_targets}
        python scripts/catyper_prepper_10.py --locus {wildcards.c} --cas_operons_file {input.crispr_positive_samples} --output_folder {params.outdir} --host_genomes_folder {params.host_genomes_folder} --mode post_hmm --hmm_targets {output.hmm_targets} --hmm_rows {input.hmm_rows} --catyper_out {output.catyper} --catyper_type "validated_new_effector" --effector_plot_data {output.plottable_effector_positions} --tmhmm_model_path {params.tmhmm_model_path} 2> {log.err} 1> {log.out}
        '''

rule concatenate_validated_new_effectors_analysis:
    input: aggregate_validated_new_effectors_analysis
    output: base_path + "/61_validated_new_effectors_analysis/validated_new_effectors_all_hmm.tsv"
    threads: thread_ultrasmall
    shell:
        """
        awk '(NR == 1) || (FNR > 1)' {input} > {output}
        """


rule concatenate_validated_new_effectors_scores:
    '''
    Concatenates info on the effector hits from cATyper analysis
    '''
    input:
        protein_to_effector = aggregate_validated_new_effectors_pte,
        effector_to_protein = aggregate_validated_new_effectors_etp
    output:
        protein_to_effector_concatenated = base_path + "/61_validated_new_effectors_analysis/val_new_effectors_protein_to_effector.tsv",
        effector_to_protein_concatenated = base_path + "/61_validated_new_effectors_analysis/val_new_effectors_to_protein.tsv",
    threads: thread_ultrasmall
    shell:
        '''
        cat {input.protein_to_effector} > {output.protein_to_effector_concatenated}
        cat {input.effector_to_protein} > {output.effector_to_protein_concatenated}
        '''

rule analyse_validated_new_effectors_scores:
    '''
    The algorithm that analyses how well each hmm profile performs in validated_new_effectors.
    '''
    input:
        pte = rules.concatenate_validated_new_effectors_scores.output.protein_to_effector_concatenated,
        etp = rules.concatenate_validated_new_effectors_scores.output.effector_to_protein_concatenated,
    output:
        effector_scores = base_path + "/61_validated_new_effectors_analysis/val_new_effectors_scores.tsv",
        effector_scores_summary = base_path + "/61_validated_new_effectors_analysis/val_new_effectors_scores_summary.tsv",
        #effector_scores_plot1 = base_path + "/11_cATyper_analysis/cATyper_effector_scores_plot1.png",
    params:
        outdir = base_path + "/61_validated_new_effectors_analysis"
    threads: thread_ultrasmall
    conda: "envs/analysis.yaml"
    shell:
        '''
        python scripts/catyper_effector_scores.py --pte {input.pte} --etp {input.etp} --output {output.effector_scores} --output_summary {output.effector_scores_summary} --outdir {params.outdir}
        touch {output.effector_scores}
        touch {output.effector_scores_summary}
        '''

rule heatmap_known_validated_effectors:
    input:
        etp_validated = rules.concatenate_validated_new_effectors_scores.output.effector_to_protein_concatenated,
        etp_known = rules.concatenate_cATyper_analysis_effector_scores.output.effector_to_protein_concatenated
    output:
        effectors_combined = base_path + "/70_validated_and_known_effectors_heatmap/etp_combined.tsv",
        effector_scores_summary = base_path + "/70_validated_and_known_effectors_heatmap/effectors.tsv",
    params:
        outdir = base_path + "/70_validated_and_known_effectors_heatmap",
        dummyout = base_path + "/70_validated_and_known_effectors_heatmap/dummyout.tsv"
    threads: thread_ultrasmall
    conda: "envs/analysis.yaml"
    shell:
        '''
        cat {input.etp_validated} {input.etp_known} > {output.effectors_combined}
        python scripts/catyper_effector_scores.py --etp {output.effectors_combined} --output {params.dummyout} --output_summary {output.effector_scores_summary} --outdir {params.outdir}
        '''


#### RING NUCLEASE RULES ####

rule ring_nucleases:
    '''
    Reuses code from validated new effectors to find ring nucleases.
    Runs hmmscan against all loci using the custom ring nuclease database as, well, the database.
    '''
    input:
        proteins = rules.cATyper_hmm_search.output.contig_proteins
    output:
        temp_rows = base_path + "/70_ring_nucleases/{c}/{c}_ring_nucleases_temp.tsv",
        hmm_rows = base_path + "/70_ring_nucleases/{c}/{c}_ring_nucleases_hmm.tsv"
    conda: "envs/hmmer.yaml"
    params:
        rn_db = ring_nuclease_db,
        temp_hmm = base_path + "/70_ring_nucleases/{c}/{c}_ring_nucleases.temp",
        evalue = "1e-8"
    log:
        out = base_path + "/70_ring_nucleases/logs/{c}/{c}_ring_nucleases.out",
        err = base_path + "/70_ring_nucleases/logs/{c}/{c}_ring_nucleases.err"
    threads: thread_ultrasmall
    shell:
        '''
        hmmscan --domtblout {params.temp_hmm} --cpu {threads} -E {params.evalue} {params.rn_db} {input.proteins} &> /dev/null
        echo "Removing commented rows" >> {log.out}
        grep -v "^#" {params.temp_hmm} > {output.temp_rows} ||:
        echo "Writing header" >> {log.out}
        echo -e 'target_name\taccession\ttlen\tquery_name\taccession\tqlen\tE-value_fullseq\tscore_fullseq\tbias_fullseq\t#_domain\tof_domain\tc-Evalue_domain\ti-Evalue_domain\tscore_domain\tbias_domain\t_hmm_from\thmm_to\t_ali_from\tali_to\tenv_from\tenv_to\tacc\tdescription' > {output.hmm_rows}
        echo "Checking if hits were found" >> {log.out}
        if [ -s {output.temp_rows} ]; then 
            echo "Hits found for {wildcards.c}" >> {log.out}
            cat {output.temp_rows} >> {output.hmm_rows}
        else
            echo "No hits found for {wildcards.c}" >> {log.out}
            touch {output.hmm_rows}
        fi
        '''

rule concatenate_validate_ring_nucleases_hmm:
    input: aggregate_ring_nucleases_hmm
    output: base_path + "/70_ring_nucleases/ring_nucleases_all_hmm.tsv"
    threads: thread_ultrasmall
    shell:
        """
        awk '(NR == 1) || (FNR > 1)' {input} > {output}
        """

rule ring_nucleases_analysis:
    '''
    Generates locus-specific information on the presence of ring nucleases.
    '''
    input:
        crispr_positive_samples = base_path + "/071_cctyper_loci/{c}/cas_operons.tsv",
        contig_proteins = base_path + "/10_cATyper_hmm/{c}/{c}_contig_proteins.faa",
        hmm_rows = base_path + "/70_ring_nucleases/{c}/{c}_ring_nucleases_hmm.tsv",
    output:
        catyper = base_path + "/71_ring_nucleases_analysis/{c}/{c}_cATyper_results.tsv",
        hmm_targets = base_path + "/71_ring_nucleases_analysis/{c}/{c}_cATyper_hmm_targets.tsv",
        effector_to_protein = base_path + "/71_ring_nucleases_analysis/{c}/{c}_effector_to_protein.tsv",
        protein_to_effector = base_path + "/71_ring_nucleases_analysis/{c}/{c}_protein_to_effector.tsv",
        plottable_effector_positions = base_path + "/71_ring_nucleases_analysis/{c}/{c}_plottable_effector_positions.tsv",
    params:
        outdir = base_path + "/71_ring_nucleases_analysis/{c}",
        host_genomes_folder = base_path + "/06_host_genomes",
        hmm_msa_folder = ring_nuclease_folder + "/profiles",
    conda: "envs/hmmer.yaml"
    log:
        out = base_path + "/71_ring_nucleases_analysis/logs/{c}.out",
        err = base_path + "/71_ring_nucleases_analysis/logs/{c}.err",
    threads: thread_ultrasmall
    shell:
        '''
        echo "Running validated effector analysis" >> {log.out}
        echo "Listing targets" >> {log.out}
        ls {params.hmm_msa_folder}/*/*.hmm > {output.hmm_targets}
        python scripts/catyper_prepper_10.py --locus {wildcards.c} --cas_operons_file {input.crispr_positive_samples} --output_folder {params.outdir} --host_genomes_folder {params.host_genomes_folder} --mode post_hmm --hmm_targets {output.hmm_targets} --hmm_rows {input.hmm_rows} --catyper_out {output.catyper} --catyper_type "ring_nucleases" --effector_plot_data {output.plottable_effector_positions} --ring_nuclease True 2> {log.err} 1> {log.out}
        '''

rule concatenate_ring_nucleases_analysis:
    input: aggregate_ring_nucleases_analysis
    output: base_path + "/71_ring_nucleases_analysis/ring_nucleases_all_hmm.tsv"
    threads: thread_ultrasmall
    shell:
        """
        awk '(NR == 1) || (FNR > 1)' {input} > {output}
        """


rule concatenate_ring_nucleases_scores:
    '''
    Concatenates info on the ring nuclease hits from the ring nuclease analysis
    '''
    input:
        protein_to_effector = aggregate_ring_nucleases_pte,
        effector_to_protein = aggregate_ring_nucleases_etp,
        ring_nuclease_concat = rules.concatenate_ring_nucleases_analysis.output
    output:
        protein_to_effector_concatenated = base_path + "/71_ring_nucleases_analysis/ring_nucleases_protein_to_effector.tsv",
        effector_to_protein_concatenated = base_path + "/71_ring_nucleases_analysis/ring_nucleases_effectors_to_protein.tsv",
    threads: thread_small
    shell:
        '''
        cat {input.protein_to_effector} > {output.protein_to_effector_concatenated}
        cat {input.effector_to_protein} > {output.effector_to_protein_concatenated}
        '''

rule analyse_ring_nucleases_scores:
    '''
    The algorithm that analyses how well each hmm profile performs in ring_nucleases.
    '''
    input:
        pte = rules.concatenate_ring_nucleases_scores.output.protein_to_effector_concatenated,
        etp = rules.concatenate_ring_nucleases_scores.output.effector_to_protein_concatenated,
    output:
        effector_scores = base_path + "/71_ring_nucleases_analysis/val_new_effectors_scores.tsv",
        effector_scores_summary = base_path + "/71_ring_nucleases_analysis/val_new_effectors_scores_summary.tsv",
        #effector_scores_plot1 = base_path + "/11_cATyper_analysis/cATyper_effector_scores_plot1.png",
    params:
        outdir = base_path + "/71_ring_nucleases_analysis"
    threads: thread_small
    conda: "envs/analysis.yaml"
    shell:
        '''
        python scripts/catyper_effector_scores.py --pte {input.pte} --etp {input.etp} --output {output.effector_scores} --output_summary {output.effector_scores_summary} --outdir {params.outdir}
        touch {output.effector_scores}
        touch {output.effector_scores_summary}
        '''

rule heatmap_ring_nucleases:
    input:
        ring_nuclease_etp = rules.concatenate_ring_nucleases_scores.output.effector_to_protein_concatenated
    output:
        #effectors_combined = base_path + "/72_ring_nucleases_heatmap/etp_combined.tsv",
        effector_scores_summary = base_path + "/72_ring_nucleases_heatmap/effectors.tsv",
    params:
        outdir = base_path + "/72_ring_nucleases_heatmap",
        dummyout = base_path + "/72_ring_nucleases_heatmap/dummyout.tsv"
    threads: thread_small
    conda: "envs/hmmer.yaml"
    log:
        out = base_path + "/72_ring_nucleases_heatmap/logs/heatmap.out",
        err = base_path + "/72_ring_nucleases_heatmap/logs/heatmap.err"
    shell:
        '''
        python scripts/catyper_effector_scores.py --etp {input.ring_nuclease_etp} --output {params.dummyout} --output_summary {output.effector_scores_summary} --outdir {params.outdir} 2> {log.err} 1> {log.out}
        '''

rule crn4_aligner:
    '''
    Deprecated. Csx14 lumped with Csx1 and 15,16,20 kept separate.

    Legacy text:
    Crn4 refers to Csx15, 16 and 20 (a, b, c). These are structurally almost identical proteins,
    but sequence-wise quite diverged. Here we align them by aa sequence
    and create tree in subsequent rule.

    We concatenate the sequences of the three proteins and align them with muscle.

    Crn1 is used as outgroup. Csx14 added for now too (21.5.2024)
    '''
    input: rules.concatenate_ring_nucleases_analysis.output #this is just a trigger to start rule
    output: base_path + "/74_crn4_align/crn4_alignment.afa",
    conda: "envs/trees.yaml"
    log:
        out = base_path + "/74_crn4_align/logs/crn4_align.out",
        err = base_path + "/74_crn4_align/logs/crn4_align.err"
    threads: thread_hogger
    params:
        concatenated_crn4_seqs = base_path + "/74_crn4_align/crn4_concatenated.faa",
        base_folder = base_path + "/71_ring_nucleases_analysis",
        csx14_renamed = base_path + "/74_crn4_align/csx14_renamed.faa",
        csx14_orig = base_path + "/74_crn4_align/csx14_orig.faa",
        csx15_renamed = base_path + "/74_crn4_align/csx15_renamed.faa",
        csx15_orig = base_path + "/74_crn4_align/csx15_orig.faa",
        csx16_renamed = base_path + "/74_crn4_align/csx16_renamed.faa",
        csx16_orig = base_path + "/74_crn4_align/csx16_orig.faa",
        csx20_renamed = base_path + "/74_crn4_align/csx20_renamed.faa",
        csx20_orig = base_path + "/74_crn4_align/csx20_orig.faa",
        crn1_renamed = base_path + "/74_crn4_align/crn1_renamed.faa",
        crn1_orig = base_path + "/74_crn4_align/crn1_orig.faa"
    shell:
        '''
        cat {params.base_folder}/*/csx20.faa > {params.csx20_orig}
        cat {params.base_folder}/*/csx14.faa > {params.csx14_orig}
        cat {params.base_folder}/*/csx15.faa > {params.csx15_orig}
        cat {params.base_folder}/*/csx16.faa > {params.csx16_orig}
        cat {params.base_folder}/*/crn1.faa > {params.crn1_orig}

        python scripts/crn4_header_renamer.py --input {params.csx20_orig} --output {params.csx20_renamed} --prefix "csx20"
        python scripts/crn4_header_renamer.py --input {params.csx14_orig} --output {params.csx14_renamed} --prefix "csx14"
        python scripts/crn4_header_renamer.py --input {params.csx15_orig} --output {params.csx15_renamed} --prefix "csx15"
        python scripts/crn4_header_renamer.py --input {params.csx16_orig} --output {params.csx16_renamed} --prefix "csx16"
        python scripts/crn4_header_renamer.py --input {params.crn1_orig} --output {params.crn1_renamed} --prefix "crn1"

        cat {params.csx14_renamed} {params.csx15_renamed} {params.csx16_renamed} {params.csx20_renamed} {params.crn1_renamed} > {params.concatenated_crn4_seqs}
        
        muscle -super5 "{params.concatenated_crn4_seqs}" -output "{output}" -threads {threads} 2> {log.err} 1> {log.out}
        '''

rule crn4_tree:
    '''
    '''
    input: rules.crn4_aligner.output
    output:
        tree1 = base_path + "/75_crn4_tree/crn4_tree1.txt",
        tree2 = base_path + "/75_crn4_tree/crn4_tree2.txt",
        tree3 = base_path + "/75_crn4_tree/crn4_tree3.txt",
        tree4 = base_path + "/75_crn4_tree/crn4_tree4.txt"
    conda: "envs/trees.yaml"
    threads: thread_hogger
    shell:
        '''
        FastTree -wag -gamma {input} > {output.tree1}
        FastTree -wag -gamma {input} > {output.tree2}
        FastTree -wag -gamma {input} > {output.tree3}
        FastTree -wag -gamma {input} > {output.tree4}
        '''

rule ring_nucleases_fusions:
    '''
    Looks at the original ring nuclease hits and finds out if a protein has multiple credible ring nuclease annotations
    to indicate a fusion between two ring nucleases.
    '''
    input:
        hmm_rows = base_path + "/70_ring_nucleases/{c}/{c}_ring_nucleases_hmm.tsv",
    output:
        ring_nuclease_fusions = base_path + "/73_ring_nuclease_fusions/{c}/{c}_ring_nuclease_fusions.tsv",
    params:
        outdir = base_path + "/73_ring_nuclease_fusions/{c}",
    conda: "envs/hmmer.yaml"
    log:
        out = base_path + "/73_ring_nuclease_fusions/logs/{c}.out",
        err = base_path + "/73_ring_nuclease_fusions/logs/{c}.err",
    threads: thread_ultrasmall
    shell:
        '''
        python scripts/ringnucleasefusionfinder.py --locus {wildcards.c} --output_folder {params.outdir} --hmm_rows {input.hmm_rows} 2> {log.err} 1> {log.out}
        '''

rule concatenate_ring_nuclease_fusions:
    input: aggregate_ring_nuclease_fusions
    output: base_path + "/73_ring_nuclease_fusions/ring_nuclease_fusions_all_hmm.tsv"
    threads: thread_ultrasmall
    shell:
        """
        awk '(NR == 1) || (FNR > 1)' {input} > {output}
        """


rule ring_fusions:
    '''
    Takes protein_to_effector data from all three search modules (known effectors, new validated effectors and ring nucleases)
    and by comparing the proteins listed in ring nucleases, creates a table
    showing to which effectors the ring nucleases are fused to
    '''
    input:
        known_effectors = rules.concatenate_cATyper_analysis_effector_scores.output.protein_to_effector_concatenated,
        validated_new_effectors = rules.concatenate_validated_new_effectors_scores.output.protein_to_effector_concatenated,
        ring_nucleases = rules.concatenate_ring_nucleases_scores.output.protein_to_effector_concatenated,
    output:
        ring_fusions = base_path + "/100_ring_fusions/ring_fusions.tsv"
    threads: thread_ultrasmall
    log:
        out = base_path + "/100_ring_fusions/logs/ring_fusions.out",
        err = base_path + "/100_ring_fusions/logs/ring_fusions.err"
    shell:
        '''
        python scripts/ring_fusions.py --known_effectors {input.known_effectors} --validated_new_effectors {input.validated_new_effectors} --ring_nucleases {input.ring_nucleases} --output {output.ring_fusions} 2> {log.err} 1> {log.out}
        '''

rule ring_nuclease_cas10_fusions:
    '''
    Runs Cas10 sequences as query against the ring nucleases database to detect any fusions
    '''
    input:
        cas10 = rules.Cas10_cluster.output.proteins
    output:
        temp_rows = base_path + "/101_ring_fusions_cas10/ring_nuclease_cas10_temp.tsv",
        hmm_rows = base_path + "/101_ring_fusions_cas10/ring_nuclease_cas10_hmm.tsv"
    conda: "envs/hmmer.yaml"
    params:
        ring_nuclease_hmm_db = ring_nuclease_db,
        temp_hmm = base_path + "/101_ring_fusions_cas10/ring_nucleases_cas10.temp",
        evalue = "1e-10"
    log:
        out = base_path + "/101_ring_fusions_cas10/logs/ring_fusions_cas10.out",
        err = base_path + "/101_ring_fusions_cas10/logs/ring_fusions_cas10.err"
    threads: thread_hogger
    shell:
        '''        
        hmmscan --domtblout {params.temp_hmm} --cpu {threads} -E {params.evalue} {params.ring_nuclease_hmm_db} {input.cas10}
        echo "Removing commented rows" >> {log.out}
        grep -v "^#" {params.temp_hmm} > {output.temp_rows} ||:
        echo -e 'target_name\taccession\ttlen\tquery_name\taccession\tqlen\tE-value_fullseq\tscore_fullseq\tbias_fullseq\t#_domain\tof_domain\tc-Evalue_domain\ti-Evalue_domain\tscore_domain\tbias_domain\t_hmm_from\thmm_to\t_ali_from\tali_to\tenv_from\tenv_to\tacc\tdescription' > {output.hmm_rows}
        cat {output.temp_rows} >> {output.hmm_rows}
        '''


rule validated_effectors_cas10_fusions:
    '''
    Runs Cas10 sequences as query against the validated effectors database to detect any fusions
    '''
    input:
        cas10 = rules.Cas10_cluster.output.proteins
    output:
        temp_rows = base_path + "/101_validated_effectors_cas10/validated_effectors_cas10_temp.tsv",
        hmm_rows = base_path + "/101_validated_effectors_cas10/validated_effectors_cas10_hmm.tsv"
    conda: "envs/hmmer.yaml"
    params:
        validated_effectors_hmm_db = validated_effectors_hmm_db,
        temp_hmm = base_path + "/101_validated_effectors_cas10/ring_nucleases_cas10.temp",
        evalue = "1e-10"
    log:
        out = base_path + "/101_validated_effectors_cas10/logs/validated_effectors_cas10.out",
        err = base_path + "/101_validated_effectors_cas10/logs/validated_effectors_cas10.err"
    threads: thread_hogger
    shell:
        '''        
        hmmscan --domtblout {params.temp_hmm} --cpu {threads} -E {params.evalue} {params.validated_effectors_hmm_db} {input.cas10}
        echo "Removing commented rows" >> {log.out}
        grep -v "^#" {params.temp_hmm} > {output.temp_rows} ||:
        echo -e 'target_name\taccession\ttlen\tquery_name\taccession\tqlen\tE-value_fullseq\tscore_fullseq\tbias_fullseq\t#_domain\tof_domain\tc-Evalue_domain\ti-Evalue_domain\tscore_domain\tbias_domain\t_hmm_from\thmm_to\t_ali_from\tali_to\tenv_from\tenv_to\tacc\tdescription' > {output.hmm_rows}
        cat {output.temp_rows} >> {output.hmm_rows}
        '''

rule known_effectors_cas10_fusions:
    '''
    Runs Cas10 sequences as query against the known effectors database to detect any fusions
    '''
    input:
        cas10 = rules.Cas10_cluster.output.proteins
    output:
        temp_rows = base_path + "/101_known_effectors_cas10/known_effectors_cas10_temp.tsv",
        hmm_rows = base_path + "/101_known_effectors_cas10/known_effectors_cas10_hmm.tsv"
    conda: "envs/hmmer.yaml"
    params:
        known_effectors_db = hmm_database_folder + "/" + hmm_database_file,
        temp_hmm = base_path + "/101_known_effectors_cas10/ring_nucleases_cas10.temp",
        evalue = "1e-10"
    log:
        out = base_path + "/101_known_effectors_cas10/logs/known_effectors_cas10.out",
        err = base_path + "/101_known_effectors_cas10/logs/known_effectors_cas10.err"
    threads: thread_hogger
    shell:
        '''        
        hmmscan --domtblout {params.temp_hmm} --cpu {threads} -E {params.evalue} {params.known_effectors_db} {input.cas10}
        echo "Removing commented rows" >> {log.out}
        grep -v "^#" {params.temp_hmm} > {output.temp_rows} ||:
        echo -e 'target_name\taccession\ttlen\tquery_name\taccession\tqlen\tE-value_fullseq\tscore_fullseq\tbias_fullseq\t#_domain\tof_domain\tc-Evalue_domain\ti-Evalue_domain\tscore_domain\tbias_domain\t_hmm_from\thmm_to\t_ali_from\tali_to\tenv_from\tenv_to\tacc\tdescription' > {output.hmm_rows}
        cat {output.temp_rows} >> {output.hmm_rows}
        '''

rule mastercombiner:
    '''
    Added 4.9.2023 to incorporate merging and data handling previously done in R in here
    '''
    input:
        concat_taxInfo = rules.concat_taxInfo.output,
        catyper = rules.concatenate_cATyper_analysis.output,
        type_iii_info = rules.combine_GGDD_HMM_to_mastertable.output.final_info_table,
        clustered_unknowns_info = rules.concatenate_unknowns_locus_info.output.info,
        temperature_data = temperature_data,
        validated_new_effectors = rules.concatenate_validated_new_effectors_analysis.output,
        ring_fusions = rules.ring_fusions.output.ring_fusions,
        ring_nucleases = rules.concatenate_ring_nucleases_analysis.output
    output:
        final_info_table = base_path + "/mastertable_v2.tsv",
        group4_IDs = base_path + "/group4/group4_IDs.txt"
    threads: thread_small
    run:
        import pandas as pd
        mastertable = pd.read_csv(str(input.type_iii_info), sep = "\t")
        taxInfo = pd.read_csv(str(input.concat_taxInfo), sep = "\t")
        catyper = pd.read_csv(str(input.catyper), sep = "\t")
        clustered_unknowns_info = pd.read_csv(str(input.clustered_unknowns_info), sep = "\t")
        validated_new_effectors_df = pd.read_csv(str(input.validated_new_effectors), sep = "\t")
        ring_nucleases_df = pd.read_csv(str(input.ring_nucleases), sep = "\t")

        #### RING FUSION PREP ####
        ring_fusions = pd.read_csv(str(input.ring_fusions), sep = "\t")
        #in ring fusions only retain columns locus, ring_nuclease, fusion_components, fusion_protein
        ring_fusions = ring_fusions[["locus", "ring_nuclease", "fusion_components", "fusion_protein"]]

        #### KNOWN EFFECTORS PREP ####
        #in catyper, if any of the following columns are True, set has_known_effector column to True: ca3	ca4	ca5	ca6	sam-amp
        catyper["has_known_effector"] = False
        catyper.loc[catyper[["ca3", "ca4", "ca5", "ca6", "sam-amp"]].any(axis = 1), "has_known_effector"] = True

        #### VALIDATED NEW EFFECTOR PREP ####
#       validated_new_effectors_df = validated_new_effectors_df[["can3", "tirsaved", "cam2", "mem-01", "mem2", "locus"]] #drop cOA columns
        validated_new_effectors_df["has_validated_new_effector"] = False
        #if any column except locus is True, set has_validated_new_effector to True
        validated_new_effectors_df.loc[validated_new_effectors_df[["can3", "tirsaved", "cam3", "cam2"]].any(axis = 1), "has_validated_new_effector"] = True

        #rename column no_effectors to no_validated_new_effectors
        validated_new_effectors_df = validated_new_effectors_df.rename(columns = {"no_effectors": "no_validated_new_effectors"})

        #rename columns ca3	ca4	ca5	ca6 to NE_ca3 NE_ca4 NE_ca5 NE_ca6
        validated_new_effectors_df = validated_new_effectors_df.rename(columns = {"ca3": "NE_ca3", "ca4": "NE_ca4", "ca5": "NE_ca5", "ca6": "NE_ca6", "sam-amp": "NE_sam-amp"})

        #### RING NUCLEASE PREP ####
        ring_nucleases_df["has_ring_nuclease"] = False
        #if any of the following are true, set has_ring_nuclease to True
        ring_nucleases_df.loc[ring_nucleases_df[["ae1-1", "ae1-2", "crn1", "crn2", "crn3", "csx15", "csx16", "csx20", "unk01", "phrogRN", "proteaseRN", "solosavedRN"]].any(axis = 1), "has_ring_nuclease"] = True

        #drop ca4, ca4, ca6 sam-amp and unk from ring_nucleases
        ring_nucleases_df = ring_nucleases_df.drop(columns = ["ca3", "ca4", "ca5", "ca6", "sam-amp", "unk", "mem", "no_effectors"])

        #### MERGE ####
        mastertable = pd.merge(mastertable, taxInfo, on = "Sample", how = "left")
        mastertable = pd.merge(mastertable, catyper, left_on = "Locus", right_on = "locus", how = "left")
        mastertable = pd.merge(mastertable, clustered_unknowns_info, left_on = "Locus", right_on = "locus_id", how = "left")
        mastertable = pd.merge(mastertable, validated_new_effectors_df, left_on = "Locus", right_on = "locus", how = "left", suffixes = ("", "new_eff"))
        mastertable = pd.merge(mastertable, ring_nucleases_df, left_on = "Locus", right_on = "locus", how = "left", suffixes = ("", "rn"))
        mastertable = pd.merge(mastertable, ring_fusions, left_on = "Locus", right_on = "locus", how = "left", suffixes = ("", "rn_fusions"))

        #create columns con_ca3, con_ca4, con_ca5, con_ca6, con_sam_amp. Con stands for consensus
        mastertable["con_ca3"] = False
        mastertable["con_ca4"] = False
        mastertable["con_ca5"] = False
        mastertable["con_ca6"] = False
        mastertable["con_sam-amp"] = False
        #print columsn from mastertable
        #print(mastertable.columns)
        #if ca3 or NE_ca3 is True, set con_ca3 to True. Same for all signal molecules
        mastertable.loc[(mastertable["ca3"] == True) | (mastertable["NE_ca3"] == True), "con_ca3"] = True
        mastertable.loc[(mastertable["ca4"] == True) | (mastertable["NE_ca4"] == True), "con_ca4"] = True
        mastertable.loc[(mastertable["ca5"] == True) | (mastertable["NE_ca5"] == True), "con_ca5"] = True
        mastertable.loc[(mastertable["ca6"] == True) | (mastertable["NE_ca6"] == True), "con_ca6"] = True
        mastertable.loc[(mastertable["sam-amp"] == True) |(mastertable["NE_sam-amp"]), "con_sam-amp"] = True

        mastertable["multiple_signals"] = False
        #if more than one of the con columns is True, set multiple_signals to True
        mastertable.loc[(mastertable[["con_ca3", "con_ca4", "con_ca5", "con_ca6", "con_sam-amp"]].sum(axis = 1) > 1), "multiple_signals"] = True

        mastertable["effector_count_known_new_sum"] = 0
        #sum the number of effectors in the known and validated new effectors columns. Effectors to sum are can3, tirsaved, cam2, cam3, cam2, saved-chat, nucc, calpl, cami1, can1-2, csx1, csx23 csm6-ca6, cora
        mastertable["effector_count_known_new_sum"] = mastertable[["can3", "tirsaved", "cam3", "cam2", "saved-chat", "nucc", "calpl", "cami1", "can1-2", "csx1", "csx23", "csm6-ca6", "cora", "cam1"]].sum(axis = 1)

        #for mastertable rows in which fusion_protein is NaN, make it False
        mastertable["fusion_protein"] = mastertable["fusion_protein"].fillna(False)

        #open the temperature data file
        temperature_data = pd.read_csv(str(input.temperature_data), sep = ",")
        #group by column "genus" and calculate mean temperature for each genus
        mean_temperatures = temperature_data.groupby("genus")["Topt_ave"].mean().dropna()
        print(mean_temperatures)
        #merge mean_temperature with mastertable
        mastertable = pd.merge(mastertable, mean_temperatures, on = "genus", how = "left")
        #rename topt_ave to "temperature"
        mastertable = mastertable.rename(columns = {"Topt_ave": "mean_temperature"})

        #Tidy the mastertable
        #drop the following columns (list)
        droppable_columns = ["sample_x", "Unnamed: 0", "val_x", "locus_x", "sample_y", "val_y", "locus_y", "Unnamed: 0_y", "mem_y", "Unnamed: 0_x"]

        #drop the columns in a for loop with try/catch for each column in case the mastertable changes in the future
        for column in droppable_columns:
            try:
                mastertable = mastertable.drop(columns = [column])
            except:
                pass

        #create a new column multilocus_sample. This a boolean value that reports whether the same sample has other loci
        mastertable["multilocus_sample"] = False

        #create a new column multisignal_sample. This is a boolean value that reports whether the same sample has multiple signals (ca3, ca4, ca6 or sam-amp)
        mastertable["multisignal_sample"] = False

        #create a new column all_signals_sample. This is a string that contains comma separated list of all signals present in the sample
        mastertable["all_signals_sample"] = ""

        #for each sample, go through loci and list which signal molecules are true and update the multisingal_sample and all_signals_sample columns accordingly
        for sample in mastertable["Sample"].unique():
            #get all loci for the sample
            sample_loci = mastertable[mastertable["Sample"] == sample]["Locus"].unique()
            #if the number of loci is greater than 1, set multilocus_sample to True
            if len(sample_loci) > 1:
                mastertable.loc[mastertable["Sample"] == sample, "multilocus_sample"] = True

            #get all signals for the sample
            sample_signals = mastertable[mastertable["Sample"] == sample][["con_ca3", "con_ca4", "con_ca5", "con_ca6", "con_sam-amp"]].sum(axis = 1)
            #if the number of signals is greater than 1, set multisignal_sample to True
            if sample_signals.sum() > 1:
                mastertable.loc[mastertable["Sample"] == sample, "multisignal_sample"] = True

            if sample_signals.sum() > 0:
                mastertable.loc[mastertable["Sample"] == sample, "signal_sample"] = True

            #using columns ca3, ca4, ca6 and sam-amp, create a string of all signals present in the sample by checking if their boolean value is True
            all_signals = ""
            if mastertable.loc[mastertable["Sample"] == sample, "con_ca3"].any() == True:
                all_signals += "ca3,"
            if mastertable.loc[mastertable["Sample"] == sample, "con_ca4"].any() == True:
                all_signals += "ca4,"
            if mastertable.loc[mastertable["Sample"] == sample, "con_ca6"].any() == True:
                all_signals += "ca6,"
            if mastertable.loc[mastertable["Sample"] == sample, "con_sam-amp"].any() == True:
                all_signals += "sam-amp,"
            #remove the last comma from the string
            all_signals = all_signals[:-1]
            #update the all_signals_sample column with the string
            mastertable.loc[mastertable["Sample"] == sample, "all_signals_sample"] = all_signals

        #create column effector_elsewhere_in_sample
        mastertable["effector_elsewhere_in_sample"] = False
        mastertable["effector_elsewhere_in_sample_but_not_here"] = False

        #for each locus, check if the sample has another locus that generates a signal and if so, set cas10_effector_elsewhere_in_sample to True
        for locus in mastertable["Locus"].unique():
            #get the sample for the locus
            sample = mastertable[mastertable["Locus"] == locus]["Sample"].unique()[0]
            #get all loci for the sample
            sample_loci = mastertable[mastertable["Sample"] == sample]["Locus"].unique()
            #if any of the loci in the sample have a signal, set cas10_effector_elsewhere_in_sample to True
            if mastertable[mastertable["Locus"].isin(sample_loci)][["con_ca3", "con_ca4", "con_ca6", "con_sam-amp"]].sum(axis = 1).any() == True:
                mastertable.loc[mastertable["Locus"] == locus, "effector_elsewhere_in_sample"] = True

            #if any of the loci in the sample have a signal but this locus does not have a signal, set cas10_effector_elsewhere_in_sample_but_not_here to True
            #first check if current locus has signal
            if mastertable.loc[mastertable["Locus"] == locus][["con_ca3", "con_ca4", "con_ca6", "con_sam-amp"]].sum(axis = 1).any() == False:
                #if it does not have a signal, check if any other loci in the sample have a signal
                if mastertable[mastertable["Locus"].isin(sample_loci)][["con_ca3", "con_ca4", "con_ca6", "con_sam-amp"]].sum(axis = 1).any() == True:
                    mastertable.loc[mastertable["Locus"] == locus, "effector_elsewhere_in_sample_but_not_here"] = True
            

        #subset the mastertable to include samples where (GGDD_hmm_boolean is True or GGDE_hmm_boolean is True) AND HD_hmm_boolean is False
        mastertable["Cas10_HD_coord"] = pd.to_numeric(mastertable["Cas10_HD_coord"], errors='coerce')
        mastertable_group4 = mastertable[
            ((mastertable["GGDD_hmm_boolean"] == True) | (mastertable["GGDE_hmm_boolean"] == True)) & 
            #(mastertable["HD_hmm_boolean"] == False) & 
            #((mastertable["Cas10_HD_coord"] > 100) | (mastertable["Cas10_HD"] == False)) & 
            (mastertable["unknown_proteins"] == True)
        ]
        mastertable_group4["Locus"].to_csv(str(output.group4_IDs), sep = "\t", index = False)

        #drop any duplicate rows
        mastertable = mastertable.drop_duplicates()
        mastertable.to_csv(str(output.final_info_table), sep = "\t", index = False)


rule node_graph:
    '''
    This rule produces files for Gephi to view effector colocalisation
    using a node graph. It also generates upset plots from the same data.

    For Gephi, it works by first creating a table with all effector combinations (edges) with the headers
    source, target, type and weight. It then looks at each locus in the mastertable
    and upon finding co-occurrence of an effector pair, it adds 1 to the weight column.
    '''
    input:
        mastertable = rules.mastercombiner.output.final_info_table,
        effector_list = rules.heatmap_known_validated_effectors.output.effectors_combined
    output:
        edges =  base_path + "/80_node_graph/edges_all.tsv",
        nodes = base_path + "/80_node_graph/nodes_all.tsv",
    params:
        outdir = base_path + "/80_node_graph",
        edges_basename = base_path + "/80_node_graph/edges.tsv",
        nodes_basename = base_path + "/80_node_graph/nodes.tsv"
    threads: thread_small
    conda: "envs/effector_nodes.yaml"
    shell:
        '''
        python scripts/effector_nodes_and_other_viz.py -i {input.mastertable} -e {input.effector_list} -o {params.outdir} -n {params.nodes_basename} -d {params.edges_basename}
        '''

rule transmembrane_prediction:
    '''
    Using TMHMM predicts transmembrane regions for all proteins in all type III loci
    '''
    input:
        contig_proteins = rules.cATyper_hmm_search.output.contig_proteins
    threads: thread_ultrasmall
    output:
        tmhmm_results = base_path + "/90_transmembrane_prediction/{c}/{c}_tmhmm_results.tsv"
    

rule cOA_RN_explorer:
    '''
    Looks at mastertable_2 output to find loci with effectors associated with given cOA.
    Then examines small proteins in those loci that are potential cOA ring nucleases.
    '''
    input:
        mastertable = rules.mastercombiner.output.final_info_table
    output:
        ca3_loci = base_path + "/91_coa_RN_explorer/ca3_loci.txt",
        ca4_loci = base_path + "/91_coa_RN_explorer/ca4_loci.txt",
        ca6_loci = base_path + "/91_coa_RN_explorer/ca6_loci.txt",
        sam_amp_loci = base_path + "/91_coa_RN_explorer/sam-amp_loci.txt",
    params:
        outbase = base_path + "/91_coa_RN_explorer"
    threads: thread_small
    shell:
        '''
        python scripts/cOA_RN_explorer.py --mastertable {input.mastertable} --output {params.outbase} --coa ca3,ca4,ca6,sam_amp
        '''

rule cOA_small_unk_proteins:
    '''
    Uses the locus list from cOA_RN_explorer and the "unknown" protein from tule unknown_finder.
    It extracts all <150 AA unknown proteins from the loci defined by cOA_RN_explorer
    and writes them to a new .tsv file and .faa file.
    Does this for ca3, ca4, ca6 and sam-amp loci.
    '''
    input:
        cA3_loci = rules.cOA_RN_explorer.output.ca3_loci,
        cA4_loci = rules.cOA_RN_explorer.output.ca4_loci,
        cA6_loci = rules.cOA_RN_explorer.output.ca6_loci,
        sam_amp_loci = rules.cOA_RN_explorer.output.sam_amp_loci,
        all_loci = rules.mastercombiner.output.final_info_table,
        all_unknown_proteins = rules.concatenate_unknowns.output.info
    output:
        cA3_small_unk_proteins_info = base_path + "/92_cOA_small_proteins/ca3_small_unk_proteins.tsv",
        cA3_small_unk_proteins_faa = base_path + "/92_cOA_small_proteins/ca3_small_unk_proteins.faa",
        cA4_small_unk_proteins_info = base_path + "/92_cOA_small_proteins/ca4_small_unk_proteins.tsv",
        cA4_small_unk_proteins_faa = base_path + "/92_cOA_small_proteins/ca4_small_unk_proteins.faa",
        cA6_small_unk_proteins_info = base_path + "/92_cOA_small_proteins/ca6_small_unk_proteins.tsv",
        cA6_small_unk_proteins_faa = base_path + "/92_cOA_small_proteins/ca6_small_unk_proteins.faa",
        sam_amp_small_unk_proteins_info = base_path + "/92_cOA_small_proteins/sam-amp_small_unk_proteins.tsv",
        sam_amp_small_unk_proteins_faa = base_path + "/92_cOA_small_proteins/sam-amp_small_unk_proteins.faa",
        all_small_info = base_path + "/92_cOA_small_proteins/all_small_unk_proteins.tsv",
        all_small_faa = base_path + "/92_cOA_small_proteins/all_small_unk_proteins.faa",
    params:
        outfolder = base_path + "/92_cOA_small_proteins",
        length_cutoff = 150,
        cOAs = "ca3,ca4,ca6,sam_amp"
    threads: thread_ultrasmall
    shell:
        '''
        python scripts/cOA_small_unk_proteins.py --loci "{input.cA3_loci},{input.cA4_loci},{input.cA6_loci},{input.sam_amp_loci}" --all_loci {input.all_loci} --unknown_proteins {input.all_unknown_proteins} --output_folder {params.outfolder} --length_cutoff {params.length_cutoff}
        '''
rule small_unk_blaster:
    '''
    Blasts the small unknown proteins against a local blast database
    '''
    input:
        proteins = rules.cOA_small_unk_proteins.output.all_small_faa
    output:
        blast_results = base_path + "/93_small_unk_blast/small_unk_blast_results.tsv"
    threads: thread_hogger
    params:
        blast_db = local_protein_blast_db
    conda: "envs/blast.yaml"
    shell:
        '''
        blastp -query {input.proteins} -db {params.blast_db} -outfmt 6 -out {output.blast_results} -num_threads {threads}
        '''

rule small_unk_blast_filter:
    '''
    Filters the blast hits by e-value and concatenates the accession numbers into a list from header position 1
    '''
    input: rules.small_unk_blaster.output.blast_results
    output:
        info = base_path + "/93_small_unk_blast/small_unk_blast_filtered.tsv",
        accessions = base_path + "/93_small_unk_blast/small_unk_blast_accessions.txt"
    threads: thread_ultrasmall
    params:
        evalue = 1e-20
    run:
        import pandas as pd
        blast_results = pd.read_csv(str(input[0]), sep = "\t", header = None)
        blast_results = blast_results[blast_results[10] < params.evalue]
        blast_results.to_csv(str(output.info), sep = "\t", index = False, header = False)
        accessions = blast_results[1].unique()
        with open(str(output.accessions), "w") as f:
            for acc in accessions:
                f.write(acc + "\n")

rule fetch_small_unk_targets_ncbi:
    '''
    Uses the accesion list from small_unk_blast_filter to download the sequences from the NCBI db
    and concatenates them into a single .faa file. Removes any proteins longer than cutoff using seqkit program.
    '''
    input:
        accessions = rules.small_unk_blast_filter.output.accessions
    output:
        all_proteins = base_path + "/94_small_unk_ncbi/small_unk_ncbi_proteins.faa",
        under_cutoff = base_path + "/94_small_unk_ncbi/small_unk_ncbi_proteins_under_cutoff.faa",
    params:
        length_cutoff = 200,
        blast_db = local_protein_blast_db,
        all_proteins_precleaned = base_path + "/94_small_unk_ncbi/small_unk_ncbi_proteins_precleaned.faa"
    threads: thread_ultrasmall
    conda: "envs/blast.yaml"
    shell:
        '''
        blastdbcmd -db {params.blast_db} -entry_batch {input.accessions} -out {params.all_proteins_precleaned}
        python scripts/clean_up_blastdbcmd.py --in {params.all_proteins_precleaned} --out {output.all_proteins}
        seqkit seq -M {params.length_cutoff} {output.all_proteins} > {output.under_cutoff}
        '''

rule cluster_small_unk_proteins:
    '''
    Clusters the small unknown protein homologs downloaded from NCBI using CD-HIT
    '''
    input:
        proteins = rules.fetch_small_unk_targets_ncbi.output.under_cutoff
    output:
        clustered_proteins = base_path + "/95_small_unk_clustered/small_unk_clustered.faa",
        cluster_info = base_path + "/95_small_unk_clustered/small_unk_clustered.faa.clstr",
        acc = base_path + "/95_small_unk_clustered/small_unk_ncbi_proteins_under_cutoff_acc.txt", #only accessions for webflags
    threads: thread_hogger
    conda: "envs/groupChar.yaml"
    params:
        cluster_cutoff = 0.4
    shell:
        '''
        cd-hit -i {input.proteins} -o {output.clustered_proteins} -c {params.cluster_cutoff} -n 2 -d 0 -M 16000 -T {threads}
        grep ">" {output.clustered_proteins} | cut -d " " -f 1 | cut -d ">" -f 2 > {output.acc}
        '''

rule webflags_cluster_small_unk_proteins:
    '''
    DEPRECATED and replaced by wildcarded version (rule webflags_cluster_small_unk_proteins_wildcarded)
    Runs webflags on the clustered small unknown proteins.
    
    This rule assumes that the Webflags project has been cloned into the working directory.
    I.e. run git clone https://github.com/GCA-VH-lab/FlaGs2 in the working folder.
   
    When run for first time, it downloads the databases genBank.db and refSeq.db to the running folder.
    These are 1.1 gb and 162 mb respectively.
    It also downloads all genomes (.faa and .gff) that have a match for the query in the temp_genome_files folder. These files are removed after finishing.
    '''
    input:
        acc = rules.cluster_small_unk_proteins.output.acc
    output:
        done = base_path + "/96_small_unk_clustered_webflags/webflags.done"
    threads: thread_small
    params:
        out_base = base_path + "/96_small_unk_clustered_webflags/test_out",
        flanking_gene_count = 10
    conda: "envs/webflags.yaml"
    shell:
        '''
        cd FlaGs2
        python3 FlaGs2.py -p {input.acc} -o {params.out_base} -u {ncbi_email} -vb --gene {params.flanking_gene_count}
        touch {output.done}
        '''

checkpoint small_unk_protein_clustered_wildcarder:
    '''
    Generates wildcards based on the accession in the clustered small unknown proteins
    '''
    input:
        accessions = rules.cluster_small_unk_proteins.output.acc
    output: directory(base_path + "/96_small_unk_clustered_wildcards")
    threads: thread_ultrasmall
    run:
        #create output directory
        import os
        os.makedirs(str(output), exist_ok = True)
        #for every line in the input txt file, create a .txt file in the output directory named after the accession
        with open(str(input.accessions), "r") as f:
            for line in f:
                with open(str(output) + "/" + line.strip() + ".txt", "w") as out:
                    out.write(line.strip())

rule webflags_cluster_small_unk_proteins_wildcarded:
    '''
    Runs webflags on the clustered small unknown proteins.
    
    This rule assumes that the Webflags project has been cloned into the working directory.
    I.e. run git clone https://github.com/GCA-VH-lab/FlaGs2 in the working folder.
   
    When run for first time, it downloads the databases genBank.db and refSeq.db to the running folder.
    These are 1.1 gb and 162 mb respectively.
    It also downloads all genomes (.faa and .gff) that have a match for the query in the temp_genome_files folder. These files are removed after finishing.
    '''
    input:
        protein = base_path + "/96_small_unk_clustered_wildcards/{protein}.txt"
    output:
        done = base_path + "/96_small_unk_clustered_webflags_wildcards/{protein}/webflags.done"
    threads: thread_hogger #using this to prevent concurrent webflags jobs
    params:
        out_base = base_path + "/96_small_unk_clustered_webflags_wildcards/{protein}/{protein}",
        temp_genome_files = base_path + "/96_small_unk_clustered_webflags_wildcards/temp_genome_files",
        flanking_gene_count = 10
    conda: "envs/webflags.yaml"
    shell:
        '''
        mkdir -p {params.temp_genome_files}
        cd FlaGs2
        python3 FlaGs2.py -p {input.protein} -o {params.out_base} -u {ncbi_email} -vb --gene {params.flanking_gene_count}
        touch {output.done}
        '''

rule concatenate_wildcarded_webflags:
    '''
    Uses aggregator function to concatenate wildcarded webflags rule results
    '''
    input: aggregate_webflags_cluster_small_unk_proteins_wildcarded
    output: base_path + "/96_small_unk_clustered_webflags_wildcards/aggregated.done"
    threads: thread_small
    shell:
        '''
        touch {output}
        '''
        

rule small_unk_proteins_clusters_hmmer:
    '''
    Uses hmmscan against pfam to annotate the clustered small unknown protein cluster representatives
    '''
    input:
        multifasta = rules.cluster_small_unk_proteins.output.clustered_proteins
    output:
        raw_table = base_path + "/97_small_unk_hmmer/small_unk_hmmer_raw.tsv"
    params:
        out = base_path + "/97_small_unk_hmmer/small_unk_hmmer_temp.out",
        rows = base_path + "/97_small_unk_hmmer/small_unk_hmmer_rows.out",
        base_path = base_path + "/97_small_unk_hmmer",
        pfam_db = pfam_db,
        E = 1e-3
    threads: thread_hogger
    conda: "envs/hmmer.yaml"
    shell:
        '''
        hmmscan --tblout {params.out} --cpu {threads} -E {params.E} {params.pfam_db} {input.multifasta}
        grep -v "#" {params.out} > {params.rows}||:
        echo "target_name accession query_name accession E-value score bias E-value_best_domain score bias exp reg clu over env dom rep inc description_of_target" > {output.raw_table}
        cat {params.rows} >> {output.raw_table}
        '''


rule align_small_unk_protein_cluster_representatives:
    '''
    Aligns the clustered small unknown protein cluster representatives using muscle
    '''
    input:
        proteins = rules.cluster_small_unk_proteins.output.clustered_proteins
    output:
        aligned_proteins = base_path + "/98_small_unk_align/small_unk_aligned.afa",
    threads: thread_hogger
    conda: "envs/hmmer.yaml"
    log:
        out = base_path + "logs/98_small_unk_align/small_unk_align.out",
        err = base_path + "logs/98_small_unk_align/small_unk_align.err"
    shell:
        '''
        muscle -super5 "{input}" -output "{output}" -threads {threads} 2> {log.err} 1> {log.out} 
        '''

rule tree_small_unk_protein_cluster_representatives:
    '''
    Creates phylogenetic tree of the aligned small protein cluster representatives
    '''
    input: rules.align_small_unk_protein_cluster_representatives.output.aligned_proteins
    output: base_path + "/99_small_unk_tree/small_unk_tree.nwk"
    threads: thread_hogger
    conda: "envs/trees.yaml"
    shell:
        '''
        FastTree -wag -gamma {input} > {output}
        '''

rule align_small_unk_proteins_all:
    '''
    Identical to rule align_small_unk_protein_cluster_representatives but aligns all small unknown proteins,
    not just cluster representatives.
    '''
    input:
        proteins = rules.fetch_small_unk_targets_ncbi.output.under_cutoff
    output:
        aligned_proteins = base_path + "/98_small_unk_align/small_unk_aligned_all.afa",
    threads: thread_hogger
    conda: "envs/hmmer.yaml"
    log:
        out = base_path + "logs/98_small_unk_align/small_unk_align_all.out",
        err = base_path + "logs/98_small_unk_align/small_unk_align_all.err"
    shell:
        '''
        muscle -super5 "{input}" -output "{output}" -threads {threads} 2> {log.err} 1> {log.out} 
        '''

rule tree_small_unk_proteins_all:
    '''
    Identical to rule tree_small_unk_protein_cluster_representatives but creates a tree of all small unknown proteins
    '''
    input: rules.align_small_unk_proteins_all.output.aligned_proteins
    output: base_path + "/99_small_unk_tree/small_unk_tree_all.nwk"
    threads: thread_hogger
    conda: "envs/trees.yaml"
    shell:
        '''
        FastTree -wag -gamma {input} > {output}
        '''
    
rule ae1_family_aligner:
    '''
    Aligns the AE1 family proteins using muscle.
    Uses a manual input list of Ae1 families obtained through extraction from iTol.
    '''
    input:
        proteins = rules.fetch_small_unk_targets_ncbi.output.under_cutoff,
        accessions = "manual_inputs/ae1_family_accessions_all.txt" #manual input stored in the snakemake root directory
    output:
        ae1_family_proteins = base_path + "/100_ae1_family_align/ae1_family_proteins.faa",
        aligned_proteins = base_path + "/100_ae1_family_align/ae1_family_aligned.afa",
    threads: thread_hogger
    conda: "envs/hmmer.yaml"
    log:
        out = base_path + "logs/100_ae1_align/ae1_align.out",
        err = base_path + "logs/100_ae1_align/ae1_align.err"
    shell:
        '''
        python scripts/extract_accessions_from_fasta.py --proteins {input.proteins} --accessions {input.accessions} --output {output.ae1_family_proteins}
        muscle -super5 "{output.ae1_family_proteins}" -output "{output.aligned_proteins}" -threads {threads} 2> {log.err} 1> {log.out} 
        '''




rule ae1_family_tree:
    '''
    Creates a phylogenetic tree of the AE1 family proteins
    '''
    input: rules.ae1_family_aligner.output.aligned_proteins
    output: base_path + "/101_ae1_family_tree/ae1_family_tree.nwk"
    threads: thread_hogger
    conda: "envs/trees.yaml"
    shell:
        '''
        FastTree -wag -gamma {input} > {output}
        '''

rule ae1_family_collector:
    '''
    Collects Ae1 grouped into families based on manual inputs.
    These 7 output .faa files are used when creating the Ae1 hmm profiles
    '''
    input:
        proteins = rules.fetch_small_unk_targets_ncbi.output.under_cutoff,
        accessions = "manual_inputs/ae1_1.txt" #this is just a placeholder.
    output:
        done = base_path + "/102_ae1_family_collector/ae1_family_collector.done"
    threads: thread_ultrasmall
    params:
        base_name = "manual_inputs/ae1_",
        output_folder_basename = base_path + "/102_ae1_family_collector/ae1_"
    shell:
        '''
        #find all base_name files in the manual_inputs directory where base_name is extended by a number and then .txt
        for file in manual_inputs/ae1_[0-7].txt
        do
            #extract the number from the file name
            number=$(echo $file | grep -o -P '(?<=ae1_)[0-9]+')
            #extract the accessions from the proteins using the file
            python scripts/extract_accessions_from_fasta.py --proteins {input.proteins} --accessions $file --output {params.output_folder_basename}$number.faa
        done
        touch {output.done}
        '''

rule potential_other_rn_collector:
    '''
    Similar to the Ae1 collector above, creates fasta files from accessions of other potential
    RN families using the accession list files in the manual_inputs directory.
    The names of the input files are
    rn_phrogs.txt
    rn_proteases.txt
    rn_solosaveds.txt
    '''
    input:
        proteins = rules.fetch_small_unk_targets_ncbi.output.under_cutoff,
        accessions = expand("manual_inputs/rn_{name}.txt", name=["phrogs", "proteases", "solosaveds"])
    output:
        done = base_path + "/103_other_rn_collector/other_rn_collector.done"
    threads: thread_ultrasmall
    params:
        output_folder_basename = base_path + "/103_other_rn_collector/other_rn_"
    shell:
        '''
        for file in {input.accessions}; do
            name=$(echo $file | grep -o -P '(?<=rn_)[a-z]+');
            python scripts/extract_accessions_from_fasta.py --proteins {input.proteins} --accessions $file --output {params.output_folder_basename}$name.faa;
        done;
        touch {output.done}
        '''


rule groupCharacteriser:
    '''
    HMMER based databases (Pfam and CARFSAVED)
    This takes as input a list of locus names (from rule mastercombiner) and performs characterisation of all unknown proteins in those loci.
    Characterisation involves:
        - HMM search against Pfam
        - HMM search against CARF/SAVED
    '''
    input:
        group4_list = rules.mastercombiner.output.group4_IDs,
        unknown_proteins = rules.concatenate_unknowns.output.proteins
    output:
        #group4_CARFSAVED = base_path + "/group4/group4.CARFSAVED.rawresults.tsv",
        group4_pfam = base_path + "/group4/group4.all.pfam.rawresults.tsv",
        clustered_unknowns = base_path + "/group4/group4_clustered_unknowns.faa",
        individual_fastas_created = base_path + "/group4/individual_fastas/fastas_created.done"
    params:
        outdir = base_path + "/group4",
        output_base = base_path + "/group4/group4",
        output_individual_fastas = base_path + "/group4/individual_fastas",
        pfams_temp = base_path + "/group4/group4.all.pfam.rawresults.txt.temp",
        tmhmm_model = TM_path + "/TMHMM2.0.model"
    log:
        out = base_path + "/group4/logs/group4_characteriser.out",
        err = base_path + "/group4/logs/group4_characteriser.err"
    conda:
        "envs/groupChar.yaml"
    threads: thread_small
    shell:
        '''
        python3 scripts/groupCharacteriser.py --input {input.group4_list} --unknown_proteins {input.unknown_proteins} --output_basename {params.output_base} --clustered_unknowns {output.clustered_unknowns} --individual_fastas_dir {params.output_individual_fastas} --tmhmm_model {params.tmhmm_model} --outputfolder {params.outdir} --pfam_db {pfam_db} --carfsaved_db {carfsaved_db}  2> {log.err} 1> {log.out}
        touch {output.individual_fastas_created}
        echo -e "{group4_pfam_headers}" > {output.group4_pfam}
        grep -v "#" {params.pfams_temp} >> {output.group4_pfam}
        '''

rule group4_PDB:
    '''
    This rule runs HHBlits using PDB as database for the group4 unknown proteins
    '''
    input:
        multifasta = rules.groupCharacteriser.output.clustered_unknowns
    output:
        hhsuite_concat = base_path + "/group4/pdb/group4_pdb.tsv",
    params:
        outdir = base_path + "/group4/pdb",
        pdb30 = pdb30_db,
    conda: "envs/hhsuite.yaml"
    threads: thread_hogger
    log:
        out = base_path + "/group4/pdb/logs/pdb.out",
        err = base_path + "/group4/pdb/logs/pdb.err",
    shell:
        '''
        python3 scripts/hhblits_wrapper.py --input {input.multifasta} --output_basepath {params.outdir} --database {params.pdb30}
        cat {params.outdir}/hhblits/*.tsv > {output.hhsuite_concat}
        '''

rule parse_hhsuite_group4_pdb:
    '''
    Parser output from HHSuite (PDB)
    '''
    input: rules.group4_PDB.output.hhsuite_concat
    output: base_path + "/group4/pdb/group4_pdb_parsed.tsv"
    conda: "envs/hhsuite.yaml"
    params:
        database = "PDB"
    threads: thread_ultrasmall
    shell:
        '''
        python3 scripts/hhsuite_parser.py --infile {input} --outfile {output} --database {params.database}
        '''

rule group4_COGs:
    '''
    This rule runs HHBlits using COGs as database for the group4 unknown proteins
    '''
    input:
        multifasta = rules.groupCharacteriser.output.clustered_unknowns
    output:
        hhsuite_concat = base_path + "/group4/cog/group4_cog.tsv",
    params:
        outdir = base_path + "/group4/cog",
        cogs = cogs_db
    conda: "envs/hhsuite.yaml"
    threads: thread_hogger
    log:
        out = base_path + "/group4/cog/logs/cog.out",
        err = base_path + "/group4/cog/logs/cog.err",
    shell:
        '''
        python3 scripts/hhblits_wrapper.py --input {input.multifasta} --output_basepath {params.outdir} --database {params.cogs}
        cat {params.outdir}/hhblits/*.tsv > {output.hhsuite_concat}
        '''

rule parse_hhsuite_group4_cog:
    '''
    Parser output from HHSuite (PDB)
    '''
    input: rules.group4_COGs.output.hhsuite_concat
    output: base_path + "/group4/cog/group4_cog_parsed.tsv"
    conda: "envs/hhsuite.yaml"
    params:
        database = "COGs",
        mapping = cogs_db + "/cog-20.def.tab"
    threads: thread_small
    shell:
        '''
        python3 scripts/hhsuite_parser.py --infile {input} --outfile {output} --database {params.database} --mapping {params.mapping}
        '''

rule group4_commonness:
    '''
    Blasts all group4 proteins against the proteomes of a CRISPR-Cas locus.
    The idea is to measure how common these proteins are and what kind of type III loci they are associated with.
    The rule first creates a diamond database for the proteome of each locus.
    Then, it blasts each group4 protein against each locus proteome using Diamond
    '''
    input:
        group4_proteins = rules.groupCharacteriser.output.clustered_unknowns,
        locus_proteins = rules.crispr_locus_proteins.output.crispr_locus_proteins
    output:
        blast_result = base_path + "/group4_prober/{c}/{c}.blast",
        diamond_db = base_path + "/group4_prober/{c}/{c}.dmnd"
    params:
        blast_result_temp = base_path + "/group4_prober/{c}/{c}.blast.temp"
    conda:
        "envs/diamond.yaml"
    threads: thread_ultrasmall
    shell:
        '''
        diamond makedb --in {input.locus_proteins} -d {output.diamond_db} --quiet
        diamond blastp --query {input.group4_proteins} --db {output.diamond_db} --outfmt 6 --out {params.blast_result_temp} --quiet
        awk -v locus={wildcards.c} '{{print $0,"\t",locus}}' {params.blast_result_temp} > {output.blast_result}
        rm {params.blast_result_temp}
        '''

rule concatenate_group4_commonness:
    '''
    Concatenates all blast results from group4_commonness into one file
    '''
    input: aggregate_group4_blasts
    output:
        concatenated_blast_results = base_path + "/group4_prober/group4_prober.blast"
    threads: thread_ultrasmall
    shell:
        '''
        echo -e "{blast_headers_group4}" > {output}
        cat {input} >> {output}
        '''

rule analyseGroup4Hits:
    input:
        blast_results = rules.concatenate_group4_commonness.output.concatenated_blast_results,
        pfam_hits = rules.groupCharacteriser.output.group4_pfam
    output:
        results_txt = base_path + "/group4_prober/group4_prober_analysis.txt"
    params:
        outpath = base_path + "/group4_prober"
    threads: thread_ultrasmall
    conda: "envs/hmmer.yaml"
    shell:
        '''
        python3 scripts/group4HitsAnalyser.py --input_blast {input.blast_results} --input_pfam_annotations {input.pfam_hits} --output {output.results_txt} --outpath {params.outpath}
        '''

rule effector_commonness:
    '''
    This rule counts the number of different effectors across the loci
    '''
    input:
        crispr_loci = rules.mastercombiner.output.final_info_table
    output:
        effector_commonness_tsv = base_path + "/13_effector_commonness/effector_commonness_master.tsv",
        #effector_commonness_plot = base_path + "/13_effector_commonness/effector_commonness.png"
    params:
        base_path = base_path + "/13_effector_commonness/effector_commonness"
    threads: thread_ultrasmall
    shell:
        '''
        python3 scripts/effector_commonness.py --input {input.crispr_loci} --output_basepath {params.base_path}
        '''

rule cctyper_gene_locations:
    '''
    This uses CCTyper outputs to generate a table that shows each cas gene's
    coordinate on the contig and other attributes. All attributes are:
    | protein_id     | start   | end     | effector | locus             | sample           | strand | type
    '''
    input:
        cas_operons_cctyper = base_path + "/071_cctyper_loci/{c}/cas_operons.tsv",
        cas_proteins = rules.crispr_locus_proteins.output.crispr_locus_proteins,
    output:
        cctyper_gene_locations_plottable = base_path + "/14_cctyper_gene_locations/{c}/{c}_cctyper_gene_locations.tsv"
    params:
        sample_folder = base_path + "/06_host_genomes/", #this is used to find the gff and proteins files that are strain-, not locus-specific
        this_folder = base_path + "/14_cctyper_gene_locations/{c}",
        outputfolder = base_path + "/14_cctyper_gene_locations/{c}",
        cctyper_folder = base_path + "/07_cctyper"
    log:
        out = base_path + "/14_cctyper_gene_locations/logs/{c}.out",
        err = base_path + "/14_cctyper_gene_locations/logs/{c}.err"
    threads: thread_ultrasmall
    shell:
        '''
        python scripts/cctyper_gene_locations.py --locus_id {wildcards.c} --sample_folder {params.sample_folder} --this_folder {params.this_folder} --outputfolder {params.outputfolder} --cas_operons {input.cas_operons_cctyper} --cctyper_path {params.cctyper_folder} --protein_fasta {input.cas_proteins} 2> {log.err} 1> {log.out}
        '''

rule locus_visualiser:
    '''
    This rule visualises a given range of a gff file.
    Using locus as wildcard.
    GFF comes from genome.
    Coordinates for CRISPR locus (and beyond) comes from 
    Takes CRISPR coordinates from 
    '''
    input:
        #genome_gff = rules.crispr_locus_proteins.output.crispr_locus_gff,
        crispr_positive_samples = base_path + "/071_cctyper_loci/{c}/cas_operons.tsv",
        cctyper_plottable = rules.cctyper_gene_locations.output.cctyper_gene_locations_plottable,
        plottable_effector_positions = rules.validated_new_effectors_analysis.output.plottable_effector_positions,
        known_effectors = rules.cATyper_analysis.output.plottable_effector_positions,
        ring_nucleases = rules.ring_nucleases_analysis.output.plottable_effector_positions
    output:
        visualisation = base_path + "/90_locus_viz/{c}/{c}_viz.png"
    conda:
        "envs/locus_visualiser.yaml"
    params:
        outdir = base_path + "/90_locus_viz/{c}",
        host_genomes_folder = base_path + "/06_host_genomes",
        cctyper_folder = base_path + "/07_cctyper",
    log:
        out = base_path + "/90_locus_viz/logs/{c}.out",
        err = base_path + "/90_locus_viz/logs/{c}.err",
    threads: thread_ultrasmall
    shell:
        '''
        python scripts/locus_visualiser.py --locus {wildcards.c} --cas_operons_file {input.crispr_positive_samples} --output_folder {params.outdir} --host_genomes_folder {params.host_genomes_folder} --cctyper_folder {params.cctyper_folder} --validated_effectors {input.plottable_effector_positions} --cctyper_protein_table {input.cctyper_plottable} --known_effector_table {input.known_effectors} --ring_nucleases {input.ring_nucleases}  2> {log.err} 1> {log.out}
        '''

rule concatenate_locus_viz:
    input: aggregate_locus_viz
    output: base_path + "/90_locus_viz/locus_viz.done"
    threads: thread_ultrasmall
    shell:
        '''
        touch {output}
        '''


rule create_html_file:
    '''
    This creates an html file out of the mastertable and adds hyperlinks
    to the loci figures and the associated .tsv files
    '''
    input:
        mastertable = rules.mastercombiner.output.final_info_table,
    output:
        html_file = base_path + "/type_iii_mastertable.html"
    params:
        viz_dir_full = base_path + "/data",
        viz_dir_relative = "data",
    conda:
        "envs/excel_output.yaml"
    threads: thread_ultrasmall
    shell:
        '''
        python scripts/html_writer.py --mastertable {input.mastertable} --viz_dir_full {params.viz_dir_full} --viz_dir_relative {params.viz_dir_relative} --output_html {output.html_file}
        '''

rule casR_clustering:
    '''
    Creating clusters of everything annotated as CasR by CCTyper.
    Here we concatenate all CasRs, so we are not using loci-based wildcards.
    Instead we use output from concatenate_locus_viz as an indirect signal that all CasRs have been generated
    '''
    input:
        viz_done = rules.concatenate_locus_viz.output,
    output:
        clustered_CasR = base_path + "/casR_clustering/casR_clustered.faa",
        concatenated_casR = base_path + "/casR_clustering/casR_concatenated.faa",
    conda:
        "envs/trees.yaml"
    params:
        casR_wildcarded = base_path + "/14_cctyper_gene_locations/*/*casR.faa",
        cutoff = 0.4,
        n = 2,
    threads: thread_ultrasmall
    shell:
        '''
        cat {params.casR_wildcarded} > {output.concatenated_casR}
        cd-hit -i {output.concatenated_casR} -o {output.clustered_CasR} -c {params.cutoff} -n {params.n} -d 0 -M 16000 -T {threads}
        '''

rule csx19_finder:
    '''
    Csx19-labeled genes from cctyper might be ring nucleases in cA3 loci.
    This script uses a modified version of the rule unknown_finder.py to find these genes and
    extract them for further study.

    If this works OK, consider generalising this rule to any gene of interest. Now just hardcode Csx19.
    '''
    input:
        cctyper = base_path + "/071_cctyper_loci/{c}/cas_operons.tsv",
        hmm = rules.cATyper_hmm_search.output.temp_rows
    output:
        csx19_proteins = base_path + "/110_csx19/{c}/{c}_csx19.faa",
        info = base_path + "/110_csx19/{c}/{c}_csx19.tsv", #each row is a protein
        locus_info = base_path + "/110_csx19/{c}/{c}_csx19_info.tsv", #each row is a locus
    params:
        outputfolder = base_path + "/110_csx19/{c}",
        host_genomes_folder = base_path + "/06_host_genomes",
        cctyper_folder = base_path + "/07_cctyper",
        sample_folder = base_path + "/06_host_genomes/", #this is used to find the gff and proteins files that are strain-, not locus-specific
        extra_cas_db = additional_cas_proteins + "/new_effectors_additional_cas.dmnd",
    conda: "envs/gff_utils.yaml"
    log:
        out = base_path + "/110_csx19/logs/{c}.out",
        err = base_path + "/110_csx19/logs/{c}.err"
    threads: thread_ultrasmall
    shell:
        '''
        python3 scripts/csx19_finder.py --locus_id {wildcards.c} --cctyper {input.cctyper} --hmm {input.hmm} --additional_cas_db {params.extra_cas_db} --outputfolder {params.outputfolder} --info_out {output.info} --unknown_proteins_output {output.csx19_proteins} --host_genomes_folder {params.host_genomes_folder} --cctyper_path {params.cctyper_folder} --sample_folder {params.sample_folder} --output_locus_info {output.locus_info} 2> {log.err} 1> {log.out}
        touch {output.csx19_proteins}
        touch {output.info}
        '''

rule concatenate_csx19:
    input: aggregate_csx19
    output:
        proteins = base_path + "/110_csx19/csx19s.faa",
        info = base_path + "/110_csx19/csx19_info.tsv"
    threads: thread_ultrasmall
    shell:
        '''
        find '{base_path}/110_csx19' -maxdepth 2 -type f -wholename '*/*_csx19.faa' -print0 | xargs -0 cat >> {output.proteins}
        echo "locus_id\tsample\tsequence\tlength\tevalue\tposition\tid\tcctyper" > {output.info}
        find '{base_path}/110_csx19' -maxdepth 2 -type f -wholename '*/*_csx19_info.tsv' -print0 | xargs -0 cat >> {output.info}
        '''


#### Phage genomes ####
millard_fa = "/mnt/shared/scratch/vhoikkal/private/databases/millard_phages/4May2024_genomes.fa"
millard_proteins = "/mnt/shared/scratch/vhoikkal/private/databases/millard_phages/4May2024_vConTACT2_proteins.faa"
millard_metadata = "/mnt/shared/scratch/vhoikkal/private/databases/millard_phages/4May2024_data.tsv"

def aggregate_millard_phage_RN_analysis(wildcards):
    '''
    This function is used to aggregate the outputs of the millard_phage_RN_analysis rule.
    If not working, create a checkpoint to generate wildcards first.
    '''
    checkpoint_output = checkpoints.divide_millard_phages_to_folders.get(**wildcards).output[0]
    cvals = glob_wildcards(os.path.join(checkpoint_output,"{phage}/done.txt")).phage
    return expand(base_path + "/P3_hmm_analysis/{phage}/{phage}_cATyper_results.tsv", phage=cvals)


checkpoint divide_millard_phages_to_folders:
    '''
    Takes in the genome Millard genomes.fa file and creates a folder and fasta for each phage in the fasta.
    The headers are in the format >phage_name description. We only want to keep the phage_name for the folder.
    '''
    output: directory(base_path + "/P1_genomes")
    input:
        millard_fa = millard_fa, #nt multifasta of all phage genomes
        millard_proteins = millard_proteins #aa multifasta of all phage genones. Headers are >phage_name_runningnumber
    params:
        outdir = base_path + "/P1_genomes",
    threads: thread_ultrasmall
    run:
        from Bio import SeqIO
        import os
        print("Starting divide_millard_phages_to_folders")
        
        # Create a dictionary to store sequences for each phage
        phage_sequences = {}
        phage_protein_sequences = {}
        
        # Parse the millard_fa file and store sequences in the dictionary
        for record in SeqIO.parse(input.millard_fa, "fasta"):
            phage_name = record.id.split(" ")[0]
            if phage_name not in phage_sequences:
                phage_sequences[phage_name] = []
            phage_sequences[phage_name].append(record)
        
        # Parse the millard_proteins file and store sequences in the dictionary
        for record in SeqIO.parse(input.millard_proteins, "fasta"):
            phage_name = record.id.split("_")[0]
            if phage_name not in phage_protein_sequences:
                phage_protein_sequences[phage_name] = []
            phage_protein_sequences[phage_name].append(record)
        
        # Write sequences to files for each phage
        for phage_name, sequences in phage_sequences.items():
            phage_dir = f"{params.outdir}/{phage_name}"
            os.makedirs(phage_dir, exist_ok=True)
            
            # Write nucleotide sequences to file
            with open(f"{phage_dir}/{phage_name}.fna", "w") as out:
                SeqIO.write([seq for seq in sequences if seq.id.split(" ")[0] == phage_name], out, "fasta")
            
        for phage_name, sequences in phage_protein_sequences.items():
            phage_dir = f"{params.outdir}/{phage_name}"
            os.makedirs(phage_dir, exist_ok=True)

            # Write protein sequences to file
            with open(f"{phage_dir}/{phage_name}_proteins.faa", "w") as out:
                SeqIO.write([seq for seq in sequences if seq.id.split("_")[0] == phage_name], out, "fasta")
            
            # Write a "done" file in every phage folder
            with open(f"{phage_dir}/done.txt", "w") as out:
                out.write("done")

rule millard_phage_RN:
    '''
    Searches for ring nucleases in the Millard phage proteomes using HMMER and
    our pre-made RN hmm profiles similar to rule ring_nucleases
    '''
    input:
        phage_proteins = base_path + "/P1_genomes/{phage}/{phage}_proteins.faa"
    output:
        temp_rows = base_path + "/P2_hmm/{phage}/{phage}_RN_temp.tsv",
        hmm_rows = base_path + "/P2_hmm/{phage}/{phage}_RN_hmm.tsv"
    params:
        rn_db = ring_nuclease_db,
        evalue = "1e-8",
        temp_hmm = base_path + "/P2_hmm/{phage}/{phage}_RN_hmm.temp"
    conda: "envs/hmmer.yaml"
    threads: thread_ultrasmall
    log:
        out = base_path + "/P2_hmm/{phage}/logs/RN_search.out",
        err = base_path + "/P2_hmm/{phage}/logs/RN_search.err"
    shell:
        '''
        hmmscan --domtblout {params.temp_hmm} --cpu {threads} -E {params.evalue} {params.rn_db} {input.phage_proteins} &> /dev/null
        echo "Removing commented rows" >> {log.out}
        grep -v "^#" {params.temp_hmm} > {output.temp_rows} ||:
        echo "Writing header" >> {log.out}
        echo -e 'target_name\taccession\ttlen\tquery_name\taccession\tqlen\tE-value_fullseq\tscore_fullseq\tbias_fullseq\t#_domain\tof_domain\tc-Evalue_domain\ti-Evalue_domain\tscore_domain\tbias_domain\t_hmm_from\thmm_to\t_ali_from\tali_to\tenv_from\tenv_to\tacc\tdescription' > {output.hmm_rows}
        echo "Checking if hits were found" >> {log.out}
        if [ -s {output.temp_rows} ]; then 
            echo "Hits found for {wildcards.phage}" >> {log.out}
            cat {output.temp_rows} >> {output.hmm_rows}
        else
            echo "No hits found for {wildcards.phage}" >> {log.out}
            touch {output.hmm_rows}
        fi
        '''

rule millard_phage_RN_analysis:
    '''
    Analyses the ring nuclease hits in the Millard phage proteomes by reusing some code from rule ring_nucleases_analysis
    and phage-specific scripts too.
    '''
    input:
        hmm_rows = rules.millard_phage_RN.output.hmm_rows,
        proteins = base_path + "/P1_genomes/{phage}/{phage}_proteins.faa"
    output:
        catyper = base_path + "/P3_hmm_analysis/{phage}/{phage}_cATyper_results.tsv",
        hmm_targets = base_path + "/P3_hmm_analysis/{phage}/{phage}_RN_hmm_targets.tsv",
        effector_to_protein = base_path + "/P3_hmm_analysis/{phage}/{phage}_effector_to_protein.tsv",
        protein_to_effector = base_path + "/P3_hmm_analysis/{phage}/{phage}_protein_to_effector.tsv",
        plottable_effector_positions = base_path + "/P3_hmm_analysis/{phage}/{phage}_plottable_effector_positions.tsv",
    params:
        outdir = base_path + "/P3_hmm_analysis/{phage}",
        hmm_msa_folder = ring_nuclease_folder + "/profiles"
    conda: "envs/hmmer.yaml"
    threads: thread_ultrasmall
    log:
        out = base_path + "/P3_hmm_analysis/{phage}/logs/RN_analysis.out",
        err = base_path + "/P3_hmm_analysis/{phage}/logs/RN_analysis.err"
    shell:
        '''
        echo "Running ring nuclease analysis for phages" >> {log.out}
        echo "Listing targets" >> {log.out}
        ls {params.hmm_msa_folder}/*/*.hmm > {output.hmm_targets}
        python scripts/phage_RN_analyzer.py --sample {wildcards.phage} --output_folder {params.outdir} --proteins_fasta {input.proteins} --hmm_targets {output.hmm_targets} --hmm_rows {input.hmm_rows} --catyper_out {output.catyper} --catyper_type "ring_nucleases" --effector_plot_data {output.plottable_effector_positions} --ring_nuclease True 2> {log.err} 1> {log.out}
        '''


rule concatenate_millard_phage_RN_analysis:
    input:
        files = aggregate_millard_phage_RN_analysis
    output: base_path + "/P3_hmm_analysis/phages_RN_hmm.tsv"
    threads: thread_ultrasmall
    log:
        out = base_path + "/P3_hmm_analysis/logs/phages_RN_hmm.out",
        err = base_path + "/P3_hmm_analysis/logs/phages_RN_hmm.err"
    params:
        temp_file = base_path + "/P3_hmm_analysis/phages_RN_hmm.temp"
    shell:
        """
        # Get the header from the first file
        header=$(head -n 1 {input.files[0]})

        # Write the header to the output file
        echo "$header" > {output}

        # Concatenate files without header and append to the output file
        for file in {input.files}; do
            tail -n +2 "$file" >> {output}
        done
        """

rule phage_cas10:
    '''
    Searches the phage proteomes for Cas10s using the Cas10 hmm profiles
    '''
    input:
        phage_proteins = base_path + "/P1_genomes/{phage}/{phage}_proteins.faa"
    output:
        temp_rows = base_path + "/P4_cas10hmm/{phage}/{phage}_Cas10_temp.tsv",
        hmm_rows = base_path + "/P4_cas10hmm/{phage}/{phage}_Cas10_hmm.tsv"
    params:
        cas10_db = cas10_db,
        evalue = "1e-8",
        temp_hmm = base_path + "/P4_cas10hmm/{phage}/{phage}_Cas10_hmm.temp"
    conda: "envs/hmmer.yaml"
    threads: thread_ultrasmall
    log:
        out = base_path + "/P4_cas10hmm/{phage}/logs/Cas10_search.out",
        err = base_path + "/P4_cas10hmm/{phage}/logs/Cas10_search.err"
    shell:
        '''
        hmmscan --domtblout {params.temp_hmm} --cpu {threads} -E {params.evalue} {params.cas10_db} {input.phage_proteins} &> /dev/null
        echo "Removing commented rows" >> {log.out}
        grep -v "^#" {params.temp_hmm} > {output.temp_rows} ||:
        echo "Writing header" >> {log.out}
        echo -e 'target_name\taccession\ttlen\tquery_name\taccession\tqlen\tE-value_fullseq\tscore_fullseq\tbias_fullseq\t#_domain\tof_domain\tc-Evalue_domain\ti-Evalue_domain\tscore_domain\tbias_domain\t_hmm_from\thmm_to\t_ali_from\tali_to\tenv_from\tenv_to\tacc\tdescription' > {output.hmm_rows}
        echo "Checking if hits were found" >> {log.out}
        if [ -s {output.temp_rows} ]; then 
            echo "Hits found for {wildcards.phage}" >> {log.out}
            cat {output.temp_rows} >> {output.hmm_rows}
        else
            echo "No hits found for {wildcards.phage}" >> {log.out}
            touch {output.hmm_rows}
        fi
        '''

def aggregate_phage_cas10(wildcards):
    checkpoint_output = checkpoints.divide_millard_phages_to_folders.get(**wildcards).output[0]
    cvals = glob_wildcards(os.path.join(checkpoint_output,"{phage}/done.txt")).phage
    return expand(base_path + "/P4_cas10hmm/{phage}/{phage}_Cas10_hmm.tsv", phage=cvals)

rule phage_cas10_concatenate:
    input: aggregate_phage_cas10
    output: base_path + "/P4_cas10hmm/phages_Cas10_hmm.tsv"
    threads: thread_ultrasmall
    log:
        out = base_path + "/P4_cas10hmm/logs/phages_Cas10_hmm.out",
        err = base_path + "/P4_cas10hmm/logs/phages_Cas10_hmm.err"
    params:
        temp_file = base_path + "/P4_cas10hmm/phages_Cas10_hmm.temp"
    shell:
        """
        # Get the header from the first file
        header=$(head -n 1 {input.files[0]})
        
        # Concatenate files without header
        tail -q -n +2 {input.files} > {params.temp_file}

        # Add header to the concatenated file
        echo "$header" > {output}

        # Concatenate the rest of the files
        cat {params.temp_file} >> {output}

        # Remove the temporary file
        rm {params.temp_file}
        """

#finish writing analysis if any cas10s are found
# rule phage_cas10_analysis:
#     '''
#     Takes the hmm hits from rule phage_cas10 and creates simplified tables similar to
#     ring nuclease rules before
#     '''
#     input:
#         hmm_rows = rules.phage_cas10.output.hmm_rows,
#         proteins = base_path + "/P1_genomes/{phage}/{phage}_proteins.faa"
#     output:
#         catyper = base_path + "/P5_cas10_analysis/{phage}/{phage}_cATyper_results.tsv",
#         hmm_targets = base_path + "/P5_cas10_analysis/{phage}/{phage}_Cas10_hmm_targets.tsv",
#         effector_to_protein = base_path + "/P5_cas10_analysis/{phage}/{phage}_effector_to_protein.tsv",
#         protein_to_effector = base_path + "/P5_cas10_analysis/{phage}/{phage}_protein_to_effector.tsv",
#         plottable_effector_positions = base_path + "/P5_cas10_analysis/{phage}/{phage}_plottable_effector_positions.tsv",
#     params:
#         outdir = base_path + "/P5_cas10_analysis/{phage}",
#         hmm_msa_folder = cas10_folder + "/profiles"
#     conda: "envs/hmmer.yaml"
#     threads: thread_ultrasmall
#     log:
#         out = base_path + "/P5_cas10_analysis/{phage}/logs/Cas10_analysis.out",
#         err = base_path + "/P5_cas10_analysis/{phage}/logs/Cas10_analysis.err"
#     shell:
#         '''
#         echo "Running Cas10 analysis for phages" >> {log.out}
#         echo "Listing targets" >> {log.out}
#         ls {params.hmm_msa_folder}/*/*.hmm > {output.hmm_targets}
#         python scripts/phage_RN_analyzer.py --sample {wildcards.phage} --output_folder {params.outdir} --proteins_fasta {input.proteins} --hmm_targets {output.hmm_targets} --hmm_rows {input.hmm_rows} --catyper_out {output.catyper} --catyper_type "Cas10" --effector_plot_data {output.plottable_effector_positions} --ring_nuclease False 2> {log.err} 1> {log.out}
#         '''


rule spacers:
    input:
        done = rules.CRISPRCasTyper.output,
        renaming_done = rules.CRISPRCasTyper_rename.output
    params:
        fastafolder = base_path + "/07_cctyper/{j}/spacers",
        fastafile = base_path + "/07_cctyper/{j}/spacers/",
        crisprresultfolder = base_path + "/07_cctyper/{j}/"
    output:
        base_path + "/03_spacers/all/{j}_spacers.csv"
    run:
        from Bio import SeqIO
        import os.path
        from os import path
        import csv

        crispr_loci = [] #create a crispr_loci list. This is just a list of cctyper loci.
        crispr_results = str(params.crisprresultfolder + "crisprs_all.tab")
        print("Opening " + crispr_results)
        if os.path.isfile(crispr_results):
            with open(crispr_results) as csv_file:
                csv_reader = csv.reader(csv_file, delimiter='\t')
                line_count = 0
                for row in csv_reader:
                    if line_count == 0:
                        #print(f'Column names are {", ".join(row)}')
                        line_count += 1
                    else:
                        #print(", ".join(row))
                        #print(row[11])
                        if row[11] == "True":
                            #print(row[1] + " trusted")
                            crispr_loci.append(row[1])
                        line_count += 1
                #print(f'Processed {line_count} lines.')
                #print("Trusted loci: " +  str(crispr_loci))
        else:
            print("No CRISPR-Cas loci found")

        fastafolder_object = Path(params.fastafolder) #store path to spacer folder as object
        shell("touch '{base_path}/03_spacers/all/{wildcards.j}_spacers.csv'") #create spacer .csv file for every sample
        if len(crispr_loci) > 0: #check if there are any loci in the trusted list
            print("------------- Spacers found in  " + str(wildcards.j) + "!")
            spacers_dict = {}
            for filename in os.listdir(fastafolder_object): #iterate through all fastas
                filename_reduced = re.search('(.*)\.[^.]+$', filename).group(1)
                #print("Reduced locus name: " + str(filename_reduced))
                #print(crispr_loci)
                #print(filename)
                if filename_reduced in crispr_loci:
                    #print("Match for file " + str(filename) + " in trusted loci")
                    for seq_record in SeqIO.parse(str(fastafolder_object) + '/' + filename, 'fasta'): #iterate through all entries in fasta
                        spacers_dict[seq_record.id] = seq_record.seq #add entry as dictionary entry
            with open(str(output), 'w') as f:
                for key in spacers_dict.keys():
                    f.write('%s,%s\n'%(key + "@" + str(wildcards.j),spacers_dict[key]))
        else:
            print("No spacers found in " + str(wildcards.j))

def aggregate_spacers(wildcards):
    '''
    Returns all spacers from all samples
    '''
    checkpoint_output = checkpoints.expand_host_genomelist_to_wildcards_postcas10.get(**wildcards).output[0]
    print(checkpoint_output)
    ivals = glob_wildcards(os.path.join(checkpoint_output,"{j}.txt")).j
    print("Aggregating spacers from trusted locus observations")
    return expand(base_path + "/03_spacers/all/{j}_spacers.csv", j=ivals)

rule combined_spacers:
    '''
    This rule aggregates all spacer.csv files into one big .csv collection. The expand function is used for
    the concatenation, and it uses the list returned by the aggregate_IDs function.
    '''
    input: aggregate_spacers
    output: base_path + "/03_spacers/combined_spacers.csv"
    shell:
        """
        cat {input} > {output}
        """

rule spacers_blast_convert_spacers:
    """
    This converts the concatenated .csv into fasta format for blasting
    """
    input:
        rules.combined_spacers.output #the concatenated csv file
    output:
        fasta = base_path + "/03_spacers/combined_spacers.fasta",
        formatted_spacers = base_path + "/03_spacers/formatted_spacers.csv"
    threads: thread_ultrasmall
    shell:
        """
        sed 's/^/>/' {input} > {output.formatted_spacers}
        tr "," "\n" < {output.formatted_spacers} > {output.fasta}
        """

rule spacers_blast_runblast:
    '''
    Blasts spacers against the phage database.
    Also adds headers to blast result file.

    Remember to make a blast db of the phage genomes first.
    '''
    input:
        fasta = rules.spacers_blast_convert_spacers.output.fasta #the concatenated csv file
    output:
        blast = base_path + "/P5_spacer_to_phage_blast/blast.out",
    params:
        blast_db = millard_fa,
    conda: "envs/blast.yaml"
    threads: thread_hogger
    shell:
        '''
        blastn -num_threads {threads} -db {params.blast_db} -task blastn-short -query {input.fasta} -outfmt "6 qseqid qlen sseqid  stitle staxids saccver sstart send evalue length pident mismatch gapopen gaps sstrand" -evalue 0.01 > {output.blast}
        echo -e 'Query_id\tQuery_length\tSubject_id\tSubject_sci_title\tSubject_taxID\tSubject_accession_ID_version\tSubject_start\tSubject_end\tEvalue\tLength\tPercentage_identical\tMismatches\tOpen_gaps\tGaps\tSubject_strand' | cat - {output.blast} > blast.cat && mv blast.cat {output.blast}      
        '''

def aggregate_allHosts(wildcards):
    '''
    Returns all hosts after the filtering by cas10
    '''
    checkpoint_output = checkpoints.expand_host_genomelist_to_wildcards_postcas10.get(**wildcards).output[0]
    ivals = glob_wildcards(os.path.join(checkpoint_output,"{j}.txt")).j
    return expand(base_path + "/03_postfiltering_genome_wildcards/{j}.txt", j=ivals)

rule hostTable:
    '''
    Prints all bacterial genome names to file without extension. Used in spacer mapping analysis.
    This table is after filtering by cas10 ("03_postfiltering_genome_wildcards"), so all genomes are type III,
    but may or may not have spacer matches against phages.
    '''
    input: aggregate_allHosts
    output: base_path + "/allHosts.tsv"
    run:
        hostlist = str(input).split(" ")
        basenameList = []
        for h in hostlist:
            basename = Path(h).stem
            basenameList.append(basename)
        #print(input)
        with open(str(output), mode='wt', encoding='utf-8') as hostfile:
            hostfile.write('\n'.join(basenameList))


rule spacer_hits_analyzer:
    '''
    Using data from hosts, spacer, protospacer, and phage databases, this script
    creates merged tables showing phage/host interactions and network tables
    for visualisation with Gephi.
    '''
    input:
        mastertable = rules.mastercombiner.output.final_info_table, #ok
        phage_data = millard_metadata, #ok
        allHosts = rules.hostTable.output,
        spacer_blast_hits = rules.spacers_blast_runblast.output.blast,
        phage_rns = rules.concatenate_millard_phage_RN_analysis.output,
    output:
        spacer_links = base_path + "/P6_spacer_analysis/spacer_protospacer_hits.csv",
        gephi_nodes = base_path + "/P6_spacer_analysis/spacer_nodes.csv",
        gephi_edges = base_path + "/P6_spacer_analysis/spacer_edges.csv",
        ca3_phages_by_host = base_path + "/P6_spacer_analysis/ca3_phagecounts.tsv", #all phages that infect ca3 hosts
        non_ca3_phages_by_host = base_path + "/P6_spacer_analysis/ca3_phagecounts_false.tsv" #all phages that do not infect ca3 hosts
    conda:
        "envs/merge.yaml"
    params:
        outfolder = base_path + "/P6_spacer_analysis"
    shell:
        """
        python scripts/spacer_hits_analyzer.py --spacer_blast_hits {input.spacer_blast_hits} --out_nodes {output.gephi_nodes} --out_edges {output.gephi_edges} --out_merged_all {output.spacer_links} --phage_data {input.phage_data} --host_mastertable {input.mastertable} --all_hosts {input.allHosts} --cATyper {input.cATyper_hmm} --base_path {base_path} --outfolder {params.outfolder} --phage_rns {input.phage_rns}
        """


rule separate_coa_proteomes_phages:
    '''
    This rule is the first step in doing enrichment analysis with the goal
    of finding phage proteins that are enriched in phages that are targeted
    by cA3-utilising hosts.

    The division into cOA-specific phages is done already in the rule spacer_hits_analyser.
    Here, we take those lists, find phages that are uniquely targeted by
    cA3-hosts, and create a fasta file of their proteins.

    Clustering is done in the next rule.
    '''
    input:
        ca3_phages = rules.spacer_hits_analyser.output.ca3_phages_by_host, #list of phages that infect ca3-hosts
        non-ca3_phages = rules.spacer_hits_analyser.output.non_ca3_phages_by_host, #list of phages that infect non-ca3 hosts
        all_phage_proteins = millard_proteins
    output:
        ca3_proteins = base_path + "/P7_enrichment_analysis/ca3_phages_proteins.faa",
        non_ca3_proteins = base_path + "/P7_enrichment_analysis/non_ca3_phages_proteins.faa"
    threads: thread_ultrasmall
    shell:
        '''
        python scripts/separate_coa_proteomes_phages.py --ca3_phages {input.ca3_phages} --non_ca3_phages {input.non-ca3_phages} --all_phage_proteins {input.all_phage_proteins} --out_ca3_proteins {output.ca3_proteins} --out_non_ca3_proteins {output.non_ca3_proteins}
        '''

rule cluster_ca3_phagge_proteomes:
    '''
    This rule creates clusters of ca3 phage proteomes
    '''
    input:
        ca3_proteins = rules.separate_coa_proteomes_phages.output.ca3_proteins,
        non_ca3_proteins = rules.separate_coa_proteomes_phages.output.non_ca3_proteins
    output:
        ca3_clusters = base_path + "/P7_enrichment_analysis/ca3_phages_clusters.tsv",
    threads: thread_hogger
    params:
        cluster_cutoff = 0.4,
        word_length = 2
    shell:
        '''
        cd-hit -i "{input.ca3_proteins}" -o "{output.ca3_clusters}" -c {params.cluster_cutoff} -n {params.word_length} -T {threads} 2> {log.err} 1> {log.out}
        '''

rule cluster_nonca3_phage_proteomes:
    '''
    This rule creates clusters of non-ca3 phage proteomes
    '''
    input:
        ca3_protein_cluster_done = rules.cluster_ca3_phagge_proteomes.output.ca3_clusters,
        non_ca3_proteins = rules.separate_coa_proteomes_phages.output.non_ca3_proteins
    output:
        non_ca3_clusters = base_path + "/P7_enrichment_analysis/non_ca3_phages_clusters.tsv"
    threads: thread_hogger
    shell:
        '''
        cd-hit -i "{input.non_ca3_proteins}" -o "{output.ca3_clusters}" -c {params.cluster_cutoff} -n {params.word_length} -T {threads} 2> {log.err} 1> {log.out}
        '''

rule enrichment_analysis:
    '''
    This rule performs enrichment analysis on the phage proteins that are targeted by cA3 hosts.
    
    Specifically, it looks for phage proteins that are enriched in cA3-targeted phages compared to non-cA3-targeted phages by using
    the Fisher's exact test.

    The output is a table of enriched phage proteins.
    '''
    input:
        ca3_clusters = rules.cluster_ca3_phagge_proteomes.output.ca3_clusters,
        non_ca3_clusters = rules.cluster_nonca3_phage_proteomes.output.non_ca3_clusters
    output:
        enrichment_results = base_path + "/P7_enrichment_analysis/enrichment_results.tsv"
    shell:
        '''
        python scripts/enrichment_analysis.py --ca3_clusters {input.ca3_clusters} --non_ca3_clusters {input.non_ca3_clusters} --out {output.enrichment_results}
        '''


rule final:
    input:
        tax = base_path + "/06_host_genomes/taxInfo.txt",
        tree_Cas10 = rules.Cas10_tree.output,
        type_iii_info = rules.combine_GGDD_HMM_to_mastertable.output.final_info_table,
        mastercombiner_final_info_table = rules.mastercombiner.output.final_info_table,
        concat_taxInfo = rules.concat_taxInfo.output,
        catyper_hmm = rules.concatenate_cATyper_hmm.output,
        catyper_analysis = rules.concatenate_cATyper_analysis.output,
        #tree_CorA = rules.CorA_tree.output,
        #tree_CorA_unclustered = rules.CorA_tree_unclustered.output,
        clustered_unknowns = rules.cluster_unknowns.output.proteins, #each row is an unknown protein
        clustered_unknowns_info = rules.concatenate_unknowns_locus_info.output.info, #contains locus-specific info on unknowns
        cas10_HD_faa = rules.Cas10_HD_hmm_maker.output.faa,
        cas10_hd_hmm = rules.Cas10_HD_hmmer.output,
        cas10_HD_hmm_merged = rules.merge_HD_hmm_with_info.output.merged_table,
        cas10_GGDD_faa = rules.Cas10_GGDD_hmm_maker.output.faa,
        cas10_GGDD_hmm = rules.Cas10_GGDD_hmmer.output,
        cas10_GGDD_hmm_merged = rules.merge_GGDD_hmm_with_info.output.merged_table,
        groupCharacteriser = rules.groupCharacteriser.output.group4_pfam,
        aggregate_crispr_locus_proteins = rules.aggregate_crispr_locus_proteins.output,
        concatenated_blast_results = rules.concatenate_group4_commonness.output,
        #analyseGroup4Hits = rules.analyseGroup4Hits.output.results_txt,
        effector_commonness = rules.effector_commonness.output.effector_commonness_tsv,
        effector_scores_summary = rules.analyse_cATyper_effector_scores.output.effector_scores_summary,
        validated_effectors_scores_summary = rules.analyse_validated_new_effectors_scores.output.effector_scores_summary,
        concatenate_effector_wilcards = rules.concatenate_effector_wildcards.output,
        concatenate_cATyper_hmm_hhsuite = rules.concatenate_cATyper_hmm_hhsuite.output,
        parse_hhsuite = rules.concatenate_cATyper_hmm_hhsuite.output,
        parse_hhsuite_cogs = rules.concatenate_cATyper_hmm_hhsuite_cogs.output,
        #group4_pdb = rules.parse_hhsuite_group4_pdb.output,
        #group4_cog = rules.parse_hhsuite_group4_cog.output,
        concatenate_validated_new_effectors_analysis = rules.concatenate_validated_new_effectors_analysis.output,
        heatmap_known_validated_effectors = rules.heatmap_known_validated_effectors.output.effector_scores_summary,
        #node_graph = rules.node_graph.output.edges,
        concatenate_locus_viz = rules.concatenate_locus_viz.output,
        #excel = rules.create_excel_file.output.excel_file,
        html = rules.create_html_file.output.html_file,
        #R_HD = rules.R_HD.output.HD_histogram,
        casR_cluster = rules.casR_clustering.output.clustered_CasR,
        heatmap_ring_nucleases = rules.heatmap_ring_nucleases.output.effector_scores_summary,
        ring_fusions = rules.ring_fusions.output.ring_fusions,
        ring_fusions_cas10 = rules.ring_nuclease_cas10_fusions.output.hmm_rows,
        validated_effectors_cas10_fusions = rules.validated_effectors_cas10_fusions.output.hmm_rows,
        known_effectors_cas10_fusions = rules.known_effectors_cas10_fusions.output.hmm_rows,
        ring_nucleases_fusions = rules.concatenate_ring_nuclease_fusions.output,
        cOA_RN_explorer = rules.cOA_RN_explorer.output.ca3_loci,
        #cA3_small_unk_proteins_info = rules.cOA_small_unk_proteins.output.cA3_small_unk_proteins_info,
        #fetch_small_unk_targets_ncbi = rules.fetch_small_unk_targets_ncbi.output.all_proteins,
        #clustered_small_unknowns = rules.cluster_small_unk_proteins.output.clustered_proteins,
        #small_unk_proteins_clusters_hmmer = rules.small_unk_proteins_clusters_hmmer.output.raw_table,
        #tree_small_unk_protein_cluster_representatives = rules.tree_small_unk_protein_cluster_representatives.output,
        #webflags_cluster_small_unk_proteins = rules.webflags_cluster_small_unk_proteins.output.done,
        #tree_small_unk_proteins_all = rules.tree_small_unk_proteins_all.output,
        #concatenate_wildcarded_webflags = rules.concatenate_wildcarded_webflags.output,
        #ae1_family_tree = rules.ae1_family_tree.output,
        #ae1_family_collector = rules.ae1_family_collector.output.done,
        #potential_other_rn_collector = rules.potential_other_rn_collector.output.done,
        concatenate_csx19 = rules.concatenate_csx19.output.proteins,
        crn4_tree = rules.crn4_tree.output.tree1,
        #phages_merge_hmm_and_metadata = rules.phages_merge_hmm_and_metadata.output
        concatenate_millard_phage_RN_analysis = rules.concatenate_millard_phage_RN_analysis.output,
        #spacers_blast_convert_spacers = rules.spacers_blast_convert_spacers.output.fasta,
        spacers_blast = rules.spacers_blast_runblast.output.blast,
        #phage_cas10_concatenate = rules.phage_cas10_concatenate.output,
        hostTable = rules.hostTable.output,
        spacer_links = rules.spacer_hits_analyzer.output.spacer_links,
    output: base_path + "/done"
    threads: thread_small
    shell:
        '''
        touch {output}
        mkdir -p {base_path}/R
        mkdir -p {base_path}/R/effector_trees

        mkdir -p {base_path}/R/effector_hmm
        mkdir -p {base_path}/R/effector_hmm/pfam
        mkdir -p {base_path}/R/effector_hmm/pdb
        mkdir -p {base_path}/R/effector_hmm/cog

        mkdir -p {base_path}/R/group4
        mkdir -p {base_path}/R/group4/pfam
        mkdir -p {base_path}/R/group4/pdb
        mkdir -p {base_path}/R/group4/cog
        mkdir -p {base_path}/R/group4/carfsaved
        
        cp {base_path}/12_Cas10_tree/cas10_tree.txt {base_path}/R
        cp {base_path}/06_host_genomes/taxInfo.txt {base_path}/R
        cp {base_path}/09_crispr_iii_CorA/loci/type_iii_info.tsv {base_path}/R
        cp {base_path}/30_unknown_effectors/unknowns_info_loci.tsv {base_path}/R
        cp {base_path}/13_effector_commonness/effector_commonness_master.tsv {base_path}/R
        cp {base_path}/mastertable_v2.tsv {base_path}/R

        cp -r {base_path}/45_effector_tree/* {base_path}/R/effector_trees

        cp -r {base_path}/43_effector_hmmer_analysis/*/*_sorted_filtered_mapped.tsv {base_path}/R/effector_hmm/pfam
        cp -r {base_path}/42_effector_hhsuite/*/*_parsed.tsv {base_path}/R/effector_hmm/pdb
        cp -r {base_path}/42_effector_hhsuite_cogs/*/*_parsed_cogs.tsv {base_path}/R/effector_hmm/cog
        
        cp -r {base_path}/group4/group4.CARFSAVED.info.tsv {base_path}/R/group4/carfsaved
        cp -r {base_path}/group4/group4.pfam.info.tsv {base_path}/R/group4/pfam
        cp -r {base_path}/group4/cog/group4_cog_parsed.tsv {base_path}/R/group4/cog
        cp -r {base_path}/group4/pdb/group4_pdb_parsed.tsv {base_path}/R/group4/pdb
        cp -r {base_path}/group4_prober/group4_prober_analysis.txt {base_path}/R/group4
        '''