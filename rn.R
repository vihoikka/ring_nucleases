#### Ring nuclease R script ####
# 
# This script takes as input the 'R' folder output by the rn.smk pipeline.
# Copy the 'R' folder in the same folder as this script. Your "project" is now
# called R (modify the project name below and the folder name if you want to
# run the script for another run of the pipeline with different output files).
#
# NOTE: there are some rows here that are specific to a certain locus
# in a certain run of the pipeline. As the naming of the loci is stochastic,
# these will probably not work for your run of the pipeline. Such rows have
# been marked in this script.
# 

#### Libraries ####

#install if missing install.packages("UpSetR")
if (!requireNamespace("UpSetR", quietly = TRUE)) {
  install.packages("UpSetR")
}

#same for devtools::install_github("krassowski/complex-upset")
if (!requireNamespace("complex-upset", quietly = TRUE)) {
  devtools::install_github("krassowski/complex-upset")
  install.packages('ComplexUpset')
}

#same for install.packages("ggcorrplot")
if (!requireNamespace("ggcorrplot", quietly = TRUE)) {
  install.packages("ggcorrplot")
}

#same for ggtreeExtra
if (!requireNamespace("ggtreeExtra", quietly = TRUE)) {
  install.packages("ggtreeExtra")
}

#same for ggtree
if (!requireNamespace("ggtree", quietly = TRUE)) {
  install.packages("ggtree")
}

# install.packages("ggnewscale")
if (!requireNamespace("ggnewscale", quietly = TRUE)) {
  install.packages("ggnewscale")
}

# install.packages("ggalt")
if (!requireNamespace("ggalt", quietly = TRUE)) {
  install.packages("ggalt")
}

library("ggalt")
library("tidyr")
library("dplyr")
library("UpSetR")
library("ComplexUpset")
library("stringr")

library("ggplot2")
library("ggcorrplot")
library("ggtree")
library("ggtreeExtra")
library("ggnewscale")

library("treeio")
library("phytools")


#### Load data ####
project = "projectname" #modify this to the project name. The folder must be in the data folder (relative to current R script)

#define paths
info_path <- paste("data/",project,"/mastertable_v2.tsv", sep = "")
cas10_tree_path <- paste("data/",project,"/cas10_tree.txt", sep = "")
effectors <- paste("data/", project, "/effector_commonness_master.tsv", sep = "")
tax_info_path <- paste("data/", project, "/taxInfo.txt", sep = "")
phage_rn <- paste("data/", project, "/phages_RN_hmm.tsv", sep = "")
rn_outside_crispr <- paste("data/", project, "/hits_outside_crispr.tsv", sep = "")


#open files
cas10_tree <- read.newick(cas10_tree_path)
cas10_tree <- midpoint.root(cas10_tree)
effectors <- read.csv(effectors, sep="\t")
tax_info <- read.csv(tax_info_path, sep="\t") #produced by a separate rule in the Snakemake pipeline
info <- read.csv(info_path, sep="\t")
rn_outside_crispr <- read.csv(rn_outside_crispr, sep="\t")
phage_rn <- read.csv(phage_rn, sep="\t")


info <- subset(info, info$Cas10_length > 499) #this was previously set to 500, resulting in a mismatch between tree and info df

#### Settings ####

colorBlindBlack8  <- c("#404040", "#E69F00", "#56B4E9", "#009E73", 
                       "#F0E442", "#CC79A7", "#D55E00", "#CC79A7")

#make up a similar color scheme with 18 colors
colorBlindBlack18 <- c("#404040", "#E69F00", "#56B4E9", "#009E73", 
                       "#F0E442", "#CC79A7", "#D55E00", "#CC79A7",
                       "#404040", "#E69F00", "#56B4E9", "#009E73", 
                       "#F0E442", "#CC79A7", "#D55E00", "#CC79A7",
                       "#404040", "#E69F00")

#### Tidy data ####

#create folder plots
if (!dir.exists(paste0("data/",project,"/plots"))) {
  dir.create(paste0("data/",project,"/plots"))
}

#create folder plots/effectors
if (!dir.exists(paste0("data/",project,"/plots/effectors"))) {
  dir.create(paste0("data/",project,"/plots/effectors"))
}

#Reduce hybrid naming to the type III component
info <- info %>%
  mutate(Subtype = str_extract_all(Subtype, "III-\\w"), # extract all "III-letter"'s and get a list of character vectors
         Subtype = sapply(Subtype, paste0, collapse=", ")) # concatenate multiple III sequences if they exist in one string

#In cases where the subtype is "III-A, III-D" rename to III-D
info$Subtype[info$Subtype == "III-A, III-D"] <- "III-D"

#mark thermophiles based on temperature
info$thermophile <- FALSE
info$thermophile[info$mean_temperature > 45] <- TRUE

#Mark hyperthermophiles
info$hyperthermophile <- FALSE
info$hyperthermophile[info$genus == "Thermus" | info$genus == "Thermotoga" | info$genus == "Aquifex"] <- TRUE

#change variable types
info$has_known_effector <- as.logical(info$has_known_effector)
info$has_validated_new_effector <- as.logical(info$has_validated_new_effector)
info$multiple_signals <- as.logical(info$multiple_signals)

#make locus GCF_003201765.2_2 as III-D. NOTE: this is specific to a certain run of the pipeline
info$Subtype[info$Locus == "GCF_003201765.2_2"] <- "III-D"

#Rename hybrid loci
info$Hybrid <- FALSE #no hybrids is default
info$Hybrid[str_detect(info$Subtype, 'Hybrid') == TRUE] <- TRUE #mark some as hybrids
info$Subtype[str_detect(info$Subtype, 'Hybrid') == TRUE & str_detect(info$Subtype, 'III-A') == TRUE] <- "III-A"
info$Subtype[str_detect(info$Subtype, 'Hybrid') == TRUE & str_detect(info$Subtype, 'III-B') == TRUE] <- "III-B"
info$Subtype[str_detect(info$Subtype, 'Hybrid') == TRUE & str_detect(info$Subtype, 'III-C') == TRUE] <- "III-C"
info$Subtype[str_detect(info$Subtype, 'Hybrid') == TRUE & str_detect(info$Subtype, 'III-D') == TRUE] <- "III-D"

info$effector_count_known_new_sum <- as.numeric(info$effector_count_known_new_sum)

#create a boolean column version of the above if > 0
info$has_effector <- FALSE
info$has_effector[info$effector_count_known_new_sum > 0] <- TRUE

#transform to some strings to boolean
info$Cas10_HD <- as.logical(info$Cas10_HD)

info$Cas10_GGDD <- as.logical(info$Cas10_GGDD)
#info$unknown_proteins <- as.logical(info$unknown_proteins)

info$Cas10_HD_coord <- as.numeric(info$Cas10_HD_coord)

info$effector_elsewhere_in_sample_but_not_here <- as.logical(info$effector_elsewhere_in_sample_but_not_here)
info$effector_elsewhere_in_sample <- as.logical(info$effector_elsewhere_in_sample)
info$multisignal_sample <- as.logical(info$multisignal_sample)
info$GGDD_hmm_boolean <- as.logical(info$GGDD_hmm_boolean)
info$cyclase <- as.logical(info$cyclase)

#same for "ca3", "ca4", "ca6", "sam.amp"
info$ca3 <- as.logical(info$ca3)
info$ca4 <- as.logical(info$ca4)
info$ca6 <- as.logical(info$ca6)
info$sam.amp <- as.logical(info$sam.amp)
info$can3 <- as.logical(info$can3)

#create a column locus_signal that shows which signal is true in the locus
info$locus_signal <- "none"
info$locus_signal[info$ca3 == TRUE] <- "ca3"
info$locus_signal[info$ca4 == TRUE] <- "ca4"
info$locus_signal[info$ca6 == TRUE] <- "ca6"
info$locus_signal[info$sam.amp == TRUE] <- "sam.amp"

#in cases where ca6 and sam.amp are true, make locus_signal ca6
info$locus_signal[info$sam.amp == TRUE & info$can3 == TRUE] <- "ca6"

nrow(subset(info, locus_signal == "ca6"))

info_samamp_ca6 <- subset(info, sam.amp == TRUE & ca6 == TRUE)
info_cora_csm62 <- subset(info, cora == TRUE & can3 == TRUE)

# in info$signal_sample, fill empty cells with FALSE
info$signal_sample <- as.logical(info$signal_sample)
info$signal_sample[is.na(info$signal_sample)] <- FALSE
#change all columns to logical
info[,c("can3", "tirsaved", "cam1", "cam2", "cam3", "saved.chat", "nucc", "calpl", "cami1", "can1.2", "csx1", "csx23", "csm6.ca6", "cora", 'csx16', 'csx20', 'crn1', 'crn2', 'crn3', 'csx15')] <- lapply(info[,c("can3", "tirsaved", "cam1", "cam2", "cam3", "saved.chat", "nucc", "calpl", "cami1", "can1.2", "csx1", "csx23", "csm6.ca6", "cora", 'csx16', 'csx20', 'crn1', 'crn2', 'crn3', 'csx15')], as.logical)
#same for con_ca3, con_ca4, con_ca6, con_sam_amp
info[,c("con_ca3", "con_ca4", "con_ca6", "con_sam.amp")] <- lapply(info[,c("con_ca3", "con_ca4", "con_ca6", "con_sam.amp")], as.logical)

info$has_multiple_effectors <- ifelse(info$effector_count_known_new_sum > 1, TRUE, FALSE)


#set rownames by locus
rownames(info) <- info$Locus

#setup Cas10 length plot for Cas10 tree
info_for_cas10length <- select(info, Locus, Cas10_length, Cas10_HD, mean_temperature)
#info_for_cas10length <- na.omit(info_for_cas10length)

#this renaming must be done. Otherwise will get an error "Error: object 'Cas10_length' not found"
info_for_cas10length <- dplyr::rename(info_for_cas10length, Cas10_length_plot = Cas10_length)
info_for_cas10length <- dplyr::rename(info_for_cas10length, Cas10_HD_plot = Cas10_HD)
info_for_cas10length <- dplyr::rename(info_for_cas10length, mean_temperature_plot = mean_temperature)

#for rn outside crispr data, make within_crispr and within_prophage columns logical
rn_outside_crispr$within_crispr <- as.logical(rn_outside_crispr$within_crispr)
rn_outside_crispr$within_prophage <- as.logical(rn_outside_crispr$within_prophage)

#make has_ring_nuclease into logical
info$has_ring_nuclease <- as.logical(info$has_ring_nuclease)

#count the number of loci with any effector
effector_count <- sum(info$has_effector)

#count the number of ring nucleases
ring_nuclease_count <- sum(info$has_ring_nuclease)


#in info, add column ring_nuclease_count and count number of RNs
info$ring_nuclease_count <- rowSums(info[,c("crn1", "crn2", "crn3", "csx15", "csx16", "csx20")])


#### 1. UPSET EFFECTOR AND RING NUCLEASES ####
#When subsetting for all, archaea or thermophiles, change the following:
# - The subsetting function at the top (first thing in this section)
# - Signal molecule associations at the top of this section
# - The name of the plot file at the end of this section
# - Title of the plot in the plot script


upset_info <- info

#filter only Archaea
#upset_info <- subset(info, domain == "Archaea")

str(upset_info)
#upset_info <- subset(info, domain == "Archaea") #uncomment this for running only for archaea
#upset_info <- subset(info, mean_temperature >= 60) #uncomment this for running only for thermophiles)

info_multiple_loci_effectors <- upset_info %>%
  pivot_longer(cols = c(can3, tirsaved, cam1, cam2, cam3, saved.chat, tirsaved, nucc, calpl, cami1, can1.2, csx1, csx23, csm6.ca6, cora,
                        csx16, csx20, crn1, crn2, crn3, csx15, csx16, csx20, blender),#, phrogRN, proteaseRN, solosavedRN
               names_to = "Effector",
               values_to = "Presence")

str(info_multiple_loci_effectors)

crn1s <- subset(info_multiple_loci_effectors, Effector == "crn1" & Presence == TRUE)
crn1s$Subtype

#presence as factor
info_multiple_loci_effectors$Presence <- as.factor(info_multiple_loci_effectors$Presence)

# Sort out duplicates and spread the data
info_multiple_loci_effectors_upset <- info_multiple_loci_effectors %>%
  filter(Presence == TRUE) %>%  # Select only rows where effector is present
  select(Locus, Effector) %>%    # Select necessary columns
  mutate(Presence = 1) %>%      # Convert presence to 1
  group_by(Locus, Effector) %>% 
  summarise(Presence = max(Presence)) %>% # Get maximum presence value for each combination
  pivot_wider(names_from = Effector, values_from = Presence, values_fill = list(Presence = 0)) # Pivot data

#create column total in info_multiple_loci_effectors_upset
info_multiple_loci_effectors_upset$total <- rowSums(info_multiple_loci_effectors_upset[,2:ncol(info_multiple_loci_effectors_upset)])

#backup dataframe with all occurrences of effectors
info_multiple_loci_beforeDuplicateRemoval <- info_multiple_loci_effectors_upset

#remove total column
info_multiple_loci_effectors_upset <- subset(info_multiple_loci_effectors_upset, select = -c(total))
info_multiple_loci_beforeDuplicateRemoval <- subset(info_multiple_loci_beforeDuplicateRemoval, select = -c(total))

#convert tibble to df
info_multiple_loci_effectors_upset <- as.data.frame(info_multiple_loci_effectors_upset)
info_multiple_loci_beforeDuplicateRemoval <- as.data.frame(info_multiple_loci_beforeDuplicateRemoval)

# Convert your dataframe to list
info_multiple_loci_effectors_upset_list <- lapply(names(info_multiple_loci_effectors_upset)[-1], function(x){
  info_multiple_loci_effectors_upset$Locus[info_multiple_loci_effectors_upset[,x] == 1]})
info_multiple_loci_beforeDuplicateRemoval_list <- lapply(names(info_multiple_loci_beforeDuplicateRemoval)[-1], function(x){
  info_multiple_loci_beforeDuplicateRemoval$Locus[info_multiple_loci_beforeDuplicateRemoval[,x] == 1]})

# Make names of listInput as the names of Effectors
names(info_multiple_loci_effectors_upset_list) <- names(info_multiple_loci_effectors_upset)[-1]
names(info_multiple_loci_beforeDuplicateRemoval_list) <- names(info_multiple_loci_beforeDuplicateRemoval)[-1]

names(info_multiple_loci_effectors_upset)[names(info_multiple_loci_effectors_upset) == 'cam1'] <- 'Cam1'
names(info_multiple_loci_effectors_upset)[names(info_multiple_loci_effectors_upset) == 'cam2'] <- 'Cam2'
names(info_multiple_loci_effectors_upset)[names(info_multiple_loci_effectors_upset) == 'cam3'] <- 'Cam3'
names(info_multiple_loci_effectors_upset)[names(info_multiple_loci_effectors_upset) == 'can3'] <- 'Csm6-2'
names(info_multiple_loci_effectors_upset)[names(info_multiple_loci_effectors_upset) == 'csm6.ca6'] <- 'Csm6'
names(info_multiple_loci_effectors_upset)[names(info_multiple_loci_effectors_upset) == 'saved.chat'] <- 'SAVED-CHAT'
names(info_multiple_loci_effectors_upset)[names(info_multiple_loci_effectors_upset) == 'tirsaved'] <- 'TIR-SAVED'
names(info_multiple_loci_effectors_upset)[names(info_multiple_loci_effectors_upset) == 'csx1'] <- 'Csx1'
names(info_multiple_loci_effectors_upset)[names(info_multiple_loci_effectors_upset) == 'csx23'] <- 'Csx23'
names(info_multiple_loci_effectors_upset)[names(info_multiple_loci_effectors_upset) == 'calpl'] <- 'CalpL'
names(info_multiple_loci_effectors_upset)[names(info_multiple_loci_effectors_upset) == 'cami1'] <- 'Cami1'
names(info_multiple_loci_effectors_upset)[names(info_multiple_loci_effectors_upset) == 'can1.2'] <- 'Can1-2'
names(info_multiple_loci_effectors_upset)[names(info_multiple_loci_effectors_upset) == 'cora'] <- 'CorA'
names(info_multiple_loci_effectors_upset)[names(info_multiple_loci_effectors_upset) == 'nucc'] <- 'NucC'
names(info_multiple_loci_effectors_upset)[names(info_multiple_loci_effectors_upset) == 'csx16'] <- 'Csx16'
names(info_multiple_loci_effectors_upset)[names(info_multiple_loci_effectors_upset) == 'csx20'] <- 'Csx20'
names(info_multiple_loci_effectors_upset)[names(info_multiple_loci_effectors_upset) == 'crn1'] <- 'Crn1'
names(info_multiple_loci_effectors_upset)[names(info_multiple_loci_effectors_upset) == 'crn2'] <- 'Crn2'
names(info_multiple_loci_effectors_upset)[names(info_multiple_loci_effectors_upset) == 'crn3'] <- 'Crn3'
names(info_multiple_loci_effectors_upset)[names(info_multiple_loci_effectors_upset) == 'csx15'] <- 'Csx15'
names(info_multiple_loci_effectors_upset)[names(info_multiple_loci_effectors_upset) == 'blender'] <- 'Blender'

#same for info_multiple_loci_beforeDuplicateRemoval
names(info_multiple_loci_beforeDuplicateRemoval)[names(info_multiple_loci_beforeDuplicateRemoval) == 'cam1'] <- 'Cam1'
names(info_multiple_loci_beforeDuplicateRemoval)[names(info_multiple_loci_beforeDuplicateRemoval) == 'cam2'] <- 'Cam2'
names(info_multiple_loci_beforeDuplicateRemoval)[names(info_multiple_loci_beforeDuplicateRemoval) == 'cam3'] <- 'Cam3'
names(info_multiple_loci_beforeDuplicateRemoval)[names(info_multiple_loci_beforeDuplicateRemoval) == 'can3'] <- 'Can3'
names(info_multiple_loci_beforeDuplicateRemoval)[names(info_multiple_loci_beforeDuplicateRemoval) == 'csm6.ca6'] <- 'Csm6'
names(info_multiple_loci_beforeDuplicateRemoval)[names(info_multiple_loci_beforeDuplicateRemoval) == 'saved.chat'] <- 'SAVED-CHAT'
names(info_multiple_loci_beforeDuplicateRemoval)[names(info_multiple_loci_beforeDuplicateRemoval) == 'tirsaved'] <- 'TIR-SAVED'
names(info_multiple_loci_beforeDuplicateRemoval)[names(info_multiple_loci_beforeDuplicateRemoval) == 'csx1'] <- 'Csx1'
names(info_multiple_loci_beforeDuplicateRemoval)[names(info_multiple_loci_beforeDuplicateRemoval) == 'csx23'] <- 'Csx23'
names(info_multiple_loci_beforeDuplicateRemoval)[names(info_multiple_loci_beforeDuplicateRemoval) == 'calpl'] <- 'CalpL'
names(info_multiple_loci_beforeDuplicateRemoval)[names(info_multiple_loci_beforeDuplicateRemoval) == 'cami1'] <- 'Cami1'
names(info_multiple_loci_beforeDuplicateRemoval)[names(info_multiple_loci_beforeDuplicateRemoval) == 'can1.2'] <- 'Can1-2'
names(info_multiple_loci_beforeDuplicateRemoval)[names(info_multiple_loci_beforeDuplicateRemoval) == 'cora'] <- 'CorA'
names(info_multiple_loci_beforeDuplicateRemoval)[names(info_multiple_loci_beforeDuplicateRemoval) == 'nucc'] <- 'NucC'
names(info_multiple_loci_beforeDuplicateRemoval)[names(info_multiple_loci_beforeDuplicateRemoval) == 'csx16'] <- 'Csx16'
names(info_multiple_loci_beforeDuplicateRemoval)[names(info_multiple_loci_beforeDuplicateRemoval) == 'csx20'] <- 'Csx20'
names(info_multiple_loci_beforeDuplicateRemoval)[names(info_multiple_loci_beforeDuplicateRemoval) == 'crn1'] <- 'Crn1'
names(info_multiple_loci_beforeDuplicateRemoval)[names(info_multiple_loci_beforeDuplicateRemoval) == 'crn2'] <- 'Crn2'
names(info_multiple_loci_beforeDuplicateRemoval)[names(info_multiple_loci_beforeDuplicateRemoval) == 'crn3'] <- 'Crn3'
names(info_multiple_loci_beforeDuplicateRemoval)[names(info_multiple_loci_beforeDuplicateRemoval) == 'csx15'] <- 'Csx15'
names(info_multiple_loci_beforeDuplicateRemoval)[names(info_multiple_loci_beforeDuplicateRemoval) == 'blender'] <- 'Blender'

#merge the info_multiple_loci_effectors_upset and info dataframes
upset_merged_info <- merge(info_multiple_loci_effectors_upset, info, by.x = "Locus", by.y = "Locus")
upsetAll_merged_info <- merge(info_multiple_loci_beforeDuplicateRemoval, info, by.x = "Locus", by.y = "Locus")

# Signal molecule associations
#All
effector_coa_data = data.frame(
  set=names(info_multiple_loci_effectors_upset)[-1],
  cOAs=c('cA4', 'cA4', 'cA4', 'RN', 'RN', 'cA4', 'RN', 'cA3', 'SAM-AMP', 'cA6', 'RN', 'RN', 'cA6', 'cA4', 'RN', 'cA3', 'cA3', "BL", 'cA4', 'cA4', 'cA4')
)

# Archaea 
#effector_coa_data = data.frame(
#  set=names(info_multiple_loci_effectors_upset)[-1],
#  cOAs=c('cA4', 'cA4', 'cA4', 'cA4', 'SAM-AMP', 'cA4', 'cA4')
#)

#Thermophiles
#effector_coa_data = data.frame(
#  set=names(info_multiple_loci_effectors_upset)[-1],
#  cOAs=c('cA4', 'cA4', 'cA4', 'cA4', 'cA4', 'cA4', 'SAM-AMP')
#)

#png(paste("data/",project,"/plots/",project,"_upset_effRN.png", sep = ""), width = 12000, height = 8000, res = 600)
png(paste("data/",project,"/plots/",project,"_upset_effRN_archaea.png", sep = ""), width = 12000, height = 8000, res = 600)

#save as pdf
#pdf(paste("data/",project,"/plots/",project,"_upset_effRN.pdf", sep = ""), width = 10, height = 10)
upset(upset_merged_info[-1],
      names(info_multiple_loci_effectors_upset)[-1],
      min_degree=1,
      #min_size=3,
      #n_intersections=3,
      height_ratio=1,
      name = "Effector or combination",
      width_ratio = 0.15,
      set_sizes=(
        upset_set_size(
          geom=geom_bar(
            stat = "count",
            position='fill',
            aes(fill=has_multiple_effectors, x=group),
            width=0.8
          ),
          position='right'
        ) +
          scale_fill_manual(values = c("#D6D6D6", "#B3B3B3")) +
          theme(legend.position = "none")
      ),
      base_annotations=list(
        'Count'=intersection_size(
          bar_number_threshold=1,
          text=list(vjust=1.1),
          counts=FALSE,
          mapping=aes(fill=Subtype)) +
          scale_fill_manual(values = colorBlindBlack8) +
          ggtitle('Effector co-occurrences and subtype distributions') +
          annotate(
            geom='text', x=Inf, y=Inf,
            color='#595959',
            size=3,
            label=paste('Total occurrences:', nrow(info_multiple_loci_effectors_upset)),
            vjust=1, hjust=1
          ) +
          geom_text(stat='count', aes(group=intersection, label=..count.., y=..count..), nudge_y = 10, size = 3) + #use nudge_y = 3 for archaea
          theme(axis.title.y = element_text(vjust = -10),
                plot.margin = margin(0, 0, 0, 0, "cm"))
        
      ),
      stripes=upset_stripes(
        mapping=aes(color=cOAs),
        colors=c(
          'cA3'='#e0c6d8',
          'cA4'='#c1d6e8',
          'cA6'='#e8fae8',
          'SAM-AMP' = '#ffd9cc',
          'RN' = "darkgrey",
          "BL" = "white"
        ),
        data=effector_coa_data
      ),
      guides='over',
) +
  xlab('Number of loci') +
  ylab("Co-occurrence\nproportion (dark)") +
  theme(plot.title = element_text(size = 8),
        axis.text.x = element_blank())

dev.off()

#### 2. PHAGE RN ####

#convert RN columns in to logical
list_of_RNS <- c("crn1", "crn2", "crn3", "csx15", "csx16", "csx20")
#only keep the above columns plus sample in the df phage_rn
phage_rn <- subset(phage_rn, select = c("sample", list_of_RNS))
phage_rn <- phage_rn %>%
  mutate_at(vars(list_of_RNS), as.logical)

# create column boolean column for RN if any of the RNs are present
phage_rn$RN <- FALSE
phage_rn$RN[phage_rn$crn1 == TRUE | phage_rn$crn2 == TRUE | phage_rn$crn3 == TRUE | phage_rn$csx15 == TRUE | phage_rn$csx16 == TRUE | phage_rn$csx20 == TRUE] <- TRUE


#make phage_rn dataframe long
phage_rn_long <- phage_rn %>%
  pivot_longer(cols = c(crn1, crn2, crn3, csx15, csx16, csx20), names_to = "RN_name", values_to = "Presence")

phage_rn_long_true <- subset(phage_rn_long, Presence == TRUE)

#generate a table showing the presence of phage RNs
phage_rn_table <- phage_rn_long %>%
  group_by(RN_name, Presence) %>%
  summarise(count = n()) %>%
  pivot_wider(names_from = Presence, values_from = count, values_fill = list(count = 0))

#convert to dataframe
phage_rn_table <- as.data.frame(phage_rn_table)

#convert to long
phage_rn_table_long <- phage_rn_table %>%
  pivot_longer(cols = c("TRUE", "FALSE"), names_to = "Presence", values_to = "count")

#convert to dataframe
phage_rn_table_long <- as.data.frame(phage_rn_table_long)

phage_rn_table_long <- subset(phage_rn_table_long, Presence == TRUE)

#order by count
phage_rn_table_long <- phage_rn_table_long[order(phage_rn_table_long$count, decreasing = TRUE),]

#order for plotting
phage_rn_table_long$RN <- factor(phage_rn_table_long$RN, levels = phage_rn_table_long$RN)

#plot counts
phage_rn_plot <- ggplot(phage_rn_table_long, aes(x = RN, y = count)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme_minimal() +
  labs(title = "Phage RN presence", x = "RN", y = "Count") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  scale_fill_manual(values = c("#D6D6D6", "#B3B3B3"))

phage_rn_plot

#save as png
plot_path <- paste("data/",project,"/plots/phage_rn_",project,".png", sep = "")
ggsave(plot_path, plot = phage_rn_plot, width = 4, height = 4, bg = "white")

# investigate phages WITH RNs
phages_with_RNs <- subset(phage_rn, RN == TRUE)
nrow(phages_with_RNs)

#### 4. Fig 2: CAS10-cOA-RN TREE ####

RNs <- c("crn1", "crn2", "crn3", "csx15", "csx16", "csx20")
RN_colors <- c("#388dd1", "#1561a3", "#123675", "#091848", "#091848", "#091848")
offsets <- c(1, 2.5, 3.5, 4.5, 5.5, 6.5, 8)
custom_labels <- c("Crn1", "Crn2", "Crn3", "Csx15", "Csx16", "Csx20")
names(custom_labels) <- RNs
associated_colors_signals <- c("#CC79A7", "#1561a3", "#9ee493", "grey", "#ff743e")
HD_colors <- c("#4a4a4a", "#bfbfbf")
label_fontsize <- 2
offsets <- c(1, 2, 3, 4, 5, 6, 7, 8, 9)

rownames(info) <- info$Locus
cas10_tree$tip.label <- as.character(cas10_tree$tip.label)

#make the Locus column the first column in info
info <- info[, c('Locus', setdiff(names(info), 'Locus'))]

pX <- ggtree(cas10_tree, layout="circular", branch.length = "none", size=0.15, open.angle=20, aes(labels=genus)) +
               #geom_tiplab(size = 2, hjust = 0.5, vjust = 0.5, offset = 10) +
               theme(legend.position = "none")

pX <- pX %<+% info



pX <- gheatmap(pX, info[c("Subtype")], offset=-1.3, width=0.04, colnames = TRUE, colnames_angle=85, colnames_offset_y = 0, font.size = label_fontsize, color = "grey", hjust = 1, legend_title = "Subtype (inner ring)") + 
  scale_fill_manual(values = colorBlindBlack8, name = "Subtype\n(inner ring)") +
  theme(text = element_text(size = 14)) # + geom_point2(aes(subset = isTip == TRUE & domain == "Archaea")), color = "red", size = 0.1 #uncomment to label archaea

pX <- pX + geom_point2(aes(subset = (cyclase == FALSE)), color = "red", size = 0.1)

pX <- pX  + theme(legend.position=c(0.95, 0.5),
                  legend.title=element_text(size=7),
                  legend.text=element_text(size=6),
                  legend.spacing.y = unit(0.02, "cm"))



pX <- pX + new_scale_fill()
pX <- gheatmap(pX, info[c("locus_signal")], offset=0, width=lane_width, colnames = TRUE,
               colnames_angle=85, colnames_offset_y = 0, font.size = label_fontsize, color = "grey",
               hjust = 1, legend_title = "Signal", custom_column_labels = c("Signal")) +
              scale_fill_manual(values = associated_colors_signals, name = "Signal", labels = c("cA3", "cA4", "cA6", "None", "SAM-AMP"))

# Iterate over each ring nuclease
for(i in 1:length(RNs)){
  pX <- pX + new_scale_fill()
  #print current RN
  print(RNs[i])
  pX <- gheatmap(pX, info[c(RNs[i])], offset=offsets[i] + 0.3, width=lane_width, colnames = TRUE,
                 colnames_angle=85, colnames_offset_y = 0, colnames_offset_x = 0, font.size = label_fontsize, color = "grey",
                 hjust = 1, custom_column_labels = custom_labels[[RNs[i]]]) +
    scale_fill_manual(values=c("white", RN_colors[i]), guide = "none") + 
    theme(legend.position = 'none')
}


pX_opened <- open_tree(pX, 15)

pX_opened <- pX_opened  + theme(legend.position=c(0.95, 0.5),
                                legend.title=element_text(size=7),
                                legend.text=element_text(size=6),
                                legend.spacing.y = unit(0.02, "cm"))

#save pX
plot_path <- paste("data/",project,"/plots/cas10_tree_",project,".png", sep = "")
ggsave(plot_path, plot = pX_opened, width = 10, height = 10)
ggsave(plot_path, plot = pX_opened, width = 49, height = 49)

#### 5. RN OUTSIDE CRISPR ####

#remove any rows of rn_outside_crispr where diamond_target is crn5a or crn5b
rn_outside_crispr <- subset(rn_outside_crispr, diamond_target != "crn5a")
rn_outside_crispr <- subset(rn_outside_crispr, diamond_target != "crn5b")

rn_outside_crispr_filtered_outside <- rn_outside_crispr %>%
  filter(within_crispr == FALSE)

# compress using dplyr to get rid of multiple lines per RN
rn_df_compressed <- rn_outside_crispr %>%
  group_by(diamond_query) %>%
  summarize(
    within_crispr_all = any(within_crispr),
    within_prophage_all = any(within_prophage),
    across(everything(), first, .names = "first_{.col}") # get the first value of each column
  )

#filter by evalue
e_cutoff <- 1e-50
rn_df_compressed <- rn_df_compressed %>%
  filter(first_evalue < e_cutoff)

# convert the tibble into a df that shows shows how many RNs are within crispr, within prophage, both or neither
rn_plot <- ggplot(rn_df_compressed, aes(x = within_crispr_all, fill = within_prophage_all)) +
  geom_bar(stat = "count") +
  labs(title = "RN outside CRISPR", x = "Within CRISPR", y = "Count")
rn_plot

#make a table that shows the average number of ring nucleases per host (column first_host shows host)
rn_table <- rn_df_compressed %>%
  group_by(first_host) %>%
  summarize(count = n()) %>%
  arrange(desc(count))
rn_table

average_rn_per_host <- mean(rn_table$count)
average_rn_per_host

#summarise by the ring nuclease type (column first_diamond_target)
rn_summarised_by_rn <- rn_df_compressed %>%
  group_by(first_diamond_target) %>%
  summarize(count = n()) %>%
  arrange(desc(count))

#plot 
rn_plot_by_rn <- ggplot(rn_summarised_by_rn, aes(x = reorder(first_diamond_target, -count), y = count)) +
  geom_bar(stat = "identity") +
  labs(title = "RN outside CRISPR", x = "RN", y = "Count") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
rn_plot_by_rn

#### 6. NUMBERS ####
# Reporting numbers and statistics. Note that these numbers are based on filtering loci with Cas10 length above 500 aa.

#number of loci in total
total_locus_count <- nrow(info)

#loci with rn
locus_with_rn <- subset(info, has_ring_nuclease == TRUE)
locus_with_rn_count <- nrow(locus_with_rn)

#loci with rn that have more than one rn
locus_with_multiple_rn <- subset(locus_with_rn, ring_nuclease_count > 1)
locus_with_multiple_rn_count <- nrow(locus_with_multiple_rn)

#number of loci with RN
total_locus_count_rn <- nrow(subset(info, has_ring_nuclease == TRUE))

#percentage of loci with rn
percentage_locus_count_rn <- total_locus_count_rn / total_locus_count * 100

#number of loci with crn1
locus_with_crn1 <- subset(info, crn1 == TRUE)
nrow(locus_with_crn1)

#number of loci with crn2
locus_with_crn2 <- subset(info, crn2 == TRUE)
nrow(locus_with_crn2)

#number of loci with crn3
locus_with_crn3 <- subset(info, crn3 == TRUE)
nrow(locus_with_crn3)

#number of loci with csx15
locus_with_csx15 <- subset(info, csx15 == TRUE)
nrow(locus_with_csx15)

#number of loci with csx16
locus_with_csx16 <- subset(info, csx16 == TRUE)
nrow(locus_with_csx16)

#number of loci with csx20
locus_with_csx20 <- subset(info, csx20 == TRUE)
nrow(locus_with_csx20)

#number of loci with Crn1 and Csx15
locus_with_crn1_and_csx15 <- subset(info, crn1 == TRUE & csx15 == TRUE)
nrow(locus_with_crn1_and_csx15)

#create a table that shows counts for all RNs
rn_table <- data.frame(RN = c("Crn1", "Crn2", "Crn3", "Csx15", "Csx16", "Csx20"), Count = c(nrow(locus_with_crn1), nrow(locus_with_crn2), nrow(locus_with_crn3), nrow(locus_with_csx15), nrow(locus_with_csx16), nrow(locus_with_csx20)))

#create dataframe that contains percentages for having and not having RN
locus_count_df <- data.frame(Feature = c("RN", "No RN"), Count = c(total_locus_count_rn, total_locus_count - total_locus_count_rn), Percentage = c(percentage_locus_count_rn, 100 - percentage_locus_count_rn))

#loci with RN but no effector
locus_with_rn_no_effector <- subset(info, has_ring_nuclease == TRUE & effector_count_known_new_sum == 0)
locus_with_rn_no_effector_count <- nrow(locus_with_rn_no_effector)

#percentage of loci with RN but no effector
percentage_locus_with_rn_no_effector <- locus_with_rn_no_effector_count / total_locus_count_rn * 100

#loci that have an RN but no effector, but have an effector elsewhere in the genome
locus_with_rn_no_effector_but_effector_elsewhere <- subset(info, has_ring_nuclease == TRUE & effector_count_known_new_sum == 0 & effector_elsewhere_in_sample == TRUE)
nrow(locus_with_rn_no_effector_but_effector_elsewhere)

#percentage of loci that have an RN but no effector, but have an effector elsewhere in the genome
percentage_locus_with_rn_no_effector_but_effector_elsewhere <- nrow(locus_with_rn_no_effector_but_effector_elsewhere) / total_locus_count_rn * 100

locus_with_rn_no_effector_but_effector_elsewhere <- subset(info, has_ring_nuclease == TRUE & effector_count_known_new_sum == 0 & effector_elsewhere_in_sample == TRUE)

percentage_locus_with_rn_no_effector_elsewhere <- nrow(locus_with_rn_no_effector_but_effector_elsewhere) / total_locus_count_rn * 100

#loci with effectors
locus_with_effector <- subset(info, effector_count_known_new_sum > 0)
nrow(locus_with_effector)

#loci with effectors that have also RN
locus_with_effector_and_rn <- subset(locus_with_effector, has_ring_nuclease == TRUE)
nrow(locus_with_effector_and_rn)

#how many loci with effectors and rns have an effector elsewhere in the genome
locus_with_effector_and_rn_elsewhere <- subset(locus_with_effector_and_rn, effector_elsewhere_in_sample == TRUE)
nrow(locus_with_effector_and_rn_elsewhere)

#numbers and percentages of loci that have effectors and associated RNs
locus_with_effector_count <- nrow(locus_with_effector)
locus_with_effector_and_rn_count <- nrow(locus_with_effector_and_rn)
percentage_locus_with_effector_and_rn <- locus_with_effector_and_rn_count / locus_with_effector_count * 100

#loci with cA4 effectors
locus_with_cA4 <- subset(info, ca4 == TRUE)
locus_with_ca4_count <- nrow(locus_with_cA4)

#loci with cA4 effectors that have also RN
locus_with_cA4_and_rn <- subset(locus_with_cA4, has_ring_nuclease == TRUE)
nrow(locus_with_cA4_and_rn)
#numbers and percentages of loci that have cA4 effectors and associated RNs
locus_with_cA4_and_rn_count <- nrow(locus_with_cA4_and_rn)
percentage_locus_with_cA4_and_rn <- locus_with_cA4_and_rn_count / locus_with_ca4_count * 100


# create pie chart of locus_count_df
A_locus_count_pie <- ggplot(locus_count_df, aes(x = "", y = Percentage, fill = Feature)) +
  geom_bar(stat = "identity", width = 1) +
  coord_polar("y", start = 0) +
  theme_void() +
  theme(legend.position = "bottom") +
  #use a modern color palette with red and grey
  scale_fill_manual(values = c("#e74c3c", "#D6D6D6")) +
  geom_text(aes(label = paste0(round(Percentage, 1), "%")), position = position_stack(vjust = 0.5))

A_locus_count_pie

#create a similar pie with the same colors but with the percentage of loci with RN but no effector
locus_count_df_rn <- data.frame(Feature = c("RN without effector", "RN with effector"), Count = c(locus_with_rn_no_effector_count, total_locus_count_rn - locus_with_rn_no_effector_count), Percentage = c(percentage_locus_with_rn_no_effector, 100 - percentage_locus_with_rn_no_effector))

B_locus_count_pie <- ggplot(locus_count_df_rn, aes(x = "", y = Percentage, fill = Feature)) +
  geom_bar(stat = "identity", width = 1) +
  coord_polar("y", start = 0) +
  theme_void() +
  theme(legend.position = "bottom") +
  #use a modern color palette with red and grey
  scale_fill_manual(values = c("#e74c3c", "#D6D6D6")) +
  geom_text(aes(label = paste0(round(Percentage, 1), "%")), position = position_stack(vjust = 0.5))

B_locus_count_pie

#create a similar pie with the same colors but with the percentage of loci with RN but no effector, but have an effector elsewhere in the genome
locus_count_df_rn_elsewhere <- data.frame(Feature = c("RN without effector", "RN with effector elsewhere"), Count = c(locus_with_rn_no_effector_count - nrow(locus_with_rn_no_effector_but_effector_elsewhere), nrow(locus_with_rn_no_effector_but_effector_elsewhere)), Percentage = c(percentage_locus_with_rn_no_effector - percentage_locus_with_rn_no_effector_elsewhere, percentage_locus_with_rn_no_effector_elsewhere))

C_locus_count_pie <- ggplot(locus_count_df_rn_elsewhere, aes(x = "", y = Percentage, fill = Feature)) +
  geom_bar(stat = "identity", width = 1) +
  coord_polar("y", start = 0) +
  theme_void() +
  theme(legend.position = "bottom") +
  #use a modern color palette with red and grey
  scale_fill_manual(values = c("#e74c3c", "#D6D6D6")) +
  geom_text(aes(label = paste0(round(Percentage, 1), "%")), position = position_stack(vjust = 0.5))

C_locus_count_pie

#Number of loci with effector but no RN
locus_with_effector_no_rn <- subset(info, has_ring_nuclease == FALSE & effector_count_known_new_sum > 0)
locus_with_effector_no_rn_count <- nrow(locus_with_effector_no_rn)

crn1_count <- sum(info$crn1)
crn2_count <- sum(info$crn2)
crn3_count <- sum(info$crn3)
csx15_count <- sum(info$csx15)
csx16_count <- sum(info$csx16)
csx20_count <- sum(info$csx20)

#create a results table with the counts of each RN
rn_count_table <- data.frame(RN = c("Crn1", "Crn2", "Crn3", "Csx15", "Csx16", "Csx20"),
                            Count = c(crn1_count, crn2_count, crn3_count, csx15_count, csx16_count, csx20_count))
rn_count_table

#subset only loci with crn1
crn1_loci <- subset(info, crn1 == TRUE)

#get loci with csx1 and no ring nucleases
csx1_no_rn <- subset(info, has_ring_nuclease == FALSE & csx1 == TRUE)
csx1_no_rn_count <- nrow(csx1_no_rn)
csx1_no_rn_count

#get loci with csx1 and ring nucleases
csx1_rn <- subset(info, has_ring_nuclease == TRUE & csx1 == TRUE)
csx1_rn_count <- nrow(csx1_rn)
csx1_rn_count

#calculate total csx1 loci and calculate the proportion that has no ring nuclease
csx1_total <- sum(info$csx1)
csx1_no_rn_proportion <- csx1_no_rn_count / csx1_total
csx1_no_rn_proportion




#what proportion of effector instances are associated with a ring nuclease
info_only_singular_effectors <- subset(info, effector_count_known_new_sum == 1)
calpl <- subset(info_only_singular_effectors, calpl == TRUE)
can12 <- subset(info_only_singular_effectors, can1.2 == TRUE)

nrow(subset(can12, has_ring_nuclease == TRUE))
nrow(subset(can12, has_ring_nuclease == FALSE))

rn_proportion <- function(rn){
  rn_no_rn <- subset(info_only_singular_effectors, get(rn) == TRUE & has_ring_nuclease == TRUE)
  rn_no_rn_count <- nrow(rn_no_rn)
  
  rn_total <- sum(info_only_singular_effectors[,rn])
  rn_no_rn_proportion <- rn_no_rn_count / rn_total
  return(rn_no_rn_proportion)
}

# run the function for all effectors ("can3", "tirsaved", "cam1", "cam2", "cam3", "saved.chat", "nucc", "calpl", "cami1", "can1.2", "csx1", "csx23", "csm6.ca6", "cora")
effectors <- c("can3", "tirsaved", "cam1", "cam2", "cam3", "saved.chat", "nucc", "calpl", "cami1", "can1.2", "csx1", "csx23", "csm6.ca6", "cora")
effectors_only_singular_with_RNs <- c("cam1", "cami1", "can1.2", "csx1")
effectors_rn_proportions <- sapply(effectors, rn_proportion)
effectors_rn_proportions

#remove effectors with value 1 or NaN
effectors_rn_proportions <- effectors_rn_proportions[!is.nan(effectors_rn_proportions)]
effectors_rn_proportions <- effectors_rn_proportions[effectors_rn_proportions != 0]
effectors_rn_proportions

#get total counts for each effector
effectors_total_counts <- sapply(effectors_only_singular_with_RNs, function(x) sum(info_only_singular_effectors[,x]))

#bundle effector counts and proportions into a dataframe and leave out nonmatching effectors
effectors_table <- data.frame(Effectors = effectors_only_singular_with_RNs, Total = effectors_total_counts, RN_proportion = effectors_no_rn_proportions)

#create column with_rn
effectors_table$with_rn <- effectors_table$Total - effectors_table$Total * effectors_table$RN_proportion

#create column that shows with_rn/total as string
effectors_table$with_rn_string <- paste(effectors_table$with_rn, "/", effectors_table$Total, sep = "")

effectors_table_plot <- ggplot(effectors_table, aes(y = reorder(Effectors, RN_proportion), x = RN_proportion)) +
  geom_lollipop(horizontal = TRUE, point.colour = "#e74c3c", point.size = 4, line.size = 1.2, line.colour = "#2c3e50") +
  #also display counts next to lollipop
  #geom_text(aes(label = with_rn_string), hjust = -0.6, size = 3) +
  labs(
    #title = "Proportion of effectors with no RN",
    x = "Associated RN proportion",
    y = "Effector"
  ) +
  scale_x_continuous(limits = c(0, 1)) +  # Remove padding around x-axis
  theme_minimal(base_size = 14) +  # Use a minimal theme for a clean look
  theme(
    axis.text.y = element_text(color = "#34495e", hjust = 1),  # Set y-axis text color
    axis.text.x = element_text(color = "#34495e"),  # Set x-axis text color
    #plot.title = element_text(hjust = 0.5, face = "bold", color = "#2c3e50", size = 12),  # Center and bold the title
    axis.title.x = element_text(margin = margin(t = 10)),  # Add margin to x-axis title
    axis.title.y = element_text(margin = margin(r = 10)),  # Add margin to y-axis title
    panel.grid.major.y = element_blank(),  # Remove horizontal grid lines
    panel.grid.minor = element_blank()  # Remove minor grid lines
  )

effectors_table_plot

#save path
plot_path <- paste("data/",project,"/plots/effectors_no_rn_proportions_",project,".png", sep = "")
ggsave(plot_path, plot = effectors_table_plot, width = 4, height = 4, units = "in", dpi = 300, bg = "white")

#overall, how many loci do not encode a ring nuclease but have an effector
no_rn_has_effector <- subset(info, has_ring_nuclease == FALSE & has_effector == TRUE)
no_rn_has_effector_count <- nrow(no_rn_has_effector)
#proportion
no_rn_has_effector_proportion <- no_rn_has_effector_count / nrow(info)
no_rn_has_effector_proportion

#how many loci have no effector (has_known_effector == FALSE) and have a ring nuclease
no_effector_has_rn <- subset(info, has_known_effector == FALSE & has_ring_nuclease == TRUE)
no_effector_has_rn_count <- nrow(no_effector_has_rn)
no_effector_has_rn_count

#how often such effector_elsewhere_in_sample is True
no_effector_has_rn_and_effector_elsewhere <- subset(no_effector_has_rn, effector_elsewhere_in_sample_but_not_here == TRUE)
no_effector_has_rn_and_effector_elsewhere_count <- nrow(no_effector_has_rn_and_effector_elsewhere)
no_effector_has_rn_and_effector_elsewhere_count

#contrast this with same counts from loci with effector
effector_has_rn <- subset(info, has_known_effector == TRUE & has_ring_nuclease == TRUE)
effector_has_rn_count <- nrow(effector_has_rn)
effector_has_rn_count

#get all effector_has_rn that have multilocus sample
effector_has_rn_multilocus <- subset(effector_has_rn, multilocus_sample == "True")

#find entries info whose Sample column matches the Sample column in effector_has_rn_multilocus
effector_has_rn_multilocus_info <- subset(info, Sample %in% effector_has_rn_multilocus$Sample)

#from effector_has_rn_multilocus_info remove all entries whose Locus matches those in effector_has_rn_multilocus
effector_has_rn_multilocus_info_no_effector <- subset(effector_has_rn_multilocus_info, !(Locus %in% effector_has_rn_multilocus$Locus))

#then remove entries where there is no effector present (effector_count_known_new_sum == 0)
effector_has_rn_multilocus_info_no_effector_no_effector <- subset(effector_has_rn_multilocus_info_no_effector, effector_count_known_new_sum == 0)
nrow(effector_has_rn_multilocus_info_no_effector_no_effector)

effector_has_rn_and_effector_elsewhere <- subset(effector_has_rn, effector_elsewhere_in_sample_but_not_here == TRUE)
effector_has_rn_and_effector_elsewhere_count <- nrow(effector_has_rn_and_effector_elsewhere)
effector_has_rn_and_effector_elsewhere_count

#investigate how often RNs co-occur in same locus. Create new column RN_cooccurrence, which is set to TRUE if multiple RNs are found in locus
info$RN_cooccurrence <- FALSE
#calculate the number of RNs in each locus
info$RN_count <- rowSums(info[,c("crn1", "crn2", "crn3", "csx15", "csx16", "csx20")])
#set RN_cooccurrence to TRUE if RN_count is greater than 1
info$RN_cooccurrence[info$RN_count > 1] <- TRUE
rn_cooccurrence_count <- sum(info$RN_cooccurrence)

loci_with_cooccurring_RNs <- subset(info, RN_cooccurrence == TRUE)

#how often csx15 is the only RN
csx15_only <- subset(info, csx15 == TRUE & RN_count == 1)
csx15_only_count <- nrow(csx15_only)

#and how often csx15 is not the only RN
csx15_not_only <- subset(info, csx15 == TRUE & RN_count > 1)
csx15_not_only_count <- nrow(csx15_not_only)

#### Csx15 outside CRISPR ####

csx15_hmm_hits <- read.csv(paste("data/",project,"/csx15_hits.tsv", sep = ""), sep = "\t")

#keep only unique sample from column csx15_sample
csx15_hmm_hits <- csx15_hmm_hits %>%
  group_by(csx15_sample) %>%
  filter(row_number() == 1)
nrow(csx15_hmm_hits)

# merge csx15_hmm_hits with info using csx15_hmm_hits on left and info on right. Merge by column csx15_sample in csx15_hmm_hits and sample in info. Do left join
csx15_hmm_hits_merged_info <- merge(csx15_hmm_hits, info, by.x = "csx15_sample", by.y = "sample", all.x = TRUE)
nrow(csx15_hmm_hits_merged_info)
str(csx15_hmm_hits_merged_info)

#remove duplicates based on column csx15_sample.
csx15_hmm_hits_merged_info_no_duplicates <- csx15_hmm_hits_merged_info[!duplicated(csx15_hmm_hits_merged_info$csx15_sample),]
nrow(csx15_hmm_hits_merged_info_no_duplicates)

# The df now contains all csx15 samples and if one or more crispr loci were found, the column Locus will not be empty.
# We disregard the subtype of the locus, it's enough to know there is some kind of a type iii locus present.

#identify loci with csx15 and no crispr by finding any rows where Locus is Nan
genomes_with_csx15_no_crispr <- subset(csx15_hmm_hits_merged_info, is.na(Locus))
nrow(genomes_with_csx15_no_crispr)

#identify loci with csx15 and crispr by finding any rows where Locus is not Nan
genomes_with_csx15_crispr <- subset(csx15_hmm_hits_merged_info, !is.na(Locus))
nrow(genomes_with_csx15_crispr)

#bundle to a dataframe
csx15_crispr_table <- data.frame(Genomes_with_Csx15 = c("With CRISPR", "Without CRISPR"), Count = c(nrow(genomes_with_csx15_crispr), nrow(genomes_with_csx15_no_crispr)))
#add percentages and round to 2 decimal places
csx15_crispr_table$Percentage <- round((csx15_crispr_table$Count / sum(csx15_crispr_table$Count)) * 100, 2)
#add total row 
csx15_crispr_table <- rbind(csx15_crispr_table, c("Total", sum(csx15_crispr_table$Count), 100))
csx15_crispr_table               

#use a logaritmic scale for the x axis as the values vary from 0.0001 to 1e-100
csx15_evalue_plot <- ggplot(subset(csx15_hmm_hits_merged_info_no_duplicates, evalue < 1e-4), aes(x = evalue)) +
  geom_histogram(bins = 50) +
  scale_x_log10() +
  labs(title = "E-values of csx15 hits", x = "E-value", y = "Count")
csx15_evalue_plot

#subest the data on evalue 1e-20 and lower and recreate the final df
csx15_hmm_hits_merged_info_no_duplicates_filtered <- subset(csx15_hmm_hits_merged_info_no_duplicates, evalue < 1e-20)
nrow(csx15_hmm_hits_merged_info_no_duplicates_filtered)

#bundle to a dataframe
csx15_crispr_table_filtered <- data.frame(Genomes_with_Csx15 = c("With CRISPR", "Without CRISPR"), Count = c(nrow(subset(csx15_hmm_hits_merged_info_no_duplicates_filtered, !is.na(Locus))), nrow(subset(csx15_hmm_hits_merged_info_no_duplicates_filtered, is.na(Locus)))))
#add percentages and round to 2 decimal places
csx15_crispr_table_filtered$Percentage <- round((csx15_crispr_table_filtered$Count / sum(csx15_crispr_table_filtered$Count)) * 100, 2)
#add total row
csx15_crispr_table_filtered <- rbind(csx15_crispr_table_filtered, c("Total", sum(csx15_crispr_table_filtered$Count), 100))
csx15_crispr_table_filtered

#### Other ####

#Extract all type III-D loci
type_iii_d <- subset(info, Subtype == "III-D")
nrow(type_iii_d)

#create a .txt file of the column Locus with each locus on a new line in current project
write.table(type_iii_d$sample, file = paste("data/",project,"/type_iii_d_loci.txt", sep = ""), quote = FALSE, row.names = FALSE, col.names = FALSE)

