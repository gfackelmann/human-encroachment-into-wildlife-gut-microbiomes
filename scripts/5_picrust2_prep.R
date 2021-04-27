#Copyright (c) 2021 Insitute of Evolutionary Ecology and Conservation Genomics, Ulm University, Germany. See the LICENSE file at the top-level directory of this distribution (https://github.com/gfackelmann/human-encroachment-into-wildlife-gut-microbiomes/blob/main/LICENSE). 
#Author: Gloria Fackelmann (gloria.fackelmann@uni-ulm.de; https://www.researchgate.net/profile/Gloria-Fackelmann)

##############################################################################################################
#Load necessary libraries
##############################################################################################################

#R version 3.6.1
library("phyloseq") #phyloseq_1.28.0
library("dplyr") #dplyr_0.8.3
library("biomformat") #biomformat_1.12.0
library("Biostrings") #Biostrings_2.52.0
library("dada2") #dada2_1.12.1

##############################################################################################################
#Export ASV table of entire dataset as a biom file to generate unstrat table using PICRUSt2
##############################################################################################################
otu_biom_whole_dataset <- make_biom(data=as(otu_table(proec_pellet_depth_decontam),"matrix"))
write_biom(otu_biom_whole_dataset, "otu_biom_whole_dataset.biom")

##############################################################################################################
#Export representative sequences file of entire dataset as a .fna file to use in PICRUSt2
##############################################################################################################
#read the full rep seq file (downloaded and saved as fasta file from qiime2 visualization of the qza file)
sequences_fasta <- readDNAStringSet("sequences.fasta")

#arrange in a dataframe:
sequences_fasta_ASV <- names(sequences_fasta)
sequences_fasta_sequences <- paste(sequences_fasta)
sequences_fasta_df <- data.frame(sequences_fasta_ASV, sequences_fasta_sequences)

#build vector of taxa names from phyloseq object and use this to subset fasta sequences
keep_taxa_picrust <- taxa_names(proec_pellet_depth_decontam)
sequences_fasta_df_subset <- sequences_fasta_df[with(sequences_fasta_df, sequences_fasta_ASV %in% keep_taxa_picrust),]

#rename the dataframe columns to export using dada2
sequences_fasta_df_subset_names <- rename(sequences_fasta_df_subset, abundance = sequences_fasta_ASV, sequence = sequences_fasta_sequences)
uniquesToFasta(sequences_fasta_df_subset_names, fout = "rep-seqs.fna", ids = sequences_fasta_df_subset_names$abundance)

##############################################################################################################
#Divide dataset into 4 by samples (not by ASVs) in order to generate strat table using PICRUSt2
##############################################################################################################
proec_pellet_depth_decontam_subsets <- proec_pellet_depth_decontam
sample_data(proec_pellet_depth_decontam_subsets)$picrust_subsets <- c("")
sample_data(proec_pellet_depth_decontam_subsets)$picrust_subsets[1:96] <- c("group1")
sample_data(proec_pellet_depth_decontam_subsets)$picrust_subsets[97:192] <- c("group2")
sample_data(proec_pellet_depth_decontam_subsets)$picrust_subsets[193:288] <- c("group3")
sample_data(proec_pellet_depth_decontam_subsets)$picrust_subsets[289:384] <- c("group4")

proec_pellet_depth_decontam_picrust_part1 <- proec_pellet_depth_decontam_subsets %>%
  subset_samples(picrust_subsets %in% ("group1"))%>%
  prune_taxa(taxa_sums(.) > 0, .)

proec_pellet_depth_decontam_picrust_part2 <- proec_pellet_depth_decontam_subsets %>%
  subset_samples(picrust_subsets %in% ("group2"))%>%
  prune_taxa(taxa_sums(.) > 0, .)

proec_pellet_depth_decontam_picrust_part3 <- proec_pellet_depth_decontam_subsets %>%
  subset_samples(picrust_subsets %in% ("group3"))%>%
  prune_taxa(taxa_sums(.) > 0, .)

proec_pellet_depth_decontam_picrust_part4 <- proec_pellet_depth_decontam_subsets %>%
  subset_samples(picrust_subsets %in% ("group4"))%>%
  prune_taxa(taxa_sums(.) > 0, .)

##############################################################################################################
#Export ASVs tables of each subset as a biom file
##############################################################################################################
otu_biom_whole_dataset_part1 <- make_biom(data=as(otu_table(proec_pellet_depth_decontam_picrust_part1),"matrix"))
write_biom(otu_biom_whole_dataset_part1, "otu_biom_whole_dataset_part1.biom")

otu_biom_whole_dataset_part2 <- make_biom(data=as(otu_table(proec_pellet_depth_decontam_picrust_part2),"matrix"))
write_biom(otu_biom_whole_dataset_part2, "otu_biom_whole_dataset_part2.biom")

otu_biom_whole_dataset_part3 <- make_biom(data=as(otu_table(proec_pellet_depth_decontam_picrust_part3),"matrix"))
write_biom(otu_biom_whole_dataset_part3, "otu_biom_whole_dataset_part3.biom")

otu_biom_whole_dataset_part4 <- make_biom(data=as(otu_table(proec_pellet_depth_decontam_picrust_part4),"matrix"))
write_biom(otu_biom_whole_dataset_part4, "otu_biom_whole_dataset_part4.biom")

#Now use the files “otu_biom_whole_dataset.biom”, otu_biom_whole_dataset_parts 1-4, and “rep-seqs.fna” to run PICRUSt2 in the terminal.