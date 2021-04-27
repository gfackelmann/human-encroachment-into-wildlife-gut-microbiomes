#Copyright (c) 2021 Insitute of Evolutionary Ecology and Conservation Genomics, Ulm University, Germany. See the LICENSE file at the top-level directory of this distribution (https://github.com/gfackelmann/human-encroachment-into-wildlife-gut-microbiomes/blob/main/LICENSE). 
#Author: Gloria Fackelmann (gloria.fackelmann@uni-ulm.de; https://www.researchgate.net/profile/Gloria-Fackelmann)

##############################################################################################################
#Load necessary libraries
##############################################################################################################

#R version 3.6.1
library("phyloseq") #phyloseq_1.28.0
library("dplyr") #dplyr_0.8.3
library("psych") #psych_1.8.12
library("vegan") #vegan_2.5-5
source("/utils/dredge_adonis.R")
source("/utils/AICc_PERMANOVA_DREDGE.R")

##############################################################################################################
#Calculate beta diversity distances
##############################################################################################################
#phyloseq object is 'proec_beta_decontam'.

#weighted unifrac
DistW <- phyloseq::distance(proec_beta_decontam, method="wUniFrac")
#unweighted unifrac
DistUW <- phyloseq::distance(proec_beta_decontam, method="uunifrac")

##############################################################################################################
#PERMANOVA: IT model selection (Supplementary Data 6-7 & 17-18)
##############################################################################################################

SequencingDepth <- scale(as(sample_data(proec_beta_decontam), "data.frame")$SequencingDepth)
Site <- as(sample_data(proec_beta_decontam), "data.frame")$capture_site
landscape <- as(sample_data(proec_beta_decontam), "data.frame")$landscape
density <- as(sample_data(proec_beta_decontam), "data.frame")$density
extractionBatch <- as(sample_data(proec_beta_decontam), "data.frame")$extractionBatch
sex <- as(sample_data(proec_beta_decontam), "data.frame")$sex
season <- as(sample_data(proec_beta_decontam), "data.frame")$season

#weighted UniFrac without density
output_dredge_adonis_W_without_density <- dredge_adonis("DistW", c("SequencingDepth", "landscape", "landscape / Site", "extractionBatch", "season", "sex"), list(c("landscape / Site", "landscape")),"strata = Site", "proec_beta_decontam")

#weighted UniFrac with density
output_dredge_adonis_W_with_density <- dredge_adonis("DistW", c("SequencingDepth", "landscape", "landscape / Site", "extractionBatch", "season", "sex", "density"), list(c("landscape / Site", "landscape")),"strata = Site", "proec_beta_decontam")

#unweighted UniFrac without density
output_dredge_adonis_UW_without_density <- dredge_adonis("DistUW", c("SequencingDepth", "landscape", "landscape / Site", "extractionBatch", "season", "sex"), list(c("landscape / Site", "landscape")),"strata = Site", "proec_beta_decontam")

#unweighted UniFrac with density
output_dredge_adonis_UW_with_density <- dredge_adonis("DistUW", c("SequencingDepth", "landscape", "landscape / Site", "extractionBatch", "season", "sex", "density"), list(c("landscape / Site", "landscape")),"strata = Site", "proec_beta_decontam")

##############################################################################################################
#effect sizes:

##############################################################################################################
#weighted UniFrac
ordW <- ordinate(proec_beta_decontam, method = "PCoA", distance = DistW)

pcoa_axis_1 <- ordW$vectors[,'Axis.1']
pcoa_axis_2 <- ordW$vectors[,'Axis.2']

pcoa_axis_1_df <- stack(pcoa_axis_1)
pcoa_axis_2_df <- stack(pcoa_axis_2)

names(pcoa_axis_1_df) <- c("PCoA1", "sampleID")
names(pcoa_axis_2_df) <- c("PCoA2", "sampleID")

pcoa_axis_1_2_df <- merge(pcoa_axis_1_df, pcoa_axis_2_df, by="sampleID")

pcoa_axis_landscape_df <- data.frame(
                sampleID = proec_alpha_decontam$X.SampleID, 
                landscape = proec_alpha_decontam$landscape)

pcoa_axis_1_2_landscape_df <- merge(pcoa_axis_1_2_df, pcoa_axis_landscape_df, by="sampleID")

pcoa_coords_weighted_FF <- pcoa_axis_1_2_landscape_df %>% filter(landscape == "forest.fragment")
pcoa_coords_weighted_CF <- pcoa_axis_1_2_landscape_df %>% filter(landscape == "continuous.forest")
pcoa_coords_weighted_I <- pcoa_axis_1_2_landscape_df %>% filter(landscape == "island")

vector_CF_PCoA1 <- pcoa_coords_weighted_CF$PCoA1
vector_FF_PCoA1 <- pcoa_coords_weighted_FF$PCoA1
vector_I_PCoA1 <- pcoa_coords_weighted_I$PCoA1

vector_CF_PCoA2 <- pcoa_coords_weighted_CF$PCoA2
vector_FF_PCoA2 <- pcoa_coords_weighted_FF$PCoA2
vector_I_PCoA2 <- pcoa_coords_weighted_I$PCoA2

d.raw.data<-function(E.data,C.data){  
  d<-(mean(E.data)-mean(C.data))/
  sqrt(((length(E.data)-1)*var(E.data)+
  (length(C.data)-1)*var(C.data))/
  (length(E.data)+length(C.data)-2))
  names(d)<-"effect size d"
  return(d)
}

##############################################################################################################
#PCoA1
#A-C
effect_size <- d.raw.data(vector_CF_PCoA1,vector_FF_PCoA1)
psych::cohen.d.ci(d = effect_size, n1 = length(vector_FF_PCoA1), n2 = length(vector_CF_PCoA1))

#A-I
effect_size <- d.raw.data(vector_I_PCoA1,vector_FF_PCoA1)
psych::cohen.d.ci(d = effect_size, n1 = length(vector_FF_PCoA1), n2 = length(vector_I_PCoA1))

#C-I
effect_size <- d.raw.data(vector_I_PCoA1, vector_CF_PCoA1)
psych::cohen.d.ci(d = effect_size, n1 = length(vector_CF_PCoA1), n2 = length(vector_I_PCoA1))

##############################################################################################################
#PCoA2
#A-C
effect_size <- d.raw.data(vector_CF_PCoA2,vector_FF_PCoA2)
psych::cohen.d.ci(d = effect_size, n1 = length(vector_FF_PCoA2), n2 = length(vector_CF_PCoA2))

#A-I
effect_size <- d.raw.data(vector_I_PCoA2,vector_FF_PCoA2)
psych::cohen.d.ci(d = effect_size, n1 = length(vector_FF_PCoA2), n2 = length(vector_I_PCoA2))

#C-I
effect_size <- d.raw.data(vector_I_PCoA2, vector_CF_PCoA2)
psych::cohen.d.ci(d = effect_size, n1 = length(vector_CF_PCoA2), n2 = length(vector_I_PCoA2))

##############################################################################################################
#Unweighted UniFrac
ordUW <- ordinate(proec_beta_decontam, method = "PCoA", distance = DistUW)

pcoa_axis_1_uw <- ordUW$vectors[,'Axis.1']
pcoa_axis_2_uw <- ordUW$vectors[,'Axis.2']

pcoa_axis_1_uw_df <- stack(pcoa_axis_1_uw)
pcoa_axis_2_uw_df <- stack(pcoa_axis_2_uw)

names(pcoa_axis_1_uw_df) <- c("PCoA1", "sampleID")
names(pcoa_axis_2_uw_df) <- c("PCoA2", "sampleID")

pcoa_axis_1_2_uw_df <- merge(pcoa_axis_1_uw_df, pcoa_axis_2_uw_df, by="sampleID")

pcoa_axis_landscape_df <- data.frame(
                sampleID = proec_alpha_decontam$X.SampleID, 
                landscape = proec_alpha_decontam$landscape)

pcoa_axis_1_2_uw_landscape_df <- merge(pcoa_axis_1_2_uw_df, pcoa_axis_landscape_df, by="sampleID")

pcoa_coords_unweighted_FF <- pcoa_axis_1_2_uw_landscape_df %>% filter(landscape == "forest.fragment")
pcoa_coords_unweighted_CF <- pcoa_axis_1_2_uw_landscape_df %>% filter(landscape == "continuous.forest")
pcoa_coords_unweighted_I <- pcoa_axis_1_2_uw_landscape_df %>% filter(landscape == "island")

vector_CF_PCoA1_uw <- pcoa_coords_unweighted_CF$PCoA1
vector_FF_PCoA1_uw <- pcoa_coords_unweighted_FF$PCoA1
vector_I_PCoA1_uw <- pcoa_coords_unweighted_I$PCoA1

vector_CF_PCoA2_uw <- pcoa_coords_unweighted_CF$PCoA2
vector_FF_PCoA2_uw <- pcoa_coords_unweighted_FF$PCoA2
vector_I_PCoA2_uw <- pcoa_coords_unweighted_I$PCoA2

d.raw.data<-function(E.data,C.data){  
  d<-(mean(E.data)-mean(C.data))/
  sqrt(((length(E.data)-1)*var(E.data)+
  (length(C.data)-1)*var(C.data))/
  (length(E.data)+length(C.data)-2))
  names(d)<-"effect size d"
  return(d)
}

##############################################################################################################
#PCoA1
#A-C
effect_size <- d.raw.data(vector_CF_PCoA1_uw,vector_FF_PCoA1_uw)
psych::cohen.d.ci(d = effect_size, n1 = length(vector_FF_PCoA1_uw), n2 = length(vector_CF_PCoA1_uw))

#A-I
effect_size <- d.raw.data(vector_I_PCoA1_uw,vector_FF_PCoA1_uw)
psych::cohen.d.ci(d = effect_size, n1 = length(vector_FF_PCoA1_uw), n2 = length(vector_I_PCoA1_uw))

#C-I
effect_size <- d.raw.data(vector_I_PCoA1_uw, vector_CF_PCoA1_uw)
psych::cohen.d.ci(d = effect_size, n1 = length(vector_CF_PCoA1_uw), n2 = length(vector_I_PCoA1_uw))

##############################################################################################################
#PCoA2
#A-C
effect_size <- d.raw.data(vector_CF_PCoA2_uw,vector_FF_PCoA2_uw)
psych::cohen.d.ci(d = effect_size, n1 = length(vector_FF_PCoA2_uw), n2 = length(vector_CF_PCoA2_uw))

#A-I
effect_size <- d.raw.data(vector_I_PCoA2_uw,vector_FF_PCoA2_uw)
psych::cohen.d.ci(d = effect_size, n1 = length(vector_FF_PCoA2_uw), n2 = length(vector_I_PCoA2_uw))

#C-I
effect_size <- d.raw.data(vector_I_PCoA2_uw, vector_CF_PCoA2_uw)
psych::cohen.d.ci(d = effect_size, n1 = length(vector_CF_PCoA2_uw), n2 = length(vector_I_PCoA2_uw))

##############################################################################################################
#PERMANOVA: Null hypothesis significance testing (Supplementary Data 8 & 19)
##############################################################################################################

SequencingDepth <- scale(as(sample_data(proec_beta_decontam), "data.frame")$SequencingDepth)
Site <- as(sample_data(proec_beta_decontam), "data.frame")$capture_site
landscape <- as(sample_data(proec_beta_decontam), "data.frame")$landscape
density <- as(sample_data(proec_beta_decontam), "data.frame")$density
extractionBatch <- as(sample_data(proec_beta_decontam), "data.frame")$extractionBatch
sex <- as(sample_data(proec_beta_decontam), "data.frame")$sex
season <- as(sample_data(proec_beta_decontam), "data.frame")$season

#weighted UniFrac without density
PERMANOVA_strat_all_w <- adonis(DistW ~ landscape/Site + landscape + season + sex + SequencingDepth + extractionBatch, strata = Site, data = as(sample_data(proec_beta_decontam), "data.frame"), permutations=9999)

#weighted UniFrac with density
PERMANOVA_strat_all_w_density <- adonis(DistW ~ landscape/Site + landscape + season + sex + SequencingDepth + extractionBatch + density, strata = Site, data = as(sample_data(proec_beta_decontam), "data.frame"), permutations=9999)

#unweighted UniFrac without density
PERMANOVA_strat_all_uw <- adonis(DistUW ~ landscape/Site + landscape + season + sex + SequencingDepth + extractionBatch, strata = Site, data = as(sample_data(proec_beta_decontam), "data.frame"), permutations=9999)

#unweighted UniFrac with density
PERMANOVA_strat_all_uw_density <- adonis(DistUW ~ landscape/Site + landscape + season + sex + SequencingDepth + extractionBatch + density, strata = Site, data = as(sample_data(proec_beta_decontam), "data.frame"), permutations=9999)