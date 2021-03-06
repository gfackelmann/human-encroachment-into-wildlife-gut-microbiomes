#Copyright (c) 2021 Insitute of Evolutionary Ecology and Conservation Genomics, Ulm University, Germany. See the LICENSE file at the top-level directory of this distribution (https://github.com/gfackelmann/human-encroachment-into-wildlife-gut-microbiomes/blob/main/LICENSE). 
#Author: Gloria Fackelmann (gloria.fackelmann@uni-ulm.de; https://www.researchgate.net/profile/Gloria-Fackelmann)

##############################################################################################################
#Load necessary libraries
##############################################################################################################

#R version 3.6.1
library("exactRankTests") #exactRankTests_0.8-31
library("foreach") #foreach_1.5.0
source("/utils/funcs_interval_of_5.R")

##############################################################################################################
#ANCOM on ASV table (Figure 2)
##############################################################################################################

system.time({ancom_taxonomy <- ANCOM.main(
	OTUdat=proec_pellet_depth_decontam, #data frame (sample x ASVs; ASVs are columns)
    Vardat=proec_alpha_decontam,
    adjusted=FALSE,
    repeated=FALSE,
    main.var="combinedLandscape",
	longitudinal = FALSE,
    repeat.var=NULL,
    multcorr=2,
    sig=0.05,
    prev.cut=1)
})

##############################################################################################################
#ANCOM on unstrat table generated by PICRUSt2 (see PICRUSt2_script.py)
##############################################################################################################

system.time({ancom_picrust_unstrat <- ANCOM.main(
	OTUdat=path_abun_unstrat_descrip, #data frame (sample x pathways; pathways are columns; make sure cell values are numeric)
    Vardat=proec_alpha_decontam,
    adjusted=FALSE,
    repeated=FALSE,
    main.var="combinedLandscape",
	longitudinal = FALSE,
    repeat.var=NULL,
    multcorr=2,
    sig=0.05,
    prev.cut=1)
})

##############################################################################################################
#ANCOM on strat table generated by PICRUSt2, subsetted to each pathway identified in unstrat table by ANCOM
##############################################################################################################

#Before we can run ANCOM, we will subset the 4 tables according to the pathways previously identified by ANCOM (unstrat table), then merge all 4 tables and then create one table per pathway identified by ANCOM on the strat table to run in ANCOM:

#read output tables from PICRUSt2 run on the 4 subsets:
part1_strat <- read.table("/part1/pathways_out/path_abun_strat.tsv", header=TRUE, sep="")
part2_strat <- read.table("/part2/pathways_out/path_abun_strat.tsv", header=TRUE, sep="")
part3_strat <- read.table("/part3/pathways_out/path_abun_strat.tsv", header=TRUE, sep="")
part4_strat <- read.table("/part4/pathways_out/path_abun_strat.tsv", header=TRUE, sep="")

#subset each of the 4 tables according to pathway:
part1_strat_sub <- subset(part1_strat, pathway=="PWY-2941" | pathway=="TEICHOICACID-PWY" | pathway=="PWY-7210" | pathway=="PWY-1861")
part2_strat_sub <- subset(part2_strat, pathway=="PWY-2941" | pathway=="TEICHOICACID-PWY" | pathway=="PWY-7210" | pathway=="PWY-1861")
part3_strat_sub <- subset(part3_strat, pathway=="PWY-2941" | pathway=="TEICHOICACID-PWY" | pathway=="PWY-7210" | pathway=="PWY-1861")
part4_strat_sub <- subset(part4_strat, pathway=="PWY-2941" | pathway=="TEICHOICACID-PWY" | pathway=="PWY-7210" | pathway=="PWY-1861")

#merge all 4 tables together:
part1_2_strat <- merge(part1_strat_sub, part2_strat_sub, by = c("pathway", "sequence"), all = TRUE)
part1_2_3_strat <- merge(part1_2_strat, part3_strat_sub, by = c("pathway", "sequence"), all = TRUE)
part1_2_3_4_strat <- merge(part1_2_3_strat, part4_strat_sub, by = c("pathway", "sequence"), all = TRUE)

#replace NAs with zeros:
part1_2_3_4_strat[is.na(part1_2_3_4_strat)] <- 0

#create one table per pathway
part1_2_3_4_strat_PWY_2941 <- subset(part1_2_3_4_strat, pathway=="PWY-2941")
part1_2_3_4_strat_TEICHOICACID <- subset(part1_2_3_4_strat, pathway=="TEICHOICACID-PWY")
part1_2_3_4_strat_PWY_7210 <- subset(part1_2_3_4_strat, pathway=="PWY-7210")
part1_2_3_4_strat_PWY_1861 <- subset(part1_2_3_4_strat, pathway=="PWY-1861")

##############################################################################################################
#Pathway 1861 (Supplementary Figure 6)
system.time({ancom_picrust_part1_2_3_4_strat_PWY_1861 <- ANCOM.main(
	OTUdat=part1_2_3_4_strat_PWY_1861, #data frame (sample x pathways; pathways are columns; make sure cell values are numeric)
    Vardat=proec_alpha_decontam,
    adjusted=FALSE,
    repeated=FALSE,
    main.var="combinedLandscape",
	longitudinal = FALSE,
    repeat.var=NULL,
    multcorr=2,
    sig=0.05,
    prev.cut=1)
})

##############################################################################################################
#Pathway 2941 (Supplementary Figure 7)
system.time({ancom_picrust_part1_2_3_4_strat_PWY_2941 <- ANCOM.main(
	OTUdat=part1_2_3_4_strat_PWY_2941, #data frame (sample x pathways; pathways are columns; make sure cell values are numeric)
    Vardat=proec_alpha_decontam,
    adjusted=FALSE,
    repeated=FALSE,
    main.var="combinedLandscape",
	longitudinal = FALSE,
    repeat.var=NULL,
    multcorr=2,
    sig=0.05,
    prev.cut=1)
})

##############################################################################################################
#Teichoic acid pathway (Supplementary Figure 8)
system.time({ancom_picrust_part1_2_3_4_strat_TEICHOICACID <- ANCOM.main(
	OTUdat=part1_2_3_4_strat_TEICHOICACID, #data frame (sample x pathways; pathways are columns; make sure cell values are numeric)
    Vardat=proec_alpha_decontam,
    adjusted=FALSE,
    repeated=FALSE,
    main.var="combinedLandscape",
	longitudinal = FALSE,
    repeat.var=NULL,
    multcorr=2,
    sig=0.05,
    prev.cut=1)
})

##############################################################################################################
#Pathway 7210 (Supplementary Figure 9)
system.time({ancom_picrust_part1_2_3_4_strat_PWY_7210 <- ANCOM.main(
	OTUdat=part1_2_3_4_strat_PWY_7210, #data frame (sample x pathways; pathways are columns; make sure cell values are numeric)
    Vardat=proec_alpha_decontam,
    adjusted=FALSE,
    repeated=FALSE,
    main.var="combinedLandscape",
	longitudinal = FALSE,
    repeat.var=NULL,
    multcorr=2,
    sig=0.05,
    prev.cut=1)
})