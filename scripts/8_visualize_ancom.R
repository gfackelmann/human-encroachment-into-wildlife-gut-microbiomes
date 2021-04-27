#Copyright (c) 2021 Insitute of Evolutionary Ecology and Conservation Genomics, Ulm University, Germany. See the LICENSE file at the top-level directory of this distribution (https://github.com/gfackelmann/human-encroachment-into-wildlife-gut-microbiomes/blob/main/LICENSE). 
#Author: Gloria Fackelmann (gloria.fackelmann@uni-ulm.de; https://www.researchgate.net/profile/Gloria-Fackelmann)

##############################################################################################################
#Load necessary libraries
##############################################################################################################

#R version 3.6.1
library("stats") #stats_3.6.1

##############################################################################################################
#Calculate clr F stat
##############################################################################################################

#The volcano plots (Figure 2b, Supplementary Figures 6b-9b) have the W value for each ASV (output from ANCOM) on the y axis, and the F stat from the clr (centered log-ratio transformaton) for the x axis, which is not an output of this ANCOM version, so it has to be calculated separately. The clr is applied to the raw ASV table for each sample. From each value for each sample, the log of that value is subtracted from the mean log of all the values in that sample.

#add +1 to all values in the dataframe, since that is what ANCOM does as well
proec_pellet_depth_decontam <- proec_pellet_depth_decontam +1
part1_2_3_4_strat_PWY_2941 <- part1_2_3_4_strat_PWY_2941 +1
part1_2_3_4_strat_TEICHOICACID <- part1_2_3_4_strat_TEICHOICACID +1
part1_2_3_4_strat_PWY_7210 <- part1_2_3_4_strat_PWY_7210 +1
part1_2_3_4_strat_PWY_1861 <- part1_2_3_4_strat_PWY_1861 +1

#define the clr function
clr_coordinates <- function(X) {
  lX = log(X)
  lX - apply(lX, 1, mean)
} 

#apply the clr function
clr_values_coords_decontam <- clr_coordinates(proec_pellet_depth_decontam)
clr_values_coords_PWY_2941 <- clr_coordinates(part1_2_3_4_strat_PWY_2941)
clr_values_coords_TEICHOICACID <- clr_coordinates(part1_2_3_4_strat_TEICHOICACID)
clr_values_coords_PWY_7210 <- clr_coordinates(part1_2_3_4_strat_PWY_7210)
clr_values_coords_PWY_1861 <- clr_coordinates(part1_2_3_4_strat_PWY_1861)

#The following code was adapted from qiime2 which generates a volcano plot for the ANCOM test (ANCOM is implemented in the q2-composition, skbio.stats.composition.ancom) by Gloria Fackelmann in order to derive the F stat by running an anova on each ASV (not each sample) and extracing the F values for each ASV (this is how the matrix of clr values is turned into a 1 dimensional list of numbers that can then be plotted): 
#https://github.com/qiime2/q2-composition/blob/master/q2_composition/_ancom.py 
#https://github.com/qiime2/q2-composition/blob/master/q2_composition/plugin_setup.py 
#https://github.com/qiime2/q2-composition/blob/master/q2_composition/_impute.py

calculate_f_vals <- function(otu_table_with_landscape, name_of_new_df){
  
  n_otu <- dim(otu_table_with_landscape)[2]-1 #ASV table contains variable of interest (in this case landscape) in the last column
  otu_ids <- colnames(otu_table_with_landscape)
  #define the vector to store the value in
  fValues <- vector(mode="double", length=n_otu)
  est_FF <- vector(mode="double", length=n_otu)
  est_CI <- vector(mode="double", length=n_otu)
  se_FF <- vector(mode="double", length=n_otu)
  se_CI <- vector(mode="double", length=n_otu)
  #run the anova on each OTU
  for(ii in 1:(n_otu)){
    #define the summary anova function:
    command <- paste("fValTest <- summary(aov( formula=",otu_ids[ii], " ~ landscape, data = otu_table_with_landscape))")
    #call the summary anova function
    eval(parse(text=command))
    #extract the f values from the anova and store it in the vector
    fValues[ii] <- fValTest[[1]]$F[1]
    #define the anova function:
    command2 <- paste("parameters.se <- aov( formula=",otu_ids[ii], " ~ landscape, data = otu_table_with_landscape)")
    #call the anova function
    eval(parse(text=command2))
    #extract the parameter estimates and standard errors from the anova and store it in the vector
    est_FF[ii] <- coef(summary.lm(parameters.se))[1] #estimate from the forest fragment
    est_CI[ii] <- coef(summary.lm(parameters.se))[2] #estimates from continuous forest/island combined landscape
    se_FF[ii] <- coef(summary.lm(parameters.se))[3] #standard erros from the forest fragment
    se_CI[ii] <- coef(summary.lm(parameters.se))[4] #standard erros from continuous forest/island combined landscape
  }
  
  
  #create a dataframe with the F values and the corresponding OTU IDs
  name_of_new_df <- data.frame(
    otu.names = otu_ids[1:length(colnames(otu_table_with_landscape))-1], 
    fVals = fValues,
    FF_estimates = est_FF,
    CI_estimates = est_CI,
    FF_SE = se_FF,
    CI_SE = se_CI
    )
  
  return(name_of_new_df)
}

#All ASV tables contains variable of interest (in this case landscape) in the last column:
#apply 'calculate_f_vals' function to the whole dataset and to each pathway:
decontam_f_vals <- calculate_f_vals(clr_values_coords_decontam, decontam_f_vals)
PWY_2941_f_vals <- calculate_f_vals(clr_values_coords_PWY_2941, PWY_2941_f_vals)
TEICHOICACID_f_vals <- calculate_f_vals(clr_values_coords_TEICHOICACID, TEICHOICACID_f_vals)
PWY_7210_f_vals <- calculate_f_vals(clr_values_coords_PWY_7210, PWY_7210_f_vals)
PWY_1861_f_vals <- calculate_f_vals(clr_values_coords_PWY_1861, PWY_1861_f_vals)

#For Figure 2 and Supplementary Figures 6-9, the fVals, FF_estimates, CI_estimates, FF_SE, CI_SE, and W_stat for each ASV identified as differentially abundant by ANCOM can be plotted.