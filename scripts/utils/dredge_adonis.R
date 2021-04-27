#Copyright (c) 2021 Insitute of Evolutionary Ecology and Conservation Genomics, Ulm University, Germany. See the LICENSE file at the top-level directory of this distribution (https://github.com/gfackelmann/human-encroachment-into-wildlife-gut-microbiomes/blob/main/LICENSE). 
#Author: Gloria Fackelmann (gloria.fackelmann@uni-ulm.de; https://www.researchgate.net/profile/Gloria-Fackelmann)

dredge_adonis <- function(distance_matrix, predictors, exclusion, strata, phyloseq_object){

  print("Copyright (c) 2021 Insitute of Evolutionary Ecology and Conservation Genomics, Ulm University, Germany. See the LICENSE file at the top-level directory of this distribution (https://github.com/gfackelmann/human-encroachment-into-wildlife-gut-microbiomes/blob/main/LICENSE). Author: Gloria Fackelmann (gloria.fackelmann@uni-ulm.de)", quote = FALSE)

  ###############################################################################################
  #Step 0: Load the necessary libaries
  ###############################################################################################

  library("phyloseq") #phyloseq_1.28.0
  library("vegan") #vegan_2.5-5
  library("dplyr") #dplyr_0.8.3
  source("/Source_codes/AICc_PERMANOVA_DREDGE.R")

  ###############################################################################################
  #Step 1: Create a list with all the possible combinations that can be made using the predictors (explanatory variables) given
  ###############################################################################################
  
  all_possible_combinations <- list()
  for (i in 1:(2^(length(predictors)))){
  	current_combination <- c()
  	for (x in 1:(length(predictors))){ 
  		checker <- bitwAnd(i-1, (bitwShiftL(1, x-1)))     
  		if (checker > 0) { 
  			current_combination <- c(current_combination, predictors[x])
  		}
  	}
  	is_valid <- TRUE
    for (pair in exclusion){
      contains_element1 <- is.element(pair[1], current_combination)
      contains_element2 <- is.element(pair[2], current_combination)
      is_valid <- !(contains_element1 && contains_element2) && is_valid 
    }
    if (is_valid){ 
      all_possible_combinations[[length(all_possible_combinations)+1]] <- current_combination	
    } 
  }
  
  ###############################################################################################
  #Step 2: Paste all of the combinations into the adonis function and call this function to calculate one model for each combination. This step takes a while to compute.
  ###############################################################################################

  results <- list()

  for (elements in 1:length(all_possible_combinations)){
    command <- paste("results[[elements]] <- adonis(",distance_matrix, "~ ", paste(all_possible_combinations[[elements]], collapse = " + "), " ,", strata, ", as(sample_data(", phyloseq_object, "), \"data.frame\"))")
    eval(parse(text=command))
  }
  
  ###############################################################################################
  #Step 3: Extract the AICc, residuals and log likelihood for all the generated models, including the null model, which only needs to be called once on any of the combined models.
  ###############################################################################################
  
  model_residuals <- c()
  aicc_list <- c()
  loglik_list <- c()
  r_squared <- c()
  model_names <- c()

  for (model in results){
    model_residuals <- c(model_residuals, model[[1]][[5]][[length(model[[1]][[5]])-1]]) 
    aicc_list <- c(aicc_list, AICc.PERMANOVA(model)$AICc)   #accesses the AICc.PERMANOVA() function in the AICc_PERMANOVA_DREDGE.R file in the Source_codes files of this respository
    loglik_list <- c(loglik_list, AICc.PERMANOVA(model)$loglik)
  }
  
  model_residuals <- c(model_residuals, 1.0)
  aicc_list <- c(aicc_list, AICc.NULL(results[[1]])$AICc)
  loglik_list <- c(loglik_list, AICc.NULL(results[[1]])$loglik)

  for (x in all_possible_combinations){
    model_names <- c(model_names, paste(x, collapse = " + "))
  }

  for (z in model_residuals){
    r_squared <- c(r_squared, (1-z))
  }

  model_names <- c(model_names, "Null Model")
  meta_list <- list("Model" = model_names, "ResidualR^2" = model_residuals, "R^2" = r_squared, "LogLikelihood" = loglik_list, "AICc" = aicc_list)
  df_all_models <- as.data.frame(meta_list, check.names=FALSE)
  
  ###############################################################################################
  #Step 4:Create a list for each predictor that adds a plus sign if the predictor was in a combination and leaves a blank if it was not in a combination.
  ###############################################################################################
  
  all_possible_combinations <- c(all_possible_combinations, "")

  for (k in predictors){
    all_marks <- c()
    for (model in all_possible_combinations){
      current_mark <- "" 
      for (variable in model){
        if (variable == k){
          current_mark <- "+"
        }
      }
      all_marks <- c(all_marks, current_mark)
    }
    df_all_models[[k]] <- all_marks 
  }
  
  ###############################################################################################
  #Step 5: Calculate delta AICc and the AICc weight (using the relative likelihood)
  ###############################################################################################
  
  df_all_models <- df_all_models[order(df_all_models$AICc),]
  delta <- c()
  
  for (number in df_all_models$AICc){
    delta <- c(delta, number-df_all_models$AICc[1])
  }
  
  df_all_models$delta <- delta

  relative_likelihood <- c()
  for (y in delta){
    relative_likelihood <- c(relative_likelihood, (exp(y * -0.5)))
  }
  sum_relative_likelihood <- sum(relative_likelihood)

  AICc_weight <- c()
  for (k in relative_likelihood){
    AICc_weight <- c(AICc_weight, (k/sum_relative_likelihood))
  }
  
  df_all_models$AICc_weight <- AICc_weight

  ###############################################################################################
  #Step 6: Assemble the dataframe/table which will be the final output of this function
  ###############################################################################################

  df_all_models <- dplyr::select(df_all_models, AICc_weight, everything())
  df_all_models <- dplyr::select(df_all_models, delta, everything()) 
  df_all_models <- dplyr::select(df_all_models, AICc, everything())

  model_names_arranged <- c()
  model_names_arranged <- df_all_models$Model
  df_all_models$Model <- NULL

  name.width <- max(sapply(names(df_all_models), nchar))
  names(df_all_models) <- format(names(df_all_models), width = name.width, justify = "centre")
  df_all_models <- format(df_all_models, width = name.width, justify = "centre")

  df_all_models$Model <- model_names_arranged

  return(df_all_models)
}

