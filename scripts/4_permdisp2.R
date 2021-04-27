#Copyright (c) 2021 Insitute of Evolutionary Ecology and Conservation Genomics, Ulm University, Germany. See the LICENSE file at the top-level directory of this distribution (https://github.com/gfackelmann/human-encroachment-into-wildlife-gut-microbiomes/blob/main/LICENSE). 
#Author: Gloria Fackelmann (gloria.fackelmann@uni-ulm.de; https://www.researchgate.net/profile/Gloria-Fackelmann)

##############################################################################################################
#Load necessary libraries
##############################################################################################################

#R version 3.6.1
library("phyloseq") #phyloseq_1.28.0
library("vegan") #vegan_2.5-5
library("MuMIn") #MuMIn_1.43.6
library("lme4") #lme4_1.1-21
library("glmmTMB") #glmmTMB_0.2.3


##############################################################################################################
#Calculate beta diversity distances
##############################################################################################################
#phyloseq object is 'proec_beta_15'. Metadata file behind phyloseq object is 'proec_beta_15_meta'.

#calculate distances with capture sites where n > 15 individuals
#weighted unifrac
DistW_proec_beta_15 <- phyloseq::distance(proec_beta_15, method="wUniFrac")
#unweighted unifrac
DistUW_proec_beta_15 <- phyloseq::distance(proec_beta_15, method="uunifrac")

#betadisper using spatial median of sites:
#weighted UniFrac
modsW_15 <- betadisper(DistW_proec_beta_15, sample_data(proec_beta_15)$capture_site, type = "median")
proec_beta_15_meta$distance_to_centroid_w <- modsW_15$distances
#unweighted unifrac
modsUW_15 <- betadisper(DistUW_proec_beta_15, sample_data(proec_beta_15)$capture_site, type = "median")
proec_beta_15_meta$distance_to_centroid_uw <- modsUW_15$distances

##############################################################################################################
#Generalized linear mixed models: weighted UniFrac (Supplementary Data 9 & 11-12)
##############################################################################################################
#weighted unifrac without optimizer to dredge:
disp_w <- glmer(distance_to_centroid_w ~ scale(SequencingDepth) + landscape + season + sex + (1|landscape:capture_site) + (1|extractionBatch), na.action=na.fail, data=proec_beta_15_meta, family=Gamma (log))

disp_w <- update(disp_w, control=glmerControl(optimizer="bobyqa"))

#Information-theoretic approach
nullmodel <- MuMIn:::.nullFitRE(disp_w)
d_disp_w <- dredge(disp_w, extra =list(R2 = function(x) {
    r.squaredGLMM(x, null = nullmodel)["delta", ]
}))

confset.95p.d_disp_w <- get.models(d_disp_w, cumsum(weight) <= 0.95)

avgmod.95p.d_disp_w <- model.avg(confset.95p.d_disp_w)

cond_coef_table_avgmod.95p.d_disp_w <- summary(avgmod.95p.d_disp_w)$coefmat.subset

CIs_avgmod.95p.d_disp_w <- confint(avgmod.95p.d_disp_w)

#check spatial autocorrelation:

#full model without spatial autocorrelation
m1 <- glmmTMB(distance_to_centroid_w ~ scale(SequencingDepth) + landscape + season + sex + (1|landscape:capture_site) + (1|extractionBatch), na.action=na.fail, data=proec_beta_15_meta, family=Gamma(log))

#full model with spatial autocorrelation
m2 <- glmmTMB(distance_to_centroid_w ~ scale(SequencingDepth) + landscape + season + sex + (1|landscape:capture_site) + (1|extractionBatch) + exp(coordinates + 0 | capture_site), na.action=na.fail, data=proec_beta_15_meta, family=Gamma(log))

AICc(m1)
AICc(m2)

logLik(m1)
logLik(m2)

AICc_weights <- function(delta1, delta2){
  return(exp(-0.5*delta1)/(exp(-0.5*delta1)+exp(-0.5*delta2)))
  }

AICc_weights(0, (AICc(m2)-AICc(m1)))  
AICc_weights((AICc(m2)-AICc(m1)), 0)

##############################################################################################################
#Generalized linear mixed models: unweighted UniFrac (Supplementary Data 10-12)
##############################################################################################################
#unweighted unifrac without optimizer to dredge:
disp_uw <- glmer(distance_to_centroid_uw ~ scale(SequencingDepth) + landscape + season + sex + (1|landscape:capture_site) + (1|extractionBatch), na.action=na.fail, data=proec_beta_15_meta, family=Gamma (log))

disp_uw <- update(disp_uw, control=glmerControl(optimizer="bobyqa"))

#Information-theoretic approach
nullmodel <- MuMIn:::.nullFitRE(disp_uw)
d_disp_uw <- dredge(disp_uw, extra =list(R2 = function(x) {
    r.squaredGLMM(x, null = nullmodel)["delta", ]
}))

confset.95p.d_disp_uw <- get.models(d_disp_uw, cumsum(weight) <= 0.95)

avgmod.95p.d_disp_uw <- model.avg(confset.95p.d_disp_uw)

cond_coef_table_avgmod.95p.d_disp_uw <- summary(avgmod.95p.d_disp_uw)$coefmat.subset

CIs_avgmod.95p.d_disp_uw <- confint(avgmod.95p.d_disp_uw)

#check spatial autocorrelation:

#full model without spatial autocorrelation
m3 <- glmmTMB(distance_to_centroid_uw ~ scale(SequencingDepth) + landscape + season + sex + (1|landscape:capture_site) + (1|extractionBatch), na.action=na.fail, data=proec_beta_15_meta, family=Gamma(log))

#full model with spatial autocorrelation
m4 <- glmmTMB(distance_to_centroid_uw ~ scale(SequencingDepth) + landscape + season + sex + (1|landscape:capture_site) + (1|extractionBatch) + exp(coordinates + 0 | capture_site), na.action=na.fail, data=proec_beta_15_meta, family=Gamma(log))

AICc(m3)
AICc(m4)

logLik(m3)
logLik(m4)

AICc_weights <- function(delta1, delta2){
  return(exp(-0.5*delta1)/(exp(-0.5*delta1)+exp(-0.5*delta2)))
  }

AICc_weights(0, (AICc(m4)-AICc(m3)))  
AICc_weights((AICc(m4)-AICc(m3)), 0)

##############################################################################################################
#Generalized linear mixed models: weighted & unweighted UniFrac with host density (Supplementary Data 20-22)
##############################################################################################################

disp_w_density <- glmer(distance_to_centroid_w ~ scale(SequencingDepth) + landscape + season + sex + (1|landscape:capture_site) + (1|extractionBatch) + density, na.action=na.fail, data=proec_beta_15_meta, family=Gamma (log))

disp_uw_density <- glmer(distance_to_centroid_uw ~ scale(SequencingDepth) + landscape + season + sex + (1|landscape:capture_site) + (1|extractionBatch) + density, na.action=na.fail, data=proec_beta_15_meta, family=Gamma (log))

disp_w_density <- update(disp_w_density, control=glmerControl(optimizer="bobyqa"))
disp_uw_density <- update(disp_uw_density, control=glmerControl(optimizer="bobyqa"))

nullmodel <- MuMIn:::.nullFitRE(disp_w_density)
d_disp_w_density <- dredge(disp_w_density, extra =list(R2 = function(x) {
    r.squaredGLMM(x, null = nullmodel)["delta", ]
}))

nullmodel <- MuMIn:::.nullFitRE(disp_uw_density)
d_disp_uw_density <- dredge(disp_uw_density, extra =list(R2 = function(x) {
    r.squaredGLMM(x, null = nullmodel)["delta", ]
}))

confset.95p.d_disp_w_density <- get.models(d_disp_w_density, cumsum(weight) <= 0.95)
confset.95p.d_disp_uw_density <- get.models(d_disp_uw_density, cumsum(weight) <= 0.95)

avgmod.95p.d_disp_w_density <- model.avg(confset.95p.d_disp_w_density)
avgmod.95p.d_disp_uw_density <- model.avg(confset.95p.d_disp_uw_density)

cond_coef_table_avgmod.95p.d_disp_w_density <- summary(avgmod.95p.d_disp_w_density)$coefmat.subset
cond_coef_table_avgmod.95p.d_disp_uw_density <- summary(avgmod.95p.d_disp_uw_density)$coefmat.subset

CIs_avgmod.95p.d_disp_w_density <- confint(avgmod.95p.d_disp_w_density)
CIs_avgmod.95p.d_disp_uw_density <- confint(avgmod.95p.d_disp_uw_density)