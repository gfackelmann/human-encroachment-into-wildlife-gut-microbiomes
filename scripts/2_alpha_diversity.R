#Copyright (c) 2021 Insitute of Evolutionary Ecology and Conservation Genomics, Ulm University, Germany. See the LICENSE file at the top-level directory of this distribution (https://github.com/gfackelmann/human-encroachment-into-wildlife-gut-microbiomes/blob/main/LICENSE). 
#Author: Gloria Fackelmann (gloria.fackelmann@uni-ulm.de; https://www.researchgate.net/profile/Gloria-Fackelmann)

##############################################################################################################
#Load necessary libraries
##############################################################################################################

#R version 3.6.1
library("phyloseq") #phyloseq_1.28.0
library("btools") #btools_0.0.1
library("MuMIn") #MuMIn_1.43.6
library("glmmTMB") #glmmTMB_0.2.3
library("lme4") #lme4_1.1-21

##############################################################################################################
#Calculate alpha diversity metrics
##############################################################################################################
#phyloseq object is 'proec_pellet_depth_decontam'. Metadata file behind phyloseq object is 'proec_alpha_decontam'.

#Faith's PD
PD <- estimate_pd(proec_pellet_depth_decontam)
proec_alpha_decontam$NotRarePD <- PD$PD

#observed number of ASVs and Shannon Diversity
Richness <- estimate_richness(proec_pellet_depth_decontam, measures=c("Observed","Shannon"))
proec_alpha_decontam$NotRareObserved <- Richness$Observed
proec_alpha_decontam$NotRareShannon <- Richness$Shannon

#calculate sequencing depth
proec_alpha_decontam$SequencingDepth <- sample_sums(proec_pellet_depth_decontam)

##############################################################################################################
#Generalized linear mixed models: Observed number of ASVs (Supplementary Data 1, 4 & 5)
##############################################################################################################
#without density and optimizer to dredge:
m1 <- glmer.nb(NotRareObserved ~ scale(SequencingDepth) + landscape + season + sex + (1|landscape:capture_site) + (1|extractionBatch), na.action=na.fail, data=proec_alpha_decontam)

#update with optimizer
m1 <- update(m1, control=glmerControl(optimizer="bobyqa"))

#Information-theoretic approach
nullmodel <- MuMIn:::.nullFitRE(m1)
d_m1 <- dredge(m1, extra =list(R2 = function(x) {
    r.squaredGLMM(x, null = nullmodel)["delta", ]
}))

confset.95p.m1 <- get.models(d_m1, cumsum(weight) <= 0.95)

avgmod.95p.m1 <- model.avg(confset.95p.m1)

cond_coef_table_avgmod.95p.m1 <- summary(avgmod.95p.m1)$coefmat.subset

CIs_avgmod.95p.m1 <- confint(avgmod.95p.m1)

#check spatial autocorrelation:

#full model without spatial autocorrelation
m2 <- glmmTMB(NotRareObserved ~ scale(SequencingDepth) + landscape + season + sex + (1|landscape:capture_site) + (1|extractionBatch), na.action=na.fail, family = nbinom2(), data=proec_alpha_decontam)

#full model with spatial autocorrelation
m3 <- glmmTMB(NotRareObserved ~ scale(SequencingDepth) + landscape + season + sex + (1|landscape:capture_site) + (1|extractionBatch) + exp(coordinates + 0 | capture_site), na.action=na.fail, family = nbinom2(), data=proec_alpha_decontam)

AICc(m2)
AICc(m3)

logLik(m2)
logLik(m3)

AICc_weights <- function(delta1, delta2){
  return(exp(-0.5*delta1)/(exp(-0.5*delta1)+exp(-0.5*delta2)))
  }

AICc_weights(0, (AICc(m3)-AICc(m2)))
AICc_weights((AICc(m3)-AICc(m2)), 0)

##############################################################################################################
#Generalized linear mixed models: Shannon Diversity (Supplementary Data 2, 4 & 5)
##############################################################################################################
#without density and optimizer to dredge:
m4 <- glmer(NotRareShannon ~ scale(SequencingDepth) + landscape + season + sex + (1|landscape:capture_site) + (1|extractionBatch), family=Gamma(log), na.action=na.fail, data=proec_alpha_decontam)

#update with optimizer
m4 <- update(m4, control=glmerControl(optimizer="bobyqa"))

#Information-theoretic approach
nullmodel <- MuMIn:::.nullFitRE(m4)
d_m4 <- dredge(m4, extra =list(R2 = function(x) {
    r.squaredGLMM(x, null = nullmodel)["delta", ]
}))

confset.95p.m4 <- get.models(d_m4, cumsum(weight) <= 0.95)

avgmod.95p.m4 <- model.avg(confset.95p.m4)

cond_coef_table_avgmod.95p.m4 <- summary(avgmod.95p.m4)$coefmat.subset

CIs_avgmod.95p.m4 <- confint(avgmod.95p.m4)

#check spatial autocorrelation:

#full model without spatial autocorrelation
m5 <- glmmTMB(NotRareShannon ~  scale(SequencingDepth) + landscape + season + sex + (1|landscape:capture_site) + (1|extractionBatch), na.action=na.fail, family = Gamma(log), data=proec_alpha_decontam)

##full model with spatial autocorrelation
m6 <- glmmTMB(NotRareShannon ~ scale(SequencingDepth) + landscape + season + sex + (1|landscape:capture_site) + (1|extractionBatch) + exp(coordinates + 0 | capture_site), na.action=na.fail, family = Gamma(log), data=proec_alpha_decontam)

AICc(m5)
AICc(m6)

logLik(m5)
logLik(m6)

AICc_weights <- function(delta1, delta2){
  return(exp(-0.5*delta1)/(exp(-0.5*delta1)+exp(-0.5*delta2)))
  }

AICc_weights(0, (AICc(m6)-AICc(m5)))
AICc_weights((AICc(m6)-AICc(m5)), 0)

##############################################################################################################
#Generalized linear mixed models: Faith's PD (Supplementary Data 3-5)
##############################################################################################################
#without density and optimizer to dredge:
m7 <- glmer(NotRarePD ~ scale(SequencingDepth) + landscape + season + sex + (1|landscape:capture_site) + (1|extractionBatch), family=Gamma(log), na.action=na.fail, data=proec_alpha_decontam)

#update with optimizer
m7 <- update(m7, control=glmerControl(optimizer="bobyqa"))

#Information-theoretic approach
nullmodel <- MuMIn:::.nullFitRE(m7)
d_m7 <- dredge(m7, extra =list(R2 = function(x) {
    r.squaredGLMM(x, null = nullmodel)["delta", ]
}))

confset.95p.d_m7 <- get.models(d_m7, cumsum(weight) <= 0.95)

avgmod.95p.d_m7 <- model.avg(confset.95p.d_m7)

cond_coef_table_avgmod.95p.d_m7 <- summary(avgmod.95p.d_m7)$coefmat.subset

CIs_avgmod.95p.d_m7 <- confint(avgmod.95p.d_m7)

#check spatial autocorrelation:

#full model without spatial autocorrelation
m8 <- glmmTMB(NotRarePD ~ scale(SequencingDepth) + landscape + season + sex + (1|landscape:capture_site) + (1|extractionBatch), na.action=na.fail, family = Gamma(log), data=proec_alpha_decontam)

#full model with spatial autocorrelation
m9 <- glmmTMB(NotRarePD ~ scale(SequencingDepth) + landscape + season + sex + (1|landscape:capture_site) + (1|extractionBatch) + exp(coordinates + 0 | capture_site), na.action=na.fail, family = Gamma(log), data=proec_alpha_decontam)

AICc(m8)
AICc(m9)

logLik(m8)
logLik(m9)

AICc_weights <- function(delta1, delta2){
  return(exp(-0.5*delta1)/(exp(-0.5*delta1)+exp(-0.5*delta2)))
  }

AICc_weights(0, (AICc(m9)-AICc(m8)))
AICc_weights((AICc(m9)-AICc(m8)), 0)

##############################################################################################################
#Generalized linear mixed models: Observed ASVs, Shannon Diversity & Faith's PD with host density (Supplementary Data 13-16)
##############################################################################################################
#Observed ASVs
observed_density <- glmer.nb(NotRareObserved ~ scale(SequencingDepth) + landscape + season + sex + (1|landscape:capture_site) + (1|extractionBatch) + scale(density), na.action=na.fail, data=proec_alpha_decontam)

observed_density_opt <- update(observed_density, control=glmerControl(optimizer="bobyqa"))

nullmodel <- MuMIn:::.nullFitRE(observed_density_opt)
d_observed_density_opt <- dredge(observed_density_opt, extra =list(R2 = function(x) {
    r.squaredGLMM(x, null = nullmodel)["delta", ]
}))

confset.95p.d_observed_density_opt <- get.models(d_observed_density_opt, cumsum(weight) <= 0.95)
avgmod.95p.d_observed_density_opt <- model.avg(confset.95p.d_observed_density_opt)
cond_coef_table_avgmod.95p.d_observed_density_opt <- summary(avgmod.95p.d_observed_density_opt)$coefmat.subset
CIs_avgmod.95p.d_observed_density_opt <- confint(avgmod.95p.d_observed_density_opt)

#Shannnon diversity
shannon_density <- glmer(NotRareShannon ~ scale(SequencingDepth) + landscape + season + sex + (1|landscape:capture_site) + (1|extractionBatch) + scale(density), family=Gamma(log), na.action=na.fail, data=proec_alpha_decontam)

shannon_density_opt <- update(shannon_density, control=glmerControl(optimizer="bobyqa"))

nullmodel <- MuMIn:::.nullFitRE(shannon_density_opt)
d_shannon_density_opt <- dredge(shannon_density_opt, extra =list(R2 = function(x) {
    r.squaredGLMM(x, null = nullmodel)["delta", ]
}))

confset.95p._shannon_density_opt <- get.models(d_shannon_density_opt, cumsum(weight) <= 0.95)
avgmod.confset.95p._shannon_density_opt <- model.avg(confset.95p._shannon_density_opt)
cond_coef_table_avgmod.confset.95p._shannon_density_opt <- summary(avgmod.confset.95p._shannon_density_opt)$coefmat.subset
CIs_avgmod.confset.95p._shannon_density_opt <- confint(avgmod.confset.95p._shannon_density_opt)

#Faith's PD
PD_density <- glmer(NotRarePD ~ scale(SequencingDepth) + landscape + season + sex + (1|landscape:capture_site) + (1|extractionBatch) + scale(density), family=Gamma(log), na.action=na.fail, data=proec_alpha_decontam)

PD_density_opt <- update(PD_density, control=glmerControl(optimizer="bobyqa"))

nullmodel <- MuMIn:::.nullFitRE(PD_density_opt)
d_PD_density_opt <- dredge(PD_density_opt, extra =list(R2 = function(x) {
    r.squaredGLMM(x, null = nullmodel)["delta", ]
}))

confset.95p.d_PD_density_opt <- get.models(d_PD_density_opt, cumsum(weight) <= 0.95)
avgmod.95p.d_PD_density_opt <- model.avg(confset.95p.d_PD_density_opt)
cond_coef_table_avgmod.95p.d_PD_density_opt <- summary(avgmod.95p.d_PD_density_opt)$coefmat.subset
CIs_avgmod.95p.d_PD_density_opt <- confint(avgmod.95p.d_PD_density_opt)