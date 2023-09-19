#####
# 
## Step 3c: SDM corrections
#
#####

### load libraries

library(tidyverse)
library(randomForest)
library(poibin)

### wd is one level up, recommend set up Rproj or:
# setwd("../")

### selected sessionInfo()
# R version 4.1.3 (2022-03-10)
# Platform: x86_64-w64-mingw32/x64 (64-bit)
# Running under: Windows 10 x64 (build 17763)
# 
# attached base packages:
#   [1] stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#   [1] poibin_1.5         randomForest_4.7-1 forcats_0.5.1      stringr_1.4.1      dplyr_1.0.8        purrr_0.3.4       
# [7] readr_2.1.2        tidyr_1.2.0        tibble_3.1.6       ggplot2_3.3.5      tidyverse_1.3.1   

####
# 1. load data
####

# species subset: at least 5 observations. used for final model fits
# add flag for auc>0.7 cutoff
rich.flag <- "subset"
# rich.flag <- "auc07"
pres.in <- read.csv(paste0("processed_data/understory_model/understory_",rich.flag,"_presence_cover.csv")) %>%
  dplyr::select(-cover)

# richness of subset
rich.in <- pres.in %>%
  group_by(plot_id) %>%
  summarise(richness=sum(pres))

# predictors: subset with 3 climate and 3 forest
pred.in <-  read.csv("processed_data/understory_model/predictors_subset.csv") %>%
  # only final predictor set, 3 climate, 3 forest, 3 site
  dplyr::select(c(plot_id,mean_temp,summer_prec,rad,TSFdec,BA,prop_fasy,sand,whc,soil_fert))

####
# 2. prep data
####

# predictor subset
pred.std <- pred.in %>%
  mutate(across(c(mean_temp:soil_fert), ~scale(.))) %>%
  group_by(plot_id) %>%
  mutate(across(everything(), ~as.numeric(.))) %>%
  as.data.frame() %>%
  # set NAs prop_fasy to 0 (mean) as equalizing value
  mutate(across(c(prop_fasy), ~replace_na(.,0)))

# create lookup for species names, altered when transformed to columns
sp.lookup <- data.frame(species_name = unique(pres.in$species_name)) %>%
  mutate(species_alt = gsub("[_.-]","",species_name))

####
# 3. helper functions
####

### adjust probabilities using calabrese et al. 2014 approach, as implemented in Zurell et al. 2020
# code from Zurell et al. 2020 J Biogeography 47:101-113
# Calabrese correction: based on Calabrese et al. (2014) Methods in Ecology and Evolution 23: 99-112.
nLL.Calabrese <- function(par,sr,probs) {
  require(poibin)
  logit <- function(x) {
    x=ifelse(x<0.0001,0.0001,ifelse(x>0.9999,.9999,x))
    log(x/(1 - x))
  }
  invlogit <- function(x) {exp(x)/(1+exp(x))}
  
  
  bysite <- function(j) {
    logit.probs <- logit(as.numeric(probs[j,]))
    corr.probs <- invlogit( logit.probs + par[1]*sr[j] + par[2] )
    dp <- dpoibin(sr[j],as.numeric(corr.probs))
    log(ifelse(dp<.0001,.0001,dp))
  }
  - sum(sapply(seq_len(length(sr)),bysite)) 	# optim will perform minimization but we aim for maximum likelihood and thus invert
}

# also set these functions for later use
logit <- function(x) {
  x=ifelse(x<0.0001,0.0001,ifelse(x>0.9999,.9999,x))
  log(x/(1 - x))
}
invlogit <- function(x) {exp(x)/(1+exp(x))}

# probability ranking rule - adapted from the ecospat package by Zurell et al. 2020
SESAM.prr <- function (proba, sr) {
  projSR <- round(round(as.vector(sr[[1]])))
  new.prob.prr <- proba
  dataSSDM_p <- proba
  for (i in 1:nrow(proba)) {
    print(paste("test.prr, processing row ", i, sep = ""))
    SR <- projSR[i]
    if (SR > 0) {
      predcom <- dataSSDM_p[i, ]
      predcom_p <- as.matrix(dataSSDM_p[i, ])  # convert to matrix for ordering
      com <- order(predcom_p, decreasing = TRUE)
      pres <- com[1:SR]
      predcom[, pres] <- 1
      predcom[, -pres] <- 0
    }
    else {
      predcom[, ] <- 0
    }
    new.prob.prr[i, ] <- predcom
  }
  new.prob.prr
}

####
# 4. corrections with aggregated individual SDM training v. testing datasets, 20 reps of training v. test corrections/MEMs
####

### input datasets
# weighted individual SDM, 9 predictor (3 climate, 3 forest, 3 site) set
run_nbr <- "run13"
mod.name <- "mem"
mod.pred <- "9"

# training and test evaluation runs for individual SDMs
train.out <- read.csv(paste0("processed_data/understory_model/evaluation/sdm_pres_evaluation_",run_nbr,".csv"))

# MEM holdout test predictions
rich.test <- read.csv(paste0("processed_data/understory_model/evaluation/",mod.name,"_richness_",rich.flag,"_evaluation_",mod.pred,"pred.csv"))

head(train.out)

# subset in event fit sdms for more species that using in richness, model eval
train.sub <- train.out %>%
  filter(species_name %in% c(pres.in$species_name))

# do we have 150 test obs for each species?
train.sub %>%
  filter(train=="test") %>%
  group_by(species_name) %>%
  summarise(nplots=length(unique(plot_id))) %>%
  filter(nplots<150) # 6 missing observations, ok for evaluation

### create separate training and test datasets
# get averages for training and tests by species, plot totals
train.avg <- train.sub %>%
  group_by(plot_id,train,pres,species_name) %>%
  summarise(pred_mean = mean(rf_pred)) %>%
  # use altered species names
  left_join(sp.lookup, by="species_name")

# prep training, 1 row per plot
sp.train <- train.avg %>%
  filter(train=="train") %>%
  ungroup() %>%
  dplyr::select(-c(pres,train,species_name)) %>%
  pivot_wider(names_from="species_alt", values_from="pred_mean") %>%
  arrange(plot_id) 

# prep test, 1 row per plot
sp.test <- train.avg %>%
  filter(train=="test") %>%
  ungroup() %>%
  dplyr::select(-c(pres,train,species_name)) %>%
  pivot_wider(names_from="species_alt", values_from="pred_mean") %>%
  arrange(plot_id) 

### compare multiple approaches: individual SDMs, with Calabrese adjustment, with MEM constraints, with probability ranking rule applied
# go through each model fit run
sample(1:1000, 20)
seeds <- c(524, 307, 249, 910, 480, 613, 205, 628, 927, 981, 973, 945, 707, 252, 904, 24, 291, 720, 461, 869)

rich.out <- data.frame()
pred.corrprr <- data.frame()

# iterate through MEM training/test runs
for(i in 1:20) {
  print(i)
  
  # extract test plots
  rich.testrun <- rich.test %>%
    filter(run==i) %>%
    arrange(plot_id)
  
  # extract training plots
  rich.trainrun <- rich.in %>%
    filter(!plot_id %in% c(rich.testrun$plot_id)) %>%
    arrange(plot_id)
  
  # subset same plots in training and test datasets for species predictions
  sp.trainrun <- sp.train %>%
    filter(!plot_id %in% c(rich.testrun$plot_id)) %>%
    # ensure in same order
    arrange(plot_id) %>%
    dplyr::select(-plot_id)
  
  sp.testrun <- sp.test %>%
    filter(plot_id %in% c(rich.testrun$plot_id)) %>%
    arrange(plot_id) %>%
    dplyr::select(-plot_id) %>%
    # assign 0s to any NAs to prevent errors
    mutate(across(everything(),~replace_na(.,0)))
  
  # stack individual sdm probabilities and mem richness predictions
  prob.stack <- rowSums(sp.testrun) # sum test probabilities
  mem.stack <- as.array(rich.testrun$richness_pred) # sum mem predictions
  
  # estimate correction parameters based on Calabrese correction, use real richness from training data
  set.seed(seeds[i])
  adj.par.nll <- optim(par=c(0,0), fn=nLL.Calabrese, sr= rich.trainrun$richness, probs= sp.trainrun) 
  
  # correct test data probabilities using probability stack predictions
  prob.corr.probsum.df <- data.frame( apply(sp.testrun,2,FUN=function(x){invlogit(logit(x)+adj.par.nll$par[1]*prob.stack+adj.par.nll$par[2])}))
  rich.testrun$richness_pred_corr <- rowSums(prob.corr.probsum.df)
  
  # correct test data probabilities using MEM predictions
  prob.corr.mem.df <- data.frame( apply(sp.testrun,2,FUN=function(x){invlogit(logit(x)+adj.par.nll$par[1]*mem.stack+adj.par.nll$par[2])}))
  rich.testrun$richness_pred_memcorr <- rowSums(prob.corr.mem.df)
  
  # probability ranking rule, using probability sums as constraint
  prr.probsum.df <- SESAM.prr(sp.testrun, data.frame(prob.stack) ) # warning related to use of order() for data.frame, but results are ok
  rich.testrun$richness_pred_prr <- rowSums(prr.probsum.df)
  
  # probability ranking rule, using MEMs as constraint
  prr.mem.df <- SESAM.prr(sp.testrun, data.frame(mem.stack) )
  rich.testrun$richness_pred_memprr <- rowSums(prr.mem.df)
  
  # probability ranking rule, using sums of corrected probabilities as constraint
  prr.corr.probsum.df <- SESAM.prr(prob.corr.probsum.df, data.frame(rowSums(prob.corr.probsum.df)) )
  rich.testrun$richness_pred_corrprr <- rowSums(prr.corr.probsum.df)
  
  # probability ranking rule, using sums of corrected probabilities as constraint (but probabilities corrected based on MEM prediction) 
  prr.corr.mem.df <- SESAM.prr(prob.corr.mem.df, data.frame(rowSums(prob.corr.mem.df)) )
  rich.testrun$richness_pred_memcorrprr <- rowSums(prr.corr.mem.df)
  
  rich.out <- rbind(rich.out,rich.testrun)
  
  # probability ranking rule, using sums of corrected probabilities, individual species outputs
  sp.prrout <- rich.testrun %>%
    dplyr::select(c(plot_id,run)) %>%
    cbind(prr.corr.probsum.df) %>%
    pivot_longer(-c(plot_id,run), names_to="species_alt", values_to="prr_pres") %>%
    left_join(sp.lookup, by="species_alt") %>%
    dplyr::select(-species_alt) %>%
    mutate(train="test")
  
  pred.corrprr <- rbind(pred.corrprr,sp.prrout)
  
}

write.csv(rich.out, paste0("processed_data/understory_model/evaluation/corrected_richness_",rich.flag,"_evaluation_",run_nbr,".csv"), row.names=FALSE)
write.csv(pred.corrprr, paste0("processed_data/understory_model/evaluation/sdm_corrprr_",rich.flag,"_evaluation_",run_nbr,".csv"), row.names=FALSE)

####
# 5. final corrections for prediction
####

### input final fit predictions
train.out <- read.csv(paste0("processed_data/understory_model/evaluation/sdm_pres_finalFit_9pred.csv"))

head(train.out)

# subset in event fit sdms for more species that using in richness, model eval
train.sub <- train.out %>%
  filter(species_name %in% c(pres.in$species_name)) %>%
  left_join(sp.lookup, by="species_name")

# prepare richness and species
rich.trainrun <- rich.in %>%
  arrange(plot_id) 

sp.trainrun <- train.sub %>%
  ungroup() %>%
  dplyr::select(-c(pres,species_name)) %>%
  pivot_wider(names_from="species_alt", values_from="rf_pred") %>%
  arrange(plot_id) %>%
  dplyr::select(-plot_id)

set.seed(11)
adj.par.nll <- optim(par=c(0,0), fn=nLL.Calabrese, sr= rich.trainrun$richness, probs= sp.trainrun) 

# save this
save(adj.par.nll, file=paste0("processed_data/understory_model/final_fits/richness_prediction_corrections/richness_",rich.flag,"_corrections.RData"))

### double check implementation, using training data
# correct test data probabilities using probability stack predictions
prob.stack <- rowSums(sp.trainrun)
prob.corr.probsum.df <- data.frame( apply(sp.trainrun,2,FUN=function(x){invlogit(logit(x)+adj.par.nll$par[1]*prob.stack+adj.par.nll$par[2])}))
rich.trainrun$richness_pred_corr <- rowSums(prob.corr.probsum.df)

rich.trainrun %>%
  ggplot(aes(x=richness,y=richness_pred_corr)) +
  geom_point() +
  geom_smooth(method="lm") +
  theme_bw()

summary(lm(richness_pred_corr~richness,data=rich.trainrun))
