#####
# 
## Step 4: Model evaluation
#
#####

### load libraries

library(tidyverse)
library(ggpubr)
library(openxlsx)
library(caret) # conf mat
library(PresenceAbsence)

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
#   [1] PresenceAbsence_1.1.10 caret_6.0-91           lattice_0.20-45        openxlsx_4.2.5         ggpubr_0.4.0          
# [6] forcats_0.5.1          stringr_1.4.1          dplyr_1.0.8            purrr_0.3.4            readr_2.1.2           
# [11] tidyr_1.2.0            tibble_3.1.6           ggplot2_3.3.5          tidyverse_1.3.1   

###
# 1. load data
###

### responses
# presence full and subset
pres.in <- read.csv("processed_data/understory_model/understory_presence_cover.csv")

# add flag for auc>0.7 cutoff
rich.flag <- "subset"
# rich.flag <- "auc07"
pres.sub <- read.csv(paste0("processed_data/understory_model/understory_",rich.flag,"_presence_cover.csv")) %>%
  dplyr::select(-cover)

# cover and richness all species
cover.in <- pres.in %>%
  group_by(plot_id) %>%
  summarise(cover=sum(cover),
            richness=sum(pres))

# richness species subset
rich.in <- pres.sub %>%
  group_by(plot_id) %>%
  summarise(richness=sum(pres))

### lookup tables
# pft lookup
pft.in <- read.csv("processed_data/species_lookup/final_pft_categories_ms.csv") %>%
  # select columns
  dplyr::select(c(species_name,species_name_ellenberg,pft))

# ellenberg lookup
ellen.in <- read.csv("processed_data/species_lookup/species_lookup_ellenberg.csv")

####
# 2. helper functions
####

# compute ellenberg mean and wtmean licht, T, F, N for this subset of species
# functionalize
ellen_compute <- function(df.in, sp.list) {
  ellen.sub <- df.in %>%
    # subset to species list
    filter(species_name %in% c(sp.list)) %>%
    # combine with ellenberg indicator data
    # first get sp name lookup from pft table
    left_join(pft.in, by="species_name") %>%
    left_join(ellen.in, by=c("species_name_ellenberg"="species_name")) %>%
    # remove missing values, also values coded as ? (unknown) or x (insensitive)
    filter(!is.na(OriglName), !OrigValueStr %in% c("?","x")) %>%
    # assign 0B as 0
    mutate(OrigValueStr = ifelse(OrigValueStr=="0B",0,OrigValueStr)) %>%
    # all values should now be numeric
    mutate(OrigValueStr=as.numeric(OrigValueStr)) %>%
    # get 1 entry per individual
    dplyr::select(c(plot_id,species_name,cover,OriglName,OrigValueStr)) %>%
    # only present species
    filter(cover>0) %>%
    # mean and wted mean
    group_by(plot_id,OriglName) %>%
    summarise(mean=mean(OrigValueStr),
              wtmean=sum(OrigValueStr*cover)/sum(cover))

}

###
# 3. MEM evaluations
###

# test fit
mod.name <- "mem"
mod.pred <- "9"

### cover and richness test fits
cover.test <- read.csv(paste0("processed_data/understory_model/evaluation/",mod.name,"_cover_evaluation_",mod.pred,"pred.csv"))
rich.test <- read.csv(paste0("processed_data/understory_model/evaluation/",mod.name,"_richness_",rich.flag,"_evaluation_",mod.pred,"pred.csv"))

# quick plot assessment: mean fit
cover.test %>%
  group_by(plot_id) %>%
  summarise(across(c(cover,cover_pred),~mean(.))) %>%
  ggplot(aes(x=cover,y=cover_pred)) +
  geom_point() +
  geom_smooth(method="lm") +
  ylim(c(min(cover.test$cover),max(cover.test$cover))) +
  xlim(c(min(cover.test$cover),max(cover.test$cover))) +
  geom_abline(intercept=0,slope=1,col="red") +
  stat_cor(
    aes(label = paste(..rr.label..)),
    label.x.npc = "left"
  ) +
  stat_cor(
    method="spearman",cor.coef.name="rho",label.x.npc="center",p.digits=NA, label.sep=""
  ) +
  theme_bw()

rich.test %>%
  group_by(plot_id) %>%
  summarise(across(c(richness,richness_pred),~mean(.))) %>%
  ggplot(aes(x=richness,y=richness_pred)) +
  geom_point() +
  geom_smooth(method="lm") +
  ylim(c(min(rich.test$richness,rich.test$richness_pred),max(rich.test$richness))) +
  xlim(c(min(rich.test$richness,rich.test$richness_pred),max(rich.test$richness))) +
  geom_abline(intercept=0,slope=1,col="red") +
  stat_cor(
    aes(label = paste(..rr.label..)),
    label.x.npc = "left"
  ) +
  stat_cor(
    method="spearman",cor.coef.name="rho",label.x.npc="center",p.digits=NA, label.sep=""
  ) +
  theme_bw()

# for cover, use same predictors as individual SDMs for consistency
# for richness, model fit is not an improvement over the sum of individual SDMs; apply final Calabrese corrections to sum of predicted probabilities rather than fitting additional MEM

### also load final model vimp, rcurves
# final model for cover
cover.vimp <- read.csv(paste0("processed_data/understory_model/evaluation/mem_cover_finalFit_vimp_9pred.csv"))

# use IncMSE
cover.vimp %>%
  mutate(variables = factor(variables)) %>%
  ggplot(aes(x=reorder(variables,-desc(X.IncMSE)),y=X.IncMSE)) +
  geom_bar(stat='identity') +
  coord_flip() +
  theme_bw() +theme(legend.position = "none")

# final model for cover
cover.rcurves <- read.csv(paste0("processed_data/understory_model/evaluation/mem_cover_finalFit_rcurves_9pred.csv"))

# plot response curves: cover
# free y to look at effect of each predictor
cover.rcurves %>%
  ggplot(aes(x=x,y=y)) +
  facet_wrap(~variable, scales="free") +
  geom_point() +
  geom_smooth(method="loess") +
  theme_bw()

# fixed y highlights effect of most important predictors
cover.rcurves %>%
  ggplot(aes(x=x,y=y)) +
  facet_wrap(~variable, scales="free_x") +
  geom_point() +
  geom_smooth(method="loess") +
  theme_bw()

###
# 4. corrected richness evaluations
###

# run
run_nbr <- "run13"

rich.out <- read.csv(paste0("processed_data/understory_model/evaluation/corrected_richness_",rich.flag,"_evaluation_",run_nbr,".csv"))

# avg predicted richness
rich.out %>%
  group_by(plot_id,richness) %>%
  summarise(across(c(richness_pred_corr:richness_pred_memcorrprr),~mean(.))) %>%
  pivot_longer(c(richness_pred_corr:richness_pred_memcorrprr)) %>%
  ggplot(aes(x=richness,y=value)) +
  facet_wrap(~name) +
  geom_point() +
  geom_smooth(method="lm") +
  ylim(c(min(rich.out$richness),max(rich.out$richness,rich.out$richness_pred_prr))) +
  xlim(c(min(rich.out$richness),max(rich.out$richness,rich.out$richness_pred_prr))) +
  geom_abline(intercept=0,slope=1,col="red") +
  stat_cor(
    aes(label = paste(..rr.label..)),
    label.x.npc = "left"
  ) +
  stat_cor(
    method="spearman",cor.coef.name="rho",label.x.npc="center",p.digits=NA, label.sep=""
  ) +
  theme_bw()

# ggsave(paste0("figures/model_evaluation/corrected_richness_",rich.flag,"_comparison_",run_nbr,".pdf"))
# model performance similar across all
# corrected probability sums best match 1:1 line

lm(richness_pred_corrprr~richness, data=rich.out) # steeper slope than mem corrected
lm(richness_pred_memcorr~richness, data=rich.out)
lm(richness_pred_prr~richness, data=rich.out)

# choose simple corrected probability or corrected probability w/ PRR, depending on individual species results

###
# 5. individual SDM and derived values (richness, Ellenberg indicator) fits
###

# run
run_nbr <- "run13"
prr <- ""

train.out <- read.csv(paste0("processed_data/understory_model/evaluation/sdm_pres_evaluation_",run_nbr,".csv"))

head(train.out)
summary(train.out)

# subset in event fit sdms for more species that using in richness, model eval
test.sub <- train.out %>%
  filter(species_name %in% c(pres.sub$species_name)) %>%
  # test data only for evaluation
  filter(train=="test")

# do we have 150 test obs for each species?
test.sub %>%
  # filter(train=="test") %>%
  group_by(species_name) %>%
  summarise(nplots=length(unique(plot_id))) %>%
  filter(nplots<150) # 6 observations, seems ok

# get averages for training and tests by species, plot totals
test.avg <- test.sub %>%
  group_by(plot_id,pres,species_name) %>%
  summarise(pred_mean = mean(rf_pred)) 

### check systematic bias in pres/abs values based on species prevalence
# corrected in training data, visible in test data
# always the case, but tempered somewhat by weighted SDM fits
# avg species prevalence
pres.prev <- pres.in %>%
  ungroup() %>%
  group_by(species_name) %>%
  summarise(prev = sum(pres)/150) 

# look for bias
train.bias <- train.out %>%
  group_by(species_name,pres,train) %>%
  summarise(pred_mean = mean(rf_pred)) %>%
  left_join(pres.prev, by="species_name")

train.bias %>%
  ggplot(aes(x=prev,y=pred_mean, color=factor(species_name))) +
  facet_grid(train~pres) +
  geom_point() +
  theme_bw() +
  theme(legend.position = "none")

### standard assessment: obs v. pred, AUC, r2, correlation
test.confMat <- test.avg 

# obs v. pred
obs.pred <- test.confMat %>%
  ungroup() %>%
  dplyr::select(c(plot_id,pres,pred_mean)) %>%
  as.data.frame()

presence.absence.hist(obs.pred, legend.cex=1,N.bars=15, 
                      ylab="Number of Sites", xlab="Predicted Probability",
                      truncate.tallest = TRUE)

# overall model assessment
sdm::roc(x=test.confMat$pres,p=test.confMat$pred_mean)
ModelMetrics::auc(actual=test.confMat$pres,predicted=test.confMat$pred_mean)
summary(lm(pred_mean~pres, test.confMat))
cor(test.confMat$pres,test.confMat$pred_mean, method="spearman")

# model fit r2
test.confMat %>%
  ggplot(aes(x=pres,y=pred_mean)) +
  geom_point() +
  geom_abline(intercept=0,slope=1,col="red") +
  stat_cor(
    aes(label = paste(..rr.label..)),
    label.x.npc = "left"
  ) +
  stat_cor(
    method="spearman",cor.coef.name="rho",label.x.npc="center",p.digits=NA, label.sep=""
  ) +
  theme_bw()

### individual species runs: individual AUC values, stability across reps
# average individ species AUC
test.bysp <- test.confMat %>%
  group_by(species_name) %>%
  summarise(auc = ModelMetrics::auc(actual=pres,predicted=pred_mean),
            rho = cor(pres,pred_mean,method="spearman")) %>%
  left_join(pres.prev, by=c("species_name")) %>%
  mutate(pres = prev*150)

# compare AUC with prevalence, fairly stable
test.bysp %>%
  ggplot(aes(x=pres, y=auc)) +
  geom_point() +
  geom_smooth(method="loess",span=0.3) +
  theme_bw()

# write.csv(test.bysp,paste0("processed_data/understory_model/evaluation/species_evaluation_auc_",run_nbr,".csv"),row.names=FALSE)

# look at averages vs. prevalence
test.bysp %>%
  mutate(count=1) %>%
  group_by(pres) %>%
  summarise(auc=mean(auc),n=sum(count)) %>% 
  # filter(n>1) %>%
  ggplot(aes(x=pres,y=auc, size=n)) +
  geom_point() +
  theme_bw()

# # table with totals by species
# tots.sdm <- test.confMat %>%
#   group_by(species_name) %>%
#   summarise(across(c(pred_mean), ~sum(.))) %>%
#   left_join(test.bysp, by="species_name")

# stability across reps, calc coefficient of variation
test.byrun <- test.sub %>%
  group_by(species_name,run_nbr) %>%
  summarise(auc = ModelMetrics::auc(actual=pres,predicted=rf_pred),
            rho = cor(pres,rf_pred,method="spearman")) %>%
  ungroup() %>%
  group_by(species_name) %>%
  summarise(across(c(auc,rho), list(mean=mean, cv = ~sd(., na.rm=TRUE)/mean(., na.rm=TRUE)*100), na.rm=TRUE))

test.byrun %>%
  pivot_longer(c(auc_mean,auc_cv)) %>%
  ggplot(aes(x=1,y=value)) +
  facet_wrap(~name,scales="free") +
  geom_boxplot() +
  theme_bw()

### compare to plot-level values: richness, ellenberg indicators
test.rich <- test.confMat %>% 
  group_by(plot_id) %>%
  summarise(richness=sum(pres), richness_pred = sum(pred_mean)) 

# overall richness
test.rich %>%
  ggplot(aes(x=richness, y=richness_pred)) +
  geom_point() +
  geom_smooth(method="lm") +
  ylim(c(min(test.rich$richness,test.rich$richness_pred),max(test.rich$richness,test.rich$richness_pred))) +
  xlim(c(min(test.rich$richness,test.rich$richness_pred),max(test.rich$richness,test.rich$richness_pred))) +
  geom_abline(intercept=0,slope=1,col="red") +
  stat_cor(
    aes(label = paste(..rr.label..)),
    label.x.npc = "left"
  ) +
  stat_cor(
    method="spearman",cor.coef.name="rho",label.x.npc="center",p.digits=NA, label.sep=""
  ) +
  theme_bw()

# # does excluding poorer models help in predicting richness? yes...
# temp <- test.bysp %>%
#   # filter(pres>=8)
#   filter(auc>0.7)
# 
# test.confMat %>% 
#   filter(species_name %in% c(temp$species_name)) %>%
#   group_by(plot_id) %>%
#   summarise(richness=sum(pres), richness_pred = sum(pred_mean)) %>%
#   ggplot(aes(x=richness, y=richness_pred)) +
#   geom_point() +
#   geom_smooth(method="lm") +
#   geom_abline(intercept=0,slope=1,col="red") +
#   stat_cor(
#     aes(label = paste(..rr.label..)),
#     label.x.npc = "left"
#   ) +
#   stat_cor(
#     method="spearman",cor.coef.name="rho",label.x.npc="center",p.digits=NA, label.sep=""
#   ) +
#   theme_bw()

# ellenberg indicator values
# prep observed data
ellen.comp <- test.confMat %>%
  left_join(pres.in, by=c("plot_id","species_name","pres")) %>%
  ungroup() %>%
  ellen_compute(., unique(test.confMat$species_name))

# prep test predictions
ellen.test <- test.confMat %>%
  ungroup() %>%
  rename(cover=pred_mean) %>%
  ellen_compute(., unique(test.confMat$species_name)) %>%
  rename(pred_mean = wtmean) %>%
  dplyr::select(-mean) %>%
  left_join(ellen.comp, by=c("plot_id","OriglName"))

# ellen overall
ellen.test %>%
  filter(OriglName %in% c("F","Licht","N","T")) %>%
  ggplot(aes(x=mean, y=pred_mean)) +
  facet_wrap(~OriglName, scales="free") +
  geom_point() +
  geom_smooth(method="lm") +
  geom_abline(intercept=0,slope=1,col="red") +
  stat_cor(
    aes(label = paste(..rr.label..)),
    label.x.npc = "left"
  ) +
  stat_cor(
    method="spearman",cor.coef.name="rho",label.x.npc="center",p.digits=NA, label.sep=""
  ) +
  theme_bw()

write.csv(ellen.test, paste0("processed_data/understory_model/evaluation/corrected_ellenberg_",rich.flag,"_evaluation_",run_nbr,prr,".csv"),row.names=FALSE)

### evaluate impact of applying the probability ranking rule for individual species, plot-level predictions
prr <- "prr"

# prr, test data only
train.out <- read.csv(paste0("processed_data/understory_model/evaluation/sdm_corrprr_",rich.flag,"_evaluation_",run_nbr,".csv")) %>%
  left_join(pres.in, by=c("plot_id","species_name")) %>%
  dplyr::select(-cover) %>%
  rename(rf_pred=prr_pres,
         run_nbr=run)

### look at confusion matrix and thresholds for test data
confusionMatrix(data=as.factor(train.out$rf_pred),reference=as.factor(train.out$pres))

### richness and ellenberg
test.confMat <- train.out %>%
  group_by(plot_id,pres,species_name) %>%
  summarise(pred_mean = mean(rf_pred)) 

test.rich <- test.confMat %>% 
  group_by(plot_id) %>%
  summarise(richness=sum(pres), richness_pred = sum(pred_mean)) 

# overall richness
test.rich %>%
  ggplot(aes(x=richness, y=richness_pred)) +
  geom_point() +
  geom_smooth(method="lm") +
  ylim(c(min(test.rich$richness,test.rich$richness_pred),max(test.rich$richness,test.rich$richness_pred))) +
  xlim(c(min(test.rich$richness,test.rich$richness_pred),max(test.rich$richness,test.rich$richness_pred))) +
  geom_abline(intercept=0,slope=1,col="red") +
  stat_cor(
    aes(label = paste(..rr.label..)),
    label.x.npc = "left"
  ) +
  stat_cor(
    method="spearman",cor.coef.name="rho",label.x.npc="center",p.digits=NA, label.sep=""
  ) +
  theme_bw()

# ellenberg indicator values
# prep observed data
ellen.comp <- test.confMat %>%
  left_join(pres.in, by=c("plot_id","species_name","pres")) %>%
  ungroup() %>%
  ellen_compute(., unique(test.confMat$species_name))

# prep test predictions
ellen.test <- test.confMat %>%
  ungroup() %>%
  rename(cover=pred_mean) %>%
  ellen_compute(., unique(test.confMat$species_name)) %>%
  rename(pred_mean = wtmean) %>%
  dplyr::select(-mean) %>%
  left_join(ellen.comp, by=c("plot_id","OriglName"))

# ellen overall
ellen.test %>%
  filter(OriglName %in% c("F","Licht","N","T")) %>%
  ggplot(aes(x=mean, y=pred_mean)) +
  facet_wrap(~OriglName, scales="free") +
  geom_point() +
  geom_smooth(method="lm") +
  geom_abline(intercept=0,slope=1,col="red") +
  stat_cor(
    aes(label = paste(..rr.label..)),
    label.x.npc = "left"
  ) +
  stat_cor(
    method="spearman",cor.coef.name="rho",label.x.npc="center",p.digits=NA, label.sep=""
  ) +
  theme_bw()

write.csv(ellen.test, paste0("processed_data/understory_model/evaluation/corrected_ellenberg_",rich.flag,"_evaluation_",run_nbr,prr,".csv"),row.names=FALSE)

###
# 6. final fits: variable importance and response curves
###

# read in
rcurves.out <- read.csv("processed_data/understory_model/evaluation/sdm_pres_finalFit_rcurves_9pred.csv") 
vimp.out <- read.csv("processed_data/understory_model/evaluation/sdm_pres_finalFit_vimp_9pred.csv")

### assessment for core group: species with 5 or more observations
rcurves.sub <- rcurves.out %>%
  filter(species_name %in% c(pres.sub$species_name))
vimp.sub <- vimp.out %>%
  filter(species_name %in% c(pres.sub$species_name))

### variable importance
vimp.sub %>%
  rename(IncMSE = X.IncMSE) %>%
  # set negative MSE values to 0, so don't artificially deflate importance
  mutate(IncMSE = ifelse(IncMSE <0,0,IncMSE)) %>%
  group_by(variables) %>%
  summarise(IncMSE = mean(IncMSE)) %>%
  mutate(variables = factor(variables)) %>%
  ggplot(aes(x=reorder(variables,-desc(IncMSE)),y=IncMSE)) +#, fill=species_name)) +
  geom_bar(stat='identity') +
  coord_flip() +
  ggtitle("Variable importance for presence") +
  theme_bw() +
  theme(
    legend.position = "none",
    axis.title.y = element_blank())

### does model align with ecological expectations
# response curves, match up with pft categories
rcurves.summ <- rcurves.sub %>%
  left_join(pft.in, by="species_name")

# light-loving, cold-habitat species
rcurves.summ %>%
  group_by(variable,x) %>%
  filter(pft=="light-cold") %>%
  summarise(mean=mean(y), p25=quantile(y,0.25),p75=quantile(y,0.75)) %>%
  ggplot(aes(x=x,y=mean)) +
  facet_wrap(~variable, scales="free_x") +
  # geom_ribbon(aes(ymin=p25,ymax=p75),alpha=0.2) +
  geom_line(color="blue", size=1) +
  ggtitle("Avg response: Light-loving, subalpine-alpine species (56)") +
  theme_bw()

# shade-tolerant, warm-habitat species
rcurves.summ %>%
  group_by(variable,x) %>%
  filter(pft=="shade-warm") %>%
  summarise(mean=mean(y), ci_low=quantile(y,0.25),p75=quantile(y,0.75)) %>%
  ggplot(aes(x=x,y=mean)) +
  facet_wrap(~variable, scales="free_x") +
  # geom_ribbon(aes(ymin=p25,ymax=p75),alpha=0.2) +
  geom_line(color="blue", size=1) +
  ggtitle("Avg response: Shade-tolerant, submontane-montane species (54)") +
  theme_bw()

# light-loving, temperature indifferent
rcurves.summ %>%
  group_by(variable,x) %>%
  filter(pft=="light-notemp") %>%
  summarise(mean=mean(y), p25=quantile(y,0.25),p75=quantile(y,0.75)) %>%
  ggplot(aes(x=x,y=mean)) +
  facet_wrap(~variable, scales="free_x") +
  # geom_ribbon(aes(ymin=p25,ymax=p75),alpha=0.2) +
  geom_line(color="blue", size=1) +
  ggtitle("Avg response: Light-loving, temperature-indifferent species (60)") +
  theme_bw()

# light-loving, warm-habitat
rcurves.summ %>%
  group_by(variable,x) %>%
  filter(pft=="light-warm") %>%
  summarise(mean=mean(y), p25=quantile(y,0.25),p75=quantile(y,0.75)) %>%
  ggplot(aes(x=x,y=mean)) +
  facet_wrap(~variable, scales="free_x") +
  # geom_ribbon(aes(ymin=p25,ymax=p75),alpha=0.2) +
  geom_line(color="blue", size=1) +
  ggtitle("Avg response: Light-loving, submontane-montane species (31)") +
  theme_bw()

# shade-tolerant, cold-habitat
rcurves.summ %>%
  group_by(variable,x) %>%
  filter(pft=="shade-cold") %>%
  summarise(mean=mean(y), p25=quantile(y,0.25),p75=quantile(y,0.75)) %>%
  ggplot(aes(x=x,y=mean)) +
  facet_wrap(~variable, scales="free_x") +
  # geom_ribbon(aes(ymin=p25,ymax=p75),alpha=0.2) +
  geom_line(color="blue", size=1) +
  ggtitle("Avg response: Shade-tolerant, subalpine-alpine species (12)") +
  theme_bw()

# shade-tolerant, temperature indifferent
rcurves.summ %>%
  group_by(variable,x) %>%
  filter(pft=="shade-notemp") %>%
  summarise(mean=mean(y), p25=quantile(y,0.25),p75=quantile(y,0.75)) %>%
  ggplot(aes(x=x,y=mean)) +
  facet_wrap(~variable, scales="free_x") +
  # geom_ribbon(aes(ymin=p25,ymax=p75),alpha=0.2) +
  geom_line(color="blue", size=1) +
  ggtitle("Avg response: Shade-tolerant, temperature-indifferent species (35)") +
  theme_bw()
