#####
# 
## Step 10: Data analysis
#
#####

# if running on server
# write("TMPDIR='/mnt/public/temp/kbraziunas'", file=file.path("~/.Renviron")) 

### load libraries

library(tidyverse)
library(plotrix) # std.error
library(car) # qq plots
library(lme4) # lmer
library(nlme) # lme with unequal variance
library(MuMIn) # marginal and conditional r2

### wd is one level up, recommend set up Rproj or:
# setwd("../")

### selected sessionInfo()
# R version 4.1.3 (2022-03-10)
# Platform: x86_64-w64-mingw32/x64 (64-bit)
# Running under: Windows 10 x64 (build 19044)
# 
# Matrix products: default
# 
# locale:
#   [1] LC_COLLATE=English_United States.1252  LC_CTYPE=English_United States.1252    LC_MONETARY=English_United States.1252
# [4] LC_NUMERIC=C                           LC_TIME=English_United States.1252    
# 
# attached base packages:
#   [1] stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#   [1] MuMIn_1.46.0    nlme_3.1-155    lme4_1.1-28     Matrix_1.4-0    car_3.0-12      carData_3.0-5   plotrix_3.8-2  
# [8] forcats_0.5.1   stringr_1.4.1   dplyr_1.0.8     purrr_0.3.4     readr_2.1.2     tidyr_1.2.0     tibble_3.1.6   
# [15] ggplot2_3.3.5   tidyverse_1.3.1

###
# 1. Q1: Future forest floor, understory plant community change in response to 21st-century changes in climate and forests
###

### read in data
# generate for either full species subset or auc>0.7 cutoff
rich.flag <- "subset"
# rich.flag <- "auc07"

# species list
pres.sub <- read.csv(paste0("processed_data/understory_model/understory_",rich.flag,"_presence_cover.csv"))

# richness and T: hx
resp.hx <- setNames(
  list.files(paste0("processed_data/iland_sdm_predictions/historical_predictions_10m/",rich.flag,"/"), pattern="_plot_", full.names=TRUE), 
  list.files(paste0("processed_data/iland_sdm_predictions/historical_predictions_10m/",rich.flag,"/"), pattern="_plot_")) %>%
  map_dfr(read.csv, .id = NULL)

# cover: hx
cover.hx <- read.csv("processed_data/iland_sdm_predictions/historical_predictions_10m/historical_total_cover.csv")

# hx drivers
pred.hx <- read.csv("processed_data/iland_sdm_predictions/sdm_predictor_set_nospecies_notstd_10m.csv") %>%
  rename(id_model=cell_id) %>%
  dplyr::select(c(id_model,mean_temp,summer_prec,rad,BA,prop_fasy,lif))

# future responses
resp.fut <- read.csv(paste0("processed_data/iland_sdm_predictions/future_prediction_summaries/scenario_mean_cv_",rich.flag,".csv"))

# lost species
lost.fut <- read.csv("processed_data/iland_sdm_predictions/future_prediction_summaries/scenario_species_lost.csv") %>%
  filter(rich_flag==rich.flag)

# elevation lookup
elev.in <- read.csv("processed_data/iland_sdm_predictions/id_elev_lookup.csv")

# random sample ids
samp.in <- read.csv("processed_data/iland_sdm_predictions/rand_sample_points.csv")

# future responses: spatial
sp.in <- data.frame(read.csv(paste0("processed_data/iland_sdm_predictions/future_prediction_summaries/spatial_response_change_rcp85_",rich.flag,".csv")), rcp="rcp85") %>%
  rbind(data.frame(read.csv(paste0("processed_data/iland_sdm_predictions/future_prediction_summaries/spatial_response_change_rcp45_",rich.flag,".csv")), rcp="rcp45"))

# red list lookup
redlist.in <- read.csv("processed_data/species_lookup/red_list_status.csv")

### 
# average change across the full landscape, by scenario
###
# average plot-level changes in: richness, thermophilization, cover
# also average species lost, species no longer present anywhere in the landscape

# calculate changes relative to mean values
# cover only included with subset
if(rich.flag=="subset") {
  resp.change <- resp.fut %>%
    # selection
    dplyr::select(c(pred_decade:rep,mean_cover_pct,mean_richness,mean_T)) %>%
    rename(cover=mean_cover_pct,richness=mean_richness,eivTemp=mean_T) %>%
    # pivot
    pivot_longer(c(cover:eivTemp)) %>%
    # add historical values
    mutate(historical = ifelse(name=="cover",mean(cover.hx$cover_pct),
                               ifelse(name=="richness",mean(resp.hx$corr_richness),
                                      ifelse(name=="eivTemp",mean(resp.hx$corr_T),NA)))) %>%
    # calculate change
    mutate(delta_change = value-historical,
           rel_change= (delta_change/historical)*100) 
} else if(rich.flag=="auc07") {
  resp.change <- resp.fut %>%
    # selection
    dplyr::select(c(pred_decade:rep,mean_richness,mean_T)) %>%
    rename(richness=mean_richness,eivTemp=mean_T) %>%
    # pivot
    pivot_longer(c(richness:eivTemp)) %>%
    # add historical values
    mutate(historical = ifelse(name=="richness",mean(resp.hx$corr_richness),
                               ifelse(name=="eivTemp",mean(resp.hx$corr_T),NA))) %>%
    # calculate change
    mutate(delta_change = value-historical,
           rel_change= (delta_change/historical)*100)
}

# calculate species lost
land.lost <- lost.fut %>%
  filter(pred_decade != 2020) %>%
  rename(rep=pred_rep) %>%
  group_by(pred_decade,pred_gcm,pred_wind,rep) %>%
  summarise(delta_change=-sum(lost)) %>%
  mutate(name="lost",
         rel_change = (delta_change/length(unique(pres.sub$species_name)))*100) # divide by orig number of species

# percentiles across all scenarios
resp.pctile <- resp.change %>%
  # include land.lost
  dplyr::select(-c(forest_climate,pred_forest,pred_climate,historical,pred_reps,value)) %>% 
  rbind(land.lost) %>%
  group_by(pred_decade,name) %>%
  summarise(across(c(delta_change,rel_change), ~quantile(., c(0,0.05,0.25,0.5,0.75,0.95,1))),q=c("min","q5","q25","median","q75","q95","max")) %>%
  mutate(pred_gcm="all",
         pred_wind="all") %>%
  pivot_wider(names_from="q",values_from=c("delta_change","rel_change"))

# percentiles by disturbance scenario
dist.pctile <- resp.change %>%
  # include land.lost
  dplyr::select(-c(forest_climate,pred_forest,pred_climate,historical,pred_reps,value)) %>% 
  rbind(land.lost) %>%
  group_by(pred_decade,pred_wind,name) %>%
  summarise(across(c(delta_change,rel_change), ~quantile(., c(0,0.05,0.25,0.5,0.75,0.95,1))),q=c("min","q5","q25","median","q75","q95","max")) %>%
  mutate(pred_gcm="all",
         pred_wind = as.character(pred_wind)) %>%
  pivot_wider(names_from="q",values_from=c("delta_change","rel_change"))

# percentiles by each climate x forest scenario
scen.pctile <- resp.change %>%
  # include land.lost
  dplyr::select(-c(forest_climate,pred_forest,pred_climate,historical,pred_reps,value)) %>%
  rbind(land.lost) %>%
  group_by(pred_gcm,pred_wind,pred_decade,name) %>%
  summarise(across(c(delta_change,rel_change), ~quantile(., c(0,0.05,0.25,0.5,0.75,0.95,1))),q=c("min","q5","q25","median","q75","q95","max")) %>%
  mutate(pred_wind = as.character(pred_wind)) %>%
  pivot_wider(names_from="q",values_from=c("delta_change","rel_change"))

# combine and write out 
pctile.out <- resp.pctile %>%
  rbind(dist.pctile) %>%
  rbind(scen.pctile) 

write.csv(pctile.out, paste0("analysis/q1_futureforestfloor/understory_future_change_",rich.flag,".csv"),row.names=FALSE)

# quick plots
pctile.out %>%
  filter(pred_wind=="all",pred_gcm=="all") %>%
  ggplot(aes(x=factor(pred_decade))) +
  facet_wrap(~name, scales="free_y") +
  geom_hline(aes(yintercept=0),lty=2,color="#555555") +
  geom_pointrange(aes(y=delta_change_median, ymin=delta_change_q5,ymax=delta_change_q95),position=position_dodge(width=0.5)) +
  ylab("Actual change") +
  ggtitle("Change in response (5th-95th %iles)") +
  theme_bw() +
  theme(
    axis.title.x = element_blank(),
    panel.grid.major.x = element_blank()
  )

pctile.out %>%
  filter(pred_wind=="all",pred_gcm=="all") %>%
  ggplot(aes(x=factor(pred_decade))) +
  facet_wrap(~name) +
  geom_hline(aes(yintercept=0),lty=2,color="#555555") +
  geom_pointrange(aes(y=rel_change_median, ymin=rel_change_q5,ymax=rel_change_q95),position=position_dodge(width=0.5)) +
  ylab("Relative change (%)") +
  ggtitle("Change in response (5th-95th %iles)") +
  theme_bw() +
  theme(
    axis.title.x = element_blank(),
    panel.grid.major.x = element_blank()
  )

pctile.out %>%
  filter(!pred_wind=="all",pred_gcm=="all") %>%
  ggplot(aes(x=factor(pred_decade), color=factor(pred_wind))) +
  facet_wrap(~name) +
  geom_hline(aes(yintercept=0),lty=2,color="#555555") +
  geom_pointrange(aes(y=rel_change_median, ymin=rel_change_q5,ymax=rel_change_q95),position=position_dodge(width=0.5)) +
  ylab("Relative change (%)") +
  ggtitle("Change in response (5th-95th %iles)") +
  theme_bw() +
  theme(
    axis.title.x = element_blank(),
    panel.grid.major.x = element_blank()
  )

pctile.out %>%
  filter(!pred_wind=="all",!grepl("all",pred_gcm)) %>%
  mutate(scen = paste(pred_gcm,pred_wind,sep="_")) %>%
  ggplot(aes(x=factor(pred_decade), color=factor(scen))) +
  facet_wrap(~name) +
  geom_hline(aes(yintercept=0),lty=2,color="#555555") +
  geom_pointrange(aes(y=rel_change_median, ymin=rel_change_q5,ymax=rel_change_q95),position=position_dodge(width=0.5)) +
  ylab("Relative change (%)") +
  ggtitle("Change in response (5th-95th %iles)") +
  theme_bw() +
  theme(
    axis.title.x = element_blank(),
    panel.grid.major.x = element_blank()
  )

### 
# local-scale variability
###

### NB: this chunk takes a long time to run, especially calculating the pairwise Spearman's correlation. can skip below (to NB: end computationally intensive chunk) to load results

# first: assess whether changes in richness, thermophilization, and cover are consistent in direction and magnitude; assess separately for 2050 and 2100

dec.in <- 2050
dec.iland <- 30

# dec.in <- 2100
# dec.iland <- 80

# read in future plot-level values and calculate change, takes some time
cover.fut <- setNames(
  list.files("processed_data/iland_sdm_predictions/future_predictions_10m/", pattern=paste(dec.in,"total_cover",sep="_"), full.names=TRUE),
  list.files("processed_data/iland_sdm_predictions/future_predictions_10m/", pattern=paste(dec.in,"total_cover",sep="_"))) %>%
  map_dfr(read.csv, .id = "run_id") %>%
  rename(fut_cover_pct = cover_pct)

cover.names <- data.frame(run_id = unique(cover.fut$run_id),run_id2 = unique(cover.fut$run_id)) %>%
  separate(run_id2, into=c("pref1","pref2","gcm","rcp","gcm2","gcm3","pred_wind","pred_rep","pred_decade","suff1","suff2"), sep="_") %>%
  dplyr::select(c(run_id,gcm,rcp,pred_wind,pred_rep,pred_decade))

cover.spchange <- cover.fut %>%
  left_join(cover.names, by="run_id") %>%
  left_join(cover.hx, by="id_model") %>%
  mutate(cover_change = fut_cover_pct - cover_pct) %>%
  dplyr::select(c(gcm,rcp,pred_wind,pred_rep,pred_decade,id_model,cover_change))

# clean
rm(cover.fut)
gc()

resp.fut <- setNames(
  list.files(paste0("processed_data/iland_sdm_predictions/future_predictions_10m/",rich.flag,"/"), pattern=paste(dec.in,"plot",sep="_"), full.names=TRUE),
  list.files(paste0("processed_data/iland_sdm_predictions/future_predictions_10m/",rich.flag,"/"), pattern=paste(dec.in,"plot",sep="_"))) %>%
  map_dfr(read.csv, .id = "run_id") %>%
  rename(fut_richness=corr_richness,
         fut_T = corr_T) %>%
  dplyr::select(c(run_id,id_model,fut_richness,fut_T))

resp.names <- data.frame(run_id = unique(resp.fut$run_id),run_id2 = unique(resp.fut$run_id)) %>%
  separate(run_id2, into=c("pref1","pref2","gcm","rcp","gcm2","gcm3","pred_wind","pred_rep","pred_decade","suff1","suff2","suff3"), sep="_") %>%
  dplyr::select(c(run_id,gcm,rcp,pred_wind,pred_rep,pred_decade))

resp.spchange <- resp.hx %>%
  dplyr::select(c(id_model,corr_richness,corr_T)) %>%
  left_join(resp.fut, by="id_model") %>%
  left_join(resp.names, by="run_id") %>%
  mutate(richness_change = fut_richness - corr_richness,
         T_change = fut_T - corr_T) %>%
  dplyr::select(c(gcm,rcp,pred_wind,pred_rep,pred_decade,id_model,richness_change,T_change))

# clean
rm(resp.fut)
gc()

# match these up
full.spchange1 <- cover.spchange %>%
  left_join(resp.spchange, by=c("gcm","rcp","pred_wind","pred_rep","pred_decade","id_model"))
# clean up
rm(cover.spchange, resp.spchange)
gc()

# read in forest drivers
for.fut <- setNames(
  list.files("iland/output/combined_outputs/", pattern="forest_predictors", full.names=TRUE),
  list.files("iland/output/combined_outputs/",pattern="forest_predictors")) %>%
  map_dfr(read.csv, .id = "run_id") %>%
  # subset to decade of interest
  filter(year==dec.iland) %>%
  rename(fut_BA=BA,
         fut_prop_fasy = prop_fasy,
         fut_lif = lif,
         id_model=cell_id) %>%
  dplyr::select(c(run_id,id_model,fut_BA,fut_prop_fasy,fut_lif))

for.names <- data.frame(run_id = unique(for.fut$run_id),run_id2 = unique(for.fut$run_id)) %>%
  separate(run_id2, into=c("gcm","rcp","gcm2","gcm3","pred_wind","pred_rep","suff1","suff2"), sep="_") %>%
  dplyr::select(c(run_id,gcm,rcp,pred_wind,pred_rep))

for.spchange <- full.spchange1 %>%
  # first have to add 0s where values are NA
  dplyr::select(c(gcm,rcp,pred_wind,pred_rep,id_model)) %>%
  left_join(for.names, by=c("gcm","rcp","pred_wind","pred_rep")) %>%
  left_join(for.fut, by=c("run_id","id_model")) %>%
  # assign 0s to NA forest structure values
  mutate(fut_BA = replace_na(fut_BA,0)) %>%
  # assign 1 to lif when no forest
  mutate(fut_lif = replace_na(fut_lif, 1)) %>%
  left_join(pred.hx, by="id_model") %>%
  dplyr::select(-c(mean_temp,summer_prec,rad)) %>%
  mutate(BA_change = fut_BA-BA,
         lif_change = fut_lif-lif,
         prop_fasy_change = fut_prop_fasy-prop_fasy) %>%
  dplyr::select(c(gcm,rcp,pred_wind,pred_rep,id_model,BA_change,lif_change,prop_fasy_change))

# clean
rm(for.fut)
gc()

# match everything up
full.spchange <- full.spchange1 %>%
  left_join(for.spchange, by=c("gcm","rcp","pred_wind","pred_rep","id_model")) %>%
  left_join(elev.in, by="id_model")

# clean up
rm(full.spchange1, for.spchange)
gc()

# spearman's rho correlations for pairwise changes
pairwise_spearman <- data.frame(pred_decade=dec.in)
pairwise_spearman$richness_T <- cor(full.spchange$richness_change,full.spchange$T_change, method="spearman")
pairwise_spearman$richness_cover <- cor(full.spchange$richness_change,full.spchange$cover_change, method="spearman")
pairwise_spearman$cover_T <- cor(full.spchange$cover_change,full.spchange$T_change, method="spearman")
# cover and richness v. BA and lif
pairwise_spearman$richness_BA <- cor(full.spchange$richness_change,full.spchange$BA_change, method="spearman")
pairwise_spearman$richness_lif <- cor(full.spchange$richness_change,full.spchange$lif_change, method="spearman")
pairwise_spearman$cover_BA <- cor(full.spchange$cover_change,full.spchange$BA_change, method="spearman")
pairwise_spearman$cover_lif <- cor(full.spchange$cover_change,full.spchange$lif_change, method="spearman")
# all v. elevation
pairwise_spearman$richness_elev <- cor(full.spchange$richness_change,full.spchange$elev_m, method="spearman")
pairwise_spearman$T_elev <- cor(full.spchange$T_change,full.spchange$elev_m, method="spearman")
pairwise_spearman$cover_elev <- cor(full.spchange$cover_change,full.spchange$elev_m, method="spearman")
pairwise_spearman$BA_elev <- cor(full.spchange$BA_change,full.spchange$elev_m, method="spearman")
pairwise_spearman$lif_elev <- cor(full.spchange$lif_change,full.spchange$elev_m, method="spearman")

# remove nas for prop_fasy
sub.spchange <- full.spchange %>%
  filter(!is.na(prop_fasy_change))

pairwise_spearman$propFasy_elev <- cor(sub.spchange$prop_fasy_change,sub.spchange$elev_m, method="spearman")

# add n_obs
pairwise_spearman$n_obs <- dim(full.spchange)[1]
pairwise_spearman$n_obs_prop_fasy <- dim(sub.spchange)[1]

# can save subset for plotting purposes, use random sample of 1000 points
samp.spchange <- full.spchange %>%
  filter(id_model %in% c(samp.in$id_model))
# save
write.csv(samp.spchange, paste0("processed_data/iland_sdm_predictions/future_prediction_summaries/rand_sample_response_forest_",rich.flag,"_",dec.in,".csv"),row.names=FALSE)

# clean
rm(full.spchange,sub.spchange)
gc()

# also read in climate and match up with elevation, gcm column also includes rcp
clim.fut <- read.csv("processed_data/iland_clim_soils/future_predictors/sdm_predictor_set_climate_rcp45_10m.csv") %>%
  rbind(read.csv("processed_data/iland_clim_soils/future_predictors/sdm_predictor_set_climate_rcp85_10m.csv")) %>%
  dplyr::select(c(id_model,gcm,decade,mean_temp,summer_prec,rad)) %>%
  filter(decade==dec.in) %>%
  rename(fut_mean_temp=mean_temp,
         fut_summer_prec=summer_prec,
         fut_rad=rad) %>%
  left_join(pred.hx, by=c("id_model")) %>%
  dplyr::select(-c(BA,prop_fasy,lif)) %>%
  mutate(mean_temp_change = fut_mean_temp-mean_temp,
         summer_prec_change = fut_summer_prec-summer_prec,
         rad_change = fut_rad-rad) %>%
  dplyr::select(c(gcm,decade,id_model,mean_temp_change,summer_prec_change,rad_change)) %>%
  left_join(elev.in, by="id_model") 

clim.spear <- clim.fut %>%
  # subset to unique entries, climate cells are coarser resolution
  group_by(gcm,decade,mean_temp_change,summer_prec_change,rad_change,elev_m) %>%
  slice(1)

# add spearman correlations
pairwise_spearman$meanTemp_elev <- cor(clim.spear$mean_temp_change,clim.spear$elev_m, method="spearman")
pairwise_spearman$summerPrec_elev <- cor(clim.spear$summer_prec_change,clim.spear$elev_m, method="spearman")
pairwise_spearman$rad_elev <- cor(clim.spear$rad_change,clim.spear$elev_m, method="spearman")

pairwise_spearman$n_obs_clim <- dim(clim.spear)[1]

samp.climchange <- clim.fut %>%
  filter(id_model %in% c(samp.in$id_model))

# write out final correlations
write.csv(pairwise_spearman, paste0("analysis/q1_futureforestfloor/pairwise_spearman_correlations_",rich.flag,"_",dec.in,".csv"),row.names=FALSE)

# can save random sample for plotting purposes
write.csv(samp.climchange, paste0("processed_data/iland_sdm_predictions/future_prediction_summaries/rand_sample_clim_",rich.flag,"_",dec.in,".csv"),row.names=FALSE)

# clean up
rm(clim.fut)
gc()

### NB: end computationally intensive chunk for calculating Spearman's rank correlations

### combine spearman outputs from above
pairwise_spearman <- read.csv(paste0("analysis/q1_futureforestfloor/pairwise_spearman_correlations_",rich.flag,"_2050.csv")) %>%
  rbind(read.csv(paste0("analysis/q1_futureforestfloor/pairwise_spearman_correlations_",rich.flag,"_2100.csv")))

# relationships because response drivers, relationships between response and climate with elevation
# CHANGE in richness and T weakly to moderately neg correlated; richness and cover moderately pos correlated; cover and T only very weakly correlated; correlation strength of all pairs is higher in 2100 than 2050
# CHANGE in richness is weakly correlated with change in BA and lif; areas that decrease in BA or increase in lif are associated with increases in richness
# CHANGE in cover is strongly correlated with change in BA and lif; areas that decrease in BA or increase in lif are associated with increases in cover
# change in responses and forest drivers mostly very weakly correlated with elevation in 2050; however, in 2100, relationships emerge with higher elevation weakly correlated with decreasing richness, moderately correlated with increasing T, moderately correlated with dereasing lif, moderately correlated with decreasing fasy (fasy increases at lower elevations)
# change in climate drivers similar relationships with elevation at 2 time periods. change in mean temp and summer precip both weakly correlated with increasing elevation (greater changes at higher elevations)
pairwise_spearman

### look at the shape of these relationships based on random sample
samp.spchange <- read.csv(paste0("processed_data/iland_sdm_predictions/future_prediction_summaries/rand_sample_response_forest_",rich.flag,"_2050.csv")) %>%
  rbind(read.csv(paste0("processed_data/iland_sdm_predictions/future_prediction_summaries/rand_sample_response_forest_",rich.flag,"_2100.csv")))

samp.climchange <- read.csv(paste0("processed_data/iland_sdm_predictions/future_prediction_summaries/rand_sample_clim_",rich.flag,"_2050.csv")) %>%
  rbind(read.csv(paste0("processed_data/iland_sdm_predictions/future_prediction_summaries/rand_sample_clim_",rich.flag,"_2100.csv")))

# pairwise relationships
samp_pairwise <- function(x.in,y.in) {
  samp.spchange %>%
    mutate(pred_decade=factor(pred_decade)) %>%
    ggplot(aes_string(x=x.in,y=y.in, color="pred_decade",fill="pred_decade")) +
    geom_hline(aes(yintercept=0),lty=2) +
    geom_vline(aes(xintercept=0),lty=2) +
    scale_color_manual(values=c("#969696","#000000")) +
    geom_smooth(method="loess", se=FALSE) +
    theme_bw()
}

samp_pairwise("T_change","richness_change")
samp_pairwise("richness_change","cover_change")
samp_pairwise("BA_change","richness_change")
samp_pairwise("BA_change","cover_change")
samp_pairwise("lif_change","richness_change")
samp_pairwise("lif_change","cover_change")

# change in responses and forest drivers with elevation
samp.spchange %>%
  pivot_longer(c(cover_change:prop_fasy_change)) %>%
  # pivot_longer(c(T_change)) %>%
  mutate(pred_decade=factor(pred_decade)) %>%
  ggplot(aes(x=elev_m,y=value, color=pred_decade,fill=pred_decade)) +
  facet_wrap(~name, scales="free_y") +
  geom_hline(aes(yintercept=0),lty=2) +
  scale_color_manual(values=c("#969696","#000000")) +
  scale_fill_manual(values=c("#969696","#000000")) +
  geom_smooth(method="loess", se=FALSE) +
  theme_bw()

# change in climate drivers
samp.climchange %>%
  pivot_longer(c(mean_temp_change:rad_change)) %>%
  mutate(decade=factor(decade)) %>%
  ggplot(aes(x=elev_m,y=value, color=decade,fill=decade)) +
  facet_wrap(~name, scales="free_y") +
  geom_hline(aes(yintercept=0),lty=2) +
  scale_color_manual(values=c("#969696","#000000")) +
  scale_fill_manual(values=c("#969696","#000000")) +
  # geom_point(alpha=0.2) +
  # fewer obs for climate, turn on SE
  geom_smooth(method="loess", se=FALSE) +
  theme_bw()

### group richness and thermophilization responses by elevation zone
head(resp.hx) # historical
head(sp.in) # future

# hx mean by elev group
hx.elev <- resp.hx %>%
  dplyr::select(c(id_model,corr_richness,corr_T)) %>%
  left_join(elev.in, by="id_model") %>%
  mutate(elev_group = ifelse(elev_m<850,"submontane",
                             ifelse(elev_m>850 & elev_m<1400,"montane",
                                    ifelse(elev_m>1400,"subalpine",NA)))) %>%
  group_by(elev_group) %>%
  summarise(mean_richness = mean(corr_richness),
            mean_T = mean(corr_T)) %>%
  pivot_longer(c(mean_richness,mean_T), values_to="hx_mean")

# future mean by elev group
sp.elev <- sp.in %>%
  dplyr::select(c(id_model,pred_decade,mean_richness,mean_T)) %>%
  left_join(elev.in, by="id_model") %>%
  mutate(elev_group = ifelse(elev_m<850,"submontane",
                             ifelse(elev_m>850 & elev_m<1400,"montane",
                                    ifelse(elev_m>1400,"subalpine",NA)))) %>%
  group_by(elev_group,pred_decade) %>%
  summarise(mean_richness = mean(mean_richness),
            mean_T = mean(mean_T)) %>%
  pivot_longer(c(mean_richness,mean_T), values_to="fut_mean")

# difference by elevation group
sp.elev %>%
  left_join(hx.elev, by=c("elev_group","name")) %>%
  mutate(diff = fut_mean-hx_mean,
         rel_diff = diff/hx_mean * 100) 
  

###	
# Shifts in plant community composition with 21st century change
###
# look at shifts in PFTs (winners/losers)
# in life forms (winners/losers)
# look at turnover (# retained, # new, # lost)

# use average future spatial responses
head(sp.in)

sp.in %>%
  group_by(pred_decade) %>%
  summarise(turnover=mean(mean_turnover)*100)

# summarise by elevation
sp.summ <- sp.in %>%
  left_join(elev.in, by="id_model") %>%
  mutate(elev_group = ifelse(elev_m<850,"submontane",
                             ifelse(elev_m>850 & elev_m<1400,"montane",
                                    ifelse(elev_m>1400,"subalpine",NA)))) %>%
  # calculate proportions
  mutate(across(c(light_cold:shrub),~(.)/mean_richness),
         # calculate turnover, multipy by 100 to get %
         mean_turnover = mean_turnover*100,
         turnover_lost = lost/(lost+new+retained)*100,
         turnover_new = new/(lost+new+retained)*100) %>%
  # remove unused columns
  dplyr::select(-c(cv_richness,cv_T,cv_turnover,lost:retained)) %>%
  # mean and variance across space
  pivot_longer(c(mean_richness:shrub,turnover_lost,turnover_new)) %>%
  # add count
  mutate(n_obs = 1) %>%
  group_by(pred_decade,elev_group,name) %>%
  summarise(mean=mean(value), sd=sd(value), n_obs = sum(n_obs))

# check sums
sp.summ %>%
  filter(name %in% c("fern","graminoid","herb","shrub")) %>%
  ungroup() %>%
  group_by(pred_decade,elev_group) %>%
  summarise(sum=sum(mean)) %>%
  summary()
  
sp.summ %>%
  filter(name %in% c("light_cold","light_warm","light_notemp","shade_cold","shade_warm","shade_notemp")) %>%
  ungroup() %>%
  group_by(pred_decade,elev_group) %>%
  summarise(sum=sum(mean)) %>%
  summary()

sp.summ %>%
  filter(name %in% c("mean_turnover","turnover_lost","turnover_new")) %>%
  dplyr::select(-c(sd)) %>%
  ungroup() %>%
  pivot_wider(names_from="name", values_from="mean") %>%
  group_by(pred_decade,elev_group) %>%
  summarise(diff = mean_turnover - turnover_lost - turnover_new) %>%
  summary() # slight difference, ok for plotting

# add historical
sphx.summ <- resp.hx %>%
  left_join(elev.in, by="id_model") %>%
  mutate(elev_group = ifelse(elev_m<850,"submontane",
                             ifelse(elev_m>850 & elev_m<1400,"montane",
                                    ifelse(elev_m>1400,"subalpine",NA)))) %>%
  # rename
  rename(mean_richness=corr_richness) %>%
  # calculate proportions
  mutate(across(c(light_cold:shrub),~(.)/mean_richness),
         pred_decade=2020) %>%
  # remove unused columns
  dplyr::select(-c(raw_richness,raw_T,corr_T)) %>%
  # mean and variance across space
  pivot_longer(c(mean_richness:shrub)) %>%
  # add count
  mutate(n_obs = 1) %>%
  group_by(pred_decade,elev_group,name) %>%
  summarise(mean=mean(value), sd=sd(value), n_obs = sum(n_obs))

sp.out <- sp.summ %>%
  rbind(sphx.summ)

write.csv(sp.out, paste0("analysis/q1_futureforestfloor/pft_lifeform_turnover_",rich.flag,".csv"),row.names=FALSE)

### winners and losers: pfts
sp.out %>%
  filter(name %in% c("light_cold","light_warm","light_notemp","shade_cold","shade_warm","shade_notemp")) %>%
  mutate(elev_group=factor(elev_group,levels=c("submontane","montane","subalpine"))) %>%
  mutate(name=factor(name,levels=c("light_cold","shade_cold","light_notemp","shade_notemp","light_warm","shade_warm"))) %>%
  ggplot(aes(x=as.character(pred_decade), y=mean, fill=name)) +
  facet_grid(~elev_group) +
  geom_bar(stat="identity") +
  ylab("Presence proportion") +
  scale_y_continuous(expand=c(0,0)) +
  scale_fill_manual(name="PFT",
                    values=c("#a6cee3","#1f78b4","#ffff99","#b15928","#fdbf6f","#ff7f00")) +
  theme_bw() +
  theme(
    axis.title.x = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle=45,hjust=1)
  )

# mean decline
sp.out %>%
  filter(name %in% c("light_cold","light_warm","light_notemp","shade_cold","shade_warm","shade_notemp")) %>%
  ungroup() %>%
  group_by(pred_decade,name) %>%
  summarise(pft_prop=mean(mean))
0.176-0.0913 # 0.08 decline in light_cold
(0.176-0.0913)/0.176 # 48% decline

### look at winners and losers in terms of life forms
sp.out %>%
  filter(name %in% c("graminoid","herb","shrub","fern")) %>%
  mutate(elev_group=factor(elev_group,levels=c("submontane","montane","subalpine"))) %>%
  mutate(name=factor(name,levels=c("graminoid","fern","herb","shrub"))) %>%
  ggplot(aes(x=as.character(pred_decade), y=mean, fill=name)) +
  facet_grid(~elev_group) +
  geom_bar(stat="identity") +
  ylab("Presence proportion") +
  scale_y_continuous(expand=c(0,0)) +
  scale_fill_manual(name="Life form",
                    values=c("#EEEEBB","#225522","#CCDDAA","#666633")) +
  theme_bw() +
  theme(
    axis.title.x = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle=45,hjust=1)
  )


### turnover: # retained, # new, # lost
sp.out %>%
  # filter(pred_decade==2050) %>%
  filter(name %in% c("turnover_new","turnover_lost")) %>%
  mutate(elev_group=factor(elev_group,levels=c("submontane","montane","subalpine"))) %>%
  mutate(name=factor(name,levels=c("turnover_new","turnover_lost"))) %>%
  ggplot(aes(x=as.character(pred_decade), y=mean, fill=name)) +
  facet_grid(~elev_group) +
  geom_bar(stat="identity") +
  ylab("Species turnover") +
  scale_y_continuous(expand=c(0,0), limits=c(0,100)) +
  scale_fill_manual(name="name",
                    values=c("#44AA99","#225555"),
                    labels=c("New species","Lost species")) +
  theme_bw() +
  theme(
    axis.title.x = element_blank(),
    legend.title=element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major.x=element_blank(),
    axis.text.x = element_text(angle=45,hjust=1)
  )

# mean turnover: values
sp.out %>%
  filter(name=="mean_turnover")

### species lost across entire landscape, by red list status
# read in counts for species
spsums.in <- read.csv("processed_data/iland_sdm_predictions/future_prediction_summaries/scenario_species_lost.csv") %>%
  filter(rich_flag==rich.flag, pred_decade!=2020)

# already have median species loss
pctile.out %>%
  filter(name=="lost")
# median -19 species, -7.7% relative to original species pool

summary(spsums.in)

# look at loss by red list status
redlist.lost <- spsums.in %>%
  left_join(redlist.in, by="species_name") %>%
  # add prevalence and retained
  mutate(prev = pres_plots/total_plots,
         retained=ifelse(pres_plots>0,1,0)) %>%
  # avg by rep
  group_by(pred_gcm,pred_wind,pred_rep,pred_decade,RLD) %>%
  summarise(lost=sum(lost),
            retained=sum(retained),
            prev=mean(prev)) %>%
  # replace nas for hx
  mutate(lost=replace_na(lost,0)) %>% 
  # then across reps
  ungroup() %>%
  group_by(pred_decade,RLD) %>%
  summarise(across(c(lost,retained,prev), ~quantile(., c(0,0.05,0.25,0.5,0.75,0.95,1))),q=c("min","q5","q25","median","q75","q95","max")) 
 

write.csv(redlist.lost, paste0("analysis/q1_futureforestfloor/turnover_redlist_",rich.flag,".csv"),row.names=FALSE)

# prop lost by red list classification
redlist.overall <- redlist.lost %>%
  filter(q=="median",pred_decade=="2050") %>%
  mutate(total=lost+retained) %>%
  ungroup() %>%
  dplyr::select(c(RLD,total))

redlist.lost %>%
  filter(q=="median") %>%
  left_join(redlist.overall, by="RLD") %>%
  mutate(prop_lost = lost/total) 
# mostly similar to overall loss rate, 8% not listed v. 4% for endangered, 33% for rare, 5% for warning

redlist.lost %>%
  filter(q=="median") %>%
  left_join(redlist.overall, by="RLD") %>%
  # across all 3 threatened categories
  mutate(threatened = ifelse(RLD %in% c("endangered","rare","warning"), "threatened",
                             ifelse(RLD %in% c("not_on_list"),"not_threatened",NA))) %>%
  group_by(pred_decade,threatened) %>%
  summarise(lost=sum(lost),total=sum(total)) %>%
  mutate(prop_lost = lost/total) # 7% v. 8%

redlist.lost %>%
  filter(RLD %in% c("warning","rare","endangered"),q=="median") %>%
  pivot_longer(c(lost:retained)) %>%
  mutate(RLD=factor(RLD,levels=c("warning","rare","endangered"))) %>%
  mutate(name=factor(name,levels=c("retained","lost"))) %>%
  ggplot(aes(x=RLD, y=value, fill=name)) +
  facet_grid(~pred_decade) +
  geom_bar(stat="identity") +
  ylab("Number of species") +
  scale_y_continuous(expand=c(0,0)) +
  scale_x_discrete(labels=c("Nearly threatened","Extremely rare","Threatened")) +
  scale_fill_manual(name="name",
                    values=c("#BBBBBB","#225555"),
                    labels=c("Retained species","Lost species")) +
  theme_bw() +
  theme(
    axis.title.x = element_blank(),
    legend.title=element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major.x=element_blank(),
    axis.text.x = element_text(angle=45,hjust=1)
  )


###
# Q2: relative importance of direct climate versus forest change
###

### read in samples
samp.read <- setNames(
  list.files(paste0("processed_data/iland_sdm_predictions/rand_sample_predictions/",rich.flag,"/"), full.names=TRUE), 
  list.files(paste0("processed_data/iland_sdm_predictions/rand_sample_predictions/",rich.flag,"/"))) %>%
  map_dfr(read.csv, .id = "run_id")

# add cover from subset for auc>0.7
if(rich.flag=="auc07") {
  cover.read <- setNames(
    list.files(paste0("processed_data/iland_sdm_predictions/rand_sample_predictions/subset/"), full.names=TRUE), 
    list.files(paste0("processed_data/iland_sdm_predictions/rand_sample_predictions/subset/"))) %>%
    map_dfr(read.csv, .id = "run_id") %>%
    dplyr::select(c(run_id,id_model,cover_pct))
  
  samp.read <- samp.read %>%
    left_join(cover.read, by=c("run_id","id_model"))
}

samp.names <- data.frame(run_id = unique(samp.read$run_id),run_id2 = unique(samp.read$run_id)) %>%
  separate(run_id2, into=c("pref1","pref2","gcm","rcp","gcm2","gcm3","pred_wind","pred_rep","pred_decade","suff"), sep="_") %>%
  mutate(forest_climate = paste(pref1,pref2,sep="_"),
         pred_gcm = paste(gcm,rcp,gcm2,gcm3,sep="_")) %>%
  dplyr::select(c(run_id,forest_climate,pred_decade,gcm,rcp,pred_wind,pred_rep)) 

samp.all <- samp.read %>%
  left_join(samp.names, by=c("run_id")) %>%
  dplyr::select(-run_id) %>%
  # prepare scenario levels
  mutate(forest = ifelse(forest_climate %in% c("historicalForest_futureClimate","historicalForest_historicalClimate"),"historicalForest",
                         ifelse(pred_wind==0,"futureForestLowDisturbance",
                                ifelse(pred_wind==15,"futureForestHighDisturbance",NA))),
         climate = ifelse(forest_climate %in% c("futureForest_futureClimate","historicalForest_futureClimate"),"futureClimate","historicalClimate")) %>%      
  # whether sorting by forest or climate first slightly affects NLME model fit
  arrange(forest,climate,id_model) 

### quick gut check: compare relative change among climate versus forest scenarios
head(samp.all)

samp.hx <- samp.all %>%
  filter(pred_decade==2020) %>%
  pivot_longer(c(corr_richness,corr_T,cover_pct)) %>%
  group_by(name) %>%
  summarise(value=mean(value)) 

samp.all %>%
  filter(pred_decade==2050) %>%
  pivot_longer(c(corr_richness,corr_T,cover_pct)) %>%
  group_by(forest,climate,pred_rep,pred_decade,name) %>%
  summarise(value=mean(value)) %>%
  ggplot(aes(y=value, fill=factor(forest))) +
  facet_wrap(~name, scales="free") +
  geom_hline(aes(yintercept=value),data=samp.hx, lty=2) +
  geom_boxplot() +
  theme_bw()

samp.all %>%
  filter(pred_decade==2050) %>%
  pivot_longer(c(corr_richness,corr_T,cover_pct)) %>%
  group_by(forest,climate,pred_rep,pred_decade,name) %>%
  summarise(value=mean(value)) %>%
  ggplot(aes(y=value, fill=factor(climate))) +
  facet_wrap(~name, scales="free") +
  geom_hline(aes(yintercept=value),data=samp.hx, lty=2) +
  geom_boxplot() +
  theme_bw()

samp.all %>%
  filter(pred_decade==2100) %>%
  pivot_longer(c(corr_richness,corr_T,cover_pct)) %>%
  group_by(forest,climate,pred_rep,pred_decade,name) %>%
  summarise(value=mean(value)) %>%
  ggplot(aes(y=value, fill=factor(forest))) +
  facet_wrap(~name, scales="free") +
  geom_hline(aes(yintercept=value),data=samp.hx, lty=2) +
  geom_boxplot() +
  theme_bw()

samp.all %>%
  filter(pred_decade==2100) %>%
  pivot_longer(c(corr_richness,corr_T,cover_pct)) %>%
  group_by(forest,climate,pred_rep,pred_decade,name) %>%
  summarise(value=mean(value)) %>%
  ggplot(aes(y=value, fill=factor(climate))) +
  facet_wrap(~name, scales="free") +
  geom_hline(aes(yintercept=value),data=samp.hx, lty=2) +
  geom_boxplot() +
  theme_bw()

# look at range of change 
clim.imp <- samp.all %>%
  filter(pred_decade!=2020) %>%
  pivot_longer(c(corr_richness,corr_T,cover_pct)) %>%
  group_by(climate,pred_decade,name) %>%
  summarise(value=mean(value)) %>%
  left_join(samp.hx, by=c("name")) %>%
  mutate(change=value.x-value.y) %>%
  ungroup() %>%
  group_by(pred_decade,name) %>%
  summarise(max=max(change),
            min=min(change),
            range=max-min) %>%
  mutate(group="climate")

for.imp <- samp.all %>%
  filter(pred_decade!=2020) %>%
  pivot_longer(c(corr_richness,corr_T,cover_pct)) %>%
  group_by(forest,pred_decade,name) %>%
  summarise(value=mean(value)) %>%
  left_join(samp.hx, by=c("name")) %>%
  mutate(change=value.x-value.y) %>%
  ungroup() %>%
  group_by(pred_decade,name) %>%
  summarise(max=max(change),
            min=min(change),
            range=max-min) %>%
  mutate(group="forest")
  
# compare the three
clim.imp %>%
  rbind(for.imp) %>%
  dplyr::select(-c(max,min)) %>%
  pivot_wider(names_from="group",values_from="range") %>%
  mutate(rel_clim = climate/(climate+forest),
         rel_forest = forest/(climate+forest))
# climate most impt except forest in 2050 for pct cover

### fit models: first, evaluate assumptions for each response variable
# this section is used to iteratively evaluate model assumptions based on 1 randomly selected replicate (6)
# response, decade, and transformation must be commented in/out below
# choose transformations to address linearity and normality
# choose model type and parameters to address unequal variance among groups

select_rep <- "6"

var.in <- "corr_richness"
# var.in <- "corr_T"
# var.in <- "cover_pct"

decade.in <- "2050"
# decade.in <- "2100"

# subset to variable, decade, rep of interest
simp.sub <- samp.all %>%
  pivot_longer(c(corr_richness,corr_T,cover_pct)) %>%
  mutate(pred_rep=ifelse(pred_wind=="NA"|pred_decade==2020,select_rep,pred_rep)) %>%
  filter(name==var.in,
         pred_decade %in% c("2020",decade.in),
         pred_rep==select_rep)

# number of observations in each group
simp.sub %>%
  group_by(forest) %>%
  tally()

simp.sub %>%
  group_by(climate) %>%
  tally()

# start with simple model, lm

# transformations
# transform <- 1 # no transformation, best T because transformation does not improve
transform <- (-1/2) # best richness
# transform <- (1/2) # best cover

# quick check of distribution and variance among groups
hist(simp.sub$value^transform)
simp.sub %>%
  group_by(scen) %>%
  summarise(var = var(value^transform))
leveneTest(lm(value~forest*climate,data=simp.sub))

# fit simple lm to assess linearity, normality
lm.fit <- lm(value^transform~forest*climate+id_model,data=simp.sub)
plot(lm.fit,1) 
car::qqPlot(resid(lm.fit))
qqline(rnorm(1000,mean(resid(lm.fit)),sd(resid(lm.fit))),col="red") # compare to simulated fit line
boxCox(lm(value~forest+climate+id_model,data=simp.sub))
summary(lm.fit)

# fit more complex model, address unequal variance
comp.sub <- simp.sub %>%
  mutate(value=value^transform)

# fit full model with main effects and their interactions
lme.fit <- lme(fixed=value~forest*climate,
               random=list(~1|id_model),
               data=comp.sub,
               method="ML",
               weights=varIdent(form=~1|forest*climate),
               control=lmeControl(msMaxIter = 1000, msMaxEval = 1000))

plot(lme.fit) 
car::qqPlot(resid(lme.fit))
qqline(rnorm(1000,mean(resid(lme.fit)),sd(resid(lme.fit))),col="red") # compare to simulated fit line

lme.fit
summary(lme.fit)
r.squaredGLMM(lme.fit)

# sqrt transform improves cover residual normality
# 1/sqrt transform (-1/2) improves richness residual normality
# transform does not consistently improve eivTemp residual normality, has long tails regardless

# assess relative importance of forest versus climate
lme.fit1 <- lme(fixed=value~forest,
                random=list(~1|id_model),
                data=comp.sub,
                method="ML",
                weights=varIdent(form=~1|forest),
                control=lmeControl(msMaxIter = 1000, msMaxEval = 1000))

lme.fit2 <- lme(fixed=value~climate,
                random=list(~1|id_model),
                data=comp.sub,
                method="ML",
                weights=varIdent(form=~1|climate),
                control=lmeControl(msMaxIter = 1000, msMaxEval = 1000))

summary(lme.fit1)
summary(lme.fit2)
r.squaredGLMM(lme.fit1)[1]
r.squaredGLMM(lme.fit2)[1]

### relative importance calculation
# total fixed
r.squaredGLMM(lme.fit)[1]

# sum of individ fixed
r.squaredGLMM(lme.fit1)[1] + r.squaredGLMM(lme.fit2)[1]

# denominator is max of these two
denom <- max(r.squaredGLMM(lme.fit)[1],r.squaredGLMM(lme.fit1)[1] + r.squaredGLMM(lme.fit2)[1])
# rel forest
r.squaredGLMM(lme.fit1)[1]/denom
# rel climate
r.squaredGLMM(lme.fit2)[1]/denom
# shared, might be negative
r.squaredGLMM(lme.fit)[1] - (r.squaredGLMM(lme.fit1)[1] + r.squaredGLMM(lme.fit2)[1])
1-r.squaredGLMM(lme.fit1)[1]/denom-r.squaredGLMM(lme.fit2)[1]/denom # shared only positive
# contribution of random effect
r.squaredGLMM(lme.fit)

### fit final models 
# fit separate models for each replicate, otherwise data are severely unbalanced, which affects results
# functionalize to extract variation partitioning outputs

fit_full <- function(df.in, decade.in, var.in,transform.in,scen) {
  # subset to variable, decade of interest
  comp.sub <- df.in %>%
    pivot_longer(c(corr_richness,corr_T,cover_pct)) %>%
    filter(name==var.in,
           pred_decade %in% c("2020",decade.in)) %>%
    mutate(value=value^transform.in)
  
  # fit full model with main effects and their interactions
  lme.fit <- lme(fixed=value~forest*climate,
                 random=list(~1|id_model),
                 data=comp.sub,
                 method="ML",
                 weights=varIdent(form=~1|forest*climate),
                 control=lmeControl(msMaxIter = 1000, msMaxEval = 1000))
  
  # fit models for relative importance of forest versus climate
  lme.fit1 <- lme(fixed=value~forest,
                  random=list(~1|id_model),
                  data=comp.sub,
                  method="ML",
                  weights=varIdent(form=~1|forest),
                  control=lmeControl(msMaxIter = 1000, msMaxEval = 1000))
  
  lme.fit2 <- lme(fixed=value~climate,
                  random=list(~1|id_model),
                  data=comp.sub,
                  method="ML",
                  weights=varIdent(form=~1|climate),
                  control=lmeControl(msMaxIter = 1000, msMaxEval = 1000))

  ### relative importance calculation
  data.out <- data.frame(scen=scen,
                         # marg and cond r2 for full model
                         r2m = r.squaredGLMM(lme.fit)[1],
                         r2c = r.squaredGLMM(lme.fit)[2],
                         # marg r2 for forest and climate
                         r2m_for= r.squaredGLMM(lme.fit1)[1],
                         r2m_clim= r.squaredGLMM(lme.fit2)[1]) %>%
    # calc relative contribution to r2m, shared contribution is difference with total
    mutate(denom = max(r.squaredGLMM(lme.fit)[1],r.squaredGLMM(lme.fit1)[1] + r.squaredGLMM(lme.fit2)[1]),
           for_ri = r.squaredGLMM(lme.fit1)[1]/denom,
           clim_ri = r.squaredGLMM(lme.fit2)[1]/denom,
           shared_ri = 1-for_ri-clim_ri,
           randEff_propR2 = (r2c-r2m)/r2c)
  
  return(data.out)
  
}

# apply model and transformations for each rcp and rep
data.reps <- data.frame()

for(i in 1:10) {
  print(i)
  select_rep <- i
  
  # subset to variable, decade, rep of interest
  df.sub <- samp.all %>%
    mutate(pred_rep=ifelse(pred_wind=="NA"|pred_decade==2020,select_rep,pred_rep)) %>%
    filter(pred_rep==select_rep)
  
  # apply models and transformations
  data.all <- 
    rbind(fit_full(df.sub,2050,"cover_pct",(1/2),"cover_2050"),
          fit_full(df.sub,2050,"corr_richness",(-1/2),"richness_2050"),
          fit_full(df.sub,2050,"corr_T",1,"temp_2050"),
          fit_full(df.sub,2100,"cover_pct",(1/2),"cover_2100"),
          fit_full(df.sub,2100,"corr_richness",(-1/2),"richness_2100"),
          fit_full(df.sub,2100,"corr_T",1,"temp_2100")) %>%
    mutate(rep=select_rep)
  
  data.reps <- rbind(data.reps,data.all)
  
}
data.reps

write.csv(data.reps, paste0("analysis/q2_clim_forest/replicate_model_varpart_",rich.flag,".csv"),row.names=FALSE)

### summary values: mean across all reps
# mean and se across reps
rep.summ <- data.reps %>%
  dplyr::select(c(scen,r2m,r2c,for_ri,clim_ri,shared_ri,randEff_propR2)) %>%
  group_by(scen) %>%
  summarise(across(everything(),list(mean=mean,se=std.error)))
write.csv(rep.summ,paste0("analysis/q2_clim_forest/model_varpart_summary_",rich.flag,".csv"),row.names=FALSE)

# figure
data.reps %>%
  separate(scen, into=c("response","decade"),sep="_") %>%
  pivot_longer(c(for_ri,clim_ri,shared_ri)) %>%
  group_by(response,decade,name) %>%
  summarise(value=round(mean(value),3)) %>%
  mutate(response = factor(response, levels=c("richness","temp","cover")),
         name = factor(name,levels=c("for_ri","shared_ri","clim_ri"))) %>%
  ggplot(aes(x=response, y=value, fill=name)) +
  facet_grid(~decade) +
  geom_bar(stat="identity") +
  ylab("Relative importance") +
  scale_y_continuous(expand=c(0,0), limits=c(0,1.002)) +
  scale_x_discrete(labels=c("Richness","Ellen T","Cover")) +
  scale_fill_manual(name="name",
                    values=c("#EECC66","#969696","#994455"),
                    labels=c("Forest change","Shared","Climate change")) +
  theme_bw() +
  theme(
    axis.title.x = element_blank(),
    legend.title=element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major.x=element_blank(),
    axis.text.x = element_text(angle=45,hjust=1)
  )


# when is forest more impt than climate and vice versa
data.reps %>%
  mutate(forest = ifelse(for_ri>clim_ri,1,0),
         climate = ifelse(clim_ri>for_ri,1,0)) %>%
  group_by(scen) %>%
  summarise(forest = sum(forest),climate=sum(climate))
# almost always climate  
