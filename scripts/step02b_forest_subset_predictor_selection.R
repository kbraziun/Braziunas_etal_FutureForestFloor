#####
# 
## Step 2b: Select most important forest predictors with RF model
#
#####

# if running on server
# write("TMP = 'temp''", file = file.path('~/.Renviron'))

# here, taking forest structure predictors and identifying top 5 subset for building final model
# using random forest

# load libraries
library(tidyverse)
library(randomForest)

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
#   [1] randomForest_4.7-1 forcats_0.5.1      stringr_1.4.1      dplyr_1.0.8        purrr_0.3.4        readr_2.1.2       
# [7] tidyr_1.2.0        tibble_3.1.6       ggplot2_3.3.5      tidyverse_1.3.1  

####
# 1. load data
####

# response
pres.in <- read.csv("processed_data/understory_model/understory_subset_presence_cover.csv") %>%
  dplyr::select(-cover)

# predictors: forest
pred.in <- read.csv("processed_data/understory_model/predictors_fullset.csv") %>%
  # only forest predictor set
  dplyr::select(c(plot_id,TSFdec,tph,BA,mean_dbh,sd_height,dom_height,tree_richness,prop_evergreen,prop_fasy))
                  
# pft lookup
pft.in <- read.csv("processed_data/species_lookup/final_pft_categories_ms.csv") %>%
  separate(pft, into=c("light","temp"), sep="-")

####
# 2. prep data
####

# standardize all predictors
pred.std <- pred.in %>%
  mutate(across(c(TSFdec:prop_fasy), ~scale(.))) %>%
  group_by(plot_id) %>%
  mutate(across(everything(), ~as.numeric(.))) %>%
  as.data.frame() %>%
  # set NAs prop_evergreen and prop_fasy to 0 (mean) as equalizing value
  mutate(across(c(prop_evergreen,prop_fasy), ~replace_na(.,0)))

###
# 3. fit SDMs for each species, extract variable importance
###

# use this to get random numbers
# set.seed ensures same result if code rerun, choose initial value based on random number between 1 and 1000
sample(1:1000,1)

# number of predictors, used in model fitting
npred <- ncol(pred.std) - 1

### fit final sdms
# fit model to full dataset for each species
sp.list <- unique(pres.in$species_name)

# variable importance and partial response curves
vimp.out <- data.frame()
rcurves.out <- data.frame()

for(j in sp.list) {
  print(j) 
  sp.name <- j
  
  # join response and predictors
  rf.data <- pres.in %>%
    filter(species_name == sp.name) %>%
    # join with predictors
    left_join(pred.std, by="plot_id")
  
  # count of p and a
  pres.count <- rf.data %>%
    group_by(pres) %>%
    tally() %>%
    mutate(pres = ifelse(pres==1,"present","absent")) %>%
    pivot_wider(names_from="pres",values_from="n")
  
  # add weights based on relative proportion presence (weight 1) and absence (rel weight)
  pres.wt <- rf.data %>%
    summarise(wt_1=1,
              wt_0=pres.count$present/pres.count$absent) %>%
    # join with rf.data
    cbind(rf.data) %>%
    # set 0s to 1
    mutate(pres_wt = ifelse(pres==1,wt_1,wt_0)) %>%
    dplyr::select(-c(wt_1,wt_0))
  
  # # double check that overall weights are equal
  # pres.wt %>%
  #   group_by(pres) %>%
  #   summarise(pres_wt=sum(pres_wt))
  
  # subset to data for model fitting
  pres.mod <- pres.wt %>%
    dplyr::select(-c(species_name,pres_wt,plot_id))
  
  # fit rf
  set.seed(269)
  rf.mod <- randomForest(pres~.,data=pres.mod, 
                         weights=pres.wt$pres_wt,
                         ntree=1000,
                         # default nodesize for regression=5
                         nodesize=5,
                         # default mtry for regression is floor(npred/3)
                         mtry=floor(npred/3),
                         sampsize=0.7*nrow(pres.mod),
                         importance=TRUE)
  
  # get variable importance
  vimp.mod <- randomForest::importance(rf.mod) %>%
    as.data.frame() %>%
    mutate(species_name = j,
           variables=row.names(.))
  
  vimp.out <- rbind(vimp.out, vimp.mod)
  
  # get partial plot values
  rcurves.mod <- data.frame()
  for (i in seq_along(vimp.mod$variables)) {
    rcurves.var <- as.data.frame(partialPlot(rf.mod, as.data.frame(pres.mod[,-1]), vimp.mod$variables[i],plot=FALSE)) %>%
      mutate(variable = vimp.mod$variables[i])

    rcurves.mod <- rbind(rcurves.mod,rcurves.var)
  }

  rcurves.mod$species_name <- j
  rcurves.out <- rbind(rcurves.out,rcurves.mod)
  
}

###
# 4. evaluate SDMs, response curves and variable importance
###

load("processed_data/understory_model/evaluation/forest_predictor_selection/vimp_forest_predictors.RData")

head(rcurves.out)

### does model align with ecological expectations
# response curves, look based on temperature and light response
rcurves.summ <- rcurves.out %>%
  left_join(pft.in, by="species_name")

# check response curves based on light response
pft.ck <- "light"
pft.ck <- "shade"

# # fast, simple mean
# rcurves.summ %>%
#   filter(licht_cat==pft.ck) %>%
#   group_by(licht_cat,variable,x) %>%
#   summarise(y=mean(y)) %>%
#   ggplot(aes(x=x,y=y)) +
#   facet_wrap(~variable, scales="free") +
#   geom_line() +
#   theme_bw()

# loess fit to all points, slower to plot
rcurves.summ %>%
  filter(light==pft.ck) %>%
  ggplot(aes(x=x,y=y)) +
  # facet_wrap(~variable, scales="free") +
  facet_wrap(~variable, scales="free_x") +
  geom_smooth(method="loess") +
  theme_bw()
# for light-pref: higher values assoc with lower BA, lower dom height, higher mean_dbh, higher prop_evergreen, lower prop_fasy, intermediate and higher sd_height with some noise, lower tph, no clear trend richness, higher TSF; TSF dominance clear when set to "free_x"
# for shade-tol: higher values assoc with lower BA, higher dom height, lower mean_dbh, lower prop_evergreen, higher prop_fasy, higher sd_height, lower tph, lower or no effect richness, lower TSF; TSF, dom_height, and prop_fasy most dominant predictors

head(vimp.out)

# variable importance
vimp.summ  <- vimp.out %>%
  # set neg values of IncMSE to 0 to avoid depressing overall mean
  rename(IncMSE = "%IncMSE") %>%
  mutate(IncMSEPos = ifelse(IncMSE>0, IncMSE, 0)) %>%
  group_by(variables) %>%
  summarise(IncMSEPos=mean(IncMSEPos))

vimp.summ %>%
  mutate(variables = factor(variables)) %>%
  ggplot(aes(x=reorder(variables,-desc(IncMSEPos)),y=IncMSEPos)) +
  geom_bar(stat='identity') +
  coord_flip() +
  theme_bw() +theme(legend.position = "none")

save(rcurves.out,vimp.out,file=paste0("processed_data/understory_model/evaluation/forest_predictor_selection/","vimp_forest_predictors.RData"))

# most important vars: TSFdec, BA, prop_fasy.
# tph and dom_ht similar; chose dom_ht a priori, provides distinct info and differing relationships based on response curves
# next best composition is prop_evergreen
# larger set: TSF, BA, dom_height, prop_fasy, prop_evergreen
# smaller set: TSF, dom_height, prop_fasy