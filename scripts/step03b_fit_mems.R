#####
# 
## Step 3b: Fit random forest to predict percent cover of PFTs
#
#####

###
# load libraries
###

# Load relevant packages
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

###
# 1. import data
###

### responses
# total cover
cover.in <- read.csv("processed_data/understory_model/understory_presence_cover.csv") %>%
  group_by(plot_id) %>%
  summarise(cover=sum(cover))

# richness of subset
rich.in <- read.csv("processed_data/understory_model/understory_subset_presence_cover.csv") %>%
  group_by(plot_id) %>%
  summarise(richness=sum(pres))

### predictors: subset with 3 climate and 3 forest
pred.in <-  read.csv("processed_data/understory_model/predictors_subset.csv") %>%
  # only final predictor set, 3 climate, 3 forest, 3 site
  dplyr::select(c(plot_id,mean_temp,summer_prec,rad,TSFdec,BA,prop_fasy,sand,whc,soil_fert))

####
# 2. prep data
####

# standardize all predictors
pred.std <- pred.in %>%
  mutate(across(c(mean_temp:soil_fert), ~scale(.))) %>%
  group_by(plot_id) %>%
  mutate(across(everything(), ~as.numeric(.))) %>%
  as.data.frame() %>%
  # set NAs prop_fasy to 0 (mean) as equalizing value
  mutate(across(c(prop_fasy), ~replace_na(.,0)))

###
# 3. fit random forest regressions for evaluation
###

### fit macroecological models (MEMs) for total cover and individual species richness (for subset with individual SDMs)

# set predictors, cover or richness
pred.set <- pred.std
resp.mod <- cover.in
mod.name <- "cover"
# resp.mod <- rich.in
# mod.name <- "richness"
npred <- ncol(pred.set) -1

rf.data <- resp.mod %>%
  left_join(pred.set, by="plot_id")

### use same test partitioning approach as with sdms to allow explicit testing of training data
sample(1:1000,20)

seeds <- c(886, 818, 63, 705, 558, 496, 323, 1000, 33, 551, 742, 763, 244, 229, 422, 955, 66, 710, 171, 492)
test.out <- data.frame()

for(i in c(1:20)) {
  print(i)
  
  # subset to 70% of plots
  set.seed(seeds[i])
  rf.sub <- sample(unique(rf.data$plot_id), 45, replace=FALSE)
  
  # subset to data for model fitting
  rf.fit <- rf.data %>%
    filter(!plot_id %in% c(rf.sub)) %>%
    dplyr::select(-plot_id)
  
  # fit rf model
  if(mod.name=="cover") {
    set.seed(seeds[i])
    rf.mod <- randomForest(cover~.,data=rf.fit, 
                           ntree=1000,
                           # default nodesize for regression=5
                           nodesize=5,
                           # default mtry for regression=floor(npred/3)
                           mtry=floor(npred/3),
                           sampsize=nrow(rf.fit))
    
    # predict for test dataset
    test.mod <- rf.data %>%
      filter(plot_id %in% c(rf.sub))%>%
      mutate(run = i)
    
    test.mod$cover_pred <- predict(rf.mod, newdata=test.mod) 
  } else if(mod.name=="richness") {
    set.seed(seeds[i])
    rf.mod <- randomForest(richness~.,data=rf.fit, 
                           ntree=1000,
                           # default nodesize for regression=5
                           nodesize=5,
                           # default mtry for regression=floor(npred/3)
                           mtry=floor(npred/3),
                           sampsize=nrow(rf.fit))
    
    # predict for test dataset
    test.mod <- rf.data %>%
      filter(plot_id %in% c(rf.sub))%>%
      mutate(run = i)
    
    test.mod$richness_pred <- predict(rf.mod, newdata=test.mod) 
  }
  
  test.out <- rbind(test.out,test.mod)
  
}

if(mod.name=="cover") {
  test.eval <- test.out %>%
    dplyr::select(c(plot_id,cover,run,cover_pred))
} else if(mod.name=="richness") {
  test.eval <- test.out %>%
    dplyr::select(c(plot_id,richness,run,richness_pred))
}

# write out
write.csv(test.eval, paste0("processed_data/understory_model/evaluation/mem_",mod.name,"_evaluation_",npred,"pred.csv"),row.names=FALSE)

###
# 4. final MEM fit: cover only
###

### fit macroecological models (MEMs) for total cover

# set predictors, cover or richness
pred.set <- pred.std
resp.mod <- cover.in
mod.name <- "cover"
npred <- ncol(pred.set) -1

rf.data <- resp.mod %>%
  left_join(pred.set, by="plot_id")

# fit full model for use in prediction, retain previously used sample size
rf.fit <- rf.data %>%
  dplyr::select(-plot_id)

set.seed(693)
final.mod <- randomForest(cover~.,data=rf.fit, 
                          ntree=1000,
                          # default nodesize for regression=5
                          nodesize=5,
                          # default mtry for regression=floor(npred/3)
                          mtry=floor(npred/3),
                          sampsize=nrow(rf.fit),
                          importance=TRUE)

final.mod # OOB var explained aligns with tests

# save model for prediction
save(final.mod, file=paste0("processed_data/understory_model/final_fits/mem_prediction_cover/cover_final.RData"))

### evaluate partial plots and variable importance for final model
varImpPlot(final.mod)

# extract importance
imp <- importance(final.mod)

imp.out <- imp %>%
  as.data.frame() %>%
  mutate(variables = (row.names(imp)))

write.csv(imp.out, paste0("processed_data/understory_model/evaluation/mem_cover_finalFit_vimp_",npred,"pred.csv"),row.names=FALSE)

# partial dependence
impvar <- rownames(imp)[order(imp[, 1], decreasing=TRUE)]

rcurves.mod <- data.frame()

for (i in seq_along(impvar)) {
  rcurves.var <- as.data.frame(partialPlot(final.mod, as.data.frame(rf.data[,-1]), impvar[i],plot=FALSE)) %>%
    mutate(variable = impvar[i])
  
  rcurves.mod <- rbind(rcurves.mod,rcurves.var)
}

write.csv(rcurves.mod, paste0("processed_data/understory_model/evaluation/mem_cover_finalFit_rcurves_",npred,"pred.csv"),row.names=FALSE)
