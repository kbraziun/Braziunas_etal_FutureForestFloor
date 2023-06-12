#####
# 
## Step 3a: Fit SDMs to presence-absence data
#
#####

# # if running on server
# write("TMP = 'temp''", file = file.path('~/.Renviron'))

###
# load libraries
###

# libraries
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

# species subset: at least 5 observations. used for final model fits
pres.in <- read.csv("processed_data/understory_model/understory_subset_presence_cover.csv") %>%
  dplyr::select(-cover)

# predictors: subset with 3 climate and 3 forest
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
# 3. fit evaluation SDMs for each species
###

# can fit weighted or unweighted, assign predictor set, track evaluation
weighted <- "weighted"
pred.set <- pred.std
npred <- ncol(pred.set) -1
eval <- "run13"

# use this to get random numbers
# set.seed ensures same result if code rerun, choose initial value based on random number between 1 and 1000
sample(1:1000,20)

### fit evaluation sdms
# fit 20 reps, each with subset of 70% of plots to minimize overfitting
# use holdout 30% for evaluation
# final prediction will be average of model test data

seeds <- c(785,  65, 693, 842, 230, 823, 540, 897, 975, 832,389, 439, 260, 715, 357,46, 911, 922, 129, 352)
sp.list <- unique(pres.in$species_name)

for(j in sp.list) {
  print(j)
  sp.name <- j
  
  rf.data <- pres.in %>%
    filter(species_name == sp.name) %>%
    # join with predictors
    left_join(pred.set, by="plot_id") 
  
  # subset presences and absences separately
  pres.sub <- rf.data %>%
    filter(pres==1)
  
  abs.sub <- rf.data %>%
    filter(pres==0)
  
  ntot <- nrow(rf.data) # total # of plots
  
  # # check to see if all plots included in at least 1 test set
  # test.plots <- data.frame()
  
  # run 20 reps
  for(i in c(1:20)) {
    print(i)
    
    # subset to 70% each of absences and presences
    set.seed(seeds[i])
    pres.samp <- sample(pres.sub$plot_id, round(0.7 * nrow(pres.sub)), replace=FALSE)
    set.seed(seeds[i])
    abs.samp <- sample(abs.sub$plot_id, round(0.7 * nrow(abs.sub)), replace=FALSE)

    # training data: 70% subset
    rf.train <- rf.data %>%
      filter(plot_id %in% c(pres.samp,abs.samp)) %>%
      mutate(run_nbr = i)
    
    # testing data: the rest
    rf.test <- rf.data %>%
      filter(!plot_id %in% c(pres.samp,abs.samp)) %>%
      mutate(run_nbr = i)
    
    # separate process for weighted or unweighted
    if(weighted=="weighted") {
      
      # add weights based on relative proportion presence (weight 1) and absence (rel weight)
      pres.wt <- data.frame(wt_1=1,
                  wt_0=length(pres.samp)/length(abs.samp)) %>%
        # join with rf.train
        cbind(rf.train) %>%
        # set 0s to 1
        mutate(pres_wt = ifelse(pres==1,wt_1,wt_0)) %>%
        dplyr::select(-c(wt_1,wt_0))
      
      # # double check that overall weights are equal
      # pres.wt %>%
      #   group_by(pres) %>%
      #   summarise(pres_wt=sum(pres_wt))

      # subset to data for model fitting
      pres.mod <- pres.wt %>%
        dplyr::select(-c(species_name,pres_wt,plot_id,run_nbr))
      
      set.seed(seeds[i])
      rf.mod <- randomForest(pres~.,data=pres.mod,
                             weights=pres.wt$pres_wt,
                             ntree=1000,
                             # default nodesize for regression=5
                             nodesize=5,
                             # default mtry for regression=floor(npred/3)
                             mtry=floor(npred/3),
                             sampsize=nrow(pres.mod))
        
    } else if(weighted=="unweighted") {

      # subset to data for model fitting
      pres.mod <- rf.train %>%
        dplyr::select(-c(species_name,plot_id,run_nbr))

      set.seed(seeds[i])
      rf.mod <- randomForest(pres~.,data=pres.mod,
                             ntree=1000,
                             # default nodesize for regression=5
                             nodesize=5,
                             # default mtry for regression=floor(npred/3)
                             mtry=floor(npred/3),
                             sampsize=nrow(pres.mod))
      
    }

    save(rf.mod,rf.train,file=paste0("processed_data/understory_model/evaluation/sdm_evaluation_fits/",j,"_run",i,".RData"))
    
    # test.plots <- rbind(test.plots,rf.test) # use this to check to see if all plots end up in training data at least once, n=10 not quite enough, trying n=20 reps which gets very close
    
    rm(pres.samp,abs.samp,rf.train,rf.test,pres.mod,rf.mod)
    
  }
  
  rm(rf.data,pres.sub,abs.sub)
  gc()
  
}

###
# 4. load evaluation SDMs for each species and save only predictions
###

train.out <- data.frame()

for(k in c(list.files("processed_data/understory_model/evaluation/sdm_evaluation_fits_run13/", full.names = TRUE))) {
  
  print(k)
  
  load(k)

  sp_name <- as.data.frame(k) %>%
    separate(k, into=c("dir1","dir2","dir3","dir4","name"), sep="/") %>%
    separate(name, into=c("species_name","suf"), sep="_run") %>%
    separate(suf, into=c("run","suf"), sep="[.]")
  
  # add flag for training v. test
  sub.lookup <- rf.train %>%
    dplyr::select(c(plot_id,species_name,run_nbr)) %>%
    mutate(train="train")
  
  # join with full dataset
  rf.pres <- pres.in %>%
    filter(species_name %in% c(sub.lookup$species_name)) %>%
    left_join(pred.set, by="plot_id") %>%
    left_join(sub.lookup, by=c("plot_id","species_name")) %>%
    mutate(run_nbr = mean(run_nbr, na.rm=TRUE)) %>%
    # if not in training dataset, label as test
    mutate(train = ifelse(is.na(train),"test",train))
  
  # add prediction
  rf.pres$rf_pred <- predict(rf.mod, newdata=rf.pres)
  
  # reduce output
  pred.train <- rf.pres %>%
    dplyr::select(c(plot_id,pres,species_name,run_nbr,train,rf_pred))
  
  # add to dataframe
  train.out <- rbind(train.out, pred.train)
  
}

write.csv(train.out, paste0("processed_data/understory_model/evaluation/sdm_pres_evaluation_",eval,".csv"),row.names=FALSE)


###
# 5. fit final SDMs for each species
###

# use this to get random numbers
# set.seed ensures same result if code rerun, choose initial value based on random number between 1 and 1000
sample(1:1000,1)

### fit final sdms
# fit model to full dataset for each species to use for prediction
# weighted approach
# predictor set with 3 climate, 3 forest, 3 site

sp.list <- unique(pres.in$species_name)
pred.set <- pred.std3
npred <- ncol(pred.set) -1

for(j in sp.list) {
  print(j)
  sp.name <- j 
  
  rf.data <- pres.in %>%
    filter(species_name == sp.name) %>%
    # join with predictors
    left_join(pred.set, by="plot_id") 
  
  # count of p and a
  pres.count <- rf.data %>%
    group_by(pres) %>%
    tally() %>%
    mutate(pres = ifelse(pres==1,"present","absent")) %>%
    pivot_wider(names_from="pres",values_from="n")
  
  # add weights based on relative proportion presence (weight 1) and absence (rel weight)
  pres.wt <- data.frame(wt_1=1,
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
  set.seed(540)
  rf.mod <- randomForest(pres~.,data=pres.mod,
                       weights=pres.wt$pres_wt,
                       ntree=1000,
                       # default nodesize for regression=5
                       nodesize=5,
                       # default mtry for regression=floor(npred/3)
                       mtry=floor(npred/3),
                       sampsize=nrow(pres.mod),
                       importance=TRUE)
  print(rf.mod)
  
  save(rf.mod,file=paste0("processed_data/understory_model/final_fits/sdm_prediction_fits/",j,".RData"))
  
}

###
# 6. get variable importance and rcurves from final fits
###

train.out <- data.frame()
rcurves.out <- data.frame()
vimp.out <- data.frame()

for(k in c(list.files("processed_data/understory_model/final_fits/sdm_prediction_fits/", full.names = TRUE))) {
  
  print(k)
  
  load(k)
  
  sp_name <- as.data.frame(k) %>%
    separate(k, into=c("dir1","dir2","dir3","dir4","name"), sep="/") %>%
    separate(name, into=c("species_name","suf"), sep=".RD")
  
  # predict based on full model
  rf.pres <- pres.in %>%
    filter(species_name %in% c(sp_name$species_name)) %>%
    left_join(pred.set, by="plot_id") 
  
  # add prediction
  rf.pres$rf_pred <- predict(rf.mod, newdata=rf.pres)
  
  # reduce output
  pred.train <- rf.pres %>%
    dplyr::select(c(plot_id,pres,species_name,rf_pred))
  
  train.out <- rbind(train.out, pred.train)
  
  # get variable importance
  vimp.mod <- importance(rf.mod) %>%
    as.data.frame() %>%
    mutate(species_name = sp_name$species_name,
           variables=row.names(.))

  vimp.out <- rbind(vimp.out, vimp.mod)

  # get partial plot values
  rcurves.mod <- data.frame()
  for (i in seq_along(vimp.mod$variables)) {
    
    rcurves.var <- as.data.frame(partialPlot(rf.mod, as.data.frame(rf.pres[,-1]), vimp.mod$variables[i],plot=FALSE)) %>%
      mutate(variable = vimp.mod$variables[i])

    rcurves.mod <- rbind(rcurves.mod,rcurves.var)
  }

  rcurves.mod$species_name <- sp_name$species_name
  rcurves.out <- rbind(rcurves.out,rcurves.mod)
}


head(train.out)
head(rcurves.out)
head(vimp.out)

# save these
write.csv(train.out, paste0("processed_data/understory_model/evaluation/sdm_pres_finalFit_",npred,"pred.csv"),row.names=FALSE)
write.csv(rcurves.out, paste0("processed_data/understory_model/evaluation/sdm_pres_finalFit_rcurves_",npred,"pred.csv"),row.names=FALSE)
write.csv(vimp.out, paste0("processed_data/understory_model/evaluation/sdm_pres_finalFit_vimp_",npred,"pred.csv"),row.names=FALSE)

