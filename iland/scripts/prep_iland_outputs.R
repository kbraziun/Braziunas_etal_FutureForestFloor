#####
# 
## Prepare iland outputs for SDM model predictions
#
#####

### load libraries
library(tidyverse)
library(sp)
library(terra)
library(sf)
library(RSQLite)

### to update
# path for future runs as argument in
# lif locations
# sqlite output location

# running from command line
args <- commandArgs(TRUE)
path <- as.character(args[1])
sqlite <- as.character(args[2])

# path <- "future_basewind"
# sqlite <- "future_basewind_ECEARTHrcp85.sqlite"

###
# 1. import lif data 
###

# standgrid for masking
grid.in <- rast("iland/gis/standgrid.asc")
crs(grid.in) <- "epsg:31468"

# remove seedbelt and classify rest as 1
grid.noseed <- classify(grid.in, matrix(c(-Inf,0,NA,1,Inf,1),ncol=3,byrow=TRUE))

### lif at 2m height
lif <- c(rast(paste0("iland/output/",path,"/lif_30.asc")),rast(paste0("iland/output/",path,"/lif_80.asc")))
crs(lif) <- crs(grid.in)

# mask
lif.mask <- crop(lif,grid.noseed)*grid.noseed

# plot(lif.mask)

###
# 2. import iland tree outputs
###

conn=DBI::dbConnect(RSQLite::SQLite(), dbname = paste0("iland/output/",path,"/",sqlite)) # connect to the db
dbListTables(conn)

# tree output, initial filtering
beetle.in <-  tbl(conn, "barkbeetle") %>%
  collect()

wind.in <-  tbl(conn, "wind") %>%
  collect()

land.in <-  tbl(conn, "landscape") %>%
  collect()

# tree output, initial filtering
trees.in <-  tbl(conn, "tree") %>%
  dplyr::select(c(year,ru,rid,species,id,x,y,height,basalArea)) %>%
  mutate(x_crs = x + 4559008, y_crs = y+5260630) %>%
  collect()

dbDisconnect(conn) # close the connection

# summary(trees.in)
# summary(land.in)
# summary(beetle.in)
# summary(wind.in)

### write out disturbance and landscape outputs
write.csv(land.in,paste0("iland/output/combined_outputs/",path,"_landscape.csv"), row.names=FALSE)
write.csv(beetle.in,paste0("iland/output/combined_outputs/",path,"_beetle.csv"), row.names=FALSE)
write.csv(wind.in,paste0("iland/output/combined_outputs/",path,"_wind.csv"), row.names=FALSE)

###
# 3. create ids for 10m cells, match up with tree data 
###

# load copy of raster with unique ids, mask to same grid
id.rast <- rast("processed_data/iland_sdm_predictions/stand_id_raster.tif")

trees.pts <- vect(trees.in, geom=c("x_crs","y_crs"), crs=crs(id.rast))

# extract raster ids
trees.in$cell_id <- extract(id.rast, trees.pts)[,2]

rm(trees.pts)
gc()

# # how many trees per 10m cell: 1-36, mean 4.4
# trees.in %>%
#   group_by(cell_id) %>%
#   tally() %>% summary()

###
# 4. summarize predictor variables of interest at 10m resolution
###

### light

lif.pred <- values(c(lif.mask,id.rast)) %>%
  as.data.frame() %>%
  rename(cell_id=lyr.1) %>%
  filter(!is.na(lif_30)) %>%
  pivot_longer(cols=c(lif_30:lif_80),names_to="year",names_prefix="lif_",values_to="lif") %>%
  mutate(year=as.integer(year))

### forest structure
# unique(trees.in$species)
# evergreen: "piab", "abal", "pisy","psme","pini","pice","pimu"
# not evergreen: "lade","fasy","quro","acps","frex","cabe","bepe","alin","qupe","algl","casa","acca","acpl","qupu","soau","soar","coav","alvi","potr","poni","tico","tipl","ulgl","saca","rops"

# head(trees.in)

# structure
trees.struct <- trees.in %>%
  group_by(year,cell_id) %>%
  # basal area
  summarise(BA = sum(basalArea),
            # dominant height
            dom_height = quantile(height, probs=c(0.95)))

# composition
trees.comp <- trees.in %>%
  # assign identifiers to evergreen
  # no trailing NAs, good
  mutate(evergreen = ifelse(species %in% c("piab", "abal", "pisy","psme","pini","pice","pimu"),1,
                            ifelse(species %in% c("lade","fasy","quro","acps","frex","cabe","bepe","alin","qupe","algl","casa","acca","acpl","qupu","soau","soar","coav","alvi","potr","poni","tico","tipl","ulgl","saca","rops"),0,NA))) %>%
  group_by(year,cell_id,species,evergreen) %>%
  # BA by species
  summarise(BA_spec = sum(basalArea)) %>%
  # add counter for species
  mutate(n_species=1)

# fasy BA
trees.fasy <- trees.comp %>%
  filter(species=="fasy") %>%
  ungroup() %>%
  # drop unneeded columns
  dplyr::select(-c(species,evergreen,n_species))

# evergreen BA
trees.ever <- trees.comp %>%
  ungroup() %>%
  # only evergreen
  filter(evergreen==1) %>%
  group_by(year,cell_id,evergreen) %>%
  summarise(BA_spec=sum(BA_spec))

# bringing it all together
trees.agg <- trees.struct %>%
  left_join(trees.ever, by=c("year","cell_id")) %>%
  # convert to proportion
  mutate(prop_evergreen = BA_spec/BA) %>%
  # drop excess columns
  dplyr::select(-c(evergreen,BA_spec)) %>%
  # add fasy
  left_join(trees.fasy, by=c("year","cell_id")) %>%
  mutate(prop_fasy = BA_spec/BA) %>%
  dplyr::select(-c(BA_spec)) %>%
  # replace nas with 0s, because trees are present in all of these plots
  mutate(across(c("prop_evergreen","prop_fasy"),~replace_na(.,0))) %>%
  # convert BA from m2/10x10m to m2/ha
  mutate(BA=BA*100)

### climate + site
# do separately, 1 site for all, climate by GCM and year

### all together
pred.for <- trees.agg %>%
  left_join(lif.pred, by=c("year","cell_id"))

write.csv(pred.for,paste0("iland/output/combined_outputs/",path,"_forest_predictors.csv"),row.names=FALSE)
