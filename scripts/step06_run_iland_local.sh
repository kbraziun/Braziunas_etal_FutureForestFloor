#!/bin/bash
# bash script to run iland for a list of scenarios
# script should be run from the main directory, which is the
# directory containing the readme
# script takes 1 argument: number of reps
# usage: bash scripts/step06_run_iland_local.sh [reps]

#####
# catch errors
#####

set -e # terminate if non-zero exit
set -u # terminate if variable unset
set -o pipefail # terminate if pipe command fails

#####
# arguments
#####

# input arguments
reps=$1
path="/iland_exe/executable_qt512.6/" # set this path to location of iland exe
simulation_years="80"

# set filename of csv that will be used to customize climate and environment in each project file
# this script is hard-coded to specific columns in a specific order
csv_name="iland/master_file.csv"
# xml="iland/cc_wind.xml"

#####
# loop through entire csv and run iland for each site
#####

sed 1d $csv_name | while IFS=, read -r gcm rcp wind xml
do 
    for rep in $(seq 1 $reps)
    do
        echo "running gcm $gcm rcp $rcp with wind speed $wind rep $rep xml $xml"

        # run iland with timevents
        ${path}ilandc.exe iland/${xml} $simulation_years \
        system.database.climate=${gcm}.sqlite \
        model.world.timeEventsFile=timeevents_${wind}_0_${rcp}_${rep}.txt
        
        mkdir iland/output/${gcm}_${wind}_${rep}

        mv iland/output/future.sqlite iland/output/${gcm}_${wind}_${rep}
        mv iland/output/lif_30.asc iland/output/${gcm}_${wind}_${rep}
        mv iland/output/lif_80.asc iland/output/${gcm}_${wind}_${rep}

        # run R script to generate outputs
        # need to set PATH and environment variables appropriately on your computer
        Rscript.exe iland/scripts/prep_iland_outputs.R ${gcm}_${wind}_${rep} "future.sqlite"

        # remove iland output, takes up space
        rm -r iland/output/${gcm}_${wind}_${rep}
  
    done
done
