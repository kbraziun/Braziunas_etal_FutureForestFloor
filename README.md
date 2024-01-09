# readme for Braziunas et al. Future Forest Floor

## purpose

This readme gives an overview of directory structure, files, and steps for recreating
outputs and analyses associated with Braziunas et al. Future Forest Floor.

**Manuscript citation:** Braziunas, K.H., L. Geres, T. Richter, F. Glasmann, C. Senf, D. Thom, S. Seibold, and R. Seidl. 2024. Projected climate and canopy change lead to thermophilization and homogenization of forest floor vegetation in a hotspot of plant species richness. Global Change Biology 30(1):e17121.

**Dataset citation:** Braziunas, K.H., L. Geres, T. Richter, F. Glasmann, C. Senf, D. Thom, S. Seibold, and R. Seidl. 2023. Projected climate and canopy change lead to thermophilization and homogenization of forest floor vegetation in a hotspot of plant species richness, Berchtesgaden National Park, Bavaria, Germany ver 1. Environmental Data Initiative. https://doi.org/10.6073/pasta/06954751c87ea7695dd043844e2745f9

## platforms

- Operating systems and software used for development and implementation
  - OS: Windows 10
  - R version: 4.1.3

## directory overview

Directory structure and files:

- `analysis/`: Results from data analysis.
- `iland/`: Inputs and project files used to rerun forest simulations. Will be available in EDI data deposit.
- `iland_exe/`: iLand executable used in this study. Will be available in EDI data deposit.
- `processed_data/`: Any data altered from its raw format. Includes data derived from field measurements; lookup tables modified from TRY Plant Trait Database and other sources; species distribution predictor selection, model fits, and evaluation; processed iLand model outputs; and species distribution model predictions. Will be available in EDI data deposit.
- `scripts/`: R and bash scripts for reproducing results.

## scripts

Scripts are named in order (`step01_`, etc.). `.Rmd` scripts rely on external or intermediate inputs not included in this deposit. `.R` scripts can be rerun with data included in deposit.

Script descriptions:

- `step01_data_prep_cleaning.Rmd`: Prepares response variables (occurrence and percent cover) from field data on understory plants collected in Berchtesgaden National Park. Subsets list to individual species occurring in at least 5 plots or individual species with better model fits (AUC &gt; 0.7) and evaluates how well these represents full dataset. Evaluates how well Ellenberg Indicator Values represent environmental gradients and assigns individual species to Plant Functional Types. Outputs in `processed_data/species_lookup/` and `processed_data/understory_model/`.
- `step02a_predictor_selection.Rmd`: Compiles species distribution model predictors representing climate, forest, and soil conditions. Evaluates collinearity and selects final predictors. Outputs in `processed_data/gis/` and `processed_data/understory_model/`.
- `step02b_forest_subset_predictor_selection.R`: Identifies most important forest predictors for building final species distribution models.
- `step03a_fit_sdms.R`: Fits species distribution models to predict individual species presence and absence. Includes evaluation, final fits, and variable importance and response curve (partial plot) outputs. Outputs in `processed_data/understory_model/`.
- `step03b_fit_mems.R`: Fits macroecological model for understory cover (evaluation and final fit). Also fits macroecological model for total species richness, for use in evaluation comparison with sum of individual species SDM predictions. Outputs in `processed_data/understory_model/`.
- `step03c_sdm_corrections.R`: Fits and evaluates species richness corrections based on Zurell et al. 2020 J Biogeography and Calabrese et al. 2014 Methods in Ecology and Evolution. Also applies probability ranking rule. Outputs in `processed_data/understory_model/`.
- `step04_model_evaluation.R`: Performs model evaluations on holdout test data. Includes comparison of individual SDMs with macroecological model predictions, as well as with and without richness corrections. Also evaluates variable importance and partial plots by Plant Functional Type. Outputs in `processed_data/understory_model/`.
- `step05_map_understory_iLand.Rmd`: Derives predictors at 10 m resolution for the full forested area in Berchtesgaden National Park based on historical climate, soil, and initial conditions in the contemporary iLand 2020 landscape. Uses these predictors to map understory plant species and total cover throughout the landscape. Individual species predictions not included in this data deposit, but aggregate plot-level values are (e.g., species richness). Outputs in `processed_data/iland_sdm_predictions/` and `processed_data/understory_model/final_fits/predictor_scaling/`.

To recreate iLand outputs used as inputs for this script, open the `historical_2020.xml` project file in iLand, save the light influence field grid at year 0 at 2 m height:

```
Globals.gridToFile('lifc','output/lif_baseline2020.asc',2);
```

Then, run the model for 1 year to get tree outputs at year 0.

- `step06_run_iland_local.sh`: Reruns iLand simulations used in this study. Must be run from home directory. Make sure to set `path` to location of folder with iLand exe. Also automatically runs R script to pre-process outputs. Outputs will be saved in `iland/output/`.

Run this script in bash with 10 replicates:

```
bash scripts/step06_run_iland_local.sh 10
```

- `step07_eval_iland_runs.Rmd`: Evaluates climate, disturbances, and forest change in the iLand simulation runs.
- `step08_future_climate_site_predictors.Rmd`: Compiles climate and site predictors for future SDM predictions at 10 m resolution. Outputs in `processed_data/iland_clim_soils/future_predictors/`.
- `step09_map_future_understory_iland.Rmd`: Predicts future understory plant communities based on future climate and forest change scenarios. Generates `.csv` summaries of predictions and creates rasters with average predictions at 10 m resolution for focal scenarios. Takes a random sample of 1000 points and predicts understory plant communities using either future v. historical climate or future v. historical forest conditions, used for evaluation of the relative importance of these two drivers. Outputs in `processed_data/iland_sdm_predictions/`.
- `step10_data_analysis.R`: Analyzes understory plant community change to answer manuscript research questions. Outputs in `analysis/` and `processed_data/iland_sdm_predictions/future_prediction_summaries/`.
