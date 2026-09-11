
# Set Up

require(rerddap); require(tidyr); require(dplyr); library(geosphere); library(ggplot2); library(this.path)
setwd(this.path::here()); setwd('..')
library(ggplot2)

source("scripts/function_download_process_telemetry.R")
source("scripts/function_basic_survival_regression.R")
source("scripts/functions_misc.R")


# Download data

#run download_process_telemetry function on american river
download_process_telemetry(trib = 'american', season = NULL, save_dir = "data_processed", return_type = "all", censor_upstream = 5, speed_limit = 120, waterfall_plots = F)



