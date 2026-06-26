setwd(this.path::here()); setwd('..')
library(ggplot2)

source("scripts/function_download_process_telemetry.R")
source("scripts/function_basic_survival_regression.R")
source("scripts/functions_misc.R")

my_trib = 'feather'
my_season = 'spring'

# download and process telemetry data (takes a minute)
alldat = download_process_telemetry(trib = my_trib, season = my_season, return_type = 'all')
fishdat = alldat$fishdat; dat_fish_recv = alldat$dat_fish_recv
rm(alldat) #remove big object
fishdat$fish_condition =  fishdat$fish_weight/fishdat$fish_length^3
fishdat$release_location = factor(fishdat$release_location) #GAM requirement

# make network plot
network_plot = plot_river_network(data=dat_fish_recv, filter_perc = 10, save_dir = 'figures')
network_plot

### Fit GLM logistic model
out = basic_survival_regression(data = fishdat,
                                y = 'survived.trib',
                                x = c("release_location", "fish_condition"), 
                                mod_type = 'GLM')
out$res #raw model output
summary(out$res)

preds_df = out$predictions

ggplot(preds_df) + 
  geom_line(aes(x = fish_condition, y = prediction, color = release_location)) +
  theme_classic()

### Fit GAM 
out2 = basic_survival_regression(data = fishdat,
                                y = 'survived.trib',
                                x = c("release_location", "fish_condition"), 
                                mod_type = 'GAM')
out2$res #raw model output
summary(out2$res)

preds_df2 = out2$predictions

ggplot(preds_df2) + 
  geom_line(aes(x = fish_condition, y = prediction, color = release_location)) +
  theme_classic()
