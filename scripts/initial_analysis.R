setwd(this.path::here())
library(ggplot2)

source("function_download_process_telemetry.R")
source("function_basic_survival_regression.R")

# download and process telemetry data (takes a minute)
fishdat = download_process_telemetry(trib = 'feather', season = 'spring')
fishdat$fish_condition =  fishdat$fish_weight/fishdat$fish_length^3
fishdat$release_location = factor(fishdat$release_location) #GAM requirement

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
