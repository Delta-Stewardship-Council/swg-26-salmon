setwd(this.path::here()); setwd('..')
library(ggplot2)

source("scripts/function_download_process_telemetry.R")
source("scripts/function_basic_survival_regression.R")
source("scripts/functions_model_setup.R")
source("scripts/functions_misc.R")

my_trib = 'feather'
my_season = 'spring'
speed_limit = 120
censor_upstream = 3

# download and process telemetry data (takes a minute)
alldat = download_process_telemetry(trib = my_trib, season = my_season, return_type = 'all',
                                    censor_upstream = censor_upstream, speed_limit = speed_limit)

dat_fish_recv = alldat$dat_fish_recv #full detection history for all fish
fishdat = alldat$fishdat #summary of fish and tributary survival
recvdat = alldat$recvdat #full receiver info
recvs_location_summary = alldat$recvs_location_summary #summary by receiver location
recvs_general_summary = alldat$recvs_general_summary #summary by receiver general location

rm(alldat) #remove big object


# make network plot
network_plot = plot_river_network(data=dat_fish_recv, filter_perc = 10,
                                  plot_only_consecutive = T, plot_all_names = T,
                                  save_dir = 'figures')
network_plot

# create list of general receiver locations defining each reach
reaches_list = list()
reaches_list[[1]] = c("Gridley", "FR_Gridley_Rel_DS", "FR_Gridley_Rel") #start of first reach
reaches_list[[2]] = c("BoydsPump", 'FR_Boyds_Rel_Rec', 'FR_Boyds_Rel') #start of second reach
reaches_list[[3]] = c("Blw_FRConf", "FR_Verona_Bot", "SutterWestSacRiver") #start of third reach... etc.
reaches_list[[4]] = c("TowerBridge", 'I80-50_Br') 
reaches_list[[5]] = c("Sac_Rio_Vista_1", "Sac_Rio_Vista_2", "RioVistaBr")
reaches_list[[6]] = c("BeniciaW", "BeniciaE")
reaches_list[[7]] = c("GoldenGateW", "GoldenGateE")

plot_reaches(recvs_location_summary, reaches_list) #quick plot of the reaches

# check receiver coverage 
recvs_coverage = get_reach_coverage(recvdat=recvdat, fishdat=fishdat, reaches_list=reaches_list, minDays = 45)
    # for each receiver group / fish release group pair,
    # complete_coverage = T indicates that at least one receiver is present at the site from [release date, release date + X], where e.g. X = minDays
    # partial_coverage = T (with complete_coverage = F) indicates that some, but not all of the date range [release date, release date + X],
    #                       has one or more receivers. Should be manually inspected to see if the coverage is sufficient or not.
    # parital_coverage = F indicates that no receivers were present at the site anytime during [release date, release date + X]. 
    #                        This site should definitely be excluded from the survival model for that year (e.g., fixing detection prob to 0)

recvs_coverage_check = subset(recvs_coverage, partial_coverage == T & complete_coverage == F)
recvs_coverage_exclude = subset(recvs_coverage, partial_coverage == F)


# corresponding to reaches_list above, manually dictate which receiver sites should *not* be included for a specific study year
excl_reaches_studies_list = list()
excl_reaches_studies_list[[1]] = NA 
excl_reaches_studies_list[[2]] = NA
excl_reaches_studies_list[[3]] = NA 
excl_reaches_studies_list[[4]] = NA
excl_reaches_studies_list[[5]] = c("FR_Spring_2013") 
excl_reaches_studies_list[[6]] = c("FR_Spring_2015")
excl_reaches_studies_list[[7]] = c("FR_Spring_2014", "FR_Spring_2015") #maybe 2013 too? coverage began a week after fish release

# aggregate fish detection history by reach








### Basic analyses ###

## fish length vs weight
ggplot(fishdat) + geom_point(aes(x = fish_length, y = fish_weight), alpha = .3) +
  theme_classic() + labs(x = 'Fish length', y = 'Fish weight')
out_lw = lm(data = fishdat, fish_weight ~ fish_length)
summary(out_lw)
plot(out_lw)
fishdat$expected_weight = predict(out_lw, newdata = fishdat)
fishdat$weight_ratio = fishdat$fish_weight/fishdat$expected_weight

fishdat$fish_condition =  fishdat$fish_weight/(fishdat$fish_length^3)
fishdat$release_location = factor(fishdat$release_location) #GAM requirement
fishdat$release_doy = as.numeric(strftime(fishdat$fish_release_date, format = "%j"))



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
  theme_classic() + labs(x = 'Fish condition', y = 'Survival probability', color = 'Release Site')  +
  ylim(c(0, 1))

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
  theme_classic() + labs(x = 'Fish condition', y = 'Survival probability', color = 'Release Site') +
  ylim(c(0, 1))



## Fit GLM/GAM by release date
out3 = basic_survival_regression(data = fishdat,
                                 y = 'survived.trib',
                                 x = c("release_location", 'release_doy'), 
                                 mod_type = 'GLM')
summary(out3$res)
preds_df3 = out3$predictions

ggplot(preds_df3) + 
  geom_line(aes(x = release_doy, y = prediction, color = release_location)) +
  theme_classic() + labs(x = 'Release day of year', y = 'Survival probability', color = 'Release Site') +
  ylim(c(0, 1))

ggplot(subset(preds_df3, fish_length %in% c(75, 91.5, 108) & release_location == 'FR_Boyds_Rel')) + 
  geom_line(aes(x = release_doy, y = prediction, color = factor(fish_length))) +
  theme_classic() + labs(x = 'Release day of year', y = 'Survival probability', color = 'Length') +
  ylim(c(0, 1))
