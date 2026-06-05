#amer_exploratory_figures
library(tidyverse)
library(here)
library(ggridges)
library(here)
library(dataRetrieval)
library(zoo)

# Fall run catch timing by water year (ridgelines)

amer <- read_csv(here("data_processed", "amer_combined.csv")) %>%
  mutate(date = ymd(date))

# Fall run only, expand catch counts
fall_catch <- amer %>%
  filter(final.run == "Fall") %>%
  mutate(doy = yday(date)) %>%
  uncount(count)

ggplot(fall_catch, aes(x = doy, y = factor(water.year))) +
  geom_density_ridges(fill = "steelblue", alpha = 0.6, scale = 0.9) +
  scale_x_continuous(
    breaks = c(1, 32, 60, 91, 121, 152, 182, 213, 244, 274, 305, 335),
    labels = month.abb
  ) +
  scale_y_discrete(expand = expansion(add = c(0.2, 1.5))) +
  labs(x = NULL, y = "Water year",
       title = "Fall run catch timing by water year")


# fall run raw catch bar plot 

amer <- read_csv(here("data_processed", "amer_combined.csv")) %>%
  mutate(date = ymd(date))

fall_run_catch_by_month <- amer %>%
  filter(final.run == "Fall") %>%
  mutate(month = month(date, label = TRUE)) %>%
  group_by(month) %>%
  summarise(total_catch = sum(count), .groups = "drop") %>%
  ggplot(aes(month, total_catch)) +
  geom_col(fill = "steelblue") +
  scale_y_continuous(labels = scales::comma) +
  labs(x = NULL, y = "Total catch",
       title = "Total fall run catch by month, 2014–2025")

ggsave(here("figures", "fall_run_catch_by_month.png"),
       plot = fall_run_catch_by_month,
       width = 8, height = 5, dpi = 300)

# Fall run only, Average across years
average_amer <- amer %>%
  filter(final.run == "Fall") %>%
  mutate(doy = yday(date),
         week = floor(doy / 7) * 7) %>%
  group_by(water.year, week) %>%
  summarise(weekly_catch = sum(count), .groups = "drop") %>%
  group_by(week) %>%
  summarise(mean_catch = mean(weekly_catch),
            se = sd(weekly_catch) / sqrt(n()),
            .groups = "drop") %>%
  ggplot(aes(week, mean_catch)) +
  geom_ribbon(aes(ymin = pmax(mean_catch - se, 0), ymax = mean_catch + se),
              fill = "lightblue", alpha = 0.5) +
  geom_line(color = "steelblue", linewidth = 1) +
  scale_x_continuous(
    breaks = c(1, 32, 60, 91, 121, 152, 182, 213, 244, 274, 305, 335),
    labels = month.abb
  ) +
  labs(x = NULL, y = "Mean weekly catch (± SE)",
       title = "Average fall run catch timing")

ggsave(here("figures", "average_amer.png"),
       plot = average_amer,
       width = 8, height = 5, dpi = 300)



visit <- amer %>%
  group_by(trap.visit.id, water.vel) %>%
  summarise(catch = sum(count), .groups = "drop") %>%
  filter(!is.na(water.vel))

ggplot(visit, aes(water.vel, catch)) +
  geom_smooth(method = "gam", formula = y ~ s(x),
              color = "steelblue", fill = "lightblue") +
  scale_y_log10() +
  labs(x = "Water velocity (m/s)", y = "Catch per trap visit (log)",
       title = "Relationship between water velocity and catch")




# Total catch per month per water year, then average across years
library(here)

fall_monthly <- amer %>%
  filter(final.run == "Fall") %>%
  mutate(month = month(date, label = TRUE)) %>%
  group_by(water.year, month) %>%
  summarise(monthly_catch = sum(count), .groups = "drop") %>%
  group_by(month) %>%
  summarise(
    mean_catch = mean(monthly_catch),
    se = sd(monthly_catch) / sqrt(n()),
    .groups = "drop"
  )

fall_monthly_avg_plot <- ggplot(fall_monthly, aes(month, mean_catch)) +
  geom_col(fill = "steelblue") +
  geom_errorbar(aes(ymin = pmax(mean_catch - se, 0), ymax = mean_catch + se),
                width = 0.3) +
  labs(x = NULL, y = "Mean monthly catch (± SE)",
       title = "Average fall run catch by month")

ggsave(here("figures", "fall_run_avg_catch_by_month.png"),
       plot = fall_monthly_avg_plot,
       width = 8, height = 5, dpi = 300)


#discharge 

# Pull mean USGS daily discharge for the American River at Fair Oaks
flow <- readNWISdv(
  siteNumbers = "11446500",
  parameterCd = "00060",
  startDate = "2014-01-01",
  endDate = "2025-06-30"
) %>%
  renameNWISColumns() %>%
  select(date = Date, discharge_cfs = Flow)

# Filter to fall run, aggregate to trap visit, join discharge
visit_fall <- amer %>%
  filter(final.run == "Fall") %>%
  group_by(trap.visit.id, date) %>%
  summarise(catch = sum(count), .groups = "drop") %>%
  left_join(flow, by = "date") %>%
  filter(!is.na(discharge_cfs))

# 1. Simple discharge–catch relationship (fall run only)
ggplot(visit_fall, aes(discharge_cfs, catch)) +
  geom_smooth(method = "gam", formula = y ~ s(x),
              color = "steelblue", fill = "lightblue") +
  scale_x_log10() +
  scale_y_log10() +
  labs(x = "Discharge (cfs, log)", y = "Catch per visit (log)",
       title = "Fall run: discharge vs. catch")


flow <- readNWISdv(
  siteNumbers = "11446500",
  parameterCd = "00060",
  startDate = "2014-01-01",
  endDate = "2025-06-30"
) %>%
  renameNWISColumns() %>%
  select(date = Date, discharge_cfs = Flow)

# Fall run, one row per trap visit
visit_fall <- amer %>%
  filter(final.run == "Fall") %>%
  group_by(trap.visit.id, date) %>%
  summarise(catch = sum(count), .groups = "drop") %>%
  left_join(flow, by = "date") %>%
  filter(!is.na(discharge_cfs))

# Raw scatter, no transformation
ggplot(visit_fall, aes(discharge_cfs, catch)) +
  geom_point(alpha = 0.3) +
  labs(x = "Discharge (cfs)", y = "Catch per trap visit",
       title = "Fall run: raw catch vs. discharge")


#General trap info 

# Total visits
amer %>% distinct(trap.visit.id) %>% nrow()

# Visits by station
amer %>% distinct(trap.visit.id, station) %>% count(station)

# Visits by water year
amer %>% distinct(trap.visit.id, water.year) %>% count(water.year)

# Visits by run (note: a visit can catch multiple runs, so totals may exceed unique visits)
amer %>% distinct(trap.visit.id, final.run) %>% count(final.run)

# life stage catch vs discharge
visit_fall_stage <- amer %>%
  filter(final.run == "Fall") %>%
  group_by(trap.visit.id, date, life.stage) %>%
  summarise(catch = sum(count), .groups = "drop") %>%
  left_join(flow, by = "date") %>%
  filter(!is.na(discharge_cfs))

ggplot(visit_fall_stage, aes(discharge_cfs, catch)) +
  geom_smooth(method = "gam", formula = y ~ s(x)) +
  scale_x_log10() +
  scale_y_log10() +
  facet_wrap(~ life.stage)

# Fall run overall — distribution of catch across discharge

visit_fall %>%
  mutate(discharge_bin = cut(discharge_cfs,
    breaks = c(0, 500, 1000, 2000, 3000, 5000, 10000, Inf),
    labels = c("<500", "500-1000", "1000-2000", "2000-3000",
               "3000-5000", "5000-10000", ">10000"))) %>%
  group_by(discharge_bin) %>%
  summarise(
    total_catch = sum(catch),
    n_visits = n(),
    .groups = "drop"
  ) %>%
  mutate(pct_of_catch = round(100 * total_catch / sum(total_catch), 1)) %>%
  arrange(desc(pct_of_catch))

sum(amer$count)

#summary table of fall run juvenile chinook timing color coded

amer <- read_csv(here("data_processed", "amer_combined.csv")) %>%
  mutate(date = ymd(date))

fall_timing <- amer %>%
  filter(final.run == "Fall") %>%
  mutate(month = month(date, label = TRUE, abbr = FALSE)) %>%
  group_by(month) %>%
  summarise(total_catch = sum(count), .groups = "drop") %>%
  mutate(
    pct_of_catch = round(100 * total_catch / sum(total_catch), 1),
    timing = case_when(
      pct_of_catch >= 15 ~ "High",
      pct_of_catch >= 5  ~ "Medium",
      pct_of_catch >  0  ~ "Low",
      TRUE               ~ "None"
    )
  )

# fall_timing_table includes the gt() conversion
fall_timing_table <- fall_timing %>%
  gt() %>%
  tab_header(
    title = "Fall-run juvenile Chinook timing",
    subtitle = "Lower American River RST, water years 2014–2025"
  ) %>%
  cols_label(
    month = "Month",
    total_catch = "Total catch",
    pct_of_catch = "% of catch",
    timing = "Timing"
  ) %>%
  fmt_number(columns = total_catch, decimals = 0, use_seps = TRUE)

gtsave(fall_timing_table, here("figures", "fall_run_timing.png"))
