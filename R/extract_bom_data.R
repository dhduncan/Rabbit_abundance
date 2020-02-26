# get rainfall data from the Bureau of Meteorology

library(bomrang)
library(dplyr)
library(lubridate)

# arrange weather stations by distance from a point roughly in the middle of the sites
bom_stats <- sweep_for_stations(latlon = c(-35.4,142))

# select only stations within 100 km. This is odd as only three stations come back in, so there is some filter being applied. Perhaps for active stations?
bom_stats <- bom_stats %>% filter(distance < 100)

# this also finds closest station to a defined point
get_historical(latlon = c(-35.4, 142), radius = 60, type = "rain")

get_historical(stationid = 076065, type = "rain")

patche_full <- get_historical(stationid = 077033, type = "rain")

patche <- patche_full %>%
  select(year, month, rain = rainfall) %>% 
  group_by(year, month) %>% 
  summarise(rain = sum(rain)) %>% 
  mutate(Deemed_date = ymd(paste(year, month, 15, sep = " "))
  )

# trim to include only the relevant time period - though before this construct some contextual graphs of the full time series.

patche <- filter(patche, year > 1972)


