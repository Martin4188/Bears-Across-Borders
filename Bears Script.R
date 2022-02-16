library(tidyverse)
library(sf)

captures <- read_csv("data/captures.csv")

capturesFinal <- captures %>%
  filter(sex == "Hona") %>%
  group_by(id) %>%
  mutate(meanlon = mean(lon), meanlat = mean(lat)) %>%
  group_by(id) %>%
  mutate(DistanceFromMean = sqrt((lon - meanlon)^2+(lat - meanlat)^2))

map <- st_read("LanSweref99TM/Lan_Sweref99TM_region.shp", quiet = TRUE) %>% 
  mutate(Inventering = case_when(LnKod == 25 ~ "2021",
                                 LnKod %in% 22:23 ~ "2020",
                                 LnKod == 24 ~ "2019",
                                 LnKod %in% 20:21 ~ "2017"))


