library(tidyverse)
library(sf)

#Creating Tibbles

captures <- read_csv("data/captures.csv")

OutlierId <- list("BI040801 XW274 ZF-310 X12-044", "BI406977 Z15-247 +", "BI407276 Z15-546 +")

capturesFemale <- captures %>%
  filter(sex == "Hona") %>%
  group_by(id) %>%
  mutate(meanlon = mean(lon), meanlat = mean(lat), n = n()) %>%
  group_by(id) %>%
  mutate(DistanceFromMean = sqrt((lon - meanlon)^2+(lat - meanlat)^2) %>% round(0)) %>%
  select(id, sex,  meanlon, meanlat, DistanceFromMean, n) %>%
  top_n(1, DistanceFromMean) %>%
  distinct()

capturesMale <- captures %>%
  filter(sex == "Hane") %>%
  group_by(id) %>%
  mutate(meanlon = mean(lon), meanlat = mean(lat), n = n()) %>%
  group_by(id) %>%
  mutate(DistanceFromMean = sqrt((lon - meanlon)^2+(lat - meanlat)^2) %>% round(0)) %>%
  select(id, sex, meanlon, meanlat, DistanceFromMean, n) %>%
  top_n(1, DistanceFromMean) %>%
  distinct()

map <- st_read("LanSweref99TM/Lan_Sweref99TM_region.shp", quiet = TRUE) %>% 
  mutate(Inventering = case_when(LnKod == 25 ~ "2021",
                                 LnKod %in% 22:23 ~ "2020",
                                 LnKod == 24 ~ "2019",
                                 LnKod %in% 20:21 ~ "2017"))


#Region Functions

RegionRatio <- function(mu1, mu2, r){
  #Function that takes the coordinates of a midpoint and a radius to create a circle and then calculates the ratio of the circle that is contained inside each of the 4 areas or outside.
  Territory <- st_point(c(mu1, mu2)) %>%
    st_buffer(dist = r)
  
  TotalArea <- st_area(Territory)
  
  #Defining the Regions
  Region1 <- st_union(map$geometry[[16]], map$geometry[[17]])
  Region2 <- map$geometry[[20]]
  Region3 <- st_union(map$geometry[[18]], map$geometry[[19]])
  Region4 <- map$geometry[[21]]
  
  #Helper function that returns the area of the intersection of the two arguments divided by the area of the territory circle.
  IntersectionRatio <- function(Region){
    st_intersection(Territory, Region) %>%
      st_area() / TotalArea
  }
  
  #Calculating the ratios.
  R1 <- IntersectionRatio(Region1)
  R2 <- IntersectionRatio(Region2)
  R3 <- IntersectionRatio(Region3)
  R4 <- IntersectionRatio(Region4)
  Sweden <- (R1 + R2 + R3 + R4)
  Other <- 1 - Sweden
  
  list(Region1=R1, Region2=R2, Region3=R3, Region4=R4, RegionSweden=Sweden, RegionOther=Other) %>%
    as_tibble() %>%
    round(3) %>%
    return()
}

#Special Tibbles

#Calculate Radius for each sex

RadiusFemale <- capturesFemale %>%
  filter(!(id %in% OutlierId)) %>%
  arrange(desc(DistanceFromMean)) %>%
  ungroup() %>%
  select(DistanceFromMean) %>%
  head(1) %>%
  .[[1]]
  
RadiusMale <- capturesMale %>%
  filter(!(id %in% OutlierId)) %>%
  arrange(desc(DistanceFromMean)) %>%
  ungroup() %>%
  select(DistanceFromMean) %>%
  head(1) %>%
  .[[1]]

#Create Region Tibbles

RegionFemale <- capturesFemale %>%
  select(id, sex, meanlon, meanlat, n) %>%
  mutate(Ratio = RegionRatio(meanlon, meanlat, RadiusFemale)) %>%
  mutate(Region1 = Ratio[[1]], Region2 = Ratio[[2]], Region3 = Ratio[[3]], Region4 = Ratio[[4]], RegionSweden = Ratio[[5]], RegionOther = Ratio[[6]]) %>%
  select(-Ratio)

RegionMale <- capturesMale %>%
  select(id, sex, meanlon, meanlat, n) %>%
  mutate(Ratio = RegionRatio(meanlon, meanlat, RadiusMale)) %>%
  mutate(Region1 = Ratio[[1]], Region2 = Ratio[[2]], Region3 = Ratio[[3]], Region4 = Ratio[[4]], RegionSweden = Ratio[[5]], RegionOther = Ratio[[6]]) %>%
  select(-Ratio)


#MAXIMUM LIKELIHOOD

MaxLikelihoodGenerator <- function(Tibble){
  
  function(lambda){
    Tibble %>%
      mutate(Likelihood = (RegionSweden*lambda) ^ n * exp(-(RegionSweden*lambda)) / factorial(n)) %>%
      mutate(Likelihood = -log(Likelihood)) %>%
      ungroup() %>%
      select(Likelihood) %>%
      summarise(Sum = sum(Likelihood)) %>%
      .[[1]]
  } %>%
    return()
}

MaxLikelihoodGeneratorZeroTrunc <- function(Tibble){
  
  function(lambda){
    Tibble %>%
      mutate(Likelihood = (RegionSweden*lambda) ^ n / ((exp(RegionSweden*lambda)-1) * factorial(n))) %>%
      mutate(Likelihood = -log(Likelihood)) %>%
      ungroup() %>%
      select(Likelihood) %>%
      summarise(Sum = sum(Likelihood)) %>%
      .[[1]]
  } %>%
    return()
}

FUNC <- MaxLikelihoodGenerator(RegionFemale)
FUNC2 <- MaxLikelihoodGeneratorZeroTrunc(RegionFemale)
FUNC3 <- MaxLikelihoodGenerator(RegionMale)
FUNC4 <- MaxLikelihoodGeneratorZeroTrunc(RegionMale)

LambdaFemale1 <- optim(14 ,f=FUNC1, method = "BFGS", hessian = TRUE)$par[[1]]

LambdaFemale2 <- optim(14 ,f=FUNC2, method = "BFGS", hessian = TRUE)$par[[1]]

LambdaMale2 <- optim(14 ,f=FUNC3, method = "BFGS", hessian = TRUE)$par[[1]]
  
LambdaMale2 <- optim(14 ,f=FUNC4, method = "BFGS", hessian = TRUE)$par[[1]]

#Estimating numbers of bears in each region

BearEstimator <- function(Tibble){
  
  Estimator <- function(Col){
   Tibble %>%
     ungroup() %>%
     select(Col) %>%
     map(sum) %>%
     .[[1]]
  }
  
  R1 <- Estimator("Region1")
  R2 <- Estimator("Region2")
  R3 <- Estimator("Region3")
  R4 <- Estimator("Region4")
  RS <- Estimator("RegionSweden")
  RO <- Estimator("RegionOther")
  Sex <- Tibble$sex[[1]]
  
  list(sex = Sex, Region1 = R1, Region2 = R2, Region3 = R3, Region4 = R4, RegionSweden = RS, RegionOther = RO)
}

FemaleEstimate <- BearEstimator(RegionFemale) %>%
  as_tibble()
MaleEstimate <- BearEstimator(RegionMale) %>%
  as_tibble()

Estimate <- full_join(FemaleEstimate, MaleEstimate) %>%
  select(-sex) %>%
  map_dfr(sum) %>%
  mutate(sex = "Total") %>%
  full_join(FemaleEstimate) %>%
  full_join(MaleEstimate)

FemaleAverageSweden<- RegionFemale %>%
  ungroup() %>%
  select(RegionSweden) %>%
  map(mean) %>%
  .[[1]]
  
MaleAverageSweden <- RegionMale %>%
  ungroup() %>%
  select(RegionSweden) %>%
  map(mean) %>%
  .[[1]]



