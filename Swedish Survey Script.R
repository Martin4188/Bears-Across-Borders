
library(tidyverse)
library(knitr)
library(sf)
library(polyCub)
library(nngeo)




captures <- read_csv("data/captures.csv") %>%
  group_by(id, year) %>%
  rename(sample_lon = lon,
         sample_lat = lat) %>%
  mutate(center_lon = mean(sample_lon),
         center_lat = mean(sample_lat),
         n = n()) %>%
  ungroup()


Outliers <- captures %>%
  select(id, center_lat, center_lon) %>%
  distinct() %>%
  rowwise() %>%
  mutate(Inside1 = st_within(st_point(c(center_lon, center_lat)), Region1) %>% list()) %>%
  mutate(Inside2 = st_within(st_point(c(center_lon, center_lat)), Region2) %>% list()) %>%
  mutate(Inside3 = st_within(st_point(c(center_lon, center_lat)), Region3) %>% list()) %>%
  mutate(Inside4 = st_within(st_point(c(center_lon, center_lat)), Region4) %>% list()) %>%
  mutate(Inside1 = Inside1[1]) %>%
  mutate(Inside2 = Inside2[1]) %>%
  mutate(Inside3 = Inside3[1]) %>%
  mutate(Inside4 = Inside4[1]) %>%
  filter(length(Inside1) == 0,
         length(Inside2) == 0,
         length(Inside3) == 0,
         length(Inside4) == 0)

captures <- captures %>%
  filter(!(id %in% Outliers$id),
         !(id == "BI406882 Z15-152" & year == "2016"))


map <- st_read("LanSweref99TM/Lan_Sweref99TM_region.shp", quiet = TRUE) %>% 
  mutate(Inventering = case_when(LnKod == 25 ~ "2021",
                                 LnKod %in% 22:23 ~ "2020",
                                 LnKod == 24 ~ "2019",
                                 LnKod %in% 20:21 ~ "2017")) %>%
  mutate(geometry = st_remove_holes(geometry))



#Function

EstimatorFunction <- function(Year, Sex, tibble = captures){
  
  
  
  tibble <- tibble %>%
    filter(year == Year,
           sex == Sex)
  
  map <- st_read("LanSweref99TM/Lan_Sweref99TM_region.shp", quiet = TRUE) %>%
    mutate(geometry = st_remove_holes(geometry))
  
  
  RegionList <- list(st_union(map$geometry[[16]], map$geometry[[17]]), 
                     map$geometry[[20]], 
                     st_union(map$geometry[[18]], map$geometry[[19]]), 
                     map$geometry[[21]])
  
  RegionNumber <- case_when(Year == "2017" ~ 1,
                            Year == "2019" ~ 2,
                            Year == "2015" ~ 3,
                            Year == "2020" ~ 3,
                            Year == "2016" ~ 4)
  
  Region <- RegionList[[RegionNumber]]
  
  #Estimate sigma
  sigmaEstimate <- tibble %>%
    filter(n != 1) %>%
    arrange(id) %>%
    mutate(Longitude = (sample_lon - center_lon)^2,
           Latitude = (sample_lat - center_lat)^2) %>%
    pivot_longer(Longitude:Latitude, names_to = "direction", values_to = "deviation") %>%
    select(id, n, direction, deviation) %>%
    group_by(id, direction) %>%
    mutate(deviation = sum(deviation)) %>%
    distinct() %>%
    ungroup() %>%
    summarise(deviation = sum(deviation),
              "n-1" = sum(n - 1)) %>%
    mutate(s = sqrt(deviation / `n-1`)) %>%
    select(s) %>%
    .[[1]]
  
  #Function that Generates a maximum likelihood function that can be optimized.
  MaxLikelihoodGeneratorZeroTrunc <- function(Tibble){
    
    function(lambda){
      Tibble %>%
        mutate(Likelihood = (IRatio * lambda) ^ n / ((exp(IRatio * lambda) - 1) * factorial(n))) %>%
        mutate(Likelihood = -log(Likelihood)) %>%
        ungroup() %>%
        select(Likelihood) %>%
        summarise(Sum = sum(Likelihood)) %>%
        .[[1]]
    } %>%
      return()
  }
  
  
  #Function to calculate InsideRatio. MOVE SQUARE
  InsideRatioFunction <- function(mu1, mu2, sigma, region, Precision = 50){
    
    region <- region - c(mu1, mu2)
    
    #polyCub requires a zero mean normal distribution, as such i translate the 
    #middlepoint so that it becomes (0, 0) and perform the same translation on the Region.
    
    Binormal <- function (s, S = sigma)
    {
      exp(-rowSums(s^2)/2/S^2) / (2*pi*S^2)
    }
    
    Intersection <- st_point(c(0, 0)) %>% 
      st_buffer(dist = 3 * sigma) %>%
      st_intersection(region)
    
    
    
    if(st_is_empty(Intersection)){
      return(0)
    }
    
    Intersection %>%
      polyCub.SV(f = Binormal, nGQ = Precision) %>%
      return()
  }
  
  #function that adjust the center point.
  center_fit <- function(data){
    
    n <- data$n[1]
    CenterLon <- data$center_lon[1]
    CenterLat <- data$center_lat[1]
    
    if(n == 1){
      return(c(CenterLon, CenterLat))
    }
    
    
    Binormal <- function (s, S = sigma)
    {
      exp(-rowSums(s^2)/2/S^2) / (2*pi*S^2)
    }
    
    LikelihoodFunction <- function(MidPoint){
      data %>%
        ungroup() %>%
        rowwise() %>%
        mutate(likelihood = Binormal(matrix(c(sample_lon - MidPoint[1], sample_lat - MidPoint[2]), 1, 2), sigmaEstimate),
               likelihood = likelihood / InsideRatioFunction(MidPoint[1], MidPoint[2], sigmaEstimate, Region),
               likelihood = -log(likelihood)) %>%
        ungroup() %>%
        summarise(likelihood = sum(likelihood)) %>%
        select(likelihood) %>%
        .[[1]]
    }
    
    fit <-  optim(par = c(CenterLon, CenterLat), 
                  fn = LikelihoodFunction, 
                  method = "BFGS")
    
    fit$par %>%
      return()
  }
  
  
  AdjustedCenter <- tibble %>% 
    nest_by(id) %>%
    mutate(center_point = list(center_fit(data))) %>%
    select(id, center_point) %>%
    mutate(center_lon = center_point[1],
           center_lat = center_point[2]) %>%
    select(-center_point)
  
  
  OptimTibble <- tibble %>%
    select(-center_lon, -center_lat, -sample_lon, -sample_lat, -date) %>%
    full_join(AdjustedCenter, by = "id") %>%
    distinct() %>%
    rowwise() %>%
    mutate(IRatio = InsideRatioFunction(center_lon, center_lat, sigmaEstimate, Region))
  
  fit <- MaxLikelihoodGeneratorZeroTrunc(OptimTibble) %>%
    optim(lambda, ., 
          hessian = TRUE, 
          method = "L-BFGS-B",
          lower = 0.1)
  
  NObs <- OptimTibble %>%
    ungroup() %>%
    summarise(n=n()) %>%
    .[[1]]
  
  InsideRatioSum <- OptimTibble %>%
    ungroup() %>%
    summarise(sum = sum(IRatio)) %>%
    .[[1]]
  
  RemoveOutsiders <- OptimTibble %>%
    mutate(Inside = st_within(st_point(c(center_lon, center_lat)), Region) %>% list()) %>%
    mutate(Inside = Inside[1]) %>%
    filter(length(Inside) != 0) %>%
    ungroup() %>%
    summarise(n = n()) %>%
    .[[1]]
  
  
  Temp <- list(NObserved = NObs, 
               WithoutOutsiders = RemoveOutsiders,
               RatioSum = InsideRatioSum, 
               LambdaHat = fit$par,
               SigmaHat = sigmaEstimate, 
               NormalizedSigma = (sigmaEstimate * 2) / sqrt(st_area(Region))) %>%
    as_tibble()
  
  return(Temp)
  
}


SwedishSurvey <- expand_grid(sex = c("Hane", "Hona"), year = c("2015", "2016", "2017", "2019", "2020")) %>%
  rowwise() %>%
  mutate(results = EstimatorFunction(Year = year, Sex = sex)) 


SwedishSurvey1 <- expand_grid(sex = c("Hane"), year = c("2015")) %>%
  rowwise() %>%
  mutate(results = EstimatorFunction(Year = year, Sex = sex))

SwedishSurvey2 <- expand_grid(sex = c("Hane"), year = c("2016")) %>%
  rowwise() %>%
  mutate(results = EstimatorFunction(Year = year, Sex = sex))

SwedishSurvey3 <- expand_grid(sex = c("Hane"), year = c("2017")) %>%
  rowwise() %>%
  mutate(results = EstimatorFunction(Year = year, Sex = sex))

SwedishSurvey4 <- expand_grid(sex = c("Hane"), year = c("2019")) %>%
  rowwise() %>%
  mutate(results = EstimatorFunction(Year = year, Sex = sex))

SwedishSurvey5 <- expand_grid(sex = c("Hane"), year = c("2020")) %>%
  rowwise() %>%
  mutate(results = EstimatorFunction(Year = year, Sex = sex))

SwedishSurvey6 <- expand_grid(sex = c("Hona"), year = c("2015")) %>%
  rowwise() %>%
  mutate(results = EstimatorFunction(Year = year, Sex = sex))

SwedishSurvey7 <- expand_grid(sex = c("Hona"), year = c("2016")) %>%
  rowwise() %>%
  mutate(results = EstimatorFunction(Year = year, Sex = sex))

SwedishSurvey8 <- expand_grid(sex = c("Hona"), year = c("2017")) %>%
  rowwise() %>%
  mutate(results = EstimatorFunction(Year = year, Sex = sex))

SwedishSurvey9 <- expand_grid(sex = c("Hona"), year = c("2019")) %>%
  rowwise() %>%
  mutate(results = EstimatorFunction(Year = year, Sex = sex))

SwedishSurvey10 <- expand_grid(sex = c("Hona"), year = c("2020")) %>%
  rowwise() %>%
  mutate(results = EstimatorFunction(Year = year, Sex = sex))





expand_grid(sex = c("Hona"), year = c("2017")) %>%
  rowwise() %>%
  mutate(results = EstimatorFunction(Year = year, Sex = sex)) 


SwedishSurvey %>%
  unpack(results)

SwedishSurvey %>%
  unpack(results) %>%
  write_csv("data/SwedishSurveyEstimate")










#Function2

EstimatorFunctionNoSex <- function(Year, tibble = captures){
  
  
  
  tibble <- tibble %>%
    filter(year == Year)
  
  map <- st_read("LanSweref99TM/Lan_Sweref99TM_region.shp", quiet = TRUE) %>%
    mutate(geometry = st_remove_holes(geometry))
  
  
  RegionList <- list(st_union(map$geometry[[16]], map$geometry[[17]]), 
                     map$geometry[[20]], 
                     st_union(map$geometry[[18]], map$geometry[[19]]), 
                     map$geometry[[21]])
  
  RegionNumber <- case_when(Year == "2017" ~ 1,
                            Year == "2019" ~ 2,
                            Year == "2015" ~ 3,
                            Year == "2020" ~ 3,
                            Year == "2016" ~ 4)
  
  Region <- RegionList[[RegionNumber]]
  
  #Estimate sigma
  sigmaEstimate <- tibble %>%
    filter(n != 1) %>%
    arrange(id) %>%
    mutate(Longitude = (sample_lon - center_lon)^2,
           Latitude = (sample_lat - center_lat)^2) %>%
    pivot_longer(Longitude:Latitude, names_to = "direction", values_to = "deviation") %>%
    select(id, n, direction, deviation) %>%
    group_by(id, direction) %>%
    mutate(deviation = sum(deviation)) %>%
    distinct() %>%
    ungroup() %>%
    summarise(deviation = sum(deviation),
              "n-1" = sum(n - 1)) %>%
    mutate(s = sqrt(deviation / `n-1`)) %>%
    select(s) %>%
    .[[1]]
  
  #Function that Generates a maximum likelihood function that can be optimized.
  MaxLikelihoodGeneratorZeroTrunc <- function(Tibble){
    
    function(lambda){
      Tibble %>%
        mutate(Likelihood = (IRatio * lambda) ^ n / ((exp(IRatio * lambda) - 1) * factorial(n))) %>%
        mutate(Likelihood = -log(Likelihood)) %>%
        ungroup() %>%
        select(Likelihood) %>%
        summarise(Sum = sum(Likelihood)) %>%
        .[[1]]
    } %>%
      return()
  }
  
  
  #Function to calculate InsideRatio. MOVE SQUARE
  InsideRatioFunction <- function(mu1, mu2, sigma, region, Precision = 50){
    
    region <- region - c(mu1, mu2)
    
    #polyCub requires a zero mean normal distribution, as such i translate the 
    #middlepoint so that it becomes (0, 0) and perform the same translation on the Region.
    
    Binormal <- function (s, S = sigma)
    {
      exp(-rowSums(s^2)/2/S^2) / (2*pi*S^2)
    }
    
    Intersection <- st_point(c(0, 0)) %>% 
      st_buffer(dist = 3 * sigma) %>%
      st_intersection(region)
    
    
    
    if(st_is_empty(Intersection)){
      return(0)
    }
    
    Intersection %>%
      polyCub.SV(f = Binormal, nGQ = Precision) %>%
      return()
  }
  
  #function that adjust the center point.
  center_fit <- function(data){
    
    n <- data$n[1]
    CenterLon <- data$center_lon[1]
    CenterLat <- data$center_lat[1]
    
    if(n == 1){
      return(c(CenterLon, CenterLat))
    }
    
    
    Binormal <- function (s, S = sigma)
    {
      exp(-rowSums(s^2)/2/S^2) / (2*pi*S^2)
    }
    
    LikelihoodFunction <- function(MidPoint){
      data %>%
        ungroup() %>%
        rowwise() %>%
        mutate(likelihood = Binormal(matrix(c(sample_lon - MidPoint[1], sample_lat - MidPoint[2]), 1, 2), sigmaEstimate),
               likelihood = -log(likelihood)) %>%
        ungroup() %>%
        summarise(likelihood = sum(likelihood)) %>%
        mutate(likelihood = likelihood / InsideRatioFunction(MidPoint[1], MidPoint[2], sigmaEstimate, Region)) %>%
        select(likelihood) %>%
        .[[1]]
    }
    
    fit <-  optim(par = c(CenterLon, CenterLat), 
                  fn = LikelihoodFunction, 
                  method = "BFGS")
    
    fit$par %>%
      return()
  }
  
  
  AdjustedCenter <- tibble %>% 
    nest_by(id) %>%
    mutate(center_point = list(center_fit(data))) %>%
    select(id, center_point) %>%
    mutate(center_lon = center_point[1],
           center_lat = center_point[2]) %>%
    select(-center_point)
  
  
  OptimTibble <- tibble %>%
    select(-center_lon, -center_lat, -sample_lon, -sample_lat, -date) %>%
    full_join(AdjustedCenter, by = "id") %>%
    distinct() %>%
    rowwise() %>%
    mutate(IRatio = InsideRatioFunction(center_lon, center_lat, sigmaEstimate, Region))
  
  fit <- MaxLikelihoodGeneratorZeroTrunc(OptimTibble) %>%
    optim(lambda, ., 
          hessian = TRUE, 
          method = "L-BFGS-B",
          lower = 0.1)
  
  NObs <- OptimTibble %>%
    ungroup() %>%
    summarise(n=n()) %>%
    .[[1]]
  
  InsideRatioSum <- OptimTibble %>%
    ungroup() %>%
    summarise(sum = sum(IRatio)) %>%
    .[[1]]
  
  RemoveOutsiders <- OptimTibble %>%
    mutate(Inside = st_within(st_point(c(center_lon, center_lat)), Region) %>% list()) %>%
    mutate(Inside = Inside[1]) %>%
    filter(length(Inside) != 0) %>%
    ungroup() %>%
    summarise(n = n()) %>%
    .[[1]]
  
  
  Temp <- list(NObserved = NObs, 
               WithoutOutsiders = RemoveOutsiders,
               RatioSum = InsideRatioSum, 
               LambdaHat = fit$par,
               SigmaHat = sigmaEstimate, 
               NormalizedSigma = (sigmaEstimate * 2) / sqrt(st_area(Region))) %>%
    as_tibble()
  
  return(Temp)
  
}


SwedishSurveyNoSex <- expand_grid(year = c("2015", "2016", "2017", "2019", "2020")) %>%
  rowwise() %>%
  mutate(results = EstimatorFunctionNoSex(Year = year)) 












NaiveEstimate <- function(YEAR, SEX){
  
  
  captures <- read_csv("data/captures.csv") %>%
    group_by(id, year) %>%
    rename(sample_lon = lon,
           sample_lat = lat) %>%
    mutate(center_lon = mean(sample_lon),
           center_lat = mean(sample_lat),
           n = n()) %>%
    ungroup() %>%
    filter(year == YEAR,
           sex == SEX)
  
  
  
  MaxLikelihoodGeneratorZeroTrunc <- function(Tibble){
    
    function(lambda){
      Tibble %>%
        mutate(Likelihood = (lambda) ^ n / ((exp(lambda) - 1) * factorial(n))) %>%
        mutate(Likelihood = -log(Likelihood)) %>%
        ungroup() %>%
        select(Likelihood) %>%
        summarise(Sum = sum(Likelihood)) %>%
        .[[1]]
    } %>%
      return()
  }
  
  N <- captures %>%
    select(id, n) %>%
    distinct() %>%
    summarise(n=n()) %>%
    .[[1]]
  
  
  LambdaEstimate <- captures %>% 
    select(id, n) %>%
    distinct() %>%
    MaxLikelihoodGeneratorZeroTrunc() %>%
    optim(lambda, ., 
          hessian = TRUE, 
          method = "L-BFGS-B",
          lower = 0.1) %>%
    .$par
  
  list(Estimate = N / (1 - exp(-LambdaEstimate)), LambdaHat = LambdaEstimate) %>%
    as_tibble() %>%
    return()
  
  
  
}

StandardEstimate <- expand.grid(year = c("2015", "2016", "2017", "2019", "2020"), sex = c("Hane", "Hona")) %>%
  rowwise() %>%
  mutate(Estimate = NaiveEstimate(year, sex)) %>%
  unpack(Estimate) %>%
  write_csv("data/StandardEstimate")







