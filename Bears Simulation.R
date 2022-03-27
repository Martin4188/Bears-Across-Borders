library(tidyverse)
library(sf)


#Estimating lambda without paying attention to the bears territory.
NaiveSimulationFunction <- function(mu, lambda, sigma){
  
  #Simulations
  N <- rpois(1, 9 * mu) # Total population
  bears <- tibble(id = 1:N,
                  center_lon = runif(N, -3, 3),
                  center_lat = runif(N, -3, 3),
                  n_samples = rpois(N, lambda)) %>%
    filter(abs(center_lon) <= 1 + 3 * sigma, abs(center_lat) <= 1 +3 * sigma) %>% #Filter bears whose vertical and/or horizontal distance from the middle square is larger than 3 sigma.
    mutate(Inside = if_else(abs(center_lon) <= 1 & abs(center_lat) <= 1, 1, 0)) #Dummy variable for being inside the central square.
  
  NTrue <- bears %>%
    select(Inside) %>%
    summarise(Inside = sum(Inside)) %>%
    .[[1]]
  
  samples <- bears %>%
    rowwise() %>%
    mutate(samples = list(tibble(sample_lon = rnorm(n_samples, center_lon, sigma),
                                 sample_lat = rnorm(n_samples, center_lat, sigma)))) %>%
    unnest(samples) %>%
    filter(abs(sample_lon) <= 1, abs(sample_lat) <= 1) %>%
    select(-n_samples) %>%
    group_by(id) %>%
    mutate(n = n())
  
  #Function that Generates a maximum likelihood function that can be optimized.
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
  
  
  #Calculate Mean number of samples.
  
  nMean <- samples %>%
    select(id, n) %>%
    distinct() %>%
    ungroup() %>%
    summarise(Mean = mean(n)) %>%
    .[[1]]
  
  #Calculate Observed number of bears
  
  NObs <- samples %>%
    ungroup() %>%
    select(id) %>%
    distinct() %>%
    summarise(n=n()) %>%
    .[[1]]
  
  #Calculate the number of observed bears whose true midpoint lies inside the middle square.
  
  NObsTrue <- samples %>%
    select(id, Inside) %>%
    distinct() %>%
    ungroup() %>%
    summarise(Sum = sum(Inside)) %>%
    .[[1]]
  
  #Using Optim with "L-BFGS-B"
  
  fit <- samples %>% 
    select(id, n) %>%
    distinct() %>%
    MaxLikelihoodGeneratorZeroTrunc() %>%
    optim(lambda, ., 
          hessian = TRUE, 
          method = "L-BFGS-B",
          lower = 0.1)
  
  
  Results <- list(MLE = fit$par,
                  MeanEstimate = nMean,
                  Fisher = as.numeric(sqrt(1 / fit$hessian)),
                  NTrue = NTrue,
                  NObs = NObs,
                  NObsTrue = NObsTrue) %>%
    as_tibble()
  
  return(Results)
  
}


set.seed(2022)
TEST <- expand_grid(mu = 500, lambda = 2:4, sigma = 1:3/10, sim = 1:10) %>% 
  rowwise() %>%
  mutate(results = NaiveSimulationFunction(mu, lambda, sigma)) %>%
  mutate(MLE = results[[1]], 
         Mean = results[[2]], 
         Fisher = results[[3]], 
         NTrue = results[[4]], 
         NObs = results[[5]], 
         NObsTrue = results[[6]]) %>%
  select(-results)






set.seed(2022)
simulation_study <- expand_grid(mu = 500, lambda = 2:4, sigma = 1:3/10, sim = 1:1000) %>% 
  rowwise() %>%
  mutate(results = NaiveSimulationFunction(mu, lambda, sigma)) %>%
  mutate(MLE = results[[1]], 
         Mean = results[[2]], 
         Fisher = results[[3]], 
         NTrue = results[[4]], 
         NObs = results[[5]], 
         NObsTrue = results[[6]]) %>%
  select(-results)

simulation_study %>%
  write.csv("data/NaiveSimulation2")



set.seed(2022)
simulation_study_small_sigma <- expand_grid(mu = 500, lambda = 2:4, sigma = 1:3/40, sim = 1:1000) %>% 
  rowwise() %>%
  mutate(results = NaiveSimulationFunction(mu, lambda, sigma)) %>%
  mutate(MLE = results[[1]], 
         Mean = results[[2]], 
         Fisher = results[[3]], 
         NTrue = results[[4]], 
         NObs = results[[5]], 
         NObsTrue = results[[6]]) %>%
  select(-results)







set.seed(2022)
LambdaEqualsTwo <- expand_grid(mu = 500, lambda = 2, sigma = 1:3/10, sim = 1:1000) %>% 
  rowwise() %>%
  mutate(results = NaiveSimulationFunction(mu, lambda, sigma)) %>%
  mutate(MLE = results[[1]], Mean = results[[2]], Fisher = results[[3]]) %>%
  select(-results)


LambdaEqualsTwo %>%
  write.csv("data/NaiveSimulation1Lambda2")











simulation_study %>%
  group_by(lambda) %>%
  summarise(MLE = mean(MLE), Mean = mean(Mean), Fisher = mean(Fisher))



NaiveTimingFunction <- function(lambda, sigma, mu){
  start_time <- Sys.time()
  NaiveSimulationFunction(mu, lambda, sigma)
  end_time <- Sys.time()
  return(end_time - start_time)
}


NaiveTimingFunction(500, 4, 0.1)
NaiveTimingFunction(500, 4, 0.3)
NaiveTimingFunction(500, 6, 0.1)
NaiveTimingFunction(500, 6, 0.3)





















































SimulationFunction <- function(mu, lambda, sigma){
  
  N <- rpois(1, 9 * mu) # Total population
  bears <- tibble(id = 1:N,
                  center_lon = runif(N, -3, 3),
                  center_lat = runif(N, -3, 3),
                  n_samples = rpois(N, lambda))
  
  samples <- bears %>%
    rowwise() %>%
    mutate(samples = list(tibble(sample_lon = rnorm(n_samples, center_lon, sigma),
                                 sample_lat = rnorm(n_samples, center_lat, sigma)))) %>%
    unnest(samples) %>%
    filter(abs(sample_lon) <= 1, abs(sample_lat) <= 1) %>%
    select(-n_samples) %>%
    group_by(id) %>%
    mutate(n = n()) %>%
    mutate(center_lon = mean(sample_lon), center_lat = mean(sample_lat))
  
  #Create a square
  
  Square <- list(lat = c(-1, -1, 1, 1), lon = c(-1, 1, -1, 1)) %>%
    as_tibble() %>%
    st_as_sf(coords = c("lon", "lat")) %>% 
    st_bbox() %>% 
    st_as_sfc()
  
  
  #Region Functions
  
  RegionRatio <- function(mu1, mu2, r){
    #Function that takes the coordinates of a midpoint and a radius to create a circle and then calculates the ratio of the circle that is contained inside each of the 4 areas or outside.
    Territory <- st_point(c(mu1, mu2)) %>%
      st_buffer(dist = r)
    
    TotalArea <- st_area(Territory)
    
    st_intersection(Territory, Square) %>%
      st_area / TotalArea %>%
      return()
  }
  
  
  #Set the radius of all bears territories to be the largest radius from a midpoint to a corresponding spill.
  
  TerritoryRadius <- samples %>%
    mutate(DistanceFromMean = (center_lat - sample_lat)^2 + (center_lon - sample_lon)^2 %>% sqrt()) %>%
    arrange(desc(DistanceFromMean)) %>%
    ungroup() %>%
    select(DistanceFromMean) %>%
    head(1) %>%
    .[[1]]
  
  
  
  #Function that Generates a maximum likelihood function that can be optimized.
  MaxLikelihoodGeneratorZeroTrunc <- function(Tibble){
    
    function(lambda){
      Tibble %>%
        mutate(Likelihood = (Inside*lambda) ^ n / ((exp(Inside*lambda) - 1) * factorial(n))) %>%
        mutate(Likelihood = -log(Likelihood)) %>%
        ungroup() %>%
        select(Likelihood) %>%
        summarise(Sum = sum(Likelihood)) %>%
        .[[1]]
    } %>%
      return()
  }
  
  
  #Calculate Mean number of samples.
  
  nMean <- samples %>%
    ungroup() %>%
    select(n) %>%
    summarise(mean(n)) %>%
    head(1) %>%
    .[[1]]
  
  #Using Optim with "BFGS"
  
  fit <- samples %>% 
    select(id, center_lon, center_lat, n) %>%
    distinct() %>%
    mutate(Inside = RegionRatio(center_lon, center_lat, TerritoryRadius)) %>%
    MaxLikelihoodGeneratorZeroTrunc() %>%
    optim(lambda, ., 
          hessian = TRUE, 
          method = "BFGS")
  
  
  Results <- list(MLE = fit$par,
                  MeanEstimate = nMean,
                  Fisher = as.numeric(sqrt(1 / fit$hessian)))
  
  return(Results)
  
}


set.seed(2022)
simulation_study <- expand_grid(lambda = 3:5, sigma = 2:3/10, mu = 4:6*100, sim = 1:3) %>% 
  rowwise() %>%
  mutate(results = list(SimulationFunction(mu, lambda, sigma))) %>%
  mutate(MLE = results[[1]], Mean = results[[2]], Fisher = results[[3]]) %>%
  select(-results)



simulation_study

write_csv(simulation_study, "data/SimulationResults1")

test <- read_csv("data/SimulationResults1")



TimingFunction <- function(lambda, sigma, mu){
  start_time <- Sys.time()
  SimulationFunction(mu, lambda, sigma)
  end_time <- Sys.time()
  return(end_time - start_time)
}


TimingFunction(3,0.1,500)
TimingFunction(3,0.3,500)
TimingFunction(14,0.1,500)
TimingFunction(14,0.3,500)

