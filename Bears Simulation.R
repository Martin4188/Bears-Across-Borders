library(tidyverse)
library(sf)
library(polyCub)
library(VGAM)


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
  
  #Estimate sigma using pooled variance.
  sigmaEstimate <- samples %>%
    filter(n != 1) %>%
    group_by(id) %>%
    mutate(center_lon = mean(sample_lon),
           center_lat = mean(sample_lat)) %>%
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
                  NObsTrue = NObsTrue,
                  SigmaHat = sigmaEstimate) %>%
    as_tibble()
  
  return(Results)
  
}

set.seed(2022)
WithSigma1 <- expand_grid(mu = 500, lambda = 2:4, sigma = 1:3/10, sim = 1:1000) %>% 
  rowwise() %>%
  mutate(results = NaiveSimulationFunction(mu, lambda, sigma)) %>%
  mutate(MLE = results[[1]], 
         Mean = results[[2]], 
         Fisher = results[[3]], 
         NTrue = results[[4]], 
         NObs = results[[5]], 
         NObsTrue = results[[6]],
         SigmaHat = results[[7]]) %>%
  select(-results)


set.seed(2022)
WithSigma2 <- expand_grid(mu = 500, lambda = 2:4, sigma = 1:3/40, sim = 1:1000) %>% 
  rowwise() %>%
  mutate(results = NaiveSimulationFunction(mu, lambda, sigma)) %>%
  mutate(MLE = results[[1]], 
         Mean = results[[2]], 
         Fisher = results[[3]], 
         NTrue = results[[4]], 
         NObs = results[[5]], 
         NObsTrue = results[[6]],
         SigmaHat = results[[7]]) %>%
  select(-results)

set.seed(2022)
WithSigma3 <- expand_grid(mu = c(400, 650, 1000), lambda = 3, sigma = 1:4/40, sim = 1:1000) %>% 
  rowwise() %>%
  mutate(results = NaiveSimulationFunction(mu, lambda, sigma)) %>%
  mutate(MLE = results[[1]], 
         Mean = results[[2]], 
         Fisher = results[[3]], 
         NTrue = results[[4]], 
         NObs = results[[5]], 
         NObsTrue = results[[6]],
         SigmaHat = results[[7]]) %>%
  select(-results)


WithSigma1 %>%
  write_csv("data/NaiveSimulation")

WithSigma2 %>%
  write_csv("data/NaiveSimulationSmallSigma")

WithSigma3 %>%
  write_csv("data/NaiveSimulationMu")




###################################
###################################
###################################


###IMPROVEMENTS###



NaiveSimulationFunction2 <- function(mu, lambda, sigma){
  
  #Simulations
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
    mutate(n = n(),
           center_lon = mean(sample_lon),
           center_lat = mean(sample_lat))
  
  #Function that Generates a maximum likelihood function that can be optimized.
  #WITHOUT IRATIO
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
  
  
  #Function that Generates a maximum likelihood function that can be optimized.
  #Using IRATIO
  MaxLikelihoodGeneratorZeroTrunc2 <- function(Tibble){
    
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
  
  
  InsideRatioFunction <- function(mu1, mu2, sigma, Precision = 5){
    
    
    
    
    #Create a square corresponding to the area we are investigating.
    square <- list(lon = c(1, 1, -1, -1), 
                   lat = c(1, -1, -1, 1)) %>%
      as_tibble() %>%
      st_as_sf(coords = c("lon", "lat")) %>%
      st_bbox() %>%
      st_as_sfc() - c(mu1, mu2) 
    
    #polyCub requires a zero mean normal distribution, as such i translate the 
    # middlepoint so that it becomes (0, 0) and perform the same translation on the Region.
    
    Binormal <- function (s, S = sigma)
    {
      exp(-rowSums(s^2)/2/S^2) / (2*pi*S^2)
    }
    
    Intersection <- st_point(c(0, 0)) %>% 
      st_buffer(dist = 3 * sigma) %>%
      st_intersection(square)
    
    
    
    if(st_is_empty(Intersection)){
      return(0)
    }
    
    
    Intersection %>%
      polyCub.SV(f = Binormal, nGQ = Precision) %>%
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
  
  #Estimate sigma using pooled variance.
  sigmaEstimate <- samples %>%
    filter(n != 1) %>%
    group_by(id) %>%
    mutate(center_lon = mean(sample_lon),
           center_lat = mean(sample_lat)) %>%
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
  
  
  samples <- samples %>%
    mutate(IRatio = InsideRatioFunction(center_lon, center_lat, sigmaEstimate),
           IRatio = if_else(IRatio > 0.985, 1, IRatio))
  
  
  FractionalBears <- samples %>%
    ungroup() %>%
    select(id, IRatio) %>%
    distinct() %>%
    summarise(IRatio = sum(IRatio)) %>%
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
  
  #Alternate Optim
  
  fit2 <- samples %>% 
    select(id, n, IRatio) %>%
    distinct() %>%
    MaxLikelihoodGeneratorZeroTrunc2() %>%
    optim(lambda, ., 
          hessian = TRUE, 
          method = "L-BFGS-B",
          lower = 0.1)
  
  
  Results <- list(MLE = fit$par,
                  MeanEstimate = nMean,
                  Fisher = as.numeric(sqrt(1 / fit$hessian)),
                  NTrue = NTrue,
                  NObs = NObs,
                  NObsTrue = NObsTrue,
                  SigmaHat = sigmaEstimate,
                  AltMLE = fit2$par,
                  RatioSum = FractionalBears) %>%
    as_tibble()
  
  return(Results)
  
}


set.seed(2022)
RatioSimulation1 <- expand_grid(mu = 500, lambda = 2:4, sigma = 1:4/40, sim = 1:1000) %>% 
  rowwise() %>%
  mutate(results = NaiveSimulationFunction2(mu, lambda, sigma)) %>%
  mutate(MLE = results[[1]], 
         Mean = results[[2]], 
         Fisher = results[[3]], 
         NTrue = results[[4]], 
         NObs = results[[5]], 
         NObsTrue = results[[6]],
         SigmaHat = results[[7]],
         AltMLE = results[[8]],
         RatioSum = results[[9]]) %>%
  select(-results)


RatioSimulation1 %>%
  write_csv("data/RatioSimulation")


