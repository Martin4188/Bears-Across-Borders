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


WithSigma1 %>%
  write.csv("data/NaiveSimulation")

WithSigma2 %>%
  write.csv("data/NaiveSimulationSmallSigma")










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


simulation_study_small_sigma %>%
  write.csv("data/NaiveSimulationSmallSigma")



set.seed(2022)
simulation_study_400_mu <- expand_grid(mu = 400, lambda = 2:4, sigma = 1:3/10, sim = 1:1000) %>% 
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
simulation_study_600_mu <- expand_grid(mu = 600, lambda = 2:4, sigma = 1:3/10, sim = 1:1000) %>% 
  rowwise() %>%
  mutate(results = NaiveSimulationFunction(mu, lambda, sigma)) %>%
  mutate(MLE = results[[1]], 
         Mean = results[[2]], 
         Fisher = results[[3]], 
         NTrue = results[[4]], 
         NObs = results[[5]], 
         NObsTrue = results[[6]]) %>%
  select(-results)









```{r}
library(tidyverse)
library(sf)
library(polyCub)

Square <- list(lon = c(-1, -1, 1, 1), lat = c(-1, 1, -1, 1)) %>%
  as_tibble() %>%
  st_as_sf(coords = c("lon", "lat")) %>%
  st_bbox() %>%
  st_as_sfc()



Square <- list(
  list(x = c(1, 1, -1, -1),
       y = c(1, -1, -1, 1))
)

Square

f <- function (s, sigma = 0.1)
{
  exp(-rowSums(s^2)/2/sigma^2) / (2*pi*sigma^2)
}

hexagon <- st_point(c(0, 0)) %>%
  st_buffer(dist = 20)


#plotpolyf(Square, f, xlim = c(-8,8), ylim = c(-8,8))

polyCub.SV(Square, f, nGQ = 20, plot = TRUE)

```





Square <- list(lat = c(-1, -1, 1, 1), lon = c(-1, 1, -1, 1)) %>%
  as_tibble() %>%
  st_as_sf(coords = c("lon", "lat")) %>% 
  st_bbox() %>% 
  st_as_sfc()


TestCircle <- st_point(c(1, 1)) %>%
  st_buffer(dist = 1)



  #Function that takes the coordinates of a midpoint and a radius to create a circle and then calculates the ratio of the circle that is contained inside each of the 4 areas or outside.
  Territory <- st_point(c(mu1, mu2)) %>%
    st_buffer(dist = r)
  
  TotalArea <- st_area(Territory)
  
  st_intersection(Territory, Square) %>%
    st_area / TotalArea

checkintrfr(intrfr = TestCircle , f =  BiNormal)




























#HEEEEEEEEEEEEEEEEHEHEHEHEHEHHEHE







SimulationFunction2 <- function(mu, lambda, sigma){
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
    mutate(n = n()) %>%
    mutate(center_lon = mean(sample_lon),
           center_lat = mean(sample_lat)) %>%
    ungroup()
  
  
  
  NTrue <- bears %>%
    select(Inside) %>%
    summarise(Inside = sum(Inside)) %>%
    .[[1]]
  
  
  #Calculate Observed number of bears
  
  NObs <- samples %>%
    select(id) %>%
    distinct() %>%
    summarise(n=n()) %>%
    .[[1]]
  
  #Calculate the number of observed bears whose true midpoint lies inside the middle square.
  
  NObsTrue <- samples %>%
    select(id, Inside) %>%
    distinct() %>%
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
    
    st_point(c(0, 0)) %>% 
      st_buffer(dist = 3 * sigma) %>%
      st_intersection(square) %>%
      polyCub.SV(f = Binormal, nGQ = Precision) %>%
      return()
  }
  
  
  
  
  
  #Adjusting Center Points
  
  center_fit <- function(data){
    
    n <- data$n[1]
    CenterLon <- data$center_lon[1]
    CenterLat <- data$center_lat[1]
    
    if(n == 1){
      return(c(CenterLon, CenterLat))
    }
    
    if(abs(CenterLon) < 1 - 3 * sigmaEstimate & abs(CenterLat) < 1 - 3 * sigmaEstimate){
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
        mutate(likelihood = likelihood / InsideRatioFunction(MidPoint[1], MidPoint[2], sigmaEstimate)) %>%
        select(likelihood) %>%
        .[[1]]
    }
    
    fit <-  optim(par = c(CenterLon, CenterLat), 
                  fn = LikelihoodFunction, 
                  method = "L-BFGS-B", 
                  lower = c(-(1+2*sigmaEstimate), -(1+2*sigmaEstimate)),
                  upper = c( (1+2*sigmaEstimate),  (1+2*sigmaEstimate)))
    
    fit$par %>%
      return()
  }
  
  
  
  
  
  AdjustedCenter <- samples %>% 
    nest_by(id) %>%
    mutate(center_point = list(center_fit(data))) %>%
    select(id, center_point) %>%
    mutate(center_lon = center_point[1],
           center_lat = center_point[2]) %>%
    select(-center_point)
  
  
  OptimTibble <- samples %>%
    select(-center_lon, -center_lat, -sample_lon, -sample_lat) %>%
    full_join(AdjustedCenter, by = "id") %>%
    distinct() %>%
    rowwise() %>%
    mutate(IRatio = InsideRatioFunction(center_lon, center_lat, sigmaEstimate))
  
  fit <- MaxLikelihoodGeneratorZeroTrunc(OptimTibble) %>%
    optim(lambda, ., 
          hessian = TRUE, 
          method = "L-BFGS-B",
          lower = 0.1)
  
  
  #Estimating population through fractions
  
  FractionalBears <- OptimTibble %>%
    ungroup() %>%
    summarise(Iratio = sum(IRatio)) %>%
    .[[1]]
  
  #Estimating population through removal of outside bears
  
  NonOutisdeBears <- OptimTibble %>%
    filter(abs(center_lon) < 1, abs(center_lat) < 1) %>%
    ungroup() %>%
    summarise(n=n()) %>%
    .[[1]]
  
  
  Results <- list(MLE = fit$par,
                  Sigma = sigmaEstimate,
                  FractionEstimate = FractionalBears,
                  OnlyTrueBears = NonOutisdeBears,
                  NTrue = NTrue,
                  NObs = NObs,
                  NObsTrue = NObsTrue) %>%
    as_tibble()
  
  return(Results)
}

set.seed(2022)
simulation_study_Alternative <- expand_grid(mu = 500, lambda = 4, sigma = 1:3/10, sim = 1:10) %>% 
  rowwise() %>%
  mutate(results = SimulationFunction2(mu, lambda, sigma)) %>%
  mutate(MLE = results[[1]], 
         SigmaEstimate = results[[2]], 
         FractionEstimate = results[[3]], 
         OnlyTrueBears = results[[4]],
         NTrue = results[[5]], 
         NObs = results[[6]], 
         NObsTrue = results[[7]]) %>%
  select(-results)
































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

