library(tidyverse)
library(sf)


set.seed(19911001)

mu <- 500 # Mean population size in 2 x 2 square
N <- rpois(1, 9 * mu) # Total population
lambda <- 3 # Mean number of samples
sigma <- .1 # Sample spread
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
  mutate(n = n())

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
  

  #Helper function that returns the area of the intersection of the two arguments divided by the area of the territory circle.
  IntersectionRatio <- function(Region){
    st_intersection(Territory, Region) %>%
      st_area() / TotalArea
  }
  
  #Calculating the ratios.
  In <- IntersectionRatio(Square)
  Out <- 1 - In
  
  list(Inside = In, Outside = Out) %>%
    as_tibble() %>%
    round(3) %>%
    return()
}

TerritoryRadius <- samples %>%
  mutate(DistanceFromMean = (center_lat - sample_lat)^2 + (center_lon - sample_lon)^2 %>% sqrt()) %>%
  arrange(desc(DistanceFromMean)) %>%
  ungroup() %>%
  select(DistanceFromMean) %>%
  head(1) %>%
  .[[1]]


MaxLikelihoodGeneratorZeroTrunc <- function(Tibble){
  
  function(lambda){
    Tibble %>%
      mutate(Likelihood = (Inside*lambda) ^ n / ((exp(Inside*lambda)-1) * factorial(n))) %>%
      mutate(Likelihood = log(Likelihood)) %>%
      ungroup() %>%
      select(Likelihood) %>%
      map(sum) %>%
      .[[1]]
  } %>%
    return()
}


#Calculate Mean number of samples.

nMean <- samples %>%
  ungroup() %>%
  select(n) %>%
  summarise(mean(n))

#Using Optimize

samples %>% 
  select(id, center_lon, center_lat, n) %>%
  distinct() %>%
  mutate(Ratio = RegionRatio(center_lon, center_lat, TerritoryRadius)) %>%
  mutate(Inside = Ratio[[1]], Outside = Ratio[[2]]) %>%
  select(-Ratio) %>%
  MaxLikelihoodGeneratorZeroTrunc() %>%
  optimize(maximum = TRUE, lower = 0, upper = 50)


#Using Optim with "BFGS"

samples %>% 
  select(id, center_lon, center_lat, n) %>%
  distinct() %>%
  mutate(Ratio = RegionRatio(center_lon, center_lat, TerritoryRadius)) %>%
  mutate(Inside = Ratio[[1]], Outside = Ratio[[2]]) %>%
  select(-Ratio) %>%
  MaxLikelihoodGeneratorZeroTrunc() %>%
  optim(nMean, ., 
        hessian = TRUE, method = "BFGS")


TestFunction <- samples %>% 
  select(id, center_lon, center_lat, n) %>%
  distinct() %>%
  mutate(Ratio = RegionRatio(center_lon, center_lat, TerritoryRadius)) %>%
  mutate(Inside = Ratio[[1]], Outside = Ratio[[2]]) %>%
  select(-Ratio) %>%
  MaxLikelihoodGeneratorZeroTrunc()


TestFunction(3.268379)







simulation_study <- expand_grid(lambda = 1:20, sigma = 0.1, mu = 500, sim = 1:10) %>% 
  rowwise()


















SAMPLE FUNCTIONS FROM MARTIN GITHUB


rztpois <- function(n, lambda){
  # Zero-truncated Poisson random draws
  assertthat::assert_that(length(lambda) == 1, msg = "Only works for scalar lambda")
  x <- numeric()
  while (length(x) < n) {
    x <- c(x, rpois(n, lambda))
    x <- x[x > 0]
  }
  x[1:n]
}
dztpois <- function(x, lambda){
  # Zero-truncated Poisson density
  dpois(x, lambda) / (1 - dpois(0, lambda))
}
fitztpois <- function(x){
  # Zero-truncated Poisson MLE
  fit <- optim(mean(x), function(lambda) -sum(log(dztpois(x, lambda))), 
               hessian = TRUE, method = "BFGS")
  list(estimate = fit$par,
       se = as.numeric(sqrt(1 / fit$hessian))
  )
}


# Simulations


set.seed(2022) # For reproducibility
# For each of 4 values of lambda, perform 1000 simulations of n = 100 observations
sim_study <- expand_grid(lambda = 1:20, n = 100, sim = 1:1000) %>% 
  rowwise() %>% 
  mutate(data = list(tibble(x = rztpois(n, lambda))), # Simulate data
         fit = list(fitztpois(data$x)), # Fit the MLE
         results = list(tibble( # Store the results
           method = c("MLE", "mean"),
           estimate = c(fit$estimate, mean(data$x))
         )))
sim_study


# Reporting a table


sim_table <- sim_study %>% 
  select(lambda, results) %>% 
  unnest(results) %>% 
  group_by(lambda, method) %>% 
  summarise(sqrtMSE = sqrt(mean((estimate - lambda)^2)),
            se = sd(estimate),
            bias = mean(estimate - lambda),
            .groups = "drop")
sim_table


# Reporting a graph


sim_table %>%  
  pivot_longer(-(lambda:method), names_to = "metric") %>% 
  ggplot(aes(x = lambda, y = value, color = method)) + geom_line()  + geom_point() + 
  facet_wrap(~metric)

