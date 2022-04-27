library(tidyverse)
library(sf)
library(knitr)
library(polyCub)
library(VGAM)
options(dplyr.summarise.inform=F)


NaiveSimulation <- read_csv("data/NaiveSimulation")
NaiveSimulationSmallSigma <- read_csv("data/NaiveSimulationSmallSigma")
SwedishSurveyResults <- read.csv("data/SwedishSurveyEstimate")
SwedishSurveyResults2 <- read.csv("data/SwedishSurveyEstimate2")



captures <- read_csv("data/captures.csv") %>%
  group_by(id, year) %>%
  rename(sample_lon = lon,
         sample_lat = lat) %>%
  mutate(center_lon = mean(sample_lon),
         center_lat = mean(sample_lat),
         n = n()) %>%
  ungroup()





NaiveEstimate <- function(YEAR){

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
    filter(year == YEAR) %>%
    select(id, n) %>%
    distinct() %>%
    summarise(n=n()) %>%
    .[[1]]
    
  
  LambdaEstimate <- captures %>% 
    filter(year == YEAR) %>%
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

#StandardEstimate <- expand.grid(year = c("2015", "2016", "2017", "2019", "2020")) %>%
#  rowwise() %>%
#  mutate(Estimate = NaiveEstimate(year))







###FIGURES###







#Figure 1: Bias and standard error of the population estimate.
Figure1 <- NaiveSimulation %>%
  full_join(NaiveSimulationSmallSigma, 
            by = c("...1", "mu", "lambda", "sigma", "sim", "MLE", "Mean", "Fisher", "NTrue", "NObs", "NObsTrue", "SigmaHat")) %>%
  select(lambda, sigma, MLE, NObs, NTrue) %>%
  mutate(PopEstimate = NObs / (1 - exp(-MLE)),
         bias = PopEstimate - NTrue,
         relativeBias = bias / NTrue) %>%
  select(lambda, sigma, PopEstimate, bias, relativeBias) %>%
  group_by(lambda, sigma) %>%
  summarise(SE = sd(PopEstimate),
            across(where(is.numeric), mean)) %>%
  select(-PopEstimate) %>%
  rename('Standard Error' = SE,
         'Bias' = bias,
         'Bias relative to true population' = relativeBias) %>%
  pivot_longer('Standard Error':'Bias relative to true population', names_to = "metric")  %>%
  mutate(lambda = lambda %>% as.factor()) %>%
  ggplot(aes(x = sigma, y = value, color = lambda)) + geom_line()  + geom_point() + 
  facet_wrap(~metric, scales = "free_y") +
  labs(title = "Figure 1: Bias and standard error of the population estimate.", x = "Sigma", y = "")



#Figure 2: Mean number of false bears
Figure2 <- NaiveSimulation %>%
  full_join(NaiveSimulationSmallSigma, 
            by = c("...1", "mu", "lambda", "sigma", "sim", "MLE", "Mean", "Fisher", "NTrue", "NObs", "NObsTrue", "SigmaHat")) %>%
  mutate(NObsFalse = NObs - NObsTrue,
         factor = 1 / (1 - exp(- MLE)),
         Relative = NObsFalse / NObs) %>%
  group_by(lambda, sigma) %>%
  summarise(SE = sd(NObsFalse),
            across(where(is.numeric), mean)) %>%
  mutate(lambda = as.factor(lambda)) %>%
  select(lambda, sigma, NObsFalse, Relative, SE) %>%
  rename('Mean number of false bears' = NObsFalse,
         'Ratio of observed population' = Relative,
         'Standard Error' = SE) %>%
  pivot_longer('Mean number of false bears':'Standard Error', names_to = "metric") %>%
  ggplot(aes(x = sigma, y = value, color = lambda)) + geom_line()  + geom_point() + 
  facet_wrap(~metric, scales = "free_y") +
  labs(title = "Figure 2: Mean number of observed false bears", x = "Sigma", y = "")


#Figure 3 Sigma bias and standard error

Figure3 <- NaiveSimulation %>%
  full_join(NaiveSimulationSmallSigma, 
            by = c("...1", "mu", "lambda", "sigma", "sim", "MLE", "Mean", "Fisher", "NTrue", "NObs", "NObsTrue", "SigmaHat")) %>%
  select(lambda, sigma, SigmaHat) %>%
  mutate(bias = SigmaHat - sigma,
         relativeBias = bias / sigma) %>%
  group_by(lambda, sigma) %>%
  summarise(SE = sd(SigmaHat),
            across(where(is.numeric), mean)) %>%
  mutate(lambda = lambda %>% as.factor()) %>%
  rename('Bias' = bias,
         'Relative Bias' = relativeBias,
         'Standard Error' = SE) %>%
  pivot_longer(c('Standard Error', 'Bias', 'Relative Bias'), names_to = "metric") %>%
  ggplot(aes(x = sigma, y = value, color = lambda)) + geom_line()  + geom_point() + 
  facet_wrap(~metric, scales = "free_y") +
  labs(title = "Figure 3: Bias and standard error of estimate of sigma.", x = "Sigma", y = "")

#Figure 4: MLE bias and standard error
Figure4 <- NaiveSimulation %>%
  full_join(NaiveSimulationSmallSigma, 
            by = c("...1", "mu", "lambda", "sigma", "sim", "MLE", "Mean", "Fisher", "NTrue", "NObs", "NObsTrue", "SigmaHat")) %>%
  select(lambda, sigma, MLE) %>%
  mutate(bias = MLE - lambda,
         SQRTMSE = (MLE - lambda)^2,
         relativeBias = bias / lambda) %>%
  group_by(lambda, sigma) %>%
  summarise(SE = sd(MLE),
            across(where(is.numeric), mean)) %>%
  mutate(lambda = lambda %>% as.factor(),
         SQRTMSE = SQRTMSE %>% sqrt()) %>%
  rename('Bias' = bias,
         'Relative Bias' = relativeBias,
         'Standard Error' = SE) %>%
  pivot_longer(c('Standard Error', 'Bias', 'Relative Bias'), names_to = "metric") %>%
  ggplot(aes(x = sigma, y = value, color = lambda)) + geom_line()  + geom_point() + 
  facet_wrap(~metric, scales = "free_y") +
  labs(title = "Figure 4: Bias and standard error of the maximum likelihood estimate of lambda.", x = "Sigma", y = "")



#Figure 5: Bias in the multiplication factor.

Figure5 <- NaiveSimulation %>%
  full_join(NaiveSimulationSmallSigma, 
            by = c("...1", "mu", "lambda", "sigma", "sim", "MLE", "Mean", "Fisher", "NTrue", "NObs", "NObsTrue", "SigmaHat")) %>%
  select(lambda, sigma, MLE, NTrue, NObs, NObsTrue) %>%
  mutate(NObsFalse = NObs - NObsTrue,
         factor = 1 / (1 - exp(-lambda)),
         factorHat = 1 / (1 - exp(- MLE)),
         bias = factorHat - factor,
         RelativeBias = bias/factor) %>%
  group_by(lambda, sigma) %>%
  summarise(SE = sd(factorHat),
            across(where(is.numeric), mean)) %>%
  select(lambda, sigma, bias, RelativeBias, SE) %>%
  rename('Bias' = bias,
         'Relative bias' = RelativeBias,
         'Standard Error' = SE) %>%
  pivot_longer('Bias':'Standard Error', names_to = "metric") %>%
  mutate(lambda = lambda %>% as.factor()) %>%
  ggplot(aes(x = sigma, y = value, color = lambda)) + geom_line()  + geom_point() + 
  facet_wrap(~metric, scales = "free_y") +
  labs(title = "Figure 5: Bias and standard error for the multiplication factor.", x = "Sigma", y = "")


#Visualizing sigma



#FIGURE SOMETHING

NaiveSimulation %>%
  full_join(NaiveSimulationSmallSigma, 
            by = c("...1", "mu", "lambda", "sigma", "sim", "MLE", "Mean", "Fisher", "NTrue", "NObs", "NObsTrue", "SigmaHat")) %>%
  select(lambda, sigma, MLE, NTrue, NObs, NObsTrue) %>%
  mutate(NObsFalse = NObs - NObsTrue,
         FalseBearFactor = (NObsTrue + NObsFalse)/NObsTrue,
         LambdaBias = NObsTrue / (1 - exp(-MLE)),
         LambdaBias2 = (NObsTrue / (1 - exp(-MLE)) / NObs)) %>%
  group_by(lambda,sigma) %>%
  summarise(across(where(is.numeric), mean))


#Error in estimate if we only observed true bears.

FigureX <- NaiveSimulation %>%
  full_join(NaiveSimulationSmallSigma, 
            by = c("...1", "mu", "lambda", "sigma", "sim", "MLE", "Mean", "Fisher", "NTrue", "NObs", "NObsTrue", "SigmaHat")) %>%
  select(lambda, sigma, MLE, NTrue, NObs, NObsTrue) %>%
  mutate(TrueEstimate = NObsTrue / (1- exp(-lambda)),
         LambdaError = NObsTrue / (1 - exp(-MLE))) %>%
  group_by(lambda,sigma) %>%
  summarise(across(where(is.numeric), mean)) %>%
  mutate(bias = LambdaError - NTrue,
         relativeBias = bias/NTrue) %>%
  rename('Bias' = bias,
         'Relative bias' = relativeBias) %>%
  pivot_longer('Bias':'Relative bias', names_to = "metric") %>%
  mutate(lambda = lambda %>% as.factor()) %>%
  ggplot(aes(x = sigma, y = value, color = lambda)) + geom_line()  + geom_point() + 
  facet_wrap(~metric, scales = "free_y") +
  labs(title = "Figure X: Bias in the case only true bears were observed.", x = "Sigma", y = "")
  


#Function Graph



base <-
  ggplot() +
  ylim(0, 3) +
  xlim(0, 6)

FunctionGraph <- base + geom_function(fun = ~1 / (1- exp(-.x)))







##########################


#SwedishSurveyResults2 <- SwedishSurveyResults2 %>%
#  mutate(year = as.factor(year)) %>%
#  mutate(PopulationEstimate = NObserved / (1 - exp(-LambdaHat)))

#SwedishResultsTable <- StandardEstimate %>%
#  unpack(Estimate) %>%
#  full_join(SwedishSurveyResults2, by = "year") %>%
#  select(year, Estimate, PopulationEstimate, LambdaHat.x, LambdaHat.y, NormalizedSigma) %>%
#  rename('Standard Estimate' = Estimate,
#         'Alt Estimate' = PopulationEstimate,
#         'Standard Lambda' = LambdaHat.x,
#         'Alt Lambda' = LambdaHat.y,
#         'Normalized Sigma' = NormalizedSigma)
  

#SwedishSurveyResults2 %>%
#  mutate(year = as.factor(year))





sigmaVisualizer <- function(SIGMA){
  square <- list(
    list(x = c(1, 1, -1, -1),
         y = c(1, -1, -1, 1))
  )
  
  f <- function (s, sigma = SIGMA)
  {
    exp(-rowSums(s^2)/2/sigma^2) / (2*pi*sigma^2)
  }
  
  polyCub.SV(square, f, nGQ = 1, plot = TRUE)
  
  
}





