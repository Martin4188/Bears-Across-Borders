library(tidyverse)
library(sf)
library(knitr)
library(polyCub)
library(VGAM)
options(dplyr.summarise.inform=F)


#NaiveSimulation <- read_csv("data/NaiveSimulation") %>%
#  filter(sigma == 0.1)
#NaiveSimulationSmallSigma <- read_csv("data/NaiveSimulationSmallSigma")
NaiveSimulationMu <- read_csv("data/NaiveSimulationMu")
SwedishSurveyResults <- read.csv("data/SwedishSurveyEstimate")
SwedishSurveyResults2 <- read.csv("data/SwedishSurveyEstimate2")
StandardEstimate <- read_csv("data/StandardEstimate") %>%
  mutate(year = as.factor(year))
RatioSimulation <- read_csv("data/RatioSimulation")



###Creating the Swedish Results Figure###
SwedishSurveyResults <- SwedishSurveyResults %>%
  mutate(year = as.factor(year)) %>%
  mutate(PopulationEstimate = NObserved / (1 - exp(-LambdaHat)),
         PopulationEstimate2 = RatioSum / (1 - exp(-LambdaHat)))

SwedishResults <- StandardEstimate %>%
  full_join(SwedishSurveyResults, by = c("year", "sex")) %>%
  select(year, sex, Estimate, PopulationEstimate, PopulationEstimate2, LambdaHat.x, LambdaHat.y, NormalizedSigma) %>%
  rename('Standard Estimate' = Estimate,
         'Alt Estimate' = PopulationEstimate,
         'Ratio Estimate' = PopulationEstimate2,
         'Standard Lambda' = LambdaHat.x,
         'Alt Lambda' = LambdaHat.y,
         'Normalized Sigma' = NormalizedSigma) %>%
  arrange(year)

SwedishResultsTable <- SwedishResults %>%
  select(sex, year, `Standard Estimate`, `Alt Estimate`, `Ratio Estimate`) %>%
  group_by(year) %>%
  summarise(across(where(is.numeric), sum)) %>%
  mutate(sex = "Total") %>%
  full_join(SwedishResults, by = c("year", "sex", "Standard Estimate", "Alt Estimate", "Ratio Estimate")) %>%
  select(year, sex, `Standard Estimate`, `Alt Estimate`, `Ratio Estimate`, `Standard Lambda`, `Alt Lambda`, `Normalized Sigma`) %>%
  arrange(year, sex)


BiasList <- list(year = c("2015", "2015", "2016", "2016", "2017", "2017", "2019", "2019", "2020", "2020"),
                 sex = c("Hane", "Hona", "Hane", "Hona", "Hane", "Hona", "Hane", "Hona","Hane", "Hona"),
                 `Simulated Bias Estimate` = c(0.17, 0.12, 0.12, 0.07, 0.13, 0.09, 0.15, 0.10, 0.13, 0.09)) %>%
  as_tibble()


EstimatedBiasTibble <- SwedishResultsTable %>%
  select(year, sex, `Standard Estimate`,`Alt Estimate` , `Ratio Estimate`) %>%
  filter(sex != "Total") %>%
  full_join(BiasList, by = c("year", "sex")) %>%
  mutate(`Simulated Unbiased Estimate` = `Standard Estimate` / (1 + `Simulated Bias Estimate`),
         `Estimated Simulated Bias` = `Standard Estimate` - `Simulated Unbiased Estimate`)







SwedishResultsTable %>%
  select(year, sex, `Standard Estimate`, `Alt Estimate`, `Ratio Estimate`) %>%
  pivot_longer(`Standard Estimate`:`Ratio Estimate`, names_to = "metric") %>%
  ggplot(aes(y = value, x = sex, fill = factor(reorder(metric,-value)))) +
  geom_col(position = "dodge2") +
  facet_wrap(~year)

EstimatedBiasTable <- EstimatedBiasTibble %>%
  select(-`Simulated Bias Estimate`) %>%
  group_by(year) %>%
  summarise(across(where(is.numeric), sum)) %>%
  mutate(sex = "Total") %>%
  full_join(EstimatedBiasTibble) %>%
  arrange(year)


#####################################################################



###FIGURES###


#Figure 1: Bias and standard error of the population estimate.
PopulationEstimateFigure <- RatioSimulation %>%
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
  pivot_longer(c('Bias relative to true population'), names_to = "metric")  %>%
  mutate(lambda = lambda %>% as.factor()) %>%
  ggplot(aes(x = sigma, y = value, color = lambda)) + geom_line()  + geom_point() + 
  facet_wrap(~metric, scales = "free_y") +
  labs(title = "Figure 1: Mean bias of the simulated population estimate.", x = "Sigma", y = "") +
  xlim(0, 0.1) +
  ylim(0, 0.25)

#



############################################################################################################################
###Creating Figure 2###

TestPlot1 <- RatioSimulation %>%
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
  select(lambda, sigma, metric, value) %>%
  mutate(lambda = str_c("Lambda = " ,lambda, ", False Bear"))


TestPlot2 <- RatioSimulation%>%
  mutate(NObsFalse = NObs - NObsTrue,
         factor = 1 / (1 - exp(- MLE)),
         Relative = NObsFalse / NObs) %>%
  group_by(lambda, sigma) %>%
  summarise(SE = sd(NObsFalse),
            across(where(is.numeric), mean)) %>%
  mutate(lambda = as.factor(lambda)) %>%
  mutate(NObsFalse = NObsFalse * factor,
         Relative = Relative * factor) %>%
  select(lambda, sigma, NObsFalse, Relative, SE) %>%
  rename('Bias' = NObsFalse,
         'Relative bias' = Relative) %>%
  pivot_longer('Bias':'Relative bias', names_to = "metric") %>%
  select(lambda, sigma, metric, value) %>%
  mutate(lambda = str_c("Lambda = " ,lambda, ", Lambda Error"))

TestPlot3 <- full_join(TestPlot1, TestPlot2)


TrueFalseBiasFigure <- TestPlot3 %>%
  ggplot(aes(x = sigma, y = value, color = lambda)) + geom_line()  + geom_point() + 
  facet_wrap(~metric, scales = "free_y") +
  labs(title = "Figure 2: Bias divided into False bear and Lambda bias.", x = "Sigma", y = "")

##########################################################################################################################



#Figure3

ParameterBiasFigure <- RatioSimulation%>%
  mutate('Lambda' = MLE - lambda,
         'Sigma' = SigmaHat - sigma) %>%
  group_by(lambda, sigma) %>%
  summarise(across(where(is.numeric), mean)) %>%
  pivot_longer(c('Lambda', 'Sigma'), names_to = "metric") %>%
  ggplot(aes(x = sigma, y = value, color = lambda %>% as.factor())) + geom_line()  + geom_point() + 
  facet_wrap(~metric, scales = "free_y") +
  labs(title = "Figure 3: Bias in the estimation of the parameters.", x = "Sigma", y = "", color = "Lambda")


MultiplicationFactorFigure <- RatioSimulation %>%
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
  pivot_longer(c('Bias', 'Relative bias'), names_to = "metric") %>%
  mutate(lambda = lambda %>% as.factor()) %>%
  ggplot(aes(x = sigma, y = value, color = lambda)) + geom_line()  + geom_point() + 
  facet_wrap(~metric, scales = "free_y") +
  labs(title = "Figure 5: Bias for the multiplication factor.", x = "Sigma", y = "")





#Figure 4: Difference in estimate of Swedish population between the different methods.

SwedishBearPopulationEstimateFigure1 <- EstimatedBiasTable  %>%
  pivot_longer(c(`Standard Estimate`,`Alt Estimate` ,`Ratio Estimate`, `Simulated Unbiased Estimate`), names_to = "metric") %>%
  filter(metric != "Simulated Unbiased Estimate",
         year == 2016 | year == 2017 | year == 2019) %>% #Removing the simulated unbiased estimate
  ggplot(aes(y = value, x = sex, fill = factor(reorder(metric,-value)))) +
  geom_col(position = "dodge2") +
  facet_wrap(~year) +
  labs(title = "Figure 6: Estimation of the Swedish bear population in 3 different regions.", x = "", y = "", fill = "Method")


SwedishBearPopulationEstimateFigure2 <- EstimatedBiasTable  %>%
  pivot_longer(c(`Standard Estimate`,`Alt Estimate` ,`Ratio Estimate`, `Simulated Unbiased Estimate`), names_to = "metric") %>%
  filter(metric != "Simulated Unbiased Estimate",
         year == 2015 | year == 2020) %>% 
  ggplot(aes(y = value, x = sex, fill = factor(reorder(metric,-value)))) +
  geom_col(position = "dodge2") +
  facet_wrap(~year) +
  labs(title = "Figure 7: Estimation of the Swedish bear population in the same region two different years.", x = "", y = "", fill = "Method")










  



#########################################

###Relative bias does not change with the total population.

TotalPopulationRelativeBiasFigure <- NaiveSimulationMu %>%
  mutate(PopEstimate = NObs / (1 - exp(-MLE)),
         bias = PopEstimate - NTrue,
         relativeBias = bias / NTrue) %>%
  group_by(mu, lambda, sigma) %>%
  summarise(across(where(is.numeric), mean)) %>%
  select(mu, lambda, sigma, bias, relativeBias) %>%
  ggplot(aes(x = sigma, y = relativeBias, color = as_factor(mu))) + geom_line()  + geom_point() + 
  labs(title = "Figure 1: Relative bias does not change with the total population.", x = "Sigma", y = "Relative Bias")




base <-
  ggplot() +
  ylim(0, 3) +
  xlim(0, 6)

FunctionGraph <- base + geom_function(fun = ~1 / (1- exp(-.x))) +
  labs(title = "Function Graph 1: The function graph for the transformed lambda.", x = "Lambda", y = "")







##########################


#Total <- SwedishSurveyResults %>%
#  group_by(year) %>%
#  select(year, sex, NObserved, WithoutOutsiders, RatioSum) %>%
#  summarise(across(where(is.numeric), sum)) %>%
#  mutate(sex = "Total")


#SwedishSurveyResults %>%
#  select(year, sex, NObserved, WithoutOutsiders, RatioSum) %>%
#  full_join(Total, by = c("year", "sex", "NObserved", "WithoutOutsiders", "RatioSum"))


#####
#SwedishSurveyResults2 <- SwedishSurveyResults2 %>%
#  mutate(year = as.factor(year)) %>%
#  mutate(PopulationEstimate = NObserved / (1 - exp(-LambdaHat)),
#         PopulationEstimate2 = RatioSum / (1 - exp(-LambdaHat)))
#
#SwedishResultsTable <- StandardEstimate %>%
#  full_join(SwedishSurveyResults2, by = c("year", "sex")) %>%
#  select(year, Estimate, PopulationEstimate, PopulationEstimate2, LambdaHat.x, LambdaHat.y, NormalizedSigma) %>%
#  rename('Standard Estimate' = Estimate,
#         'Alt Estimate' = PopulationEstimate,
#         'Ratio Estimate' = PopulationEstimate2,
#         'Standard Lambda' = LambdaHat.x,
#         'Alt Lambda' = LambdaHat.y,
#         'Normalized Sigma' = NormalizedSigma)
  

#SwedishSurveyResults2 %>%
#  mutate(year = as.factor(year))

###################



  

#########################################################

#Visualizing sigma

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

#Comparing bias between standard and ratio estimate.

RatioEstimateBiasFigure <- RatioSimulation %>%
  select(lambda, sigma, NTrue, NObs, RatioSum, MLE, AltMLE) %>%
  mutate(Estimate = NObs / (1 - exp(-MLE)),
         altEstimate = RatioSum / (1 - exp(-AltMLE)),
         bias = (Estimate - NTrue) / NTrue,
         AltBias = (altEstimate - NTrue) / NTrue,
         lambda = lambda %>% as_factor()) %>%
  group_by(lambda, sigma) %>%
  summarise(across(where(is.numeric), mean)) %>%
  pivot_longer(c(bias, AltBias), names_to = "metric") %>%
  ggplot(aes(x = sigma, y = value, color = metric)) + geom_line()  + geom_point() + 
  facet_wrap(~lambda) +
  labs(title = "Figure 4: Comparing bias in the population estimate between the two methods.", x = "Sigma", y = "", color = "Method")
  

#Comparing bias in lambda between standard and ratio estimate.

RatioLambdaEstimateBiasFigure <- RatioSimulation %>%
  select(lambda, sigma, NTrue, NObs, RatioSum, MLE, AltMLE) %>%
  mutate(Estimate = NObs / (1 - exp(-MLE)),
         altEstimate = RatioSum / (1 - exp(-AltMLE)),
         lambdaBias = MLE - lambda,
         altLambdaBias = AltMLE - lambda,
         lambda = lambda %>% as_factor()) %>%
  group_by(lambda, sigma) %>%
  summarise(across(where(is.numeric), mean)) %>%
  pivot_longer(c(lambdaBias, altLambdaBias), names_to = "metric") %>%
  ggplot(aes(x = sigma, y = value, color = metric)) + geom_line()  + geom_point() + 
  facet_wrap(~lambda) +
  labs(title = "Figure 5: Comparing bias in the estimate of lambda between the two methods.", x = "Sigma", y = "", color = "Method")






