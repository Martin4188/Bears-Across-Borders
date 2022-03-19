Bears Across Borders
================
Martin Andersson
2022-01-29

``` r
NaiveSimulation <- read_csv("data/NaiveSimulation1")
```

    ## New names:
    ## * `` -> ...1

    ## Rows: 9000 Columns: 8

    ## -- Column specification --------------------------------------------------------
    ## Delimiter: ","
    ## dbl (8): ...1, mu, lambda, sigma, sim, MLE, Mean, Fisher

    ## 
    ## i Use `spec()` to retrieve the full column specification for this data.
    ## i Specify the column types or set `show_col_types = FALSE` to quiet this message.

``` r
NaiveSimulation %>%
  group_by(lambda, sigma) %>%
  mutate(MLEBias = abs(lambda - MLE), MeanBias = abs(lambda - Mean)) %>%
  summarise(MLE = mean(MLE), Mean = mean(Mean), Fisher = mean(Fisher), MLEBias = mean(MLEBias), MeanBias = mean(MeanBias)) %>%
  kable()
```

    ## `summarise()` has grouped output by 'lambda'. You can override using the `.groups` argument.

| lambda | sigma |      MLE |     Mean |    Fisher |   MLEBias |  MeanBias |
|-------:|------:|---------:|---------:|----------:|----------:|----------:|
|      3 |   0.1 | 2.522731 | 3.668154 | 0.0739489 | 0.4772685 | 0.6681542 |
|      3 |   0.2 | 2.147131 | 3.353964 | 0.0655808 | 0.8528692 | 0.3539639 |
|      3 |   0.3 | 1.842418 | 3.068198 | 0.0588417 | 1.1575816 | 0.0916210 |
|      4 |   0.1 | 3.307682 | 4.559242 | 0.0791524 | 0.6923185 | 0.5592425 |
|      4 |   0.2 | 2.784925 | 4.146991 | 0.0689373 | 1.2150750 | 0.1549855 |
|      4 |   0.3 | 2.375321 | 3.759442 | 0.0610733 | 1.6246787 | 0.2415091 |
|      5 |   0.1 | 4.063403 | 5.451041 | 0.0843424 | 0.9365970 | 0.4510412 |
|      5 |   0.2 | 3.384433 | 4.924717 | 0.0721564 | 1.6155667 | 0.1126374 |
|      5 |   0.3 | 2.876683 | 4.451993 | 0.0631377 | 2.1233174 | 0.5480072 |

``` r
#map %>% ggplot() + geom_sf(aes(fill = Inventering)) + theme_void() +
#  geom_point(data = captures %>% filter(id == "BI041000 ZF-108"), 
#             aes(x = lon, y = lat, size=I(2),stroke=I(0),shape=I(16))) +
#  labs(title = "Territory Size")
```

Simulation

Introduction

*Why are we doing this?*

-Bear Inventory has been performed for several years. -Hunting quotas
are set at a regional level. -Studying bias from double counting bears.
-Estimating number of bears in each of the four regions.

*What have been done (literature review)?*
</li>

-   

-kindberg2011estimating

*How will we approach the problem?*

-Simple model that can be applied to all samples. -Number of spills
poisson distributed. -Bear territory a circle of specified size.

***First Draft***

For several years now the Swedish museum of natural history has been
performing an inventory of the brown bear population in the northern
half of Sweden. The region has been divided into four parts and every
year a different part is investigated and every fifth year no
investigations are performed. The investigation is done by soliciting
hunters among others to collect stool samples from bear that they find
in the forest and send it in for DNA analysis along with information on
where the sample was found.Conventional statistical methods have been
performed on the supplied data and estimated have been made over the
years but now more specific type of analysis has been requested.

(***WRITE SOMETHING ABOUT THE KINDBERG ARTICLE***)

Brown Bears have very large areas that they wander through and pay no
heed to the borders specified by humans. As such a bear might be found
in two different areas in two different years which can cause a bear to
be counted in both areas. This is a problem since hunting quotas are not
set at a national level but a regional level and as such knowing the
exact number of bears in each area is important. As such the purpouse of
my project is to estimate the bias introduced by double counting bears
and also create my own estimate of the total numbers of bears in each
area.

As a large number of bears have only been observed a single time and a
large number of bears have not been observed at all, trying to apply a
standard model for estimating the size and shape of a bears area is
going to be difficult if not impossible. As such for the sample at hand
i will be applying a simplified model that can be applied to all bears
independant of the number of samples we have from each. I will also be
performing a simulation study in which i will be using a more realistic
model for the size and shape of the bears territories. I will be making
the assumption that the number of samples from each bear follows a
poisson distribution and combined with the model on the bears areas i
will estimate the number of bears in each region.

(***Standardize choice of words***)

-Bears Territory/Bears Area

-Area/Region

***Method***

*What tools are we going to use:*

-Rstudio

*Statistical models and assumptions.*

-Poisson Model for number of samples found from each bear which can also
be used to estimate total number of bears.

-Simple model for size of territory (circle of constant size) along with
more complicated one for simulation studies (Bivariate normal) *How are
we going to estimate parameters?*

-Maximum likelihood. *Do we need numerical methods?*

-Numerical Methods for calculating the ratio of a bears territory lies
on which side of the border.

*Results*

-   

*Discussion*

-Try to collect samples from outside Sweden as well.

*References*

-bibtex file

    ## 
    ## To cite R in publications use:
    ## 
    ##   R Core Team (2021). R: A language and environment for statistical
    ##   computing. R Foundation for Statistical Computing, Vienna, Austria.
    ##   URL https://www.R-project.org/.
    ## 
    ## A BibTeX entry for LaTeX users is
    ## 
    ##   @Manual{,
    ##     title = {R: A Language and Environment for Statistical Computing},
    ##     author = {{R Core Team}},
    ##     organization = {R Foundation for Statistical Computing},
    ##     address = {Vienna, Austria},
    ##     year = {2021},
    ##     url = {https://www.R-project.org/},
    ##   }
    ## 
    ## We have invested a lot of time and effort in creating R, please cite it
    ## when using it for data analysis. See also 'citation("pkgname")' for
    ## citing R packages.
