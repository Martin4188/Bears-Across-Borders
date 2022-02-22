Bears Across Borders
================
Martin Andersson
2022-01-29

``` r
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
  Other <- 1 - (R1 + R2 + R3 + R4)
  
  list(Region1=R1, Region2=R2, Region3=R3, Region4=R4, OtherRegion=Other) %>%
    as_tibble() %>%
    round(3) %>%
    return()
}
```

``` r
TestTibble <- capturesFemale %>%
  select(id, meanlon, meanlat) %>%
  mutate(Ratio = RegionRatio(meanlon, meanlat, 121381)) %>%
  mutate(Region1 = Ratio[[1]], Region2 = Ratio[[2]], Region3 = Ratio[[3]], Region4 = Ratio[[4]], RegionOther = Ratio[[5]]) %>%
  select(-Ratio)

TestTibble %>% 
  arrange(desc(RegionOther))
```

    ## # A tibble: 1,806 x 8
    ## # Groups:   id [1,806]
    ##    id                  meanlon  meanlat Region1 Region2 Region3 Region4 RegionOther
    ##    <chr>                 <dbl>    <dbl>   <dbl>   <dbl>   <dbl>   <dbl>       <dbl>
    ##  1 BI408946 S17-008    381306  6742271    0.373   0       0.011   0           0.616
    ##  2 BI408939 S17-001 +  382195. 6740997    0.375   0       0.01    0           0.615
    ##  3 BI407756 BD16-236   851500. 7560285    0       0       0       0.431       0.569
    ##  4 BI407674 BD16-149 + 856734. 7533838    0       0       0       0.48        0.52 
    ##  5 BI409386 X17-273    605948. 6739078.   0.484   0       0       0           0.516
    ##  6 BI407553 BD16-020 + 858771  7528040    0       0       0       0.485       0.515
    ##  7 BI409115 X17-002 +  602058. 6733611.   0.489   0       0       0           0.511
    ##  8 BI405021 NT117 +    439818  7166068.   0       0.153   0.338   0           0.51 
    ##  9 BI405768 HE165      383940  6788807    0.404   0       0.093   0           0.503
    ## 10 BI409059 W17-111    364428. 6864046.   0.227   0       0.277   0           0.496
    ## # ... with 1,796 more rows

``` r
map %>% ggplot() + geom_sf(aes(fill = Inventering)) + theme_void() +
  geom_point(data = captures %>% filter(id == "BI041000 ZF-108"), 
             aes(x = lon, y = lat, size=I(2),stroke=I(0),shape=I(16))) +
  labs(title = "Territory Size")
```

![](Bears-Across-Borders_files/figure-gfm/unnamed-chunk-5-1.png)<!-- -->

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
