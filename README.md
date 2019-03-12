<!-- README.md is generated from README.Rmd. Please edit that file -->
paleoTS
=======

paleoTS allows the user to analyze paleontological time-series implementing many common models that are considered by paleontologists when looking at trait changes in a species over time, including random walks, directional change, Ornstein-Uhlenbeck models of adaptation, and stasis. In addition, more complex models are available that allow for punctuated change, and for shifts in dynamics within a time-series.

Example
-------

As a simple example, we simulate directional change, and then fit three common models used in the literature to these data:

``` r
 library(paleoTS)
 x <- sim.GRW(ns = 20, ms = 1, vs = 0.3)
 fit3models(x)
#> 
#> Comparing 3 models [n = 20, method = Joint]
#> 
#>             logL K      AICc    dAICc Akaike.wt
#> GRW    -15.10736 3  37.71472  0.00000         1
#> URW    -28.45691 2  61.61970 23.90498         0
#> Stasis -61.93319 2 128.57227 90.85754         0
```
