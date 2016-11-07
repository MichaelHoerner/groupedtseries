
<!-- README.md is generated from README.Rmd. Please edit that file -->
Analysis grouped time series
============================

Load package

``` r
library(groupedtseries)
#> Loading required package: tseries
```

Example

``` r
library(datasets)
data("EuStockMarkets")

m <- grouped_arma(EuStockMarkets[, 1])  
summary(m)
#> 
#> Call:
#> grouped_arma(x = EuStockMarkets[, 1])
#> 
#> Model:
#> ARMA(1,1)
#> 
#> Residuals:
#>       Min        1Q    Median        3Q       Max 
#> -229.8672  -11.7200   -0.6021   12.4986  156.8212 
#> 
#> Coefficient(s):
#>             Estimate  Std. Error  t value Pr(>|t|)    
#> ar1        1.0013518   0.0006936 1443.695   <2e-16 ***
#> ma1       -0.0024755   0.0236573   -0.105    0.917    
#> intercept -1.3485286   1.9080786   -0.707    0.480    
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#> Fit:
#> sigma^2 estimated as 1054,  Conditional Sum-of-Squares = 1958190,  AIC = 18230.56
```
