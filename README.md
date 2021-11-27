
<!-- README.md is generated from README.Rmd. Please edit that file -->
### Overview

**MIIVsem** is an R package for estimating structural equation models using model-implied instrumental variables.

Version 0.5.2 includes the following features:

-   Estimation of latent variable and simultaneous equation models.
-   Model-implied and traditional instrumental variable estimation.
-   Equation level specification tests.
-   Efficient computation from covariance matrix input.
-   Polychoric instrumental variable estimation for endogenous categorical variables.
-   Impose and test within- and across-equation parameter restrictions.
-   Bootstrap standard errors.
-   Variance and covariance parameter estimation.

### Installation

In R you can install MIIVsem from CRAN as follows:

``` r
install.packages("MIIVsem")
```

### Usage

MIIVsem uses a subset of the model syntax employed by [lavaan](https://lavaan.ugent.be/) (Rosseel, 2012) for model specification. The following model syntax operators are currently supported:

| Operators |                                                                    |
|-----------|--------------------------------------------------------------------|
| =~        | Used for expressing measurement relations, read as 'measured by.'  |
| ~         | Used For expressing regression relations, read as 'regressed on.'  |
| ~~        | For specifying variances and covariances, read as 'covaries with.' |
| \*        | For assigning equality or numerical constraints.                   |

### Model Syntax

**Example using Syntax Operators**

In the model below, `L1 =~ Z1 + Z2 + Z3` indicates the latent variable L1 is measured by 3 indicators, `Z1`, `Z2`, and `Z3`. Likewise, `L2` is measured by 3 indicators, `Z4`, `Z5`, and `Z6`. The statement `L1 ~ L2` specifies latent variable `L1` is regressed on latent variable `L2`. `Z1 ~~ Z2` indicates the error of `Z2` is allowed to covary with the error of `Z3`. The label `LA3` prepended to `Z3` and `Z6` in the measurement model equations constrains the factor loadings for `Z3` and `Z6` to equality.

``` r
model <- '
   L1 =~ Z1 + Z2 + LA3*Z3
   L2 =~ Z4 + Z5 + LA3*Z6
   L1  ~ L2
   Z2 ~~ Z3
'  
```

**Scaling Indicators**

Following the lavaan model syntax, latent variables are defined using the `=~` operator. For first order factors, the scaling indicator chosen is the first observed variable on the RHS of an equation. For the model below `Z1` would be chosen as the scaling indicator for `L1` and `Z4` would be chosen as the scaling indicator for `L2`.

``` r
model <- '
   L1 =~ Z1 + Z2 + Z3
   L2 =~ Z4 + Z5 + Z6
'
```

**Equality Constraints and Parameter Restrictions**

Within- and across-equation equality constraints on the factor loading and regression coefficients can be imposed directly in the model syntax. To specify equality constraints between different parameters equivalent labels should be prepended to the variable name using the `*` operator. For example, we could constrain the factor loadings for two non-scaling indicators of latent factor `L1` to equality using the following model syntax.

``` r
model <- '
   L1 =~ Z1 + LA2*Z2 + LA2*Z3
   L2 =~ Z4 + Z5 + Z6
'
```

Researchers also can constrain the factor loading and regression coefficients to specific numeric values in a similar fashion. Below we constrain the regression coefficient of `L1` on `L2` to `1`.

``` r
model <- '
   L1 =~ Z1 + Z2 + Z3
   L2 =~ Z4 + Z5 + Z6
   L3 =~ Z7 + Z8 + Z9
   L1  ~ 1*L2 + L3
'
```

**Higher Order Factor Model**

In the model below, the scaling indicator for the higher-order factor `H1` is taken to be `Z1`, the scaling indicator that would have been assigned to the first lower-order factor `L1`. The intercepts for lower-order latent variables are set to zero, by default

``` r
model <- '
      H1 =~ L1 + L2 + L3
      L1 =~ Z1 + Z2 + Z3
      L2 =~ Z4 + Z5 + Z6
      L3 =~ Z7 + Z8 + Z9
   '
```

**Model Defaults**

In addition to those relationships specified in the model syntax MIIVsem will automatically include the intercepts of any observed or latent endogenous variable. The intercepts for any scaling indicators and lower-order latent variables are set to zero. Covariances among exogenous latent and observed variables are included by default. Where appropriate the covariances of the errors of latent and observed dependent variables are also included in the model specification. These defaults correspond to those used by lavaan and `auto = TRUE`, except that endogenous latent variable intercepts are estimated by default, and the intercepts of scaling indicators are fixed to zero.

### Getting Started

**MIIV Search**

Researchers typically search for instrumental variables external to the model. The key property of valid instruments is that they are uncorrelated with equation error. The MIIV approach proposed in Bollen (1996) finds instruments among observed variables already in the model. Here, the model specification itself implies which observed variables are uncorrelated with the equation disturbance.

Using the industrialization-democracy example from Bollen (1989) we illustrate the MIIV Search:

``` r
library(MIIVsem)

model <- '

    Eta1 =~ y1 + y2  + y3  + y4  
    Eta2 =~ y5 + y6  + y7  + y8    
    Xi1  =~ x1 + x2 + x3 

    Eta1 ~ Xi1  
    Eta2 ~ Xi1 
    Eta2 ~ Eta1 

    y1   ~~ y5
    y2   ~~ y4
    y2   ~~ y6
    y3   ~~ y7
    y4   ~~ y8
    y6   ~~ y8 
  '
```

``` r
miivs(model)
#> Model Equation Information 
#> 
#>  LHS RHS    MIIVs                             
#>  y2  y1     y3, y7, y8, x1, x2, x3            
#>  y3  y1     y2, y4, y6, y8, x1, x2, x3        
#>  y4  y1     y3, y6, y7, x1, x2, x3            
#>  y6  y5     y3, y4, y7, x1, x2, x3            
#>  y7  y5     y2, y4, y6, y8, x1, x2, x3        
#>  y8  y5     y2, y3, y7, x1, x2, x3            
#>  x2  x1     y1, y2, y3, y4, y6, y5, y7, y8, x3
#>  x3  x1     y1, y2, y3, y4, y6, y5, y7, y8, x2
#>  y1  x1     x2, x3                            
#>  y5  x1, y1 y2, y3, y4, x2, x3
```

**MIIV-2SLS Estimation**

We can also estimate the industrialization-democracy model using MIIV-2SLS:

``` r
miive(model, bollen1989a)
#> MIIVsem (0.5.2) results 
#> 
#> Number of observations                                                     75
#> Number of equations                                                        10
#> Estimator                                                           MIIV-2SLS
#> Standard Errors                                                      standard
#> Missing                                                              listwise
#> 
#> 
#> Parameter Estimates:
#> 
#> 
#> STRUCTURAL COEFFICIENTS:
#>                    Estimate  Std.Err  z-value  P(>|z|)   Sargan   df   P(Chi)
#>   Eta1 =~                                                                    
#>     y1                1.000                                                  
#>     y2                1.139    0.179    6.371    0.000    8.409    5    0.135
#>     y3                0.969    0.140    6.924    0.000    5.874    6    0.437
#>     y4                1.210    0.139    8.713    0.000    4.276    5    0.510
#>   Eta2 =~                                                                    
#>     y5                1.000                                                  
#>     y6                1.051    0.165    6.377    0.000    8.712    5    0.121
#>     y7                1.180    0.151    7.814    0.000    9.538    6    0.146
#>     y8                1.203    0.154    7.798    0.000    2.795    5    0.731
#>   Xi1 =~                                                                     
#>     x1                1.000                                                  
#>     x2                2.078    0.128   16.171    0.000    8.301    8    0.405
#>     x3                1.751    0.149   11.782    0.000    8.738    8    0.365
#>                                                                              
#>   Eta1 ~                                                                     
#>     Xi1               1.261    0.426    2.962    0.003    0.503    1    0.478
#>   Eta2 ~                                                                     
#>     Xi1               1.123    0.312    3.598    0.000    0.801    3    0.849
#>     Eta1              0.724    0.101    7.140    0.000                       
#> 
#> INTERCEPTS:
#>                    Estimate  Std.Err  z-value  P(>|z|)   
#>     Eta1             -0.909    2.170   -0.419    0.675   
#>     Eta2             -4.499    1.424   -3.160    0.002   
#>     x1                0.000                              
#>     x2               -5.711    0.654   -8.727    0.000   
#>     x3               -5.292    0.758   -6.985    0.000   
#>     y1                0.000                              
#>     y2               -1.969    1.044   -1.886    0.059   
#>     y3                1.265    0.814    1.553    0.120   
#>     y4               -2.160    0.814   -2.654    0.008   
#>     y5                0.000                              
#>     y6               -2.418    0.909   -2.659    0.008   
#>     y7                0.135    0.830    0.163    0.870   
#>     y8               -2.137    0.853   -2.505    0.012
```

### Additional Features

**Estimation from Sample Moments**

``` r
sample.cov  <- cov(bollen1989a)
sample.mean <- colMeans(bollen1989a)
sample.nobs <- nrow(bollen1989a)

miive(model, sample.cov = sample.cov, sample.mean = sample.mean, sample.nobs = sample.nobs)
```

**Bootstrap Standard Errors** (Version 0.5.2)

``` r
microbenchmark::microbenchmark(
  fit <- miive(model, bollen1989a, se = "boot", bootstrap = 100L),
  times = 100L
)
#> Unit: milliseconds
#>                                                             expr      min
#>  fit <- miive(model, bollen1989a, se = "boot", bootstrap = 100L) 193.4475
#>        lq     mean   median      uq      max neval
#>  203.1075 209.0176 207.8951 212.088 340.4677   100
```

**Categorical Endogenous Variables (Bollen & Maydeu-Olivares (2007))**

``` r
model <- ' 
    female.access =~ access1 + access2 + access3 
    male.access   =~ access4 + access5 + access6 
'

miive(model, bollen1996, ordered = c("access1", "access2","access3", "access4", "access5", "access6"))
```

### Replication of Textbook Results

Following [Henningsen and Hamann (2007)](https://www.jstatsoft.org/article/view/v023i04) we replicate textbook results from

**Klein's Model I** (Greene, 2003, p.381)

``` r
data("KleinI", package = "systemfit")

model <- '
  consump  ~ corpProf + corpProfLag + wages
  invest   ~ corpProf + corpProfLag + capitalLag
  privWage ~ gnp + gnpLag + trend
'

instruments <- '
  consump  ~ govExp + taxes + govWage + trend + capitalLag + corpProfLag + gnpLag
  invest   ~ govExp + taxes + govWage + trend + capitalLag + corpProfLag + gnpLag
  privWage ~ govExp + taxes + govWage + trend + capitalLag + corpProfLag + gnpLag
'

fit <- miive(model, KleinI, instruments, miiv.check = FALSE)
estimatesTable(fit)
#>         lhs op         rhs         est         se          z       pvalue
#> 1   consump  ~    corpProf  0.01730221 0.11804941  0.1465675 8.834734e-01
#> 2   consump  ~ corpProfLag  0.21623404 0.10726796  2.0158306 4.381770e-02
#> 3   consump  ~       wages  0.81018270 0.04024971 20.1289055 4.119943e-90
#> 4    invest  ~    corpProf  0.15022182 0.17322929  0.8671849 3.858407e-01
#> 5    invest  ~ corpProfLag  0.61594358 0.16278539  3.7837767 1.544664e-04
#> 6    invest  ~  capitalLag -0.15778764 0.03612624 -4.3676741 1.255767e-05
#> 7  privWage  ~         gnp  0.43885907 0.03563192 12.3164596 7.386658e-35
#> 8  privWage  ~      gnpLag  0.14667382 0.03883613  3.7767360 1.588970e-04
#> 9  privWage  ~       trend  0.13039569 0.02914098  4.4746500 7.653659e-06
#> 10  consump ~1             16.55475577 1.32079242 12.5339573 4.867207e-36
#> 11   invest ~1             20.27820894 7.54270590  2.6884528 7.178398e-03
#> 12 privWage ~1              1.50029689 1.14778020  1.3071291 1.911689e-01
```
