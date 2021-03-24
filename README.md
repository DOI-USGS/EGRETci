# EGRETci <img src="man/figures/egret-02.png" alt="EGRETci" height="150px" align="right" />

[![R build
status](https://github.com/USGS-R/EGRETci/workflows/R-CMD-check/badge.svg)](https://github.com/USGS-R/EGRETci/actions)
[![CRAN version](http://www.r-pkg.org/badges/version/EGRETci)](https://cran.r-project.org/package=EGRETci)
[![](http://cranlogs.r-pkg.org/badges/EGRETci)](https://cran.r-project.org/package=EGRETci)
[![](http://cranlogs.r-pkg.org/badges/grand-total/EGRETci)](https://cran.r-project.org/package=EGRETci)
[![codecov](https://codecov.io/gh/USGS-R/EGRETci/branch/master/graph/badge.svg)](https://codecov.io/gh/USGS-R/EGRETci)
[![status](https://img.shields.io/badge/USGS-Research-blue.svg)](https://owi.usgs.gov/R/packages.html#research)


This package **EGRETci** implements a set of approaches to the analysis of uncertainty associated with WRTDS trend analysis as implemented in the **EGRET** package.

See: [http://authors.elsevier.com/sd/article/S1364815215300220](http://authors.elsevier.com/sd/article/S1364815215300220) for more details.

The **EGRETci** package is designed to carry out three specific types of tasks.

1)  Evaluate a water quality trend over a specific span of years and produce a variety of tabular results.  This is done with a short workflow involving the functions: `trendSetUp` and `wBT`.  The results come in three forms: 1) console output, which shows the bootstrap replicate process as it is underway and the results when it has finished, 2) a text file that shows the results of the bootstrap analysis (a subset of what is included in the console output), and 3) a set of outputs in a named list called eBoot.  The contents of eBoot are described below.  

2)  Plot histograms of values for the trend magnitudes, expressed in percent change over the specified period, for flow-normalized concentration and flow-normalized flux.  This is done with the function `plotHistogramTrend`.  It depends on outputs contained in eBoot.  Note that there are a number of custom outputs similar to these histograms that can be developed from the contents of eBoot (for example, what is the likelihood that the flow normalized flux decreased by more than 2000 kg/year over the trend period).  Such analyses would require a small amount of script writing by the user.

3)  Plot confidence bands around the computed trends in flow-normalized concentration and flow-normalized flux.  This is done using a function called `ciCalculations` and then, using the output from that function running two functions that produce the confidence band graphics for concentration and flux respectively (`plotConcHistBoot`, and `plotFluxHistBoot`).

## How to cite EGRETci:

```
citation(package = "EGRETci")

To cite EGRETci in publications, please use:

  Hirsch, R.M., Archfield, S.A., De Cicco, L.A., "A bootstrap method for
  estimating uncertainty of water quality trends", Environmental
  Modelling & Software, Vol 73, Nov 2015, p 148-166. doi:
  j.envsoft.2015.07.017

A BibTeX entry for LaTeX users is

  @Article{,
    author = {Robert M. Hirsch and Stacey A. Archfield and Laura A. {De Cicco}},
    title = {A bootstrap method for estimating uncertainty of water quality trends},
    journal = {Journal of Environmental Modelling and Software},
    year = {2015},
  }

```


### Reporting bugs

Please consider reporting bugs and asking questions on the Issues page:
[https://github.com/USGS-R/EGRETci/issues](https://github.com/USGS-R/EGRET/issues)

Follow `@USGS_R` on Twitter for updates on USGS R packages:

[![Twitter Follow](https://img.shields.io/twitter/follow/USGS_R.svg?style=social&label=Follow%20USGS_R)](https://twitter.com/USGS_R)

### Subscribe

Please email questions, comments, and feedback to: 
egret_comments@usgs.gov

Additionally, to subscribe to an email list concerning updates to these R packages, please send a request to egret_comments@usgs.gov.

### Code of Conduct

We want to encourage a warm, welcoming, and safe environment for contributing to this project. See the [code of conduct](https://github.com/USGS-R/EGRETci/blob/master/CONDUCT.md) for more information.


## Model Archive

When using the `WRTDS` model (and corresponding confidence intervals), it is important to be able to reproduce the results in the future. The following version of R and package dependencies were used most recently to pass the embedded tests within this package. There is no guarantee of reproducible results using future versions of R or updated versions of package dependencies; however, we will make diligent efforts to test and update future modeling environments.

```
devtools::session_info()
Session info ----------------------------------------------------------------------
 setting  value                       
 version  R version 3.4.1 (2017-06-30)
 system   x86_64, mingw32             
 ui       RStudio (1.0.143)           
 language (EN)                        
 collate  English_United States.1252  
 tz       America/Chicago             
 date     2017-07-27                  

Packages --------------------------------------------------------------------------
 package       * version    date       source        
 assertthat      0.2.0      2017-04-11 CRAN (R 3.4.0)
 base          * 3.4.1      2017-06-30 local         
 bindr           0.1        2016-11-13 CRAN (R 3.4.0)
 bindrcpp        0.2        2017-06-17 CRAN (R 3.4.0)
 binom           1.1-1      2014-01-02 CRAN (R 3.4.0)
 codetools       0.2-15     2016-10-05 CRAN (R 3.4.1)
 compiler        3.4.1      2017-06-30 local         
 curl            2.8.1      2017-07-21 CRAN (R 3.4.1)
 dataRetrieval   2.7.2      2017-05-23 CRAN (R 3.4.1)
 datasets      * 3.4.1      2017-06-30 local         
 devtools        1.13.2     2017-06-02 CRAN (R 3.4.1)
 digest          0.6.12     2017-01-27 CRAN (R 3.4.0)
 dotCall64       0.9-04     2016-10-07 CRAN (R 3.4.1)
 dplyr           0.7.2      2017-07-20 CRAN (R 3.4.1)
 EGRET           2.6.0      2017-07-27 CRAN (R 3.4.1)    
 EGRETci       * 1.0.5      2017-07-24 CRAN (R 3.4.1)      
 fields          9.0        2017-06-06 CRAN (R 3.4.0)
 foreach         1.4.3      2015-10-13 CRAN (R 3.4.0)
 glue            1.1.1      2017-06-21 CRAN (R 3.4.0)
 graphics      * 3.4.1      2017-06-30 local         
 grDevices     * 3.4.1      2017-06-30 local         
 grid            3.4.1      2017-06-30 local         
 hms             0.3        2016-11-22 CRAN (R 3.4.0)
 httr            1.2.1      2016-07-03 CRAN (R 3.3.1)
 iterators       1.0.8      2015-10-13 CRAN (R 3.4.0)
 jsonlite        1.5        2017-06-01 CRAN (R 3.4.0)
 lattice         0.20-35    2017-03-25 CRAN (R 3.4.1)
 lubridate       1.6.0      2016-09-13 CRAN (R 3.4.0)
 magrittr        1.5        2014-11-22 CRAN (R 3.4.0)
 maps            3.2.0      2017-06-08 CRAN (R 3.4.0)
 Matrix          1.2-10     2017-05-03 CRAN (R 3.4.1)
 memoise         1.1.0      2017-04-21 CRAN (R 3.4.0)
 methods       * 3.4.1      2017-06-30 local         
 parallel        3.4.1      2017-06-30 local         
 pkgconfig       2.0.1      2017-03-21 CRAN (R 3.4.0)
 plyr            1.8.4      2016-06-08 CRAN (R 3.3.1)
 R6              2.2.2      2017-06-17 CRAN (R 3.4.0)
 Rcpp            0.12.12    2017-07-15 CRAN (R 3.4.1)
 readr           1.1.1      2017-05-16 CRAN (R 3.4.0)
 reshape2        1.4.2      2016-10-22 CRAN (R 3.4.0)
 rlang           0.1.1      2017-05-18 CRAN (R 3.4.1)
 spam            2.1-1      2017-07-03 CRAN (R 3.4.1)
 splines         3.4.1      2017-06-30 local         
 stats         * 3.4.1      2017-06-30 local         
 stringi         1.1.5      2017-04-07 CRAN (R 3.4.0)
 stringr         1.2.0      2017-02-18 CRAN (R 3.4.0)
 survival        2.41-3     2017-04-04 CRAN (R 3.4.0)
 tibble          1.3.3      2017-05-28 CRAN (R 3.4.0)
 tools           3.4.1      2017-06-30 local         
 utils         * 3.4.1      2017-06-30 local         
 withr           1.0.2      2016-06-20 CRAN (R 3.3.1)
 xml2            1.1.1      2017-01-24 CRAN (R 3.4.0)

```

## Disclaimer

This software has been approved for release by the U.S. Geological Survey (USGS). Although the software has been subjected to rigorous review, the USGS reserves the right to update the software as needed pursuant to further analysis and review. No warranty, expressed or implied, is made by the USGS or the U.S. Government as to the functionality of the software and related material nor shall the fact of release constitute any such warranty. Furthermore, the software is released on condition that neither the USGS nor the U.S. Government shall be held liable for any damages resulting from its authorized or unauthorized use.
