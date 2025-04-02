# EGRETci <img src="man/figures/egret-02.png" alt="EGRET" style="width:90px;height:auto;" align="right" class="logo" />

This package **EGRETci** implements a set of approaches to the analysis
of uncertainty associated with WRTDS trend analysis as implemented in
the **EGRET** package.

See: <https://doi.org/10.1016/j.envsoft.2015.07.017> for more details.

The **EGRETci** package is designed to carry out three specific types of
tasks.

1.  Evaluate a water quality trend over a specific span of years and
    produce a variety of tabular results. This is done with a short
    workflow involving the functions: `trendSetUp` and `wBT`. The
    results come in three forms: 1) console output, which shows the
    bootstrap replicate process as it is underway and the results when
    it has finished, 2) a text file that shows the results of the
    bootstrap analysis (a subset of what is included in the console
    output), and 3) a set of outputs in a named list called eBoot. The
    contents of eBoot are described below.

2.  Plot histograms of values for the trend magnitudes, expressed in
    percent change over the specified period, for flow-normalized
    concentration and flow-normalized flux. This is done with the
    function `plotHistogramTrend`. It depends on outputs contained in
    eBoot. Note that there are a number of custom outputs similar to
    these histograms that can be developed from the contents of eBoot
    (for example, what is the likelihood that the flow normalized flux
    decreased by more than 2000 kg/year over the trend period). Such
    analyses would require a small amount of script writing by the user.

3.  Plot confidence bands around the computed trends in flow-normalized
    concentration and flow-normalized flux. This is done using a
    function called `ciCalculations` and then, using the output from
    that function running two functions that produce the confidence band
    graphics for concentration and flux respectively
    (`plotConcHistBoot`, and `plotFluxHistBoot`).

## How to cite EGRETci:

``` r
citation(package = "EGRETci")
#> To cite EGRETci in publications, please use:
#> 
#>   Hirsch, R.M., De Cicco, L.A., Archfield, S.A., Murphy, J.C., 2023,
#>   Exploration and Graphics for RivEr Trends (EGRET) Uncertainty and
#>   Confidence Intervals, version 2.0.5, doi:10.5066/P9CC9JEX
#> 
#> A BibTeX entry for LaTeX users is
#> 
#>   @Manual{,
#>     author = {Robert M. Hirsch and Laura A. {De Cicco} and Stacey A. Archfield and Jennifer C. Murphy},
#>     title = {EGRETci},
#>     publisher = {U.S. Geological Survey},
#>     year = {2023},
#>     url = {https://pubs.usgs.gov/tm/04/a10/},
#>   }
```

### Reporting bugs

Please consider reporting bugs and asking questions on the Issues page:
[https://github.com/DOI-USGS/EGRETci/issues](https://github.com/DOI-USGS/EGRET/issues)

### Code of Conduct

We want to encourage a warm, welcoming, and safe environment for
contributing to this project. See the [code of
conduct](https://github.com/DOI-USGS/EGRETci/blob/main/CONDUCT.md) for
more information.

# Disclaimer

This software is preliminary or provisional and is subject to revision.
It is being provided to meet the need for timely best science. The
software has not received final approval by the U.S. Geological Survey
(USGS). No warranty, expressed or implied, is made by the USGS or the
U.S. Government as to the functionality of the software and related
material nor shall the fact of release constitute any such warranty. The
software is provided on the condition that neither the USGS nor the U.S.
Government shall be held liable for any damages resulting from the
authorized or unauthorized use of the software.
