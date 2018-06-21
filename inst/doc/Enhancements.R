## ----setup, echo = FALSE, message=FALSE----------------------------------
library(EGRET)
library(EGRETci)
library(lubridate)
library(dplyr)
library(knitr)

knitr::opts_chunk$set(echo = TRUE, 
                      warning = FALSE, 
                      message = FALSE, 
                      eval = nzchar(Sys.getenv("EGRET_eval")),
                      fig.width=7, fig.height=7)


## ---- echo = FALSE-------------------------------------------------------
load("pairResults2.RData")
load("Chop.OPbase.RData")

## ------------------------------------------------------------------------
bootPairOut2 <- runPairsBoot(eList, pairResults2, nBoot = 10)

## ----eval=FALSE----------------------------------------------------------
#  plotHistogramTrend(eList, eBoot, caseSetUp = NA,
#                     flux = TRUE, xMin = NA, xMax = NA,
#                     xStep = NA, printTitle = TRUE,
#                     cex.main = 1.1, cex.axis = 1.1,
#                     cex.lab = 1.1, col.fill = "grey", ...)

## ---- echo = FALSE-------------------------------------------------------
load("bootPairOut2.RData")

## ------------------------------------------------------------------------
plotHistogramTrend(eList,bootPairOut2, caseSetUp = NA)
plotHistogramTrend(eList,bootPairOut2, caseSetUp = NA, flux = FALSE)

## ----eval=FALSE----------------------------------------------------------
#  CIAnnualResults <- ciCalculations(eList, startSeed = 494817,
#                                    verbose = TRUE, ...)

## ---- echo = FALSE-------------------------------------------------------
load("eListOutChop.RData")

## ---- eval = FALSE-------------------------------------------------------
#  CIAnnualResults <- ciCalculations(eListOut, verbose = FALSE, nBoot = 100, blockLength = 200, widthCI = 90)

## ------------------------------------------------------------------------
plotConcHistBoot(eListOut, CIAnnualResults)
plotFluxHistBoot(eListOut, CIAnnualResults)

## ---- echo = FALSE-------------------------------------------------------
rm(list = ls())
load("Chop.OPbase.RData")
load("groupResults.RData")

## ---- eval = FALSE-------------------------------------------------------
#  bootGroupsOut <- runGroupsBoot(eList, groupResults, nBoot = 100)
#  # Laura: I can't seem to do the next step,
#  # what I want to do here is insert the output from this run
#  # which I have sitting in a txt file called "bootGroupsOutput"
#  # I'd like it inserted here, just don't know how to do it.

## ---- echo = FALSE-------------------------------------------------------
load("bootGroupsOut.RData")

## ------------------------------------------------------------------------
plotGroupsHistogramTrend <- function (eList, eBoot, flux = TRUE, xMin = NA, xMax = NA, 
    xStep = NA, printTitle = TRUE, cex.main = 1.1, cex.axis = 1.1, 
    cex.lab = 1.1, col.fill = "grey", ...) 
{
    periodName <- EGRET::setSeasonLabel(data.frame(PeriodStart = eList$INFO$paStart, 
        PeriodLong = eList$INFO$paLong))
    if ("runGroups" %in% names(attributes(eList)) | "segmentInfo" %in% 
        names(attributes(eList$INFO))) {
        periodName <- paste(periodName, "*")
    }
    if (flux) {
        change <- 100 * eBoot$bootOut$estF/eBoot$bootOut$baseFlux
        reps <- eBoot$pFlux
        xlabel <- "Flux trend, in %"
        titleWord <- "Flux: "
    }
    else {
        change <- 100 * eBoot$bootOut$estC/eBoot$bootOut$baseConc
        reps <- eBoot$pConc
        xlabel <- "Concentration trend, in %"
        titleWord <- "Concentration: "
    }
    group1firstYear <- attr(eBoot,"group1firstYear")
    group1lastYear <- attr(eBoot,"group1lastYear")
    group2firstYear <- attr(eBoot,"group2firstYear")
    group2lastYear <- attr(eBoot,"group2lastYear")
    periodWords <- paste(group2firstYear, "to", group2lastYear,
                         "minus",group1firstYear, "to",group1lastYear, sep=" ")
    titleToPrint <- ifelse(printTitle, paste("Trend magnitude in", 
        eList$INFO$paramShortName, "\nFlow Normalized", titleWord, 
        periodWords, "\n", eList$INFO$shortName, periodName), 
        "")
    minReps <- min(reps, na.rm = TRUE)
    maxReps <- max(reps, na.rm = TRUE)
    xMin <- if (is.na(xMin)) 
        min(-10, minReps)
    else xMin
    xMax <- if (is.na(xMax)) 
        max(10, maxReps)
    else xMax
    xStep <- if (is.na(xStep)) 
        (xMax - xMin)/10
    else xStep
    xSeq <- seq(xMin, xMax, xStep)
    hist(reps, breaks = xSeq, yaxs = "i", xaxs = "i", axes = FALSE, 
        ylab = "", main = titleToPrint, freq = FALSE, xlab = xlabel, 
        col = col.fill, cex.main = cex.main * 0.8, cex.lab = cex.lab, 
        ...)
    abline(v = change, lwd = 3, lty = 2)
    abline(v = 0, lwd = 3)
    box()
    axis(1, tcl = 0.5, labels = TRUE, cex.axis = cex.axis)
    axis(2, tcl = 0.5, labels = TRUE, las = 1, cex.axis = cex.axis)
    title(ylab = "Density", line = 4.5, cex.lab = cex.lab)
    axis(3, tcl = 0.5, labels = FALSE)
    axis(4, tcl = 0.5, labels = FALSE)
}

plotGroupsHistogramTrend(eList, bootGroupsOut, xMin = -30, xMax = 40, xStep = 10)
plotGroupsHistogramTrend(eList, bootGroupsOut, flux=FALSE, xMin = -30, xMax = 40, xStep = 10)

