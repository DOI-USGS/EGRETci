#' Graph of annual concentration, flow normalized concentration, 
#' and confidence bands for flow normalized concentrations
#' 
#' Uses the output of \code{\link[EGRET]{modelEstimation}} in the EGRET package (results in the named 
#' list eList), and the data frame CIAnnualResults (produced by the function ciCalculations in the EGRETci package 
#' using scripts described in the EGRETci vignette) to produce a graph of annual 
#' concentration, flow normalized concentration, and confidence bands for 
#' flow-normalized concentrations.  In addition to the arguments listed below, 
#' it will accept any additional arguments that are listed for the EGRET function 
#' \code{\link[EGRET]{plotConcHist}}.
#'
#' @param eList named list with at least the Daily, Sample, and INFO dataframes. Created from the EGRET package, after running \code{\link[EGRET]{modelEstimation}}.
#' @param CIAnnualResults data frame generated from ciBands (includes nBoot, probs, and blockLength attributes).
#' @param yearStart numeric is the calendar year containing the first estimated annual value to be plotted, default is NA (which allows it to be set automatically by the data).
#' @param yearEnd numeric is the calendar year just after the last estimated annual value to be plotted, default is NA (which allows it to be set automatically by the data).
#' @param plotFlowNorm logical variable if TRUE flow normalized concentration line is plotted, if FALSE not plotted, default is TRUE. 
#' @param col.pred character color of line for flow-normalized concentration and for the confidence limits, default is "green". 
#' @param concMax numeric specifying the maximum value to be used on the vertical axis, default is NA (which allows it to be set automatically by the data).
#' @param plotAnnual logical variable if \code{TRUE}, annual mean concentration points from WRTDS output are plotted, if \code{FALSE} not plotted. 
#' @param plotGenConc logical variable. If \code{TRUE}, annual mean concentration points from WRTDS_K output are plotted, if \code{FALSE} not plotted.
#' @param cex numeric value giving the amount by which plotting symbols should be magnified, default = 0.8.
#' @param cex.axis numeric value of magnification to be used for axis annotation relative to the current setting of cex, default = 1.1.
#' @param lwd numeric magnification of line width, default = 2.
#' @param col color of annual mean points on plot, see ?par 'Color Specification', default = "black".
#' @param col.gen color of annual mean points for WRTDS_K output on plot, see ?par 'Color Specification', default = "red". 
#' @param printTitle logical print title of the plot, default = TRUE.
#' @param cex.main numeric value of magnification to be used for plot title, default = 1.1.
#' @param customPar logical defaults to FALSE. If TRUE, par() should be set by user before calling this function 
#' (for example, adjusting margins with par(mar=c(5,5,5,5))). If customPar FALSE, EGRETci chooses the best margins.
#' @param \dots graphical parameters
#' @export
#' @importFrom graphics title
#' @importFrom graphics lines
#' @examples
#' library(EGRET)
#' eList <- Choptank_eList
#' CIAnnualResults <- Choptank_CIAnnualResults
#' plotConcHistBoot(eList, CIAnnualResults)
#' plotConcHistBoot(eList, CIAnnualResults, yearStart=1990, yearEnd=2002)
#' # Very long-running function:
#' \dontrun{
#' CIAnnualResults <- ciCalculations(eList, nBoot = 100, blockLength = 200)
#' plotConcHistBoot(eList, CIAnnualResults)
#' }
plotConcHistBoot <- function (eList, CIAnnualResults, yearStart = NA, yearEnd = NA, 
                              plotFlowNorm = TRUE, col.pred = "green", concMax = NA,
                              plotAnnual = TRUE, plotGenConc = FALSE,
                              cex = 0.8, cex.axis = 1.1,  lwd = 2, 
                              col = "black", col.gen = "red", customPar = FALSE, 
                              printTitle = TRUE, cex.main = 1.1, ...){
  
  nBoot <- attr(CIAnnualResults, "nBoot")
  blockLength <- attr(CIAnnualResults, "blockLength")
  probs <- attr(CIAnnualResults, "probs")
  
  widthCI <- (max(probs) - min(probs))*100
  
  if(plotGenConc){
    if(!all((c("GenFlux","GenConc") %in% names(eList$Daily)))){
      stop("This option requires running WRTDSKalman on eList")
    }
    
  } 
  
  localAnnualResults <- EGRET::setupYears(paStart = eList$INFO$paStart, paLong = eList$INFO$paLong,
                                          localDaily = eList$Daily)
  periodName <- EGRET::setSeasonLabel(localAnnualResults)
  if("runSeries" %in% names(attributes(eList)) |
     "segmentInfo" %in% names(attributes(eList$INFO))){
    periodName <- paste(periodName, "*")
  }
  
  title3 <- paste(widthCI,"% CI on FN Concentration, Replicates =",nBoot,"Block=",blockLength,"days")
  
  title <- paste(eList$INFO$shortName, " ", eList$INFO$paramShortName, 
                 "\n", periodName, "\n",title3)
  
  dataStart <- min(eList$Sample$DecYear, na.rm = TRUE)
  dataStartPad <- dataStart - 0.5
  
  if(is.na(yearStart)){
    yearStart <- dataStartPad
  } else {
    yearStart <- max(yearStart, dataStartPad)
    
  }
  
  dataEnd <- max(eList$Sample$DecYear, na.rm = TRUE)
  dataEndPad <- dataEnd + 0.5
  
  if(is.na(yearEnd)){
    yearEnd <- dataEndPad
  } else {
    yearEnd <- min(yearEnd, dataEndPad)
  }
  
  if(is.na(concMax)){
    numYears <- length(localAnnualResults$DecYear)
    
    subAnnualResults <- localAnnualResults[localAnnualResults$DecYear>=yearStart & localAnnualResults$DecYear <= yearEnd,]
    
    annConc <- subAnnualResults$Conc
    if(plotGenConc){
      concMax <- 1.05*max(c(CIAnnualResults$FNConcHigh,
                            annConc,
                            subAnnualResults$GenConc), na.rm=TRUE)
    } else {
      concMax <- 1.05*max(c(CIAnnualResults$FNConcHigh,
                            annConc), na.rm=TRUE)      
    }

  }
  
  EGRET::plotConcHist(eList, yearStart = yearStart, yearEnd = yearEnd,
                      col.pred=col.pred, printTitle=FALSE, cex.axis = cex.axis,
                      cex = cex, cex.main = cex.main, col.gen = col.gen,
                      plotGenConc = plotGenConc, plotAnnual = plotAnnual,
                      plotFlowNorm = plotFlowNorm, concMax = concMax, 
                      lwd = lwd, customPar = customPar, ...)
  if(printTitle) {
    title(main=title, cex.main=cex.main)
  }
  
  if(!is.na(yearStart)){
    CIAnnualResults <- CIAnnualResults[CIAnnualResults$Year >= yearStart, ]
  }
  
  if(!is.na(yearEnd)){
    CIAnnualResults <- CIAnnualResults[CIAnnualResults$Year <= yearEnd, ]
  }
  
  lines(CIAnnualResults$Year, CIAnnualResults$FNConcLow,lty=2,col=col.pred)
  lines(CIAnnualResults$Year, CIAnnualResults$FNConcHigh, lty=2,col=col.pred)
  
}

#' Graph of annual flux, flow normalized flux, and confidence bands for flow normalized flux
#'
#' Uses the output of \code{\link[EGRET]{modelEstimation}} in the EGRET package (results in the named list eList), 
#' and the data frame CIAnnualResults (produced by EGRETci package using scripts described in 
#' the vignette) to produce a graph of annual flux, flow normalized flux, and confidence bands 
#' for flow-normalized flux. In addition to the arguments listed below, it will accept any 
#' additional arguments that are listed for the EGRET function \code{\link[EGRET]{plotFluxHist}}.
#'
#' @param eList named list with at least the Daily, Sample, and INFO dataframes. Created from the EGRET package, after running \code{\link[EGRET]{modelEstimation}}.
#' @param CIAnnualResults data frame from ciBands (needs nBoot, probs, and blockLength attributes).
#' @param yearStart numeric is the calendar year containing the first estimated annual value to be plotted, default is NA (which allows it to be set automatically by the data).
#' @param yearEnd numeric is the calendar year just after the last estimated annual value to be plotted, default is NA (which allows it to be set automatically by the data).
#' @param plotFlowNorm logical variable if TRUE flow normalized flux line is plotted, if FALSE not plotted, default is TRUE. 
#' @param col.pred character color of line for flow-normalized flux and for the confidence limits, default is "green".
#' @param fluxUnit integer representing entry in pre-defined fluxUnit class array. \code{\link[EGRET]{printFluxUnitCheatSheet}}
#' @param fluxMax numeric specifying the maximum value to be used on the vertical axis, default is NA (which allows it to be set automatically by the data), uses units specificed by fluxUnit. 
#' @param plotAnnual logical variable if \code{TRUE}, annual mean flux points from WRTDS output are plotted, if \code{FALSE} not plotted. 
#' @param plotGenFlux logical variable. If \code{TRUE}, annual mean flux points from WRTDS_K output are plotted, if \code{FALSE} not plotted. 
#' @param cex numeric value giving the amount by which plotting symbols should be magnified, default = 0.8.
#' @param cex.axis numeric magnification to be used for axis annotation relative to the current setting of cex, default = 1.1.
#' @param lwd numeric magnification of line width, default = 2.
#' @param col color of annual mean points on plot, see ?par 'Color Specification', default = "black".
#' @param col.gen color of annual mean points for WRTDS_K output on plot, see ?par 'Color Specification', default = "red".
#' @param printTitle logical print title of the plot, default = TRUE.
#' @param cex.main numeric value of magnification to be used for main titles relative to the current setting of cex, default = 1.1. 
#' @param customPar logical defaults to FALSE. If TRUE, par() should be set by user before calling this function 
#' (for example, adjusting margins with par(mar=c(5,5,5,5))). If customPar FALSE, EGRET chooses the best margins.


#' @param cex.main numeric title scale
#' @param \dots graphical parameters
#' @export
#' @importFrom EGRET fluxConst
#' @importFrom graphics lines
#' @importFrom graphics title
#' @examples
#' library(EGRET)
#' eList <- Choptank_eList 
#' CIAnnualResults <- Choptank_CIAnnualResults
#' plotFluxHistBoot(eList, CIAnnualResults, fluxUnit=5)
#' 
#' \dontrun{
#' CIAnnualResults <- ciCalculations(eList, nBoot = 100, blockLength = 200)
#' plotFluxHistBoot(eList, CIAnnualResults, fluxUnit=5)
#' }
plotFluxHistBoot <- function (eList, CIAnnualResults, 
                              yearStart = NA, yearEnd = NA,
                              fluxUnit = 9, fluxMax = NA,
                              plotFlowNorm = TRUE, col.pred = "green",
                              plotAnnual = TRUE, plotGenFlux = FALSE,
                              cex = 0.8, cex.axis = 1.1,  lwd = 2, 
                              col = "black", col.gen = "red", cex.main = 1.1,
                              printTitle = TRUE, customPar = FALSE, ...){
  
  nBoot <- attr(CIAnnualResults, "nBoot")
  blockLength <- attr(CIAnnualResults, "blockLength")
  probs <- attr(CIAnnualResults, "probs")
  
  widthCI <- (max(probs) - min(probs))*100
  
  localAnnualResults <- EGRET::setupYears(paStart = eList$INFO$paStart, paLong = eList$INFO$paLong,
                                          localDaily = eList$Daily)
  periodName <- EGRET::setSeasonLabel(localAnnualResults)
  if("runSeries" %in% names(attributes(eList)) |
     "segmentInfo" %in% names(attributes(eList$INFO))){
    periodName <- paste(periodName, "*")
  }
  
  if(plotGenFlux){
    if(!all((c("GenFlux","GenConc") %in% names(eList$Daily)))){
      stop("This option requires running WRTDS_K on eList")
    }
  } 
  
  title3 <- paste(widthCI,"% CI on FN Flux, Replicates =",nBoot,", Block=",blockLength,"days")
  
  title <- paste(eList$INFO$shortName, " ", eList$INFO$paramShortName, 
                 "\n", periodName, "\n",title3)
  
  if (is.numeric(fluxUnit)) {
    fluxUnit <- fluxConst[shortCode = fluxUnit][[1]]
  } else if (is.character(fluxUnit)) {
    fluxUnit <- fluxConst[fluxUnit][[1]]
  }
  unitFactorReturn <- fluxUnit@unitFactor
  
  dataStart <- min(eList$Sample$DecYear, na.rm = TRUE)
  dataStartPad <- dataStart - 0.5
  
  if(is.na(yearStart)){
    yearStart <- dataStartPad
  } else {
    yearStart <- max(yearStart, dataStartPad)
    
  }
  
  dataEnd <- max(eList$Sample$DecYear, na.rm = TRUE)
  dataEndPad <- dataEnd + 0.5
  
  if(is.na(yearEnd)){
    yearEnd <- dataEndPad
  } else {
    yearEnd <- min(yearEnd, dataEndPad)
  }
  
  if(is.na(fluxMax)){
    numYears <- length(localAnnualResults$DecYear)
    
    subAnnualResults <- localAnnualResults[localAnnualResults$DecYear>=yearStart & localAnnualResults$DecYear <= yearEnd,]
    
    annFlux <- unitFactorReturn*subAnnualResults$Flux
    
    if(plotGenFlux){
      fluxMax <- 1.05*max(c(CIAnnualResults$FNFluxHigh*unitFactorReturn,
                            annFlux,
                            unitFactorReturn*subAnnualResults$GenFlux), na.rm=TRUE)
    } else {
      fluxMax <- 1.05*max(c(CIAnnualResults$FNFluxHigh*unitFactorReturn,annFlux), na.rm=TRUE)
    }
    
  }
  
  EGRET::plotFluxHist(eList, yearStart = yearStart, yearEnd = yearEnd,
                      fluxUnit=fluxUnit, col.pred=col.pred,fluxMax=fluxMax,
                      printTitle=FALSE, cex.axis = cex.axis,
                      cex = cex, cex.main = cex.main, col.gen = col.gen,
                      plotGenFlux = plotGenFlux, plotAnnual = plotAnnual,
                      plotFlowNorm = plotFlowNorm, 
                      lwd = lwd, customPar = customPar,...)
  if (printTitle) {
    title(main=title, cex.main=cex.main)
  }
  
  if(!is.na(yearStart)){
    CIAnnualResults <- CIAnnualResults[CIAnnualResults$Year >= yearStart, ]
  }
  
  if(!is.na(yearEnd)){
    CIAnnualResults <- CIAnnualResults[CIAnnualResults$Year <= yearEnd, ]
  }
  
  lines(CIAnnualResults$Year, CIAnnualResults$FNFluxLow*unitFactorReturn,
        lty=2,col=col.pred)
  lines(CIAnnualResults$Year, CIAnnualResults$FNFluxHigh*unitFactorReturn, 
        lty=2,col=col.pred)
  
}

#' Single confidence interval bootstrap run
#'
#' One bootstrap run used in calculating confidence interval bands.
#'
#' @param eList named list with at least the Daily, Sample, and INFO dataframes. Created from the EGRET package, after running \code{\link[EGRET]{modelEstimation}}.
#' @param blockLength integer default value is 200.
#' @param startSeed setSeed value. Defaults to 494817. This is used to make repeatable output.
#' @param verbose logical specifying whether or not to display progress message.
#' @param jitterOn logical, if TRUE, adds "jitter" to the data in an attempt to avoid some numerical problems.  Default = FALSE.  See Details below.
#' @param V numeric a multiplier for addition of jitter to the data, default = 0.2.
#' @export
#' @details
#' In some situations numerical problems are encountered in the bootstrap process, resulting in highly unreasonable spikes in the confidence intervals.
#' The use of "jitter" can often prevent these problems, but should only be used when it is clearly needed.
#' It adds a small amount of random "jitter" to the explanatory variables of the WRTDS model.  The V parameter sets the scale of variation in the log discharge values.
#' The standard deviation of the added jitter is V * standard deviation of Log Q.
#' The default for V is 0.2.  Larger values should generally be avoided, and smaller values may be sufficient.
#' @examples
#' library(EGRET)
#' eList <- Choptank_eList
#' \dontrun{
#' annualResults <- bootAnnual(eList)
#' }
bootAnnual <- function(eList, blockLength = 200, startSeed = 494817,
                       verbose = FALSE, jitterOn = FALSE, V = 0.2){
  Sample <- eList$Sample
  Daily <- eList$Daily
  INFO <- eList$INFO
  
  if(is.null(INFO$edgeAdjust)){
    INFO$edgeAdjust <- FALSE
  }
  
  paStart <- 10
  paLong <- 12
  
  if(!is.null(INFO$paLong)){
    paLong <- INFO$paLong
  }  
  if(!is.null(INFO$paStart)){
    paStart <- INFO$paStart
  }
  
  bootSample <- blockSample(localSample = Sample, blockLength = blockLength, startSeed = startSeed)
  
  if(jitterOn) bootSample <- EGRET::jitterSam(bootSample, V = V)
  
  eListBoot <- EGRET::as.egret(INFO,Daily,bootSample,NA)
  
  if(isTRUE("runSeries" %in% names(attributes(eList)) && attr(eList, "runSeries"))){
    #Indicates runSeries was run
    
    seriesEList <- EGRET::runSeries(eList = eListBoot,
                                    windowSide = INFO$windowSide,
                                    surfaceStart = INFO$surfaceStart,
                                    surfaceEnd = INFO$surfaceEnd,
                                    flowBreak = INFO$flowBreak,
                                    Q1EndDate = INFO$Q1EndDate,
                                    QStartDate = INFO$QStartDate,
                                    QEndDate = INFO$QEndDate,
                                    wall = INFO$wall, 
                                    oldSurface = FALSE,
                                    sample1EndDate = INFO$sample1EndDate,
                                    sampleStartDate = INFO$sampleStartDate,
                                    sampleEndDate = INFO$sampleEndDate,
                                    paStart = INFO$paStart,
                                    paLong = INFO$paLong,
                                    minNumObs = INFO$minNumObs,
                                    minNumUncen = INFO$minNumUncen,
                                    windowY = INFO$windowY,
                                    windowQ = INFO$windowQ,
                                    windowS = INFO$windowS,
                                    edgeAdjust = INFO$edgeAdjust,
                                    verbose = verbose)
    Daily1 <- seriesEList$Daily
  } else {
    surfaces1 <- EGRET::estSurfaces(eListBoot, 
                                    windowY = eList$INFO$windowY, 
                                    windowQ = eList$INFO$windowQ, 
                                    windowS = eList$INFO$windowS,
                                    minNumObs = eList$INFO$minNumObs, 
                                    minNumUncen = eList$INFO$minNumUncen, 
                                    edgeAdjust = eListBoot$INFO$edgeAdjust,
                                    verbose = verbose)
    seriesEList <- EGRET::as.egret(INFO, Daily, bootSample, surfaces1)
    Daily1 <- EGRET::estDailyFromSurfaces(seriesEList)
  }
  
  annualResults1 <- EGRET::setupYears(Daily1, paStart=paStart, paLong=paLong)
  annualResults1$year <- as.integer(annualResults1$DecYear)
  annualResults <- annualResults1[,c("year","FNConc","FNFlux")]
  
  attr(annualResults, "blockLength") <- blockLength
  return(annualResults)
}

#' Confidence Interval Band Calculations
#'
#' @description 
#' Computes confidence intervals for Flow-Normalized Concentration 
#' and Flow-Normalized Flux for a WRTDS model.  
#'
#' @param repAnnualResults named list returned from bootstrapping process.
#' 
#' @param eList named list with at least the Daily, Sample, 
#' and INFO dataframes. Created from the EGRET package, after running 
#' \code{\link[EGRET]{modelEstimation}}.
#' 
#' @param probs numeric vector low and high confidence interval frequencies, 
#' default = c(0.05, 0.95) (which results in a 90\% confidence interval).
#' 
#' @export
#' 
#' @importFrom stats quantile
#' 
#' @examples
#' library(EGRET)
#' eList <- Choptank_eList
#' nBoot <- 100
#' blockLength <- 200
#' \dontrun{
#' 
#' repAnnualResults <- vector(mode = "list", length = nBoot)
#' for(n in 1:nBoot){
#'    annualResults <- bootAnnual(eList, blockLength, startSeed = n) 
#'    repAnnualResults[[n]] <- annualResults
#' }
#' 
#' CIAnnualResults <- ciBands(eList, repAnnualResults)
#' plotConcHistBoot(eList, CIAnnualResults)
#' }
#' 
ciBands <- function(eList, repAnnualResults, probs = c(0.05, 0.95)){
  
  if(length(probs) != 2){
    stop("Please provide only lower and upper limit in the probs argument")
  }
  
  paStart <- 10
  paLong <- 12
  
  INFO <- eList$INFO
  
  if(!is.null(INFO$paLong)){
    paLong <- INFO$paLong
  }
  
  if(!is.null(INFO$paStart)){
    paStart <- INFO$paStart
  }
  
  AnnualResults <- EGRET::setupYears(eList$Daily, paLong = paLong, paStart=paStart)
  AnnualResults$year <- as.integer(AnnualResults$DecYear)
  names(AnnualResults)[which(names(AnnualResults) %in% c("FNFlux"))] <- c("FNFlux_1")
  names(AnnualResults)[which(names(AnnualResults) %in% c("FNConc"))] <- c("FNConc_1")
  
  AnnualResults <- AnnualResults[c("year", "DecYear", "FNFlux_1", "FNConc_1")]
  nBoot <- length(repAnnualResults)
  numYears <- nrow(repAnnualResults[[1]])
  yearStart <- repAnnualResults[[1]][1,1]
  blockLength <- attr(repAnnualResults[[1]], "blockLength")
  
  manyAnnualResults <- array(NA, dim=c(numYears,2,nBoot))
  
  for (i in 1:nBoot){
    cat(i, "\n")
    df_1 <- repAnnualResults[[i]]
    df_1 <- merge(df_1, 
                  AnnualResults, by = "year", all = TRUE)
    
    manyAnnualResults[,1,i] <- 2*log(df_1$FNConc_1) - log(df_1$FNConc)
    manyAnnualResults[,2,i] <- 2*log(df_1$FNFlux_1) - log(df_1$FNFlux)
  }
  
  CIAnnualResults <- data.frame(matrix(ncol = 5, nrow = numYears))
  names(CIAnnualResults) <- c("Year","FNConcLow","FNConcHigh","FNFluxLow","FNFluxHigh")
  
  for(iYear in 1:numYears) {
    quantConc <- quantile(manyAnnualResults[iYear,1,1:nBoot],prob=probs,type=6,na.rm = TRUE)
    quantFlux <- quantile(manyAnnualResults[iYear,2,1:nBoot],prob=probs,type=6,na.rm = TRUE)
    
    CIAnnualResults$Year[iYear] <- AnnualResults$DecYear[iYear]
    CIAnnualResults$FNConcLow[iYear] <- exp(quantConc[1])
    CIAnnualResults$FNConcHigh[iYear] <- exp(quantConc[2])
    CIAnnualResults$FNFluxLow[iYear] <- exp(quantFlux[1])
    CIAnnualResults$FNFluxHigh[iYear] <- exp(quantFlux[2])
  }
  
  attr(CIAnnualResults, "nBoot") <- nBoot
  attr(CIAnnualResults, "probs") <- probs
  attr(CIAnnualResults, "blockLength") <- blockLength
  
  return(CIAnnualResults)
}

#' plotHistogramTrend
#'
#' Produces a histogram of trend results from bootstrap process.  The histogram shows the trend results expressed as percentage change between the first year (or first period) 
#' and the second year (or second period).  It shows the zero line (no trend) and also shows the WRTDS 
#' estimate of the trend in percent.  It is based on the output of either wBT or 
#' runPairsBoot.
#'
#' @param eList named list with at least the Daily, Sample, and INFO dataframes. Created from the EGRET package, after running \code{\link[EGRET]{modelEstimation}}.
#' @param eBoot named list. Returned from \code{\link{wBT}} or from \code{\link{runPairsBoot}}.
#' @param caseSetUp data frame. Returned from \code{\link{trendSetUp}}, or if \code{\link{runPairsBoot}} was used, need to specify caseSetUp = NA.
#' @param flux logical if TRUE, plots flux results, if FALSE plots concentration results.
#' @param xMin minimum bin value for histogram, it is good to have the xMin and xMax arguments straddle zero, default is NA (value set from the data). 
#' @param xMax maximum bin value for histogram, default is NA (value set from the data).
#' @param xStep step size, typically multiples of 10 or 20, default is NA (value set from the data).
#' @param printTitle logical if TRUE, plot includes title.
#' @param cex.main numeric magnification of font size for title, default is 1.1.
#' @param cex.axis numeric magnification of font size for axis, default is 1.1.
#' @param cex.lab numeric magnification of font size for axis labels, default is 1.1. 
#' @param col.fill character fill color for histogram, default is "grey".
#' @param \dots base R graphical parameters that can be passed to the hist function
#' @export
#' @details
#' For any given set of results (from eBoot) it is best to run it first with the arguments
#' xMin = NA, xMax = NA, and xStep = NA.  Then, observing the range the histogram covers
#' it can be run again with values of these three arguments selected by the user to provide
#' for a more readable version of the histogram.
#' @importFrom graphics hist
#' @importFrom graphics abline
#' @importFrom graphics box
#' @importFrom graphics axis
#' @examples
#' library(EGRET)
#' eList <- Choptank_eList
#' eBoot <- Choptank_eBoot
#' caseSetUp <- Choptank_caseSetUp
#' plotHistogramTrend(eList, eBoot, caseSetUp, flux = FALSE)
#' 
#' \dontrun{
#' # Using wBT:	
#' caseSetUp <- trendSetUp(eList)
#' eBoot <- wBT(eList,caseSetUp)
#' plotHistogramTrend(eList, eBoot, caseSetUp,  
#'                    flux = FALSE, xMin = -20, xMax = 60, xStep = 5)
#' plotHistogramTrend(eList, eBoot, caseSetUp, 
#'                    flux = TRUE, xMin = -20, xMax = 60, xStep = 5)
#'    
#' # Using runPairs followed by runPairsBoot:
#' year1 <- 1985
#' year2 <- 2009          
#' pairOut_2 <- runPairs(eList, year1, year2, windowSide = 7)
#' boot_pair_out <- runPairsBoot(eList, pairOut_2, nBoot = 10)
#' 
#' plotHistogramTrend(eList, boot_pair_out, caseSetUp = NA, 
#'                    flux = TRUE, xMin = -20, xMax = 60, xStep = 5)          
#' }
plotHistogramTrend <- function (eList, eBoot, caseSetUp, 
                                flux = TRUE, xMin = NA, xMax = NA, xStep = NA, 
                                printTitle=TRUE, cex.main=1.1, cex.axis = 1.1, cex.lab = 1.1, col.fill="grey",...){
  
  
  if(all(c("paStart", "paLong") %in% names(attributes(eBoot)))){
    paStart <- attr(eBoot, "paStart")
    paLong <- attr(eBoot, "paLong")
  } else {
    paStart <- eList$INFO$paStart
    paLong <- eList$INFO$paLong
  }
  
  periodName <- EGRET::setSeasonLabel(data.frame(PeriodStart = paStart, 
                                                 PeriodLong = paLong))
  
  if(any(c("yearPair","group1firstYear") %in% names(attributes(eBoot))) |  "segmentInfo" %in% names(attributes(eList$INFO))){
    periodName <- paste(periodName, "*")
  }
  
  if (flux) {
    change <- 100 * eBoot$bootOut$estF/eBoot$bootOut$baseFlux
    reps <- eBoot$pFlux
    xlabel <- "Flux trend, in %"
    titleWord <- "Flux"
  } else {
    change <- 100 * eBoot$bootOut$estC/eBoot$bootOut$baseConc
    reps <- eBoot$pConc
    xlabel <- "Concentration trend, in %"
    titleWord <- "Concentration"
  }
  
  if(!("group2firstYear" %in% names(attributes(eBoot)))){
    if(all(is.na(caseSetUp))){
      if("year1" %in% names(attributes(eBoot))){
        year1 <- attr(eBoot, "year1")
      }
      if("year2" %in% names(attributes(eBoot))){
        year2 <- attr(eBoot, "year2")
      }    
      
    } else {
      year1 <- caseSetUp$year1
      year2 <- caseSetUp$year2
    }
    
    if(any(is.na(c(year1,year2)))){
      stop("Provide caseSetUp information")
    }
    
    titleToPrint <- ifelse(printTitle, paste("Trend magnitude in", 
                                             eList$INFO$paramShortName, "\nFlow Normalized", titleWord, 
                                             year1, "to", year2, "\n", eList$INFO$shortName, 
                                             periodName), "")
    
  } else {
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
    
  }
  
  minReps <- min(reps, na.rm = TRUE)
  maxReps <- max(reps, na.rm = TRUE)
  xMin <- ifelse(is.na(xMin), min(-10, minReps),xMin)
  xMax <- ifelse(is.na(xMax), max(10, maxReps), xMax)
  xStep <- ifelse(is.na(xStep), (xMax - xMin)/10, xStep)
  xSeq <- seq(xMin, xMax, xStep)
  yLim <- c(0,1.04*max(hist(reps, breaks = xSeq, plot = FALSE)$density, na.rm = TRUE))
  
  hist(reps, breaks = xSeq, axes = FALSE, ylab = "", yaxs = "i", xaxs = "i", 
       main = titleToPrint, freq = FALSE, xlab = xlabel, col = col.fill, 
       cex.main = cex.main, cex.lab = cex.lab, ylim = yLim, ...)
  abline(v = change, lwd = 3, lty = 2)
  abline(v = 0, lwd = 3)
  box()
  axis(1, tcl = 0.5, labels = TRUE, cex.axis = cex.axis)
  axis(2, tcl = 0.5, labels = TRUE, las = 1, cex.axis = cex.axis)
  title(ylab = "Density", line = 4.5, cex.lab = cex.lab)
  axis(3, tcl = 0.5, labels = FALSE)
  axis(4, tcl = 0.5, labels = FALSE)
  
}

#' ciCalculations
#'
#' Function to calculate confidence bands for flow normalized concentration or 
#' flow normalized flux.It returns the data frame CIAnnualResults, which is used
#' as input to the functions \code{plotConcHistBoot}, and
#' \code{plotFluxHistBoot} which produce the graphical output. 
#'
#' @param eList named list with at least the Daily, Sample, and INFO dataframes. Created from the EGRET package, after running \code{\link[EGRET]{modelEstimation}}.
#' @param startSeed setSeed value. Defaults to 494817. This is used to make repeatable output.
#' @param verbose logical specifying whether or not to display progress message, default = TRUE
#' @param jitterOn logical, if TRUE, adds "jitter" to the data in an attempt to avoid some numerical problems.  Default = FALSE.  See Details below.
#' @param V numeric a multiplier for addition of jitter to the data, default = 0.2.  See Details below.  
#' @param nBoot number of times the bootstrap resampling and model estimating is done.
#' Default is 100, but that will take a long time. Testing should initially be done using
#' a smaller number like 10.
#' @param blockLength integer size of subset, expressed in days.  200 days has 
#' been found to be a good choice.
#' @param widthCI numeric, the width of the confidence intervals. 0.9 means the 
#' confidence intervals will be calculated with 90\%.
#' 
#' @export
#' 
#' @return CIAnnualResults a data frame with the following columns:
#' \tabular{ll}{
#' Year \tab mean decYear value for the year being reported \cr
#' FNConcLow \tab the lower confidence limit for flow normalized concentration, in mg/L \cr
#' FNConcHigh \tab  the upper confidence limit for flow normalized concentration, in mg/L \cr
#' FNFluxLow \tab  the lower confidence limit for flow normalized flux, in kg/day \cr
#' FNFluxLow \tab  the lower confidence limit for flow normalized flux, in kg/day \cr
#' }
#' 
#' @details
#' In some situations numerical problems are encountered in the bootstrap process, resulting in highly unreasonable spikes in the confidence intervals.
#' The use of "jitter" can often prevent these problems, but should only be used when it is clearly needed.
#' It adds a small amount of random "jitter" to the explanatory variables of the WRTDS model.  The V parameter sets the scale of variation in the log discharge values.
#' The standard deviation of the added jitter is V * standard deviation of Log Q.
#' The default for V is 0.2.  Larger values should generally be avoided, and smaller values may be sufficient.
#'
#' Argument values suggested.  
#' To test the code nBoot = 10 is sufficient, but for meaningful results nBoot = 100 or even nBoot = 500 are more appropriate.
#' blockLength = 200
#' widthCI = 90 (90\% confidence interval)
#' 
#' @examples
#' library(EGRET)
#' eList <- Choptank_eList
#' \dontrun{
#' CIAnnualResults <- ciCalculations(eList,
#'                                   nBoot = 10)
#' plotConcHistBoot(eList, CIAnnualResults)
#' 
#' # run in batch mode, using non-stationary flow normalization
#' # In this example nBoot is set very small, useful for an initial trial run.
#' # A meaningful application would use nBoot values such as 100 or even 500. 
#' seriesOut_2 <- runSeries(eList, windowSide = 11)
#' CIAnnualResults <- ciCalculations(seriesOut_2, 
#'                      nBoot = 10,
#'                      blockLength = 200,
#'                      widthCI = 90)
#'                      
#'  plotConcHistBoot(seriesOut_2, CIAnnualResults)
#' 
#' }
ciCalculations <- function(eList, 
                            startSeed = 494817,
                            verbose = TRUE,
                            jitterOn = FALSE, 
                            V = 0.2,
                            nBoot = 100,
                            blockLength = 200,
                            widthCI = 90){
  
  INFO <- eList$INFO

  ciLower <- (50-(widthCI/2))/100
  ciUpper <- (50+(widthCI/2))/100
  probs <- c(ciLower,ciUpper)
  
  repAnnualResults <- vector(mode = "list", length = nBoot)
  
  if(isTRUE("runSeries" %in% names(attributes(eList)) && attr(eList, "runSeries"))){
    #Indicates runSeries was run
    cat("\nRunning the EGRET runSeries function to have that as a baseline for the Confidence Bands\n")
    
    eList <- EGRET::runSeries(eList = eList,
                              windowSide = INFO$windowSide,
                              surfaceStart = INFO$surfaceStart,
                              surfaceEnd = INFO$surfaceEnd,
                              flowBreak = INFO$flowBreak,
                              Q1EndDate = INFO$Q1EndDate,
                              QStartDate = INFO$QStartDate,
                              QEndDate = INFO$QEndDate,
                              wall = INFO$wall,
                              oldSurface = TRUE,
                              sample1EndDate = INFO$sample1EndDate,
                              sampleStartDate = INFO$sampleStartDate,
                              sampleEndDate = INFO$sampleEndDate,
                              paStart = INFO$paStart,
                              paLong = INFO$paLong,
                              minNumObs = INFO$minNumObs,
                              minNumUncen = INFO$minNumUncen,
                              windowY = INFO$windowY,
                              windowQ = INFO$windowQ,
                              windowS = INFO$windowS,
                              edgeAdjust = INFO$edgeAdjust, verbose = verbose)
  } else {
    cat("\nRunning the EGRET modelEstimation function first to have that as a baseline for the Confidence Bands")
    
    eList <- EGRET::modelEstimation(eList, windowY = eList$INFO$windowY,
                                    windowQ = eList$INFO$windowQ,
                                    windowS = eList$INFO$windowS,
                                    minNumObs = eList$INFO$minNumObs,
                                    minNumUncen = eList$INFO$minNumUncen,
                                    verbose = verbose)
    
  }
  
  for(n in 1:nBoot){
    cat(n, "\n")
    repAnnualResults[[n]] <- bootAnnual(eList, blockLength, startSeed+n, verbose = verbose, 
                                        jitterOn = jitterOn, V = V)
  }
  
  CIAnnualResults <- ciBands(eList, repAnnualResults, probs)
  
  return(CIAnnualResults)
  
}
