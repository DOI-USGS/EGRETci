#' Graph of annual concentration, flow normalized concentration, 
#' and confidence bands for flow normalized concentrations
#' 
#' Uses the output of \code{\link[EGRET]{modelEstimation}} in the EGRET package (results in the named 
#' list eList), and the data frame CIAnnualResults (produced by EGRETci package 
#' using scripts described in the vignette) to produce a graph of annual 
#' concentration, flow normalized concentration, and confidence bands for 
#' flow-normalized concentrations.  In addition to the arguments listed below, 
#' it will accept any additional arguments that are listed for the EGRET function 
#' \code{\link[EGRET]{plotConcHist}}.
#'
#' @param eList named list with at least the Daily, Sample, and INFO dataframes. Created from the EGRET package, after running \code{\link[EGRET]{modelEstimation}}.
#' @param CIAnnualResults data frame generated from ciBands (includes nBoot, probs, and blockLength attributes)
#' @param yearStart numeric is the calendar year containing the first estimated annual value to be plotted, default is NA (which allows it to be set automatically by the data)
#' @param yearEnd numeric is the calendar year just after the last estimated annual value to be plotted, default is NA (which allows it to be set automatically by the data)
#' @param plotFlowNorm logical variable if TRUE flow normalized line is plotted, if FALSE not plotted 
#' @param col.pred character prediction color
#' @param concMax number specifying the maximum value to be used on the vertical axis, default is NA (which allows it to be set automatically by the data)
#' @param printTitle logical
#' @param cex.main numeric title scale
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
#' \dontrun{
#' CIAnnualResults <- ciCalculations(eList, nBoot = 100, blockLength = 200)
#' plotConcHistBoot(eList, CIAnnualResults)
#' }
plotConcHistBoot <- function (eList, CIAnnualResults, yearStart = NA, yearEnd = NA, 
                              plotFlowNorm=TRUE, col.pred="green", concMax = NA,
                              printTitle=TRUE, cex.main=1.1, ...){
  
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
  
  title3 <- paste(widthCI,"% CI on FN Concentration, Replicates =",nBoot,"Block=",blockLength,"days")

  title <- paste(eList$INFO$shortName, " ", eList$INFO$paramShortName, 
                   "\n", periodName, "\n",title3)
  
  if(is.na(concMax)){
    numYears <- length(localAnnualResults$DecYear)

    subAnnualResults <- localAnnualResults[localAnnualResults$DecYear>=yearStart & localAnnualResults$DecYear <= yearEnd,]
    
    if(is.na(yearStart)){
      yearStart <- min(localAnnualResults$DecYear[!is.na(localAnnualResults$FNConc)], na.rm = TRUE)
    }
    
    if(is.na(yearEnd)){
      yearEnd <- max(localAnnualResults$DecYear[!is.na(localAnnualResults$FNConc)], na.rm = TRUE)
    }
    
    annConc <- subAnnualResults$Conc
    concMax <- 1.05*max(c(CIAnnualResults$FNConcHigh,annConc), na.rm=TRUE)
  }
  
  EGRET::plotConcHist(eList, yearStart = yearStart, yearEnd = yearEnd,
               col.pred=col.pred, printTitle=FALSE, 
               plotFlowNorm = plotFlowNorm, concMax = concMax, ...)
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
#' @param CIAnnualResults data frame from ciBands (needs nBoot, probs, and blockLength attributes)
#' @param yearStart numeric is the calendar year containing the first estimated annual value to be plotted, default is NA (which allows it to be set automatically by the data)
#' @param yearEnd numeric is the calendar year just after the last estimated annual value to be plotted, default is NA (which allows it to be set automatically by the data)
#' @param fluxUnit number representing entry in pre-defined fluxUnit class array. \code{\link{printFluxUnitCheatSheet}}
#' @param fluxMax number specifying the maximum value to be used on the vertical axis, default is NA (which allows it to be set automatically by the data)
#' @param plotFlowNorm logical variable if TRUE flow normalized line is plotted, if FALSE not plotted 
#' @param col.pred character prediction color
#' @param printTitle logical
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
                              yearStart=NA, yearEnd=NA,
                              plotFlowNorm=TRUE, fluxUnit = 9, fluxMax=NA,
                              col.pred="green", printTitle=TRUE, 
                              cex.main=1.1, ...){
  
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
  
  title3 <- paste(widthCI,"% CI on FN Flux, Replicates =",nBoot,", Block=",blockLength,"days")
  
  title <- paste(eList$INFO$shortName, " ", eList$INFO$paramShortName, 
                 "\n", periodName, "\n",title3)
  
  if (is.numeric(fluxUnit)) {
    fluxUnit <- fluxConst[shortCode = fluxUnit][[1]]
  } else if (is.character(fluxUnit)) {
    fluxUnit <- fluxConst[fluxUnit][[1]]
  }
  unitFactorReturn <- fluxUnit@unitFactor
  
  if(is.na(fluxMax)){
    numYears <- length(localAnnualResults$DecYear)

    yearStart <- if(is.na(yearStart)) trunc(min(localAnnualResults$DecYear[!is.na(localAnnualResults$FNFlux)],na.rm = TRUE)) else yearStart
    yearEnd <- if(is.na(yearEnd)) trunc(max(localAnnualResults$DecYear[!is.na(localAnnualResults$FNFlux)],na.rm = TRUE))+1 else yearEnd
    
    subAnnualResults <- localAnnualResults[localAnnualResults$DecYear>=yearStart & localAnnualResults$DecYear <= yearEnd,]
    
    annFlux <- unitFactorReturn*subAnnualResults$Flux
    
    fluxMax <- 1.05*max(c(CIAnnualResults$FNFluxHigh*unitFactorReturn,annFlux), na.rm=TRUE)
  }
  
  EGRET::plotFluxHist(eList, yearStart = yearStart, yearEnd = yearEnd,
               fluxUnit=fluxUnit, col.pred=col.pred,fluxMax=fluxMax,
               plotFlowNorm = plotFlowNorm, printTitle=FALSE,...)
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

#' bootAnnual
#'
#' One bootstrap run.
#'
#' @param eList named list with at least the Daily, Sample, and INFO dataframes. Created from the EGRET package, after running \code{\link[EGRET]{modelEstimation}}.
#' @param blockLength integer suggested value is 200
#' @param startSeed setSeed value. Defaults to 494817. This is used to make repeatable output.
#' @param verbose logical specifying whether or not to display progress message
#' @export
#' @examples
#' library(EGRET)
#' eList <- Choptank_eList
#' \dontrun{
#' annualResults <- bootAnnual(eList)
#' }
bootAnnual <- function(eList, blockLength=200, startSeed = 494817, verbose = FALSE){
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

#' ciBands
#'
#' Computes confidence intervals for Flow-Normalized Concentration 
#' and Flow-Normalized Flux for a WRTDS model.  
#'
#' @param repAnnualResults named list returned from bootstrapping process
#' @param eList named list with at least the Daily, Sample, and INFO dataframes. Created from the EGRET package, after running \code{\link[EGRET]{modelEstimation}}.
#' @param probs vector high and low confidence interval percentages
#' @export
#' @importFrom stats quantile
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
#' 
#' }
ciBands <- function(eList, repAnnualResults, probs=c(0.05,0.95)){

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
  
  nBoot <- length(repAnnualResults)
  numYears <- nrow(repAnnualResults[[1]])
  yearStart <- repAnnualResults[[1]][1,1]
  blockLength <- attr(repAnnualResults[[1]], "blockLength")
  
  manyAnnualResults <- array(NA, dim=c(numYears,2,nBoot))
  for (i in 1:nBoot){
    manyAnnualResults[,1,i] <- 2*log(AnnualResults$FNConc) - log(repAnnualResults[[i]]$FNConc)
    manyAnnualResults[,2,i] <- 2*log(AnnualResults$FNFlux) - log(repAnnualResults[[i]]$FNFlux)
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
#' Histogram of trend.
#'
#' @param eList named list with at least the Daily, Sample, and INFO dataframes. Created from the EGRET package, after running \code{\link[EGRET]{modelEstimation}}.
#' @param eBoot named list. Returned from \code{\link{wBT}}.
#' @param caseSetUp data frame. Returned from \code{\link{trendSetUp}}.
#' @param flux logical if TRUE, plots flux results, if FALSE plots concentration
#' @param xMin minimum bin value, it is good to have the xMin and xMax arguments straddle zero. 
#' @param xMax maximum bin value
#' @param xStep step size, should probably be multiples of 10 or 20
#' @param printTitle logical if TRUE, includes title
#' @param cex.main numeric title font size
#' @param cex.axis numeric axis font size
#' @param cex.lab numeric label font size
#' @param col.fill character fill color
#' @param \dots base R graphical parameters that can be passed to the hist function
#' @export
#' @importFrom graphics hist
#' @importFrom graphics abline
#' @importFrom graphics box
#' @importFrom graphics axis
#' @examples
#' library(EGRET)
#' eList <- Choptank_eList
#' eBoot <- Choptank_eBoot
#' caseSetUp <- Choptank_caseSetUp
#' plotHistogramTrend(eList, eBoot, caseSetUp, flux=FALSE)
#' 
#' \dontrun{
#' caseSetUp <- trendSetUp(eList)
#' eBoot <- wBT(eList,caseSetUp)
#' plotHistogramTrend(eList, eBoot, caseSetUp,  
#'                    flux=FALSE, xMin = -20, xMax = 60, xStep = 5)
#' plotHistogramTrend(eList, eBoot, caseSetUp, 
#'                    flux=TRUE, xMin = -20, xMax = 60, xStep = 5)
#'    
#' # Using runPairs:
#' year1 <- 1985
#' year2 <- 2009          
#' pairOut_2 <- runPairs(eList, year1, year2, windowSide = 7)
#' boot_pair_out <- runPairsBoot(eList, pairOut_2, nBoot = 10)
#' 
#' plotHistogramTrend(eList, boot_pair_out,caseSetUp=NA, 
#'                    flux=TRUE, xMin = -20, xMax = 60, xStep = 5)          
#' }
plotHistogramTrend <- function (eList, eBoot, caseSetUp, 
                                flux = TRUE, xMin = NA, xMax = NA, xStep = NA, 
                                printTitle=TRUE, cex.main=1.1, cex.axis = 1.1, cex.lab = 1.1, col.fill="grey",...){
  
  periodName <- EGRET::setSeasonLabel(data.frame(PeriodStart = eList$INFO$paStart, 
                                          PeriodLong = eList$INFO$paLong))
  
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
  
  hist(reps, breaks = xSeq, yaxs = "i", xaxs = "i", axes = FALSE, ylab = "",
       main = titleToPrint, freq = FALSE, xlab = xlabel, col = col.fill, 
       cex.main = cex.main, cex.lab = cex.lab, ...)
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
#' Interactive function to calculate WRTDS confidence bands
#'
#' @param eList named list with at least the Daily, Sample, and INFO dataframes. Created from the EGRET package, after running \code{\link[EGRET]{modelEstimation}}.
#' @param startSeed setSeed value. Defaults to 494817. This is used to make repeatable output.
#' @param verbose logical specifying whether or not to display progress message
#' @param \dots optionally include nBoot, blockLength, or widthCI
#' @export
#' @examples
#' library(EGRET)
#' eList <- Choptank_eList
#' \dontrun{
#' CIAnnualResults <- ciCalculations(eList)
#' 
#' seriesOut_2 <- runSeries(eList, windowSide = 7)
#' CIAnnualResults <- ciCalculations(seriesOut_2, 
#'                      nBoot = 10,
#'                      blockLength = 200,
#'                      widthCI = 90)
#'                      
#'  plotConcHistBoot(seriesOut_2, CIAnnualResults)
#' 
#' }
ciCalculations <- function (eList, 
                            startSeed = 494817,
                            verbose = TRUE,
                            ...){
  
  matchReturn <- list(...)
  
  INFO <- eList$INFO
  
  if(!is.null(matchReturn$nBoot)){
    nBoot <- matchReturn$nBoot
  } else {
    message("Enter nBoot, the number of bootstrap replicates to be used, typically 100")
    nBoot <- as.numeric(readline())
    cat("nBoot = ",nBoot," this is the number of replicates that will be run\n")
  }
  
  if(!is.null(matchReturn$blockLength)){
    blockLength <- matchReturn$blockLength
  } else {
    message("Enter blockLength, in days, typically 200 is a good choice")
    blockLength <- as.numeric(readline())
  }
  
  if(!is.null(matchReturn$widthCI)){
    widthCI <- matchReturn$widthCI
  } else {
    message("Enter confidence interval, for example 90 represents a 90% confidence interval,")
    message("the low and high returns are 5 and 95 % respectively")
    widthCI <-  as.numeric(readline())
  }
  
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
    repAnnualResults[[n]] <- bootAnnual(eList, blockLength, startSeed+n, verbose = verbose)
  }
  
  CIAnnualResults <- ciBands(eList, repAnnualResults, probs)
  
  return(CIAnnualResults)

}
