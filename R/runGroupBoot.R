
#' The bootstrap uncertainty analysis for runGroups results
#' 
#' This function that does the uncertainty analysis for determining the change 
#' between two groups of years.  The process is virtually 
#' identical to what is used for \code{\link{runPairsBoot}}.  
#' 
#' @param eList named list with at least the Daily, Sample, and INFO dataframes
#' @param groupResults data frame returned from \code{\link[EGRET]{runGroups}}
#' @param nBoot the maximum number of bootstrap replicates to be used, typically 100
#' @param blockLength days, typically 200 is a good choice
#' @param startSeed setSeed value. Defaults to 494817. This is used to make repeatable output.
#' @export
#' @return eBoot, a named list with bootOut,wordsOut,xConc,xFlux values. bootOut is a data frame with the results
#' of the bootstrapping tests. wordsOut is a character vector describing the results.
#' xConc, xFlux are vectors of length iBoot, of the change in flow normalized concentration or flux 
#' computed by each bootstrap replicate (mg/L). pConc and pFlux are vectors of length iBoot, of the change 
#' in flow normalized concentration or flux computed from each bootstrap replicate expressed as % change. 
#' @seealso \code{\link{runPairsBoot}}, \code{\link[EGRET]{runGroups}}
#' @examples 
#' library(EGRET)
#' eList <- Choptank_eList
#' 
#' \dontrun{
#' groupResults <- runGroups(eList, 
#'                           group1firstYear = 1995, 
#'                           group1lastYear = 2004, 
#'                           group2firstYear = 2005, 
#'                           group2lastYear = 2014, 
#'                           windowSide = 7, wall = TRUE, 
#'                           sample1EndDate = "2004-10-30", 
#'                           paStart = 4, paLong = 2, 
#'                           verbose = FALSE)
#' 
#' boot_group_out <- runGroupsBoot(eList, groupResults)
#' 
#' plotHistogramTrend(eList, boot_group_out, caseSetUp=NA)
#' }
runGroupsBoot <- function (eList, groupResults, nBoot = 100, 
                  startSeed = 494817, blockLength = 200){
  interactive <- FALSE
  localINFO <- eList$INFO
  localDaily <- eList$Daily
  localSample <- eList$Sample
  firstDayDaily <- min(localDaily$Date, na.rm = TRUE)
  lastDayDaily <- max(localDaily$Date, na.rm = TRUE)
  firstDaySample <- min(localSample$Date, na.rm = TRUE)
  lastDaySample <- max(localSample$Date, na.rm = TRUE)
  prob = c(0.025, 0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95, 0.975)
  
  words <- function(z) {
    out <- if (z) 
      "Reject Ho"
    else "Do Not Reject Ho"
    return(out)
  }
  bootOut <- as.data.frame(matrix(ncol = 27, nrow = 1))
  
  colnames(bootOut) <- c("rejectC", "pValC", "estC", "lowC90", 
                         "upC90", "lowC50", "upC50", "lowC95", "upC95", "likeCUp", 
                         "likeCDown", "rejectF", "pValF", "estF", "lowF90", "upF90", 
                         "lowF50", "upF50", "lowF95", "upF95", "likeFUp", "likeFDown", 
                         "baseConc", "baseFlux", "nBoot", "startSeed", "blockLength")
  
  paStart <- attr(groupResults, "groupInfo")[["paStart"]]
  paLong <- attr(groupResults, "groupInfo")[["paLong"]]
  group1firstYear <- attr(groupResults, "groupInfo")[["group1firstYear"]]
  group1lastYear <- attr(groupResults, "groupInfo")[["group1lastYear"]]
  group2firstYear <- attr(groupResults, "groupInfo")[["group2firstYear"]]
  group2lastYear <- attr(groupResults, "groupInfo")[["group2lastYear"]]
  sample1StartDate <- attr(groupResults, "SampleBlocks")[["sample1StartDate"]]
  sample1EndDate <- attr(groupResults, "SampleBlocks")[["sample1EndDate"]]
  sample2StartDate <- attr(groupResults, "SampleBlocks")[["sample2StartDate"]]
  sample2EndDate <- attr(groupResults, "SampleBlocks")[["sample1EndDate"]]
  sampleStartDate <- sample1StartDate
  sampleEndDate <- sample2EndDate
  dateInfo <- attr(groupResults, "dateInfo")
  sample1StartDate <- attr(groupResults, "SampleBlocks")[["sample1StartDate"]]
  sample1EndDate <- attr(groupResults, "SampleBlocks")[["sample1EndDate"]]
  sample2StartDate <- attr(groupResults, "SampleBlocks")[["sample2StartDate"]]
  sample2EndDate <- attr(groupResults, "SampleBlocks")[["sample2EndDate"]]
  surfaceStart <- attr(groupResults, "SampleBlocks")[["surfaceStart"]]
  surfaceEnd <- attr(groupResults, "SampleBlocks")[["surfaceEnd"]]
  minNumObs <- attr(groupResults, "Other")[["minNumObs"]]
  minNumUncen <- attr(groupResults, "Other")[["minNumUncen"]]
  windowY <- attr(groupResults, "Other")[["windowY"]]
  windowQ <- attr(groupResults, "Other")[["windowQ"]]
  windowS <- attr(groupResults, "Other")[["windowS"]]
  wall <- attr(groupResults, "Other")[["wall"]]
  edgeAdjust <- attr(groupResults, "Other")[["edgeAdjust"]]
  
  xConc <- rep(NA, nBoot)
  xFlux <- rep(NA, nBoot)
  pConc <- rep(NA, nBoot)
  pFlux <- rep(NA, nBoot)
  
  regDeltaConc <- groupResults$x22[1] - groupResults$x11[1]
  estC <- regDeltaConc
  baseConc <- groupResults$x11[1]
  regDeltaConcPct <- (regDeltaConc/baseConc) * 100
  LConcDiff <- log(groupResults$x22[1]) - log(groupResults$x11[1])
  regDeltaFlux <- (groupResults$x22[2] - groupResults$x11[2])
  estF <- regDeltaFlux
  baseFlux <- groupResults$x11[2]
  regDeltaFluxPct <- (regDeltaFlux/baseFlux) * 100
  LFluxDiff <- log(groupResults$x22[2]) - log(groupResults$x11[2])
  fcc <- format(regDeltaConc, digits = 3, width = 7)
  ffc <- format(regDeltaFlux, digits = 3, width = 8)
                       
  nBootGood <- 0
  
  for (iBoot in 1:(2 * nBoot)) {
    bootSample <- blockSample(localSample = localSample, 
                              blockLength = blockLength, startSeed = startSeed + iBoot)
    eListBoot <- suppressMessages(EGRET::as.egret(localINFO, localDaily, bootSample,NA))

    if(wall) {
      possibleError <- tryCatch(surfaces <- suppressMessages(EGRET::stitch(eListBoot, surfaceStart = surfaceStart, surfaceEnd = surfaceEnd,
                         sample1StartDate = sample1StartDate, sample1EndDate = sample1EndDate,
                         sample2StartDate = sample2StartDate, sample2EndDate = sample2EndDate,
                         windowY = windowY, windowQ = windowQ, windowS = windowS,
                         minNumObs = minNumObs, minNumUncen = minNumUncen, edgeAdjust = edgeAdjust)), 
                         error = function(e) e)
    } else {
      possibleError <- tryCatch(surfaces <- EGRET::estSurfaces(eListBoot, surfaceStart = surfaceStart, surfaceEnd = surfaceEnd,
                              windowY = windowY, windowQ = windowQ, windowS = windowS,
                              minNumObs = minNumObs, minNumUncen = minNumUncen, edgeAdjust = edgeAdjust),
                              error = function(e) e)
    }
    if (!inherits(possibleError, "error") ) {
      eListS <- suppressMessages(EGRET::as.egret(eListBoot$INFO, eListBoot$Daily, eListBoot$Sample, surfaces))
      eListOut <- suppressMessages(EGRET::flexFN(eListS, dateInfo, flowNormStartCol = "flowNormStart", 
                         flowNormEndCol = "flowNormEnd", flowStartCol = "flowStart", 
                         flowEndCol = "flowEnd"))
      eListOut$INFO$wall <- wall
      eListOut$INFO$surfaceStart <- surfaceStart
      eListOut$INFO$surfaceEnd <- surfaceEnd
      DailyFlex <- eListOut$Daily
      annFlex <- EGRET::setupYears(DailyFlex, paLong = paLong, paStart = paStart)
      annFlex$year <- floor(annFlex$DecYear + (annFlex$PeriodLong / 12) * 0.5)
      annFlex1 <- annFlex[annFlex$DecYear >= group1firstYear & annFlex$DecYear <= group1lastYear,]
      annFlex2 <- annFlex[annFlex$DecYear >= group2firstYear & annFlex$DecYear <= group2lastYear,]
      
      #  pairResults are in 10^6 kg/year, when we get to the bootstrap results
      #  Converting them all to 10^6 kg/year units
      c11 <- mean(annFlex1$FNConc, na.rm = TRUE)
      f11 <- mean(annFlex1$FNFlux, na.rm = TRUE) * 0.00036525
      c22 <- mean(annFlex2$FNConc, na.rm = TRUE)
      f22 <- mean(annFlex2$FNFlux, na.rm = TRUE) * 0.00036525
      xConc_here <- (2 * regDeltaConc) - (c22 - c11)
      xFlux_here <- (2 * regDeltaFlux) - (f22 - f11)
      if (!is.na(xConc_here) & !is.na(xFlux_here)) {
        nBootGood <- nBootGood + 1
        xConc[nBootGood] <- xConc_here
        xFlux[nBootGood] <- xFlux_here
        LConc <- (2 * LConcDiff) - (log(c22) - log(c11))
        pConc[nBootGood] <- (100 * exp(LConc)) - 100
        LFlux <- (2 * LFluxDiff) - (log(f22) - log(f11))
        pFlux[nBootGood] <- (100 * exp(LFlux)) - 100
        cat("\n iBoot, xConc and xFlux", nBootGood, xConc[nBootGood], 
            xFlux[nBootGood])
        if (nBootGood >= nBoot) {
          (break)()
        }
      }
    } else {
      stop(possibleError3, "\n", possibleError4)
    }
  }

  if (iBoot == 2 * nBoot) {
    message(iBoot, " iterations were run. They only achieved ", 
            nBootGood, " sucessful runs.")
  } else if (iBoot > nBoot) {
    message("It took ", iBoot, " iterations to achieve ", 
            nBoot, " sucessful runs.")
  }
  quantConc <- quantile(xConc, prob, type = 6, na.rm = TRUE)
  lowConc <- quantConc[["5%"]]
  highConc <- quantConc[["95%"]]
  quantFlux <- quantile(xFlux, prob, type = 6, na.rm = TRUE)
  lowFlux <- quantFlux[["5%"]]
  highFlux <- quantFlux[["95%"]]
  rejectC <- lowConc * highConc > 0
  rejectF <- lowFlux * highFlux > 0
  cat("\n\n  ", eList$INFO$shortName, "\n  ", eList$INFO$paramShortName)
  periodName <- EGRET::setSeasonLabelByUser(paStart, paLong)
  cat("\n  ", periodName, "\n")
  cat("\n Change estimates for\n average of", group2firstYear," through",group2lastYear,
      " minus average of", group1firstYear," through", group1lastYear, "\n")
  if (wall) 
    cat("\n Sample data set was partitioned with a wall at ", 
        as.character(sample1EndDate), "\n\n")
  cat("\n\nShould we reject Ho that Flow Normalized Concentration Trend = 0 ?", 
      words(rejectC))
  fquantConc <- format(quantConc, digits = 3, width = 8)
  cat("\n best estimate of change in concentration is", fcc, 
      "mg/L\n  Lower and Upper 90% CIs", fquantConc[["5%"]], fquantConc[["95%"]])
  lowC <- quantConc[["5%"]]
  upC <- quantConc[["95%"]]
  cat("\n also 95% CIs", fquantConc[["2.5%"]], fquantConc[["97.5%"]], "\n and 50% CIs", 
      fquantConc[["25%"]], fquantConc[["75%"]])
  lowC50 <- quantConc[["25%"]]
  upC50 <- quantConc[["75%"]]
  lowC95 <- quantConc[["2.5%"]]
  upC95 <- quantConc[["97.5%"]]
  pValC <- pVal(xConc)
  cat("\n approximate two-sided p-value for Conc", format(pValC, 
                                                          digits = 2, width = 9))
  xConc <- as.numeric(na.omit(xConc))
  nBootGood <- length(xConc)
  posX <- ifelse(xConc > 0, 1, 0)
  posXConc <- sum(posX)
  if (posXConc == 0 | posXConc == nBootGood) 
    cat("\n* Note p-value should be considered to be < stated value")
  likeCUp <- (posXConc + 0.5)/(nBootGood + 1)
  likeCDown <- 1 - likeCUp
  cat("\n Likelihood that Flow Normalized Concentration is trending up =", 
      format(likeCUp, digits = 3), " is trending down =", format(likeCDown, 
                                                                 digits = 3))
  if (nBootGood < nBoot) 
    cat("\n The number of good replicates in the bootstrap was ", 
        nBootGood, " out of the ", nBoot, "total")
  cat("\n\nShould we reject Ho that Flow Normalized Flux Trend = 0 ?", 
      words(rejectF))
  fquantFlux <- format(quantFlux, digits = 3, width = 8)
  cat("\n best estimate of change in flux is", ffc, "10^6 kg/year\n  Lower and Upper 90% CIs", 
      fquantFlux[["5%"]], fquantFlux[["95%"]])
  lowF <- quantFlux[["5%"]]
  upF <- quantFlux[["95%"]]
  cat("\n also 95% CIs", fquantFlux[["2.5%"]], fquantFlux[["97.5%"]], "\n and 50% CIs", 
      fquantFlux[["25%"]], fquantFlux[["75%"]])
  lowF50 <- quantFlux[["25%"]]
  upF50 <- quantFlux[["75%"]]
  lowF95 <- quantFlux[["2.5%"]]
  upF95 <- quantFlux[["97.5%"]]
  pValF <- pVal(xFlux)
  cat("\n approximate two-sided p-value for Flux", format(pValF, 
                                                          digits = 2, width = 9))
  xFlux <- as.numeric(na.omit(xFlux))
  nBootGood <- length(xFlux)
  posX <- ifelse(xFlux > 0, 1, 0)
  posXFlux <- sum(posX)
  if (posXFlux == 0 | posXFlux == nBootGood) 
    cat("\n* Note p-value should be considered to be < stated value")
  likeFUp <- (posXFlux + 0.5)/(nBootGood + 1)
  likeFDown <- 1 - likeFUp
  cat("\n Likelihood that Flow Normalized Flux is trending up =", 
      format(likeFUp, digits = 3), " is trending down =", format(likeFDown, 
                                                                 digits = 3))
  if (nBootGood < nBoot) 
    cat("\n The number of good replicates in the bootstrap was ", 
        nBootGood, " out of the ", nBoot, "total")
  bootOut <- data.frame(rejectC, pValC, estC, lowC, upC, lowC50, 
                        upC50, lowC95, upC95, likeCUp, likeCDown, rejectF, pValF, 
                        estF, lowF, upF, lowF50, upF50, lowF95, upF95, likeFUp, 
                        likeFDown, baseConc, baseFlux, nBoot, startSeed, blockLength, 
                        nBootGood)
  likeList <- c(likeCUp, likeCDown, likeFUp, likeFDown)
  wordsOut <- wordLike(likeList)
  cat("\n\n", format(wordsOut[1], width = 30), "\n", format(wordsOut[3], 
                                                            width = 30))
  cat("\n", format(wordsOut[2], width = 30), "\n", format(wordsOut[4], 
                                                          width = 30))
  pConc <- as.numeric(na.omit(pConc))
  pFlux <- as.numeric(na.omit(pFlux))
  groupBootOut <- list(bootOut = bootOut, wordsOut = wordsOut, 
                       xConc = xConc, xFlux = xFlux, pConc = pConc, pFlux = pFlux, 
                       startSeed = startSeed)
  attr(groupBootOut, "group1firstYear") <- group1firstYear
  attr(groupBootOut, "group1lastYear") <- group1lastYear
  attr(groupBootOut, "group2firstYear") <- group2firstYear
  attr(groupBootOut, "group2lastYear") <- group2lastYear
  return(groupBootOut)
}
