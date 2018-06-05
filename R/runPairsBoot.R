
#' runPairsBoot
#' 
#' runPairsBoot
#' 
#' @param eList named list with at least the Daily, Sample, and INFO dataframes
#' @param pairResults data frame returned from \code{EGRET::runPairs}
#' @param nBoot the maximum number of bootstrap replicates to be used, typically 100
#' @param blockLength days, typically 200 is a good choice
#' @param startSeed setSeed value. Defaults to 494817. This is used to make repeatable output.
#' @export
#' @examples 
#' library(EGRET)
#' eList <- Choptank_eList
#' year1 <- 1985
#' year2 <- 2009
#' 
#' \dontrun{
#' pairOut_2 <- runPairs(eList, year1, year2, windowSide = 7)
#' 
#' boot_pair_out <- runPairsBoot(eList, pairOut_2)
#' 
#' plotHistogramTrend(eList, boot_pair_out, caseSetUp=NA)
#' }
runPairsBoot <- function(eList, pairResults, 
                         nBoot=100, startSeed = 494817, 
                         blockLength = 200) {

  interactive <- FALSE
  localINFO <- eList$INFO
  localDaily <- eList$Daily
  localSample <- eList$Sample
  
  firstDayDaily <- localDaily$Date[1]
  lastDayDaily <- localDaily$Date[length(localDaily$Date)]
  firstDaySample <- localSample$Date[1]
  lastDaySample <- localSample$Date[length(localSample$Date)]
  
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
                         "baseConc", "baseFlux", "nBoot","startSeed","blockLength")
  # get the parameters out of the attributes of pairResults
  paStart <- attr(pairResults, "yearPair")[["paStart"]]
  paLong <- attr(pairResults, "yearPair")[["paLong"]]
  year1 <- attr(pairResults, "yearPair")[["year1"]]
  year2 <- attr(pairResults, "yearPair")[["year2"]]

  startEnd1 <- EGRET::startEnd(paStart, paLong, year1)
  startEnd2 <- EGRET::startEnd(paStart, paLong, year2)

  
  if(startEnd2$startDate > range(eList$Sample$Date)[2]){
    stop("year2 is outside the Sample range")
  }
  
  if(startEnd1$endDate < range(eList$Sample$Date)[1]){
    stop("year1 is outside the Sample range")
  }
    
  start1 <- as.Date(startEnd1[["startDate"]])
  end1 <- as.Date(startEnd1[["endDate"]])
  start2 <- as.Date(startEnd2[["startDate"]])
  end2 <- as.Date(startEnd2[["endDate"]])

  dateInfo <- attr(pairResults, "dateInfo")
  
  sample1StartDate <- attr(pairResults, "SampleBlocks")[["sample1StartDate"]]
  sample1EndDate <- attr(pairResults, "SampleBlocks")[["sample1EndDate"]]
  sample2StartDate <- attr(pairResults, "SampleBlocks")[["sample2StartDate"]]
  sample2EndDate <- attr(pairResults, "SampleBlocks")[["sample2EndDate"]]
  
  minNumObs <- attr(pairResults, "Other")[["minNumObs"]]
  minNumUncen <- attr(pairResults, "Other")[["minNumUncen"]]
  windowY <- attr(pairResults, "Other")[["windowY"]]
  windowQ <- attr(pairResults, "Other")[["windowQ"]]
  windowY <- attr(pairResults, "Other")[["windowY"]]
  wall <- attr(pairResults, "Other")[["wall"]]
  edgeAdjust <- attr(pairResults, "Other")[["edgeAdjust"]]
  
  xConc <- rep(NA, nBoot)
  xFlux <- rep(NA, nBoot)
  pConc <- rep(NA, nBoot)
  pFlux <- rep(NA, nBoot)
  
  #
  regDeltaConc <- pairResults[1,7] - pairResults[1,5]
  estC <- regDeltaConc
  baseConc <- pairResults[1,5]
  regDeltaConcPct <- (regDeltaConc/baseConc) * 100
  LConcDiff <- log(pairResults[1,7]) - log(pairResults[1,5])
  #  pairResults are in 10^6 kg/year, when we get to the bootstrap results
  #  we will need to convert them all to 10^6 kg/year units
  regDeltaFlux <- (pairResults[2,7] - pairResults[2,5]) 
  estF <- regDeltaFlux
  baseFlux <- pairResults[2,5] 
  regDeltaFluxPct <- (regDeltaFlux/baseFlux) * 100
  LFluxDiff <- log(pairResults[2,7]) - log(pairResults[2,5])
  fcc <- format(regDeltaConc, digits = 3, width = 7)
  ffc <- format(regDeltaFlux, digits = 3, width = 8)
  Daily1 <- localDaily[localDaily$Date >= as.Date(dateInfo$flowNormStart[1]) & localDaily$Date <= 
                         as.Date(dateInfo$flowNormEnd[1]), ]
  Daily2 <- localDaily[localDaily$Date >= as.Date(dateInfo$flowNormStart[2]) & localDaily$Date <= 
                         as.Date(dateInfo$flowNormEnd[2]), ]
  # start inserting the setup stuff
  
  # bootstrap loop starts here
  nBootGood <- 0
  for (iBoot in 1:(2*nBoot)){

    bootSample <- blockSample(localSample = localSample, blockLength = blockLength, startSeed = startSeed + iBoot)
    eListBoot <- EGRET::as.egret(localINFO, localDaily, bootSample, NA)
    
    Sample1 <- bootSample[bootSample$Date >= sample1StartDate &
                            bootSample$Date <= sample1EndDate,]
    
    possibleError3 <- tryCatch( surfaces1 <- EGRET::estSurfaces(eListBoot, surfaceStart = start1, surfaceEnd = end1,
                                                         localSample = Sample1, minNumObs = minNumObs, minNumUncen = minNumUncen,
                                                         verbose = FALSE), error = function(e) e)
    
    Sample2 <- bootSample[bootSample$Date >= sample2StartDate &
                            bootSample$Date <= sample2EndDate,]
    
    possibleError4 <- tryCatch( surfaces2 <- EGRET::estSurfaces(eListBoot, surfaceStart = start2, surfaceEnd = end2,
                                                         localSample = Sample2, minNumObs = minNumObs, minNumUncen = minNumUncen,
                                                         verbose = FALSE), error = function(e) e)
    if (!inherits(possibleError3, "error") & 
        !inherits(possibleError4, "error")) {
  # note that all the flux calculations inside the bootstrap loop are in kg/day units    
      DailyRS1FD1 <- EGRET::estDailyFromSurfaces(eListBoot, localsurfaces = surfaces1, localDaily = Daily1)
      annualFlex <- EGRET::setupYears(DailyRS1FD1, paLong = paLong, paStart = paStart)
      c11 <- mean(annualFlex$FNConc, na.rm = TRUE)
      f11 <- mean(annualFlex$FNFlux, na.rm = TRUE) * 0.00036525
      
      DailyRS2FD2 <- EGRET::estDailyFromSurfaces(eListBoot, localsurfaces = surfaces2, localDaily = Daily2)
      annualFlex <- EGRET::setupYears(DailyRS2FD2, paLong = paLong, paStart = paStart)
      c22 <- mean(annualFlex$FNConc, na.rm = TRUE)
      f22 <- mean(annualFlex$FNFlux, na.rm = TRUE) * 0.00036525
      
      

      xConc_here <- (2 * regDeltaConc) - (c22 - c11)
      xFlux_here <- (2 * regDeltaFlux) - (f22 - f11)
      
      if(!is.na(xConc_here) & !is.na(xFlux_here)){
        nBootGood <- nBootGood + 1
      
        xConc[nBootGood] <- xConc_here
        xFlux[nBootGood] <- xFlux_here
        LConc <- (2 * LConcDiff) - (log(c22) - log(c11))
        pConc[nBootGood] <- (100 * exp(LConc)) - 100
        LFlux <- (2 * LFluxDiff) - (log(f22) - log(f11))
        pFlux[nBootGood] <- (100 * exp(LFlux)) - 100
        cat("\n iBoot, xConc and xFlux",nBootGood, xConc[nBootGood], xFlux[nBootGood])
    #  end of bootstrap replicates loop
        cat(nBootGood, "\n")
        if(nBootGood >= nBoot) {
          break()
        }
      }
    } else {
      stop(possibleError3, "\n", possibleError4)
    }
  }

  if(iBoot == 2*nBoot){
    message(iBoot, " iterations were run. They only achieved ", nBootGood, " sucessful runs.")
  } else if (iBoot > nBoot){
    message("It took ", iBoot, " iterations to achieve ", nBoot, " sucessful runs.")
  }
  # now summarize the bootstrap outputs
  quantConc <- quantile(xConc, prob, type = 6, na.rm = TRUE)
  lowConc <- quantConc[2]
  highConc <- quantConc[8]
  quantFlux <- quantile(xFlux, prob, type = 6, na.rm = TRUE)
  lowFlux <- quantFlux[2]
  highFlux <- quantFlux[8]
  rejectC <- lowConc * highConc > 0
  rejectF <- lowFlux * highFlux > 0
  cat("\n  ", eList$INFO$shortName, "\n  ", eList$INFO$paramShortName)
  periodName <- EGRET::setSeasonLabelByUser(paStart, paLong)
  cat("\n  ", periodName, "\n")
  if(wall) cat("\n Sample data set was partitioned with a wall at ", as.character(sample1EndDate), "\n\n")
  cat("\n\nShould we reject Ho that Flow Normalized Concentration Trend = 0 ?", 
      words(rejectC))
  fquantConc <- format(quantConc, digits = 3, width = 8)
  cat("\n best estimate of change in concentration is", fcc, "mg/L\n  Lower and Upper 90% CIs", 
      fquantConc[2], fquantConc[8])
  lowC <- quantConc[2]
  upC <- quantConc[8]
  cat("\n also 95% CIs", fquantConc[1], fquantConc[9], 
      "\n and 50% CIs", fquantConc[4], fquantConc[6])
  lowC50 <- quantConc[4]
  upC50 <- quantConc[6]
  lowC95 <- quantConc[1]
  upC95 <- quantConc[9]
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
      format(likeCUp, digits = 3), " is trending down =", 
      format(likeCDown, digits = 3))
  if(nBootGood < nBoot) cat("\n The number of good replicates in the bootstrap was ", nBootGood,
                            " out of the ", nBoot, "total")
  # end of Concentration summary
  cat("\n\nShould we reject Ho that Flow Normalized Flux Trend = 0 ?", 
      words(rejectF))
  fquantFlux <- format(quantFlux, digits = 3, width = 8)
  cat("\n best estimate of change in flux is", ffc, "10^6 kg/year\n  Lower and Upper 90% CIs", 
      fquantFlux[2], fquantFlux[8])
  lowF <- quantFlux[2]
  upF <- quantFlux[8]
  cat("\n also 95% CIs", fquantFlux[1], fquantFlux[9], 
      "\n and 50% CIs", fquantFlux[4], fquantFlux[6])
  lowF50 <- quantFlux[4]
  upF50 <- quantFlux[6]
  lowF95 <- quantFlux[1]
  upF95 <- quantFlux[9]
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
  cat("\n Likelihood that Flow Normalized Concentration is trending up =", 
      format(likeFUp, digits = 3), " is trending down =", 
      format(likeFDown, digits = 3))
  if(nBootGood < nBoot) cat("\n The number of good replicates in the bootstrap was ", nBootGood,
                            " out of the ", nBoot, "total")
  bootOut <- data.frame(rejectC, pValC, estC, lowC, upC, 
                        lowC50, upC50, lowC95, upC95, likeCUp, likeCDown, 
                        rejectF, pValF, estF, lowF, upF, lowF50, upF50, lowF95, 
                        upF95, likeFUp, likeFDown, baseConc, baseFlux, nBoot, 
                        startSeed, blockLength, nBootGood)
  likeList <- c(likeCUp, likeCDown, likeFUp, likeFDown)
  wordsOut <- wordLike(likeList)
  cat("\n\n", format(wordsOut[1], width = 30), "\n", format(wordsOut[3], 
                                                               width = 30))
  cat("\n", format(wordsOut[2], width = 30), "\n", format(wordsOut[4], 
                                                            width = 30))
  pConc <- as.numeric(na.omit(pConc))
  pFlux <- as.numeric(na.omit(pFlux))
  pairsBootOut <- list(bootOut = bootOut, wordsOut = wordsOut, 
                xConc = xConc, xFlux = xFlux, pConc = pConc, pFlux = pFlux,
                startSeed = startSeed)
  attr(pairsBootOut, "year1") <- year1
  attr(pairsBootOut, "year2") <- year2
  return(pairsBootOut)
}
