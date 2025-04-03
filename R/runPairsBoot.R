#' The bootstrap uncertainty analysis for runPairs results
#' 
#' The function that does the uncertainty analysis for determining the change between any 
#' pair of years.  It is very similar to the \code{\link{wBT}} function that runs the WRTDS 
#' bootstrap test.  It differs from \code{\link{wBT}} in that it runs a specific number of 
#' bootstrap replicates, unlike the \code{\link{wBT}} approach that will stop running replicates 
#' based on the status of the test statistics along the way.  Also, this code can be used with
#' generalized flow normalization, which handles non-stationary discharge, 
#' whereas \code{\link{wBT}} does not. 
#' 
#' @param eList named list with at least the Daily, Sample, and INFO dataframes
#' @param pairResults data frame returned from \code{\link[EGRET]{runPairs}}
#' @param nBoot the maximum number of bootstrap replicates to be used, typically 100
#' @param blockLength integer size of subset, expressed in days.  200 days has been found to be a good choice.
#' @param startSeed setSeed value. Defaults to 494817. This is used to make repeatable output.
#' @param jitterOn logical, if TRUE, adds "jitter" to the data in an attempt to avoid some numerical problems.  Default = FALSE.  See Details below.
#' @param V numeric a multiplier for addition of jitter to the data, default = 0.2.
#' @param run.parallel logical to run bootstrapping in parallel or not
#' @export
#' @details
#' In some situations numerical problems are encountered in the bootstrap process, resulting in highly unreasonable spikes in the confidence intervals.
#' The use of "jitter" can often prevent these problems, but should only be used when it is clearly needed.
#' It adds a small amount of random "jitter" to the explanatory variables of the WRTDS model.  The V parameter sets the scale of variation in the log discharge values.
#' The standard deviation of the added jitter is V * standard deviation of Log Q.
#' The default for V is 0.2.  Larger values should generally be avoided, and smaller values may be sufficient.
#' @return eBoot, a named list with bootOut, wordsOut, xConc, xFlux, pConc, pFlux values.
#' \itemize{
#'   \item{bootOut is a data frame with the results of the bootstrap test. }
#'   \item{wordsOut is a character vector describing the results.}
#'   \item{xConc and xFlux are vectors of length iBoot, of the change in flow normalized concentration
#'    and flow normalized flux computed from each of the bootstrap replicates. }
#'   \item{pConc and pFlux are vectors of length iBoot, of the change in flow normalized concentration
#'    or flow normalized flux computed from each of the bootstrap replicates expressed as \% change.}
#' }
#' @seealso \code{\link{runGroupsBoot}}, \code{\link[EGRET]{runPairs}}
#' @examples 
#' eList <- EGRET::Choptank_eList
#' year1 <- 1985
#' year2 <- 2009
#' 
#' \dontrun{
#' pairOut_2 <- EGRET::runPairs(eList, 
#'                              year1, year2, 
#'                              windowSide = 7)
#' 
#' # For good analysis, bump up nBoot to about 100:
#' boot_pair_out <- runPairsBoot(eList, pairOut_2, nBoot = 5)
#' 
#' plotHistogramTrend(eList, boot_pair_out, caseSetUp = NA)
#' }
runPairsBoot <- function(eList, pairResults, 
                         nBoot = 100, startSeed = 494817, 
                         blockLength = 200,
                         jitterOn = FALSE, V = 0.2,
                         run.parallel = FALSE) {

  check_pair_dates(eList, pairResults)

  # bootstrap loop starts here
  nBootGood <- 0
  xConc <- c()
  xFlux <- c()
  pConc <- c()
  pFlux <- c()
  
  if(run.parallel){
    `%dopar%` <- foreach::`%dopar%`
    boot_list_out <- foreach::foreach(iBoot = 1:ceiling(1.25*nBoot), 
                                          .packages=c('EGRETci', 'EGRET')) %dopar% {
                                            boot_list <- boot_pairs_run(iBoot, startSeed,
                                                                        eList, pairResults, 
                                                                        jitterOn, V,
                                                                        blockLength)
                                          }
    
    xConc <- unlist(sapply(boot_list_out, function(x) x[["xConc"]]))
    xFlux <- unlist(sapply(boot_list_out, function(x) x[["xFlux"]]))
    pConc <- unlist(sapply(boot_list_out, function(x) x[["pConc"]]))
    pFlux <- unlist(sapply(boot_list_out, function(x) x[["pFlux"]]))
    
    if(length(xConc) < nBoot){
      # do the last set NOT in parallel because there's a lot of overhead,
      # potentially for just a few runs:
      iStart <- nBoot + 1
      iEnd <- 2*nBoot - length(xConc)
      message("Running ", iStart-iEnd, " in series")
      for(i in iStart:iEnd){
        boot_list <- boot_pairs_run(iBoot, startSeed,
                                    eList, pairResults, 
                                    jitterOn, V,
                                    blockLength)
        
        if(!is.null(boot_list$xConc) & !is.null(boot_list$xFlux)){
          if(nBootGood >= nBoot) {
            break()
          }
          nBootGood <- nBootGood + 1
        }
        xConc <- c(xConc, boot_list$xConc)
        xFlux <- c(xFlux, boot_list$xFlux)
        pConc <- c(pConc, boot_list$pConc)
        pFlux <- c(pFlux, boot_list$pFlux)
      }
      iBoot <- length(pFlux)
    } else {
      xConc <- xConc[1:nBoot]
      xFlux <- xFlux[1:nBoot]
      pConc <- pConc[1:nBoot]
      pFlux <- pFlux[1:nBoot]
      iBoot <- nBoot
    }
    
  } else {
    for (iBoot in 1:(2*nBoot)){
      boot_list <- boot_pairs_run(iBoot, startSeed,
                                  eList, pairResults, 
                                  jitterOn, V,
                                  blockLength)
      
      if(!is.null(boot_list$xConc) & !is.null(boot_list$xFlux)){
        if(nBootGood >= nBoot) {
          break()
        }
        nBootGood <- nBootGood + 1
      }
      xConc <- c(xConc, boot_list$xConc)
      xFlux <- c(xFlux, boot_list$xFlux)
      pConc <- c(pConc, boot_list$pConc)
      pFlux <- c(pFlux, boot_list$pFlux)
    }  
  }
  
  if(iBoot == 2*nBoot){
    message(iBoot, " iterations were run. They only achieved ", nBootGood, " sucessful runs.")
  } else if (iBoot > nBoot){
    message("It took ", iBoot, " iterations to achieve ", nBoot, " sucessful runs.")
  }
  
  pairsBootOut <- calc_boot_out(xConc, xFlux, pConc, pFlux,
                                pairResults, nBoot, startSeed,
                                blockLength, nBootGood)
  
  attr(pairsBootOut, "year1") <- attr(pairResults, "yearPair")[["year1"]]
  attr(pairsBootOut, "year2") <- attr(pairResults, "yearPair")[["year2"]]
  
  attr(pairResults, "paStart") <- attr(pairResults, "yearPair")[["paStart"]]
  attr(pairResults, "paLong") <- attr(pairResults, "yearPair")[["paLong"]]

  pair_boot_message(eList, 
                    pairResults,
                    pairsBootOut,
                    nBootGood, nBoot,
                    type = "pair")
  
  return(pairsBootOut)
}

boot_pairs_run <- function(iBoot, startSeed,
                           eList, pairResults, 
                           jitterOn, V,
                           blockLength){
  
  localINFO <- eList$INFO
  localDaily <- eList$Daily
  localSample <- eList$Sample
  
  dateInfo <- attr(pairResults, "dateInfo")

  Daily1 <- localDaily[localDaily$Date >= as.Date(dateInfo$flowNormStart[1]) & localDaily$Date <= 
                         as.Date(dateInfo$flowNormEnd[1]), ]
  Daily2 <- localDaily[localDaily$Date >= as.Date(dateInfo$flowNormStart[2]) & localDaily$Date <= 
                         as.Date(dateInfo$flowNormEnd[2]), ]
  
  bootSample <- blockSample(localSample = localSample, 
                            blockLength = blockLength, 
                            startSeed = startSeed + iBoot)
  
  if(jitterOn) bootSample <- EGRET::jitterSam(bootSample, V = V)
  
  eListBoot <- suppressMessages(EGRET::as.egret(localINFO, 
                                                localDaily,
                                                bootSample, NA))
  
  startEnd1 <- EGRET::startEnd(attr(pairResults, "yearPair")[["paStart"]],
                               attr(pairResults, "yearPair")[["paLong"]], 
                               attr(pairResults, "yearPair")[["year1"]])
  startEnd2 <- EGRET::startEnd(attr(pairResults, "yearPair")[["paStart"]], 
                               attr(pairResults, "yearPair")[["paLong"]],
                               attr(pairResults, "yearPair")[["year2"]])
  
  Sample1 <- bootSample[bootSample$Date >= attr(pairResults, "SampleBlocks")[["sample1StartDate"]] &
                          bootSample$Date <= attr(pairResults, "SampleBlocks")[["sample1EndDate"]],]
  
  possibleError3 <- tryCatch( 
    surfaces1 <- suppressMessages(EGRET::estSurfaces(eListBoot, 
                                                     surfaceStart = as.Date(startEnd1[["startDate"]]), 
                                                     surfaceEnd = as.Date(startEnd1[["endDate"]]),
                                                     edgeAdjust = attr(pairResults, "Other")[["edgeAdjust"]],
                                                     localSample = Sample1, 
                                                     minNumObs = attr(pairResults, "Other")[["minNumObs"]], 
                                                     minNumUncen = attr(pairResults, "Other")[["minNumUncen"]],
                                                     verbose = FALSE)),
    error = function(e) e)
  
  Sample2 <- bootSample[bootSample$Date >= attr(pairResults, "SampleBlocks")[["sample2StartDate"]] &
                          bootSample$Date <= attr(pairResults, "SampleBlocks")[["sample2EndDate"]], ]
  
  possibleError4 <- tryCatch(
    surfaces2 <- suppressMessages(EGRET::estSurfaces(eListBoot, 
                                                     surfaceStart = as.Date(startEnd2[["startDate"]]), 
                                                     surfaceEnd = as.Date(startEnd2[["endDate"]]),
                                                     edgeAdjust = attr(pairResults, "Other")[["edgeAdjust"]],
                                                     localSample = Sample2, 
                                                     minNumObs = attr(pairResults, "Other")[["minNumObs"]],
                                                     minNumUncen = attr(pairResults, "Other")[["minNumUncen"]],
                                                     verbose = FALSE)),
    error = function(e) e)
  
  if (!inherits(possibleError3, "error") & 
      !inherits(possibleError4, "error")) {
    # note that all the flux calculations inside the bootstrap loop are in kg/day units    
    DailyRS1FD1 <- EGRET::estDailyFromSurfaces(eListBoot, 
                                               localsurfaces = surfaces1, 
                                               localDaily = Daily1)
    annualFlex <- EGRET::setupYears(DailyRS1FD1, 
                                    paLong = attr(pairResults, "yearPair")[["paLong"]], 
                                    paStart = attr(pairResults, "yearPair")[["paStart"]])
    
    #  runPairs are in 10^6 kg/year, when we get to the bootstrap results
    #  Converting them all to 10^6 kg/year units
    c11 <- mean(annualFlex$FNConc, na.rm = TRUE)
    f11 <- mean(annualFlex$FNFlux, na.rm = TRUE) * 0.00036525
    
    DailyRS2FD2 <- EGRET::estDailyFromSurfaces(eListBoot, 
                                               localsurfaces = surfaces2, 
                                               localDaily = Daily2)
    annualFlex <- EGRET::setupYears(DailyRS2FD2, 
                                    paLong = attr(pairResults, "yearPair")[["paLong"]], 
                                    paStart = attr(pairResults, "yearPair")[["paStart"]])
    
    c22 <- mean(annualFlex$FNConc, na.rm = TRUE)
    f22 <- mean(annualFlex$FNFlux, na.rm = TRUE) * 0.00036525
    
    regDeltaConc <- pairResults$x22[1] - pairResults$x11[1]
    regDeltaFlux <- pairResults$x22[2] - pairResults$x11[2]
    LFluxDiff <- log(pairResults$x22[2]) - log(pairResults$x11[2])
    LConcDiff <- log(pairResults$x22[1]) - log(pairResults$x11[1])
    
    xConc <- (2 * regDeltaConc) - (c22 - c11)
    xFlux <- (2 * regDeltaFlux) - (f22 - f11)
    
    if(!is.na(xConc) & !is.na(xFlux)){

      LConc <- (2 * LConcDiff) - (log(c22) - log(c11))
      pConc <- (100 * exp(LConc)) - 100
      LFlux <- (2 * LFluxDiff) - (log(f22) - log(f11))
      pFlux <- (100 * exp(LFlux)) - 100
      message("\n iBoot, xConc and xFlux ",iBoot, ": ", 
              round(xConc, digits = 4), " ", 
              round(xFlux, digits = 4))
    }
    
    return_list <- list(xConc = xConc,
                        xFlux = xFlux,
                        pConc = pConc,
                        pFlux = pFlux)
  } else {
    return_list <- list(xConc = NULL,
                        xFlux = NULL,
                        pConc = NULL,
                        pFlux = NULL)    
  }
  
  return(return_list)
}

check_pair_dates <- function(eList, pairResults){
  localDaily <- eList$Daily
  localSample <- eList$Sample
  
  firstDayDaily <- min(localDaily$Date, na.rm = TRUE)
  lastDayDaily <- max(localDaily$Date, na.rm = TRUE)
  firstDaySample <- min(localSample$Date, na.rm = TRUE)
  lastDaySample <- max(localSample$Date, na.rm = TRUE)
  
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
}

calc_boot_out <- function(xConc,
                          xFlux,
                          pConc,
                          pFlux,
                          type_results, nBoot, startSeed,
                          blockLength, nBootGood){
  
  prob = c(0.025, 0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95, 0.975)
  quantConc <- quantile(xConc, prob, type = 6, na.rm = TRUE)
  lowConc <- quantConc[["5%"]]
  highConc <- quantConc[["95%"]]
  quantFlux <- quantile(xFlux, prob, type = 6, na.rm = TRUE)
  lowFlux <- quantFlux[["5%"]]
  highFlux <- quantFlux[["95%"]]
  rejectC <- lowConc * highConc > 0
  rejectF <- lowFlux * highFlux > 0
  
  fquantConc <- format(quantConc, digits = 3, width = 8)
  
  lowC <- quantConc[["5%"]]
  upC <- quantConc[["95%"]]
  lowC50 <- quantConc[["25%"]]
  upC50 <- quantConc[["75%"]]
  lowC95 <- quantConc[["2.5%"]]
  upC95 <- quantConc[["97.5%"]]
  pValC <- pVal(xConc)
  
  xConc <- as.numeric(na.omit(xConc))
  nBootGood <- length(xConc)
  
  posX <- ifelse(xConc > 0, 1, 0)
  posXConc <- sum(posX)
  
  likeCUp <- (posXConc + 0.5)/(nBootGood + 1)
  likeCDown <- 1 - likeCUp
  
  lowF <- quantFlux[["5%"]]
  upF <- quantFlux[["95%"]]
  lowF50 <- quantFlux[["25%"]]
  upF50 <- quantFlux[["75%"]]
  lowF95 <- quantFlux[["2.5%"]]
  upF95 <- quantFlux[["97.5%"]]
  
  pValF <- pVal(xFlux)
  
  xFlux <- as.numeric(na.omit(xFlux))
  nBootGood <- length(xFlux)
  
  posX <- ifelse(xFlux > 0, 1, 0)
  posXFlux <- sum(posX)
  
  likeFUp <- (posXFlux + 0.5)/(nBootGood + 1)
  likeFDown <- 1 - likeFUp
  
  estC <- type_results$x22[1] - type_results$x11[1]
  estF <- type_results$x22[2] - type_results$x11[2]
  
  baseConc <- type_results$x11[1]
  baseFlux <- type_results$x11[2]
  
  bootOut <- data.frame(rejectC, pValC, estC, lowC, upC, 
                        lowC50, upC50, lowC95, upC95, likeCUp, likeCDown, 
                        rejectF, pValF, estF, lowF, upF, lowF50, upF50, lowF95, 
                        upF95, likeFUp, likeFDown, baseConc, baseFlux, nBoot, 
                        startSeed, blockLength, nBootGood)
  
  likeList <- c(likeCUp, likeCDown, likeFUp, likeFDown)
  wordsOut <- wordLike(likeList)
  
  pConc <- as.numeric(na.omit(pConc))
  pFlux <- as.numeric(na.omit(pFlux))
  
  boot_out <- list(bootOut = bootOut,
                   wordsOut = wordsOut, 
                   xConc = xConc, 
                   xFlux = xFlux, 
                   pConc = pConc, 
                   pFlux = pFlux)
  boot_out
}

#' Print message of bootstrapping results
#' 
#' Uses the results of either group or pair results and prints
#' of the standard message.
#' 
#' @export
#' @param eList named list with at least the Daily, Sample, and INFO dataframes
#' @param type_results List returned from either runPairs or runGroup.
#' @param bootOut List returned from either runPairs or runGroupBoot.
#' @param nBootGood Number of good results.
#' @param nBoot Number of requested bootstraps.
#' @param type Character can be "pair" or "group"
#' 
pair_boot_message <- function(eList, 
                              type_results,
                              bootOut,
                              nBootGood, nBoot,
                              type = "pair"){
  
  match.arg(type, choices = c("pair", "group"))
  
  wordsOut <- bootOut[["wordsOut"]]
  
  regDeltaConc <- type_results$x22[1] - type_results$x11[1]
  estC <- regDeltaConc
  baseConc <- type_results$x11[1]
  regDeltaConcPct <- (regDeltaConc/baseConc) * 100
  LConcDiff <- log(type_results$x22[1]) - log(type_results$x11[1])
  #  type_results are in 10^6 kg/year, when we get to the bootstrap results
  #  we will need to convert them all to 10^6 kg/year units
  regDeltaFlux <- type_results$x22[2] - type_results$x11[2]
  estF <- regDeltaFlux
  baseFlux <- type_results$x11[2] 
  regDeltaFluxPct <- (regDeltaFlux/baseFlux) * 100
  LFluxDiff <- log(type_results$x22[2]) - log(type_results$x11[2])
  fcc <- format(regDeltaConc, digits = 3, width = 7)
  ffc <- format(regDeltaFlux, digits = 3, width = 8)
  
  paStart <- attr(type_results, "paStart")
  paLong <- attr(type_results, "paLong")
  
  cat("\n  ", eList$INFO$shortName, "\n  ", eList$INFO$paramShortName)
  periodName <- EGRET::setSeasonLabelByUser(paStart, paLong)
  cat("\n  ", periodName, "\n")
  if(type == "pair"){
    cat("\n Change estimates are for ", attr(type_results, "yearPair")[["year2"]],
        " minus ", attr(type_results, "yearPair")[["year1"]])    
  } else if(type == "group") {
    cat("\n Change estimates for\n average of", attr(type_results, "groupInfo")[["group2firstYear"]],
        " through", attr(type_results, "groupInfo")[["group2lastYear"]],
        " minus average of", attr(type_results, "groupInfo")[["group1firstYear"]],
        " through", attr(type_results, "groupInfo")[["group1lastYear"]], "\n")
  }

  if(attr(type_results, "Other")[["wall"]]) cat("\n Sample data set was partitioned with a wall at ",
                                               as.character(attr(type_results, "SampleBlocks")[["sample1EndDate"]]), "\n\n")
  cat("\n\nShould we reject Ho that Flow Normalized Concentration Trend = 0 ?", 
      words(bootOut$bootOut$rejectC))
  
  cat("\n best estimate of change in concentration is", fcc, "mg/L\n  Lower and Upper 90% CIs", 
      bootOut$bootOut$lowC, bootOut$bootOut$upC)
  
  cat("\n also 95% CIs", bootOut$bootOut$lowC95, bootOut$bootOut$upC95, 
      "\n and 50% CIs", bootOut$bootOut$lowC50, bootOut$bootOut$upC50)
  
  cat("\n approximate two-sided p-value for Conc", format(bootOut$bootOut$pValC, 
                                                          digits = 2, width = 9))
  posX <- ifelse(bootOut$xConc > 0, 1, 0)
  posXConc <- sum(posX)
  if (posXConc == 0 | posXConc == nBootGood) 
    cat("\n* Note p-value should be considered to be < stated value")
  
  cat("\n Likelihood that Flow Normalized Concentration is trending up =", 
      format(bootOut$bootOut$likeCUp, digits = 3), " is trending down =", 
      format(bootOut$bootOut$likeCDown, digits = 3))
  posX <- ifelse(bootOut$xFlux > 0, 1, 0)
  posXFlux <- sum(posX)
  if (posXFlux == 0 | posXFlux == nBootGood) 
    cat("\n* Note p-value should be considered to be < stated value")
  
  if(nBootGood < nBoot) cat("\n The number of good replicates in the bootstrap was ", nBootGood,
                            " out of the ", nBoot, "total")
  # end of Concentration summary
  cat("\n\nShould we reject Ho that Flow Normalized Flux Trend = 0 ?", 
      words(bootOut$bootOut$rejectF))
  
  cat("\n best estimate of change in flux is", ffc, "10^6 kg/year\n  Lower and Upper 90% CIs", 
      bootOut$bootOut$lowF, bootOut$bootOut$upF)
  cat("\n also 95% CIs", bootOut$bootOut$lowF95, bootOut$bootOut$upF95, 
      "\n and 50% CIs", bootOut$bootOut$lowF50, bootOut$bootOut$upF50)
  
  cat("\n approximate two-sided p-value for Flux", format(bootOut$bootOut$pValF, 
                                                          digits = 2, width = 9))
  
  if (posXFlux == 0 | posXFlux == nBootGood) 
    cat("\n* Note p-value should be considered to be < stated value")
  
  cat("\n Likelihood that Flow Normalized Flux is trending up =", 
      format(bootOut$bootOut$likeFUp, digits = 3), " is trending down =", 
      format(bootOut$bootOut$likeFDown, digits = 3))
  
  if(nBootGood < nBoot) cat("\n The number of good replicates in the bootstrap was ", nBootGood,
                            " out of the ", nBoot, "total")
  
  cat("\n\n", format(bootOut$wordsOut[1], width = 30), "\n", 
      format(bootOut$wordsOut[3], width = 30))
  
  cat("\n", format(bootOut$wordsOut[2], width = 30), "\n",
      format(bootOut$wordsOut[4], width = 30))
  
}


