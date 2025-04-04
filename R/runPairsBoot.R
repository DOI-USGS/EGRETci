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
#' plotHistogramTrend(eList, boot_pair_out)
#' 
#' boot_message(eList, pairOut_2, boot_pair_out)
#' }
runPairsBoot <- function(eList, pairResults, 
                         nBoot = 100, startSeed = 494817, 
                         blockLength = 200,
                         jitterOn = FALSE, V = 0.2,
                         run.parallel = FALSE) {

  check_pair_dates(eList, pairResults)

  boot_return <- run_bootstraps(eList = eList, 
                                type_results = pairResults, 
                                jitterOn = jitterOn, V = V,
                                blockLength = blockLength, 
                                startSeed = startSeed,
                                nBoot = nBoot,
                                type = "pair",
                                run.parallel = run.parallel)

  pairsBootOut <- calc_boot_out(boot_list_return = boot_return,
                                type_results = pairResults, 
                                nBoot = nBoot, 
                                startSeed = startSeed,
                                blockLength = blockLength)
  
  attr(pairsBootOut, "year1") <- attr(pairResults, "yearPair")[["year1"]]
  attr(pairsBootOut, "year2") <- attr(pairResults, "yearPair")[["year2"]]
  
  attr(pairResults, "paStart") <- attr(pairResults, "yearPair")[["paStart"]]
  attr(pairResults, "paLong") <- attr(pairResults, "yearPair")[["paLong"]]

  boot_message(eList, 
               pairResults,
               pairsBootOut,
               type = "pair")
  
  return(pairsBootOut)
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


#' Calculations from boot return
#' 
#' @param boot_list_return The object that is returned from \code{run_bootstraps}.
#' @param type_results data frame returned from either \code{\link[EGRET]{runGroups}}
#' or \code{\link[EGRET]{runPairs}} depending on context.
#' @param nBoot the maximum number of bootstrap replicates to be used, typically 100
#' @param blockLength integer size of subset, expressed in days.  200 days has been found to be a good choice.
#' @param startSeed sets the random seed value. This is used to make repeatable output.
#' @keywords internal  
#' @export 
#' 
calc_boot_out <- function(boot_list_return,
                          type_results,
                          nBoot, 
                          startSeed,
                          blockLength){
  
  prob = c(0.025, 0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95, 0.975)
  
  xConc <- boot_list_return$xConc
  xFlux <- boot_list_return$xFlux
  pConc <- boot_list_return$pConc
  pFlux <- boot_list_return$pFlux
  
  quantConc <- quantile(xConc, prob, type = 6, na.rm = TRUE)
  quantFlux <- quantile(xFlux, prob, type = 6, na.rm = TRUE)

  lowConc <- quantConc[["5%"]]
  highConc <- quantConc[["95%"]]
  
  lowFlux <- quantFlux[["5%"]]
  highFlux <- quantFlux[["95%"]]
  
  rejectC <- lowConc * highConc > 0
  rejectF <- lowFlux * highFlux > 0

  nBootGood <- length(boot_list_return$xConc)
  
  posX <- ifelse(xConc > 0, 1, 0)
  posXConc <- sum(posX)
  
  likeCUp <- (posXConc + 0.5)/(nBootGood + 1)
  likeCDown <- 1 - likeCUp

  posX <- ifelse(xFlux > 0, 1, 0)
  posXFlux <- sum(posX)
  
  likeFUp <- (posXFlux + 0.5)/(nBootGood + 1)
  likeFDown <- 1 - likeFUp
  
  bootOut <- data.frame(rejectC = rejectC, 
                        pValC = pVal(boot_list_return$xConc), 
                        estC = type_results$x22[1] - type_results$x11[1], 
                        lowC = quantConc[["5%"]],
                        upC = quantConc[["95%"]], 
                        lowC50 = quantConc[["25%"]],
                        upC50 = quantConc[["75%"]],
                        lowC95 = quantConc[["2.5%"]],
                        upC95 = quantConc[["97.5%"]],
                        likeCUp = likeCUp,
                        likeCDown = likeCDown, 
                        rejectF = rejectF,
                        pValF =  pVal(xFlux),
                        estF = type_results$x22[2] - type_results$x11[2],
                        lowF = quantFlux[["5%"]],
                        upF = quantFlux[["95%"]],
                        lowF50 = quantFlux[["25%"]],
                        upF50 = quantFlux[["75%"]],
                        lowF95 = quantFlux[["2.5%"]], 
                        upF95 = quantFlux[["97.5%"]], 
                        likeFUp = likeFUp,
                        likeFDown = likeFDown,
                        baseConc = type_results$x11[1],
                        baseFlux = type_results$x11[2],
                        nBoot = nBoot, 
                        startSeed = startSeed, 
                        blockLength = blockLength,
                        nBootGood = nBootGood)
  
  likeList <- c(likeCUp, likeCDown, likeFUp, likeFDown)
  wordsOut <- wordLike(likeList)
  
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
#' @param type_results data frame returned from either \code{\link[EGRET]{runGroups}}
#' or \code{\link[EGRET]{runPairs}} depending on context.
#' @param bootOut List returned from either \code{runPairsBoot} or \code{runGroupBoot}.
#' @param type Character can be "pair" or "group".
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
#' boot_message(eList,
#'              pairOut_2,
#'              boot_pair_out)
#' }
boot_message <- function(eList, 
                         type_results,
                         bootOut,
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
  
  paStart <- eList$INFO$paStart
  paLong <- eList$INFO$paLong
  
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
  if (posXConc == 0 | posXConc == bootOut$bootOut$nBootGood) 
    cat("\n* Note p-value should be considered to be < stated value")
  
  cat("\n Likelihood that Flow Normalized Concentration is trending up =", 
      format(bootOut$bootOut$likeCUp, digits = 3), " is trending down =", 
      format(bootOut$bootOut$likeCDown, digits = 3))
  posX <- ifelse(bootOut$xFlux > 0, 1, 0)
  posXFlux <- sum(posX)

  # end of Concentration summary
  cat("\n\nShould we reject Ho that Flow Normalized Flux Trend = 0 ?", 
      words(bootOut$bootOut$rejectF))
  
  cat("\n best estimate of change in flux is", ffc, "10^6 kg/year\n  Lower and Upper 90% CIs", 
      bootOut$bootOut$lowF, bootOut$bootOut$upF)
  cat("\n also 95% CIs", bootOut$bootOut$lowF95, bootOut$bootOut$upF95, 
      "\n and 50% CIs", bootOut$bootOut$lowF50, bootOut$bootOut$upF50)
  
  cat("\n approximate two-sided p-value for Flux", format(bootOut$bootOut$pValF, 
                                                          digits = 2, width = 9))
  
  if (posXFlux == 0 | posXFlux == bootOut$bootOut$nBootGood) 
    cat("\n* Note p-value should be considered to be < stated value")
  
  cat("\n Likelihood that Flow Normalized Flux is trending up =", 
      format(bootOut$bootOut$likeFUp, digits = 3), " is trending down =", 
      format(bootOut$bootOut$likeFDown, digits = 3))
  
  if(bootOut$bootOut$nBootGood < bootOut$bootOut$nBoot) {
    cat("\n The number of good replicates in the bootstrap was ", 
        bootOut$bootOut$nBootGood,
        " out of the ", 
        bootOut$bootOut$nBoot, "total")
  }
  
  cat("\n\n", format(bootOut$wordsOut[1], width = 30), "\n", 
      format(bootOut$wordsOut[3], width = 30))
  
  cat("\n", format(bootOut$wordsOut[2], width = 30), "\n",
      format(bootOut$wordsOut[4], width = 30))
  
}


