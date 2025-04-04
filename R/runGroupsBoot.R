
#' The bootstrap uncertainty analysis for runGroups results
#' 
#' This function that does the uncertainty analysis for determining the change 
#' between two groups of years.  The process is virtually 
#' identical to what is used for \code{\link{runPairsBoot}} which looks at a change
#' between a pair of years.  
#' 
#' @param eList named list with at least the Daily, Sample, and INFO dataframes
#' @param groupResults data frame returned from \code{\link[EGRET]{runGroups}}
#' @param nBoot the maximum number of bootstrap replicates to be used, typically 100
#' @param blockLength integer size of subset, expressed in days.  200 days has been found to be a good choice.
#' @param startSeed setSeed value. Defaults to 494817. This is used to make repeatable output.
#' @param jitterOn logical, if TRUE, adds "jitter" to the data in an attempt to avoid some numerical problems.
#'   Default = FALSE.  See Details below.
#' @param V numeric a multiplier for addition of jitter to the data, default = 0.2.
#' @param run.parallel logical to run bootstrapping in parallel or not
#' @return eBoot, a named list with bootOut, wordsOut, xConc, xFlux, pConc, pFlux values.
#' \itemize{
#'   \item{bootOut is a data frame with the results of the bootstrap test.}
#'   \item{wordsOut is a character vector describing the results.}
#'   \item{xConc and xFlux are vectors of length iBoot, of the change in flow normalized concentration
#'    and flow normalized flux computed from each of the bootstrap replicates.}
#'   \item{pConc and pFlux are vectors of length iBoot, of the change in flow normalized concentration
#'    or flow normalized flux computed from each of the bootstrap replicates expressed as \% change.}
#'}
#' @seealso \code{\link{runPairsBoot}}, \code{\link[EGRET]{runGroups}}
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
#' # For good analysis, bump up nBoot to about 100:
#' boot_group_out <- runGroupsBoot(eList, groupResults, nBoot = 3)
#' 
#' }
runGroupsBoot <- function (eList, groupResults, nBoot = 100, 
                           startSeed = 494817, blockLength = 200,
                           jitterOn = FALSE, V = 0.2,
                           run.parallel = FALSE){

  boot_return <- run_bootstraps(eList = eList, 
                                type_results = groupResults, 
                                jitterOn = jitterOn, V = V,
                                blockLength = blockLength, 
                                startSeed = startSeed,
                                nBoot = nBoot,
                                type = "group",
                                run.parallel = run.parallel)
  
  groupBootOut <- calc_boot_out(boot_list_return = boot_return,
                                type_results = groupResults, 
                                nBoot = nBoot, 
                                startSeed = startSeed,
                                blockLength = blockLength, 
                                nBootGood = length(boot_return$xConc))
  
  attr(groupResults, "paStart") <- attr(groupResults, "groupInfo")[["paStart"]]
  attr(groupResults, "paLong") <- attr(groupResults, "groupInfo")[["paLong"]]
  
  boot_message(eList = eList,
               type_results = groupResults,
               bootOut = groupBootOut,
               nBootGood = length(boot_return$xConc),
               nBoot = nBoot,
               type = "group")
  
  return(groupBootOut)
}

single_boot_run <- function(iBoot, startSeed,
                           eList, type_results, 
                           jitterOn, V,
                           blockLength, type){
  
  match.arg(type, choices = c("pair", "group"))
  
  bootSample <- blockSample(localSample = eList$Sample, 
                            blockLength = blockLength, 
                            startSeed = startSeed + iBoot)
  
  if(jitterOn) bootSample <- EGRET::jitterSam(bootSample, V = V)
  
  eListBoot <- suppressMessages(EGRET::as.egret(eList$INFO,
                                                eList$Daily,
                                                bootSample,NA))
  dateInfo <- attr(type_results, "dateInfo")
  
  if(type == "group"){
    
    if(attr(type_results, "Other")[["wall"]]) {
      possibleError <- tryCatch(surfaces <- suppressMessages(EGRET::stitch(eListBoot, 
                                                                           surfaceStart = attr(type_results, "SampleBlocks")[["surfaceStart"]], 
                                                                           surfaceEnd =  attr(type_results, "SampleBlocks")[["surfaceEnd"]],
                                                                           sample1StartDate = attr(type_results, "SampleBlocks")[["sample1StartDate"]],
                                                                           sample1EndDate = attr(type_results, "SampleBlocks")[["sample1EndDate"]],
                                                                           sample2StartDate = attr(type_results, "SampleBlocks")[["sample2StartDate"]], 
                                                                           sample2EndDate = attr(type_results, "SampleBlocks")[["sample2EndDate"]],
                                                                           windowY = attr(type_results, "Other")[["windowY"]], 
                                                                           windowQ = attr(type_results, "Other")[["windowQ"]], 
                                                                           windowS = attr(type_results, "Other")[["windowS"]],
                                                                           minNumObs = attr(type_results, "Other")[["minNumObs"]], 
                                                                           minNumUncen =  attr(type_results, "Other")[["minNumUncen"]], 
                                                                           edgeAdjust =  attr(type_results, "Other")[["edgeAdjust"]])), 
                                error = function(e) e)
    } else {
      possibleError <- tryCatch(surfaces <- EGRET::estSurfaces(eListBoot, 
                                                               surfaceStart = attr(type_results, "SampleBlocks")[["surfaceStart"]], 
                                                               surfaceEnd =  attr(type_results, "SampleBlocks")[["surfaceEnd"]],
                                                               windowY = attr(type_results, "Other")[["windowY"]],
                                                               windowQ = attr(type_results, "Other")[["windowQ"]],
                                                               windowS = attr(type_results, "Other")[["windowS"]],
                                                               minNumObs = attr(type_results, "Other")[["minNumObs"]],
                                                               minNumUncen =  attr(type_results, "Other")[["minNumUncen"]],
                                                               edgeAdjust =  attr(type_results, "Other")[["edgeAdjust"]]),
                                error = function(e) e)
    }
    
    if (!inherits(possibleError, "error") ) {
      eListS <- suppressMessages(EGRET::as.egret(eListBoot$INFO, 
                                                 eListBoot$Daily,
                                                 eListBoot$Sample, surfaces))
      eListOut <- suppressMessages(EGRET::flexFN(eListS, dateInfo,
                                                 flowNormStartCol = "flowNormStart", 
                                                 flowNormEndCol = "flowNormEnd", flowStartCol = "flowStart", 
                                                 flowEndCol = "flowEnd"))
      
      eListOut$INFO$wall <- attr(type_results, "Other")[["wall"]]
      eListOut$INFO$surfaceStart <- attr(type_results, "SampleBlocks")[["surfaceStart"]]
      eListOut$INFO$surfaceEnd <- attr(type_results, "SampleBlocks")[["surfaceEnd"]]
      DailyFlex1 <- eListOut$Daily
      annFlex <- EGRET::setupYears(DailyFlex1, 
                                   paLong =  attr(type_results, "groupInfo")[["paLong"]],
                                   paStart =  attr(type_results, "groupInfo")[["paStart"]])
      annFlex$year <- floor(annFlex$DecYear + (annFlex$PeriodLong / 12) * 0.5)
      annFlex1 <- annFlex[annFlex$DecYear >= attr(type_results, "groupInfo")[["group1firstYear"]] &
                            annFlex$DecYear <= attr(type_results, "groupInfo")[["group1lastYear"]],]
      annFlex2 <- annFlex[annFlex$DecYear >= attr(type_results, "groupInfo")[["group2firstYear"]] & 
                            annFlex$DecYear <= attr(type_results, "groupInfo")[["group2lastYear"]],]
      
    } else {
      return(list(xConc = NULL,
                  xFlux = NULL,
                  pConc = NULL,
                  pFlux = NULL))
    }    
  } else if (type == "pair") {
    localDaily <- eList$Daily
    Daily1 <- localDaily[localDaily$Date >= as.Date(dateInfo$flowNormStart[1]) & 
                           localDaily$Date <= as.Date(dateInfo$flowNormEnd[1]), ]
    Daily2 <- localDaily[localDaily$Date >= as.Date(dateInfo$flowNormStart[2]) & 
                           localDaily$Date <= as.Date(dateInfo$flowNormEnd[2]), ]
    
    startEnd1 <- EGRET::startEnd(attr(type_results, "yearPair")[["paStart"]],
                                 attr(type_results, "yearPair")[["paLong"]], 
                                 attr(type_results, "yearPair")[["year1"]])
    startEnd2 <- EGRET::startEnd(attr(type_results, "yearPair")[["paStart"]], 
                                 attr(type_results, "yearPair")[["paLong"]],
                                 attr(type_results, "yearPair")[["year2"]])
    
    Sample1 <- bootSample[bootSample$Date >= attr(type_results, "SampleBlocks")[["sample1StartDate"]] &
                            bootSample$Date <= attr(type_results, "SampleBlocks")[["sample1EndDate"]],]
    Sample2 <- bootSample[bootSample$Date >= attr(type_results, "SampleBlocks")[["sample2StartDate"]] &
                            bootSample$Date <= attr(type_results, "SampleBlocks")[["sample2EndDate"]], ]
    
    possibleError3 <- tryCatch(
      surfaces1 <- suppressMessages(EGRET::estSurfaces(eListBoot, 
                                                       surfaceStart = as.Date(startEnd1[["startDate"]]), 
                                                       surfaceEnd = as.Date(startEnd1[["endDate"]]),
                                                       edgeAdjust = attr(type_results, "Other")[["edgeAdjust"]],
                                                       localSample = Sample1, 
                                                       minNumObs = attr(type_results, "Other")[["minNumObs"]], 
                                                       minNumUncen = attr(type_results, "Other")[["minNumUncen"]],
                                                       verbose = FALSE)),
      error = function(e) e)
    
    possibleError4 <- tryCatch(
      surfaces2 <- suppressMessages(EGRET::estSurfaces(eListBoot, 
                                                       surfaceStart = as.Date(startEnd2[["startDate"]]), 
                                                       surfaceEnd = as.Date(startEnd2[["endDate"]]),
                                                       edgeAdjust = attr(type_results, "Other")[["edgeAdjust"]],
                                                       localSample = Sample2, 
                                                       minNumObs = attr(type_results, "Other")[["minNumObs"]],
                                                       minNumUncen = attr(type_results, "Other")[["minNumUncen"]],
                                                       verbose = FALSE)),
      error = function(e) e)
    if (!inherits(possibleError3, "error") & 
        !inherits(possibleError4, "error")) {
      # note that all the flux calculations inside the bootstrap loop are in kg/day units    
      DailyRS1FD1 <- EGRET::estDailyFromSurfaces(eListBoot, 
                                                 localsurfaces = surfaces1, 
                                                 localDaily = Daily1)
      annFlex1 <- EGRET::setupYears(DailyRS1FD1, 
                                      paLong = attr(type_results, "yearPair")[["paLong"]], 
                                      paStart = attr(type_results, "yearPair")[["paStart"]])
      DailyRS2FD2 <- EGRET::estDailyFromSurfaces(eListBoot, 
                                                 localsurfaces = surfaces2, 
                                                 localDaily = Daily2)
      annFlex2 <- EGRET::setupYears(DailyRS2FD2, 
                                      paLong = attr(type_results, "yearPair")[["paLong"]], 
                                      paStart = attr(type_results, "yearPair")[["paStart"]])
    } else {
      return(list(xConc = NULL,
                  xFlux = NULL,
                  pConc = NULL,
                  pFlux = NULL))
    }

  }

  #  results are in 10^6 kg/year, when we get to the bootstrap results
  #  Converting them all to 10^6 kg/year units
  c11 <- mean(annFlex1$FNConc, na.rm = TRUE)
  f11 <- mean(annFlex1$FNFlux, na.rm = TRUE) * 0.00036525
  c22 <- mean(annFlex2$FNConc, na.rm = TRUE)
  f22 <- mean(annFlex2$FNFlux, na.rm = TRUE) * 0.00036525
  
  regDeltaConc <- type_results$x22[1] - type_results$x11[1]
  regDeltaFlux <- type_results$x22[2] - type_results$x11[2]
  LFluxDiff <- log(type_results$x22[2]) - log(type_results$x11[2])
  LConcDiff <- log(type_results$x22[1]) - log(type_results$x11[1])
  
  xConc <- (2 * regDeltaConc) - (c22 - c11)
  xFlux <- (2 * regDeltaFlux) - (f22 - f11)
  if (!is.na(xConc) & !is.na(xFlux)) {

    LConc <- (2 * LConcDiff) - (log(c22) - log(c11))
    pConc<- (100 * exp(LConc)) - 100
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

  
  return(return_list)
  
}

words <- function(z) {
  out <- if (z) 
    "Reject Ho"
  else "Do Not Reject Ho"
  return(out)
}

# bootOut <- as.data.frame(matrix(ncol = 27, nrow = 1))
# colnames(bootOut) <- c("rejectC", "pValC", "estC", "lowC90", 
#                        "upC90", "lowC50", "upC50", "lowC95", "upC95", "likeCUp", 
#                        "likeCDown", "rejectF", "pValF", "estF", "lowF90", "upF90", 
#                        "lowF50", "upF50", "lowF95", "upF95", "likeFUp", "likeFDown", 
#                        "baseConc", "baseFlux", "nBoot","startSeed","blockLength")
# 
# save(bootOut, file = "R/sysdata.rda", compress = "xz")


