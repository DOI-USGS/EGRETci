
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
#' plotHistogramTrend(eList, boot_group_out, caseSetUp=NA)
#' }
runGroupsBoot <- function (eList, groupResults, nBoot = 100, 
                           startSeed = 494817, blockLength = 200,
                           jitterOn = FALSE, V = 0.2,
                           run.parallel = FALSE){

  xConc <- rep(NA, nBoot)
  xFlux <- rep(NA, nBoot)
  pConc <- rep(NA, nBoot)
  pFlux <- rep(NA, nBoot)
  
  nBootGood <- 0
  
  if(run.parallel){
    `%dopar%` <- foreach::`%dopar%`
    boot_list_out <- foreach::foreach(iBoot = 1:ceiling(1.25*nBoot), 
                                      .packages=c('EGRETci', 'EGRET')) %dopar% {
                                        boot_list <- boot_group_run(iBoot, startSeed,
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
        boot_list <- boot_group_run(iBoot, startSeed,
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
  
    for (iBoot in 1:(2 * nBoot)) {
      boot_list <- boot_group_run(iBoot, startSeed,
                                  eList, groupResults, 
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
  
  if (iBoot == 2 * nBoot) {
    message(iBoot, " iterations were run. They only achieved ", 
            nBootGood, " sucessful runs.")
  } else if (iBoot > nBoot) {
    message("It took ", iBoot, " iterations to achieve ", 
            nBoot, " sucessful runs.")
  }
  
  groupBootOut <- calc_boot_out(xConc, xFlux, pConc, pFlux,
                                groupResults, nBoot, startSeed,
                                blockLength, nBootGood)
  
  attr(groupResults, "paStart") <- attr(groupResults, "groupInfo")[["paStart"]]
  attr(groupResults, "paLong") <- attr(groupResults, "groupInfo")[["paLong"]]
  
  pair_boot_message(eList, 
                    groupResults,
                    groupBootOut,
                    nBootGood, nBoot,
                    type = "group")
  
  return(groupBootOut)
}

boot_group_run <- function(iBoot, startSeed,
                           eList, groupResults, 
                           jitterOn, V,
                           blockLength){
  
  bootSample <- blockSample(localSample = eList$Sample, 
                            blockLength = blockLength, 
                            startSeed = startSeed + iBoot)
  
  if(jitterOn) bootSample <- EGRET::jitterSam(bootSample, V = V)
  
  eListBoot <- suppressMessages(EGRET::as.egret(eList$INFO,
                                                eList$Daily,
                                                bootSample,NA))
  dateInfo <- attr(groupResults, "dateInfo")
  
  if(attr(groupResults, "Other")[["wall"]]) {
    possibleError <- tryCatch(surfaces <- suppressMessages(EGRET::stitch(eListBoot, 
                                                                         surfaceStart = attr(groupResults, "SampleBlocks")[["surfaceStart"]], 
                                                                         surfaceEnd =  attr(groupResults, "SampleBlocks")[["surfaceEnd"]],
                                                                         sample1StartDate = attr(groupResults, "SampleBlocks")[["sample1StartDate"]],
                                                                         sample1EndDate = attr(groupResults, "SampleBlocks")[["sample1EndDate"]],
                                                                         sample2StartDate = attr(groupResults, "SampleBlocks")[["sample2StartDate"]], 
                                                                         sample2EndDate = attr(groupResults, "SampleBlocks")[["sample2EndDate"]],
                                                                         windowY = attr(groupResults, "Other")[["windowY"]], 
                                                                         windowQ = attr(groupResults, "Other")[["windowQ"]], 
                                                                         windowS = attr(groupResults, "Other")[["windowS"]],
                                                                         minNumObs = attr(groupResults, "Other")[["minNumObs"]], 
                                                                         minNumUncen =  attr(groupResults, "Other")[["minNumUncen"]], 
                                                                         edgeAdjust =  attr(groupResults, "Other")[["edgeAdjust"]])), 
                              error = function(e) e)
  } else {
    possibleError <- tryCatch(surfaces <- EGRET::estSurfaces(eListBoot, 
                                                             surfaceStart = attr(groupResults, "SampleBlocks")[["surfaceStart"]], 
                                                             surfaceEnd =  attr(groupResults, "SampleBlocks")[["surfaceEnd"]],
                                                             windowY = attr(groupResults, "Other")[["windowY"]],
                                                             windowQ = attr(groupResults, "Other")[["windowQ"]],
                                                             windowS = attr(groupResults, "Other")[["windowS"]],
                                                             minNumObs = attr(groupResults, "Other")[["minNumObs"]],
                                                             minNumUncen =  attr(groupResults, "Other")[["minNumUncen"]],
                                                             edgeAdjust =  attr(groupResults, "Other")[["edgeAdjust"]]),
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
    eListOut$INFO$wall <- attr(groupResults, "Other")[["wall"]]
    eListOut$INFO$surfaceStart <- attr(groupResults, "SampleBlocks")[["surfaceStart"]]
    eListOut$INFO$surfaceEnd <- attr(groupResults, "SampleBlocks")[["surfaceEnd"]]
    DailyFlex <- eListOut$Daily
    annFlex <- EGRET::setupYears(DailyFlex, 
                                 paLong =  attr(groupResults, "groupInfo")[["paLong"]],
                                 paStart =  attr(groupResults, "groupInfo")[["paStart"]])
    annFlex$year <- floor(annFlex$DecYear + (annFlex$PeriodLong / 12) * 0.5)
    annFlex1 <- annFlex[annFlex$DecYear >= attr(groupResults, "groupInfo")[["group1firstYear"]] &
                          annFlex$DecYear <= attr(groupResults, "groupInfo")[["group1lastYear"]],]
    annFlex2 <- annFlex[annFlex$DecYear >= attr(groupResults, "groupInfo")[["group2firstYear"]] & 
                          annFlex$DecYear <= attr(groupResults, "groupInfo")[["group2lastYear"]],]
    
    #  pairResults are in 10^6 kg/year, when we get to the bootstrap results
    #  Converting them all to 10^6 kg/year units
    c11 <- mean(annFlex1$FNConc, na.rm = TRUE)
    f11 <- mean(annFlex1$FNFlux, na.rm = TRUE) * 0.00036525
    c22 <- mean(annFlex2$FNConc, na.rm = TRUE)
    f22 <- mean(annFlex2$FNFlux, na.rm = TRUE) * 0.00036525
    
    regDeltaConc <- groupResults$x22[1] - groupResults$x11[1]
    regDeltaFlux <- groupResults$x22[2] - groupResults$x11[2]
    LFluxDiff <- log(groupResults$x22[2]) - log(groupResults$x11[2])
    LConcDiff <- log(groupResults$x22[1]) - log(groupResults$x11[1])
    
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
  } else {
    return_list <- list(xConc = NULL,
                        xFlux = NULL,
                        pConc = NULL,
                        pFlux = NULL)
  }
  
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


