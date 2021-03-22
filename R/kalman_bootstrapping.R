#' WRTDSKalman Bootstrapping
#' 
#' Function to get multiple bootstrap replicates at a daily time step using the WRTDS_K
#' method.  It is done by doing bootstrap resampling of the original Sample data frame.
#' The number of these replicate samples that are created is called nBoot and in each
#' case the WRTDS model is estimated.  Then, for each of these models, there are
#' nKalman time series of daily values computed, using all of the sample values in
#' the original Sample data frame.  The total number of replicates of the complete
#' process is nBoot * nKalman.  For example we might generate 500 replicates by setting
#' nBoot = 20 and nKalman = 25.
#' 
#' @param eList is the data with a fitted model already done. Note that the eList$Sample 
#' may have multiple values on a given day and it can also have censored values.
#' @param nBoot number of times the bootstrap resampling and model estimating is done.
#' @param nKalman number of different realizations of the daily time series for each re-estimated model.
#' @param rho numeric the lag one autocorrelation. Default is 0.9.
#' @param setSeed value. Defaults is \code{NA}, which will not specify a randomized seed.
#' This can be used to make repeatable output.
#' @param jitterOn logical, if TRUE, adds "jitter" to the data in an attempt to avoid some numerical problems.  Default = FALSE.  See Details below.
#' @param V numeric a multiplier for addition of jitter to the data, default = 0.2.  See Details below.  
#' @export
#' @details
#' In some situations numerical problems are encountered in the bootstrap process, resulting in highly unreasonable spikes in the confidence intervals.
#' The use of "jitter" can often prevent these problems, but should only be used when it is clearly needed.
#' It adds a small amount of random "jitter" to the explanatory variables of the WRTDS model.  The V parameter sets the scale of variation in the log discharge values.
#' The standard deviation of the added jitter is V * standard deviation of Log Q.
#' The default for V is 0.2.  Larger values should generally be avoided, and smaller values may be sufficient.
#' @return
#' dailyBootOut a matrix of daily flux values (in kg/day).  
#' The number of columns of the matrix is the number of replicates produced
#' which is nBoot * nKalman
#' The number of rows is the number of days in the record.  
#' The set of days simulated is the same set of days that are in the eList$Daily data frame. 
#' @rdname kalman
#' @examples 
#' 
#' eList <- EGRET::Choptank_eList
#' \donttest{
#' dailyBootOut <- genDailyBoot(eList, nBoot = 2, nKalman = 2)
#' }
#' 
genDailyBoot <- function(eList, nBoot = 10, nKalman = 10, 
                         rho = 0.9, setSeed = NA,
                         jitterOn = FALSE, V = 0.2) {
  
  # to see if the results are roughly reproducible, seed value can be changed
  #
  nTotalReps <- nBoot * nKalman
  localDaily <- eList$Daily
  localSample <- eList$Sample
  localINFO <- eList$INFO
  nDaily <- length(localDaily$Date)
  if(!is.na(setSeed)){
    set.seed(setSeed)
  }
  
  if(!all(c("windowY", "windowQ", "windowS",
           "minNumObs", "minNumUncen") %in%
         names(localINFO))){
    stop("Run EGRET::setUpEstimation on eList before running genDailyBoot")
  }
  
  dailyBootOut <- matrix(data = NA, nrow = nDaily, ncol = nTotalReps)
  for(iBoot in 1: nBoot){
    message("Boot: ", iBoot)
    bootSample <- blockSample(localSample, 200)
    
    if(jitterOn) bootSample <- jitterSam(bootSample, V = V)
    
    eListBoot <- EGRET::as.egret(localINFO, localDaily, bootSample)
    surfaces1 <- EGRET::estSurfaces(eListBoot, verbose = FALSE,
                                    windowY = localINFO$windowY, 
                                    windowQ = localINFO$windowQ, 
                                    windowS = localINFO$windowS, 
                                    minNumObs = localINFO$minNumObs, 
                                    minNumUncen = localINFO$minNumUncen, 
                                    edgeAdjust = ifelse(is.null(localINFO$edgeAdjust),
                                                        TRUE, localINFO$edgeAdjust))
    eListBoot <- EGRET::as.egret(localINFO, localDaily, localSample, 
                                 surfaces1)
    message("made surfaces from boot sample", iBoot,"replicate")
    # note that we use the surfaces object made with the bootstrap sample
    # but the Kalman filtering is done with the full Sample set
    for(iKalman in 1:nKalman) {
      message("     - Kalman index: ", iKalman)
      eListK <- EGRET::WRTDSKalman(eListBoot, rho = rho, verbose = FALSE,
                                   niter = 1, seed = setSeed + iKalman)
      
      iter <- ((iBoot-1) * nKalman) + iKalman 
      dailyBootOut[,iter] <- eListK$Daily$GenFlux
    }
  }
  
  return(dailyBootOut)
}


