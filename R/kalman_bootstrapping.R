#' WRTDSKalman Bootstrapping
#' 
#' Funtion to get
#' @param eList is the data with a fitted model already done. Note that the eList$Sample 
#' may have multiple values on a given day and it can also have censored values.
#' @param nBoot number of times the bootstrap resampling and model estimating is done
#' @param nKalman number of different realizations for each re-estimated model
#' @param rho numeric the lag one autocorrelation. Default is 0.9.
#' @param niter number of iterations. Default is 200. 
#' @param setSeed value. Defaults is \code{NA}, which will not specify a randomized seed.
#' This can be used to make repeatable output.
#' @export
#' @rdname kalman
#' @examples 
#' 
#' eList <- EGRET::Choptank_eList
#' dailyBoot <- genDailyBoot(eList, nBoot = 2, nKalman = 2)
#' 
#' 
genDailyBoot <- function(eList, nBoot = 10, nKalman = 10, 
                         rho = 0.9, niter = 200, setSeed = NA) {

  # to see if the results are roughly reproducible,seed value can be changed
  #
  nTotalReps <- nBoot * nKalman
  localDaily <- eList$Daily
  localSample <- eList$Sample
  localINFO <- eList$INFO
  nDaily <- length(localDaily$Date)
  if(!is.na(setSeed)){
    set.seed(setSeed)
  }
  dailyBootOut <- matrix(data = NA, nrow = nDaily, ncol = nTotalReps)
  for(iBoot in 1: nBoot){
    cat("Boot: ", iBoot, "\n")
    bootSample <- blockSample(localSample, 200)
    eListBoot <- EGRET::as.egret(localINFO, localDaily, bootSample)
    surfaces1 <- EGRET::estSurfaces(eListBoot, verbose = FALSE)
    eListBoot <- EGRET::as.egret(localINFO, localDaily, localSample, 
                                 surfaces1)
    cat("\n made surfaces from boot sample", iBoot,"replicate")
    # note that we use the surfaces object made with the bootstrap sample
    # but the Kalman filtering is done with the full Sample set
    for(iKalman in 1:nKalman) {
      cat("     - Kalman index: ", iKalman, "\n")
      eListK <- EGRET::WRTDSKalman(eListBoot, rho = rho, niter = niter)
      
      iter <- ((iBoot-1) * nKalman) + iKalman 
      dailyBootOut[,iter] <- eListK$Daily$GenFlux
    }
  }
  
  return(dailyBootOut)
}


#' make Months
#' 
#' Month statistics
#' 
#' @param dailyBootOut data frame returned from \code{\link{genDailyBoot}}
#' @param eList is the data with a fitted model already done. Note that the eList$Sample 
#' may have multiple values on a given day and it can also have censored values.
#' @importFrom stats aggregate
#' @export
#' @examples 
#' eList <- EGRET::Choptank_eList
#' dailyBoot <- genDailyBoot(eList, nBoot = 2, nKalman = 2)
#' monthPcts <- makeMonths(dailyBoot, eList)
makeMonths <- function(dailyBootOut, eList){
  nDaily <- length(dailyBootOut[,1])
  nIter <- length(dailyBootOut[1,])
  firstMonth <- eList$Daily$MonthSeq[1]
  lastMonth <- eList$Daily$MonthSeq[nDaily]
  nMonths <- lastMonth - firstMonth + 1
  monthsConc <- matrix(data = NA, nrow = nMonths, ncol = nIter)
  attr(monthsConc, "firstMonth") <- firstMonth
  monthsFlux <- monthsConc
  attr(monthsFlux, "firstMonth") <- firstMonth
  for(iter in 1:nIter){
    monthsFlux[,iter] <- aggregate(dailyBootOut[,iter], by = list(eList$Daily$MonthSeq), mean)$x
    monthsConc[,iter] <- aggregate(dailyBootOut[,iter]/(eList$Daily$Q * 86.4), by = list(eList$Daily$MonthSeq), mean)$x
  }
  # now do the percentiles by month
  probs <- c(0.01, 0.025, 0.05, 0.1, 0.25, 0.5, 0.75, 0.90, 0.95, 0.975, 0.99)
  monthFluxPcts <- matrix(data = NA, nrow = nMonths, ncol = 11)
  monthConcPcts <- matrix(data = NA, nrow = nMonths, ncol = 11)
  for(iMonth in 1:nMonths){
    monthFluxPcts[iMonth,] <- quantile(monthsFlux[iMonth,], probs = probs, type = 6, na.rm = TRUE)
    monthConcPcts[iMonth,] <- quantile(monthsConc[iMonth,], probs = probs, type = 6, na.rm = TRUE)
  }
  monthPcts <- list(monthFluxPcts = monthFluxPcts, 
                    monthConcPcts = monthConcPcts,
                    firstMonth = firstMonth)
  return(monthPcts)
}


#' monthSeqToDec
#' 
#' Convert a sequence of month integers into their decimal years.
#' 
#' @export
#' @param monthSeq integer vector of months. Month 1 is considered Jan. 1850.
#' @examples 
#' 
#' months <- 1558:1600
#' monthSeqToDec(months)
monthSeqToDec <- function(monthSeq){
  # function that turns a monthSeq value into a decYear value
  dec <- floor((monthSeq - 1)/12) + 
    1850 + (1/24) +
    ((monthSeq-1)%%12/12)
  return(dec)
}

