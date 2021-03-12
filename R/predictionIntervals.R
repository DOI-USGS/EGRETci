#' Make Monthly Prediction Intervals
#' 
#' Month statistics using WRTDSKalman bootstrapping approach.  The input to this
#' function is the dailyBootOut matrix which contains nReplicate sets of 
#' daily flux values for the period of interest.
#' The results are in the form of quantiles of concentration and of flux 
#' for each of these months.
#' 
#' @param dailyBootOut data frame returned from \code{\link{genDailyBoot}}
#' @param eList is the data with a fitted model already done. Note that the eList$Sample 
#' may have multiple values on a given day and it can also have censored values.
#' @importFrom stats aggregate
#' @return a list of 2 data frames, one for average concentration, in mg/L
#' and one for flux (in kg/day)
#' In each data frame the first column is monthSeq that corresponds to the months
#' in the "MonthSeq" column in the eList$Daily data frame. The remaining columns are
#' quantiles of the flux or concentration (depending on the data frame).
#' @export
#' @examples 
#' eList <- EGRET::Choptank_eList
#' dailyBoot <- genDailyBoot(eList, nBoot = 2, nKalman = 2)
#' monthPcts <- makeMonthPI(dailyBoot, eList)
#' head(monthPcts[["month_flux_results"]])
#' head(monthPcts[["month_conc_results"]])
#' 
makeMonthPI <- function(dailyBootOut, eList){
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
  
  
  lastMonth <- firstMonth + length(monthFluxPcts[,1]) - 1
  monthSeqVec <- seq(firstMonth, lastMonth)
  
  month_flux_results <- data.frame(monthSeq = monthSeqVec,
                                   p1 = monthFluxPcts[,1],
                                   p2.5 = monthFluxPcts[,2],
                                   p5 = monthFluxPcts[,3],
                                   p10 = monthFluxPcts[,4],
                                   p25 = monthFluxPcts[,5],
                                   p50 = monthFluxPcts[,6],
                                   p75 = monthFluxPcts[,7],
                                   p90 = monthFluxPcts[,8],
                                   p95 = monthFluxPcts[,9],
                                   p97.5 = monthFluxPcts[,10],
                                   p99 = monthFluxPcts[,11])
  
  month_conc_results <- data.frame(monthSeq = monthSeqVec,
                                   p1 = monthConcPcts[,1],
                                   p2.5 = monthConcPcts[,2],
                                   p5 = monthConcPcts[,3],
                                   p10 = monthConcPcts[,4],
                                   p25 = monthConcPcts[,5],
                                   p50 = monthConcPcts[,6],
                                   p75 = monthConcPcts[,7],
                                   p90 = monthConcPcts[,8],
                                   p95 = monthConcPcts[,9],
                                   p97.5 = monthConcPcts[,10],
                                   p99 = monthConcPcts[,11])
  
  return(list(month_flux_results = month_flux_results,
              month_conc_results = month_conc_results))
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

