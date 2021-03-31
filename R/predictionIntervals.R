#' Make Monthly Prediction Intervals
#' 
#' Month statistics using WRTDSKalman bootstrapping approach.  The input to this
#' function is the dailyBootOut matrix which contains nReplicate sets of 
#' daily flux values for the period of interest.
#' The results are in the form of quantiles of concentration and of flux 
#' for each of these months.
#' 
#' @param dailyBootOut data frame returned from \code{\link{genDailyBoot}}
#' @param fluxUnit number representing entry in pre-defined fluxUnit class array. 
#' \code{\link[EGRET]{printFluxUnitCheatSheet}}
#' @param eList named list with at least the Daily, Sample, and INFO 
#' dataframes. Created from the EGRET package, after running \code{\link[EGRET]{modelEstimation}}.
#' @importFrom stats aggregate
#' @return a list of 2 data frames, one for average concentration, in mg/L
#' and one for flux (unit depends on fluxUnit argument)
#' In each data frame the first column is monthSeq that corresponds to the months
#' in the "MonthSeq" column in the eList$Daily data frame. The remaining columns are
#' quantiles of the flux or concentration (depending on the data frame).
#' @export
#' @examples 
#' eList <- EGRET::Choptank_eList
#' # This example is only based on 4 iterations
#' # Actual prediction intervals should be calculated on
#' # a much larger number of iterations (several hundred). 
#' dailyBoot <- Choptank_dailyBootOut
#' monthPcts <- makeMonthPI(dailyBoot, eList)
#' head(monthPcts[["flux"]])
#' head(monthPcts[["conc"]])
#' 
makeMonthPI <- function(dailyBootOut, eList, fluxUnit = 3){
  nDaily <- length(dailyBootOut[,1])
  nIter <- length(dailyBootOut[1,])
  firstMonth <- eList$Daily$MonthSeq[1]
  lastMonth <- eList$Daily$MonthSeq[nDaily]
  nMonths <- lastMonth - firstMonth + 1
  monthsConc <- matrix(data = NA, nrow = nMonths, ncol = nIter)
  attr(monthsConc, "firstMonth") <- firstMonth
  monthsFlux <- monthsConc
  attr(monthsFlux, "firstMonth") <- firstMonth
  
  if (is.numeric(fluxUnit)){
    fluxUnit <- EGRET::fluxConst[shortCode=fluxUnit][[1]]    
  } else if (is.character(fluxUnit)){
    fluxUnit <- EGRET::fluxConst[fluxUnit][[1]]
  }
  for(iter in 1:nIter){
    monthsFlux[,iter] <- aggregate(dailyBootOut[,iter], by = list(eList$Daily$MonthSeq), mean)$x * fluxUnit@unitFactor
    aveFlow <- aggregate(eList$Daily$Q, by = list(eList$Daily$MonthSeq), mean)$x
    monthsConc[,iter] <- monthsFlux[,iter] / (86.4 * aveFlow)
  }
  
  DecYear <- aggregate(eList$Daily$DecYear, by = list(eList$Daily$MonthSeq), mean)$x
  
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
                                   DecYear = DecYear,
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
                                   DecYear = DecYear,
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
  
  return(list(flux = month_flux_results,
              conc = month_conc_results))
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

#' Make Annual Prediction Intervals
#' 
#' This function takes the output from \code{\link{genDailyBoot}} and 
#' calculates the quantiles for an annual (based on paStart/paLong) aggregation. 
#' This means that the function can be used for seasons.
#' 
#' @param dailyBootOut data frame returned from \code{\link{genDailyBoot}}
#' @param fluxUnit number representing entry in pre-defined fluxUnit class array. 
#' \code{\link[EGRET]{printFluxUnitCheatSheet}}
#' @param eList named list with at least the Daily, Sample, and INFO 
#' dataframes. Created from the EGRET package, after running \code{\link[EGRET]{modelEstimation}}.
#' @param paLong numeric integer specifying the length of the period of analysis, in months, 1<=paLong<=12, default is 12
#' @param paStart numeric integer specifying the starting month for the period of analysis, 1<=paStart<=12, default is 10 
#' @importFrom stats aggregate
#' 
#' @return a list of 2 data frames, one for average concentration, in mg/L
#' and one for flux (unit depends on fluxUnit argument)
#' In each data frame the first column is DecYear. The remaining columns are
#' quantiles of the flux or concentration (depending on the data frame).
#' @export
#' @examples 
#' eList <- EGRET::Choptank_eList
#' # This example is only based on 4 iterations
#' # Actual prediction intervals should be calculated on
#' # a much larger number of iterations (several hundred).  
#' dailyBoot <- Choptank_dailyBootOut
#' annualPcts <- makeAnnualPI(dailyBoot, eList)
#' head(annualPcts[["flux"]])
#' head(annualPcts[["conc"]])
#' 
makeAnnualPI <- function(dailyBootOut, eList, 
                        paLong = 12, paStart = 10, fluxUnit = 3){
  
  nIter <- length(dailyBootOut[1,])

  localDaily <- eList$Daily
  Annual_results <- data.frame()
  
  if (is.numeric(fluxUnit)){
    fluxUnit <- EGRET::fluxConst[shortCode=fluxUnit][[1]]    
  } else if (is.character(fluxUnit)){
    fluxUnit <- EGRET::fluxConst[fluxUnit][[1]]
  }
  Annual_results <- data.frame()
  
  for(i in 1:nIter){
    localDaily$GenFlux <- dailyBootOut[,i]
    
    localDaily$GenConc <- localDaily$GenFlux /(localDaily$Q *86.40)
    
    AnnualResults_i <- EGRET::setupYears(localDaily, 
                                  paLong = paLong, paStart = paStart)
    AnnualResults_i <- AnnualResults_i[,c("DecYear", "GenFlux", "GenConc")]
    Annual_results <- rbind(Annual_results, AnnualResults_i)
  }
  
  probs <- c(0.01, 0.025, 0.05, 0.1, 0.25, 0.5, 0.75, 0.90, 0.95, 0.975, 0.99)
  
  years <- unique(Annual_results$DecYear)
  yearFluxPcts <- matrix(data = NA, nrow = length(years), ncol = length(probs))
  yearConcPcts <- matrix(data = NA, nrow = length(years), ncol = length(probs))
  
  for(year in years){
    i <- which(years %in% year)
    flux_conc_year <- Annual_results[Annual_results$DecYear == year,] 
    yearFluxPcts[i,] <- quantile(flux_conc_year$GenFlux, probs = probs, type = 6, na.rm = TRUE) * fluxUnit@unitFactor
    yearConcPcts[i,] <- quantile(flux_conc_year$GenConc, probs = probs, type = 6, na.rm = TRUE)
  }
  
  year_flux_results <- data.frame(DecYear = years,
                                  p1 = yearFluxPcts[,1],
                                  p2.5 = yearFluxPcts[,2],
                                  p5 = yearFluxPcts[,3],
                                  p10 = yearFluxPcts[,4],
                                  p25 = yearFluxPcts[,5],
                                  p50 = yearFluxPcts[,6],
                                  p75 = yearFluxPcts[,7],
                                  p90 = yearFluxPcts[,8],
                                  p95 = yearFluxPcts[,9],
                                  p97.5 = yearFluxPcts[,10],
                                  p99 = yearFluxPcts[,11])
  year_conc_results <- data.frame(DecYear = years,
                                  p1 = yearConcPcts[,1],
                                  p2.5 = yearConcPcts[,2],
                                  p5 = yearConcPcts[,3],
                                  p10 = yearConcPcts[,4],
                                  p25 = yearConcPcts[,5],
                                  p50 = yearConcPcts[,6],
                                  p75 = yearConcPcts[,7],
                                  p90 = yearConcPcts[,8],
                                  p95 = yearConcPcts[,9],
                                  p97.5 = yearConcPcts[,10],
                                  p99 = yearConcPcts[,11])
  return(list(flux = year_flux_results,
              conc = year_conc_results))
}


#' Make Daily Prediction Intervals
#' 
#' This function takes the output from \code{\link{genDailyBoot}} and 
#' calculates the quantiles for a daily aggregation.
#' 
#' @param dailyBootOut data frame returned from \code{\link{genDailyBoot}}
#' @param fluxUnit number representing entry in pre-defined fluxUnit class array. 
#' \code{\link[EGRET]{printFluxUnitCheatSheet}}
#' @param eList named list with at least the Daily, Sample, and INFO 
#' dataframes. Created from the EGRET package, after running \code{\link[EGRET]{modelEstimation}}.
#' @importFrom stats aggregate
#' 
#' @return a list of 2 data frames, one for average concentration, in mg/L
#' and one for flux (unit depends on fluxUnit argument)
#' In each data frame the first column is Date. The remaining columns are
#' quantiles of the flux or concentration (depending on the data frame).
#' @export
#' @examples 
#' eList <- EGRET::Choptank_eList
#' # This example is only based on 4 iterations
#' # Actual prediction intervals should be calculated on
#' # a much larger number of iterations (several hundred).
#' dailyBoot <- Choptank_dailyBootOut
#' dailyPcts <- makeDailyPI(dailyBoot, eList)
#' head(dailyPcts[["flux"]])
#' head(dailyPcts[["conc"]])
#' 
makeDailyPI <- function(dailyBootOut, eList, fluxUnit = 3){
  
  nDaily <- length(dailyBootOut[,1])
  nIter <- length(dailyBootOut[1,])
  
  if (is.numeric(fluxUnit)){
    fluxUnit <- EGRET::fluxConst[shortCode=fluxUnit][[1]]    
  } else if (is.character(fluxUnit)){
    fluxUnit <- EGRET::fluxConst[fluxUnit][[1]]
  }

  # now do the percentiles by month
  probs <- c(0.01, 0.025, 0.05, 0.1, 0.25, 0.5, 0.75, 0.90, 0.95, 0.975, 0.99)
  dailyFluxPcts <- matrix(data = NA, nrow = nDaily, ncol = 11)
  dailyConcPcts <- matrix(data = NA, nrow = nDaily, ncol = 11)
  for(iDay in 1:nDaily){
    dailyFluxPcts[iDay,] <- quantile(dailyBootOut[iDay,], probs = probs, type = 6, na.rm = TRUE)
    dailyConcPcts[iDay,] <- quantile(dailyBootOut[iDay,]/ (86.4 * eList$Daily$Q[iDay]), probs = probs, type = 6, na.rm = TRUE)
  }

  daily_flux_results <- data.frame(Date = eList$Daily$Date,
                                   DecYear = eList$Daily$DecYear,
                                   p1 = dailyFluxPcts[,1],
                                   p2.5 = dailyFluxPcts[,2],
                                   p5 = dailyFluxPcts[,3],
                                   p10 = dailyFluxPcts[,4],
                                   p25 = dailyFluxPcts[,5],
                                   p50 = dailyFluxPcts[,6],
                                   p75 = dailyFluxPcts[,7],
                                   p90 = dailyFluxPcts[,8],
                                   p95 = dailyFluxPcts[,9],
                                   p97.5 = dailyFluxPcts[,10],
                                   p99 = dailyFluxPcts[,11])
  
  daily_conc_results <- data.frame(Date = eList$Daily$Date,
                                   DecYear = eList$Daily$DecYear,
                                   p1 = dailyConcPcts[,1],
                                   p2.5 = dailyConcPcts[,2],
                                   p5 = dailyConcPcts[,3],
                                   p10 = dailyConcPcts[,4],
                                   p25 = dailyConcPcts[,5],
                                   p50 = dailyConcPcts[,6],
                                   p75 = dailyConcPcts[,7],
                                   p90 = dailyConcPcts[,8],
                                   p95 = dailyConcPcts[,9],
                                   p97.5 = dailyConcPcts[,10],
                                   p99 = dailyConcPcts[,11])
  
  return(list(flux = daily_flux_results,
              conc = daily_conc_results))
}

