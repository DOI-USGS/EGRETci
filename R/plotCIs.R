#' plotConcHistBoot
#'
#' plotConcHistBoot
#'
#' @param eList named list
#' @param CIAnnualResults data frame from ciBands (needs nBoot, probs, and blockLength attributes)
#' @param plotFlowNorm logical
#' @param col.pred character prediction color
#' @param printTitle logical
#' @param cex.main numeric title scale
#' @export
#' @import EGRET
#' @examples
#' library(EGRET)
#' eList <- Choptank_eList
#' nBoot <- 100
#' blockLength <- 200
#' \dontrun{
#' repAnnualResults <- vector(mode = "list", length = nBoot)
#' for(n in 1:nBoot){
#'    repAnnualResults[[n]] <- bootAnnual(eList, blockLength)
#' }
#' 
#' CIAnnualResults <- ciBands(eList, repAnnualResults)
#' plotConcHistBoot(eList, CIAnnualResults)
#' }
plotConcHistBoot <- function (eList, CIAnnualResults, 
                              plotFlowNorm=TRUE, col.pred="green", printTitle=TRUE, cex.main=1.1, ...){
  
  nBoot <- attr(CIAnnualResults, "nBoot")
  blockLength <- attr(CIAnnualResults, "blockLength")
  probs <- attr(CIAnnualResults, "probs")
  widthCI <- (max(probs) - min(probs))*100
  
  localAnnualResults <- setupYears(paStart = eList$INFO$paStart, paLong = eList$INFO$paLong,
                                   localDaily = eList$Daily)
  periodName <- setSeasonLabel(localAnnualResults)
  title3 <- paste(widthCI,"% CI on FN Concentration, Replicates =",nBoot,"Block=",blockLength,"days")
  title <- ""
  if (printTitle) {
    title <- paste(eList$INFO$shortName, " ", eList$INFO$paramShortName, 
                   "\n", periodName, "\n",title3)
  }
  
  plotConcHist(eList,  col.pred=col.pred, printTitle=FALSE,...)
  title(main=title, cex.main=cex.main)
  lines(CIAnnualResults$Year, CIAnnualResults$FNConcLow,lty=2,col=col.pred)
  lines(CIAnnualResults$Year, CIAnnualResults$FNConcHigh, lty=2,col=col.pred)
  
  if (plotFlowNorm) {
    lines(CIAnnualResults$Year, CIAnnualResults$FNConcLow,lty=2,col=col.pred)
    lines(CIAnnualResults$Year, CIAnnualResults$FNConcHigh, lty=2,col=col.pred)
  }
  
}

#' plotFluxHistBoot
#'
#' plotFluxHistBoot
#'
#' @param eList named list
#' @param CIAnnualResults data frame from ciBands (needs nBoot, probs, and blockLength attributes)
#' @param fluxUnit number representing entry in pre-defined fluxUnit class array. \code{\link{printFluxUnitCheatSheet}}
#' @param plotFlowNorm logical
#' @param col.pred character prediction color
#' @param printTitle logical
#' @param cex.main numeric title scale
#' @export
#' @import EGRET
#' @examples
#' library(EGRET)
#' eList <- Choptank_eList
#' nBoot <- 100
#' blockLength <- 200
#' \dontrun{
#' repAnnualResults <- vector(mode = "list", length = nBoot)
#' for(n in 1:nBoot){
#'    repAnnualResults[[n]] <- bootAnnual(eList, blockLength)
#' }
#' 
#' CIAnnualResults <- ciBands(eList, repAnnualResults)
#' plotFluxHistBoot(eList, CIAnnualResults)
#' }
plotFluxHistBoot <- function (eList, CIAnnualResults, 
                              plotFlowNorm=TRUE, fluxUnit = 9, 
                              col.pred="green", printTitle=TRUE, cex.main=1.1, ...){
  
  localAnnualResults <- setupYears(paStart = eList$INFO$paStart, paLong = eList$INFO$paLong,
                                   localDaily = eList$Daily)
  periodName <- setSeasonLabel(localAnnualResults)
  title3 <- paste(widthCI,"% CI on FN Flux, Replicates =",nBoot,", Block=",blockLength,"days")
  title <- ""
  if (printTitle) {
    title <- paste(eList$INFO$shortName, " ", eList$INFO$paramShortName, 
                   "\n", periodName, "\n",title3)
  }
  
  if (is.numeric(fluxUnit)) {
    fluxUnit <- fluxConst[shortCode = fluxUnit][[1]]
  } else if (is.character(fluxUnit)) {
    fluxUnit <- fluxConst[fluxUnit][[1]]
  }
  unitFactorReturn <- fluxUnit@unitFactor
  
  plotFluxHist(eList,  col.pred=col.pred, printTitle=FALSE,...)
  title(main=title, cex.main=cex.main)
  lines(CIAnnualResults$Year, CIAnnualResults$FNFluxLow,lty=2,col=col.pred)
  lines(CIAnnualResults$Year, CIAnnualResults$FNFluxHigh, lty=2,col=col.pred)
  
  if (plotFlowNorm) {
    lines(CIAnnualResults$Year, CIAnnualResults$FNFluxLow * unitFactorReturn,lty=2,col=col.pred)
    lines(CIAnnualResults$Year, CIAnnualResults$FNFluxHigh * unitFactorReturn, lty=2,col=col.pred)
  }
  
}


#' saveCB
#'
#' saveCB
#'
#' @param eList named list with at least the Daily, Sample, and INFO dataframes
#' @export
#' @import EGRET
#' @examples
#' library(EGRET)
#' eList <- Choptank_eList
#' \dontrun{
#' saveCB(eList)
#' }
saveCB<-function(eList){ 
  INFO <-EGRET::getInfo(eList)
  saveName <- paste0(INFO$staAbbrev,".",INFO$constitAbbrev,".CB.RData")
  save.image(file = saveName)
  message("Saved to: ",getwd(),"/",saveName)
}

#' bootAnnual
#'
#' bootAnnual One bootstrap run.
#'
#' @param eList named list with at least the Daily, Sample, and INFO dataframes
#' @param blockLength integer suggested value is 200
#' @export
#' @import EGRET
#' @examples
#' library(EGRET)
#' eList <- Choptank_eList
#' \dontrun{
#' annualResults <- bootAnnual(eList)
#' }
bootAnnual <- function(eList, blockLength=200){
  Sample <- eList$Sample
  Daily <- eList$Daily
  INFO <- eList$INFO
  
  bootSample <- EGRETci::blockSample(Sample, blockLength)
  eListBoot <- EGRET::as.egret(INFO,Daily,bootSample,NA)
  surfaces1<-EGRET::estSurfaces(eListBoot)
  eListBoot<-EGRET::as.egret(INFO,Daily,bootSample,surfaces1)
  Daily1<-EGRET::estDailyFromSurfaces(eListBoot)
  annualResults1 <- EGRET::setupYears(Daily1)
  annualResults1$year <- as.integer(annualResults1$DecYear)
  annualResults <- annualResults1[,c("year","FNConc","FNFlux")]
  
  attr(annualResults, "blockLength") <- blockLength
  return(annualResults)
}

#' ciBands
#'
#' ciBands
#'
#' @param repAnnualResults named list returned from bootstrapping process
#' @param AnnualResults data frame
#' @param probs vector high and low confidence interval percentages
#' @export
#' @import EGRET
#' @examples
#' library(EGRET)
#' eList <- Choptank_eList
#' nBoot <- 100
#' blockLength <- 200
#' \dontrun{
#' AnnualResults <- setupYears(eList$Daily)
#' repAnnualResults <- vector(mode = "list", length = nBoot)
#' for(n = 1:nBoot){
#'    annualResults <- bootAnnual(eList, blockLength) 
#'    repAnnualResults[[n]] <- bootAnnual(eList, blockLength)
#' }
#' 
#' CIAnnualResults <- ciBands(eList, repAnnualResults)
#' }
ciBands <- function(eList, repAnnualResults, probs=c(0.05,0.95)){

  if(length(probs) != 2){
    stop("Please provide only lower and upper limit in the probs argument")
  }

  AnnualResults <- setupYears(eList$Daily)
  
  nBoot <- length(repAnnualResults)
  numYears <- nrow(repAnnualResults[[1]])
  yearStart <- repAnnualResults[[1]][1,1]
  blockLength <- attr(repAnnualResults[[1]], "blockLength")
  
  manyAnnualResults <- array(NA, dim=c(numYears,2,nBoot))
  for (i in 1:nBoot){
    manyAnnualResults[,1,i] <- 2*AnnualResults$FNConc - repAnnualResults[[i]]$FNConc
    manyAnnualResults[,2,i] <- 2*AnnualResults$FNFlux - repAnnualResults[[i]]$FNFlux
  }
  
  CIAnnualResults <- data.frame(matrix(ncol = 5, nrow = numYears))
  names(CIAnnualResults) <- c("Year","FNConcLow","FNConcHigh","FNFluxLow","FNFluxHigh")
  
  for(iYear in 1:numYears) {
    quantConc <- quantile(manyAnnualResults[iYear,1,1:nBoot],prob=probs,type=6)
    quantFlux <- quantile(manyAnnualResults[iYear,2,1:nBoot],prob=probs,type=6)
    
    CIAnnualResults$Year[iYear] <- yearStart + iYear - 1
    CIAnnualResults$FNConcLow[iYear] <- quantConc[1]
    CIAnnualResults$FNConcHigh[iYear] <- quantConc[2]
    CIAnnualResults$FNFluxLow[iYear] <- quantFlux[1]
    CIAnnualResults$FNFluxHigh[iYear] <- quantFlux[2]
  }
  
  attr(CIAnnualResults, "nBoot") <- nBoot
  attr(CIAnnualResults, "probs") <- probs
  attr(CIAnnualResults, "blockLength") <- blockLength
  
  return(CIAnnualResults)
}
