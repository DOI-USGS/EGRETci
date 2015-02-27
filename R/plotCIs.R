plotConcHistBoot <- function (eList, yearStart = NA, yearEnd = NA, 
                              concMax = NA, printTitle = TRUE, tinyPlot = FALSE, plotFlowNorm = TRUE, 
                              cex = 0.8, cex.axis = 1.1, cex.main = 1.1, lwd = 2, col = "black", 
                              col.pred = "green", customPar = FALSE, ...){
  
  localDaily <- eList$Daily
  localINFO <- eList$INFO
  if (sum(c("paStart", "paLong") %in% names(localINFO)) ==  2) {
    paLong <- localINFO$paLong
    paStart <- localINFO$paStart
  } else {
    paLong <- 12
    paStart <- 10
  }
  localAnnualResults <- setupYears(paStart = paStart, paLong = paLong, 
                                   localDaily = localDaily)
  periodName <- setSeasonLabel(localAnnualResults = localAnnualResults)
  title3 <- paste(widthCI,"% CI on FN Concentration, Replicates =",nBoot,"Block=",blockLength,"days")
  title <- ""
  if (printTitle) {
    title <- paste(localINFO$shortName, " ", localINFO$paramShortName, 
                   "\n", periodName, "\n",title3)
  }
  
  xInfo <- generalAxis(x = localAnnualResults$DecYear, minVal = yearStart, 
                       maxVal = yearEnd, padPercent = 0, tinyPlot = tinyPlot)
  
  combinedY <- c(localAnnualResults$Conc, localAnnualResults$FNConc[localAnnualResults$DecYear > 
                                                                      xInfo$bottom & localAnnualResults$DecYear < xInfo$top])
  yInfo <- generalAxis(x = combinedY, minVal = 0, maxVal = concMax, 
                       padPercent = 5, tinyPlot = tinyPlot)
  
  genericEGRETDotPlot(x = localAnnualResults$DecYear, y = localAnnualResults$Conc, 
                      xTicks = xInfo$ticks, yTicks = yInfo$ticks, xDate = TRUE, 
                      xlim = c(xInfo$bottom, xInfo$top), ylim = c(yInfo$bottom, 
                                                                  yInfo$top), ylab = "Concentration in mg/L", col = col, 
                      cex = cex, plotTitle = title, cex.axis = cex.axis, cex.main = cex.main, 
                      tinyPlot = tinyPlot, customPar = customPar, ...)
  if (plotFlowNorm) {
    with(localAnnualResults, lines(DecYear[DecYear > xInfo$bottom & 
                                             DecYear < xInfo$top], FNConc[DecYear > xInfo$bottom & 
                                                                            DecYear < xInfo$top], col = col.pred, lwd = lwd))
    lines(CIAnnualResults$Year, CIAnnualResults$FNConcLow,lty=2,col="green")
    lines(CIAnnualResults$Year, CIAnnualResults$FNConcHigh, lty=2,col="green")
  }
}

#
#


plotFluxHistBoot <- function (eList, yearStart = NA, yearEnd = NA, fluxUnit = 9, fluxMax = NA, printTitle = TRUE, plotFlowNorm = TRUE, 
                              tinyPlot = FALSE, col = "black", col.pred = "green", cex = 0.8, 
                              cex.axis = 1.1, cex.main = 1.1, lwd = 2, customPar = FALSE, 
                              ...) {
  
  localDaily <- eList$Daily
  localINFO <- eList$INFO
  
  if (sum(c("paStart", "paLong") %in% names(localINFO)) == 
        2) {
    paLong <- localINFO$paLong
    paStart <- localINFO$paStart
  } else {
    paLong <- 12
    paStart <- 10
  }
  
  localAnnualResults <- setupYears(paStart = paStart, paLong = paLong, 
                                   localDaily = localDaily)
  if (is.numeric(fluxUnit)) {
    fluxUnit <- fluxConst[shortCode = fluxUnit][[1]]
  } else if (is.character(fluxUnit)) {
    fluxUnit <- fluxConst[fluxUnit][[1]]
  }
  
  unitFactorReturn <- fluxUnit@unitFactor
  ylabel <- paste("Flux in ", fluxUnit@unitName, sep = "")
  numYears <- length(localAnnualResults$DecYear)
  
  if (is.na(yearStart)) {
    yearStart <-  trunc(localAnnualResults$DecYear[1])
  }
  
  if (is.na(yearEnd)) {
    yearEnd <- trunc(localAnnualResults$DecYear[numYears]) + 1
  }
  
  subAnnualResults <- subset(localAnnualResults, DecYear >= 
                               yearStart)
  subAnnualResults <- subset(subAnnualResults, DecYear <= yearEnd)
  annFlux <- unitFactorReturn * subAnnualResults$Flux
  fnFlux <- unitFactorReturn * subAnnualResults$FNFlux
  periodName <- setSeasonLabel(localAnnualResults = localAnnualResults)
  title3 <- paste(widthCI,"% CI on FN Flux, Replicates =",nBoot,", Block=",blockLength,"days")
  
  title <- ""
  if (printTitle) {
    title <- paste(localINFO$shortName, " ", localINFO$paramShortName, 
                   "\n", periodName, "\n",title3)
  }
  
  xInfo <- generalAxis(x = subAnnualResults$DecYear, minVal = yearStart, 
                       maxVal = yearEnd, padPercent = 0, tinyPlot = tinyPlot)
  combinedY <- c(annFlux, fnFlux)
  yInfo <- generalAxis(x = combinedY, minVal = 0, maxVal = fluxMax, 
                       padPercent = 5, tinyPlot = tinyPlot)
  genericEGRETDotPlot(x = subAnnualResults$DecYear, y = annFlux, 
                      xTicks = xInfo$ticks, yTicks = yInfo$ticks, xDate = TRUE, 
                      xlim = c(xInfo$bottom, xInfo$top), ylim = c(0, yInfo$top), 
                      col = col, ylab = ylabel, plotTitle = title, customPar = customPar, 
                      cex = cex, cex.axis = cex.axis, cex.main = cex.main, 
                      tinyPlot = tinyPlot, ...)
  if (plotFlowNorm) {
    lines(subAnnualResults$DecYear, fnFlux, col = col.pred, 
          lwd = lwd)
    lines(CIAnnualResults$Year, CIAnnualResults$FNFluxLow * unitFactorReturn,lty=2,col="green")
    lines(CIAnnualResults$Year, CIAnnualResults$FNFluxHigh * unitFactorReturn, lty=2,col="green")
  }
}
#
#
saveCB<-function(eList){ 
  INFO <-getInfo(eList)
  saveName <- paste(INFO$staAbbrev,".",INFO$constitAbbrev,".CB.RData",sep = "")
  save.image(file = saveName)
}