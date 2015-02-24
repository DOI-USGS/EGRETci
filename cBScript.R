# this file is called cBScript.R  (which stands for "Confidence Band script")
# this scripts assumes that you have loaded a workspace 
# that workspace must contain an object called eList (which contains INFO, Daily, and Sample)
#  it doesn't matter if there is a surfaces object inside of eList
#   you will need to have installed EGRET 2.1.1 or higher
#   if this is the first time you have done this you may need to install these other packages:
#      foreach, doParallel, iterators 
# it is best to run this in the terminal, If you run it from the console
#  it will make it impossible to do anything else with R while it is running
#  and will slow down your computer for other activities as well
#  it can take a long time to run (say 20 minutes)
#  the main thing that determines how low it will take is the number of replicates (it will prompt you for a choice)
#  I've found that even with 30 replicates I get pretty reasonable looking outputs but you might want to
#  go up to something like 120 if you really want to nail things down very well
# 
# After it has run, then you can execute the two graphical functions
# and have them produce pdfs of the confidence regions and you can save the bootstrap output 
#  you might want to save it if you are thinking you might want to redo the graphics later on
#
#        ##############  Workflow  ######################
#    total workflow inside a Terminal window might look like this 
#    (by the way on the Mac we say Terminal window, I'm not sure what this is called on a Windows computer)
#
#     Once R has started up, do setwd(".....") that sets up the working directory for all this work
#      load("... .RData")  this is the workspace in which the eList resides.  
#     note that if it is in the working directory you called for in setwd, you just need the file name
#     then do source("cBScript.R")  or full path name if the file isn't in this directory
#     respond to the prompts you will get
#
#      then there will be a long time when nothing happens, just be patient (could be 20 minutes)
#   When it is done, you may get a warning message, don't worry about it.
#   You might want to save the results as a workspace.  This command will do it (saving in the working directory)
#
#   saveCB(eList)
#
#    the saved file name will be just like the regular workspace file names but will have ".CB." before "RData"
#
#      Now set up to do two plots  -- note xxx is just a name to help you track the file
#       in the calls to the functions plotConcHistBoot and plotFluxHistBoot you can add any of the
#       other arguments that are specified in the plotConcHist or plotFluxHist functions respectively
#
#      pdf("plotConcHistBootxxx.pdf",height=6,width=8)
#      plotConcHistBoot(eList)
#      dev.off()
#      pdf("plotFluxHistBootxxx.pdf",height=6,width=8)
#      plotFluxHistBoot(eList)
#      dev.off()
#         Now you can go to the directory and open them up and view them.  
#         If you want to adjust some arguments of the plot routines or the graphics sizes
#         Just go back to the terminal and do the three lines of commands adding in whatever
#             arguments you want 
#         It will overwrite the older versions of the graphic
#
#  If you have saved the output and want to work on the figures later, you can do this:
#   copy the two functions shown near the end of this document, plotConcHistBoot, and plotFluxHistBoot into your console
#   do a <carriage return> and then load the workspace that you created (the one with .CB. in it)
#   then just give the commands plotConcHistBoot(eList) and plotFluxHistBoot(eList) to create your plots
#
#
library(foreach)
library(doParallel)
library(iterators)
library(EGRET)
# at this point make sure that the data workspace is loaded
rm(list=ls()[!(ls() %in% c("eList"))])
# will do the regular model estimation, just to make sure it
#   corresponds to the bootstrap estimates that are about to run
#
#   for this next line you will need to change the full pathname for your file structure
source("~/Dropbox/WBT/wBTCode.R")
#     
eList <- modelEstimation(eList)
surfaces <- eList$surfaces
message("\nEnter nBoot, the number of bootstrap replicates, at least 30, 100 is a good number")
nBoot <- as.numeric(readline())
message("\nEnter blockLength (in days), suggested value is 200")
blockLength <- as.numeric(readline())
message("\nEnter Confidence Interval percentage (suggested entry is 90)")
widthCI <- as.numeric(readline())
message("\nIf you want one core to stay out of the parallel processing enter 1, otherwise, 0")
coreOut <- as.numeric(readline())
#
INFO <- eList$INFO
paStart <- 10
paLong <- 12
ciLower <- (50-(widthCI/2))/100
ciUpper <- (50+(widthCI/2))/100
INFO$paStart <- paStart
INFO$paLong <- paLong
INFO$nBoot <- nBoot
INFO$blockLength <- blockLength
INFO$widthCI <- widthCI
INFO$windowY <- 7
eList <- as.egret(INFO,eList$Daily,eList$Sample,eList$surfaces)

probs <- c(ciLower,0.5,ciUpper)

windowY <- INFO$windowY
windowQ <- INFO$windowQ
windowS <- INFO$windowS
minNumObs <- INFO$minNumObs
minNumUncen <- INFO$minNumUncen
AnnualResults <- setupYears(eList$Daily, paLong=12, paStart=10)
numYears <- length(AnnualResults$DecYear)
yearStart <- AnnualResults$DecYear[1]
Daily <- eList$Daily
Sample <- eList$Sample
INFO <- eList$INFO
surfaces <- eList$surfaces

nCores <- detectCores() - coreOut
cl <- makeCluster(nCores)
registerDoParallel(cl)

repAnnualResults <- foreach(n = 1:nBoot,.packages=c('EGRET'),
                   .export="blockSample") %dopar% {
  bootSample <- blockSample(Sample, INFO$blockLength)
  eListBoot <- as.egret(INFO,Daily,bootSample,NA)
  surfaces1<-estSurfaces(eListBoot)
  eListBoot<-as.egret(INFO,Daily,bootSample,surfaces1)
  Daily1<-estDailyFromSurfaces(eListBoot)
  annualResults1 <- setupYears(Daily1)
  annualResults1$year <- as.integer(annualResults1$DecYear)
  return(annualResults1[,c("year","FNConc","FNFlux")])                
}

stopCluster(cl)               



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
  CIAnnualResults$FNConcHigh[iYear] <- quantConc[3]
  CIAnnualResults$FNFluxLow[iYear] <- quantFlux[1]
  CIAnnualResults$FNFluxHigh[iYear] <- quantFlux[3]
}
#
#
plotConcHistBoot <- function (eList, yearStart = NA, yearEnd = NA, 
    concMax = NA, printTitle = TRUE, tinyPlot = FALSE, plotFlowNorm = TRUE, 
    cex = 0.8, cex.axis = 1.1, cex.main = 1.1, lwd = 2, col = "black", 
    col.pred = "green", customPar = FALSE, ...) 
{
    localDaily <- eList$Daily
    localINFO <- eList$INFO
    if (sum(c("paStart", "paLong") %in% names(localINFO)) == 
        2) {
        paLong <- localINFO$paLong
        paStart <- localINFO$paStart
    }
    else {
        paLong <- 12
        paStart <- 10
    }
    localAnnualResults <- setupYears(paStart = paStart, paLong = paLong, 
        localDaily = localDaily)
    periodName <- setSeasonLabel(localAnnualResults = localAnnualResults)
    title3 <- paste(widthCI,"% CI on FN Concentration, Replicates =",nBoot,"Block=",blockLength,"days")
    title <- if (printTitle) 
        paste(localINFO$shortName, " ", localINFO$paramShortName, 
            "\n", periodName, "\n",title3)
    else ""
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
    if (plotFlowNorm) 
        with(localAnnualResults, lines(DecYear[DecYear > xInfo$bottom & 
            DecYear < xInfo$top], FNConc[DecYear > xInfo$bottom & 
            DecYear < xInfo$top], col = col.pred, lwd = lwd))
            lines(CIAnnualResults$Year, CIAnnualResults$FNConcLow,lty=2,col="green")
            lines(CIAnnualResults$Year, CIAnnualResults$FNConcHigh, lty=2,col="green")

}

#
#
plotFluxHistBoot <- function (eList, yearStart = NA, yearEnd = NA, fluxUnit = 9, fluxMax = NA, printTitle = TRUE, plotFlowNorm = TRUE, 
    tinyPlot = FALSE, col = "black", col.pred = "green", cex = 0.8, 
    cex.axis = 1.1, cex.main = 1.1, lwd = 2, customPar = FALSE, 
    ...) 
{
    localDaily <- eList$Daily
    localINFO <- eList$INFO
    if (sum(c("paStart", "paLong") %in% names(localINFO)) == 
        2) {
        paLong <- localINFO$paLong
        paStart <- localINFO$paStart
    }
    else {
        paLong <- 12
        paStart <- 10
    }
    localAnnualResults <- setupYears(paStart = paStart, paLong = paLong, 
        localDaily = localDaily)
    if (is.numeric(fluxUnit)) {
        fluxUnit <- fluxConst[shortCode = fluxUnit][[1]]
    }
    else if (is.character(fluxUnit)) {
        fluxUnit <- fluxConst[fluxUnit][[1]]
    }
    unitFactorReturn <- fluxUnit@unitFactor
    ylabel <- paste("Flux in ", fluxUnit@unitName, sep = "")
    numYears <- length(localAnnualResults$DecYear)
    yearStart <- if (is.na(yearStart)) 
        trunc(localAnnualResults$DecYear[1])
    else yearStart
    yearEnd <- if (is.na(yearEnd)) 
        trunc(localAnnualResults$DecYear[numYears]) + 1
    else yearEnd
    subAnnualResults <- subset(localAnnualResults, DecYear >= 
        yearStart)
    subAnnualResults <- subset(subAnnualResults, DecYear <= yearEnd)
    annFlux <- unitFactorReturn * subAnnualResults$Flux
    fnFlux <- unitFactorReturn * subAnnualResults$FNFlux
    periodName <- setSeasonLabel(localAnnualResults = localAnnualResults)
    title3 <- paste(widthCI,"% CI on FN Flux, Replicates =",nBoot,", Block=",blockLength,"days")
    title <- if (printTitle) 
        paste(localINFO$shortName, " ", localINFO$paramShortName, 
            "\n", periodName, "\n",title3)
    else ""
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
    if (plotFlowNorm) 
        lines(subAnnualResults$DecYear, fnFlux, col = col.pred, 
            lwd = lwd)
            lines(CIAnnualResults$Year, CIAnnualResults$FNFluxLow * unitFactorReturn,lty=2,col="green")
            lines(CIAnnualResults$Year, CIAnnualResults$FNFluxHigh * unitFactorReturn, lty=2,col="green")
}
#
#
saveCB<-function(eList)
{ 
	INFO <-getInfo(eList)
	saveName <- paste(INFO$staAbbrev,".",INFO$constitAbbrev,".CB.RData",sep = "")
	save.image(file = saveName)
	}