## ----setup, include=FALSE---------------------------------
library(xtable)
options(continue=" ")
options(width=60)
library(knitr)
library(EGRET)
library(EGRETci)


## ----eval=FALSE, echo=TRUE--------------------------------
#  library(EGRET)
#  library(EGRETci)
#  eList <- Choptank_eList
#  
#  #Interactive function to set up trend analysis:
#  caseSetUp <- trendSetUp(eList)
#  
#  #Currently, only water-year calculations are supported
#  eList <- setPA(eList)
#  eList <- setForBoot(eList)
#  
#  eBoot <- wBT(eList,caseSetUp,
#               saveOutput = TRUE, fileName = "outputText.txt")
#  
#  #Interactive save output function:
#  saveEGRETci(eList, eBoot)
#  
#  

## ----, fig.height=6, fig.width=6--------------------------
library(EGRET)
library(EGRETci)

# Example data included in package:
eList <- Choptank_eList # Example data from EGRET package
eBoot <- Choptank_eBoot
caseSetUp <- Choptank_caseSetUp

#Concentration:
plotHistogramTrend(eBoot, caseSetUp, eList, 
                   flux=FALSE, xSeq = seq(-20,60,5))

#Flux
plotHistogramTrend(eBoot, caseSetUp, eList,
                   flux=TRUE, xSeq = seq(-20,60,5))


## ----, histExampleCombo, fig.width=7, fig.height=4--------
par(mfrow=c(1,2))
plotHistogramTrend(eBoot, caseSetUp, eList, flux=FALSE,
                   printTitle=FALSE, ylim=c(0,0.07), xSeq=seq(-10,70,10))
plotHistogramTrend(eBoot, caseSetUp, eList, flux=TRUE,
                   printTitle=FALSE, ylim=c(0,0.07), xSeq=seq(-10,70,10))
par(mfrow=c(1,1))                   

## ----, eval=FALSE-----------------------------------------
#  library(EGRET)
#  library(EGRETci)
#  
#  eList <- Choptank_eList
#  
#  CIAnnualResults <- ciCalculations(eList)
#  
#  save(eList,CIAnnualResults, file="CIAnnualResults.RData")
#  
#  

## ----, eval=FALSE-----------------------------------------
#  library(foreach)
#  library(doParallel)
#  library(iterators)
#  library(EGRET)
#  library(EGRETci)
#  
#  eList <- Choptank_eList
#  
#  nBoot <- 100
#  blockLength <- 200
#  coreOut <- 1 #Number of cores to leave out of processing tasks
#  
#  widthCI <- 90
#  ciLower <- (50-(widthCI/2))/100
#  ciUpper <- (50+(widthCI/2))/100
#  probs <- c(ciLower,ciUpper)
#  
#  nCores <- detectCores() - coreOut
#  cl <- makeCluster(nCores)
#  registerDoParallel(cl)
#  repAnnual <- foreach(n = 1:nBoot,.packages=c('EGRETci')) %dopar% {
#     annualResults <- bootAnnual(eList, blockLength)
#  }
#  stopCluster(cl)
#  
#  # save(repAnnualResults, file="repAnnual.RData")
#  
#  CIAnnualResults <- ciBands(eList, repAnnual, probs)
#  save(eList,CIAnnualResults, file="CIAnnualResults.RData")
#  

## ----, fig.height=5, fig.width=7--------------------------
eList <- Choptank_eList

CIAnnualResults <- Choptank_CIAnnualResults

plotConcHistBoot(eList, CIAnnualResults)

plotFluxHistBoot(eList, CIAnnualResults)


