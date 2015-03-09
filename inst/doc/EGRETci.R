## ----openLibrary, echo=FALSE, message=FALSE---------------
library(xtable)
options(continue=" ")
options(width=60)
library(knitr)
library(EGRET)
library(EGRETci)

## ----include=TRUE ,echo=FALSE,eval=TRUE-------------------
opts_chunk$set(highlight=TRUE, tidy=TRUE, keep.space=TRUE, keep.blank.space=FALSE, keep.comment=TRUE, concordance=TRUE,tidy=FALSE,comment="")

knit_hooks$set(inline = function(x) {
   if (is.numeric(x)) round(x, 3)})
knit_hooks$set(crop = hook_pdfcrop)

bold.colHeaders <- function(x) {
  x <- gsub("\\^(\\d)","$\\^\\1$",x)
  x <- gsub("\\%","\\\\%",x)
  x <- gsub("\\_"," ",x)
  returnX <- paste("\\multicolumn{1}{c}{\\textbf{\\textsf{", x, "}}}", sep = "")
}

addSpace <- function(x) ifelse(x != "1", "[5pt]","")


## ----workflowFlowHistory, echo=TRUE,eval=FALSE------------
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

## ----histExample, echo=TRUE, eval=TRUE, fig.cap="Concentration histogram"----
library(EGRET)
library(EGRETci)

# Example data included in package:
eList <- Choptank_eList # Example data from EGRET package
eBoot <- Choptank_eBoot
caseSetUp <- Choptank_caseSetUp

#Concentration:
plotHistogramTrend(eBoot, caseSetUp, eList, 
                   flux=FALSE, xSeq = seq(-20,60,5))


## ----histExampleFlux, echo=TRUE, eval=TRUE, fig.cap="Flux histogram"----
#Flux
plotHistogramTrend(eBoot, caseSetUp, eList,
                   flux=TRUE, xSeq = seq(-20,60,5))

## ----histExampleCombo, echo=TRUE, eval=TRUE, fig.cap="Combo example",  fig.width=7, fig.height=4----
par(mfrow=c(1,2))
plotHistogramTrend(eBoot, caseSetUp, eList, flux=FALSE,
                   printTitle=FALSE, ylim=c(0,0.07))
plotHistogramTrend(eBoot, caseSetUp, eList, flux=TRUE,
                   printTitle=FALSE, ylim=c(0,0.07))

## ----workflowBaseR, echo=TRUE, eval=FALSE-----------------
#  library(EGRET)
#  library(EGRETci)
#  
#  eList <- Choptank_eList
#  
#  CIAnnualResults <- ciCalculations(eList)
#  
#  save(CIAnnualResults, file="CIAnnualResults.RData")
#  

## ----workflowFlowHistoryBands, echo=TRUE,eval=FALSE-------
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
#  coreOut <- 2 #Number of cores to leave out of processing tasks
#  
#  widthCI <- 90
#  ciLower <- (50-(widthCI/2))/100
#  ciUpper <- (50+(widthCI/2))/100
#  probs <- c(ciLower,ciUpper)
#  
#  nCores <- detectCores() - coreOut
#  cl <- makeCluster(nCores)
#  registerDoParallel(cl)
#  repAnnualResults <- foreach(n = 1:nBoot,.packages=c('EGRETci')) %dopar% {
#     annualResults <- bootAnnual(eList, blockLength)
#  }
#  stopCluster(cl)
#  
#  # save(repAnnualResults, file="repAnnualResults.RData")
#  
#  CIAnnualResults <- ciBands(eList, repAnnualResults, probs)
#  save(CIAnnualResults, file="CIAnnualResults.RData")

## ----plots1, echo=TRUE,eval=TRUE, fig.cap="plotConcHistBoot"----
#Load package data:
eList <- Choptank_eList
CIAnnualResults <- Choptank_CIAnnualResults

plotConcHistBoot(eList, CIAnnualResults)

## ----plots2, echo=TRUE,eval=TRUE, fig.cap="plotFluxHistBoot"----

plotFluxHistBoot(eList, CIAnnualResults)

