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
#  #Interactive function to set up trend analysis:
#  caseSetUp <- trendSetUp(eList)
#  eList <- setPA(eList)
#  eList <- setForBoot(eList)
#  
#  eBoot <- wBT(eList,caseSetUp,
#               saveOutput = TRUE, fileName = "outputText.txt")
#  
#  #Save output
#  saveEGRETci(eList, eBoot)

## ----workflowBaseR, echo=TRUE, eval=FALSE-----------------
#  
#  eList <- Choptank_eList
#  
#  nBoot <- 100
#  blockLength <- 200
#  
#  widthCI <- 90
#  ciLower <- (50-(widthCI/2))/100
#  ciUpper <- (50+(widthCI/2))/100
#  probs <- c(ciLower,ciUpper)
#  
#  repAnnualResults <- vector(mode = "list", length = nBoot)
#  for(n in 1:nBoot){
#      repAnnualResults[[n]] <- bootAnnual(eList, blockLength)
#  }
#  
#  CIAnnualResults <- ciBands(eList, repAnnualResults, probs)
#  
#  

## ----workflowFlowHistoryBands, echo=TRUE,eval=FALSE-------
#  library(foreach)
#  library(doParallel)
#  library(iterators)
#  
#  nBoot <- 100
#  blockLength <- 200
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
#  CIAnnualResults <- ciBands(eList, repAnnualResults, probs)
#  

