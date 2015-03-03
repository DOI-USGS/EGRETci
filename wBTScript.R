library(EGRET)
library(EGRETci)
eList <- Choptank_eList
caseSetUp <- trendSetUp(eList)
eList <- setPA(eList)
eList <- setForBoot(eList)

eBoot <- wBT(eList,caseSetUp, 
             saveOutput = TRUE, fileName = "outputText.txt")

#Save output
saveEGRETci(eList, eBoot)

#TODO:
# note this thing only works if the pa is water year:
# makeTwoYearsResults