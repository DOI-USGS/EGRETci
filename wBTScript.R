# this code is written to be interactive  (version of 2015-02-26)
# it assumes that the user has supplied an object called eList
# eList is an egret object and it contains (in this order) INFO, Daily, and Sample 
#    surfaces is not needed, but it is ok if it exists in eList
#  
#  If you have a workspace from a version of EGRET from before 2014-11-01
#       and your Sample data frame has the discharge values in it 
#  all you need to do prior to running it is give the command
#      eList<-as.egret(INFO,Daily,Sample)
#  If you are creating the data sets from scratch, you just use the standard workflow to create eList
#  That means you have to get your workflow to this point:  eList <- mergeReport(INFO,Daily,Sample)
#  There are various dependencies that you will need: EGRET, binom and survival
#  The code is written here using edgeAdjust=TRUE and with windowY=7, 
#    but you can override the defaults
#  by modifying the script near the bottom at call to function: "setForBoot"
#  There are some other parameters for the bootstrap process that are hard coded here
#  You can change them in the code below - they are nBoot (which is Mmax in the paper)
#              I've set it to 100, 
#        also confStop (which is 1 - alphap in the paper), I've set it to 0.7
#     
#  The code will produce a very verbose output, but it helps to see what is happening
#  The tail end of the output is intended to be a bit better formatted
#  when it finishes it will produce and save a set of outputs
#    they are caseSetUp, bootOut, wordsOut, xConc, xFlux, and INFO 
#    you can learn about them by reading the comments at the top of wBTCode.R
#
#      Before the first time you run this script you need to take care of two things
#     1. Make sure that you have installed from CRAN, the latest versions of three packages:
#        EGRET, binom, and survival
#     2. Edit this script so that the line that starts with the word "source"
#           shows the actual full name of the file wBTCode.R 
#    
#   
#    Each time you run this script you need to do two things once you have started R
#    1.  Give the command setwd("..name of directory..") this should be a directory where
#           you want to accumulate bootstrap results, 
#           could be the same as where you have the workspaces you will use
#    2. load the workspace that you want to work on, make sure it contains eList
#           load("...name of workspace...")
#      
#
#
library(EGRET)
library(EGRETci)
eList <- Choptank_eList
caseSetUp <- trendSetUp(eList)
eList <- setPA(eList)
eList <- setForBoot(eList)
eBoot <- wBT(eList,caseSetUp)

#Save output
saveEGRETci(eList, eBoot)
