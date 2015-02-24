# this code is written to be interactive
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
library(binom)
library(survival)
source("~/Dropbox/WBT/wBTCode.R")
prob<-c(0.025, 0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95, 0.975)
numSamples <- length(eList$Sample$Date)
cat("\n Sample set runs from",eList$Sample$DecYear[1]," to",eList$Sample$DecYear[numSamples])
message("\nEnter first water year of trend period\n")
year1 <- as.numeric(readline())
message("Enter last water year of trend period\n")
year2 <- as.numeric(readline())
yearData1 <- trunc(eList$Sample$DecYear[1]+0.25)
yearData2 <- trunc(eList$Sample$DecYear[numSamples]+0.25)
nBoot <- 100  # if you want to make this flexible you can uncomment the next two lines
# message("Enter nBoot the largest number of boot replicates allowed, typically 100\n")
# nBoot <- as.numeric(readline())
cat("\nnBoot = ",nBoot," this is the maximum number of replicates that will be run\n")
message("Enter Mmin (minimum number of replicates), between 9 and nBoot, values of 39 or greater produce more accurate CIs\n")
bootBreak <- as.numeric(readline())
bootBreak <- if(bootBreak>nBoot) nBoot else bootBreak
message("Enter blockLength, in days, typically 200 is a good choice\n")
blockLength <- as.numeric(readline())
confStop <- 0.7  # testing suggests that confStop = 0.7 is good
# it is the confidence level required when checking to see if we can be confident that
#  p is really below 0.1 or really above 0.1
# message("Enter confidence interval for stopping, confStop, try 0.7\n")
# confStop <- as.numeric(readline())
message("Enter a filename for output (it will go in the working directory)\n")
fileName<-readline()
fullName<-paste(fileName,".RData",sep="")
cat("\n\n",eList$INFO$shortName,"  ",eList$INFO$paramShortName)
cat("\n\n  Bootstrap process, for change from Water Year",year1,"to Water Year",year2 )
cat("\n                   data set runs from WaterYear",yearData1, "to Water Year", yearData2)
cat("\n  Bootstrap block length in days", blockLength)
cat("\n  bootBreak is",bootBreak," confStop is",confStop)
calStart <- yearData1 - 1
countConcReject <- 0
countFluxReject <- 0
countFACReject <- 0
caseSetUp <- data.frame(year1,yearData1,year2,yearData2,numSamples,nBoot,bootBreak,blockLength,confStop)

		
		eList <- setPA(eList)
		eList <- setForBoot(eList,windowY = 7, windowQ = 2, windowS = 0.5, edgeAdjust=TRUE)
#
		eBoot <- wBT(eList,caseSetUp)
#	
localINFO <- eList$INFO
bootOut <- eBoot$bootOut
wordsOut <- eBoot$wordsOut
xConc <- eBoot$xConc
xFlux <- eBoot$xFlux									
save(caseSetUp,bootOut,wordsOut,xConc,xFlux,localINFO,file=fullName)
