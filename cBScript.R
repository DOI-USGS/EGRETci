#     version of 2015-02-26
# this file is called cBScript.R  (which stands for "Confidence Band script")
# this scripts assumes that you have loaded a workspace 
# that workspace must contain an object called eList (which contains INFO, Daily, and Sample)
#  it doesn't matter if there is a surfaces object inside of eList
#   you will need to have installed EGRET 2.1.1 or higher
#   you will also need to change the line in the script that specifies the location of the wBTCode.R file on your computer
#   Don't move that line up higher in the script - it needs to go in that place in the workflow
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
# and have them produce pdfs of the confidence bands and you can save the bootstrap output 
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
#       Of course, if you have saved the file workspace as *****.CB.RData then you can
#       do what you want in the console, just load that RData file.  It will contain the plotConcHistBoot 
#       and plotFluxHistBoot and you can just call them with with the argument eList, or you can add other arguments as well
#
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
source("wBTCode.R")
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
