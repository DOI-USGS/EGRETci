#   version of 2015-02-26
# script for making histograms of trends (expressed in percent)
# must have loaded the *.RData file that has the results of wBTScript
#  it is strongly suggested that if you run this code, you should have done
#  at least 50 bootstrap replicates in wBTScript  (set Mmin to 50 or better still 100)
concChange<-100*bootOut$estC/bootOut$baseConc
cReps <- 100*xConc/bootOut$baseConc
cat("\n The minimum and maximum percent changes are",min(cReps),max(cReps))
cat("\n for your histogram specify the min, max and step size\n it is a good idea to have the min and max straddle zero\n")
message("Enter minimum")
xmin <- as.numeric(readline())
message("Enter maximum")
xmax <- as.numeric(readline())
message("Enter step size")
xstep <- as.numeric(readline())
breaks <- seq(xmin,xmax,xstep)
title<-paste("Histogram of trend in",INFO$paramShortName,"\nFlow Normalized Concentration:",caseSetUp$year1," to ",caseSetUp$year2,"\n",INFO$shortName)
hist(cReps,breaks=breaks,xaxs="i",main=title,freq=FALSE,xlab="Trend in Percent",col="grey")
abline(v=concChange,lwd=3)
abline(v=0,lwd=2,lty=2)
message("now is the time to save your graph, when you are ready to proceed enter any character")
junk <- readline()
rm(junk)
fluxChange<-100*bootOut$estF/bootOut$baseFlux
fReps <- 100*xFlux/bootOut$baseFlux
cat("\n The minimum and maximum percent changes are",min(fReps),max(fReps))
cat("\n for your histogram specify the min, max and step size\n it is a good idea to have the min and max straddle zero\n")
message("Enter minimum")
xmin <- as.numeric(readline())
message("Enter maximum")
xmax <- as.numeric(readline())
message("Enter step size")
xstep <- as.numeric(readline())
breaks <- seq(xmin,xmax,xstep)
title<-paste("Histogram of trend in",INFO$paramShortName,"\nFlow Normalized Flux:",caseSetUp$year1," to ",caseSetUp$year2,"\n",INFO$shortName)
hist(fReps,breaks=breaks,xaxs="i",main=title,freq=FALSE,xlab="Trend in Percent",col="grey")
abline(v=fluxChange,lwd=3)
abline(v=0,lwd=2,lty=2)

  