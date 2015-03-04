library(foreach)
library(doParallel)
library(iterators)


coreOut <- 1
nBoot <- 30
widthCI <- 90
blockLength <- 200

ciLower <- (50-(widthCI/2))/100
ciUpper <- (50+(widthCI/2))/100

probs <- c(ciLower,ciUpper)

nCores <- detectCores() - coreOut
cl <- makeCluster(nCores)
registerDoParallel(cl)
repAnnualResults <- foreach(n = 1:nBoot,.packages=c('EGRETci')) %dopar% {
   annualResults <- bootAnnual(eList, blockLength)  
}
stopCluster(cl)               

CIAnnualResults <- ciBands(eList, repAnnualResults, probs)

plotConcHistBoot(eList, CIAnnualResults)
plotFluxHistBoot(eList, CIAnnualResults)

#
#
