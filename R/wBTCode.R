# This code carries out the WRTDS bootstrap test (WBT) (2015-02-26)
#   as described in the draft manuscript on that subject by Hirsch, Archfield, and De Cicco
#   dated February 12, 2015 with one update on the words used to describe likelihoods 
#     Note that this version will only do tests on the Water Year
#     Later on, when fully packaged it will do any Period of Analysis
#
#     the arguments are two ojects (eList and caseSetUp):
#   eList (which is just a standard EGRET eList, it must contain
#     at least the INFO, Daily and Sample data frames,  it doesn't need surfaces)
#   caseSetUp which define the details of how the WBT will be run
#       these are the values in caseSetUp
#         year1 = the water year that is the start of the trend period (an integer)
#         yearData1 = the water year that is the start of the data set (an integer)
#         year2 = the water year that is the end of the trend period (an integer)
#         yearData2 = the water year that is the end of the data set (an integer)
#         numSamples = number of samples in eList$Sample 
#         nBoot = maximum number of replicates that will be done (called Mmax in paper)
#         bootBreak = the minimum number of replicates to be done (called Mmin in paper)
#         blockLength = length of blocks for bootstrap (called B in the paper)
#         confStop = 1 - alphap, the width of the confidence interval used in adaptive
#              stopping rule (paper suggests alphap of 0.3 so confStop would be 0.7)
#  
#   What is returned is a list (called eBoot) of four objects, they are
#   bootOut (which is a whole bunch of results from the finished bootstrap estimation)
#       the names of these variables is prety self-explanatory except perhaps for the last three
#       baseConc = the FN Concentration in the first water year of the trend period (mg/L) 
#           these are needed to do percentage changes in FN Concentration
#       baseFlux = the FN Flux in the first water year of the trend period (10^6 kg/yr)
#           these are needed to do percentage changes in FN Flux
#       iBoot = the number of bootstrap replicates actually run
#   wordsOut (which is a vector of 4 phrases describing the likelihood of up or down trends) 
#	xConc  (the bootstrap replicates for change in concentration) a vector of length iBoot
#       the units are mg/L
#   xFlux  (the bootstrap replicates for change in flux) a vector of length iBoot
#       the units are 10^6 kg/yr







wBT<-function(eList,caseSetUp, ){
	localINFO <- eList$INFO
	localDaily <- eList$Daily
	localSample <- eList$Sample
	prob<-c(0.025, 0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95, 0.975)
	bootOut <- as.data.frame(matrix(ncol=25,nrow=1))
  colnames(bootOut) <- c("rejectC","pValC","estC","lowC","upC",
                         "lowC50","upC50","lowC95","upC95","likeCUp",
                         "likeCDown","rejectF","pValF","estF","lowF",
                         "upF","lowF50","upF50","lowF95","upF95","likeFUp",
                         "likeFDown","baseConc","baseFlux","iBoot")
	year1 <- caseSetUp$year1
	yearData1 <- caseSetUp$yearData1
	year2 <- caseSetUp$year2
	yearData2 <- caseSetUp$yearData2
	numSamples <- caseSetUp$numSamples
	nBoot <- caseSetUp$nBoot
	bootBreak <- caseSetUp$bootBreak
	blockLength <-caseSetUp$blockLength
	confStop <- caseSetUp$confStop
  xConc<-rep(NA,nBoot)
  xFlux<-rep(NA,nBoot)
  posXConc<-0
  posXFlux<-0
  surfaces1 <- estSliceSurfacesSimpleAlt(eList,year1)
  surfaces2 <- estSliceSurfacesSimpleAlt(eList,year2)
  combo <- makeCombo(surfaces1,surfaces2)
  eListCombo <- as.egret(localINFO,localDaily,localSample,combo)
  res <- makeTwoYearsResults(eListCombo,year1,year2)
  regDeltaConc <- res[2] - res[1]
  estC <- regDeltaConc
  baseConc <- res[1]
  regDeltaConcPct <- (regDeltaConc/baseConc) * 100
  regDeltaFlux <- (res[4] - res[3]) * 0.00036525
  estF <- regDeltaFlux
  baseFlux <- res[3] * 0.00036525
  regDeltaFluxPct <- (regDeltaFlux/baseFlux) * 100
  fcc <- format(regDeltaConc, digits = 3, width = 7)
  ffc <- format(regDeltaFlux, digits = 4, width = 8)
	cat("\n\n WRTDS estimated concentration change is",fcc," mg/L" )
	cat("\n WRTDS estimated flux change is        ",ffc," 10^6 kg/yr")
	cat("\n value is bootstrap replicate result (deltack or deltafk in paper)")
	cat("\n nPos is cumulative number of positive trends")
	cat("\n post_p is posterior mean estimate of probability of a positive trend")
	cat("\n Lower and Upper are estimates of the 90% CI values for magnitude of trend")
	cat("\n\n      rep              Concentration             |              Flux")
  cat("\n          value     nPos post_p   Lower   Upper  |     value   nPos  post_p    Lower   Upper")
	for (iBoot in 1:nBoot){
		bootSample<-blockSample(localSample=localSample,blockLength=blockLength)
		eListBoot <- as.egret(localINFO,localDaily,bootSample,NA)
 		surfaces1 <- estSliceSurfacesSimpleAlt(eListBoot,year1)
		surfaces2 <- estSliceSurfacesSimpleAlt(eListBoot,year2)
		combo <- makeCombo(surfaces1,surfaces2)
		eListBoot <- as.egret(localINFO,localDaily,bootSample,combo)
		res <- makeTwoYearsResults(eListBoot,year1,year2)
#			cat("\nres",res[1],res[2],res[3],res[4])
		xConc[iBoot] <- (2 * regDeltaConc) - (res[2] - res[1])
		xFlux[iBoot] <- (2 * regDeltaFlux) - ((res[4] - res[3])*0.00036525)
#			cat("\nxConc[iBoot],iBoot",xConc[iBoot],iBoot)
		posXConc <- ifelse(xConc[iBoot]>0,posXConc+1,posXConc)
		binomIntConc<-binom.bayes(posXConc,iBoot,confStop,"central")
		belowConc <- ifelse(binomIntConc$upper<0.05,1,0)
		aboveConc <- ifelse(binomIntConc$lower>0.95,1,0)
		midConc <- ifelse(binomIntConc$lower>0.05 & binomIntConc$upper<0.95,1,0)
		posXFlux <- ifelse(xFlux[iBoot]>0,posXFlux+1,posXFlux)
		binomIntFlux<-binom.bayes(posXFlux,iBoot,confStop,"central")
		belowFlux <- ifelse(binomIntFlux$upper<0.05,1,0)
		aboveFlux <- ifelse(binomIntFlux$lower>0.95,1,0)
		midFlux <- ifelse(binomIntFlux$lower>0.05 & binomIntFlux$upper<0.95,1,0)
#cat("\niBoot",iBoot,belowConc,aboveConc,midConc,belowFlux,aboveFlux,midFlux)
#	cat("\niBoot,xConc[1],xConc[iBoot]",iBoot,xConc[1],xConc[iBoot])			
		quantConc<-quantile(xConc[1:iBoot],prob,type=6)
		lowConc <- quantConc[2]
		highConc <- quantConc[8]
		quantFlux<-quantile(xFlux[1:iBoot],prob,type=6)
		lowFlux <- quantFlux[2]
		highFlux <- quantFlux[8]
		prints<-c(format(iBoot,digits=3,width=7),format(xConc[iBoot],digits=3,width=7),format(posXConc,digits=3,width=5),format(binomIntConc$mean,digits=3,width=7),format(quantConc[2],digits=3,width=7),format(quantConc[8],digits=3,width=7),"  |  ",format(xFlux[iBoot],digits=4,width=8),format(posXFlux,digits=3,width=5),format(binomIntFlux$mean,digits=3,width=7),format(quantFlux[2],digits=4,width=8),format(quantFlux[8],digits=4,width=8))
		cat("\n",prints)
		test1 <- if(belowConc + aboveConc + midConc > 0.5 & belowFlux + aboveFlux + midFlux > 0.5 & iBoot >= bootBreak & iBoot > 30) 1 else 0
		test2 <- if(midConc > 0.5 & midFlux > 0.5 & iBoot >= bootBreak & iBoot <= 30) 1 else 0
		if(test1 + test2 > 0.5) break
	}

	rejectC <- lowConc * highConc > 0
	rejectF <- lowFlux * highFlux > 0
			
	cat("\n\nShould we reject Ho that Flow Normalized Concentration Trend = 0 ?", words(rejectC))
	fquantConc <- format(quantConc,digits=3,width=8)
	cat("\n best estimate is",fcc, "mg/L\n  Lower and Upper 90% CIs",fquantConc[2],fquantConc[8])
	lowC <- quantConc[2]
	upC <- quantConc[8]
	cat("\n also 95% CIs",fquantConc[1],fquantConc[9],"\n and 50% CIs",fquantConc[4],fquantConc[6])
	lowC50 <- quantConc[4]
	upC50 <- quantConc[6]
	lowC95 <- quantConc[1]
	upC95 <- quantConc[9]
	p <- binomIntConc$mean
	pValC <- 2*(min(p,(1-p)))
	cat("\n approximate two-sided p-value for Conc",format(pValC,digits=2,width=9))
	if(posXConc==0|posXConc==iBoot) cat("\n* Note p-value should be considered to be < stated value")
	likeCUp <- (posXConc + 0.5) / (iBoot + 1)
	likeCDown <- 1 - likeCUp
	cat("\n Likelihood that Flow Normalized Concentration is trending up =",format(likeCUp,digits=3,width=10)," is trending down =",format(likeCDown,digits=3,width=10))
	cat("\n\nShould we reject Ho that Flow Normalized Flux Trend = 0 ?",words(rejectF))
	fquantFlux <- format(quantFlux,digits=3,width=8)
	cat("\n best estimate is",ffc, "10^6 kg/year\n  Lower and Upper 90% CIs",fquantFlux[2],fquantFlux[8])
	lowF <- quantFlux[2]
	upF <- quantFlux[8]
	cat("\n also 95% CIs",fquantFlux[1],fquantFlux[9],"\n and 50% CIs",fquantFlux[4],fquantFlux[6])
	lowF50 <- quantFlux[4]
	upF50 <- quantFlux[6]
	lowF95 <- quantFlux[1]
	upF95 <- quantFlux[9]
	p <- binomIntFlux$mean
	pValF <- 2*(min(p,(1-p)))
	cat("\n approximate two-sided p-value for Flux",format(pValF,digits=2,width=9))
	if(posXFlux==0|posXFlux==iBoot) cat("\n* Note p-value should be considered to be < stated value")
	
	likeFUp <- (posXFlux + 0.5) / (iBoot + 1)
	likeFDown <- 1 - likeFUp

  cat("\n Likelihood that Flow Normalized Flux is trending up =",format(likeFUp,digits=3,width=10)," is trending down=",format(likeFDown,digits=3,width=10))			
  bootOut<-data.frame(rejectC,pValC,estC,lowC,upC,lowC50,upC50,lowC95,upC95,likeCUp,likeCDown,rejectF,pValF,estF,lowF,upF,lowF50,upF50,lowF95,upF95,likeFUp,likeFDown,baseConc,baseFlux,iBoot)
  likeList <- c(likeCUp,likeCDown,likeFUp,likeFDown)
  wordsOut <- wordLike(likeList)
  cat("\n\n",format(wordsOut[1],width=30),"\n",format(wordsOut[3],width=30))
  cat("\n",format(wordsOut[2],width=30),"\n",format(wordsOut[4],width=30))
  xConc <- xConc[1:iBoot]
  xFlux <- xFlux[1:iBoot]
  eBoot <- list(bootOut=bootOut,wordsOut=wordsOut,xConc=xConc,xFlux=xFlux)
  return(eBoot)						
}
#
#
#
#
estSliceSurfacesSimpleAlt<-function(eList,year){
  localINFO <- eList$INFO
  localSample <- eList$Sample
  localDaily <- eList$Daily
  windowY <- localINFO$windowY
  windowS <- localINFO$windowS
  windowQ <- localINFO$windowQ
  edgeAdjust <- localINFO$edgeAdjust
  originalColumns <- names(localSample)
  minNumObs<-localINFO$minNumObs
  minNumUncen<-localINFO$minNumUncen
  bottomLogQ<-localINFO$bottomLogQ
  stepLogQ <- localINFO$stepLogQ
  topLogQ<-bottomLogQ + 13 * stepLogQ
  vectorLogQ<-seq(bottomLogQ,topLogQ,stepLogQ)
  nVectorLogQ <- localINFO$nVectorLogQ
  stepYear<- localINFO$stepYear
  bottomYear<-localINFO$bottomYear
  nVectorYear <- localINFO$nVectorYear
  topYear<-bottomYear + (nVectorYear - 1)* stepYear 
  vectorYear <-seq(bottomYear,topYear,stepYear)
  surfaces<-array(NA,dim=c(14,length(vectorYear),3))
  if(!is.null(localINFO$paStart) & !is.null(localINFO$paLong)){
    vectorIndex <- paVector(year,localINFO$paStart,localINFO$paLong,vectorYear)
    vectorIndex <- c(vectorIndex[1]-1,vectorIndex,vectorIndex[length(vectorIndex)]+1)
  } else {
    # Water year
    vectorIndex <- paVector(year,10,12,vectorYear)
    vectorIndex <- c(vectorIndex[1]-1,vectorIndex,vectorIndex[length(vectorIndex)]+1)
  }
  
  vectorYear <- vectorYear[vectorIndex]
  nVectorYear<-length(vectorYear)
  estPtLogQ<-rep(vectorLogQ,nVectorYear)
  estPtYear<-rep(vectorYear,each=14)
  
#  colToKeep <- c("ConcLow","ConcHigh","Uncen","DecYear","SinDY","CosDY","LogQ")
#  cat("\ncolToKeep",colToKeep)
#  localSampleMin <- localSample[,which(originalColumns %in% colToKeep)]
#  cat("\ngot past localSampleMin")
  numDays <- localINFO$numDays
  DecLow <- localINFO$DecLow
  DecHigh <- localINFO$DecHigh
  resultSurvReg<-runSurvReg(estPtYear,estPtLogQ,numDays,DecLow,DecHigh, localSample,windowY,windowQ,windowS,minNumObs,minNumUncen,interactive=FALSE,edgeAdjust)
  
  for(iQ in 1:14) {
    for(iY in 1:length(vectorIndex)){ 
      k<-(iY-1)*14+iQ
      surfaces[iQ,vectorIndex[iY],]<-resultSurvReg[k,]
    }
  }
  
  return(surfaces)
}
#
#
#
#
paVector <- function(year,paStart,paLong, vectorYear){
  
  nDaysInMonth <- c(31,28,31,30,31,30,31,31,30,31,30,31)
  
  Jan1Lo <- as.POSIXct(paste(year, "-01-01 00:00:00",sep=""),format="%Y-%m-%d %H:%M:%S")
  Jan1Hi <- as.POSIXct(paste(year+1, "-01-01 00:00:00",sep=""),format="%Y-%m-%d %H:%M:%S")
  
  if (paStart + paLong > 12){
    # Crosses January      
    Lo <- as.POSIXct(paste(year-1,"-", paStart,"-01 00:00:00",sep=""),format="%Y-%m-%d %H:%M:%S")
    Hi <- as.POSIXct(paste(year,"-", (paStart+paLong-13),"-",nDaysInMonth[paStart+paLong-13], " 23:59:59",sep=""),format="%Y-%m-%d %H:%M:%S")
    
    Jan1HiMinus <- as.POSIXct(paste(year-1, "-01-01 00:00:00",sep=""),format="%Y-%m-%d %H:%M:%S")
    
    yearsToAdd <- c(year-1,year)
    
    low <- c(Jan1HiMinus, Jan1Lo)
    hi <- c(Jan1Lo, Jan1Hi)
    
  } else {
    # Completely comprises a calendar year
    
    Lo <- as.POSIXct(paste(year,"-", paStart,"-01 00:00:00",sep=""),format="%Y-%m-%d %H:%M:%S")
    Hi <- as.POSIXct(paste(year,"-", paStart+paLong-1,"-",nDaysInMonth[paStart+paLong-1], " 23:59:59",sep=""),format="%Y-%m-%d %H:%M:%S")
    
    low <- c(Jan1Lo, Jan1Lo)
    hi <- c(Jan1Hi, Jan1Hi)
    
    yearsToAdd <- c(year,year)
    
  }
  
  dates <- c(Lo, Hi)
  
  decimal <- as.numeric(difftime(dates, low, units = "secs"))
  nonzero <- decimal != 0
  decimal[nonzero] <- decimal[nonzero]/as.numeric(difftime(hi, low, units = "secs"))
  decimal <- yearsToAdd + decimal
  
  vectorIndex <- which(vectorYear >= decimal[1] & vectorYear <= decimal[2])
  
  return(vectorIndex)
}
#
#

#
makeCombo <- function (surfaces1,surfaces2) {
	surfaces1[is.na(surfaces1)]<-0
	surfaces2[is.na(surfaces2)]<-0
	combo <- surfaces1 + surfaces2
	combo[combo == 0] <- NA
	return(combo)
}
#
#
makeTwoYearsResults <- function(eList,year1,year2){
# note this thing only works if the pa is water year
	returnDaily <- estDailyFromSurfaces(eList)
	bootAnnRes<-setupYears(localDaily=returnDaily)
	AnnBase <- bootAnnRes[1,1]
	index1 <- year1 - trunc(AnnBase) + 1
	index2 <- year2 - trunc(AnnBase) + 1
	twoYearsResults <- c(bootAnnRes$FNConc[index1],bootAnnRes$FNConc[index2],bootAnnRes$FNFlux[index1],bootAnnRes$FNFlux[index2])
	return(twoYearsResults)
}
#
#
#
#
setForBoot<-function (eList,windowY,windowQ,windowS,edgeAdjust) {
#  does the setup functions usually done by modelEstimation
	localINFO <- eList$INFO
	localDaily <- eList$Daily
localSample <- eList$Sample
  numDays <- length(localDaily$DecYear)
  DecLow <- localDaily$DecYear[1]
  DecHigh <- localDaily$DecYear[numDays]
  numSamples <- length(localSample$Julian)
  surfaceIndexParameters <- surfaceIndex(localDaily)
  localINFO$bottomLogQ <- surfaceIndexParameters[1]
  localINFO$stepLogQ <- surfaceIndexParameters[2]
  localINFO$nVectorLogQ <- surfaceIndexParameters[3]
  localINFO$bottomYear <- surfaceIndexParameters[4]
  localINFO$stepYear <- surfaceIndexParameters[5]
  localINFO$nVectorYear <- surfaceIndexParameters[6]
  localINFO$windowY <- windowY
  localINFO$windowQ <- 2
  localINFO$windowS <- 0.5
  localINFO$minNumObs <- min(100,numSamples-20)
  localINFO$minNumUncen <- min(50,numSamples-20)
  localINFO$numDays <- numDays
  localINFO$DecLow <- DecLow
  localINFO$DecHigh <- DecHigh
  localINFO$edgeAdjust <- edgeAdjust
  eList$INFO <- localINFO
  return(eList)
}
#
#
#
#
blockSample <- function(localSample, blockLength){
  numSamples <- length(localSample$Julian)
  dayOne <- localSample$Julian[1]
  newSample <- data.frame()
  firstJulian <- localSample$Julian[1] - blockLength + 1
  lastJulian <- localSample$Julian[numSamples]
  possibleStarts <- seq(firstJulian,lastJulian)
  while(nrow(newSample) <= nrow(localSample)){
    randomDate <- sample(possibleStarts, 1)
    blockStart <- max(randomDate,dayOne)
    blockEnd <- min(lastJulian,randomDate+blockLength-1)
    oneYear <- subset(localSample, localSample$Julian >= blockStart & 
                        localSample$Julian < blockEnd)
    newSample <- rbind(oneYear, newSample)
        
  }
  newSample <- newSample[-c((nrow(localSample)+1):nrow(newSample)),]
  newSample <- newSample[order(newSample$Julian),]
  return(newSample)
}
#
#
wordLike <- function(likeList){
	firstPart <- c("Upward trend in concentration is","Downward trend in concentration is","Upward trend in flux is","Downward trend in flux is")
	secondPart <- c("highly unlikely","very unlikely","unlikely","about as likely as not","likely","very likely","highly likely")
	breaks<-c(0,0.05,0.1,0.33,0.66,0.9,0.95,1)
	levelLike <- cut(likeList,breaks=breaks,labels=FALSE)
	wordLikeFour <- paste(firstPart,secondPart[levelLike],sep=" ")
	return(wordLikeFour)
}
#
#
words <- function(z){
	out <- if(z) "Reject Ho" else "Do Not Reject Ho"
	return(out)	
}
