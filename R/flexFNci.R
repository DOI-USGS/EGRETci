#' Flexible flow-normalization confidence intervals
#' 
#' Code for flexFNci
#' 
#' @param eList named list with at least the Daily, Sample, and INFO dataframes. Created from the EGRET package, after running \code{\link[EGRET]{modelEstimation}}.
#' @param nBoot integer maximum number of replicates (called Mmax in paper)
#' @param rs0cy is the calendar year that we want to be the first period for RS. 
#' When we want water years, these are set to the year prior to the desired water year
#' @param rs1cy is the calendar year that we want to be the second period for RS
#' @param repSeed setSeed value
#' @export
#' @examples 
#' library(EGRET)
#' eList <- Choptank_eList
#' rs0cy <- 1985
#' rs1cy <- 2000
#' \dontrun{
#' flexFNlist <- flexFNci(eList, rs0cy, rs1cy)
#' }
flexFNci <- function(eList, rs0cy, rs1cy, nBoot = 100, 
                     repSeed = 38772){
  
  Daily <- getDaily(eList)
  Sample <- getSample(eList)
  INFO <- getInfo(eList)

  surfaceIndexParameters <- surfaceIndex(Daily)
  bottomLogQ <- surfaceIndexParameters[1]
  stepLogQ <- surfaceIndexParameters[2]

  nObs <- length(Sample$ConcHigh)
  
  minNumObs <- if(nObs > 2000) 400 else 100
  minNumUncen <- if(nObs > 2000) 200 else 50
  
  windowY <- INFO$windowY#7
  windowQ <- INFO$windowQ#2
  windowS <- INFO$windowS#0.5
  edgeAdjust <- FALSE
  
  nDaysFull <- length(eList$Daily$MonthSeq)
  fullStart <- eList$Daily$MonthSeq[1]
  fullEnd <- eList$Daily$MonthSeq[nDaysFull]
  fullLength <- (fullEnd - fullStart + 1) / 12
  p0StartMonth <- 10
  p0StartYear <- rs0cy
  p0StartMonthSeq <- ((p0StartYear - 1850) * 12) + p0StartMonth
  paStart <- p0StartMonth
  paLong <- 12
  nPlus1 <- nBoot + 1
  # note the z vectors have units of kg/day
  z00 <- rep(NA,nPlus1)
  z10 <- rep(NA,nPlus1)
  z01 <- rep(NA,nPlus1)
  z11 <- rep(NA,nPlus1)
  blockLength <- 200
  
  set.seed(repSeed)
  # p1 is the first period of the surfaces
  # next we get the period for flow normalization
  windowFlex <- 15
  if(fullLength < windowFlex) stop
  half <- (windowFlex - 1) / 2
  half <- trunc(half)
  q0StartMonthSeqTemp <- p0StartMonthSeq - (half * 12)
  q0StartMonthSeq <- max(fullStart,q0StartMonthSeqTemp)
  q0EndMonthSeq <- q0StartMonthSeq + (windowFlex * 12) - 1
  q0Daily <- subset(eList$Daily,MonthSeq >= q0StartMonthSeq & MonthSeq <= q0EndMonthSeq)
  # now do it again for the later period
  p1StartMonth <- p0StartMonth
  p1StartYear <- rs1cy
  p1StartMonthSeq <- ((p1StartYear - 1850) * 12) + p1StartMonth
  q1StartMonthSeqTempA <- p1StartMonthSeq - (half * 12) # start based window from middle
  q1StartMonthSeqTempB <- fullEnd - (windowFlex * 12) + 1 # start based on end 
  q1StartMonthSeq <- min(q1StartMonthSeqTempA, q1StartMonthSeqTempB)
  q1EndMonthSeq <- q1StartMonthSeq + (windowFlex * 12) - 1
  q1Daily <- subset(eList$Daily,MonthSeq >= q1StartMonthSeq & MonthSeq <= q1EndMonthSeq)
  # done with creating the two flow-normalizing Daily data frames
  # next we estimate the suface for the two calendar years that contain the first
  # period
  # first period RS0
  # we also set up the FD0 data set rawDaily0
  localSample <- eList$Sample
  localINFO <- eList$INFO
  localINFO$paLong <- paLong
  localINFO$paStart <- paStart
  localINFO0 <- localINFO
  localINFO0$nVectorYear <- 33 # this is 33 because it is 2 years of 16 slices, plus 1
  localINFO0$bottomYear <- p0StartYear
  localINFO$stepYear <- 1/16
  p0SurfaceStartDate <- paste(p0StartYear,"-01-01",sep = "")
  p0SurfaceEndDate <- paste(p0StartYear + 2, "-01-01", sep = "")
  rawDaily0 <- subset(eList$Daily, Date >= p0SurfaceStartDate & Date < p0SurfaceEndDate)
  eList0 <- as.egret(localINFO0, rawDaily0, localSample)
  eList0 <- setUpEstimation(eList0, windowY = windowY, windowQ = windowQ, windowS = windowS,
                            minNumObs = minNumObs, minNumUncen = minNumUncen, edgeAdjust =
                              edgeAdjust, verbose = FALSE)
  eList0$INFO$bottomLogQ <- bottomLogQ
  eList0$INFO$stepLogQ <- stepLogQ
  surfaces0 <- estSurfaces(eList0, windowY=7, windowQ=2, windowS=0.5, minNumObs=minNumObs,
                           minNumUncen=minNumUncen, edgeAdjust=edgeAdjust, verbose = FALSE)
  # now the surface for the later period, RS1
  # we also set up the second FD period, FD1
  localINFO1 <- localINFO0
  localINFO1$bottomYear <- p1StartYear
  p1SurfaceStartDate <- paste(p1StartYear,"-01-01",sep = "")
  p1SurfaceEndDate <- paste(p1StartYear + 2, "-01-01", sep = "")
  rawDaily1 <- subset(eList$Daily, Date >= p1SurfaceStartDate & Date < p1SurfaceEndDate)
  eList1 <- as.egret(localINFO1, rawDaily1, localSample)
  eList1 <- setUpEstimation(eList1, windowY = windowY, windowQ = windowQ, windowS = windowS,
                            minNumObs = minNumObs, minNumUncen = minNumUncen, edgeAdjust =
                              edgeAdjust, verbose = FALSE)
  eList1$INFO$bottomLogQ <- bottomLogQ
  eList1$INFO$stepLogQ <- stepLogQ
  surfaces1 <- estSurfaces(eList1, windowY=7, windowQ=2, windowS=0.5, minNumObs=minNumObs,
                                  minNumUncen=minNumUncen, edgeAdjust=edgeAdjust, verbose = FALSE)
  # now do the estimation for RS0 and FD0
  eList0 <- as.egret(eList0$INFO, q0Daily, localSample, surfaces0)
  returnDaily0 <- estDailyFromSurfaces(eList0)
  bootAnnRes0 <- setupYears(localDaily = returnDaily0, paStart = paStart, paLong = paLong)
  bootAnnRes0 <- na.omit(bootAnnRes0)
  z00[nPlus1] <- bootAnnRes0$FNFlux[1]
  #
  # now do results for RS1 and FD1
  eList1 <- as.egret(eList1$INFO, q1Daily, localSample, surfaces1)
  returnDaily1 <- estDailyFromSurfaces(eList1)
  bootAnnRes1 <- setupYears(localDaily = returnDaily1, paStart = paStart, paLong = paLong)
  bootAnnRes1 <- na.omit(bootAnnRes1)
  z11[nPlus1] <- bootAnnRes1$FNFlux[1]
  # now do results for RS0 and FD1
  eList01 <- as.egret(eList1$INFO, q1Daily, localSample, surfaces0)
  returnDaily01 <- estDailyFromSurfaces(eList01)
  bootAnnRes01 <- setupYears(localDaily = returnDaily01, paStart = paStart, paLong = paLong)
  bootAnnRes01 <- na.omit(bootAnnRes01)
  z01[nPlus1] <- bootAnnRes01$FNFlux[1]
  #  finally do the results for RS1 and FD0
  eList10 <- as.egret(eList0$INFO, q0Daily, localSample, surfaces1)
  returnDaily10 <- estDailyFromSurfaces(eList10)
  bootAnnRes10 <- setupYears(localDaily = returnDaily10, paStart = paStart, paLong = paLong)
  bootAnnRes10 <- na.omit(bootAnnRes10)
  z10[nPlus1] <- bootAnnRes10$FNFlux[1]
  #
  # summary of the changes before bootstrap
  zs00 <- z00[nPlus1]
  zs01 <- z01[nPlus1]
  zs10 <- z10[nPlus1]
  zs11 <- z11[nPlus1]
  area <- eList$INFO$drainSqKm
  yieldDelta <- (zs11 - zs00) * 365.25 / area
  RSpart <- 0.5 * (zs10 - zs00 + zs11 - zs01) * 365.25 / area
  FDpart <- 0.5 * (zs01 - zs00 + zs11 - zs10) * 365.25 / area
  pctChange <- 100.0 * (zs11 - zs00) / zs00
  cat("\n change in yield", yieldDelta,"kg / km^2 / year", 
      "\n yield change due to change in response surface",RSpart,
      "\n yield change due to change in flow distribution", FDpart,
      "\n percentage change in yield", pctChange)
  #
  # now the bootstrap replicates
  #bootReps <- as.data.frame(matrix(ncol = 4, nrow = nBoot))
  #colnames(bootReps) <- c("Conc1","Conc2","Flux1","Flux2")
  for(iBoot in 1: nBoot) {
    bootSample <- blockSample(localSample = localSample, blockLength = blockLength)
    # first period
    eList0 <- as.egret(localINFO0, rawDaily0, bootSample)
    eList0 <- setUpEstimation(eList0, windowY = windowY, windowQ = windowQ, windowS = windowS,
                              minNumObs = minNumObs, minNumUncen = minNumUncen, edgeAdjust =
                                edgeAdjust, verbose = FALSE)
    eList0$INFO$bottomLogQ <- bottomLogQ
    eList0$INFO$stepLogQ <- stepLogQ
    surfaces0 <- estSurfaces(eList0, windowY=7, windowQ=2, windowS=0.5, minNumObs=minNumObs,
                                    minNumUncen=minNumUncen, edgeAdjust=edgeAdjust, verbose = FALSE)
    # now the surface for the later period, RS1
    # we also set up the second FD period, FD1
    eList1 <- as.egret(localINFO1, rawDaily1, bootSample)
    eList1 <- setUpEstimation(eList1, windowY = windowY, windowQ = windowQ, windowS = windowS,
                              minNumObs = minNumObs, minNumUncen = minNumUncen, edgeAdjust =
                                edgeAdjust, verbose = FALSE)
    eList1$INFO$bottomLogQ <- bottomLogQ
    eList1$INFO$stepLogQ <- stepLogQ
    surfaces1 <- estSurfaces(eList1, windowY=7, windowQ=2, windowS=0.5, minNumObs=minNumObs,
                                    minNumUncen=minNumUncen, edgeAdjust=edgeAdjust, verbose = FALSE)
    # now do the estimation for RS0 and FD0
    eList0 <- as.egret(eList0$INFO, q0Daily, bootSample, surfaces0)
    returnDaily0 <- estDailyFromSurfaces(eList0)
    bootAnnRes0 <- setupYears(localDaily = returnDaily0, paStart = paStart, paLong = paLong)
    bootAnnRes0 <- na.omit(bootAnnRes0)
    z00[iBoot] <- bootAnnRes0$FNFlux[1]
    #
    # now do results for RS1 and FD1
    eList1 <- as.egret(eList1$INFO, q1Daily, bootSample, surfaces1)
    returnDaily1 <- estDailyFromSurfaces(eList1)
    bootAnnRes1 <- setupYears(localDaily = returnDaily1, paStart = paStart, paLong = paLong)
    bootAnnRes1 <- na.omit(bootAnnRes1)
    z11[iBoot] <- bootAnnRes1$FNFlux[1]
    # now do the results for RS0 and FD1
    eList01 <- as.egret(eList01$INFO, q1Daily, bootSample, surfaces0)
    returnDaily01 <- estDailyFromSurfaces(eList01)
    bootAnnRes01 <- setupYears(localDaily = returnDaily01, paStart = paStart, paLong = paLong)
    bootAnnRes01 <- na.omit(bootAnnRes01)
    z01[iBoot] <- bootAnnRes01$FNFlux[1]
    #  finally do the results for RS1 and FD0
    eList10 <- as.egret(eList10$INFO, q0Daily, bootSample, surfaces1)
    returnDaily10 <- estDailyFromSurfaces(eList10)
    bootAnnRes10 <- setupYears(localDaily = returnDaily10, paStart = paStart, paLong = paLong)
    bootAnnRes10 <- na.omit(bootAnnRes10)
    z10[iBoot] <- bootAnnRes10$FNFlux[1]  
    cat("\n boot rep ",iBoot," done   ")
    cat(z00[iBoot],z11[iBoot],z01[iBoot],z10[iBoot])
  }
  deltaBoth <- rep(0,nBoot)
  deltaRS <- rep(0,nBoot)
  deltaFD <- rep(0,nBoot)
  for(i in 1:nBoot){
    deltaBoth[i] <- 2 * (zs11 - zs00) - z11[i] + z00[i]
    rsPart1 <- 2 * (zs10 - zs00) - z10[i] + z00[i]
    rsPart2 <- 2 * (zs11 - zs01) - z11[i] + z01[i]
    deltaRS[i] <- 0.5 * (rsPart1 + rsPart2)
    fdPart1 <- 2 * (zs01 - zs00) - z01[i] + z00[i]
    fdPart2 <- 2 * (zs11 - zs10) - z11[i] + z10[i]
    deltaFD[i] <- 0.5 * (fdPart1 + fdPart2)
  }
  zBoth <- (sort(deltaBoth)) * 365.25 / area
  cat("\n deltaBoth in kg/km^2/yr\n")
  print(quantile(zBoth,probs = c(0.05, 0.95), type = 6))
  CIBoth <- quantile(zBoth,probs = c(0.05, 0.95), type = 6)
  isPos <- ifelse(zBoth > 0, 1, 0)
  nPos <- sum(isPos)
  cat("\n nPos and likeUp", nPos, nPos / nBoot)
  zRS <- (sort(deltaRS)) * 365.25 / area
  cat("\n\n deltaRS in kg/km^2/yr\n")
  print(quantile(zRS,probs = c(0.05, 0.95), type = 6))
  CIRS <- quantile(zRS,probs = c(0.05, 0.95), type = 6)
  zFD <- (sort(deltaFD)) * 365.25 / area
  cat("\n\n deltaFD in kg/km^2/yr\n")
  print(quantile(zFD,probs = c(0.05, 0.95), type = 6))
  CIFD <- quantile(zFD,probs = c(0.05, 0.95), type = 6)
  flexBoot <- list(dfStart = dfStart, dfEnd = dfEnd, INFO = INFO, rs0cy = rs0cy, rs1cy = rs1cy,
                   windowY = windowY, windowQ = windowQ, windowS = windowS, minNumObs = minNumObs,
                   minNumUncen = minNumUncen, edgeAdjust = edgeAdjust, p0StartMonth = p0StartMonth,
                   paLong = paLong, nBoot = nBoot, blockLength = blockLength, repSeed = repSeed,
                   bottomLogQ = bottomLogQ, stepLogQ = stepLogQ, yieldDelta = yieldDelta,
                   RSpart = RSpart, FDpart = FDpart,  z00 = z00, 
                   z01 = z01, z10 = z10, z11 = z11, zBoth = zBoth,
                   zRS = zRS, zFD = zFD, CIBoth = CIBoth, CIRS = CIRS,
                   CIFD = CIFD, pctChange = pctChange, nPos = nPos,
                   site = site, siteName = siteName, abbrevS = INFO$staAbbrev,
                   abbrevC = INFO$constitAbbrev)
  # fileName <- paste(siteName,".",shortConstit,".RData",sep="")
  # save(flexBoot, file = fileName)
  return(flexBoot)
  
}