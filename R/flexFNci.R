#' Flexible flow-normalization confidence intervals
#' 
#' Code for flexFNci
#' 
#' @param eList named list with at least the Daily, Sample, and INFO dataframes. Created from the EGRET package, after running \code{\link[EGRET]{modelEstimation}}.
#' @param nBoot integer maximum number of replicates (called Mmax in paper)
#' @param rs0cy is the calendar year that we want to be the first period for RS - response surface. 
#' When we want water years, these are set to the year prior to the desired water year
#' @param rs1cy is the calendar year that we want to be the second period for RS
#' @param blockLength integer size of subset.
#' @param repSeed setSeed value
#' @param windowFlex interger how big of window
#' @param run.parallel logical to run bootstrapping in parallel or not
#' @param runCI logical to run the confident interval calculations or not
#' @importFrom EGRET getInfo
#' @importFrom EGRET getDaily
#' @importFrom EGRET getSample
#' @importFrom EGRET setUpEstimation
#' @importFrom parallel detectCores
#' @importFrom foreach foreach
#' @importFrom foreach %dopar%
#' @export
#' @examples 
#' library(EGRET)
#' eList <- Choptank_eList
#' rs0cy <- 1985
#' rs1cy <- 2000
#' \dontrun{
#' flexFNlist <- flexFNci(eList, rs0cy, rs1cy, runCI = FALSE)
#' }
flexFNci <- function(eList, rs0cy, rs1cy, windowFlex = 15,
                     run.parallel = TRUE, 
                     nBoot = 100, blockLength = 200,
                     repSeed = 1000, runCI = FALSE){

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
  
  edgeAdjust <- FALSE#there's a question about this...
  
  fullLength <- (diff(range(eList$Daily$MonthSeq)) + 1) / 12
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
  
  # the c vectors need to think about units!
  c00 <- rep(NA,nPlus1)
  c10 <- rep(NA,nPlus1)
  c01 <- rep(NA,nPlus1)
  c11 <- rep(NA,nPlus1)
  
  # p1 is the first period of the surfaces
  # next we get the period for flow normalization

  if(fullLength < windowFlex) stop
  
  half <- (windowFlex - 1) / 2
  half <- trunc(half)
  
  q0StartMonthSeqTemp <- p0StartMonthSeq - (half * 12)
  q0StartMonthSeq <- max(fullStart,q0StartMonthSeqTemp)
  q0EndMonthSeq <- q0StartMonthSeq + (windowFlex * 12) - 1
  q0Daily <- Daily[Daily$MonthSeq >= q0StartMonthSeq & Daily$MonthSeq <= q0EndMonthSeq,]
  # now do it again for the later period
  p1StartMonth <- p0StartMonth
  p1StartYear <- rs1cy
  p1StartMonthSeq <- ((p1StartYear - 1850) * 12) + p1StartMonth
  q1StartMonthSeqTempA <- p1StartMonthSeq - (half * 12) # start based window from middle
  q1StartMonthSeqTempB <- fullEnd - (windowFlex * 12) + 1 # start based on end 
  q1StartMonthSeq <- min(q1StartMonthSeqTempA, q1StartMonthSeqTempB)
  q1EndMonthSeq <- q1StartMonthSeq + (windowFlex * 12) - 1
  q1Daily <- Daily[Daily$MonthSeq >= q1StartMonthSeq & Daily$MonthSeq <= q1EndMonthSeq,]
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
  localINFO$bottomLogQ <- bottomLogQ
  localINFO$stepLogQ <- stepLogQ
  
  p0SurfaceStartDate <- as.Date(paste0(p0StartYear,"-01-01"))
  p0SurfaceEndDate <- as.Date(paste0(p0StartYear + 2, "-01-01"))
  rawDaily0 <- Daily[Daily$Date >= p0SurfaceStartDate & Daily$Date < p0SurfaceEndDate,]
  
  
  # now the surface for the later period, RS1
  # we also set up the second FD period, FD1
  localINFO1 <- localINFO0
  localINFO1$bottomYear <- p1StartYear
  p1SurfaceStartDate <- as.Date(paste0(p1StartYear,"-01-01"))
  p1SurfaceEndDate <- as.Date(paste0(p1StartYear + 2, "-01-01"))
  rawDaily1 <- Daily[Daily$Date >= p1SurfaceStartDate & Daily$Date < p1SurfaceEndDate,]
  
  info_list <- list(localINFO0, localINFO1)
  daily_list <- list(rawDaily0, rawDaily1)
  surface_list <- list()
  
  for(i in 1:2){
    # Note...I did put this in a foreach parallel loop and it actually was slower
    surface_list[[i]] <- setup_and_estSurfaces(info_list[[i]],daily_list[[i]], localSample,
                                                windowY, windowQ, windowS,
                                                minNumObs, minNumUncen, edgeAdjust)
  }

  # Do these in parallel? or at least a loop:
  # now do the estimation for RS0 and FD0
  eList0 <- as.egret(localINFO0, q0Daily, localSample, surface_list[[1]])
  returnDaily0 <- estDailyFromSurfaces(eList0)
  bootAnnRes0 <- setupYears(localDaily = returnDaily0, paStart = paStart, paLong = paLong)
  bootAnnRes0 <- na.omit(bootAnnRes0)

  #
  # now do results for RS1 and FD1
  eList1 <- as.egret(localINFO1, q1Daily, localSample, surface_list[[2]])
  returnDaily1 <- estDailyFromSurfaces(eList1)
  bootAnnRes1 <- setupYears(localDaily = returnDaily1, paStart = paStart, paLong = paLong)
  bootAnnRes1 <- na.omit(bootAnnRes1)

  #
  # now do results for RS0 and FD1
  eList01 <- as.egret(localINFO1, q1Daily, localSample, surface_list[[1]])
  returnDaily01 <- estDailyFromSurfaces(eList01)
  bootAnnRes01 <- setupYears(localDaily = returnDaily01, paStart = paStart, paLong = paLong)
  bootAnnRes01 <- na.omit(bootAnnRes01)

  #
  #  finally do the results for RS1 and FD0
  eList10 <- as.egret(localINFO0, q0Daily, localSample, surface_list[[2]])
  returnDaily10 <- estDailyFromSurfaces(eList10)
  bootAnnRes10 <- setupYears(localDaily = returnDaily10, paStart = paStart, paLong = paLong)
  bootAnnRes10 <- na.omit(bootAnnRes10)

  #
  # summary of the changes before bootstrap
  zs00 <- bootAnnRes0$FNFlux[1]
  zs01 <- bootAnnRes01$FNFlux[1]
  zs10 <- bootAnnRes10$FNFlux[1]
  zs11 <- bootAnnRes1$FNFlux[1]
  
  cs00 <- bootAnnRes0$FNConc[1]
  cs01 <- bootAnnRes01$FNConc[1]
  cs10 <- bootAnnRes10$FNConc[1]
  cs11 <- bootAnnRes1$FNConc[1]
  
  area <- eList$INFO$drainSqKm
  yieldDelta <- (zs11 - zs00) * 365.25 / area
  RSpart <- 0.5 * (zs10 - zs00 + zs11 - zs01) * 365.25 / area
  FDpart <- 0.5 * (zs01 - zs00 + zs11 - zs10) * 365.25 / area
  pctChange <- 100.0 * (zs11 - zs00) / zs00
  cat("\n change in yield", yieldDelta,"kg / km^2 / year", 
      "\n yield change due to change in response surface",RSpart,
      "\n yield change due to change in flow distribution", FDpart,
      "\n percentage change in yield", pctChange)
###################################
  cChange <- cs11 - cs00
  cRSpart <- 0.5 * (cs10 - cs00 + cs11 - cs01)
  cFDpart <- 0.5 * (cs01 - cs00 + cs11 - cs10)
  cpctChange <- 100 * (cs11 - cs00)/cs00
  cat("\n change in conc", cChange, "mg/L", 
      "\n conc change due to change in response surface", 
      cRSpart, "\n conc change due to change in flow distribution", 
      cFDpart, "\n percentage change in conc", cpctChange)
  
  zChange <- (zs11 - zs00) * 365.25
  zRSpart <- 0.5 * (zs10 - zs00 + zs11 - zs01) * 365.25
  zFDpart <- 0.5 * (zs01 - zs00 + zs11 - zs10) * 365.25
  zpctChange <- 100 * (zs11 - zs00)/zs00
  cat("\n change in flux", zChange, "kg/year", 
      "\n flux change due to change in response surface", 
      zRSpart, "\n flux change due to change in flow distribution", 
      zFDpart, "\n percentage change in flux", zpctChange)
  changes <- data.frame(type=c("conc", "flux"),
                        change=c(cChange, zChange), pctChange=c(cpctChange, zpctChange),
                        RSpart=c(cRSpart, zRSpart), FDpart=c(cFDpart, zFDpart))
#############################################
  localINFO0$bottomLogQ <- bottomLogQ
  localINFO0$stepLogQ <- stepLogQ
  localINFO1$bottomLogQ <- bottomLogQ
  localINFO1$stepLogQ <- stepLogQ
  
  if(runCI) {
  
    # now the bootstrap replicates  
    if(run.parallel ){
  
      Z00s <- foreach(n = 1:nBoot, 
                       .packages=c('EGRETci','EGRET')) %dopar% {
                         zList <- getZ00s(n+repSeed, localSample, 
                                          rawDaily0, rawDaily1,
                                          q0Daily, q1Daily,
                                          localINFO0,localINFO1, 
                                          paStart,paLong,
                                          eList01,eList10,
                                          blockLength, 
                                          windowY, windowQ, windowS,
                                          minNumObs, minNumUncen, edgeAdjust)
                       }
  
      z00 <- sapply(Z00s, function(x) x[["z00"]])
      z11 <- sapply(Z00s, function(x) x[["z11"]])
      z01 <- sapply(Z00s, function(x) x[["z01"]])
      z10 <- sapply(Z00s, function(x) x[["z10"]])
      
      c00 <- sapply(Z00s, function(x) x[["c00"]])
      c11 <- sapply(Z00s, function(x) x[["c11"]])
      c01 <- sapply(Z00s, function(x) x[["c01"]])
      c10 <- sapply(Z00s, function(x) x[["c10"]])
      
    } else {
      
      for(iBoot in 1:nBoot) {
        z00s <- getZ00s(iBoot+repSeed,localSample, 
                        rawDaily0, rawDaily1,
                        q0Daily, q1Daily,
                        localINFO0, localINFO1,
                        paStart,paLong,
                        eList01,eList10,
                        blockLength, 
                        windowY, windowQ, windowS,
                        minNumObs, minNumUncen, edgeAdjust) 
        z00[iBoot] <- z00s[["z00"]]
        z11[iBoot] <- z00s[["z11"]]
        z01[iBoot] <- z00s[["z01"]]
        z10[iBoot] <- z00s[["z10"]]
        
        c00[iBoot] <- z00s[["c00"]]
        c11[iBoot] <- z00s[["c11"]]
        c01[iBoot] <- z00s[["c01"]]
        c10[iBoot] <- z00s[["c10"]]
        cat("\n boot rep ",iBoot," done   ")
        cat(z00[iBoot],z11[iBoot],z01[iBoot],z10[iBoot])
      }
    }
    
    z00[nPlus1] <- bootAnnRes0$FNFlux[1]
    z11[nPlus1] <- bootAnnRes1$FNFlux[1]
    z01[nPlus1] <- bootAnnRes01$FNFlux[1]
    z10[nPlus1] <- bootAnnRes10$FNFlux[1]
    c00[nPlus1] <- bootAnnRes0$FNConc[1]
    c11[nPlus1] <- bootAnnRes1$FNConc[1]
    c01[nPlus1] <- bootAnnRes01$FNConc[1]
    c10[nPlus1] <- bootAnnRes10$FNConc[1]
    
    deltaBoth <- 2 * (zs11 - zs00) - z11 + z00
    rsPart1 <- 2 * (zs10 - zs00) - z10 + z00
    rsPart2 <- 2 * (zs11 - zs01) - z11 + z01
    deltaRS <- 0.5 * (rsPart1 + rsPart2)
    fdPart1 <- 2 * (zs01 - zs00) - z01 + z00
    fdPart2 <- 2 * (zs11 - zs10) - z11 + z10
    deltaFD <- 0.5 * (fdPart1 + fdPart2)
    
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
    
    flexBoot <- list(rs0cy = rs0cy, rs1cy = rs1cy,
                     windowY = windowY, windowQ = windowQ, windowS = windowS, minNumObs = minNumObs,
                     minNumUncen = minNumUncen, edgeAdjust = edgeAdjust, p0StartMonth = p0StartMonth,
                     paLong = paLong, nBoot = nBoot, blockLength = blockLength, repSeed = repSeed,
                     bottomLogQ = bottomLogQ, stepLogQ = stepLogQ, yieldDelta = yieldDelta,
                     RSpart = RSpart, FDpart = FDpart,  z00 = z00, 
                     z01 = z01, z10 = z10, z11 = z11, zBoth = zBoth,
                     zRS = zRS, zFD = zFD, CIBoth = CIBoth, CIRS = CIRS,
                     CIFD = CIFD, pctChange = pctChange, nPos = nPos,
                     c01 = c01, c10 = c10, c11 = c11, c00 = c00,
                     site = eList$INFO$site.no, siteName = eList$INFO$station.nm, abbrevS = INFO$staAbbrev,
                     abbrevC = INFO$constitAbbrev)
  
    
  } else {
    flexBoot <- list(rs0cy = rs0cy, rs1cy = rs1cy, windowY = windowY, 
                     windowQ = windowQ, windowS = windowS, minNumObs = minNumObs, 
                     minNumUncen = minNumUncen, edgeAdjust = edgeAdjust, p0StartMonth = p0StartMonth, 
                     paLong = paLong, nBoot = nBoot, blockLength = blockLength, repSeed = repSeed, 
                     bottomLogQ = bottomLogQ, stepLogQ = stepLogQ, cChange = cChange, cRSpart = cRSpart, 
                     cFDpart = cFDpart, cpctChange = cpctChange, zChange = zChange, zRSpart = zRSpart, 
                     zFDpart = zFDpart, zpctChange = zpctChange,site = eList$INFO$site.no, 
                     c11=cs11, c00=cs00, c10=cs10, c01=cs01, z11=zs11, z00=zs00, z10=zs10, z01=zs01, 
                     siteName = eList$INFO$station.nm, abbrevS = INFO$staAbbrev, abbrevC = INFO$constitAbbrev)
      }
  
  return(flexBoot)
}

#' @keywords internal
#' @importFrom EGRET as.egret
#' @importFrom EGRET setUpEstimation
#' @importFrom EGRET estSurfaces
setup_and_estSurfaces <- function(localINFO0, rawDaily0, localSample,
                      windowY, windowQ, windowS,
                      minNumObs, minNumUncen, edgeAdjust){
  eList0 <- as.egret(localINFO0, rawDaily0, localSample)
  eList0 <- setUpEstimation(eList0, windowY = windowY, windowQ = windowQ, windowS = windowS,
                            minNumObs = minNumObs, minNumUncen = minNumUncen, edgeAdjust =
                              edgeAdjust, verbose = FALSE)
  
  surfaces0 <- estSurfaces(eList0, windowY=windowY, windowQ=windowQ, windowS=windowS, minNumObs=minNumObs,
                           minNumUncen=minNumUncen, edgeAdjust=edgeAdjust, verbose = FALSE)
  return(surfaces0)
  
}

#' @keywords internal
#' @importFrom EGRET as.egret
getZ00s <- function(iBoot,localSample, 
                    rawDaily0, rawDaily1,
                    q0Daily, q1Daily,
                    localINFO0, localINFO1,
                    paStart,paLong,
                    eList01,eList10,
                    blockLength, 
                    windowY, windowQ, windowS,
                    minNumObs, minNumUncen, edgeAdjust){
  set.seed(iBoot)
  
  bootSample <- blockSample(localSample = localSample, blockLength = blockLength)
  # first period
  eList0 <- as.egret(localINFO0, rawDaily0, bootSample)
  eList0 <- setUpEstimation(eList0, windowY = windowY, windowQ = windowQ, windowS = windowS,
                            minNumObs = minNumObs, minNumUncen = minNumUncen, edgeAdjust =
                              edgeAdjust, verbose = FALSE)
  
  surfaces0 <- estSurfaces(eList0, windowY=windowY, windowQ=windowQ, windowS=windowS, minNumObs=minNumObs,
                           minNumUncen=minNumUncen, edgeAdjust=edgeAdjust, verbose = FALSE)
  
  # now the surface for the later period, RS1
  # we also set up the second FD period, FD1
  eList1 <- as.egret(localINFO1, rawDaily1, bootSample)
  eList1 <- setUpEstimation(eList1, windowY = windowY, windowQ = windowQ, windowS = windowS,
                            minNumObs = minNumObs, minNumUncen = minNumUncen, edgeAdjust =
                              edgeAdjust, verbose = FALSE)
  
  surfaces1 <- estSurfaces(eList1, windowY=windowY, windowQ=windowQ, windowS=windowS, minNumObs=minNumObs,
                           minNumUncen=minNumUncen, edgeAdjust=edgeAdjust, verbose = FALSE)
  
  # now do the estimation for RS0 and FD0
  eList0 <- as.egret(eList0$INFO, q0Daily, bootSample, surfaces0)
  returnDaily0 <- estDailyFromSurfaces(eList0)
  bootAnnRes0 <- setupYears(localDaily = returnDaily0, paStart = paStart, paLong = paLong)
  bootAnnRes0 <- na.omit(bootAnnRes0)
  z00 <- bootAnnRes0$FNFlux[1]
  c00 <- bootAnnRes0$FNConc[1]
  #
  # now do results for RS1 and FD1
  eList1 <- as.egret(eList1$INFO, q1Daily, bootSample, surfaces1)
  returnDaily1 <- estDailyFromSurfaces(eList1)
  bootAnnRes1 <- setupYears(localDaily = returnDaily1, paStart = paStart, paLong = paLong)
  bootAnnRes1 <- na.omit(bootAnnRes1)
  z11 <- bootAnnRes1$FNFlux[1]
  c11 <- bootAnnRes1$FNConc[1]
  
  # now do the results for RS0 and FD1
  eList01 <- as.egret(eList01$INFO, q1Daily, bootSample, surfaces0)
  returnDaily01 <- estDailyFromSurfaces(eList01)
  bootAnnRes01 <- setupYears(localDaily = returnDaily01, paStart = paStart, paLong = paLong)
  bootAnnRes01 <- na.omit(bootAnnRes01)
  z01 <- bootAnnRes01$FNFlux[1]
  c01 <- bootAnnRes01$FNConc[1]
  
  #  finally do the results for RS1 and FD0
  eList10 <- as.egret(eList10$INFO, q0Daily, bootSample, surfaces1)
  returnDaily10 <- estDailyFromSurfaces(eList10)
  bootAnnRes10 <- setupYears(localDaily = returnDaily10, paStart = paStart, paLong = paLong)
  bootAnnRes10 <- na.omit(bootAnnRes10)
  z10 <- bootAnnRes10$FNFlux[1]  
  c10 <- bootAnnRes10$FNConc[1] 
  
  return(list(z00=z00,z11=z11,z01=z01,z10=z10,
              c00=c00,c11=c11,c01=c01,c10=c10))
}
