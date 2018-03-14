context("Two-year trend")
library(EGRET)
eList <- Choptank_eList

test_that("trendSetUp",{

  caseSetUp <- trendSetUp(eList, 
                          year1=1990, 
                          year2=2012, 
                          nBoot = 50, 
                          bootBreak = 39, 
                          blockLength = 200)
  df <- data.frame(year1 = 1990,
                   yearData1 = 1980,
                   year2 = 2012,
                   yearData2 = 2011,
                   numSamples = 606,
                   nBoot = 50,
                   bootBreak = 39,
                   blockLength = 200,
                   confStop = 0.7)
  
  expect_equal(caseSetUp, df)
  expect_error(  caseSetUp <- trendSetUp(eList, 
                                         year1=1970, 
                                         year2=2012, 
                                         nBoot = 50, 
                                         bootBreak = 39, 
                                         blockLength = 200))
  
  expect_error(  caseSetUp <- trendSetUp(eList, 
                                         year1=1980, 
                                         year2=2013, 
                                         nBoot = 50, 
                                         bootBreak = 39, 
                                         blockLength = 200))
  
})

test_that("setForBoot", {
  
  
  INFO <- eList$INFO
  
  eList <- setForBoot(eList, caseSetUp)
  
  INFO2 <- eList$INFO
  
  expect_gt(ncol(INFO2),ncol(INFO))
  
  expect_true(all(c("DecLow","DecHigh") %in% names(INFO2)))
  
})


test_that("wordLike", {
  likeList <- c(0.01, 0.5, 0.55, 0.99)
  Trends <- wordLike(likeList)
  
  expect_equal(Trends, c("Upward trend in concentration is highly unlikely",         
                         "Downward trend in concentration is about as likely as not",
                         "Upward trend in flux is about as likely as not",           
                         "Downward trend in flux is highly likely"))
  
})

test_that("blockSample", {
  Sample <- eList$Sample
  set.seed(1)
  bsReturn <- blockSample(Sample, 25)
  
  expect_equal(bsReturn$ConcLow[1], 0.62)
  expect_equal(bsReturn$Date[1], as.Date("1979-10-24"))
  
})

test_that("pVal", {
  s <- c(0.01, 0.5, 0.55, 0.99)
  pValue <- pVal(s)
  expect_equal(pValue, 0.4)
})

test_that("makeTwoYearsResults", {
  testthat::skip_on_cran()
  twoResultsWaterYear <- makeTwoYearsResults(eList, 1985, 2005)
  expect_equal(floor(twoResultsWaterYear[1:2]), c(1,0))
  
})

test_that("makeCombo", {

  surfaces1 <- c(1,2,3)
  surfaces2 <- c(4, NA, 5)
  surfaces <- makeCombo(surfaces1, surfaces2)
  expect_equal(surfaces, c(5,2,8))
})

test_that("paVector", {
  
  year <- 2000
  paStart <- 10
  paLong <- 12
  vectorYear <- c(seq(1999,2001,0.0833))
  paIndexWaterYear <- paVector(year, paStart, paLong, vectorYear)
  expect_equal(paIndexWaterYear, 10:21)
  
  paStart <- 11
  paLong <- 3
  paIndexWinter <- paVector(year, paStart, paLong, vectorYear)
  expect_equal(paIndexWinter, 11:13)
  
  paStart <- 6
  paLong <- 3
  paIndexSummer <- paVector(year, paStart, paLong, vectorYear)
  expect_equal(paIndexSummer, 18:20)
  
  paStart <- 10
  paLong <- 3
  paIndexLate <- paVector(year, paStart, paLong, vectorYear)
  expect_equal(paIndexLate, 22:24)

  paCalendarYear <- paVector(year, 1, 12, vectorYear)
  expect_equal(paCalendarYear, 14:24)
  
})

test_that("estSliceSurfacesSimpleAlt", {
  testthat::skip_on_cran()
  eList <- Choptank_eList
  
  caseSetUp <- trendSetUp(eList, 
                          year1=1990, 
                          year2=2012, 
                          nBoot = 50, 
                          bootBreak = 39, 
                          blockLength = 200)
  eList <- setForBoot(eList,caseSetUp)
  
  surfaces <- estSliceSurfacesSimpleAlt(eList, 1990)
  expect_equal(surfaces[1:14,1,3], as.numeric(rep(NA, 14)))
  expect_equal(surfaces[1,173,1], 0.16541093)
})

test_that("wBT", {
  
  testthat::skip_on_cran()
  
  eList <- Choptank_eList
  
  caseSetUp <- trendSetUp(eList, 
                          year1=1985, 
                          year2=2005, 
                          nBoot = 5, 
                          bootBreak = 39, 
                          blockLength = 200)
  eList <- setForBoot(eList,caseSetUp)
  
  eBoot <- wBT(eList,caseSetUp)
  
  bootOut <- eBoot$bootOut
  expect_true(bootOut$rejectC)
  
  expect_equal(signif(bootOut$lowC, digits = 6), 0.298427)
  expect_equal(signif(bootOut$likeCUp, digits = 6), 0.916667)
  
  expect_true(bootOut$rejectF)
  expect_equal(eBoot$wordsOut, c("Upward trend in concentration is very likely" , 
                                 "Downward trend in concentration is very unlikely",
                                 "Upward trend in flux is very likely",
                                 "Downward trend in flux is very unlikely"))
  
  
  expect_equal(signif(eBoot$xConc, digits = 2), c(0.31,0.35,0.30,0.34,0.31))
  expect_equal(signif(eBoot$pFlux, digits = 2), c(18,30,31,21,18))
  expect_equal(signif(eBoot$xFlux, digits = 2), c(0.022,0.034,0.034,0.025,0.021))
  expect_equal(signif(eBoot$pConc, digits = 2), c(30,35,30,33,31))
  
})


