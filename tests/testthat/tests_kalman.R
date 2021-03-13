context("Kalman tests")

test_that("Kalman Bootstrapping", {
  eList <- EGRET::Choptank_eList
  dailyBootOut <- genDailyBoot(eList, nBoot = 2, jitterOn = TRUE,
                               nKalman = 2, setSeed = 1)
  
  expect_equal(nrow(dailyBootOut), nrow(eList$Daily))
  expect_equal(ncol(dailyBootOut), 4)
  
  expect_equal(round(dailyBootOut[1,1], digits = 2), 126.77)
  expect_equal(round(dailyBootOut[2,2], digits = 2), 152.92)
  
})


test_that("Monthly PI", {
  eList <- EGRET::Choptank_eList
  dailyBoot <- Choptank_dailyBootOut
  month_PI <- makeMonthPI(dailyBoot, eList)
  
  
  
  
})