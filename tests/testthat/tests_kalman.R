context("Kalman tests")

test_that("Kalman Bootstrapping", {
  eList <- EGRET::Choptank_eList
  dailyBootOut <- genDailyBoot(eList, nBoot = 2, jitterOn = TRUE,
                               nKalman = 2, setSeed = 1)
  
  expect_equal(nrow(dailyBootOut), nrow(eList$Daily))
  expect_equal(ncol(dailyBootOut), 4)

  expect_equal(round(dailyBootOut[1,1], digits = 2), 161.57)
  expect_equal(round(dailyBootOut[2,2], digits = 2), 131.96)
  
})


test_that("Monthly PI", {
  eList <- EGRET::Choptank_eList
  dailyBoot <- Choptank_dailyBootOut
  month_PI <- makeMonthPI(dailyBoot, eList)
  
  expect_true(all(names(month_PI) %in% c("flux", "conc")))  
  df_flux <- month_PI$flux
  df_conc <- month_PI$conc
  expect_equal(nrow(df_flux), length(unique(eList$Daily$MonthSeq)))
  
  expect_true(all(names(df_flux) %in% c("monthSeq", "p1", "p2.5",
                                        "p5", "p10", "p25", "p50",
                                        "p75", "p90", "p95", "p97.5",
                                        "p99")))
  
  expect_true(all(names(df_conc) %in% c("monthSeq", "p1", "p2.5",
                                        "p5", "p10", "p25", "p50",
                                        "p75", "p90", "p95", "p97.5",
                                        "p99")))
  
})

test_that("Daily PI", {
  eList <- EGRET::Choptank_eList
  dailyBoot <- Choptank_dailyBootOut
  daily_PI <- makeDailyPI(dailyBoot, eList)
  
  expect_true(all(names(daily_PI) %in% c("flux", "conc")))  
  df_flux <- daily_PI$flux
  df_conc <- daily_PI$conc
  expect_equal(nrow(df_flux), nrow(eList$Daily))
  
  expect_true(all(names(df_flux) %in% c("Date", "p1", "p2.5",
                                        "p5", "p10", "p25", "p50",
                                        "p75", "p90", "p95", "p97.5",
                                        "p99")))
  
  expect_true(all(names(df_conc) %in% c("Date", "p1", "p2.5",
                                        "p5", "p10", "p25", "p50",
                                        "p75", "p90", "p95", "p97.5",
                                        "p99")))
  
})

test_that("Annual PI", {
  eList <- EGRET::Choptank_eList
  dailyBoot <- Choptank_dailyBootOut
  annual_PI <- makeAnnualPI(dailyBoot, eList)
  
  expect_true(all(names(annual_PI) %in% c("flux", "conc")))  
  df_flux <- annual_PI$flux
  df_conc <- annual_PI$conc
  
  expect_equal(nrow(df_flux), nrow(EGRET::setupYears(eList$Daily)))
  
  expect_true(all(names(df_flux) %in% c("DecYear", "p1", "p2.5",
                                        "p5", "p10", "p25", "p50",
                                        "p75", "p90", "p95", "p97.5",
                                        "p99")))
  
  expect_true(all(names(df_conc) %in% c("DecYear", "p1", "p2.5",
                                        "p5", "p10", "p25", "p50",
                                        "p75", "p90", "p95", "p97.5",
                                        "p99")))
  
})