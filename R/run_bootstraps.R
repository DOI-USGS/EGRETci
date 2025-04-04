#' Run bootstrapping routine
#' 
#' Depending on input arguments, this function will run
#' a complete bootstrapping routine either in series or parallel.
#' 
#' When in parallel, if there are unsuccessful runs, they will be made
#' up in series after the parallel runs are done.#' 
#' 
#' @param eList named list with at least the Daily, Sample, and INFO dataframes
#' @param type_results data frame returned from either \code{\link[EGRET]{runGroups}}
#' or \code{\link[EGRET]{runPairs}} depending on context.
#' @param jitterOn logical, if TRUE, adds "jitter" to the data in an attempt to avoid some numerical problems.
#'   Default = FALSE.  See Details below.
#' @param V numeric a multiplier for addition of jitter to the data, default = 0.2.
#' @param nBoot the maximum number of bootstrap replicates to be used, typically 100
#' @param blockLength integer size of subset, expressed in days.  200 days has been found to be a good choice.
#' @param startSeed sets the random seed value. This is used to make repeatable output.
#' @param type Character can be "pair" or "group".
#' @param run.parallel logical to run bootstrapping in parallel or not
#' 
#' @export
#' @keywords internal
#' 
run_bootstraps <- function(eList, 
                           type_results, 
                           jitterOn, V,
                           blockLength, 
                           startSeed,
                           nBoot,
                           type,
                           run.parallel){
  
  match.arg(type, choices = c("pair", "group"))
  
  xConc <- rep(NA, nBoot)
  xFlux <- rep(NA, nBoot)
  pConc <- rep(NA, nBoot)
  pFlux <- rep(NA, nBoot)
  
  nBootGood <- 0
  
  if(run.parallel){
    `%dopar%` <- foreach::`%dopar%`
    boot_list_out <- foreach::foreach(iBoot = 1:nBoot, 
                                      .packages=c('EGRETci', 'EGRET')) %dopar% {
                                        boot_list <- single_boot_run(iBoot = iBoot, 
                                                                    startSeed = startSeed,
                                                                    eList = eList, 
                                                                    type_results = type_results, 
                                                                    jitterOn = jitterOn, 
                                                                    V = V,
                                                                    blockLength = blockLength,
                                                                    type = type)
                                      }
    
    xConc <- unlist(sapply(boot_list_out, function(x) x[["xConc"]]))
    xFlux <- unlist(sapply(boot_list_out, function(x) x[["xFlux"]]))
    pConc <- unlist(sapply(boot_list_out, function(x) x[["pConc"]]))
    pFlux <- unlist(sapply(boot_list_out, function(x) x[["pFlux"]]))
    
    if(length(xConc) < nBoot){
      # do the last sets NOT in parallel because there's a lot of overhead,
      # potentially for just a few runs:
      iStart <- nBoot + 1
      iEnd <- 2*nBoot - length(xConc)
      message("Running ", iStart-iEnd, " in series")
      for(i in iStart:iEnd){
        boot_list <- single_boot_run(iBoot = iBoot, 
                                    startSeed = startSeed,
                                    eList = eList, 
                                    type_results = type_results, 
                                    jitterOn = jitterOn, 
                                    V = V,
                                    blockLength = blockLength,
                                    type = type)
        
        if(!is.null(boot_list$xConc) & !is.null(boot_list$xFlux)){
          if(nBootGood >= nBoot) {
            break()
          }
          nBootGood <- nBootGood + 1
        }
        xConc <- c(xConc, boot_list$xConc)
        xFlux <- c(xFlux, boot_list$xFlux)
        pConc <- c(pConc, boot_list$pConc)
        pFlux <- c(pFlux, boot_list$pFlux)
      }
      iBoot <- length(pFlux)
    } else {
      xConc <- xConc[1:nBoot]
      xFlux <- xFlux[1:nBoot]
      pConc <- pConc[1:nBoot]
      pFlux <- pFlux[1:nBoot]
      iBoot <- nBoot
    }
    
  } else {
    
    for (iBoot in 1:(2 * nBoot)) {
      boot_list <- single_boot_run(iBoot = iBoot, 
                                  startSeed = startSeed,
                                  eList = eList, 
                                  type_results = type_results, 
                                  jitterOn = jitterOn, 
                                  V = V,
                                  blockLength = blockLength,
                                  type = type)
      
      if(!is.null(boot_list$xConc) & !is.null(boot_list$xFlux)){
        if(nBootGood >= nBoot) {
          break()
        }
        nBootGood <- nBootGood + 1
      }
      xConc <- c(xConc, boot_list$xConc)
      xFlux <- c(xFlux, boot_list$xFlux)
      pConc <- c(pConc, boot_list$pConc)
      pFlux <- c(pFlux, boot_list$pFlux)
      
    }
  }
  
  if (iBoot == 2 * nBoot) {
    message(iBoot, " iterations were run. They only achieved ", 
            nBootGood, " sucessful runs.")
  } else if (iBoot > nBoot) {
    message("It took ", iBoot, " iterations to achieve ", 
            nBoot, " sucessful runs.")
  }

  boot_list_return <- list(xConc = as.numeric(na.omit(xConc)),
                           xFlux = as.numeric(na.omit(xFlux)),
                           pConc = as.numeric(na.omit(pConc)),
                           pFlux = as.numeric(na.omit(pFlux)))

  return(boot_list_return)
  
}