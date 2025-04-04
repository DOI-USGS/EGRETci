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
    boot_list_out <- foreach::foreach(iBoot = 1:ceiling(1.25*nBoot), 
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
      # do the last set NOT in parallel because there's a lot of overhead,
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