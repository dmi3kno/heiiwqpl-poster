tar_target(birdsdata, {
  set.seed(1234)
  
  ## Generate simulated data
  ## data.fn() is defined in bpa-code.txt, available at
  ## http://www.vogelwarte.ch/de/projekte/publikationen/bpa/complete-code-and-data-files-of-the-book.html
  p610 <- read.table("data/p610.txt", header = TRUE)
  y <- p610[,5:9]                         # Grab counts
  y[y > 1] <- 1                           # Convert to det-nondetections
  ever.observed <- apply(y, 1, max)
  wt <- p610$bm[ever.observed == 1]       # Body mass
  yy <- as.matrix(y[ever.observed == 1,]) # Detection histories
  dimnames(yy) <- NULL
  
  ## Augment both data sets
  nz <- 150
  yaug <- rbind(yy, array(0, dim = c(nz, ncol(yy))))
  logwt3 <- c(log(wt^(1/3)), rep(NA, nz))
  
  ## Bundle data
  bsize <- logwt3[1:nrow(yy)]
  list(y = yaug,
       bsize = bsize - mean(bsize),
       M = nrow(yaug),
       `T` = ncol(yaug),
        C = nrow(yy),
        prior_sd_upper = 3)
})
