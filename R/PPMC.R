# a function to calculate pecentage of extreme o-values for item-pair odds ratios
ppmc <- function(stan_model, respMatrix, n.sim = 6000) {
  library(Rlab)
  library(MCMCpack)

  mod1 <- stan_model
  x <- as.matrix(mod1)
  PImatCols <- grep(names(mod1), pattern = "PImat", value = T)
  VcCols <- grep(names(mod1), pattern = "Vc", value = T)

  n.sim = 6000
  n.iter = mod1@stan_args[[1]]$iter
  n <- nrow(x)
  Vc <- x[(n.iter/2 + 1): nrow(x), VcCols]
  PImat <- x[(n.iter/2 + 1): nrow(x), PImatCols]

  # number of observation
  i <- dim(respMatrix)[1]
  # number of item
  j <- dim(respMatrix)[2]

  # simulate response matrix based on PImat and Vc
  SimRespMatrix = array(data = NA, dim = c(i, j, n.sim))
  for (sim in 1:n.sim) {
    iter <- sample(1:n.iter, 1, replace = TRUE)
    for (item in 1:j) {
      # extract Predicted Probability Matrix
      PImat.iter <- PImat[iter,] %>% matrix(9, 8)
      # draw from dirichlet distribution
      Vc.iter <- rdirichlet(1, rep(1, 8))
      p[item] <- PImat.iter[item,] %*% t(Vc.iter)
      # generating simulated data
      SimRespMatrix[, item, sim] <- rbern(1000, p[item])
    }
  }

  # Item-Pair Odds Ratios Matrix
  oddsratio <- function(binarydata) {
    or <- rep(NA, ncol(binarydata) * (ncol(binarydata) - 1) / 2)
    num <- 0
    for (i in 1:(ncol(binarydata)-1)) {
      for (j in (i+1):ncol(binarydata)) {
        num <- num + 1
        mat <- table(binarydata[,i], binarydata[,j])
        or[num] <- (mat[1,1] * mat[2,2]) / (mat[1,2] * mat[2,1])
      }
    }
    or
  }

  # calculate odd ratios for simulated data.
  t <- sapply(1:n.sim, function(x) oddsratio(SimRespMatrix[,,x]))
  # odds ratio vector for original response matrix
  data_or <- oddsratio(respMatrix)
  perc_or <- rep(NA, length(data_or))
  for (i in 1:length(data_or)) {
    corECDF <- ecdf(t[i,])
    perc_or[i] <- corECDF(data_or[i])
  }

  # percentage of posterior predicted extreme p-values greater than .80 or less than .20
  p_extreme <- sum(perc_or > 0.8 | perc_or < 0.2) / length(data_or)
  p_extreme
}
