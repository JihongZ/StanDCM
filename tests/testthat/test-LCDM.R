library(testthat)
context("test LCDM and DINA")

#test_check("StanDCM")


#load("../../data/data.RData")

# test_that("runing LCDM model", {
#   StanLCDM.run(Qmatrix,respMatrix,iter=20,init.list='cdm',chains = 3)
# })

# test_that("test Loo fit model comparision", {
#   mod1<-StanLCDM.run(Qmatrix,respMatrix,iter=100,init.list='cdm',chains = 3)
#   StanDCM::StanLCDM.loofit(estimatedMod)
# })

# test_that("test DINO script", {
#   StanDINO.script(Qmatrix)
# })


test_that("run LCDM model",{
  Qmatrix2 <- cbind(Qmatrix,rep(1,9))
  Qmatrix2[1,1]<-0
  mod2 <- StanLCDM.run(Qmatrix2,respMatrix, iter=100,init.list='cdm')

})

test_that("run NCRUM model",{
  Qmatrix2 <- cbind(Qmatrix,rep(1,9))
  Qmatrix2[1,1]<-0
  StanNCRUM.run(Qmatrix2,respMatrix, iter=100,init.list='cdm')
})

test_that("run CRUM model",{
  Qmatrix2 <- cbind(Qmatrix,rep(1,9))
  Qmatrix2[1,1]<-0
  StanCRUM.run(Qmatrix2,respMatrix, iter=100,init.list='cdm')
})

test_that("run ORDM model",{
  Qmatrix2 <- cbind(Qmatrix,rep(1,9))
  Qmatrix2[1,1]<-0
  StanORDM.run(Qmatrix2,respMatrix, iter=100,init.list='cdm')
})



test_that("ddd",{
  # Item-Pair Odds Ratios Matrix
  # oddsratio <- function(binarydata) {
  #   or <- rep(NA, ncol(binarydata) * (ncol(binarydata) - 1) / 2)
  #   num <- 0
  #   for (i in 1:(ncol(binarydata)-1)) {
  #     for (j in (i+1):ncol(binarydata)) {
  #       num <- num + 1
  #       mat <- table(binarydata[,i], binarydata[,j])
  #       or[num] <- (mat[1,1] * mat[2,2]) / (mat[1,2] * mat[2,1])
  #     }
  #   }
  #   or
  # }

  # calculate odd ratios for simulated data.
  # t <- sapply(1:n.sim, function(x) oddsratio(SimRespMatrix[,,x]))
  # odds ratio vector for original response matrix
  # data_or <- oddsratio(respMatrix)
  # perc_or <- rep(NA, length(data_or))
  # for (i in 1:length(data_or)) {
  #   corECDF <- ecdf(t[i,])
  #   perc_or[i] <- corECDF(data_or[i])
  # }


  # percentage of posterior predicted extreme p-values greater than .80 or less than .20
  #p_extreme <- sum(perc_or > 0.8 | perc_or < 0.2) / length(data_or)
  #p_extreme

})


test_that("test StanPrior.update",{
  Example
  priorUpdate.matrix<-matrix(0,2,2)
  priorUpdate.matrix[1,1]<-"l1_0~normal(0,15)"
  priorUpdate.matrix[2,1]<-"l4_0~normal(0,15)"
  priorUpdate.matrix[1,2]<-"l1_0~normal(2,10)"
  priorUpdate.matrix[2,2]<-"l4_0~normal(2,10)"
  StanPrior.update(priorUpdate.matrix,"./CRUM_uninf.stan")
})

