library(testthat)
context("test LCDM and DINA")


test_that("run LCDM model",{
  Qmatrix2 <- cbind(Qmatrix,rep(1,9))
  Qmatrix2[1,1]<-0
  mod2 <- StanLCDM.run(Qmatrix2,respMatrix, iter=100,init.list='cdm')

})


