library(testthat)
context("test LCDM and DINA")


test_that("run LCDM model",{
  Qmatrix2 <- cbind(Qmatrix,rep(1,9))
  Qmatrix2[1,1]<-0
  mod2 <- StanLCDM.run(Qmatrix, respMatrix, iter= 1000, init.list='cdm')
  mod.gdina <- GDINA:::GDINA(dat = respMatrix, Q = Qmatrix, model = "logitGDINA")

  summary(mod2)
  total <- extract(mod2)
  colnames(total)
  N <- nrow(respMatrix)
  l1_11 <- rstan::extract(mod2, "l1_11")

  ## Intercepts
  l1_0 <- rstan::extract(mod2, "l1_0")$l1_0[1000]
  l1_11 <- rstan::extract(mod2, "l1_11")$l1_11[1000]
  l2_0 <- rstan::extract(mod2, "l2_0")$l2_0[1000]
  l2_11 <- rstan::extract(mod2, "l2_11")$l2_11[1000]
  l3_0 <- rstan::extract(mod2, "l3_0")$l3_0[1000]
  l3_12 <- rstan::extract(mod2, "l3_12")$l3_12[1000]
  l9_0 <- rstan::extract(mod2, "l9_0")$l9_0[1000]
  l9_11 <- rstan::extract(mod2, "l9_11")$l9_11[1000]
  l9_13 <- rstan::extract(mod2, "l9_13")$l9_13[1000]
  l9_213 <- rstan::extract(mod2, "l9_213")$l9_213[1000]

  as.numeric(x$Vc[300,])
  coef(mod.gdina, "delta")$`Item 1`
  print(paste0("l1_0::",l1_0, "; l1_11: ",l1_11))
  coef(mod.gdina, "delta")$`Item 2`
  print(paste0("l2_0::",l2_0, "; l2_11: ",l2_11))
  coef(mod.gdina, "delta")$`Item 3`
  print(paste0("l3_0::",l3_0, "; l3_12: ",l3_12))

  coef(mod.gdina, "delta")$`Item 9`
  print(paste0("l9_0::",l9_0, "; l9_11: ",l9_11,
               "; l9_13: ", l9_13, "; l9_213: ", l9_213))
})


