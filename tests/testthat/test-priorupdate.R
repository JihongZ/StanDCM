context("test-priorupdate")

test_that("check prior update", {
  Qmatrix2 <- cbind(Qmatrix,rep(1,9))
  Qmatrix2[1,1]<-0
  StanPrior.update(priorUpdate.matrix = Qmatrix2, save.name = "Update_DINO", script.path = "DINO_uninf.stan")
  StanPrior.update(priorUpdate.matrix = Qmatrix2, save.name = "Update_CRUM", script.path = "CRUM_uninf.stan")
  StanPrior.update(priorUpdate.matrix = Qmatrix2, save.name = "Update_DINA", script.path = "DINA_uninf.stan")
  StanPrior.update(priorUpdate.matrix = Qmatrix2, save.name = "Update_DINO_uninf_multiG", script.path = "DINO_uninf_multiG.stan")
  StanPrior.update(priorUpdate.matrix = Qmatrix2, save.name = "Update_LCDM", script.path = "LCDM_uninf.stan")
  StanPrior.update(priorUpdate.matrix = Qmatrix2, save.name = "Update_NCRUM", script.path = "NCRUM_uninf.stan")
  StanPrior.update(priorUpdate.matrix = Qmatrix2, save.name = "Update_ORDM", script.path = "ORDM_uninf.stan")
})

