#' @title A function to generate leave-one-out cross-validation for Stan Model
#'
#' @description
#' The StanLCDM.loofit Function to automate Stan code geneartion for LCDMs with binary resposnes
#'
#' @param Qmatrix the Q-matrix specified for the LCDM
#' @param savepath save the .stan file to somewhere; the default path is getwd()
#' @param savename name the .stan
#' @return a. stan file saved at the specified path
#'
#' @author {Zhehan Jiang, University of Alabama, \email{zjiang17@@ua.edu}}
#'
#' @export
#loading needed packages

# load("data.RData")
# Qmatrix<-cbind(Qmatrix,rep(1,9)); Qmatrix[1,1]<-0
# # dim(respMatrix)
# misspecifiedQmatrix <- Qmatrix
# misspecifiedQmatrix[1:6,] <- 1-Qmatrix[1:6,]
#
# # misspecifiedQmatrix[1,3] = 0
# #

# estimatedMod_misspecified<-StanLCDM.run(misspecifiedQmatrix,respMatrix, iter=100,init.list='cdm',chains = 3)
#
# loo1 <- loo(estimatedMod, pars = "log_lik", save_psis = TRUE)
# loo2 <- loo(estimatedMod_misspecified, pars = "log_lik", save_psis = TRUE)
# print(loo1)
# print(loo2)
#
# save(Qmatrix, misspecifiedQmatrix, respMatrix, file = "R/data.RData")
# save(estimatedMod, estimatedMod_misspecified, file = "R/modelfit.RData")
#
# loo::compare(loo1, loo2)

StanDCM.fit <- function(x, pars ="log_lik", cores = 2, save_psis =TRUE) {
  loo1 <- loo(x = x, pars = pars, save_psis = save_psis, cores = cores)
  loo1
}
