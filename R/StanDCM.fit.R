#' @title A function to generate leave-one-out cross-validation for Stan Model
#'
#' @description
#' The StanLCDM.loofit Function to automate Stan code geneartion for LCDMs with binary resposnes
#'
#' @param Qmatrix the Q-matrix specified for the LCDM
#' @param savepath save the .stan file to somewhere; the default path is getwd()
#' @param savename name the .stan
#' @return a. stan file saved at the specified path
#' @importFrom rstan extract
#' @author {Zhehan Jiang, University of Alabama, \email{zjiang17@@ua.edu}}
#'
#' @export

StanDCM.fit <- function(x, pars ="log_lik", cores = 2, save_psis =TRUE) {
  loo1 <- loo(x = x, pars = pars, save_psis = save_psis, cores = cores)
  loo1
}
