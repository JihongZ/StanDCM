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
#' @author {Zhehan Jiang, University of Alabama, \email{zjiang17@@ua.edu}
#'
#' @export
#loading needed packages
#load("D:\\Dropbox\\Stan\\R\\data.RData")


StanLCDM.loofit <- function(x, pars ="log_lik", cores, save_psis) {

}
