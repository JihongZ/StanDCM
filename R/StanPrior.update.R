#' @title Generate Stan code and Run the estimation for LCDM
#'
#' @description
#' The StanPrior.update Function to automate update prior
#'
#' @param priorUpdate.matrix the Q-matrix specified for the LCDM
#' @param script.path path to save the script file
#' @param save.path save the .stan file to somewhere; the default path is getwd()
#' @param save.name name the .stan
#' @return a. stan file saved at the specified path
#'
#' @author {Zhehan Jiang, University of Alabama, \email{zjiang17@@ua.edu}}
#'
#' @export
#loading needed packages
#load("D:\\Dropbox\\Stan\\R\\data.RData") ;Qmatrix<-cbind(Qmatrix,rep(1,9));Qmatrix[1,1]<-0

StanPrior.update <- function(priorUpdate.matrix, script.path, save.path = getwd(), save.name = NA) {
  Install.package("plyr")
  Install.package("stringr")
  if (is.na(save.name)) {
    saveUpdate.name <- paste(format(Sys.time(), "%a%b%Y"), "_priorUpdate.stan", sep = "")
  } else {
    saveUpdate.name <- paste(paste(save.name, sep = "/"), ".stan", sep = "")
  }
  readStan.script <- readLines(script.path)
  for (i in 1:nrow(priorUpdate.matrix)) {
    readStan.script[which(str_remove_all(str_remove_all(readStan.script, " "), ";") == priorUpdate.matrix[i, 1])] <- paste(
      "   ",
      priorUpdate.matrix[i, 2],
      ";//UpdatedPrior "
    )
  }
  filename <- paste(paste(save.path, saveUpdate.name, sep = "/"), sep = "")
  writeLines(readStan.script, filename)
}


