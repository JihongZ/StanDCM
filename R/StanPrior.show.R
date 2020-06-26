#' @title Generate Stan code and Run the estimation for LCDM
#'
#' @description
#' The StanPrior.show Function to automate Stan code geneartion for LCDMs with binary resposnes
#'
#' @param script.path the path of .stan file
#' @return a. stan file saved at the specified path
#'
#' @author {Zhehan Jiang, University of Alabama, \email{zjiang17@@ua.edu}}
#' @importFrom plyr failwith id summarize count desc mutate arrange rename summarise
#' @import stringr CDM
#' @export

StanPrior.show<-function(script.path){
  # Install.package("plyr")
  # Install.package('stringr')
  readStan.script<- readLines(script.path)
  readStan.script<-readStan.script[grep("Prior",readStan.script):grep("Likelihood",readStan.script)]
  readStan.script<-readStan.script[grep("~",readStan.script)]
  readStan.script<-str_remove_all(readStan.script," ")
  readStan.script<-str_remove_all(readStan.script,";")
  readStan.script
}


