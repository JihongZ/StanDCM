#' @title Install and load packages.
#'
#' @description
#' \code{Install.package} allows other functions to install and load packages.
#'
#' @param needed_packages the packages to be installed and loaded.
#'
#' @author {Zhehan Jiang, University of Alabama, \email{zjiang17@@ua.edu}}
#'
#' @examples
#' \dontrun{
#' Install.package("ggplot2")
#' }
#'

Install.package = function(needed_packages){
  for (i in 1:length(needed_packages)){
    haspackage = require(needed_packages[i], character.only = TRUE)
    if (haspackage==FALSE){
      install.packages(needed_packages[i])
      require(needed_packages[i], character.only = TRUE)
    }
  }
}
