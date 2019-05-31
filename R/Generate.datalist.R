#' @title Generate Stan code and Run the estimation for LCDM
#'
#' @description
#' The StanLCDM.script Function to automate Stan code geneartion for LCDMs with binary resposnes
#'
#' @param Qmatrix the Q-matrix specified for the LCDM
#' @param savepath save the .stan file to somewhere; the default path is getwd()
#' @param savename name the .stan
#' @return a. stan file saved at the specified path
#'
#' @author {Zhehan Jiang, University of Alabama, \email{zjiang17@@ua.edu}}
#' @export
#loading needed packages
#load("D:\\Dropbox\\Stan\\R\\data.RData") ;Qmatrix<-cbind(Qmatrix,rep(1,9));Qmatrix[1,1]<-0


Generate.datalist<-function(Qmatrix,response.matrix){
  if(sum(is.na(response.matrix))!=0){'Stop!The response dataset contains missing value(s)'}else{
    if( length((unique(unlist(unique(c(response.matrix))))))>2 ){
      list(Y=response.matrix, Na=ncol(Qmatrix),
           Np = nrow(response.matrix),
           Ni = ncol(response.matrix),Nc=2^(ncol(Qmatrix)),
           Ns = length((unique(unlist(unique(c(response.matrix)))))))
    }else{
      list(Y=response.matrix, Na=ncol(Qmatrix),
           Np = nrow(response.matrix),
           Ni = ncol(response.matrix),Nc=2^(ncol(Qmatrix)))
    }
  }
}
