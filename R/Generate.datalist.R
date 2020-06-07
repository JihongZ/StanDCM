#' @title Generate and prepare the data list used for this Packgae
#'
#' @description
#' \code{Generate.datalist} is used to generate a list for estimation.
#'
#' @param Qmatrix the Q-matrix specified for the LCDM
#' @param response.matrix The response matrix without any missing values. Each cell represents the response of
#' each item for each person.
#' @param GroupID default is NA. Used for multi-group LDCM.
#' @return a list including Y (response matrix), Np (number of observation) and Ni (number of iterm)
#'
#' @author {Zhehan Jiang, University of Alabama, \email{zjiang17@@ua.edu}}
#' @export
#' @examples
#' \dontrun{
#' Qmatrix<-cbind(Qmatrix,rep(1,9));Qmatrix[1,1]<-0
#' }





Generate.datalist<-function(Qmatrix,response.matrix,GroupID=NA){
  if(sum(is.na(response.matrix))!=0){
    stop('Stop!The response dataset contains missing value(s)')
    }
  else{
    if( length((unique(unlist(unique(c(response.matrix))))))>2 ){
      Generate.dataList<-list(Y=response.matrix, Na=ncol(Qmatrix),
                              Np = nrow(response.matrix),
                              Ni = ncol(response.matrix),Nc=2^(ncol(Qmatrix)),
                              Ns = length((unique(unlist(unique(c(response.matrix)))))))
    }else{
      Generate.dataList<-list(Y=response.matrix, Na=ncol(Qmatrix),
                              Np = nrow(response.matrix),
                              Ni = ncol(response.matrix),Nc=2^(ncol(Qmatrix)))
    }
  }
  if(!is.na(GroupID)[1]){Generate.dataList$GroupID=GroupID}
  return(Generate.dataList)
}
