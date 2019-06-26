#' @title DCM calibration under the DINA model via Stan
#'
#'
#' @description
#' \code{StanDINA} uses Stan program to calibrate the deterministic inputs, noisy and gate model for dichotomous responses, and its
#' extension
#'
#' In addition, users are allowed to specify design matrix and link function for each item, and distinct models may be used
#' in a single test for different items. The attributions can be either dichomous or polytomous.

#' @usage
#' StanDINA.run<-function(Qmatrix,response.matrix,script.path=NA,save.path=getwd(),save.name="DINA_uninf",iter=1000,warmup = floor(iter/2),
#' chain.num=3,init.list='random',control.list=NA)
#'
#' @param Qmatrix A required matrix
#' @param response.matrix save the .stan file to somewhere; the default path is getwd()
#' @param script.path save the .stan file to somewhere; the default path is getwd()
#' @param save.name the name of saved stan file
#' @param iter number of iteration of MCMC estimation. defalts to 1000.
#' @param chain.num number of MCMC chain.num.
#' @param init.list the initial values. 'random' or 'CDM'
#' @param control.list the controlled parameters
#'
#' @return StanDINA returens an object of class StanDINA. Methods for StanDINA objects include
#' \code{\link{extract}} for extract for extracting various components, \code{\link{coef}} for
#' extracting strctural parameters.
#'
#'
#' @author {Zhehan Jiang, University of Alabama, \email{zjiang17@@ua.edu} \cr Jihong Zhang, University of Iowa, \email{jihong-zhang@uiowa.edu}}
#'
#' @export
#' @examples
#' \dontrun{
#' #----------- DINO model-----------#
#' mod1<-StanDINO.run(Qmatrix, respMatrix, iter=20, init.list='cdm', chain.num = 3)
#' summary(mod1)
#' }

StanDINO.run<-function(Qmatrix,response.matrix,script.path=NA,save.path=getwd(),save.name="DINO_uninf",iter=1000,warmup = 0,
                       chain.num=3,init.list='random',control.list=NA){
  rstan.detect<-tryCatch(library("rstan"),error=function(e){"rstan is not loaded properly. See https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started for details."})
  if(length(rstan.detect)==1){
    stop()
  }
  Cdm.init<-F
  if(init.list=='cdm'){
    Cdm.init<-T
    Install.package(c("CDM","stringr"))
    trueParmName<-Parm.name(Qmatrix=Qmatrix)$parm.name
    Classp.exp1<-Parm.name(Qmatrix=Qmatrix)$class.expression
    mod1<-gdina( data =respMatrix, q.matrix = Qmatrix , maxit=700,link = "logit",progress=F)
    CDMresult<-as.data.frame(coef(mod1))
    library(stringr)
    CDM.parm.name<-paste(paste(paste('l',CDMresult[,3],sep=''),'_',sep=''),str_count(CDMresult$partype.attr,"Attr"),sep='')
    CDM.parm.name<-paste(CDM.parm.name,
                         unlist(lapply(strsplit(unlist(lapply(strsplit(CDMresult$partype.attr, 'Attr', fixed=FALSE),function(x){paste(x,collapse="")})),'-'),function(x){paste(x,collapse="")})),
                         sep='')
    CDM.parm.est<-CDMresult$est
    parm.ini<-round(CDM.parm.est[match(trueParmName,CDM.parm.name)],4)
    CDM.prop.est<-mod1$attribute.patt
    prop.ini<-CDM.prop.est[match(Classp.exp1,rownames(CDM.prop.est)),1]
    inilist1<-paste('list(',paste(noquote(paste(noquote(unlist(list(
      paste(paste('Vc=c(',paste((prop.ini),collapse=','),')',collapse=','))))
    ))),collapse=',') ,')',collapse='')

    inilist1<-eval(parse(n =2000000 ,text=inilist1))
    for( i in 2:chain.num){
      temp.text<-paste('inilist',i,"<-inilist1",sep='')
      eval(parse(text=(temp.text)))
    }
    temp.text<-paste('init.list<-list(',paste(paste('inilist',1:chain.num,sep=''),collapse = ","),')',sep='')
    eval(parse(text=(temp.text)))
  }
  data.list<-Generate.datalist(Qmatrix,response.matrix)

  if(is.na(control.list)){control.list<-list(adapt_delta=0.82)}
  if(is.na(script.path)==T){
    options(warn=-1)
    StanDINA.script(Qmatrix,save.path=save.path,save.name=save.name)
    script.path<-paste(paste(save.path,save.name,sep='/'),'.stan',sep='')
    options(warn=0)
    compiled_model<-stan_model(script.path)
  }else{
    compiled_model<-stan_model(script.path)
  }
  if(Cdm.init==T){
    estimated_model<-tryCatch(sampling(compiled_model,
                                       data = data.list,
                                       iter = iter,
                                       init = init.list,
                                       warmup = warmup,
                                       chains=chain.num,
                                       control=control.list),
                              error=function(e){"The estimation process is terminated with errors"})
  }else{
    estimated_model<-tryCatch(sampling(compiled_model,
                                       data = data.list,
                                       iter = iter,
                                       init = init.list,
                                       warmup = warmup,
                                       chains=chain.num,
                                       control=control.list),
                              error=function(e){"The estimation process is terminated with errors"})

  }

  estimated_model
}
