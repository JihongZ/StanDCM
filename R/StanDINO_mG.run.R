#' @title Generate Stan code and Run the estimation for ORDM
#'
#' @description
#' The StanLCDM.script Function to automate Stan code geneartion for LCDMs with binary resposnes
#'
#' @param Qmatrix the Q-matrix specified for the LCDM
#' @param save.path save the .stan file to somewhere; the default path is getwd()
#' @param save.name name the .stan
#' @return a. stan file saved at the specified path
#'
#' @author {Zhehan Jiang, University of Alabama, \email{zjiang17@@ua.edu}}
#'
#' @export
#loading needed packages
#load("D:\\Dropbox\\Stan\\R\\Data")

StanDINO_mG.run<-function(Qmatrix,
                          response.matrix,
                          GroupID,
                          fixeditem.vector=NA,
                          class.equal=T,
                          script.path=NA,save.path=getwd(),save.name="DINO_uninf_multiG",
                          iter=1000,warmup = 0,
                          chain.num=3,init.list='random',control.list=NA){
  group.num<-length(factor(GroupID))
  rstan.detect<-tryCatch(library("rstan"),error=function(e){"rstan is not loaded properly. See https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started for details."})
  if(length(rstan.detect)==1){
    break
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
    ))),collapse=',')  ,')',collapse='')

    IniList1<-NULL
    if(!class.equal){
      temp.inilist1<-eval(parse(n =2000000 ,text=inilist1))
      eval(parse(n =2000000 ,text=paste(paste('IniList1$Vc_g',1:group.num,sep=''),"<-temp.inilist1$Vc",sep='') ))
      inilist1<-IniList1
    }else{
      inilist1<-eval(parse(n =2000000 ,text=inilist1))}


    for( i in 2:chain.num){
      temp.text<-paste('inilist',i,"<-inilist1",sep='')
      eval(parse(text=(temp.text)))
    }
    temp.text<-paste('init.list<-list(',paste(paste('inilist',1:chain.num,sep=''),collapse = ","),')',sep='')
    eval(parse(text=(temp.text)))
  }
  data.list<-Generate.datalist(Qmatrix,response.matrix,GroupID)

  if(is.na(control.list)){control.list<-list(adapt_delta=0.82)}
  if(is.na(script.path)==T){
    options(warn=-1)
    #Need to update script
    StanDINO_mG.script(Qmatrix=Qmatrix,
                       group.num=group.num,
                       fixeditem.vector=fixeditem.vector,
                       class.equal=class.equal,
                       save.path=save.path,save.name=save.name)
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

