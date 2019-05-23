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
#' @author {Zhehan Jiang, University of Alabama, \email{zjiang17@@ua.edu}
#'
#' @export


StanLCDM.script<-function(Qmatrix,savepath=getwd(),savename="LCDM_uninf"){
  
  #Packages detection
  install.packages.auto <- function(x) { 
    x <- as.character(substitute(x)) 
    if(isTRUE(x %in% .packages(all.available=TRUE))) { 
      eval(parse(text = sprintf("require(\"%s\")", x)))
    } else { 
      #update.packages(ask= FALSE) #update installed packages.
      eval(parse(text = sprintf("install.packages(\"%s\", dependencies = TRUE)", x)))
    }
    if(isTRUE(x %in% .packages(all.available=TRUE))) { 
      eval(parse(text = sprintf("require(\"%s\")", x)))
    } else {
      source("http://bioconductor.org/biocLite.R")
      #biocLite(character(), ask=FALSE) #update installed packages.
      eval(parse(text = sprintf("biocLite(\"%s\")", x)))
      eval(parse(text = sprintf("require(\"%s\")", x)))
    }
  }
  
  #Load packages
  install.packages.auto("plyr")
  install.packages.auto('stringr')

  
  nc<-ncol(Qmatrix)
  nr<-nrow(Qmatrix)
  temp.table.col<-unique(apply(combn(rep(c(0,1),nc),nc),2,function(x){paste(x,collapse = "")}))
  temp.table.col<-temp.table.col[order(temp.table.col)]
  temp.table<-matrix(0,nr,length(temp.table.col))
  colnames(temp.table)<-temp.table.col
  rownames(temp.table)<-paste('item',c(1:nr),sep='')
  temp.table<-as.data.frame(temp.table)
  for (i in 1:nr){
    temp.table[i,]<-paste('l',i,'_0',sep='')
  }
  intercept<-temp.table[,1]
  
  #Generate attribute combinations
  comb.generator<-function(x.vector){
    if(length(x.vector)>1){
      temp.attr<-x.vector
      temp.attr.sav<-NULL
      for(i in 1:length(temp.attr)){
        temp.1<-combn(temp.attr,i)
        temp.2<-apply(temp.1,2,function(x){paste(x,collapse = "")})
        temp.attr.sav<-c(temp.attr.sav,temp.2)
      }
    }
    if(length(x.vector)==1){temp.attr.sav<-x.vector}
    temp.attr.sav
  }
  #vectors needed for combination.generator
  Item.load.id<-list()
  for ( i in 1:nr){
    Item.load.id[[i]]<-grep('1',Qmatrix[i,])}

  Attr.load.id<-list()
  attr.load.id<-matrix(0,length(temp.table.col),nc)
  for ( i in 1:length(temp.table.col)){
    attr.load.id[i,]<-unlist(strsplit(temp.table.col[i],split=''))
    Attr.load.id[[i]]<-grep('1',attr.load.id[i,])
  }
  
  #Generate Combination for both Item.load and Attr.load
  Item.Comb<-list()
  for ( i in 1:nr){
    Item.Comb[[i]]<-comb.generator(Item.load.id[[i]])
  }
  Attr.Comb<-list()
  for ( i in 2:length(temp.table.col)){
    Attr.Comb[[1]]<-0
    Attr.Comb[[i]]<-comb.generator(Attr.load.id[[i]])
  }
  constraints.list<-list()
  nway.inter.list<-list()
  for(i in 1:nr){
    for(a in 2:length(temp.table.col)){
      ifzero<-as.numeric(paste(Item.Comb[[i]][Item.Comb[[i]]%in%(Attr.Comb[[a]])],collapse=''))
      if((!is.na(ifzero))){
        temp.table[i,a]<-paste(c(temp.table[i,a],
                                 paste("S","l",i,"_",nchar(Item.Comb[[i]][Item.Comb[[i]]%in%(Attr.Comb[[a]])]),Item.Comb[[i]][Item.Comb[[i]]%in%(Attr.Comb[[a]])],sep='',collapse='')
        ),collapse='')
        if(a==length(temp.table.col)){
          nway.inter.list[[i]]<-nchar(Item.Comb[[i]][Item.Comb[[i]]%in%(Attr.Comb[[a]])])
          constraints.list[[i]]<-paste("l",i,"_",nchar(Item.Comb[[i]][Item.Comb[[i]]%in%(Attr.Comb[[a]])]),Item.Comb[[i]][Item.Comb[[i]]%in%(Attr.Comb[[a]])],sep='')
        }
      }
    }
  }
  
  #Create Lambda Table
  Lamda.Table<-temp.table
  for(i in 1:nr){
    for(a in 1:length(Lamda.Table)){
      t.ref<-unique(as.character(Lamda.Table[i,]))
      pos<-c(1:length(t.ref))[Lamda.Table[i,a]==t.ref]
      temp.table[i,a]<-paste("t",i,"_",pos,sep='')}}
  
  #Generate LCDM specification
  out<-list() 
  out[[1]]<-Lamda.Table
  out[[2]]<-temp.table
  out[[3]]<-constraints.list
  out[[4]]<-nway.inter.list
  out[[5]]<-intercept
  OUTPUT<-out
  nclass<-ncol(OUTPUT[[1]]);Nc<-nclass

  #Produce kernel expressions across items and attributes 
  Kernel.exp<-OUTPUT[[1]] 
  for (i in 1:nrow(OUTPUT[[1]])){
    for ( j in 1:ncol(OUTPUT[[1]])){
      if(sum(grep('S',OUTPUT[[1]][i,j]))!=0){Kernel.exp[i,j]<-gsub('S','+',OUTPUT[[1]][i,j])}
    }
  }


  #Monotonicity constraint in terms of the interaction terms of the item effects
  Constrain.List1<-NULL
  name.inter<-unlist(OUTPUT[[3]])[unlist(OUTPUT[[4]])>=2]
  numway.inter<-unlist(OUTPUT[[4]])[unlist(OUTPUT[[4]])>=2]
  subname.inter<-substr((unlist(OUTPUT[[3]])[unlist(OUTPUT[[4]])>=2]), (nchar(unlist(OUTPUT[[3]])[unlist(OUTPUT[[4]])>=2])-unlist(OUTPUT[[4]])[unlist(OUTPUT[[4]])>=2]+1),
                      nchar(unlist(OUTPUT[[3]])[unlist(OUTPUT[[4]])>=2]))

  for (inter in 1: length(name.inter)){
    temp.nw<-numway.inter[inter]
    temp.nm<-name.inter[inter]
    temp.subnm<-strsplit(subname.inter[inter],split='')[[1]]
    temp.sel<-paste(unlist(strsplit(temp.nm,split = '_'))[1],"_",(1:(temp.nw-1)),sep='')
    first.sel<-unlist(OUTPUT[[3]])[grep(paste((temp.sel),collapse="|"),unlist(OUTPUT[[3]]))]
    second.sel<-sub(".*_.", "", first.sel)
    for (sel in 1:length(temp.subnm)){
      SEL<-second.sel[sel]
      Constrain.List1<-rbind(
        paste(temp.nm,">-(0", paste("+",first.sel[grep(SEL,second.sel)],
                                  sep='',collapse=''),")",sep=''),Constrain.List1)
    }
  }
  Constrain.List1<-as.character(Constrain.List1)

  itemParmName<-c(unlist(OUTPUT[[3]])[unlist(OUTPUT[[4]])==1],unlist(OUTPUT[[3]])[unlist(OUTPUT[[4]])==2],
                  unlist(OUTPUT[[3]])[unlist(OUTPUT[[4]])==3],
                  unlist(OUTPUT[[3]])[unlist(OUTPUT[[4]])==4],
                  unlist(OUTPUT[[3]])[unlist(OUTPUT[[4]])==5],
                  unlist(OUTPUT[[3]])[unlist(OUTPUT[[4]])==6],
                  unlist(OUTPUT[[3]])[unlist(OUTPUT[[4]])==7],
                  unlist(OUTPUT[[3]])[unlist(OUTPUT[[4]])==8],
                  unlist(OUTPUT[[3]])[unlist(OUTPUT[[4]])==9],
                  unlist(OUTPUT[[3]])[unlist(OUTPUT[[4]])==10],
                  unlist(OUTPUT[[3]])[unlist(OUTPUT[[4]])==11],OUTPUT[[5]]) 
  numMainEffect<-length(unlist(OUTPUT[[3]])[unlist(OUTPUT[[4]])==1])
  Constrain.List<-paste('  real<lower=0>',itemParmName[1:numMainEffect],';\n ')
  Unconstrain.List<-paste('  real',itemParmName[-(1:numMainEffect)],';\n ')
  Reparm<-as.data.frame(matrix(0,nr,nclass))
  

  #Produce Stan code for PImat parameter
  for(loopi in 1:nr){
    for( loopc in 1:nclass){
      Reparm[loopi,loopc]<-paste('  PImat[',loopi,',',loopc,']=inv_logit(',paste(Kernel.exp[loopi,loopc]),');\n',sep='')
    }
  }  

  Modelcontainer<-paste('   vector[Nc] contributionsC;\n','    vector[Ni] contributionsI;\n\n',sep='')
  Parmprior<-paste(c(paste('   //Prior\n'),paste('   ',itemParmName,'~normal(0,15)',';\n',sep=''),paste('   Vc~dirichlet(rep_vector(2.0, Nc));',sep='')))
  #Likelihood Stan code
  Likelihood<-'
  \n
  //Likelihood
  for (iterp in 1:Np){
    for (iterc in 1:Nc){
      for (iteri in 1:Ni){
        if (Y[iterp,iteri] == 1)
          contributionsI[iteri]=bernoulli_lpmf(1|PImat[iteri,iterc]);
        else
          contributionsI[iteri]=bernoulli_lpmf(0|PImat[iteri,iterc]); 
      }
      contributionsC[iterc]=log(Vc[iterc])+sum(contributionsI);
    }
  target+=log_sum_exp(contributionsC);    
  }
  '


  #Data Specification
  data.spec<-'
data{
  int Np;
  int Ni;
  int Nc;
  matrix[Np, Ni] Y;
}
  '

  #Parameter Specification
  parm.spec<-paste(c('
parameters{
  simplex[Nc] Vc;\n ',paste0(Constrain.List),paste0(Unconstrain.List),
'}\n'),collapse='')

  #Reparameter Specification
  transparm.spec<-paste(c('
  transformed parameters{
  matrix[Ni, Nc] PImat;\n',
  paste0(unlist(Reparm)),'}\n'),collapse='')

  #Model Specification
  model.spec<-paste(c('\nmodel {\n',paste(c(Modelcontainer,Parmprior,Likelihood),sep=''),'\n}',sep=''))

  #Generated Quantities Specification
  generatedQuantities.spec<-'
  \n
generated quantities {  
  
  vector[Np] log_lik[Nc, Ni];
  matrix[Np,Nc] contributionsPC;
  vector[Ni] contributionsI;
  //Posterior
  for (iterp in 1:Np){
      for (iterc in 1:Nc){
        for (iteri in 1:Ni){
          if (Y[iterp,iteri] == 1)
            contributionsI[iteri]=bernoulli_lpmf(1|PImat[iteri,iterc]);
          else
            contributionsI[iteri]=bernoulli_lpmf(0|PImat[iteri,iterc]);
          log_lik[iterp,iterc,iteri]=contributionsI[iteri];
        }
        contributionsPC[iterp,iterc]=prod(exp(contributionsI));
      }
    }
}
  '
  
  sink(file=paste(paste(savepath,savename,sep='\\'),'.stan',sep=''),append=FALSE)
    cat(
    paste(c('   ',
    data.spec,parm.spec,transparm.spec,model.spec,generatedQuantities.spec)
  ))
  sink(NULL)
}

StanLCDM.run<-function(Qmatrix,script.path=NA,savepath=getwd(),savename="LCDM_uninf"){
  rstan.detect<-tryCatch(library("rstan"),error=function(e){"rstan is not loaded properly. See https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started for details."})
  if(startsWith(rstan.detect,"rstan")){
    break
  }
  if(is.na(script.path)==T){
    StanLCDM.script(Qmatrix,savepath=savepath,savename=savename)
    script.path<-paste(paste(savepath,savename,sep='\\'),'.stan',sep='')
    compiled_model<-stan_model(script.path)
  }else{
    compiled_model<-stan_model(script.path)
  }
  
}
