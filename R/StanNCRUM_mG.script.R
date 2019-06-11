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

StanNCRUM_mG.script<-function(Qmatrix,
                              group.num,
                              fixeditem.vector=NA,
                              class.equal=T,
                              save.path=getwd(),save.name="NCRUM_uninf_multiG"){

  #Load packages
  Install.package("plyr")
  Install.package('stringr')

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
  Kernel.exp.detect<-OUTPUT[[1]] #052719updates
  Kernel.exp.NCRUM<-OUTPUT[[1]] #052719updates
  for (i in 1:nrow(OUTPUT[[1]])){
    for ( j in 1:ncol(OUTPUT[[1]])){
      if(sum(grep('S',OUTPUT[[1]][i,j]))!=0){Kernel.exp[i,j]<-gsub('S','+',OUTPUT[[1]][i,j])
      Kernel.exp.detect[i,j]<-NA} #052719updates
    }
  }
  for (i in 1:nrow(OUTPUT[[1]])){ #052719updates
    theClosestEffect<-which(is.na(Kernel.exp.detect[i,]))[1] #052719updates
    useToReplaceLonger<-Kernel.exp[i,theClosestEffect] #052719updates
    Kernel.exp.NCRUM[i,is.na(Kernel.exp.detect[i,])]<-useToReplaceLonger #052719updates
  } #052719updates

  #Monotonicity constraint in terms of the interaction terms of the item effects
  Constrain.List1<-NULL
  name.inter<-unlist(OUTPUT[[3]])[unlist(OUTPUT[[4]])>=2]
  numway.inter<-unlist(OUTPUT[[4]])[unlist(OUTPUT[[4]])>=2]
  subname.inter<-substr((unlist(OUTPUT[[3]])[unlist(OUTPUT[[4]])>=2]), (nchar(unlist(OUTPUT[[3]])[unlist(OUTPUT[[4]])>=2])-unlist(OUTPUT[[4]])[unlist(OUTPUT[[4]])>=2]+1),
                        nchar(unlist(OUTPUT[[3]])[unlist(OUTPUT[[4]])>=2]))

  if(length(name.inter)!=0){
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
  }else{
    Constrain.List1<-NULL
  }

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


  #########053019update:create pistar and rstar parameters###########################
  pistarParm<-rep(0,nr)
  for(loopi in 1:nr){
    pistarParm[loopi]<-paste('  piParm[',loopi,']=inv_logit(',paste(Kernel.exp[loopi,nclass]),');\n',sep='')
  }

  piParm<-str_replace_all(intercept,"l","pi");piParm<-str_replace_all(piParm,"_0","")
  rParm<-c(itemParmName[1:numMainEffect],name.inter)
  rParm<-str_replace_all(rParm,"l","r")
  Kernel.exp.NCRUM<-Kernel.exp
  for(i in 1:Nc){Kernel.exp.NCRUM[,i]<-piParm}
  #for(i in 1:nr){
  #  for (j in 1:Nc){
  #   Kernel.exp.NCRUM.inv[i,j]<-str_replace_all(Kernel.exp[i,j],
  #                                             intercept[i],
  #                                            piParm[i])
  #}
  #}
  #############################################################
  ###########060719update:The rParm update#####################

  allrParm<-NULL
  for(i in 1:nc){
    allrParm<-cbind(allrParm,paste(paste("r",1:nr,"_1",sep=''),i,sep='')
    )
  }
  allrParm[!t(apply(allrParm ,1,function(x){x%in%rParm}))]<-1
  noneOneallrParm<-allrParm[allrParm!='1']
  noneOneAllrParm<-noneOneallrParm
  noneOneAllrParm<-paste(' real',noneOneAllrParm, ';\n',sep=' ')

  for(i in 1:Nc){
    for(j in 1:nr){
      Kernel.exp.NCRUM[j,i]<-paste(paste(allrParm[j,which(strsplit(colnames(Kernel.exp)[i],'')[[1]]=='0')],"*",collapse ='',sep=''),
                                   Kernel.exp.NCRUM[j,i],sep='')
    }
  }
  Kernel.exp.NCRUM[,Nc]<-str_replace_all(c(Kernel.exp.NCRUM[,Nc]), "[*]", '')

  piTransf<-Kernel.exp.NCRUM[,Nc]
  piTransfwithoutName<-piTransf
  for(loopi in 1:nr){
    piTransf[loopi]<-paste('  piParm[',loopi,']=inv_logit(',paste(Kernel.exp[loopi,Nc]),');\n',sep='')
    piTransfwithoutName[loopi]<-paste('inv_logit(',paste(Kernel.exp[loopi,Nc]),')',sep='')
  }

  rTransf<-NULL
  for(loopi in 1:nr){
    for(loopc in 1:Nc){
      if(str_count(Kernel.exp.NCRUM[loopi,loopc],"\\*")==1){
        if(strsplit(Kernel.exp.NCRUM[loopi,loopc],"\\*")[[1]][1]=='1'){next}
        temp_r<-noneOneallrParm[(grep(strsplit(Kernel.exp.NCRUM[loopi,loopc],"\\*")[[1]][1],noneOneallrParm))]
        if(temp_r%in%noneOneallrParm){
          rTransf<-c(rTransf,paste(" ",temp_r,"=",paste( paste('inv_logit(',Kernel.exp[loopi,loopc],")",sep=''),
                                                         "/",
                                                         piTransfwithoutName[loopi])
          ))
        }
      }
    }
  }
  rTransf<-paste(rTransf,');\n',sep='')
  rTransf<-str_replace_all(rTransf,"\\)\\);","\\);")

  #############################################################
  #######060319update: Multiple Group##########################
  #############################################################
  if(!is.na(fixeditem.vector)[1]){
    fixedItem.vector<-c(1:nr)[-fixeditem.vector]
  }else{fixedItem.vector<-c(1:nr)}

  #############################################################
  #######060719update: Multiple Group##########################
  #############################################################
  Kernel.exp.NCRUM<-Kernel.exp
  for(loopi in 1:nr){
    for( loopc in 1:nclass){
      Reparm[loopi,loopc]<-paste('  PImat[',loopi,',',loopc,']=inv_logit(',paste(Kernel.exp[loopi,loopc]),');\n',sep='')
    }
  }

  #############################################################
  #######060719update: Multiple Group End######################
  #############################################################

  Modelcontainer<-paste('   vector[Nc] contributionsC;\n','    vector[Ni] contributionsI;\n\n',sep='')
  Parmprior<-paste(c(paste('   //Prior\n'),paste('   ',itemParmName,'~normal(0,15)',';\n',sep=''),paste('   Vc~dirichlet(rep_vector(2.0, Nc));',sep='')))  #############################################################


  Kernel.exp.NCRUM.groupName<-paste("Kernel.exp.NCRUM_g",c(1:group.num),sep='')
  for(i in 1:group.num){
    tempfill.Kernel.exp.NCRUM<-Kernel.exp.NCRUM
    temp.Kernel.exp.NCRUM<-Kernel.exp.NCRUM[fixedItem.vector,]
    for(j in 1:nrow(temp.Kernel.exp.NCRUM)){
      for(z in 1:ncol(temp.Kernel.exp.NCRUM)){
        temp.Kernel.exp.NCRUM[j,z]<-paste(temp.Kernel.exp.NCRUM[j,z],'_g',i,sep='')
      }
    }
    for(j in 1:nrow(temp.Kernel.exp.NCRUM)){
      for(z in 1:ncol(temp.Kernel.exp.NCRUM)){
        temp.Kernel.exp.NCRUM[j,z]<-str_replace_all(temp.Kernel.exp.NCRUM[j,z],"\\+",paste("_g",i,"+",sep=''))
      }
    }
    tempfill.Kernel.exp.NCRUM[fixedItem.vector,]<-temp.Kernel.exp.NCRUM
    assign(Kernel.exp.NCRUM.groupName[i],tempfill.Kernel.exp.NCRUM)
  }
  Kernel.exp.NCRUM.list<-list()
  for(i in 1:group.num){Kernel.exp.NCRUM.list[[i]]<-eval(parse(text=paste("Kernel.exp.NCRUM_g",i,sep='')))}
  #############################################################
  ##########060719update: Multiple Group End###################
  #############################################################

  #############################################################
  ##########060719update: Multiple Group#######################
  #############################################################

  PImat.groupName<-paste("PImat_g",c(1:group.num),sep='')
  Reparm.multigroup<-array(0,dim = c(nr,nclass,group.num))
  #Produce Stan code for PImat parameter
  for(loopi in 1:nr){
    for( loopc in 1:nclass){
      for (loopg in 1:group.num){
        Reparm.multigroup[loopi,loopc,loopg]<-paste('  PImat_g',loopg,'[',loopi,',',loopc,']=inv_logit(',paste(Kernel.exp.NCRUM.list[[loopg]][loopi,loopc]),');\n',sep='')

      }
    }
  }

  #############################################################
  ##########060319update: Multiple Group ######################
  #############################################################
  piTransf.multigroup<-array(0,dim = c(nr,1,group.num))
  rTransf.multigroup<-array(0,dim = c(length(rTransf),1,group.num))
  for(i in 1:group.num){
    tempfill.piTransf<-piTransf
    tempfill.rTransf<-rTransf

    tempfill.piTransf<-str_replace_all(tempfill.piTransf,"piTransf",paste("piParm_g",i,sep=''))
    #tempfill.rTransf<-str_replace_all(tempfill.rTransf,"rTransf",paste("rTransf_g",i,sep=''))
    for (j in 1:length(noneOneallrParm)){
      for( z in 1:length(tempfill.rTransf)){
        tempfill.rTransf[z]<-str_replace_all(tempfill.rTransf[z],noneOneallrParm[j],paste(noneOneallrParm[j],"_g",i,sep=''))

      }
    }

    temp.piTransf<-piTransf
    fixedItem.vector_rTransf<-NULL;
    for(j in fixedItem.vector){
      fixedItem.vector_rTransf<-c(fixedItem.vector_rTransf,
                                  grep(j,unlist(lapply(strsplit(unlist(lapply(strsplit(rTransf,"="),function(x){x[[1]]})),"_"),function(x){x[[1]]}))))
    }

    temp.rTransf<-rTransf[fixedItem.vector_rTransf]
    #060819update:a crtical step to handle piTransf lx_0 doesnt have group but some main effects have grouped difference
    Kernel.exp.fill<-Kernel.exp
    Kernel.exp.fill<-unique(Kernel.exp.fill)
    for (w in 1:nr){
      for( j in 1:nclass){
        for (z in 1:nr){
          #Kernel.exp.fill[w,j]<-str_replace_all(Kernel.exp.fill[w,j],paste(intercept[z],"\\+",sep=''),"")
          Kernel.exp.fill[w,j]<-str_replace_all(Kernel.exp.fill[w,j],intercept[z],"")

        }
      }
    }
    Kernel.exp.fill<-strsplit(paste(unlist(c(Kernel.exp.fill[fixedItem.vector,])),collapse="+"),"\\+")[[1]]
    Kernel.exp.fill<-Kernel.exp.fill[Kernel.exp.fill!='']
    Kernel.exp.fill<-unique(Kernel.exp.fill)
    Kernel.exp.fill<-c(Kernel.exp.fill,intercept[fixedItem.vector])
    for(j in 1:length(temp.piTransf)){
      temp.piTransf[j]<-str_replace_all(temp.piTransf[j],"piParm",paste("piParm_g",i,"",sep=''))
    }
    for(z in 1:length(Kernel.exp.fill)){
      for(j in 1:length(temp.piTransf)){
        temp.piTransf[j]<-str_replace_all(temp.piTransf[j],Kernel.exp.fill[z],paste(Kernel.exp.fill[z]
                                                                                    ,"_g",i,sep=''))
      }
    }



    for(j in 1:length(temp.rTransf)){
      #    temp.piTransf<-piTransf
      temp.rTransf[j]<-str_replace_all(temp.rTransf[j],"\\)",paste("_g",i,")",sep=''))
      temp.rTransf[j]<-str_replace_all(temp.rTransf[j],"\\+",paste("_g",i,"+",sep=''))
    }


    for (j in 1:length(noneOneallrParm)){
      for( z in 1:length(temp.rTransf)){
        temp.rTransf[z]<-str_replace_all(temp.rTransf[z],noneOneallrParm[j],paste(noneOneallrParm[j],"_g",i,sep=''))

      }
    }


    tempfill.piTransf[]<-temp.piTransf
    tempfill.rTransf[fixedItem.vector_rTransf]<-temp.rTransf

    piTransf.multigroup[,,i]<-tempfill.piTransf
    rTransf.multigroup[,,i]<-tempfill.rTransf

  }

  #############################################################
  ##########060319update: Multiple Group ######################
  #############################################################

  intercept.multigroup<-array(0,dim = c(nr,1,group.num))
  mainEff.multigroup<-array(0,dim = c(numMainEffect,1,group.num))
  interaction.multigroup<-array(0,dim = c(length(name.inter),1,group.num))
  #Group Invariant Parameter Name
  fixedParmName<-NULL
  if(!is.na(fixeditem.vector)[1]){
    for(i in fixeditem.vector){
      fixedParmName<-c(fixedParmName,out[[3]][[i]])
    }
    fixedParmName<-c(fixedParmName,out[[5]][fixeditem.vector])
  }
  #Group Variant Parameter Name
  freeParmName<-c(out[[5]],itemParmName)[!c(out[[5]],itemParmName)%in%fixedParmName]

  for(i in 1:group.num){
    tempfill.intercept<-out[[5]]
    tempfill.mainEff<-itemParmName[1:numMainEffect]
    tempfill.interaction<-name.inter

    temp.intercept<-tempfill.intercept[fixedItem.vector]
    temp.mainEff<-tempfill.mainEff[!tempfill.mainEff%in%fixedParmName]
    temp.interaction<-tempfill.interaction[!tempfill.interaction%in%fixedParmName]

    for(j in 1:length(temp.intercept)){
      temp.intercept[j]<-paste(temp.intercept[j],"_g",i,sep='')
    }

    for(j in 1:length(temp.mainEff)){
      temp.mainEff[j]<-paste(temp.mainEff[j],"_g",i,sep='')
    }

    for(j in 1:length(temp.interaction)){
      temp.interaction[j]<-paste(temp.interaction[j],"_g",i,sep='')
    }

    tempfill.intercept[fixedItem.vector]<-temp.intercept
    tempfill.mainEff[!tempfill.mainEff%in%fixedParmName]<-temp.mainEff
    tempfill.interaction[!tempfill.interaction%in%fixedParmName]<-temp.interaction

    intercept.multigroup[,,i]<-tempfill.intercept
    mainEff.multigroup[,,i]<-tempfill.mainEff
    interaction.multigroup[,,i]<-tempfill.interaction
  }

  Constrain.List<-paste('  real<lower=0>',unique(mainEff.multigroup),';\n ')
  Constrain.List<-c(Constrain.List,paste('  real<lower=0>',unique(interaction.multigroup),';\n '))
  Unconstrain.List<-paste('  real',unique(c(intercept.multigroup)),';\n ')
  #############################################################
  ##########060319update: Multiple Group End###################
  #############################################################


  Modelcontainer<-paste('   vector[Nc] contributionsC;\n','    vector[Ni] contributionsI;\n\n',sep='')
  Parmprior<-paste(c(paste('   //Prior\n'),paste('   ',itemParmName,'~normal(0,15)',';\n',sep=''),paste('   Vc~dirichlet(rep_vector(2.0, Nc));',sep='')))
  update.Parmprior<-Parmprior
  fix.Parmprior<-NULL
  update.Parmprior<-update.Parmprior[update.Parmprior!='']

  #############################################################
  #####060319update:Change from update.Parmprior&fix     ######
  #############################################################
  update.Parmprior.multiGroup<-NULL
  for(i in 1:length(update.Parmprior)){
    for (j in 1:length(freeParmName)){
      for (z in 1:group.num){
        if(sum(grepl(freeParmName[j],update.Parmprior[i]))>=1){
          temp.update.Parmprior<-str_replace_all(update.Parmprior[i],freeParmName[j],paste(freeParmName[j],"_g",z,sep=''))
          update.Parmprior.multiGroup<-c(update.Parmprior.multiGroup,
                                         temp.update.Parmprior)
        }
      }
    }
  }
  update.Parmprior.multiGroup<-unique(update.Parmprior.multiGroup)
  update.Parmprior.multiGroup<-c("   //Prior\n",paste('   ',fixedParmName,'~normal(0,15)',';\n',sep=''),update.Parmprior.multiGroup )
  if(class.equal){
    update.Parmprior.multiGroup<-c(update.Parmprior.multiGroup,paste('   Vc~dirichlet(rep_vector(2.0, Nc));',sep='') )
  }else{
    for(i in 1:group.num){
      update.Parmprior.multiGroup<-c(update.Parmprior.multiGroup,
                                     paste('   Vc_g',i,'~dirichlet(rep_vector(2.0, Nc));\n',sep='') )
    }
  }

  ##therefore we can use: fix.Parmprior,update.Parmprior
  #############################################################
  #############################################################

  #############################################################
  #####060419update:Likelihood Add PImat_g#####################
  #############################################################
  PImat.likelihood1<-NULL
  PImat.likelihood0<-NULL
  Vc.likelihood<-NULL
  for(loopg in 1:group.num){
    temp.PImat.likelihood1<-paste(paste('          if (GroupID[iterp]==',loopg,')',sep=''),
                                  paste('            contributionsI[iteri]=bernoulli_lpmf(1|PImat_g',loopg,'[iteri,iterc]);',sep=''),
                                  sep='\n')
    PImat.likelihood1<-paste(PImat.likelihood1,temp.PImat.likelihood1,sep='\n')
    temp.PImat.likelihood0<-paste(paste('          if (GroupID[iterp]==',loopg,')',sep=''),
                                  paste('            contributionsI[iteri]=bernoulli_lpmf(0|PImat_g',loopg,'[iteri,iterc]);',sep=''),
                                  sep='\n')
    PImat.likelihood0<-paste(PImat.likelihood0,temp.PImat.likelihood0,sep='\n')
  }
  if(!class.equal){
    for(loopg in 1:group.num){
      temp.Vc.likelihood<-paste(paste('       if (GroupID[iterp]==',loopg,')',sep=''),
                                paste('         contributionsC[iterc]=log(Vc_g',loopg,'[iterc])+sum(contributionsI);',sep=''),
                                sep='\n')
      Vc.likelihood<-paste(Vc.likelihood,temp.Vc.likelihood,sep='\n')
    }
  }

  #Likelihood Stan code
  if(class.equal){
    Likelihood<-paste('
                      \n
                      //Likelihood
                      for (iterp in 1:Np){
                      for (iterc in 1:Nc){
                      for (iteri in 1:Ni){
                      if (Y[iterp,iteri] == 1)'
                      ,PImat.likelihood1,'\n',
                      '        else'
                      ,PImat.likelihood0,
                      '}
                      contributionsC[iterc]=log(Vc[iterc])+sum(contributionsI);
                      }
                      target+=log_sum_exp(contributionsC);
                      }
                      ',sep='')}else{
                        Likelihood<-paste('
                                          \n
                                          //Likelihood
                                          for (iterp in 1:Np){
                                          for (iterc in 1:Nc){
                                          for (iteri in 1:Ni){
                                          if (Y[iterp,iteri] == 1)'
                                          ,PImat.likelihood1,'\n',
                                          '        else'
                                          ,PImat.likelihood0,
                                          '}\n',
                                          Vc.likelihood
                                          ,'


                                          }
                                          target+=log_sum_exp(contributionsC);
                                          }
                                          ',sep='')

                      }
  #############################################################
  #####060419update:Likelihood Add PImat_g  end################
  #############################################################

  #Data Specification
  data.spec<-'
  data{
  int Np;
  int Ni;
  int Nc;
  matrix[Np, Ni] Y;
  vector[Np] GroupID;
  }
  '


  #############################################################
  #####060419update: Stan script##############################
  #############################################################
  Constrain.List<-unique(Constrain.List);Unconstrain.List<-unique(Unconstrain.List)
  #Parameter Specification
  if(class.equal){parm.spec<-paste(c('
                                     parameters{
                                     simplex[Nc] Vc;\n ',paste0(Constrain.List),paste0(Unconstrain.List),
                                     '}\n'),collapse='')}else{
                                       parm.spec<-paste(c('
                                                          parameters{\n ',
                                                          paste(paste('   simplex[Nc] Vc_g',1:group.num, ";",sep=''),"\n"),
                                                          paste0(Constrain.List),paste0(Unconstrain.List),
                                                          '}\n'),collapse='')


                                     }
  #To set rParm to real value
  rParm.real<-
    paste('real', str_replace_all(unique(unlist(lapply(strsplit(rTransf.multigroup,"="),function(x){x[[1]]})))," ",""),';\n',sep=' ')

  #Reparameter Specification
  tranrTransf.spec<-paste(c('
                            transformed parameters{

                            ',
                            paste('  matrix[Ni, Nc] PImat_g',1:group.num,';\n',sep=''),
                            paste('  vector[Ni] piParm_g',1:group.num,';\n',sep=''),
                            rParm.real,

                            unique(c(piTransf.multigroup)), #060419update
                            unique(c(rTransf.multigroup)), #060419update
                            paste0(c(Reparm.multigroup)),'}\n'),collapse='')

  #Model Specification update052619
  model.spec<-paste(c('\nmodel {\n',paste(c(Modelcontainer,update.Parmprior.multiGroup,Likelihood),sep=''),'\n}',sep=''))
  model.spec<-model.spec[!startsWith(str_remove_all(model.spec," "),"~")]

  #Generated Quantities Specification
  IC.generatedquantities<-NULL
  if(!class.equal){
    for(loopg in 1:group.num){
      temp.IC.generatedquantities<-paste(paste('       if (GroupID[iterp]==',loopg,')',sep=''),
                                         paste('         contributionsIC[iteri,iterc]=log(Vc_g',loopg,'[iterc])+contributionsI[iteri];',sep=''),
                                         sep='\n')
      IC.generatedquantities<-paste(IC.generatedquantities,temp.IC.generatedquantities,sep='\n')
    }
  }


  if(class.equal){
    generatedQuantities.spec<-paste('
                                    \n
                                    generated quantities {
                                    vector[Ni] log_lik[Np];
                                    vector[Ni] contributionsI;
                                    matrix[Ni,Nc] contributionsIC;
                                    //Posterior
                                    for (iterp in 1:Np){
                                    for (iteri in 1:Ni){
                                    for (iterc in 1:Nc){
                                    if (Y[iterp,iteri] == 1)'
                                    ,PImat.likelihood1,'\n',
                                    '        else'
                                    ,PImat.likelihood0,
                                    '
                                    contributionsIC[iteri,iterc]=log(Vc[iterc])+contributionsI[iteri];
                                    }
                                    log_lik[iterp,iteri]=log_sum_exp(contributionsIC[iteri,]);
                                    }
                                    }
                                    }
                                    ',sep='')}else{
                                      generatedQuantities.spec <- paste( '
                                                                         \n
                                                                         generated quantities {
                                                                         vector[Ni] log_lik[Np];
                                                                         vector[Ni] contributionsI;
                                                                         matrix[Ni,Nc] contributionsIC;
                                                                         //Posterior
                                                                         for (iterp in 1:Np){
                                                                         for (iteri in 1:Ni){
                                                                         for (iterc in 1:Nc){
                                                                         if (Y[iterp,iteri] == 1)',
                                                                         PImat.likelihood1,'\n',
                                                                         '        else',
                                                                         PImat.likelihood0,'\n',
                                                                         "   ",IC.generatedquantities,'\n
                                                                         }\n',
                                                                         '     log_lik[iterp,iteri]=log_sum_exp(contributionsIC[iteri,]);
                                                                         }
                                                                         }
                                                                         }
',
                                                                         sep = ''
                                      )
                                    }




  if (.Platform$OS.type == "unix") {
    filename = paste(paste(save.path,save.name,sep='/'),'.stan',sep='')
  }else{
    filename = paste(paste(save.path,save.name,sep='\\'),'.stan',sep='')
  }

  sink(file=filename, append=FALSE)
  cat(
    paste(c('   ',
            data.spec,parm.spec,tranrTransf.spec,model.spec,generatedQuantities.spec)
    ))
  sink(NULL)

}

