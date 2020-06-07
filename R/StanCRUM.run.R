#' @title Generate Stan code and Run the estimation for ORDM
#'
#' @description
#' The StanCRUM.run Function to automate Stan code generation for LCDMs with binary responses
#' @author {Zhehan Jiang, University of Alabama, \email{zjiang17@@ua.edu}}
#'
#' @param Qmatrix the Q-matrix specified for the LCDM
#' @param response.matrix the response matrix
#' @param script.path the path to save the .stan file to somewhere; the default path is getwd()
#' @param save.path the path to save the .stan file to somewhere; the default path is getwd()
#' @param save.name name the .stan
#' @param iter name the .stan
#' @param warmup number of warmup iterations
#' @param chain.num number of MCMC chains. Default is 3.
#' @param init.list initial variables. Default is random. Other options include cdm
#' @param control.list the control values
#'
#' @return a. stan file saved at the specified path
#' @importFrom GDINA GDINA
#' @importFrom stringr str_count
#' @importFrom rstan extract
#' @importFrom stats coef
#' @export
#' @examples
#' \dontrun{
#' #----------- CRUM model-----------#
#' mod1<-StanCRUM.run(Qmatrix, response.matrix, iter=2000, init.list='cdm', chain.num = 3)
#' summary(mod1)
#' }

StanCRUM.run <- function(Qmatrix, response.matrix, script.path = NA, save.path = getwd(), save.name = "CRUM_uninf", iter = 1000, warmup = 0,
                         chain.num = 3, init.list = "random", control.list = NA) {
  if(!sum(installed.packages()%in%"rstan")){
    stop("rstan is not loaded properly. See https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started for details.")
  }
  Cdm.init <- F
  if (init.list == "cdm") {
    Cdm.init <- T
    Install.package(c("CDM", "stringr"))
    trueParmName <- Parm.name(Qmatrix = Qmatrix)$parm.name
    Classp.exp1 <- Parm.name(Qmatrix = Qmatrix)$class.expression
    mod1 <- GDINA::GDINA(data = response.matrix, q.matrix = Qmatrix, maxit = 700, link = "logit", progress = F)
    CDMresult <- as.data.frame(coef(mod1))

    CDM.parm.name <- paste(paste(paste("l", CDMresult[, 3], sep = ""), "_", sep = ""), str_count(CDMresult$partype.attr, "Attr"), sep = "")
    CDM.parm.name <- paste(CDM.parm.name,
      unlist(lapply(strsplit(unlist(lapply(strsplit(CDMresult$partype.attr, "Attr", fixed = FALSE), function(x) {
        paste(x, collapse = "")
      })), "-"), function(x) {
        paste(x, collapse = "")
      })),
      sep = ""
    )
    CDM.parm.est <- CDMresult$est
    parm.ini <- round(CDM.parm.est[match(trueParmName, CDM.parm.name)], 4)
    CDM.prop.est <- mod1$attribute.patt
    prop.ini <- CDM.prop.est[match(Classp.exp1, rownames(CDM.prop.est)), 1]
    inilist1 <- paste("list(", paste(noquote(paste(noquote(unlist(list(
      paste(paste("Vc=c(", paste((prop.ini), collapse = ","), ")", collapse = ","))
    ))))), collapse = ","), ")", collapse = "")

    inilist1 <- eval(parse(n = 2000000, text = inilist1))
    for (i in 2:chain.num) {
      temp.text <- paste("inilist", i, "<-inilist1", sep = "")
      eval(parse(text = (temp.text)))
    }
    temp.text <- paste("init.list<-list(", paste(paste("inilist", 1:chain.num, sep = ""), collapse = ","), ")", sep = "")
    eval(parse(text = (temp.text)))
  }
  data.list <- Generate.datalist(Qmatrix, response.matrix)

  if (is.na(control.list)) {
    control.list <- list(adapt_delta = 0.82)
  }
  if (is.na(script.path) == T) {
    options(warn = -1)
    StanDINA.script(Qmatrix, save.path = save.path, save.name = save.name)
    script.path <- paste(paste(save.path, save.name, sep = "/"), ".stan", sep = "")
    options(warn = 0)
    compiled_model <- stan_model(script.path)
  } else {
    compiled_model <- stan_model(script.path)
  }
  if (Cdm.init == T) {
    estimated_model <- tryCatch(sampling(compiled_model,
      data = data.list,
      iter = iter,
      init = init.list,
      warmup = warmup,
      chains = chain.num,
      control = control.list
    ),
    error = function(e) {
      "The estimation process is terminated with errors"
    }
    )
  } else {
    estimated_model <- tryCatch(sampling(compiled_model,
      data = data.list,
      iter = iter,
      init = init.list,
      warmup = warmup,
      chains = chain.num,
      control = control.list
    ),
    error = function(e) {
      "The estimation process is terminated with errors"
    }
    )
  }

  estimated_model
}

