#' Function that start estimation and prediction tasks
#'
#' @param fixed_X0.models fixed effects in the submodel for the baseline level of processes
#' @param fixed_DeltaX.models a two-sided linear formula object for specifying the response outcomes (one the left part of ~ symbol) 
#' and the covariates with fixed-effects (on the right part of ~ symbol) 
#' in the submodel for change over time of latent processes
#' @param randoms_X0.models random effects in the submodel for the baseline level of processes
#' @param randoms_DeltaX.models random effects in the submodel for change over time of latent processes
#' @param mod_trans.model model for elements of the temporal transition matrix, which captures 
#' the temporal influences between latent processes
#' @param DeltaT indicates the discretization step
#' @param outcomes indicates names of the outcomes
#' @param nD number of the latent processes
#' @param mapping.to.LP indicates which outcome measured which latent process, it is a mapping table between
#'  outcomes and latents processes
#' @param link indicates link used to transform outcome
#' @param knots indicates position of knots used to transform outcomes 
#' @param subject indicates the name of the covariate representing the grouping structure
#' @param data indicates the data frame containing all the variables for estimating the model
#' @param Time indicates the name of the covariate representing the time
#' @param Survdata Dataset for event and statusEvent
#' @param basehaz baseline hazard type (Weibull or splines)
#' @param knots_surv knots for splines (baseline hazard function)
#' @param assoc association type between outcomes and times-to-events
#' @param truncation boolean indicating if left truncation or not
#' @param fixed.survival.models fixed effects in the submodel for the event hazards
#' @param makepred indicates if predictions in the real scales of outcomes have to be done
#' @param MCnr number of replicates  to compute the predictions in the real scales of the outcomes
#' @param type_int type of Monte Carlo integration method to use
#' @param paras.ini initial values for parameters, default values is NULL
#' @param indexparaFixeUser position of parameters to be constrained
#' @param paraFixeUser values associated to the index of parameters to be constrained
#' @param maxiter maximum iteration
#' @param univarmaxiter maximum iteration for estimating univariate model 
#' @param nproc number of processor to be used for running this package, default value is 1
#' @param epsa threshold for the convergence criterion on the parameters, default value is 1.e-4
#' @param epsb threshold for the convergence criterion on the likelihood, default value is 1.e-4
#' @param epsd threshold for the convergence criterion on the derivatives, default value is 1.e-3
#' @param print.info  to print information during the liklihood optimization, default value is FALSE 
#' @param \dots optional parameters
#'
#' @return CInLPN2 object
CInLPN2.default <- function(fixed_X0.models, fixed_DeltaX.models, randoms_X0.models, randoms_DeltaX.models, mod_trans.model, 
                            DeltaT, outcomes, nD, mapping.to.LP, link, knots=NULL, subject, data, Time, 
                            Survdata = NULL, basehaz = NULL, knots_surv=NULL, assoc = 0, truncation = FALSE, fixed.survival.models = NULL, 
                            makepred, MCnr, type_int = NULL, sequence = NULL, ind_seq_i = NULL, nmes = NULL, cholesky= FALSE,
                            paras.ini= NULL, indexparaFixeUser, paraFixeUser, maxiter, zitr, ide, univarmaxiter, nproc = 1, 
                            epsa =0.0001, epsb = 0.0001, epsd= 0.001, print.info = FALSE, TimeDiscretization = TRUE, 
                            Tentry = NULL, Event = NULL, StatusEvent = NULL, ...)
{
  cl <- match.call()

  ################### discretization of the data with discretisation value given by the user ##########################
  #

  if(TimeDiscretization){#PB with covariates !
    data <- TimeDiscretization(rdata=data, subject = subject, 
                               fixed_X0.models = fixed_X0.models,
                               randoms_X0.models = randoms_X0.models, fixed_DeltaX.models = fixed_DeltaX.models, 
                               randoms_DeltaX.models = randoms_DeltaX.models, mod_trans.model = mod_trans.model,
                               outcomes = outcomes, #predictors = predictors, 
                               Time = Time, Delta = DeltaT)
  }
  ################### created formated data ##########################
  
  if(requireNamespace("survival", quietly = TRUE)){
    data_F <- DataFormat(data=data, subject = subject, fixed_X0.models = fixed_X0.models,
                         randoms_X0.models = randoms_X0.models, fixed_DeltaX.models = fixed_DeltaX.models, 
                         randoms_DeltaX.models = randoms_DeltaX.models, mod_trans.model = mod_trans.model, 
                         outcomes = outcomes, nD = nD, link=link, knots = knots, zitr= zitr, ide = ide, 
                         Time = Time, Survdata = Survdata, basehaz = basehaz, fixed.survival.models =fixed.survival.models, DeltaT=DeltaT, assoc = assoc, truncation = truncation)
    
  }else{
    stop("Need package survival to work, Please install it.")
  }

  K <- data_F$K #  number of markers
  vec_ncol_x0n <- data_F$vec_ncol_x0n # number of parameters on initial level of processes
  n_col_x <- ncol(data_F$x) # number of parameters on processes slope
  nb_RE <- data_F$nb_RE # number of random effects on the processes slope
  L <- ncol(data_F$modA_mat)
  ncolMod.MatrixY <- ncol(data_F$Mod.MatrixY)
  
  #definition quasi-random sequence
  if(!is.null(type_int)){
    if(type_int == "sobol"){
      sequence  <- randtoolbox::sobol(n = MCnr, dim = nb_RE, scrambling = 1, normal = TRUE, init=T)
    }else if(type_int == "halton"){
      sequence  <- randtoolbox::halton(n = MCnr, dim = i, normal = TRUE, init=T) 
    }else if(type_int == "torus"){
      sequence  <- randtoolbox::torus(n = MCnr, dim = i,normal = TRUE, init=T) 
    }
  }
  
  
  assocT <- NULL
  if(!is.null(assoc)){
    assocT <- ifelse(assoc==0, "r.intercept",ifelse(assoc==1, "r.slope",ifelse(assoc==2, "r.intercept/slope",ifelse(
      assoc==3, "c.value",ifelse(assoc==4, "c.slope","c.value/slope")
    ))))
  }
  
  # ### creation of arguments:  Initialising parameters
  if(K>1 & is.null(paras.ini)){
    
    paras.ini <- f_paras.ini(data = data, outcomes = outcomes, mapped.to.LP = mapping.to.LP, fixed_X0.models = fixed_X0.models, fixed_DeltaX.models = fixed_DeltaX.models,  
                             randoms_DeltaX.models = randoms_DeltaX.models, nb_RE = nb_RE, mod_trans.model = mod_trans.model, 
                             subject = subject, MCnr = MCnr, type_int = type_int,
                             Survdata = Survdata, basehaz = basehaz, Time = Time, link = link, knots = knots, zitr = zitr, ide = ide,
                             DeltaT = DeltaT, maxiter = univarmaxiter, epsd = epsd, nproc = nproc, print.info = print.info, 
                             TimeDiscretization = TimeDiscretization, fixed.survival.models = fixed.survival.models, 
                             Tentry = Tentry, Event = Event, StatusEvent = StatusEvent, assocT = assocT, truncation = truncation)
  }
  npara_k <- sapply(outcomes, function(x) length(grep(x, names(data.frame(data_F$Mod.MatrixY)))))
  
  paras <- Parametre(K=K, nD = nD, vec_ncol_x0n, n_col_x, nb_RE, indexparaFixeUser = indexparaFixeUser, 
                     paraFixeUser = paraFixeUser, L = L, ncolMod.MatrixY = ncolMod.MatrixY, paras.ini=paras.ini, 
                     link = link, npara_k = npara_k, 
                     Survdata = Survdata, basehaz = basehaz, knots_surv = knots_surv, assoc = assoc, truncation = truncation,
                     data = data, outcomes = outcomes, df= data_F$df, nE = data_F$nE, np_surv = data_F$np_surv, fixed.survival.models =fixed.survival.models,)
  if_link <- rep(0,K)
  for(k in 1:K){
    if(!link[k] %in%c("linear","thresholds")){
      if_link[k] <- 1
    }else if(link[k]=="thresholds"){
      if_link[k] <- 2
    } 
  }
  paras$npara_k <- npara_k
  
  #add zitr, ide  dans estim(). What about knots? What in dataF
  if(any(link=="thresholds")|| !is.null(Survdata) || !is.null(type_int)){
    #  nmes <- c()
    #  for (i in 1:length(unique(data$id))){
    #    for(k in 1:K){
    #      ind = which(names(data)==outcomes[k])
    #      nmes <- c(nmes, length(which(data$id==unique(data$id)[i] & !is.na(data[,ind]))))
    #    }
    #  }
    # nmes <- unique(nmes)
    # nmes <- nmes[order(nmes)]
    # 
    # sequence <- matrix(NA, MCnr*length(nmes), max(nmes))
    # 
    # for(j in 1:length(nmes)){
    #   if (type_int == "sobol") {
    #     Seq <- randtoolbox::sobol(MCnr, dim = nmes[j], normal = TRUE, scrambling = 1)
    #   } else if (type_int == "halton") {
    #     Seq <- randtoolbox::halton(MCnr, dim = nmes[j], normal = TRUE)
    #   }
    #   sequence[((j-1)*MCnr+1):(j*MCnr), 1:nmes[j]] <- Seq
    # }
    
    paras$sequence <- sequence
    paras$type_int <- ifelse(type_int=="halton",1,ifelse(type_int=="sobol",2,ifelse(type_int=="torus",3,ifelse(type_int=="MC",-1,0))))
    paras$ind_seq_i <- ind_seq_i
  }else{
    paras$type_int <- -1
    paras$sequence <- matrix(0,nD,1)
    paras$ind_seq_i <- 0
  }
  
  if(!is.null(Survdata)){
    paras$basehaz <- basehaz
    if(typeof(as.matrix(Survdata))!="double")
      stop("Remove the factor labels in data so that the matrix is all numeric.")
  }else{
    Survdata<-0
  }
  
  # estimation
  est <- CInLPN2.estim(K = K, nD = nD, mapping.to.LP = mapping.to.LP, data = data_F, if_link = if_link, cholesky = cholesky,
                       DeltaT = DeltaT, MCnr = MCnr, nmes = nmes, data_surv = Survdata, 
                       paras = paras, maxiter = maxiter, nproc = nproc, epsa = epsa, epsb = epsb,
                       epsd = epsd, print.info = print.info)
  
  
  res <- list(conv = est$istop, v = est$v, best = est$b, ca = est$ca, cb = est$cb, rdm = est$rdm, 
              #niter = est$iter, 
              niter = est$ni, coefficients = est$coefficients, posfix = est$posfix)
  
  #   # fitted value and correlation matrix
  #   m <- length(data$tau)
  res$fitted.values <- NULL
  if(makepred & res$conv==1){
    ptm2 <- proc.time() 
    cat("Computation of the predictions for the observed marker(s) \n")
    cat(paste("based on Newton Raphson algorithm (criterion: eps=1.e-9) and Monte Carlo integration with N=",MCnr),"replicates \n")
    #        cat(    ,"and  for search of inverse of a transformation function \n ")
    #    cat(paste("Number on replicates for MC integration : N=",MCnr),"\n")
    #    cat("Convergence criterion for the seach of inverse : eps=1.e-9 \n")
    cat("Execution may take some time \n ")
    
    #================== predictions
    
    res$Marginal_Predict <- data_F$id_and_Time
    res$SubjectSpecific_Predict <- data_F$id_and_Time
    col <- colnames(res$Marginal_Predict)
    # colSS <- colnames(res$SubjectSpecific_Predict)
    if(requireNamespace("splines2", quietly = TRUE)){
      Predict <- pred(K = K, nD = nD, mapping = mapping.to.LP, paras = res$coefficients,
                      m_is= data_F$m_i, Mod_MatrixY = data_F$Mod.MatrixY, df= data_F$df,
                      x = data_F$x, z = data_F$z, q = data_F$q, cholesky = cholesky, nb_paraD = data_F$nb_paraD, x0 = data_F$x0, z0 = data_F$z0,
                      q0 = data_F$q0, if_link = if_link, tau = data_F$tau,
                      tau_is=data_F$tau_is, modA_mat = data_F$modA_mat, DeltaT=DeltaT, 
                      MCnr = MCnr, minY=data_F$minY, maxY=data_F$maxY, knots=data_F$knots, data_F$degree, epsPred = 1.e-9)
    }else{
      stop("Need package MASS to work, Please install it.")
    }
    kk <- 1
    for(k in 1: K){
      res$Marginal_Predict <- cbind(res$Marginal_Predict,data_F$Y[,k],Predict[,kk], 
                                    (data_F$Y[,k]-Predict[,kk]),Predict[,(kk+1):(kk+3)])
      
      res$SubjectSpecific_Predict <- cbind(res$SubjectSpecific_Predict,data_F$Y[,k],Predict[,(kk+4)], 
                                           (data_F$Y[,k]-Predict[,(kk+4)]),Predict[,c((kk+1),(kk+5),(kk+6))])
      
      col <- c(col,outcomes[k],paste(outcomes[k], "Pred", sep="."), paste(outcomes[k], "Res", sep="."),
               paste(outcomes[k], "tr", sep="-"), paste(outcomes[k], "tr.Pred", sep="-"),
               paste(outcomes[k], "tr.Res", sep="-"))
      kk <- kk+7
    }
    colnames(res$Marginal_Predict) <- col
    colnames(res$SubjectSpecific_Predict) <- col
    
    p.time2 <- proc.time() - ptm2
    cat("Prediction computation took:", p.time2[1], "\n")
  }
  
  ## Compute number of parameters  per component of the processes network
  i1 <- 0
  res$length_para_mu0 <- length(which(res$posfix[(i1+1):(i1+ncol(data_F$x0))]==0))
  i1 <- i1 + ncol(data_F$x0)
  res$length_para_mu <- length(which(res$posfix[(i1+1):(i1+ncol(data_F$x))]==0))
  i1 <- i1 + ncol(data_F$x)
  res$length_para_RE <- length(which(res$posfix[(i1+1):(i1+data_F$nb_paraD)]==0))
  res$length_para_trY = ncol(data_F$Mod.MatrixY)
  
  ###names of estimates parameters 
  #colname transition matrix
  a <- NULL
  
  for(i in 1:nD){
    a <- c(a, paste("a",paste(i,1:nD,sep=""), sep="_"))
  }
  Col.matA <- NULL
  colmodA_mat <- colnames(data_F$modA_mat)
  for(i in 1: length(a)){
    Col.matA <- c(Col.matA, paste(a[i],colmodA_mat, sep="."))
  }
  # colname RE
  L <- NULL
  #   for(i in 1:K){
  npRE <- data_F$nb_paraD
  if(cholesky){
    L <- paste("Chol.", 1:npRE,sep="")
  }else{
    L <- paste("Rho.", 1:npRE,sep="")
  }
  
  # colname  measurement error
  sig <- paste("sigma",outcomes, sep=".")
  # Measurement model parameters
  nb_para_trY <- ncol(data_F$Mod.MatrixY)
  ParaTransformY <- colnames((data_F$Mod.MatrixY))
  
  ##
  ##output related to the statistical model
  res$ns <- data_F$nb_subject
  res$N <- data_F$nb_obs
  res$nD <- data_F$nD
  res$K <- data_F$K
  # est$best : only estimates of non constraint parameters
  res$linkstype <- link
  res$linknodes <- data_F$knots
  res$nb_paraD <- data_F$nb_paraD
  res$tau = data_F$tau
  res$outcomes <- outcomes
  res$covariates <- data_F$all.pred
  res$Time <- Time 
  res$DeltaT <- DeltaT
  res$colnames <- c(colnames(data_F$x0), colnames(data_F$x), L, Col.matA, sig,ParaTransformY)
  if(!is.null(Survdata)){
    param_survie <-c()
    if(data_F$nE>0){
      #   if(paras$basehaz=="Weibull"){
      #     param_survie <-c("shape", "scale")
      #   }else{
      #     message("rename survival parameters if not Weibull")
      #   }
      #   Xsurv <- as.matrix(model.matrix(as.formula(paste("",fixed.survival.models, sep="~")),data=Survdata)[,-1])
      #   param_survie <- c(param_survie, names(Xsurv))
      # }else if(data_F$nE==2){
      param_survie<-c()
      for(ij in 1:data_F$nE){
        if(paras$basehaz=="Weibull"){
          param_survie <- c(param_survie, paste("shape",ij,sep="."), paste("scale",ij,sep="."))
        }else{
          message("rename survival parameters if not Weibull")
        }
        
        covsurv <- unique(as.vector(strsplit(fixed.survival.models,"[|*+]")[[1]]))
        if(length(which(covsurv!="1"))>0){
          Xsurv <- as.data.frame(model.matrix(as.formula(paste("",fixed.survival.models[ij], sep="~")),data=Survdata)[,-1])
          for(ip in 1:length(names(Xsurv))){
            param_survie <- c(param_survie, paste(names(Xsurv)[ip],ij,sep="."))
          }   
        }

        
        if(assoc==0){
          name_assoc <- "r.inter"
        }else if(assoc==1){
          name_assoc <- "r.slope"
        }else if(assoc==2){
          name_assoc <- c("r.inter", "r.slope")
        }else if(assoc==3){
          name_assoc <- "c.value"
        }else if(assoc==4){
          name_assoc <- "c.slope"
        }else if(assoc==5){
          name_assoc <- c("c.value", "c.slope")
        }
        if(nD>1)
          name_assoc <- paste(name_assoc, 1:nD,sep=".")
        param_survie <- c(param_survie, name_assoc)
      }
    }
    res$colnames <- c(res$colnames, param_survie)
  }
  
  nb_RE <- data_F$q0+data_F$q
  
  if(cholesky){
    chol <- matrix(0, nb_RE*nD, nb_RE*nD)
    chol[upper.tri(chol, diag = T)] <- res$coefficients[grep("Chol", res$colnames)]
    res$varcov <- t(chol)%*%chol
    
  }else{
    prmea <- matrix(0, nb_RE*nD, nb_RE*nD)
    prmea[upper.tri(prmea,diag=TRUE)] <- res$coefficients[grep("Rho", res$colnames)]
    prmea <- t(prmea)
    prmea[upper.tri(prmea,diag=TRUE)] <- res$coefficients[grep("Rho", res$colnames)]
    
    sea <- abs(diag(prmea))
    
    corr <- (exp(prmea)-1)/(exp(prmea)+1)
    diag(corr) <- 1
    covea <- sweep(corr,1,sea,"*")
    covea <- sweep(covea,2,sea,"*")
    res$varcov <- covea
  }

  res$coefficients <- as.matrix(res$coefficients)
  rownames(res$coefficients) <- res$colnames
  colnames(res$coefficients) <- "Coef."
  
  
  res$loglik <- est$fn.value
  res$AIC <- -2*est$fn.value + 2*length(est$b)
  res$BIC <- -2*est$fn.value + log(res$N)*length(est$b)
  
  
  ##output related to the iteration process
  # res$istop: binary indicator of if convergence is reached
  # res$iter: indicates the number of iterations
  # res$ca : convergence criteria related to the stability of the parameters
  # res$cb : convergence criteria related to the stability of the likelihood
  # res$rdm : convergence criteria related to the inversibility of the Hessian matrix
  
  
  
  ##output related big data like predictions, stocked matrix
  # res$fitted.values
  res$modA_mat = data_F$modA_mat
  res$cholesky <- cholesky
  res$call <- match.call()
  class(res) <- 'CInLPN2'
  res
}
