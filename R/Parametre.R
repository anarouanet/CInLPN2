#' Function to initialize parameters k in multivariate CInLPN2 model
#'
#' @param K number of the markers
#' @param nD number of the latent processes
#' @param vec_ncol_x0n vector of number of columns of model.matrix for baseline's submodel
#' @param n_col_x number of overall columns of model.matrix for change's submodel
#' @param nb_RE number of random effects
#' @param stochErr indicates if the structural model contain stochastique components
#' @param indexparaFixeUser position of parameters to be constrained
#' @param paraFixeUser values associated to the index of parameters to be constrained
#' @param L number of columns of model.matrix for temporal infuences model
#' @param paras.ini initial values for parameters, default values is NULL
#' @param ncolMod.MatrixY vector of number of columns of model.matrix for transformation submodel
#'
#' @return a list
#' 
#' 
Parametre <- function(K, nD, vec_ncol_x0n, n_col_x, nb_RE, stochErr=FALSE, indexparaFixeUser =NULL,
                      paraFixeUser=NULL, L = 1, paras.ini, ncolMod.MatrixY, link, npara_k, 
                      Survdata = NULL, basehaz = NULL, knots_surv = NULL, assoc = NULL, truncation = F,
                      data, outcomes, df, nE = 0, np_surv = 0, fixed.survival.models = NULL, interactionY.survival.models = NULL, nYsurv = 0){
  cl <- match.call()
  #   require(MASS)
  #initialisation des parametres
  # L = number of parameters for each coefficient a of matrix A
  # K = number of outcomes
  #======================================================================================

  nb_paraD = nb_RE*(nb_RE+1)/2
  indexparaFixeForIden <- NULL
  # if user not specified initial parameters
  
  if(is.null(basehaz))
    basehaz<-"Weibull" #not to have NULL value in C++ code
  
  if(is.null(paras.ini)){
    p <- 0 # position in the initialize parameters
    cpt1 <- 0 # compteur pour tous les paras
    cpt2<-0 # compteur de boucle
    #alpha_mu0
    alpha_mu0 <-rep(.1,sum(vec_ncol_x0n))
    p <- p+ sum(vec_ncol_x0n)
    index_paraFixe_mu0_constraint <- NULL
    for(n in 1:nD){
      alpha_mu0[(cpt2+1)] <- 0
      cpt2 <- cpt2 + vec_ncol_x0n[n]
      cpt1 <- cpt1 + vec_ncol_x0n[n]
    }
    
    #alpha_mu
    alpha_mu <-rep(.3,n_col_x)
    p <- p+n_col_x
    cpt1 <- cpt1 + n_col_x

    alpha_D <-rep(.1,nb_paraD)
    to_nrow <- nb_RE
    i_alpha_D <- 0
    index_paraFixeDconstraint <- NULL
    for(n in 1:nD){
      alpha_D[i_alpha_D+1] <- 1
      i_alpha_D <- i_alpha_D + to_nrow
      cpt1 <- cpt1 + to_nrow
      to_nrow <- to_nrow -1
    }
    p <- p + nb_paraD
    # para of transition matrix vec_alpha_ij
    vec_alpha_ij <- rep(0.4, L*nD*nD)
    cpt1 <- cpt1 + L*nD*nD
    p <- p + L*nD*nD
    #paraB
    paraB <- NULL
    if(stochErr==TRUE){
      paraB <- rep(.15,nD)
      cpt1 <- cpt1 + nD
      p <- p + nD
    }
    #paraSig
    paraSig <- rep(.5, K)
    cpt1 <- cpt1 + K
    p <- p + K
    ### parameters of the link function
    ParaTransformY <- c()
    for(k in 1:K){
      if(link[k]=='thresholds'){
        ParaTransformY <- c (ParaTransformY, c(2*qunif(0.98)*(-median(data[,outcomes[k]]) + min(data[,outcomes[k]]) +1)/(length(unique(data[,outcomes[k]]))-2),
        rep(sqrt(2*qunif(0.98)/(length(unique(data[,outcomes[k]]))-2)), length(unique(data[,outcomes[k]]))-2)))
        #ParaTransformY <- c(min(data[,outcomes[k]]), rep(1, length(unique(data[,outcomes[k]]))-2))
      }else{
        ParaTransformY <- c(ParaTransformY, rep(1, ncolMod.MatrixY))
      }
    }

    cpt1 <- cpt1 + ncolMod.MatrixY
    p <- p + ncolMod.MatrixY
    
    para_surv <- NULL
    para_basehaz <- NULL
    knots_surv <- c(0,0) # changer !!
    #Survival
    if(!is.null(Survdata)){
      cov_surv <- unique(unlist(strsplit(fixed.survival.models,"[+*]")))
      
      if(nE==2){
        tmat <- trans.comprisk(2, names = c("event-free", "event1", "event2"))
        
        Survdata$stat1 <- as.numeric(Survdata$StatusEvent == 1)
        Survdata$stat2 <- as.numeric(Survdata$StatusEvent == 2)
        Survdatalong <- msprep(time = c(NA, "Event", "Event"), status = c(NA, "stat1", "stat2"), data = Survdata, keep = cov_surv, trans = tmat)
        dim0<-dim(Survdatalong)[2]
        Survdatalong <- expand.covs(Survdatalong, cov_surv)
        cov_survexpanded <- paste(names(Survdatalong)[(dim0+1):dim(Survdatalong)[2]],collapse='+')
        form <- as.formula(paste("Surv(time, status) ~ ",cov_survexpanded, "+ factor(trans)+ strata(trans)", sep=""))
        mod_surv <- survreg(form, data = Survdatalong, dist="weibull")

        #   1/(survreg's scale)  =    rweibull shape a
        #   exp(survreg's intercept) = rweibull scale b
        # (a/b) (x/b)^(a-1) shape a, scale b
        para_basehaz <- c(1/(mod_surv$coefficients[["(Intercept)"]]), exp(mod_surv$icoef[2]),1/(mod_surv$coefficients[["(Intercept)"]]+mod_surv$coefficients[["factor(trans)2"]]), exp(mod_surv$icoef[3]))
        n1 <- 1+np_surv[1]-1+ifelse(assoc%in%c(0, 1, 3, 4),1,2)*nD
        para_surv <- c(mod_surv$coefficients[2:(1+np_surv[1]-1)], rep(0, ifelse(assoc%in%c(0, 1, 3, 4),1,2)*nD), 
                       mod_surv$coefficients[(n1):(n1-1+np_surv[1]-1)], rep(0, ifelse(assoc%in%c(0, 1, 3, 4),1,2)*nD))

      }else if (nE==1){
        
        form <- as.formula(paste("Surv(Event, StatusEvent) ~ ",cov_surv, sep=""))
        mod_surv <- survreg(form, data = Survdata, dist="weibull")
        para_basehaz <- c(1/(mod_surv$coefficients[["(Intercept)"]]), exp(mod_surv$icoef[2]))
        para_surv <- c(mod_surv$coefficients[2:(1+np_surv-1)], rep(0, ifelse(assoc%in%c(0, 1, 3, 4),1,2)*nD))
      }
      p <- p + length(para_basehaz) + length(para_surv)
      np_baz <- length(para_basehaz)/nE
    }

    if(!is.null(interactionY.survival.models)){
      browser()
    }
  }

  # if user specified initial parameters
  if(!is.null(paras.ini)){
    p <- 0 # position in the initialize parameters
    cpt1 <-0 # counter for parameterd
    cpt2<-0 # loop counter
    #alpha_mu0
    alpha_mu0 <- paras.ini[(p+1):(p+sum(vec_ncol_x0n))]
    p <- p+ sum(vec_ncol_x0n)
    index_paraFixe_mu0_constraint <-NULL
    for(n in 1:nD){
      #alpha_mu0[(cpt2+1)] <- 0
      cpt2 <- cpt2 + vec_ncol_x0n[n]
      cpt1 <- cpt1 + vec_ncol_x0n[n]
    }
    paraFixe_mu0_constraint <- rep(1,nD)
    #alpha_mu
    alpha_mu <- paras.ini[(p+1):(p+n_col_x)]
    p <- p+n_col_x
    cpt1 <- cpt1 + n_col_x
    #alpha_D parameters for cholesky of all random effects
    alpha_D <- paras.ini[(p+1):(p+nb_paraD)]
    to_nrow <- nb_RE
    i_alpha_D <- 0
    index_paraFixeDconstraint <- NULL

    for(n in 1:nD){
      #if(link[n] != "thresholds")
      #alpha_D[i_alpha_D+1] <- 1
      i_alpha_D <- i_alpha_D + to_nrow
      cpt1 <- cpt1 + to_nrow
      to_nrow <- to_nrow -1
    }
    p <- p+nb_paraD
    paraFixeDconstraint <- rep(1,nD)
    # para of transition matrix vec_alpha_ij
    vec_alpha_ij <- paras.ini[(p+1):(p + L*nD*nD)]
    p <- p + L*nD*nD
    cpt1 <- cpt1 + L*nD*nD
    # paraB
    paraB <- NULL
    if(stochErr==TRUE){
      
      paraB <- paras.ini[(p+1):(p + nD)]
      p <- p + nD
      cpt1 <- cpt1 + nD
    }
    #paraSig
    paraSig <- paras.ini[(p+1):(p + K)]
    p <- p + K
    cpt1 <- cpt1 + K

    ### para of link function
    ParaTransformY <- paras.ini[(p+1):(p + ncolMod.MatrixY)]
    i_para <- 0
     for(k in 1:K){
       if(link[k]=="linear" & ParaTransformY[i_para+2]==0){
         stop('Second parameter for linear link function cannot be set at 0 (variance)')
       }
       i_para <- i_para + npara_k[k]
    }
    
    cpt1 <- cpt1 + ncolMod.MatrixY
    p <- p + ncolMod.MatrixY
    
    
    #Survival
    para_surv <- NULL
    para_basehaz <- NULL
    knots_surv <- c(0,0) # changer !!
    if(!is.null(Survdata)){
      # if(nE ==1){
      #   np_surv <- dim(Survdata)[2]-3 + ifelse(assoc%in%c(0, 1, 3, 4),1,2)
      # }else{
      #   np_surv <- dim(Survdata)[2]-3 + ifelse(assoc%in%c(0, 1, 3, 4),1,2)
      # }
      np_baz <- ifelse(basehaz=="Weibull",2, 0)# changer 0!!
      for (jj in 1:nE){
        para_basehaz <- c(para_basehaz, paras.ini[(p+1) : (p + np_baz)])  
        p <- p + np_baz  # change here?
      #}
      #for (jj in 1:nE){
        para_surv <- c(para_surv, paras.ini[(p + 1 ) : (p + np_surv[jj])]) 
        p <- p + np_surv[jj] # change here?
      }
      if(basehaz=="Splines") cat('add number of parameters for splines in p and para_surv')
      if(basehaz=="Splines") cat('Define knots_surv para_basehaz')
      
    }
    
    #if(length(paras.ini) != (p + sum(df)))
    #  stop("The length of paras.ini is not correct.")
  }

  #final vector of initial parameters
  paras <- c(alpha_mu0, alpha_mu, alpha_D, vec_alpha_ij,  paraB, paraSig, ParaTransformY)
  t1 <- 0
  t2 <- 0

  
  if(nE>0){
    for(jj in 1:nE){
      paras <- c(paras, para_basehaz[(t1+1) : (t1 + np_baz)]) # change 0!!
      t1 <- t1 + np_baz
      paras <- c(paras, para_surv[(t2 + 1) : (t2 + np_surv[jj])]) # change 0!!
      t2 <- t2 + np_surv[jj]
    }
  }

  # if(nE>0){
  #   for(jj in 1:nE){
  #     paras <- c(paras, para_basehaz[(t1+1) : (t1 + np_baz)]) # change 0!!
  #     t1 <- t1 + np_baz
  #   }
  # 
  #   for(jj in 1:nE){
  #     paras <- c(paras, para_surv[(t2 + 1) : (t2 + np_surv[jj])]) # change 0!!
  #     t2 <- t2 + np_surv[jj]
  #   }
  # }

  if(!is.null(paras.ini)){
    if(length(paras) != p || length(paras.ini) != p ){
      message("The length of paras.ini is not correct.")
      browser()
      stop("The length of paras.ini is not correct.") 
    }
  }else{
    if(length(paras) != p )
      stop("The length of paras.ini is not correct.")
  }

  #initialisation
  #   paraOpt <- paras
  posfix <- rep(0,length(paras)) # 0 = non fixe 1 = fixe # initialisation
  # constraining of parameters==============
  indexFixe <- indexparaFixeForIden

  if(!is.null(indexparaFixeUser)){
    if(length(indexparaFixeUser) != length(paraFixeUser)){
      stop("The length of paraFixe does not correspond with the length of indexparaFixe")
    }
    indexFixe <- sort(unique(c(indexFixe,indexparaFixeUser)))
  }
  paraFixe <- rep(NA, length(posfix))
  if(!is.null(paraFixeUser)){
    paraFixe[c(indexparaFixeUser)]<- paraFixeUser
  }
  paraFixe[index_paraFixe_mu0_constraint]<- rep(0,K)
  paraFixe[index_paraFixeDconstraint]<- rep(1,K)
  if(sum(!is.na(paraFixe))==0){
    paraFixe = -1
    paraOpt <- paras
  }else{
    paraFixe <- paraFixe[!is.na(paraFixe)]
    posfix[indexFixe] <- 1 # fixation des paras d'indexes dans indexparaFixe
    paras[indexFixe] <- paraFixe
    paraOpt <- paras[-indexFixe]
  }
  
  return(list(para = paras, paraOpt = paraOpt, paraFixe = paraFixe, posfix = posfix, L = L, basehaz = basehaz, knots_surv = knots_surv, np_surv = np_surv, assoc = assoc, truncation = truncation ))
}


#' Initialisation of parameters
#'
#' @param data indicates the data frame containing all the variables for estimating the model
#' @param outcomes names of the outcomes
#' @param mapped.to.LP indicates which outcome measured which latent process, it is a mapping table between
#'  outcomes and latents processes
#' @param fixed_X0.models fixed effects in the submodel for the baseline level of processes
#' @param fixed_DeltaX.models a two-sided linear formula object for specifying the response outcomes (one the left part of ~ symbol) 
#' and the covariates with fixed-effects (on the right part of ~ symbol) 
#' @param randoms_DeltaX.models random effects in the submodel for change over time of latent processes
#' @param randoms_X0.models random effects in the submodel for the baseline level of processes
#' @param nb_RE number of random effects
#' @param mod_trans.model model for elements of the temporal transition matrix, which captures 
#' the temporal influences between latent processes
#' @param subject indicates the name of the covariate representing the grouping structure
#' @param Time indicates the name of the covariate representing the time
#' @param link indicates link used to transform outcome
#' @param knots indicates position of knots used to transform outcomes 
#' @param DeltaT indicates the discretization step
#' @param maxiter maximum iteration
#' @param epsa threshold for the convergence criterion on the parameters, default value is 1.e-4
#' @param epsb threshold for the convergence criterion on the likelihood, default value is 1.e-4
#' @param epsd threshold for the convergence criterion on the derivatives, default value is 1.e-3
#' @param nproc number of processor to be used for running this package
#' @param print.info  to print information during the liklihood optimization, default value is FALSE
#'
#' @return a list

f_paras.ini <- function(data, outcomes, mapped.to.LP, fixed_X0.models, fixed_DeltaX.models, randoms_DeltaX.models, 
                        randoms_X0.models, nb_RE, mod_trans.model, subject, 
                        MCnr = NULL, type_int = NULL,
                        Survdata = NULL, basehaz = NULL,
                        Time, link, knots, zitr = NULL, ide = NULL, DeltaT, maxiter = 25, epsa = .0001, epsb = .0001,
                        epsd = .0001, nproc = 1, print.info = TRUE, 
                        TimeDiscretization = FALSE, fixed.survival.models = NULL, Tentry = NULL, Event = NULL, StatusEvent = NULL,
                        assocT = NULL, truncation = NULL)
{
  cl <- match.call()
  
  K <- length(unlist(outcomes))
  nD <- length(fixed_X0.models)
  paras.ini <-  NULL
  para.fixed_X0 <- NULL
  para.fixed_DeltaX <- NULL
  para.RE <- NULL
  para.trans <- NULL
  para.Sig <- NULL
  ParaTransformY <- NULL
  paras.ini <- list()

  fixed.survival <- NULL
  if(!is.null(fixed.survival.models)){
    fixed.survival <- paste("~", fixed.survival.models[1], sep="")
    for(j in 2:length(fixed.survival.models)){
      fixed.survival <- paste(fixed.survival, fixed.survival.models[j], sep="|")
    }
    fixed.survival<-as.formula(fixed.survival)
  }
  
  for(k in 1 : K){
    data <- data[!is.na(data[,outcomes[k]]),]
    fixed_DeltaX <- as.formula(paste(outcomes[k],"~",fixed_DeltaX.models[mapped.to.LP[k]], sep=""))
    fixed_X0 <- as.formula(paste("~", fixed_X0.models[mapped.to.LP[k]], sep=""))
    n_col_x_k <- ncol(model.matrix(object = fixed_DeltaX,data = data ))
    n_col_x0_k <- ncol(model.matrix(object = fixed_X0,data = data ))
    mod_randoms_DeltaX <- as.formula(paste("~", randoms_DeltaX.models[mapped.to.LP[k]], sep=""))
    mod_trans <- as.formula(paste("~", mod_trans.model, sep=" "))
    indexparaFixeUser <- c(1,(n_col_x0_k+n_col_x_k+1))
    paraFixeUser <- c(0,1)  #intercept latent process = 0
    #link_k <- ifelse(link[k]=="thresholds", "linear",link[k])

    structural.model <- list(fixed.LP0 = fixed_X0,
                             fixed.DeltaLP = fixed_DeltaX,
                             random.DeltaLP = mod_randoms_DeltaX, 
                             trans.matrix = mod_trans, 
                             delta.time = DeltaT,
                             fixed.survival = fixed.survival)
    measurement.model <- list(link.functions = list(links = link[k], knots = knots[k]))
    
    parameters = list(Fixed.para.index = indexparaFixeUser, Fixed.para.values = paraFixeUser)
    
    option = list(nproc = nproc, print.info = print.info, maxiter = maxiter, MCnr = MCnr, type_int = type_int,
                  basehaz = basehaz, assocT = assocT, truncation = truncation)

    mod <- CInLPN2:::CInLPN2(structural.model = structural.model, measurement.model = measurement.model, parameters = parameters,
                               option = option, Time = Time, subject = subject, data = data, links = links[k], 
                               Tentry = Tentry, Event = Event, StatusEvent = StatusEvent, TimeDiscretization = TimeDiscretization)

    L <- ncol(mod$modA_mat) 

    ## compute number of parameter per component
    coefficients <- as.vector(mod$coefficients)
    i1 <- 0
    if(k==1 | (k>1 && (mapped.to.LP[k-1]!= mapped.to.LP[k]))){
      para.fixed_X0 <- c(para.fixed_X0, coefficients[(i1+1):(i1+n_col_x0_k)])
      i1 <- i1+n_col_x0_k
      para.fixed_DeltaX <- c(para.fixed_DeltaX, coefficients[(i1+1):(i1+n_col_x_k)])
      i1 <- i1 + n_col_x_k + mod$nb_paraD
    }
    else{
      i1 <- i1 + n_col_x0_k + mod$nb_paraD
    }
    
    #eta_k =(H(k)+H(k+1))/2
    

    if(k==1){
      para.trans <- c(para.trans, coefficients[(i1+1):(i1+L)])
    }
    if(k!=1 && (mapped.to.LP[k-1]!= mapped.to.LP[k]) ){ # 
      para.trans <- c(para.trans, rep(0,nD*L), coefficients[(i1+1):(i1+L)])
    }
    i1 <- i1+L
    
    para.Sig <- c(para.Sig,  mod$b[(i1+1):(i1+1)])
    i1 <- i1+1
    
    
    ParaTransformY <- c(ParaTransformY,  coefficients[(i1+1):(i1+mod$length_para_trY)])

    i1 <- i1+mod$length_para_trY
  }
  para.RE <- rep(0.1, (nb_RE*(nb_RE+1)/2))

  paras.ini <- c(para.fixed_X0, para.fixed_DeltaX, para.RE, para.trans, para.Sig, ParaTransformY) 

  
  #Survival sub model
  if(!is.null(Survdata)){
    if(basehaz == "Weibull"){
      browser()
      mod_S <- coxph(Surv(Event, StatusEvent) ~ 1, data=lung)
      mod_S <- survreg(Surv(Event, StatusEvent, type='left') ~ 1,
                       data=Survdata, dist='weibull')
      cat("check: type= left is left truncation? add type_surv + X_surv.")
      }else{
        cat("ToDo survreg with splines")
      }

    paras.ini <- c(paras.ini, mod_S$coefficients)
  }
  return(paras.ini)
}

