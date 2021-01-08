#' @name .getweights
#' @noRd
.getweights <- function(W,Data,OE.weights=NULL,Wex.restr=NULL,variable.list=NULL){
  cN<-names(Data) # ord
  nn<-lapply(Data,colnames)
  exo=which(table(unlist(nn))==1) #exo variables are those that are only in one country model
  endo<-which(!table(unlist(nn))==1)
  exo.countries<-names(which(sapply(lapply(nn,function(y) names(exo)%in%y),any))) #gives you a vector of exo countries
  OE.flag <- FALSE
  
  W.sets<-length(W)
  if(is.null(variable.list)){
    variable.list<-list();variable.list$vars<-names(endo)
  }
  
  #------------------------------------- additional checks -----------------------------------------------------#
  if(is.list(W)&&length(W)>1){
    if(is.null(variable.list)){
      stop("You have submitted more than 1 weight matrix but not specified the according variable sets.")
    }
    if(length(W)!=length(variable.list)){
      stop("You have submitted more than 1 weight matrix but not the same number of variable sets.")
    }
    if(!all(names(endo)%in%unlist(variable.list))){
      stop("You have submitted more than 1 weight matrix but some of the variables are not assigned to a weight matrix by the variable list")
    }
  }
  if(!is.null(OE.weights)){
    OE.flag <- TRUE
    OE.sets <- length(OE.weights)
    OE.cN   <- names(OE.weights)
    OE.vars <- lapply(OE.weights,function(l)l$variables)
    OE <- list()
    for(kk in 1:OE.sets){
      OE[[OE.cN[kk]]] <- Data[[OE.cN[kk]]]
      Data[[OE.cN[kk]]] <- NULL
    }
    if(!all(unlist(lapply(OE.vars,function(l)l%in%c(names(exo),names(endo)))))){
      stop("Please specify for the additional entities variables which are also contained in the data. Please respecify.")
    }
    OE.exo <- lapply(OE.weights,function(l)l$exo)
    for(kk in 1:OE.sets){
      if(is.null(OE.exo[[kk]])) OE.exo[[kk]] <- colnames(OE[[kk]])
    }
    if(!any(sapply(1:OE.sets,function(oo)any(OE.exo[[oo]]%in%OE.vars[[oo]])))){
      stop("Please specify the exogenous variables also in the element 'variables'. Please respecify.")
    }
    OE.weights <- lapply(OE.weights,function(l)l$weights)
    if(any(unlist(lapply(OE.weights,function(l)is.null(names(l))||!all(names(l)%in%names(Data)))))){
      stop("Either you have not provided names attached to the weights of other entities or the ones you are provided are not contained in the country data. Please respecify.")
    }
  }
  #------------------------------------- build W matrix -----------------------------------------------------#
  gW<-list()
  # exclusion specification 
  if(!is.null(Wex.restr)){
    er.idx<-sapply(nn,function(x) Wex.restr %in% x)
    er.idx<-which(!er.idx)
  }else{
    er.idx<-c()
  }
  
  #make sure that W and Data names are in the same order
  cnames <- names(Data)
  W      <- lapply(W,function(x) x[cnames,cnames])
  
  xglobal <- c()
  for(jj in 1:length(Data)){
    temp <- Data[[jj]];class(temp) <- "numeric"
    temp <- temp[,colSums(is.na(temp))<nrow(temp)]
    colnames(temp) <- paste(cnames[jj],colnames(temp),sep=".")
    xglobal <- cbind(xglobal,temp)
  }
  
  Max.char<-max(nchar(colnames(xglobal)))
  cnt.char<-max(nchar(cnames))
  # for each country
  for(cc in 1:length(cnames)){
    xglobal <- xglobal[,!duplicated(colnames(xglobal))]
    #creates a dynamic list of variables
    varnames <- substr(colnames(xglobal),cnt.char+2,Max.char); varnames <- varnames[!duplicated(varnames)] 
    
    if(length(er.idx)>0){
      if(cc%in%er.idx){
        varnames <- varnames
      }else{
        varnames <- varnames[-charmatch(Wex.restr,varnames)]
      }
    }
    # names of endo variables
    endnames <- unlist(lapply(strsplit(colnames(xglobal)[grepl(cnames[[cc]],colnames(xglobal))],".",fixed=TRUE),
                              function(l)l[2]))
    Wnew <- matrix(0,length(varnames),ncol(xglobal));colnames(Wnew) <- colnames(xglobal);rownames(Wnew) <- varnames
    
    #----------------------here we specify the part for the weakly exogenous variabes-------------------------------------#
    # loop over all variables
    for (kk in 1:nrow(Wnew)){
      # gives you name of countries with this specific variable
      var.cntry.indic <- substr(colnames(Wnew),4,Max.char)==rownames(Wnew)[kk]
      cntry.indic     <- substr(names(Wnew[kk,var.cntry.indic]),1,2)
      if(length(cntry.indic)>0){
        # this selects the weight matrix according to the variable we want to weight (financial, real, etc)
        wghts <- W[[which(sapply(variable.list, function(x) rownames(Wnew)[kk] %in% x))]][cc,]
        # this gives the variable name (e.g., y)
        Wnew[kk,var.cntry.indic] <- wghts[cntry.indic]
      }else{# in case we have exo variables
        # this includes the exo variables in all country models and the exo countries
        Wnew[kk,var.cntry.indic] <- 1
        # in case we look at an exo country, where the variable is endog. defined, we have to set the exo variable to zero
        if(cnames[[cc]]%in%exo.countries&&paste(cnames[cc],rownames(Wnew)[kk],sep=".")%in%colnames(Wnew)){
          Wnew[kk,var.cntry.indic] <- 0
        }
        
      }
    }
    #----------------------here we specify the part for the endogenous variabes-----------------------------------------#
    endoW <- matrix(0,length(endnames),ncol(Wnew))
    endonr <- xglobal[,substr(colnames(xglobal),1,2)==cnames[cc]];rownames(endoW) <- colnames(endonr)
    colnames(endoW) <- colnames(Wnew)
    namesW <- colnames(endoW)
    namesNr <- colnames(endonr)
    
    for (j in 1:nrow(endoW)){
      for (i in 1:length(namesW)){
        if(namesNr[[j]]==namesW[[i]]){
          endoW[j,i]=1
        }
      }
    }
    WfinNR <- rbind(endoW,Wnew)
    
    WfinNR <- WfinNR[!(rowSums(abs(WfinNR)) == 0),]
    gW[[cc]]<-apply(WfinNR,2,function(x) x/(rowSums(WfinNR)))
  }
  names(gW)<-cnames
  #----------------------here we specify the part for extra weights-----------------------------------------#
  if(OE.flag){
    for(kk in 1:OE.sets){
      temp    <- OE[[kk]];class(temp) <- "numeric"
      xglobal <- cbind(xglobal,temp);OE.x<-ncol(OE[[kk]]); OEnames <- colnames(OE[[kk]])
      colnames(xglobal)[(ncol(xglobal)-OE.x+1):ncol(xglobal)] <- paste(OE.cN[kk],".",OEnames,sep="")
      
      for(i in 1:length(cnames)){
        aux<-gW[[cnames[i]]];ii<-nrow(aux)
        aux<-rbind(cbind(aux,matrix(0,ncol=OE.x,nrow=nrow(aux))),matrix(0,nrow=length(OE.exo[[kk]]),ncol=c(ncol(aux)+OE.x)))
        rownames(aux)[(ii+1):(ii+length(OE.exo[[kk]]))]<-OE.exo[[kk]]
        colnames(aux)[(ncol(aux)-OE.x+1):ncol(aux)]<-paste(OE.cN[kk],OEnames,sep=".")
        if(length(OE.exo[[kk]])>1){
          diag(aux[OE.exo[[kk]],paste(OE.cN[kk],".",OE.exo[[kk]],sep="")])<-1
        }else{
          aux[OE.exo[[kk]],paste(OE.cN[kk],".",OE.exo[[kk]],sep="")]<-1
        }
        gW[[cnames[i]]]<-aux
      }
      # this creates the W matrix for the other entity model
      if(!is.null(OE.vars[[kk]])){
        Wnew <- matrix(0,length(OE.vars[[kk]]),ncol(xglobal))
        colnames(Wnew) <- colnames(xglobal)
        rownames(Wnew) <- c(paste(OE.cN[kk],".",OE.vars[[kk]][!OE.vars[[kk]]%in%names(endo)],sep=""),
                            OE.vars[[kk]][OE.vars[[kk]]%in%names(endo)])
        if(OE.x>1){
          diag(Wnew[paste(OE.cN[kk],".",OEnames,sep=""),paste(OE.cN[kk],".",OEnames,sep="")])<-1
        }else{
          Wnew[paste(OE.cN[kk],".",OEnames,sep=""),paste(OE.cN[kk],".",OEnames,sep="")]<-1
        }
        vars <- OE.vars[[kk]][OE.vars[[kk]]%in%names(endo)]
        for(i in 1:length(vars)){
          Wnew[vars[i],paste(names(OE.weights[[kk]]),".",vars[i],sep="")] <- OE.weights[[kk]]
          Wnew[vars[i],]
        }
        # this creates the part if there are more than one other entities
        if(OE.sets>1 && kk < OE.sets){
          aux <- Wnew
          xx <- lapply(OE[(kk+1):OE.sets],function(l)ncol(l))
          xx <- sum(unlist(xx[!OE.cN%in%OE.cN[(kk+1):OE.sets]]))
          names <- c()
          for(kkk in (kk+1):OE.sets){
            names <- c(names,paste(OE.cN[kkk],".",colnames(OE[[kkk]]),sep=""))
          }
          Wnew <- cbind(aux,matrix(0,ncol=xx,nrow=nrow(aux)))
          colnames(Wnew)[(ncol(Wnew)-xx+1):ncol(Wnew)] <- names
        }
        gW[[(length(gW)+1)]]  <- Wnew
        names(gW)[length(gW)] <- OE.cN[kk]
      }
    }
  }
  #----------------------- return everything to main function ---------------------------------------------#
  gW<-gW[cN]
  return(list(gW=gW,bigx=xglobal,exo=exo,exo.countries=exo.countries,endo=endo))
}

#' @name .get_V
#' @noRd
.get_V <- function(k=k,M=M,Mstar=Mstar,p=p,a_bar_1,a_bar_2,a_bar_3,a_bar_4,sigma_sq,sigma_wex,trend=FALSE){
  V_i <- matrix(0,k,M)
  # endogenous part
  for(i in 1:M){
    for(pp in 1:p){
      for(j in 1:M){
        if(i==j){
          #V_i[j+M*(pp-1),i] <- a_bar_1/(pp^2) ######
          V_i[j+M*(pp-1),i] <- (a_bar_1/pp)^2
        }else{
          #V_i[j+M*(pp-1),i] <- (a_bar_2 * sigma_sq[i])/(pp^2*sigma_sq[j]) #####
          V_i[j+M*(pp-1),i] <- (a_bar_2/pp)^2 * (sigma_sq[i]/sigma_sq[j])
        }
      }
    }
  }
  # exogenous part
  for(i in 1:M){
    for(pp in 0:p){
      for(j in 1:Mstar){
        #V_i[M*p+pp*Mstar+j,i] <- a_bar_4 * sigma_sq[i]/(sigma_wex[j]*(pp+1)) #####
        V_i[M*p+pp*Mstar+j,i] <- (a_bar_4/(pp+1))^2 * (sigma_sq[i]/sigma_wex[j])
      }
    }
  }
  # deterministics
  for(i in 1:M){
    if(trend){
      V_i[(k-1):k,i] <- a_bar_3 * sigma_sq[i]
    }else{
      V_i[k,i] <- a_bar_3 * sigma_sq[i]
    }
  }
  return(V_i)
}

#' @name .bernoulli
#' @importFrom stats runif
#' @noRd
.bernoulli <- function(p){
  u <- runif(1)
  if (u<p){
    x=0
  }else{
    x=1
  }
  return(x)
}

#' @name .atau_post
#' @importFrom stats dgamma dexp
#' @noRd
.atau_post <- function(atau,lambda2,thetas,k,rat=1){
  logpost <- sum(dgamma(thetas,atau,(atau*lambda2/2),log=TRUE))+dexp(atau,rate=rat,log=TRUE)
  return(logpost)
}

#' @name .BVAR_linear_wrapper
#' @noRd
#' @importFrom utils capture.output
.BVAR_linear_wrapper <- function(cc, cN, xglobal, gW, prior, plag, draws, burnin, trend, SV, thin, default_hyperpara, Ex){
  Yraw <- xglobal[,substr(colnames(xglobal),1,2)==cN[cc],drop=FALSE]
  W    <- gW[[cc]]
  Exraw <- NULL
  if(!is.null(Ex)) if(cN[cc]%in%names(Ex)) Exraw <- Ex[[cN[cc]]]
  all  <- t(W%*%t(xglobal))
  Wraw <- all[,(ncol(Yraw)+1):ncol(all),drop=FALSE]
  class(Yraw) <- class(Wraw) <- "numeric"
  prior_in <- ifelse(prior=="MN",1,ifelse(prior=="SSVS",2,3))
  if(default_hyperpara[["a_log"]]){
    default_hyperpara["a_start"] <- 1/log(ncol(Yraw))
  }
  invisible(capture.output(bvar<-BVAR_linear(Y_in=Yraw,W_in=Wraw,p_in=plag,draws_in=draws,burnin_in=burnin,cons_in=TRUE,trend_in=trend,sv_in=SV,
                                             thin_in=thin,prior_in=prior_in,hyperparam_in=default_hyperpara,Ex_in=Exraw), type="message"))
  if(is(bvar,"try-error")){
    bvar<-.BVAR_linear_R(Y_in=Yraw,W_in=Wraw,p_in=plag,draws_in=draws,burnin_in=burnin,cons_in=TRUE,trend_in=trend,sv_in=SV,thin_in=thin,prior_in=prior_in,hyperparam_in=default_hyperpara,Ex_in=Exraw)
  }
  #------------------------------------------------ get data ----------------------------------------#
  Y <- bvar$Y; colnames(Y) <- colnames(Yraw); X <- bvar$X
  M <- ncol(Y); Mstar <- ncol(Wraw); bigT <- nrow(Y); K <- ncol(X)
  if(!is.null(Exraw)) Mex <- ncol(Exraw)
  xnames <- c(paste(rep("Ylag",M),rep(seq(1,plag),each=M),sep=""),rep("Wex",Mstar),
              paste(rep("Wexlag",Mstar),rep(seq(1,plag),each=Mstar),sep=""))
  if(!is.null(Exraw)) xnames <- c(xnames,paste(rep("Tex",Mex)))
  xnames <- c(xnames,"cons")
  if(trend) xnames <- c(xnames,"trend")
  colnames(X) <- xnames
  #-----------------------------------------get containers ------------------------------------------#
  A_store <- bvar$A_store; dimnames(A_store)[[2]] <- colnames(X); dimnames(A_store)[[3]] <- colnames(Y)
  # splitting up stores
  dims          <- dimnames(A_store)[[2]]
  a0store       <- adrop(A_store[,which(dims=="cons"),,drop=FALSE],drop=2)
  a1store <- Exstore <- NULL
  if(trend){
    a1store     <- adrop(A_store[,which(dims=="trend"),,drop=FALSE],drop=2)
  }
  if(!is.null(Exraw)){
    Exstore     <- A_store[,which(dims=="Tex"),,drop=FALSE]
  }
  Lambda0store  <- A_store[,which(dims=="Wex"),,drop=FALSE]
  Lambdastore   <- NULL
  Phistore    <- NULL
  for(jj in 1:plag){
    Lambdastore[[jj]] <- A_store[,which(dims==paste("Wexlag",jj,sep="")),,drop=FALSE]
    Phistore[[jj]]  <- A_store[,which(dims==paste("Ylag",jj,sep="")),,drop=FALSE]
  }
  SIGMA_store <- array(NA, c(draws/thin,bigT,M,M)); dimnames(SIGMA_store) <- list(NULL,NULL,colnames(Y),colnames(Y))
  L_store <- bvar$L_store
  for(irep in 1:(draws/thin)){
    for(tt in 1:bigT){
      if(M>1){
        SIGMA_store[irep,tt,,] <- L_store[irep,,]%*%diag(exp(bvar$Sv_store[irep,tt,]))%*%t(L_store[irep,,])
      }else{
        SIGMA_store[irep,tt,,] <- L_store[irep,,]%*%exp(bvar$Sv_store[irep,tt,])%*%t(L_store[irep,,])
      }
    }
  }
  SIGMAmed_store <- apply(SIGMA_store,c(1,3,4),median)
  theta_store   <- bvar$theta_store; dimnames(theta_store)[[2]] <- colnames(X); dimnames(theta_store)[[3]] <- colnames(Y)
  if(SV){
    vola_store  <- bvar$Sv_store; dimnames(vola_store) <- list(NULL,NULL,colnames(Y))
    pars_store  <- bvar$pars_store
    vola_post   <- apply(vola_store,c(2,3),median)
    pars_post   <- apply(pars_store,c(2,3),median)
  }else{
    vola_store  <- bvar$Sv_store; pars_store <- NULL;
    vola_post   <- apply(vola_store,c(2,3),median); pars_post <- NULL
  }
  res_store     <- bvar$res_store; dimnames(res_store) <- list(NULL,NULL,colnames(Y))
  # MN
  if(prior=="MN"){
    shrink_store  <- bvar$shrink_store; dimnames(shrink_store) <- list(NULL,c("shrink1","shrink2","shrink4"))
    shrink_post   <- apply(shrink_store,2,median)
  }else{
    shrink_store  <- shrink_post <- NULL
  }
  # SSVS
  if(prior=="SSVS"){
    gamma_store <- bvar$gamma_store; dimnames(gamma_store) <- list(NULL,colnames(X),colnames(Y))
    omega_store <- bvar$omega_store; dimnames(omega_store) <- list(NULL,colnames(Y),colnames(Y))
    PIP         <- apply(gamma_store,c(2,3),mean)
    PIP_omega   <- apply(omega_store,c(2,3),mean)
  }else{
    gamma_store <- omega_store <- PIP <- PIP_omega <- NULL
  }
  # NG
  if(prior=="NG"){
    lambda2_store <- bvar$lambda2_store
    tau_store     <- bvar$tau_store
    dimnames(lambda2_store) <- list(NULL,paste("lag",0:plag,sep="_"),c("endogenous","weakly exogenous","covariance"))
    dimnames(lambda2_store) <- list(NULL,paste("lag",0:plag,sep="_"),c("endogenous","weakly exogenous","covariance"))
    lambda2_post  <- apply(lambda2_store,c(2,3),median)
    tau_post      <- apply(tau_store,c(2,3),median)
  }else{
    lambda2_store <- tau_store <- lambda2_post <- tau_post <- NA
  }
  store <- list(A_store=A_store,a0store=a0store,a1store=a1store,Lambda0store=Lambda0store,Lambdastore=Lambdastore,Phistore=Phistore,Exstore=Exstore,SIGMA_store=SIGMA_store,SIGMAmed_store=SIGMAmed_store,L_store=L_store,theta_store=theta_store,vola_store=vola_store,pars_store=pars_store,res_store=res_store,shrink_store=shrink_store,gamma_store=gamma_store,omega_store=omega_store,lambda2_store=lambda2_store,tau_store=tau_store)
  #------------------------------------ compute posteriors -------------------------------------------#
  A_post      <- apply(A_store,c(2,3),median)
  SIGMA_post  <- apply(SIGMA_store,c(2,3,4),median)
  S_post      <- apply(SIGMA_post,c(2,3),mean)
  Sig         <- S_post/(bigT-K)
  theta_post  <- apply(theta_store,c(2,3),median)
  res_post    <- apply(res_store,c(2,3),median)
  # splitting up posteriors
  a0post      <- A_post[which(dims=="cons"),,drop=FALSE]
  a1post <- Expost <- NULL
  if(trend){
    a1post    <- A_post[which(dims=="trend"),,drop=FALSE]
  }
  if(!is.null(Exraw)){
    Expost    <- A_post[which(dims=="Tex"),,drop=FALSE]
  }
  Lambda0post <- A_post[which(dims=="Wex"),,drop=FALSE]
  Lambdapost  <- NULL
  Phipost     <- NULL
  for(jj in 1:plag){
    Lambdapost <- rbind(Lambdapost,A_post[which(dims==paste("Wexlag",jj,sep="")),,drop=FALSE])
    Phipost    <- rbind(Phipost,A_post[which(dims==paste("Ylag",jj,sep="")),,drop=FALSE])
  }
  post <- list(A_post=A_post,a0post=a0post,a1post=a1post,Lambda0post=Lambda0post,Lambdapost=Lambdapost,Phipost=Phipost,Expost=Expost,SIGMA_post=SIGMA_post,S_post=S_post,Sig=Sig,theta_post=theta_post,vola_post=vola_post,pars_post=pars_post,res_post=res_post,shrink_post=shrink_post,PIP=PIP,PIP_omega=PIP_omega,lambda2_post=lambda2_post,tau_post=tau_post)
  return(list(Y=Y,X=X,W=W,store=store,post=post))
}

#' @name .BVAR_linear_R
#' @importFrom stochvol svsample_fast_cpp specify_priors get_default_fast_sv
#' @importFrom MASS ginv mvrnorm
#' @importFrom methods is
#' @importFrom stats rnorm rgamma runif dnorm
#' @export
#' @noRd
.BVAR_linear_R <- function(Y_in,W_in,p_in,draws_in,burnin_in,cons_in,trend_in,sv_in,thin_in,quiet_in,prior_in,hyperparam_in,Ex_in){
  #----------------------------------------INPUTS----------------------------------------------------#
  Yraw  <- Y_in
  p     <- p_in
  Traw  <- nrow(Yraw)
  M     <- ncol(Yraw)
  K     <- M*p
  Ylag  <- .mlag(Yraw,p)
  nameslags <- NULL
  for (ii in 1:p) nameslags <- c(nameslags,rep(paste("Ylag",ii,sep=""),M))
  colnames(Ylag) <- nameslags
  
  exo   <- FALSE; Mstar <- 0; Kstar <- 0; wexnameslags <- NULL; Wraw <- NULL; Wexlag <- NULL
  if(!is.null(W_in)){
    Wraw  <- W_in; Mstar <- ncol(Wraw); Kstar <- Mstar*(p+1)
    exo <- TRUE
    Wexlag <- .mlag(Wraw,p)
    colnames(Wraw) <- rep("Wex",Mstar)
    for (ii in 1:p) wexnameslags <- c(wexnameslags,rep(paste("Wexlag",ii,sep=""),Mstar))
    colnames(Wexlag) <- wexnameslags
  }
  
  texo <- FALSE; Mex <- 0; Exraw <- NULL
  if(!is.null(Ex_in)){
    Exraw <- Ex_in; Mex <- ncol(Exraw)
    texo <- TRUE
    colnames(Exraw) <- rep("Tex",Mex)
  }
  
  X <- cbind(Ylag,Wraw,Wexlag,Exraw)
  X <- X[(p+1):nrow(X),,drop=FALSE]
  Y <- Yraw[(p+1):Traw,,drop=FALSE]
  bigT  <- nrow(X)
  
  cons  <- cons_in
  if(cons){
    X <- cbind(X,1)
    colnames(X)[ncol(X)] <- "cons"
  }
  trend <- trend_in
  if(trend){
    X <- cbind(X,seq(1,bigT))
    colnames(X)[ncol(X)] <- "trend"
  }
  
  k     <- ncol(X)
  n <- k*M
  v <- (M*(M-1))/2
  #---------------------------------------------------------------------------------------------------------
  # HYPERPARAMETERS
  #---------------------------------------------------------------------------------------------------------
  hyperpara <- hyperparam_in
  prior     <- prior_in
  sv        <- sv_in
  prmean    <- hyperpara$prmean
  a_1       <- hyperpara$a_1
  b_1       <- hyperpara$b_1
  crit_eig  <- hyperpara$crit_eig
  Bsigma    <- hyperpara$Bsigma
  a0        <- hyperpara$a0
  b0        <- hyperpara$b0
  bmu       <- hyperpara$bmu
  Bmu       <- hyperpara$Bmu
  # prior == 1: MN
  shrink1   <- hyperpara$shrink1
  shrink2   <- hyperpara$shrink2
  shrink3   <- hyperpara$shrink3
  shrink4   <- hyperpara$shrink4
  # prior == 2: SSVS
  tau00     <- hyperpara$tau0
  tau11     <- hyperpara$tau1
  p_i       <- hyperpara$p_i
  kappa0    <- hyperpara$kappa0
  kappa1    <- hyperpara$kappa1
  q_ij      <- hyperpara$q_ij
  # prior == 3: NG
  d_lambda  <- hyperpara$d_lambda
  e_lambda  <- hyperpara$e_lambda
  a_start   <- hyperpara$a_start
  sample_A  <- hyperpara$sample_A
  #---------------------------------------------------------------------------------------------------------
  # OLS Quantitites
  #---------------------------------------------------------------------------------------------------------
  XtXinv <- try(solve(crossprod(X)),silent=TRUE)
  if(is(XtXinv,"try-error")) XtXinv <- ginv(crossprod(X))
  A_OLS  <- XtXinv%*%(t(X)%*%Y)
  E_OLS  <- Y - X%*%A_OLS
  #a_OLS <- as.vector(A_OLS)
  #SSE  <-  t((Y - X%*%A_OLS))%*%(Y - X%*%A_OLS)
  SIGMA_OLS  <- crossprod(E_OLS)/(bigT-k)
  #IXY  <-   kronecker(diag(M),(t(X)%*%Y))
  #---------------------------------------------------------------------------------------------------------
  # Initial Values
  #---------------------------------------------------------------------------------------------------------
  A_draw <- A_OLS
  SIGMA  <- array(SIGMA_OLS, c(M,M,bigT))
  Em     <- Em_str <- E_OLS
  L_draw <- diag(M)
  #---------------------------------------------------------------------------------------------------------
  # PRIORS
  #---------------------------------------------------------------------------------------------------------
  # Priors on VAR coefs
  #-----------------------------
  # prior mean
  A_prior <- matrix(0,k,M)
  diag(A_prior) <- prmean
  a_prior  <-  as.vector(A_prior)
  # prior variance
  theta <- matrix(10,k,M)
  
  # MN stuff
  accept1 <- 0
  accept2 <- 0
  accept4 <- 0
  scale1  <- .43
  scale2  <- .43
  scale4  <- .43
  sigma_sq  <- matrix(0,M,1) #vector which stores the residual variance
  for (i in 1:M){
    Ylag_i        <- .mlag(Yraw[,i],p)
    Ylag_i        <- Ylag_i[(p+1):nrow(Ylag_i),,drop=FALSE]
    Y_i           <- Yraw[(p+1):nrow(Yraw),i,drop=FALSE]
    Ylag_i        <- cbind(Ylag_i,seq(1,nrow(Y_i)))
    alpha_i       <- solve(crossprod(Ylag_i))%*%crossprod(Ylag_i,Y_i)
    sigma_sq[i,1] <- (1/(nrow(Y_i)-p-1))*t(Y_i-Ylag_i%*%alpha_i)%*%(Y_i-Ylag_i%*%alpha_i)
  }
  sigma_wex <- matrix(0,Mstar,1)
  for (j in 1:Mstar){
    Ywex_i <- .mlag(Wraw[,j],p)
    Ywex_i <- Ywex_i[(p+1):Traw,]
    Yw_i   <- Wraw[(p+1):Traw,j,drop=FALSE]
    Ywex_i <- cbind(Ywex_i,seq(1,nrow(Yw_i)))
    alpha_w <- solve(crossprod(Ywex_i))%*%t(Ywex_i)%*%Yw_i
    sigma_wex[j,1] <- (1/(nrow(Yw_i)-p-1))*t(Yw_i-Ywex_i%*%alpha_w)%*%(Yw_i-Ywex_i%*%alpha_w)
  }
  if(prior==1){
    theta <- .get_V(k=k,M=M,Mstar=Mstar,p=p,a_bar_1=shrink1,a_bar_2=shrink2,a_bar_3=shrink3,
                    a_bar_4=shrink4,sigma_sq=sigma_sq,sigma_wex=sigma_wex,trend=trend)
  }
  
  # SSVS stuff
  gamma  <-  matrix(1,k,M)
  sigma_alpha  <-  sqrt(diag(kronecker(SIGMA_OLS,XtXinv)))
  tau0 <- matrix(NA, k, M); tau1 <- matrix(NA, k, M)
  ii <- 1
  for(mm in 1:M){
    for(kk in 1:k){
      tau0[kk,mm] <- tau00*sigma_alpha[ii]
      tau1[kk,mm] <- tau11*sigma_alpha[ii]
      ii <- ii+1
    }
  }
  # NG stuff
  lambda2_A    <- matrix(0.01,p+1,2)
  A_tau        <- matrix(a_start,p+1,2)
  colnames(A_tau) <- colnames(lambda2_A) <- c("endo","exo")
  rownames(A_tau) <- rownames(lambda2_A) <- paste("lag.",seq(0,p),sep="")
  A_tuning     <- matrix(.43,p+1,2)
  A_accept     <- matrix(0,p+1,2)
  lambda2_A[1,1] <- A_tau[1,1] <- A_tuning[1,1] <- A_accept[1,1] <- NA
  #------------------------------------
  # Priors on coefs in H matrix of VCV
  #------------------------------------
  # prior mean
  l_prior <- matrix(0,M,M)
  
  # prior variance
  L_prior <- matrix(kappa1,M,M)
  L_prior[upper.tri(L_prior)] <- 0; diag(L_prior) <- 0
  
  # SSVS
  omega <- matrix(1,M,M)
  omega[upper.tri(omega)] <- 0; diag(omega) <- 0
  
  # NG
  lambda2_L <- 0.01
  L_tau     <- a_start
  L_accept  <- 0
  L_tuning  <- .43
  #------------------------------------
  # SV quantities
  #------------------------------------
  Sv_draw <- matrix(-3,bigT,M)
  pars_var <- matrix(c(-3,.9,.2,-3),4,M,dimnames=list(c("mu","phi","sigma","latent0"),NULL))
  Sv_priors <- specify_priors(mu=sv_normal(mean=bmu, sd=Bmu), phi=sv_beta(a0,b0), sigma2=sv_gamma(shape=0.5,rate=1/(2*Bsigma)))
  #---------------------------------------------------------------------------------------------------------
  # SAMPLER MISCELLANEOUS
  #---------------------------------------------------------------------------------------------------------
  nsave <- draws_in
  nburn <- burnin_in
  ntot  <- nsave+nburn
  
  # thinning
  thin         <- thin_in
  count        <- 0
  thindraws    <- nsave/thin
  thin.draws   <- seq(nburn+1,ntot,by=thin)
  #---------------------------------------------------------------------------------------------------------
  # STORAGES
  #---------------------------------------------------------------------------------------------------------
  A_store      <- array(NA,c(thindraws,k,M))
  L_store      <- array(NA,c(thindraws,M,M))
  res_store    <- array(NA,c(thindraws,bigT,M))
  # SV
  Sv_store     <- array(NA,c(thindraws,bigT,M))
  pars_store   <- array(NA,c(thindraws,4,M))
  # MN
  shrink_store <- array(NA,c(thindraws,3))
  # SSVS
  gamma_store  <- array(NA,c(thindraws,k,M))
  omega_store  <- array(NA,c(thindraws,M,M))
  # NG
  theta_store  <- array(NA,c(thindraws,k,M))
  lambda2_store<- array(NA,c(thindraws,p+1,3))
  tau_store    <- array(NA,c(thindraws,p+1,3))
  #---------------------------------------------------------------------------------------------------------
  # MCMC LOOP
  #---------------------------------------------------------------------------------------------------------
  for (irep in 1:ntot){
    #----------------------------------------------------------------------------
    # Step 1: Sample coefficients
    for (mm in 1:M){
      if (mm==1){
        Y.i <- Y[,mm]*exp(-0.5*Sv_draw[,mm])
        X.i <- X*exp(-0.5*Sv_draw[,mm])
        
        V_post <- try(chol2inv(chol(crossprod(X.i)+diag(1/theta[,mm]))),silent=TRUE)
        if (is(V_post,"try-error")) V_post <- ginv(crossprod(X.i)+diag(1/theta[,mm]))
        A_post <- V_post%*%(crossprod(X.i,Y.i)+diag(1/theta[,mm])%*%A_prior[,mm])
        
        A.draw.i <- try(A_post+t(chol(V_post))%*%rnorm(ncol(X.i)),silent=TRUE)
        if (is(A.draw.i,"try-error")) A.draw.i <- mvrnorm(1,A_post,V_post)
        A_draw[,mm] <- A.draw.i
        Em[,mm] <-  Em_str[,mm] <- Y[,mm]-X%*%A.draw.i
      }else{
        Y.i <- Y[,mm]*exp(-0.5*Sv_draw[,mm])
        X.i <- cbind(X,Em[,1:(mm-1)])*exp(-0.5*Sv_draw[,mm])
        
        V_post <- try(chol2inv(chol((crossprod(X.i)+diag(1/c(theta[,mm],L_prior[mm,1:(mm-1)]))))),silent=TRUE)
        if (is(V_post,"try-error")) V_post <- ginv((crossprod(X.i)+diag(1/c(theta[,mm],L_prior[mm,1:(mm-1)]))))
        A_post <- V_post%*%(crossprod(X.i,Y.i)+diag(1/c(theta[,mm],L_prior[mm,1:(mm-1)]))%*%c(A_prior[,mm],l_prior[mm,1:(mm-1)]))
        
        A.draw.i <- try(A_post+t(chol(V_post))%*%rnorm(ncol(X.i)),silent=TRUE)
        if (is(A.draw.i,"try-error")) A.draw.i <- mvrnorm(1,A_post,V_post)
        
        A_draw[,mm] <- A.draw.i[1:ncol(X)]
        Em[,mm] <- Y[,mm]-X%*%A.draw.i[1:ncol(X)]
        Em_str[,mm] <- Y[,mm]-X%*%A.draw.i[1:ncol(X)]-Em[,1:(mm-1),drop=FALSE]%*%A.draw.i[(ncol(X)+1):ncol(X.i),drop=FALSE]
        L_draw[mm,1:(mm-1)] <- A.draw.i[(ncol(X)+1):ncol(X.i)]
      }
    }
    rownames(A_draw) <- colnames(X)
    #----------------------------------------------------------------------------
    # Step 2: different shrinkage prior setups
    # MN
    if(prior==1){
      #Step for the first shrinkage parameter (own lags)
      shrink1.prop <- exp(rnorm(1,0,scale1))*shrink1
      theta1.prop   <- .get_V(k=k,M=M,Mstar=Mstar,p=p,a_bar_1=shrink1.prop,a_bar_2=shrink2,a_bar_3=shrink3,
                              a_bar_4=shrink4,sigma_sq=sigma_sq,sigma_wex=sigma_wex)
      post1.prop<-sum(dnorm(as.vector(A_draw),as.vector(A_prior),sqrt(as.vector(theta1.prop)),log=TRUE))+dgamma(shrink1.prop,0.01,0.01,log=TRUE)
      post1.prop<-post1.prop+log(shrink1.prop) # correction term
      post1 <- sum(dnorm(as.vector(A_draw),as.vector(A_prior),sqrt(as.vector(theta)),log=TRUE))+dgamma(shrink1,0.01,0.01,log=TRUE)
      post1 <- post1+log(shrink1) # correction term
      if ((post1.prop-post1)>log(runif(1,0,1))){
        shrink1 <- shrink1.prop
        theta   <- theta1.prop
        accept1 <- accept1+1
      }
      
      #Step for the second shrinkage parameter (cross equation)
      shrink2.prop <- exp(rnorm(1,0,scale2))*shrink2
      theta2.prop   <- .get_V(k=k,M=M,Mstar=Mstar,p=p,a_bar_1=shrink1,a_bar_2=shrink2.prop,a_bar_3=shrink3,
                              a_bar_4=shrink4,sigma_sq=sigma_sq,sigma_wex=sigma_wex)
      post2.prop <- sum(dnorm(as.vector(A_draw),as.vector(A_prior),sqrt(as.vector(theta2.prop)),log=TRUE))+dgamma(shrink2.prop,0.01,0.01,log=TRUE)
      post2.prop <- post2.prop + log(shrink2.prop) # correction term
      post2 <- sum(dnorm(as.vector(A_draw),as.vector(A_prior),sqrt(as.vector(theta)),log=TRUE))+dgamma(shrink2,0.01,0.01,log=TRUE)
      post2 <- post2 + log(shrink2) # correction term
      if ((post2.prop-post2)>log(runif(1,0,1))){
        shrink2 <- shrink2.prop
        theta   <- theta2.prop
        accept2 <- accept2+1
      }
      
      #Step for the final shrinkage parameter (weakly exogenous)
      shrink4.prop <- exp(rnorm(1,0,scale4))*shrink4
      theta4.prop   <- .get_V(k=k,M=M,Mstar=Mstar,p=p,a_bar_1=shrink1,a_bar_2=shrink2,a_bar_3=shrink3,
                              a_bar_4=shrink4.prop,sigma_sq=sigma_sq,sigma_wex=sigma_wex)
      post4.prop <- sum(dnorm(as.vector(A_draw),as.vector(A_prior),sqrt(as.vector(theta4.prop)),log=TRUE))+dgamma(shrink4.prop,0.01,0.01,log=TRUE)
      post4.prop <- post4.prop + log(shrink4.prop)
      post4 <- sum(dnorm(as.vector(A_draw),as.vector(A_prior),sqrt(as.vector(theta)),log=TRUE))+dgamma(shrink4,0.01,0.01,log=TRUE)
      post4 <- post4 + log(shrink4)
      if ((post4.prop-post4)>log(runif(1,0,1))){
        shrink4  <- shrink4.prop
        theta    <- theta4.prop
        accept4  <- accept4+1
      }
      
      if (irep<(0.5*nburn)){
        if ((accept1/irep)<0.15) scale1 <- 0.99*scale1
        if ((accept1/irep)>0.3)  scale1 <- 1.01*scale1
        if ((accept2/irep)<0.15) scale2 <- 0.99*scale2
        if ((accept2/irep)>0.3)  scale2 <- 1.01*scale2
        if ((accept4/irep)<0.15) scale4 <- 0.99*scale4
        if ((accept4/irep)>0.3)  scale4 <- 1.01*scale4
      }
    }
    # SSVS
    if(prior==2){
      for(mm in 1:M){
        for(kk in 1:k){
          u_i1  <-  dnorm(A_draw[kk,mm],A_prior[kk,mm],tau0[kk,mm]) * p_i
          u_i2  <-  dnorm(A_draw[kk,mm],A_prior[kk,mm],tau1[kk,mm]) * (1-p_i)
          gst  <-  u_i1/(u_i1 + u_i2)
          if(gst=="NaN") gst <- 0
          gamma[kk,mm]  <-  .bernoulli(gst)
          gamma[is.na(gamma)] <- 1
          if (gamma[kk,mm] == 0){
            theta[kk,mm]  <-  tau0[kk,mm]^2
          }else if (gamma[kk,mm] == 1){
            theta[kk,mm]  <-  tau1[kk,mm]^2
          }
        }
      }
      for(mm in 2:M){
        for(ii in 1:(mm-1)){
          u_ij1  <-  dnorm(L_draw[mm,ii],l_prior[mm,ii],kappa0) * q_ij
          u_ij2  <-  dnorm(L_draw[mm,ii],l_prior[mm,ii],kappa1) * (1-q_ij)
          ost  <-  u_ij1/(u_ij1 + u_ij2)
          if(is.na(ost)) ost <- 1
          omega[mm,ii] <-  .bernoulli(ost)
          if (is.na(omega[mm,ii])) omega[mm,ii] <- 1
          if(omega[mm,ii]==1){
            L_prior[mm,ii] <- kappa1^2
          }else{
            L_prior[mm,ii] <- kappa0^2
          }
        }
      }
    }
    # NG
    if(prior==3){
      # Normal-Gamma for Covariances
      lambda2_L    <- rgamma(1,d_lambda+L_tau*v,e_lambda+L_tau/2*sum(L_prior[lower.tri(L_prior)]))
      #Step VI: Sample the prior scaling factors for covariances from GIG
      for(mm in 2:M){
        for(ii in 1:(mm-1)){
          L_prior[mm,ii] <- do_rgig1(lambda=L_tau-0.5, chi=(L_draw[mm,ii]-l_prior[mm,ii])^2, psi=L_tau*lambda2_L)
        }
      }
      if(sample_A){
        #Sample L_tau through a simple RWMH step
        L_tau_prop       <- exp(rnorm(1,0,L_tuning))*L_tau
        post_L_tau_prop  <- .atau_post(atau=L_tau_prop, thetas=L_prior[lower.tri(L_prior)], k=v, lambda2=lambda2_L)
        post_L_tau_old   <- .atau_post(atau=L_tau,      thetas=L_prior[lower.tri(L_prior)], k=v, lambda2=lambda2_L)
        post.diff    <- post_L_tau_prop-post_L_tau_old
        post.diff    <- ifelse(is.nan(post.diff),-Inf,post.diff)
        if (post.diff > log(runif(1,0,1))){
          L_tau      <- L_tau_prop
          L_accept   <- L_accept+1
        }
        if (irep<(0.5*nburn)){
          if ((L_accept/irep)>0.3)  L_tuning <- 1.01*L_tuning
          if ((L_accept/irep)<0.15) L_tuning <- 0.99*L_tuning
        }
      }
      # Norml-Gamma for weakly exogenous
      for (ss in 0:p){
        if (ss==0) slct.i <- which(rownames(A_draw)=="Wex") else slct.i <- which(rownames(A_draw)==paste("Wexlag",ss,sep=""))
        A.lag.star  <- A_draw[slct.i,,drop=FALSE]
        A.lag.prior <- A_prior[slct.i,,drop=FALSE]
        theta.lag   <- theta[slct.i,,drop=FALSE]
        
        MMstar <- ncol(A.lag.star)*nrow(A.lag.star)
        if (ss==0){
          lambda2_A[ss+1,2] <- rgamma(1,d_lambda+A_tau[ss+1,2]*MMstar,e_lambda+A_tau[ss+1,2]/2*sum(theta.lag))
        }else{
          lambda2_A[ss+1,2] <- rgamma(1,d_lambda+A_tau[ss+1,2]*MMstar,e_lambda+A_tau[ss+1,2]/2*prod(lambda2_A[1:ss,2])*sum(theta.lag))
        }
        for (jj in 1:Mstar){
          for (ii in 1:M){
            theta.lag[jj,ii] <- do_rgig1(lambda=A_tau[ss+1,2]-0.5,
                                         chi=(A.lag.star[jj,ii]-A.lag.prior[jj,ii])^2,
                                         psi=A_tau[ss+1,2]*prod(lambda2_A[1:(ss+1),2]))
          }
        }
        theta[slct.i,] <- theta.lag
        theta[theta<1e-7] <- 1e-7
        
        if(sample_A){
          #Sample a_tau through a simple RWMH step (on-line tuning of the MH scaling within the first 50% of the burn-in phase)
          A_tau_prop       <- exp(rnorm(1,0,A_tuning[ss+1,2]))*A_tau[ss+1,2]
          post_A_tau_prop  <- .atau_post(atau=A_tau_prop,    thetas=as.vector(theta.lag),lambda2 = prod(lambda2_A[1:(ss+1),2]))
          post_A_tau_old   <- .atau_post(atau=A_tau[ss+1,2], thetas=as.vector(theta.lag),lambda2 = prod(lambda2_A[1:(ss+1),2]))
          post.diff        <- post_A_tau_prop-post_A_tau_old
          post.diff        <- ifelse(is.nan(post.diff),-Inf,post.diff)
          if (post.diff > log(runif(1,0,1))){
            A_tau[ss+1,2]    <- A_tau_prop
            A_accept[ss+1,2] <- A_accept[ss+1,2]+1
          }
          if (irep<(0.5*nburn)){
            if ((A_accept[ss+1,2]/irep)>0.3)  A_tuning[ss+1,2] <- 1.01*A_tuning[ss+1,2]
            if ((A_accept[ss+1,2]/irep)<0.15) A_tuning[ss+1,2] <- 0.99*A_tuning[ss+1,2]
          }
        }
      }
      # Normal-Gamma for endogenous variables
      for (ss in 1:p){
        slct.i    <- which(rownames(A_draw)==paste("Ylag",ss,sep=""))
        A.lag     <- A_draw[slct.i,,drop=FALSE]
        A.prior   <- A_prior[slct.i,,drop=FALSE]
        theta.lag <- theta[slct.i,,drop=FALSE]
        
        M.end <- ncol(A.lag.star)*nrow(A.lag.star)
        if (ss==1){
          lambda2_A[ss+1,1] <- rgamma(1,d_lambda+A_tau[ss+1,1]*M.end,e_lambda+A_tau[ss+1,1]/2*sum(theta.lag))
        }else{
          lambda2_A[ss+1,1] <- rgamma(1,d_lambda+A_tau[ss+1,1]*M.end,e_lambda+A_tau[ss+1,1]/2*prod(lambda2_A[2:(ss+1),1])*sum(theta.lag))
        }
        for (jj in 1:M){
          for (ii in 1:M){
            theta.lag[jj,ii] <- do_rgig1(lambda=A_tau[ss+1,1]-0.5,
                                         chi=(A.lag[jj,ii]-A.prior[jj,ii])^2,
                                         psi=A_tau[ss+1,1]*prod(lambda2_A[2:(ss+1),1]))
          }
        }
        theta[slct.i,] <- theta.lag
        theta[theta<1e-8] <- 1e-8
        #TO BE MODIFIED
        if (sample_A){
          #Sample a_tau through a simple RWMH step (on-line tuning of the MH scaling within the first 50% of the burn-in phase)
          A_tau_prop <- exp(rnorm(1,0,A_tuning[ss+1,1]))*A_tau[ss+1,1]
          post_A_tau_prop <- .atau_post(atau=A_tau_prop,    thetas=as.vector(theta.lag),lambda2=prod(lambda2_A[2:(ss+1),1]))
          post_A_tau_old  <- .atau_post(atau=A_tau[ss+1,1], thetas=as.vector(theta.lag),lambda2=prod(lambda2_A[2:(ss+1),1]))
          post.diff <- post_A_tau_prop-post_A_tau_old
          post.diff <- ifelse(is.nan(post.diff),-Inf,post.diff)
          
          if (post.diff > log(runif(1,0,1))){
            A_tau[ss+1,1] <- A_tau_prop
            A_accept[ss+1,1] <- A_accept[ss+1,1]+1
          }
          if (irep<(0.5*nburn)){
            if ((A_accept[ss+1,1]/irep)>0.3)  A_tuning[ss+1,1] <- 1.01*A_tuning[ss+1,1]
            if ((A_accept[ss+1,1]/irep)<0.15) A_tuning[ss+1,1] <- 0.99*A_tuning[ss+1,1]
          }
        }
      }
    }
    #----------------------------------------------------------------------------
    # Step 3: Sample variances
    if (sv){
      for (mm in 1:M){
        para   <- as.list(pars_var[,mm])
        para$nu = Inf; para$rho=0; para$beta<-0
        svdraw <- svsample_fast_cpp(y=Em_str[,mm], draws=1, burnin=0, designmatrix=matrix(NA_real_), 
                                    priorspec=Sv_priors, thinpara=1, thinlatent=1, keeptime="all", 
                                    startpara=para, startlatent=Sv_draw[,mm], 
                                    keeptau=FALSE, print_settings=list(quiet=TRUE, n_chains=1, chain=1), 
                                    correct_model_misspecification=FALSE, interweave=TRUE, myoffset=0, 
                                    fast_sv=get_default_fast_sv())
        para$mu       <- svdraw$para[1,"mu"]
        para$phi      <- svdraw$para[1,"phi"]
        para$sigma    <- svdraw$para[1,"sigma"]
        para$latent0  <- svdraw$latent0[1,"h_0"]
        pars_var[,mm] <- unlist(para[c("mu","phi","sigma","latent0")])
        Sv_draw[,mm]  <- svdraw$latent[1,]
      }
    }else{
      for (jj in 1:M){
        S_1 <- a_1+bigT/2
        S_2 <- b_1+crossprod(Em_str[,jj])/2
        
        sig_eta <- 1/rgamma(1,S_1,S_2)
        Sv_draw[,jj] <- log(sig_eta)
      }
    }
    #----------------------------------------------------------------------------
    # Step 4: store draws
    if(irep %in% thin.draws){
      count <- count+1
      A_store[count,,]   <- A_draw
      L_store[count,,]   <- L_draw
      res_store[count,,] <- Y-X%*%A_draw
      # SV
      Sv_store[count,,]   <- Sv_draw
      pars_store[count,,] <- pars_var
      # MN
      shrink_store[count,] <- c(shrink1,shrink2,shrink4)
      # SSVS
      gamma_store[count,,] <- gamma
      omega_store[count,,] <- omega
      # NG
      theta_store[count,,]     <- theta
      lambda2_store[count,1,3] <- lambda2_L
      lambda2_store[count,1:(p+1),1:2] <- lambda2_A
      tau_store[count,1,3]             <- L_tau
      tau_store[count,1:(p+1),1:2]     <- A_tau
    }
  }
  #---------------------------------------------------------------------------------------------------------
  # END ESTIMATION
  #---------------------------------------------------------------------------------------------------------
  dimnames(A_store)=list(NULL,colnames(X),colnames(A_OLS))
  
  ret <- list(Y=Y,X=X,A_store=A_store,L_store=L_store,Sv_store=Sv_store,shrink_store=shrink_store,gamma_store=gamma_store,omega_store=omega_store,theta_store=theta_store,lambda2_store=lambda2_store,tau_store=tau_store,pars_store=pars_store,res_store=res_store)
  
  return(ret)
}

#' @name .gvar.stacking
#' @importFrom abind adrop
#' @importFrom Matrix bdiag
#' @importFrom stats median
#' @importFrom utils txtProgressBar setTxtProgressBar
#' @noRd
.gvar.stacking<-function(xglobal,plag,globalpost,draws,thin,trend,eigen=FALSE,trim=NULL){
  # initialize objects here
  bigT <- nrow(xglobal) 
  bigK <- ncol(xglobal)
  cN   <- names(globalpost)
  
  thindraws <- draws/thin
  F.eigen   <- numeric(thindraws)
  trim.info <- "No trimming"
  
  A_large     <- array(NA, dim=c(thindraws,bigK,(bigK*plag+1+ifelse(trend,1,0))))
  S_large     <- array(NA, dim=c(thindraws,bigK,bigK))
  Ginv_large  <- array(NA, dim=c(thindraws,bigK,bigK))
  F_large     <- array(NA, dim=c(thindraws,bigK,bigK,plag))
  dimnames(S_large)[[2]]<-dimnames(S_large)[[3]]<-dimnames(Ginv_large)[[2]]<-dimnames(Ginv_large)[[3]]<-dimnames(A_large)[[2]]<-colnames(xglobal)
  
  pb <- txtProgressBar(min = 0, max = thindraws, style = 3)
  for (irep in 1:thindraws){
    a0     <- NULL
    a1     <- NULL
    G      <- NULL
    x      <- NULL
    S_post <- list()
    
    for (cc in 1:length(cN)){
      VAR <- globalpost[[cc]]
      W   <- VAR$W
      A   <- cbind(diag(ncol(VAR$Y)),-t(adrop(VAR$store$Lambda0store[irep,,,drop=FALSE],drop=1)))
      
      for(pp in 1:plag){
        assign(paste("B",pp,sep=""),cbind(t(adrop(VAR$store$Phistore[[pp]][irep,,,drop=FALSE],drop=1)),
                                    t(adrop(VAR$store$Lambdastore[[pp]][irep,,,drop=FALSE],drop=1))))
        if(cc==1) assign(paste("H",pp,sep=""), get(paste("B",pp,sep=""))%*%W)
        if(cc>1)  assign(paste("H",pp,sep=""), rbind(get(paste("H",pp,sep="")),get(paste("B",pp,sep=""))%*%W))
      }
      G            <- rbind(G,A%*%W)
      a0           <- rbind(a0,t(VAR$store$a0store[irep,,drop=FALSE]))
      if(trend) a1 <- rbind(a1,t(VAR$store$a1store[irep,,drop=FALSE]))
      S_post[[cc]] <- apply(adrop(VAR$store$SIGMA_store[irep,,,,drop=FALSE],drop=1),c(2,3),median)
    }
    G.inv  <- solve(G)
    S_large[irep,,] <- as.matrix(bdiag(S_post))
    b0     <- G.inv%*%a0
    if(trend) b1 <- G.inv%*%a1 else b1 <- NULL
    Ginv_large[irep,,] <- G.inv
    
    ALPHA <- NULL
    for (kk in 1:plag){
      assign(paste("F",kk,sep=""),G.inv%*%get(paste("H",kk,sep="")))
      F_large[irep,,,kk] <- get(paste("F",kk,sep=""))
      ALPHA <- cbind(ALPHA,F_large[irep,,,kk])
    }
    
    ALPHA <- cbind(ALPHA,b0,b1)
    A_large[irep,,]<-ALPHA
    
    if(eigen){
      MM  <- .get_companion(ALPHA,c(ncol(xglobal),ifelse(trend,2,1),plag))$MM
      aux <- suppressWarnings(eigen(MM[1:(bigK*plag),1:(bigK*plag)]))
      F.eigen[irep] <- max(abs(Re(aux$values)))
    }
    # if(stats){
    #   X_large         <- cbind(.mlag(xglobal,plag),1)
    #   if(trend) X_large <- cbind(X_large,seq(1:nrow(X_large)))
    #   Y_large         <- xglobal[(plag+1):nrow(xglobal),]
    #   X_large         <- X_large[(plag+1):nrow(X_large),]
    #   globalLik[irep] <- .globalLik(Y=Y_large,X=X_large,Sig=G.inv%*%S_large[irep,,]%*%t(G.inv),ALPHA=ALPHA,bigT=bigT-plag)
    # }
    setTxtProgressBar(pb, irep)
  }
  
  # kick out in-stable draws
  if(!is.null(trim)){
    if(trim==TRUE) trim <- 1.05
    idx<-which(F.eigen<trim)
    
    F_large     <- F_large[idx,,,,drop=FALSE]
    S_large     <- S_large[idx,,,drop=FALSE]
    Ginv_large  <- Ginv_large[idx,,,drop=FALSE]
    A_large     <- A_large[idx,,,drop=FALSE]
    F.eigen     <- F.eigen[idx]
    
    if(length(idx)<10){
      stop("Less than 10 stable draws have been found. Please re-estimate the model.")
    }
    
    trim.info <- round((length(idx)/thindraws)*100,2)
    trim.info <- paste("Trimming leads to ",length(idx) ," (",trim.info,"%) stable draws out of ",thindraws," total draws.",sep="")
  }
  
  results<-list(S_large=S_large,F_large=F_large,Ginv_large=Ginv_large,A_large=A_large,F.eigen=F.eigen,trim.info=trim.info)
  return(results)
}

#' @name .gvar.stacking.wrapper
#' @importFrom stats median
#' @noRd
.gvar.stacking.wrapper<-function(xglobal,plag,globalpost,draws,thin,trend,eigen,trim,verbose){
  bigT      <- nrow(xglobal)
  bigK      <- ncol(xglobal)
  cN        <- names(globalpost)
  thindraws <- draws/thin
  F_large   <- array(NA, dim=c(thindraws,bigK,bigK,plag))
  trim.info <- "No trimming"
  
  ## call Rcpp
  out <- gvar_stacking(xglobal_in=xglobal, plag_in=plag, globalpost_in=globalpost, draws_in=draws,
                       thin_in=thin, trend_in=trend, eigen_in=TRUE, verbose_in=verbose)
  A_large    <- out$A_large
  for(pp in 1:plag){
    F_large[,,,pp] <- out$F_large[,,((bigK*(pp-1))+1):(bigK*pp),drop=FALSE]
  }
  S_large    <- out$S_large
  Ginv_large <- out$Ginv_large
  F.eigen    <- out$F_eigen
  dimnames(S_large)[[2]]<-dimnames(S_large)[[3]]<-dimnames(Ginv_large)[[2]]<-dimnames(Ginv_large)[[3]]<-dimnames(A_large)[[2]]<-colnames(xglobal)
  names <- c(paste(rep(colnames(xglobal),plag),".",rep(seq(1,plag),each=bigK),sep=""),"cons")
  if(trend) names <- c(names,"trend")
  dimnames(A_large)[[3]]<-names
  
  # kick out in-stable draws
  if(eigen){
    idx<-which(F.eigen<trim)
    
    F_large     <- F_large[idx,,,,drop=FALSE]
    S_large     <- S_large[idx,,,drop=FALSE]
    Ginv_large  <- Ginv_large[idx,,,drop=FALSE]
    A_large     <- A_large[idx,,,drop=FALSE]
    F.eigen     <- F.eigen[idx]
    
    if(length(idx)<10){
      stop("Less than 10 stable draws have been found. Please re-estimate the model.")
    }
    
    trim.info <- round((length(idx)/thindraws)*100,2)
    trim.info <- paste("Trimming leads to ",length(idx) ," (",trim.info,"%) stable draws out of ",thindraws," total draws",sep="")
  }
  
  results<-list(S_large=S_large,F_large=F_large,Ginv_large=Ginv_large,A_large=A_large,F.eigen=F.eigen,trim.info=trim.info)
  return(results)
}

#' @name .get_companion
#' @noRd
.get_companion <- function(Beta_,varndxv){
  nn <- varndxv[[1]] # anzahl variablen
  nd <- varndxv[[2]] # anzahl deterministics
  nl <- varndxv[[3]] # anzahl lags
  
  nkk <- nn*nl+nd
  
  Jm <- matrix(0,nkk,nn)
  Jm[1:nn,1:nn] <- diag(nn)
  if (nd>0){
    MM <- rbind(Beta_,cbind(diag((nl-1)*nn), matrix(0,(nl-1)*nn,nn+nd)),cbind(matrix(0,nd,nn*nl),diag(nd)))
  }else{
    MM <- rbind(Beta_,cbind(diag((nl-1)*nn),matrix(0,(nl-1)*nn,nn)))
  }
  
  return(list(MM=MM,Jm=Jm))
}

#' @name .globalLik
#' @noRd
.globalLik <- function(Y,X,Sig,ALPHA,bigT){
  PLS  <- sum(dmvnrm_arma_fast(Y,X%*%t(ALPHA),Sig,TRUE))
  return(PLS)
}

#' @name .avg.shrink
#' @noRd
.avg.shrink <- function(country.shrink,prior){
  cN <- names(country.shrink)
  N  <- length(cN)
  colNames <- lapply(country.shrink,colnames)
  varNames <- lapply(country.shrink,rownames)
  for(cc in 1:N) {
    colNames[[cc]] <- gsub(paste(cN[cc],".",sep=""),"",colNames[[cc]])
    varNames[[cc]] <- gsub(paste(cN[cc],".",sep=""),"",varNames[[cc]])
  }
  colNames <- unique(unlist(colNames))
  varNames <- unique(unlist(varNames))
  shrink <- array(NA,dim=c(length(varNames),length(colNames),N))
  dimnames(shrink)[[1]] <- varNames; dimnames(shrink)[[2]] <- colNames; dimnames(shrink)[[3]] <- cN
  for(cc in 1:N){
    aux <- country.shrink[[cc]];rownames(aux)<-gsub(paste(cN[cc],".",sep=""),"",rownames(aux));colnames(aux)<-gsub(paste(cN[cc],".",sep=""),"",colnames(aux))
    for(z in 1:ncol(aux)){
      shrink[rownames(aux),colnames(aux)[z],cc] <- aux[,z]
    }
  }
  avg.shrink <- apply(shrink,c(1,2),function(x) mean(x,na.rm=TRUE))
  idx1       <- apply(shrink,3,function(x) which(colSums(x,na.rm=TRUE)==0))
  idx2       <- apply(shrink,3,function(x) which(rowSums(x,na.rm=TRUE)==0))
  shrink2    <- list()
  for(cc in 1:N){
    shrink2[[cc]] <- shrink[,,cc]
    idx1.cc       <- ifelse(length(idx1)>0,idx1[[cc]],integer(0))
    idx2.cc       <- ifelse(length(idx2)>0,idx2[[cc]],integer(0))
    if(length(idx1.cc)>0){
      shrink2[[cc]] <- shrink2[[cc]][,-idx1.cc,drop=FALSE]
    }
    if(length(idx2.cc)>0){
      shrink2[[cc]] <- shrink2[[cc]][-idx2.cc,,drop=FALSE]
    }
  }
  names(shrink2) <- cN
  return(list(PIP.cc=shrink2,PIP.avg=avg.shrink))
  
}

#' @name .construct.arglist
#' @noRd
.construct.arglist = function (funobj, envir = NULL){
  namedlist = formals(funobj)
  argnames = names(namedlist)
  if (!is.environment(envir))
    envir = sys.frame(-1)
  for (argn in 1:length(namedlist)) {
    testval = as.logical(try(exists(argnames[argn], envir = envir),
                             silent = TRUE))
    if (is.na(testval))
      testval = FALSE
    if (testval) {
      testout = try(get(argnames[argn], envir = envir),silent = TRUE)
      if (is.null(testout)) {
        namedlist[[argn]] = "list(NULL)blabla"
      } else {
        namedlist[[argn]] = testout
      }
    }
  }
  namedlist = lapply(namedlist,function(x) if (any(x=="list(NULL)blabla")) NULL else x)
  lapply(namedlist, function(l) if(any(l=="list(NULL)blabla")){NULL}else{l})
  return(namedlist)
}

#' @name .mlag
#' @noRd
.mlag <- function(X,lag){
  p <- lag
  X <- as.matrix(X)
  Traw <- nrow(X)
  N <- ncol(X)
  Xlag <- matrix(0,Traw,p*N)
  for (ii in 1:p){
    Xlag[(p+1):Traw,(N*(ii-1)+1):(N*ii)] <- X[(p+1-ii):(Traw-ii),(1:N)]
  }
  colnames(Xlag) <- paste(colnames(X),".lag",rep(seq(p),each=N),sep="")
  return(Xlag)  
}

#' @name .timelabel
#' @noRd
.timelabel <- function(time){
  time   <- as.character(time)
  years  <- regmatches(time,regexpr("^[0-9]{4}",time))
  freq   <- sum(years%in%unique(years)[2])
  xlabel <- gsub("-[0-9]{2}$","",time)
  if(freq==12){
    xlabel<-gsub("-01","-Jan",xlabel,fixed=TRUE)
    xlabel<-gsub("-02","-Feb",xlabel,fixed=TRUE)
    xlabel<-gsub("-03","-Mar",xlabel,fixed=TRUE)
    xlabel<-gsub("-04","-Apr",xlabel,fixed=TRUE)
    xlabel<-gsub("-05","-May",xlabel,fixed=TRUE)
    xlabel<-gsub("-06","-Jun",xlabel,fixed=TRUE)
    xlabel<-gsub("-07","-Jul",xlabel,fixed=TRUE)
    xlabel<-gsub("-08","-Aug",xlabel,fixed=TRUE)
    xlabel<-gsub("-09","-Sep",xlabel,fixed=TRUE)
    xlabel<-gsub("-10","-Oct",xlabel,fixed=TRUE)
    xlabel<-gsub("-11","-Nov",xlabel,fixed=TRUE)
    xlabel<-gsub("-12","-Dec",xlabel,fixed=TRUE)
  }
  if(freq==4){
    xlabel<-gsub("-01"," Q1",xlabel,fixed=TRUE)
    xlabel<-gsub("-02"," Q1",xlabel,fixed=TRUE)
    xlabel<-gsub("-03"," Q1",xlabel,fixed=TRUE)
    xlabel<-gsub("-04"," Q2",xlabel,fixed=TRUE)
    xlabel<-gsub("-05"," Q2",xlabel,fixed=TRUE)
    xlabel<-gsub("-06"," Q2",xlabel,fixed=TRUE)
    xlabel<-gsub("-07"," Q3",xlabel,fixed=TRUE)
    xlabel<-gsub("-08"," Q3",xlabel,fixed=TRUE)
    xlabel<-gsub("-09"," Q3",xlabel,fixed=TRUE)
    xlabel<-gsub("-10"," Q4",xlabel,fixed=TRUE)
    xlabel<-gsub("-11"," Q4",xlabel,fixed=TRUE)
    xlabel<-gsub("-12"," Q4",xlabel,fixed=TRUE)
  }
  return(xlabel)
}

#' @name .irf.sign.zero
#' @importFrom abind abind
#' @importFrom MASS Null
#' @importFrom stats rnorm
#' @noRd
.irf.sign.zero <- function(x,plag,nhor,Ginv,Fmat,Smat,shock,sign.constr,Global,MaxTries,shock.nr,...){
  bigT     <- nrow(x)
  bigK     <- ncol(x)
  varNames <- colnames(x)
  cN       <- unique(sapply(strsplit(varNames,".",fixed=TRUE),function(x)x[1]))
  vars     <- unique(sapply(strsplit(varNames,".",fixed=TRUE),function(x)x[2]))
  N        <- length(cN)
  # shock details
  shock.cN  <- unlist(lapply(sign.constr,function(x)unlist(strsplit(x$shock,".",fixed=TRUE))[1]))
  restr.cN  <- unique(unlist(lapply(sign.constr,function(y)lapply(strsplit(y$restrictions,".",fixed=TRUE),function(x)x[1]))))
  
  P0G <- diag(bigK); colnames(P0G) <- rownames(P0G) <- varNames
  #calculates cholesky factor for blocks of gcov where shock occurrs
  pure.restr.cN <- restr.cN[!(restr.cN%in%shock.cN)]
  P0l<-list();ii<-1
  for(cc in 1:N){
    idx <- which(grepl(cN[cc],varNames))
    if(cN[cc] %in% shock.cN){
      P0             <- try(t(chol(Smat[idx,idx,drop=FALSE])),silent=TRUE)
      if(is(P0,"try-error")) suppressWarnings(P0 <- t(chol(Smat[idx,idx,drop=FALSE],pivot=TRUE)))
      P0l[[ii]]      <- P0
      names(P0l)[ii] <- cN[cc]
      P0G[idx,idx]   <- P0
      Smat[idx,idx]  <- diag(ncol(P0)) 
      ii<-ii+1
    }else if(cN[cc] %in% pure.restr.cN){
      P0l[[ii]]      <- Smat[idx,idx,drop=FALSE]
      names(P0l)[ii] <- cN[cc]
      ii<-ii+1
    }
  }
  ccode <- names(P0l) # now it has right ordering
  
  # create dynamic multiplier
  PHIx <- array(0,c(bigK,bigK,plag+nhor+1)); dimnames(PHIx)[[1]] <- dimnames(PHIx)[[2]] <- varNames
  PHIx[,,plag+1]  <-  diag(bigK)
  for (ihor in (plag+2):(plag+nhor+1)){
    acc = matrix(0,bigK,bigK)
    for (pp in 1:plag){
      acc  <-  acc + Fmat[,,pp]%*%PHIx[,,ihor-pp]
    }
    PHIx[,,ihor]  <-  acc
  }
  PHI  <-  PHIx[,,(plag+1):(plag+nhor+1)]
  #-----------------------------------------------------------------------------
  k              <- sum(unlist(lapply(P0l,ncol)))
  sign.horizon   <- unique(unlist(lapply(sign.constr, function(l) l$rest.horz)))-1 # zero impact is coded as 1
  sign.horizon   <- sort(sign.horizon, decreasing=FALSE)
  sign.shockvars <- unlist(lapply(sign.constr, function(l) l$shock))
  varnames       <- as.character(unlist(lapply(P0l,function(x) dimnames(x)[[1]])))
  H.restr        <- length(sign.horizon)
  N.restr        <- k*H.restr
  S.cube         <- array(NA, c(N.restr, N.restr, k)) # sign restrictions
  Z.cube         <- array(NA, c(N.restr, N.restr, k)) # zero restrictions
  dimnames(S.cube)[[1]] <- dimnames(Z.cube)[[1]] <-
    dimnames(S.cube)[[2]] <- dimnames(Z.cube)[[2]] <- paste(rep(varnames,H.restr),".",
                                                            rep(sign.horizon,each=k),sep="")
  dimnames(S.cube)[[3]] <- dimnames(Z.cube)[[3]] <- varnames
  for(vv in 1:length(varnames)){
    S.temp <- matrix(0, N.restr, N.restr)
    Z.temp <- matrix(0, N.restr, N.restr)
    if(varnames[vv]%in%sign.shockvars){
      sign.shock <- sign.constr[[which(sign.shockvars == varnames[vv])]]$shock
      sign.restr <- sign.constr[[which(sign.shockvars == varnames[vv])]]$restrictions
      sign.signs <- sign.constr[[which(sign.shockvars == varnames[vv])]]$sign
      sign.horiz <- sign.constr[[which(sign.shockvars == varnames[vv])]]$rest.horz-1
      
      sign.restr <- c(sign.shock,sign.restr) ## append positive shock on shock variable
      
      s.point <- which(sign.signs=="<"|sign.signs==">")
      z.point <- seq(1,length(sign.signs))[-s.point]
      
      if(length(s.point)>0){
        for(ss in 1:length(s.point)){
          grp <- which(sign.horiz[s.point[ss]] == sign.horizon)
          col <- seq(which(sign.restr[s.point[ss]]==varnames),k*grp,by=k)
          for(ii in 1:length(col)){
            S.temp[col[ii],col[ii]] <- ifelse(sign.signs[s.point[ss]]=="<",-1,1)
          }
        }
      }
      if(length(z.point)>0){
        for(zz in 1:length(z.point)){
          if(sign.signs[z.point[zz]]=="0"){
            grp <- which(sign.horiz[z.point[zz]] == sign.horizon)
            row <- (grp-1)*k+which(sign.restr[z.point[zz]]==varnames)
            Z.temp[row,row] <- 1
          }else{ # take row from above
            grp <- which(sign.horiz[z.point[zz]] == sign.horizon)
            col <- (grp-1)*k+which(sign.restr[z.point[zz]]==varnames)
            Z.temp[row,col] <- as.numeric(sign.signs[z.point[zz]])
          }
        }
      }
    }
    S.cube[,,vv] <- S.temp
    Z.cube[,,vv] <- Z.temp
  }
  
  no.zero.restr <- ifelse(base::sum(abs(Z.cube))>0,FALSE,TRUE)
  shock.order   <- rep(NA, k)
  search.Znum   <- apply(Z.cube, 3, function(x) base::sum(abs(x)))
  search.Snum   <- apply(S.cube, 3, function(x) base::sum(abs(x)))
  for(mm in 1:k){
    if(!no.zero.restr){
      max.Z <- which(search.Znum==max(search.Znum))
      if(length(max.Z)==1){
        shock.order[mm] <- max.Z
      }else{
        shock.order[mm] <- sample(max.Z,1)
      }
    } else {
      shock.order[mm] <- mm
    }
    search.Znum[shock.order[mm]] <- -1
    search.Snum[shock.order[mm]] <- -1
  }
  shock.order <- varnames[shock.order]
  
  get.cols  <- matrix(0,bigK,k); colnames(get.cols)<- varnames; rownames(get.cols)<- colnames(x)
  vars_list <- lapply(P0l,function(x) dimnames(x)[[1]])
  vars_len  <- lapply(vars_list, length)
  start_col <- 1
  for(kk in 1:length(vars_list)) {
    position  <- which(colnames(x) %in% vars_list[[kk]])
    end_col   <- start_col+length(vars_list[[kk]])-1
    get.cols[position,start_col:end_col] <- diag(length(vars_list[[kk]]))
    start_col <- start_col+length(vars_list[[kk]])
  }
  irf.restr         <- matrix(NA, N.restr, k)
  invGSigma_u       <- Ginv%*%P0G%*%Smat
  for(hh in 1:H.restr){
    # ARRW: Definition 1
    if(sign.horizon[hh]!=Inf) irf.hh<-t(get.cols)%*%PHI[,,sign.horizon[hh]+1]%*%invGSigma_u%*%get.cols
    # ARRW: Definition 2
    #if(sign.horizon[hh]==Inf) irf.hh <- solve(A0-A0%*%Cm[1:M,]%*%do.call("rbind",rep(list(diag(M)),p)))
    irf.restr[((hh-1)*k+1):(k*hh),1:k] <- irf.hh
  }
  colnames(irf.restr) <- varnames
  rownames(irf.restr) <- paste(rep(varnames,H.restr),".",
                              rep(sign.horizon,each=k),sep="")
  
  Z.cube <- Z.cube[,,shock.order]
  
  # draw rotation matrix here
  icounter <- 0
  condall <- 0
  max.counter <- MaxTries
  impresp<-Q_bar<-NA
  while(condall == 0 && icounter < max.counter){
    signCheck <- matrix(NA, k, 1)
    randMat <- matrix(rnorm(k*k,0,1),k,k)
    Q <- matrix(0, k, k)
    if(no.zero.restr){
      Q <- qr(randMat)
      Q <- qr.Q(Q)
    }else{
      for(mm in 1:k){
        Z.temp <- Z.cube[,,mm]
        Z.temp <- Z.temp[rowSums(abs(Z.temp))!=0,,drop=F]
        if(nrow(Z.temp)==0){
          Z.temp <- matrix(0, 1, N.restr)
        }
        if(all(Z.temp==0) && mm>1){
          R <- c()
        }else{
          R <- Z.temp%*%irf.restr
        }
        if(mm > 1){R <- rbind(R, t(Q[,(1:(mm-1)), drop=FALSE]))}
        
        NU  <- Null(t(R))
        x_j <- randMat[,mm,drop =FALSE]
        
        q_j <- NU%*%(t(NU)%*%x_j/sqrt(as.numeric(crossprod(t(NU)%*%x_j))))
        Q[,mm] <- q_j
      }
    }
    
    colnames(Q) <- shock.order; rownames(Q) <- varnames
    Q <- Q[,varnames]
    Q_bar <- Q%*%diag(((diag(Q)>0)-(diag(Q)<0)))
    rotS <- diag(ncol(x));colnames(rotS) <- rownames(rotS) <- substr(colnames(x),1,2)
    R_bar <- matrix(0,k,k)
    startkkk<-1
    for(cc in 1:length(ccode)){
      endkkk <- startkkk+length(vars_list[[cc]])-1
      rotS[rownames(rotS)==ccode[cc],rownames(rotS)==ccode[cc]] <- Q_bar[startkkk:endkkk,startkkk:endkkk]
      R_bar[startkkk:endkkk,startkkk:endkkk] <- Q_bar[startkkk:endkkk,startkkk:endkkk]
      startkkk <- startkkk+length(vars_list[[cc]])
    }
    
    irf.check <- irf.restr%*%R_bar
    colnames(irf.check) <- varnames
    rownames(irf.check) <- paste(rep(varnames,H.restr),".",rep(sign.horizon,each=k),sep="")
    
    for(ss in 1:k){
      STemp <- S.cube[,,ss]
      if(Global){
        var_global <- strsplit(dimnames(S.cube)[[3]][ss],".",fixed=TRUE)[[1]][2]
        IrfCheckTemp <- irf.check[,grep(paste(var_global,"$",sep=""),colnames(irf.check)),drop=FALSE]
        IrfCheckTemp <- matrix(rowSums(IrfCheckTemp),nrow=N.restr,1)
      }else{
        IrfCheckTemp <- irf.check[,ss,drop = FALSE]
      }
      signCheckVec <- matrix(NA, N.restr, 1)
      rownames(signCheckVec) <- paste(rep(varnames,H.restr),".",rep(sign.horizon,each=k),sep="")
      for(kk in 1:N.restr){
        STempRow <- STemp[kk,]
        emptyCheck <- sum(STempRow)
        if(emptyCheck == 0){
          signCheckVec[kk,1] <- 1;
        }else{
          signCheckVec[kk,1] <- as.numeric(STempRow%*%IrfCheckTemp)
        }
      }
      if(sum(abs(STemp))>0){
        shocknr <- which(unlist(lapply(sign.constr,function(l){
          l$shock==dimnames(S.cube)[3][[1]][ss]})))
        constraints <- sign.constr[[shocknr]]$constr
        if(any(constraints!=1)){
          restrictions <- c(sign.constr[[shocknr]]$shock,sign.constr[[shocknr]]$restrictions)
          cntr <- unlist(lapply(strsplit(restrictions,".",fixed=TRUE),function(l) l[1]))
          vars_constr <- unlist(lapply(strsplit(restrictions,".",fixed=TRUE),function(l) l[2]))
          vars_unique <- unique(vars_constr)
          for(kk in 1:length(vars_unique)){
            pointer <- vars_constr%in%vars_unique[kk]
            if(length(pointer)!=1 & any(constraints[pointer]<1)){
              constrainttt  <- unique(constraints[vars_constr%in%vars_unique[kk]])
              cntr_unique <- unique(cntr[pointer])
              restrictionss <- paste(restrictions[pointer],".",
                                    rep(sign.horizon,length(cntr_unique)),sep="")
              v    <- signCheckVec[rownames(signCheckVec)%in%restrictionss,]
              perz <- sum(sign(v)==1)/length(v)
              if(perz>=constrainttt){
                signCheckVec[rownames(signCheckVec)%in%restrictionss,] <- 1
              }
            }
          }
        }
      }
      signCheck[ss,] <- prod((signCheckVec > 0)*(signCheckVec > 0))
    }
    condall <- prod(signCheck)
    icounter <- icounter + 1
  }
  
  st_impulses <- array(NA,c(bigK,bigK,nhor+1));dimnames(st_impulses)[[1]] <- dimnames(st_impulses)[[2]] <- colnames(x)
  Cmhat <- invGSigma_u%*%rotS
  for (ii in 1:(nhor+1)){
    st_impulses[,,ii]  <-   as.matrix((PHI[,,ii]%*%Cmhat))
  }
  
  shock_slct <- sapply(sign.constr, function(l) l$shock)
  if(Global){ # GLOBAL SHOCK
    if(shock.nr>1){
      st_impulses2<-NULL
      for(i in 1:shock.nr){
        st_impulses2<-abind(st_impulses2,apply(st_impulses[,shock_slct[[i]],],c(1,3),sum),along=3)
      }
      #st_impulses <- aperm(st_impulses2,c(1,3,2))
      st_impulses<-st_impulses2
    }else{
      st_impulses <- apply(st_impulses[,shock_slct,],c(1,3),sum)
    }
  }else{
    st_impulses<-st_impulses[,shock_slct,,drop=FALSE]
    st_impulses<-aperm(st_impulses,c(1,3,2)) # re-indicize
  }
  
  if(icounter==MaxTries){
    st_impulses <- R_bar <- rotS <- NA
  }
  # end rotation matrix loop ----------------------------------------------------------------------------
  return(list(impl=st_impulses,rot=rotS,icounter=icounter,gcov=NA))
}

#' @name .irf.chol
#' @noRd
.irf.chol <- function(x,plag,nhor,Ginv,Fmat,Smat,shock,...){
  bigT     <- nrow(x)
  bigK     <- ncol(x)
  varNames <- colnames(x)
  cN       <- unique(sapply(strsplit(varNames,".",fixed=TRUE),function(x)x[1]))
  vars     <- unique(sapply(strsplit(varNames,".",fixed=TRUE),function(x)x[2]))
  N        <- length(cN)
  # shock details
  shock.var <- shock[["var"]]
  shock.cN  <- shock[["cN"]]
  
  P0G <- diag(bigK); colnames(P0G) <- rownames(P0G) <- varNames
  for(cc in 1:N){
    if(cN[cc] %in% shock.cN){
      idx           <- grep(cN[cc],varNames)
      P0G[idx,idx]  <- t(chol(Smat[idx,idx,drop=FALSE])) # calculate local cholesky factor of gcov
      Smat[idx,idx] <- diag(length(idx)) #set vcv matrix to identity for coutnry where shock occurs
    }
  }
  
  # create dynamic multiplier
  PHIx <- array(0,c(bigK,bigK,plag+nhor+1)); dimnames(PHIx)[[1]] <- dimnames(PHIx)[[2]] <- varNames
  PHIx[,,plag+1]  <-  diag(bigK)
  for (ihor in (plag+2):(plag+nhor+1)){
    acc = matrix(0,bigK,bigK)
    for (pp in 1:plag){
      acc  <-  acc + Fmat[,,pp]%*%PHIx[,,ihor-pp]
    }
    PHIx[,,ihor]  <-  acc
  }
  PHI  <-  PHIx[,,(plag+1):(plag+nhor+1)]
  
  # create selection matrix
  eslct <- matrix(0,bigK,1);rownames(eslct) <- varNames
  select_shocks <- paste(shock.cN,shock.var,sep=".")
  eslct[select_shocks,] <- 1
  
  # compute shock
  invGSigma_u  <-  Ginv%*%P0G%*%Smat
  cons<-as.numeric(1/sqrt(t(eslct)%*%Smat%*%eslct))
  
  # computing impulse response function
  irfa  <- array(0,c(bigK,nhor+1)); dimnames(irfa)[[1]]  <- varNames
  for (ihor in 1:(nhor+1)){
    irfa[,ihor]  <-   as.matrix((PHI[,,ihor]%*%invGSigma_u%*%eslct))*cons
  }
  
  out <- list(impl=irfa,rot=NULL)
  return(out)
}

#' @name .irf.girf
#' @noRd
.irf.girf <- function(x,plag,nhor,Ginv,Fmat,Smat,shock, ...){
  bigT     <- nrow(x)
  bigK     <- ncol(x)
  varNames <- colnames(x)
  cN       <- unique(sapply(strsplit(varNames,".",fixed=TRUE),function(x)x[1]))
  vars     <- unique(sapply(strsplit(varNames,".",fixed=TRUE),function(x)x[2]))
  N        <- length(cN)
  # shock details
  shock.var <- shock$var
  shock.cN  <- shock$cN
  
  # create dynamic multiplier
  PHIx <- array(0,c(bigK,bigK,plag+nhor+1)); dimnames(PHIx)[[1]] <- dimnames(PHIx)[[2]] <- varNames
  PHIx[,,plag+1]  <-  diag(bigK)
  for (ihor in (plag+2):(plag+nhor+1)){
    acc = matrix(0,bigK,bigK)
    for (pp in 1:plag){
      acc  <-  acc + Fmat[,,pp]%*%PHIx[,,ihor-pp]
    }
    PHIx[,,ihor]  <-  acc
  }
  PHI  <-  PHIx[,,(plag+1):(plag+nhor+1)]
  
  # create selection matrix
  eslct <- matrix(0,bigK,1);rownames(eslct) <- varNames
  select_shocks <- paste(shock.cN,shock.var,sep=".")
  eslct[select_shocks,] <- 1
  
  invGSigma_u  <-  Ginv%*%Smat
  cons<-as.numeric(1/sqrt(t(eslct)%*%Smat%*%eslct))
  # computing impulse response function
  irfa  <- array(0,c(bigK,nhor+1)); dimnames(irfa)[[1]]  <- varNames
  for (ihor in 1:(nhor+1)){
    irfa[,ihor]  <-   as.matrix((PHI[,,ihor]%*%invGSigma_u%*%eslct))*cons
  }
  
  return(list(impl=irfa,rot=NULL))
}

#' @name .irf.gir.sims
#' @noRd
.irf.girf.sims <- function(invG,lF,gcov,x,horizon=40,...){
  cN <- unique(substr(colnames(x),1,2))
  N <- length(cN)
  Cm <- as.matrix(gcov)
  K <- ncol(x)
  p<-dim(lF)[[3]] # number of lags
  
  invGSigma_u  <-  invG%*%Cm;rownames(invGSigma_u)<-colnames(invGSigma_u)<-rownames(x)
  
  impls.girf<- .impulsdtrf(lF,invGSigma_u,horizon)
  
  cons.girf<-(1/diag(invGSigma_u))*sqrt(diag(invGSigma_u))
  # irfa is a K responses times K shocks times horizon array
  irfa<-impls.girf*array(matrix(cons.girf,K,K,byrow=TRUE),dim=c(K,K,(horizon)))
  
  return(list(impl=irfa,gcov=Cm))
}

#' @name .mk_fevd.sims
#' @noRd
.mk_fevd.sims <- function(irfa){
  ny <- dim(irfa)[[1]];nH <- dim(irfa)[[3]]
  
  fevda <- apply(irfa*irfa,c(1,2),cumsum);
  fevda <- aperm(fevda,c(2,3,1))
  accm <- matrix(0,ny,ny)
  for (ih in 1:nH){
    accm <- accm+irfa[,,ih]%*%t(irfa[,,ih])
    denm <- matrix((diag(accm)),ny,ny)
    fevda[,,ih]=fevda[,,ih]/denm
  }
  return(fevda)
}

#' @name .impulsdtrf
#' @noRd
.impulsdtrf <- function(B,smat,nstep)
  ### By:             As emerges from rfvar, neqn x nvar x lags array of rf VAR coefficients.
  ### smat:           nshock x nvar matrix of initial shock vectors.  To produce "orthogonalized
  ###                 impulse responses" it should have the property that crossprod(t(smat))=sigma,
  ###                 where sigma is the Var(u(t)) matrix and u(t) is the rf residual vector.  One
  ###                 way to get such a smat is to set smat=t(chol(sigma)).  To get the smat
  ###                 corresponding to a different ordering, use
  ###                 smat = t(chol(P %*% Sigma %*% t(P)) %*% P), where P is a permutation matrix.
  ###                 To get impulse responses for a structural VAR in the form A(L)y=eps, with
  ###                 Var(eps)=I, use B(L)=-A_0^(-1)A_+(L) (where A_+ is the coefficients on strictly
  ###                 positive powers of L in A), smat=A_0^(-1).
  ###                 In general, though, it is not required that smat be invertible.
### response:       nvar x nshocks x nstep array of impulse responses.
###
### Code written by Christopher Sims,mat based on 6/03 matlab code.  This version 3/27/04.
### Added dimension labeling, 8/02/04.
{
  
  neq <- dim(B)[1]
  nvar <- dim(B)[2]
  lags <- dim(B)[3]
  dimnB <- dimnames(B)
  if(dim(smat)[2] != dim(B)[2]) stop("B and smat conflict on # of variables")
  response <- array(0,dim=c(neq,nvar,nstep+lags-1));
  response[ , , lags] <- smat
  response <- aperm(response, c(1,3,2))
  irhs <- 1:(lags*nvar)
  ilhs <- lags * nvar + (1:nvar)
  response <- matrix(response, ncol=neq)
  B <- B[, , seq(from=lags, to=1, by=-1)]  #reverse time index to allow matrix mult instead of loop
  B <- matrix(B,nrow=nvar)
  for (it in 1:(nstep-1)) {
    response[ilhs, ] <- B %*% response[irhs, ]
    irhs <- irhs + nvar
    ilhs <- ilhs + nvar
  }
  dim(response) <- c(nvar, nstep + lags - 1, nvar)
  #drop the zero initial conditions; array in usual format
  if(lags>1){
    response<-response[,-(1:(lags-1)),]
  }
  response <- aperm(response, c(1, 3, 2))
  dimnames(response) <- list(dimnB[[1]], dimnames(smat)[[2]], NULL)
  ## dimnames(response)[2] <- dimnames(smat)[1]
  ## dimnames(response)[1] <- dimnames(B)[2]
  return(response)
}

#' @name .divisors
#' @noRd
.divisors <- function (n,div) {
  div <- round(div)
  for(dd in div:1){
    if(n%%div==0) break else div<-div-1
  }
  return(div)
}
