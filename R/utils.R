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
  #make sure that W and Data names are in the same order
  cnames <- names(Data)
  W      <- lapply(W,function(x) x[cnames,cnames])
  
  xglobal <- c()
  for(jj in 1:length(Data)){
    pretemp <- Data[[jj]];class(pretemp) <- "numeric"
    temp <- as.matrix(pretemp[,colSums(is.na(pretemp))<nrow(pretemp)])
    colnames(temp) <- paste(cnames[jj],colnames(pretemp),sep=".")
    xglobal <- cbind(xglobal,temp)
  }
  
  Max.char<-max(nchar(colnames(xglobal)))
  cnt.char<-max(nchar(cnames))
  # for each country
  gW<-list()
  for(cc in 1:length(cnames)){
    xglobal <- xglobal[,!duplicated(colnames(xglobal))]
    #creates a dynamic list of variables
    varnames <- substr(colnames(xglobal),cnt.char+2,Max.char); varnames <- varnames[!duplicated(varnames)] 
    
    if(!is.null(Wex.restr)){
      varnames = varnames[-charmatch(Wex.restr,varnames)]
    }
    
    # names of endo variables
    endnames <- unlist(lapply(strsplit(colnames(xglobal)[grepl(paste0("^",cnames[[cc]]),colnames(xglobal))],".",fixed=TRUE),
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
        Wnew[kk,var.cntry.indic] = 1
        # if all zero then set zero
        if(all(W$W[cc,] == 0)) Wnew[kk,var.cntry.indic] = 0
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
    
    # WfinNR <- WfinNR[!(rowSums(abs(WfinNR)) == 0),]  # only chang MB 14/01/23 !!!! 
    #gW[[cc]]<-apply(WfinNR,2,function(x) x/rowSums(WfinNR))
    
    zero_rows <- apply(Wnew,1,sum) == 0
    if(sum(zero_rows) != nrow(Wnew)){
      Wnew <- Wnew[!zero_rows,]
      Wnew <- apply(Wnew,2,function(x)x/rowSums(Wnew))
    }
    
    WfinNR <- rbind(endoW,Wnew)
    gW[[cc]] <- WfinNR
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
        Wnew           = matrix(0,length(OE.vars[[kk]]),ncol(xglobal))
        colnames(Wnew) = colnames(xglobal)
        rownames(Wnew) = c(paste(OE.cN[kk],".",OE.vars[[kk]][!OE.vars[[kk]]%in%names(endo)],sep=""),
                           OE.vars[[kk]][OE.vars[[kk]]%in%names(endo)])
        if(OE.x>1){
          diag(Wnew[paste(OE.cN[kk],".",OEnames,sep=""),paste(OE.cN[kk],".",OEnames,sep="")])<-1
        }else{
          Wnew[paste(OE.cN[kk],".",OEnames,sep=""),paste(OE.cN[kk],".",OEnames,sep="")]<-1
        }
        vars <- OE.vars[[kk]][OE.vars[[kk]]%in%names(endo)]
        if(length(vars)>0){
          for(i in 1:length(vars)){
            Wnew[vars[i],paste(names(OE.weights[[kk]]),".",vars[i],sep="")] <- OE.weights[[kk]]
            Wnew[vars[i],]
          }
        }
        # this creates the part if there are more than one other entities
        if(OE.sets>1 && kk < OE.sets){
          aux <- Wnew
          xx <- lapply(OE[(kk+1):OE.sets],function(l)ncol(l))
          xx <- sum(unlist(xx))
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
.get_V <- function(k=k,M=M,Mstar,plag,plagstar,lambda1,lambda2,lambda3,lambda4,sigma_sq,sigma_wex,trend=FALSE,wexo=TRUE){
  V_i <- matrix(0,k,M)
  # endogenous part
  for(i in 1:M){
    for(pp in 1:plag){
      for(j in 1:M){
        if(i==j){
          #V_i[j+M*(pp-1),i] <- a_bar_1/(pp^2) ######
          V_i[j+M*(pp-1),i] <- (lambda1/pp)^2
        }else{
          #V_i[j+M*(pp-1),i] <- (a_bar_2 * sigma_sq[i])/(pp^2*sigma_sq[j]) #####
          V_i[j+M*(pp-1),i] <- (lambda1*lambda2/pp)^2 * (sigma_sq[i]/sigma_sq[j])
        }
      }
    }
  }
  # exogenous part
  if(wexo){
    for(i in 1:M){
      for(pp in 0:plagstar){
        for(j in 1:Mstar){
          #V_i[M*p+pp*Mstar+j,i] <- a_bar_4 * sigma_sq[i]/(sigma_wex[j]*(pp+1)) #####
          V_i[M*plag+pp*Mstar+j,i] <- (lambda1*lambda3/(pp+1))^2 * (sigma_sq[i]/sigma_wex[j])
        }
      }
    }
  }
  # deterministics
  for(i in 1:M){
    if(trend){
      V_i[(k-1):k,i] <- lambda4 * sigma_sq[i]
    }else{
      V_i[k,i] <- lambda4 * sigma_sq[i]
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

#' @name .get_nrc
#' @noRd
.get_nrc <- function(k){
  if(k==1) return(c(1,1))
  if(k==2) return(c(2,1))
  if(k%in%c(3,4)) return(c(2,2))
  if(k%in%c(5,6)) return(c(3,2))
  if(k>6) return(c(3,3))
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
.BVAR_linear_wrapper <- function(cc, cN, xglobal, gW, prior, lags, draws, burnin, trend, SV, thin, default_hyperpara, Ex, use_R, setting_store){
  Yraw  = xglobal[,substr(colnames(xglobal),1,2)==cN[cc],drop=FALSE]; class(Yraw) = "numeric"
  W     = gW[[cc]]
  Exraw = matrix(NA_real_)
  if(!is.null(Ex)) if(cN[cc]%in%names(Ex)) Exraw <- Ex[[cN[cc]]]
  all         = t(W%*%t(xglobal))
  if(ncol(Yraw) == ncol(all)){
    Wraw = NULL
    Mstar = 0
  }else{
    Wraw  = all[,(ncol(Yraw)+1):ncol(all),drop=FALSE]; class(Wraw) = "numeric"
    Mstar = ncol(Wraw)
  }
  
  if(all(Wraw==0)){ # case of no exogenous variables -- always uses R version (not implemented in Rcpp)
    Wraw  = NULL
    wexo  = FALSE
  }else{
    wexo  = TRUE
  }
  default_hyperpara$Mstar <- Mstar
  prior_in <- ifelse(prior=="MN",1,ifelse(prior=="SSVS",2,ifelse(prior=="NG",3,4)))
  if(default_hyperpara[["tau_log"]]){
    if(ncol(Yraw)>1) default_hyperpara["tau_theta"] <- 1/log(ncol(Yraw))
  }
  # estimation
  if(!use_R){
    # Rcpp::sourceCpp("./src/BVAR_linear.cpp")
    invisible(capture.output(bvar<-try(BVAR_linear(Yraw,Wraw,Exraw,lags,as.integer(draws),as.integer(burnin),
                                                   as.integer(thin),TRUE,trend,SV,as.integer(prior_in),
                                                   default_hyperpara,setting_store)), type="message"))
  }else{
    bvar <- structure("message",class=c("try-error","character"))
  }
  if(is(bvar,"try-error")){
    # Rcpp::sourceCpp("./src/do_rgig1.cpp")
    bvar<-try(.BVAR_linear_R(Yraw,Wraw,Exraw,lags,draws,burnin,thin,TRUE,trend,SV,
                             prior_in,default_hyperpara,TRUE,setting_store), silent=TRUE)
  }
  # error handling
  if(inherits(bvar,"try-error")){
    message("\nBGVAR incurred an error when estimating the models.\nSee original error message:")
    message(paste0("Error occured in countrymodel: ", cN[cc],". Please check."))
    message("Error in detail: \n")
    message(bvar)
    #stop()
  }
  #------------------------------------------------ get data ----------------------------------------#
  Y <- bvar$Y; colnames(Y) <- colnames(Yraw); X <- bvar$X
  M <- ncol(Y); bigT <- nrow(Y); K <- ncol(X)
  plag <- lags[1]; plagstar <- lags[2]; pmax <- max(lags)
  if(!any(is.na(Exraw))) Mex <- ncol(Exraw)
  if(wexo){
    xnames <- c(paste(rep("Ylag",M),rep(seq(1,plag),each=M),sep=""),rep("Wex",Mstar),
                paste(rep("Wexlag",Mstar),rep(seq(1,plagstar),each=Mstar),sep=""))
    if(!any(is.na(Exraw))) xnames <- c(xnames,paste(rep("Tex",Mex)))
    xnames <- c(xnames,"cons")
    if(trend) xnames <- c(xnames,"trend")
    xnames_end <- xnames
  }else{
    xnames <- c(paste0(rep("Ylag",M),rep(seq(1,plag),each=M),sep=""))
    if(!any(is.na(Exraw))) xnames <- c(xnames,paste(rep("Tex",Mex)))
    xnames <- c(xnames,"cons")
    if(trend) xnames <- c(xnames,"trend")
    
    xnames_end <- c(paste(rep("Ylag",M),rep(seq(1,plag),each=M),sep=""),rep("Wex",Mstar),
                    paste(rep("Wexlag",Mstar),rep(seq(1,plagstar),each=Mstar),sep=""))
    if(!any(is.na(Exraw))) xnames_end <- c(xnames_end,paste(rep("Tex",Mex)))
    xnames_end <- c(xnames_end,"cons")
    if(trend) xnames_end <- c(xnames_end,"trend")
  }
  colnames(X) <- xnames
  #-----------------------------------------get containers ------------------------------------------#
  A_store <- bvar$A_store; dimnames(A_store)[[1]] <- xnames_end; dimnames(A_store)[[2]] <- colnames(Y)
  # splitting up stores
  dims          <- dimnames(A_store)[[1]]
  a0store       <- adrop(A_store[which(dims=="cons"),,,drop=FALSE],drop=1)
  a1store <- Exstore <- NULL
  if(trend){
    a1store     <- adrop(A_store[which(dims=="trend"),,,drop=FALSE],drop=1)
  }
  if(!any(is.na(Exraw))){
    Exstore     <- A_store[which(dims=="Tex"),,,drop=FALSE]
  }
  Lambda0store  <- A_store[which(dims=="Wex"),,,drop=FALSE]
  Lambdastore   <- NULL
  Phistore      <- NULL
  for(pp in 1:pmax){
    if(pp %in% seq(plag)){
      Phistore[[pp]] <- A_store[which(dims==paste("Ylag",pp,sep="")),,,drop=FALSE]
    }else{
      Phistore[[pp]] <- array(0, c(M, M, draws/thin), 
                              dimnames=list(rep(paste0("Ylag",pp),M),colnames(Y),NULL))
    }
    if(pp %in% seq(plagstar)){
      Lambdastore[[pp]] <- A_store[which(dims==paste("Wexlag",pp,sep="")),,,drop=FALSE]
    }else{
      Lambdastore[[pp]] <- array(0, c(Mstar, M, draws/thin),
                                 dimnames=list(rep(paste0("Wexlag",pp),Mstar),colnames(Y),NULL))
    }
  }
  SIGMA_store <- array(NA, c(bigT,M,M,draws/thin)); dimnames(SIGMA_store) <- list(NULL,colnames(Y),colnames(Y),NULL)
  L_store <- bvar$L_store
  for(irep in 1:(draws/thin)){
    for(tt in 1:bigT){
      if(M>1){
        SIGMA_store[tt,,,irep] <- L_store[,,irep]%*%diag(exp(bvar$Sv_store[tt,,irep]))%*%t(L_store[,,irep])
      }else{
        SIGMA_store[tt,,,irep] <- L_store[,,irep]%*%exp(bvar$Sv_store[tt,,irep])%*%t(L_store[,,irep])
      }
    }
  }
  SIGMAmed_store <- apply(SIGMA_store, c(2,3,4), median)
  res_store     <- bvar$res_store; dimnames(res_store) <- list(NULL,colnames(Y),NULL)
  if(SV){
    vola_store  <- bvar$Sv_store; dimnames(vola_store) <- list(NULL,colnames(Y),NULL)
    vola_post   <- apply(vola_store,c(1,2),median)
  }else{
    vola_store  <- bvar$Sv_store; 
    vola_post   <- apply(vola_store,c(1,2),median)
  }
  # additional stuff
  if(SV & setting_store$vola_pars){
    pars_store  <- bvar$pars_store
    pars_post   <- apply(pars_store,c(1,2),median)
  }else{
    pars_store <- pars_post <- NULL 
  }
  # MN
  if(prior=="MN" & setting_store$shrink_MN){
    lambda_store  <- bvar$MN$lambda_store; dimnames(lambda_store) <- list(c("lambda1","lambda2","lambda4"),NULL,NULL)
    lambda_post   <- apply(lambda_store,c(1,2),median)
  }else{
    lambda_store  <- lambda_post <- NULL
  }
  # SSVS
  if(prior=="SSVS" & setting_store$shrink_SSVS){
    gamma_store <- bvar$SSVS$gamma_store; dimnames(gamma_store) <- list(xnames_end,colnames(Y),NULL)
    omega_store <- bvar$SSVS$omega_store; dimnames(omega_store) <- list(colnames(Y),colnames(Y),NULL)
    PIP         <- apply(gamma_store,c(1,2),mean)
    PIP_omega   <- apply(omega_store,c(1,2),mean)
  }else{
    gamma_store <- omega_store <- PIP <- PIP_omega <- NULL
  }
  # NG
  if(prior=="NG" & setting_store$shrink_NG){
    theta_store   <- bvar$NG$theta_store; dimnames(theta_store)[[1]] <- colnames(X); dimnames(theta_store)[[2]] <- colnames(Y)
    lambda2_store <- bvar$NG$lambda2_store
    tau_store     <- bvar$NG$tau_store
    dimnames(lambda2_store) <- list(paste("lag",0:plag,sep="_"),c("endogenous","weakly exogenous","covariance"),NULL)
    dimnames(lambda2_store) <- list(paste("lag",0:plag,sep="_"),c("endogenous","weakly exogenous","covariance"),NULL)
    theta_post  <- apply(theta_store,c(1,2),median)
    lambda2_post  <- apply(lambda2_store,c(1,2),median)
    tau_post      <- apply(tau_store,c(1,2),median)
  }else{
    theta_store <- lambda2_store <- tau_store <- theta_post <- lambda2_post <- tau_post <- NULL
  }
  # HS
  lambda_A_endo_store <- lambda_A_exo_store <- lambda_L_store <- NULL
  nu_A_endo_store     <- nu_A_exo_store     <- nu_L_store     <- NULL
  tau_A_endo_store    <- tau_A_exo_store    <- tau_L_store    <- NULL
  zeta_A_endo_store   <- zeta_A_exo_store   <- zeta_L_store   <- NULL
  lambda_A_endo_post  <- lambda_A_exo_post  <- lambda_L_post  <- NULL
  nu_A_endo_post      <- nu_A_exo_post      <- nu_L_post      <- NULL
  tau_A_endo_post     <- tau_A_exo_post     <- tau_L_post     <- NULL
  zeta_A_endo_post    <- zeta_A_exo_post    <- zeta_L_post    <- NULL
  if(prior=="HS" & setting_store$shrink_HS){
    lambda_A_endo_store <- bvar$HS$lambda_A_endo_store
    lambda_A_exo_store  <- bvar$HS$lambda_A_exo_store
    lambda_L_store      <- bvar$HS$lambda_L_store
    nu_A_endo_store     <- bvar$HS$nu_A_endo_store
    nu_A_exo_store      <- bvar$HS$nu_A_exo_store
    nu_L_store          <- bvar$HS$nu_L_store
    tau_A_endo_store    <- bvar$HS$tau_A_endo_store
    tau_A_exo_store     <- bvar$HS$tau_A_exo_store
    tau_L_store         <- bvar$HS$tau_L_store
    zeta_A_endo_store   <- bvar$HS$zeta_A_endo_store
    zeta_A_exo_store    <- bvar$HS$zeta_A_exo_store
    zeta_L_store        <- bvar$HS$zeta_L_store
    
    lambda_A_endo_post  <- apply(lambda_A_endo_store, 1, median)
    lambda_A_exo_post   <- apply(lambda_A_exo_store, 1, median) 
    lambda_L_post       <- apply(lambda_L_store, 1, median)
    nu_A_endo_post      <- apply(nu_A_endo_store, 1, median)
    nu_A_exo_post       <- apply(nu_A_exo_store, 1, median)
    nu_L_post           <- apply(nu_L_store, 1, median)
    tau_A_endo_post     <- apply(tau_A_endo_store, 1, median)
    tau_A_exo_post      <- apply(tau_A_exo_store, 1, median)
    tau_L_post          <- apply(tau_L_store, 1, median)
    zeta_A_endo_post    <- apply(zeta_A_endo_store, 1, median)
    zeta_A_exo_post     <- apply(zeta_A_exo_store, 1, median)
    zeta_L_post         <- apply(zeta_L_store, 1, median)
  }
  #------------------------------------ compute posteriors -------------------------------------------#
  A_post      <- apply(A_store, c(1,2), median)
  L_post      <- apply(L_store, c(1,2), median)
  SIGMA_post  <- apply(SIGMA_store,c(1,2,3),median)
  S_post      <- apply(SIGMA_post,c(1,2),mean)
  Sig         <- S_post/(bigT-K)
  res_post    <- apply(res_store,c(1,2),median)
  # splitting up posteriors
  a0post      <- A_post[which(dims=="cons"),,drop=FALSE]
  a1post <- Expost <- NULL
  if(trend){
    a1post    <- A_post[which(dims=="trend"),,drop=FALSE]
  }
  if(!any(is.na(Exraw))){
    Expost    <- A_post[which(dims=="Tex"),,drop=FALSE]
  }
  Lambda0post <- A_post[which(dims=="Wex"),,drop=FALSE]
  Lambdapost  <- NULL
  Phipost     <- NULL
  for(pp in 1:pmax){
    if(pp %in% seq(plag)){
      Phipost <- rbind(Phipost,A_post[which(dims==paste("Ylag",pp,sep="")),,drop=FALSE])
    }else{
      Phipost <- rbind(Phipost, matrix(0, M, M, dimnames=list(rep(paste0("Ylag",pp),M),colnames(Y))))
    }
    if(pp %in% seq(plagstar)){
      Lambdapost <- rbind(Lambdapost,A_post[which(dims==paste("Wexlag",pp,sep="")),,drop=FALSE])
    }else{
      Lambdapost <- rbind(Lambdapost, matrix(0, Mstar, M, dimnames=list(rep(paste0("Wexlag",pp),Mstar),colnames(Y))))
    }
  }
  post <- list(A_post=A_post,a0post=a0post,a1post=a1post,Lambda0post=Lambda0post,Lambdapost=Lambdapost,
               Phipost=Phipost,Expost=Expost,S_post=S_post,Sig=Sig,theta_post=theta_post,L_post=L_post,
               SIGMA_post=SIGMA_post,
               vola_post=vola_post,pars_post=pars_post,res_post=res_post,lambda_post=lambda_post,
               PIP=PIP,PIP_omega=PIP_omega,lambda2_post=lambda2_post,tau_post=tau_post,
               lambda_A_endo_post=lambda_A_endo_post,lambda_A_exo_post=lambda_A_exo_post,lambda_L_post=lambda_L_post,
               nu_A_endo_post=nu_A_endo_post,nu_A_exo_post=nu_A_exo_post,nu_L_post=nu_L_post,
               tau_A_endo_post=tau_A_endo_post,tau_A_exo_post=tau_A_exo_post,tau_L_post=tau_L_post,
               zeta_A_endo_post=zeta_A_endo_post,zeta_A_exo_post=zeta_A_exo_post,zeta_L_post=zeta_L_post)
  
  store <- list(a0store=a0store,a1store=a1store,Lambda0store=Lambda0store,Lambdastore=Lambdastore,
                Phistore=Phistore,Exstore=Exstore,SIGMAmed_store=SIGMAmed_store,
                L_store=L_store,theta_store=theta_store,vola_store=vola_store,pars_store=pars_store,
                res_store=res_store,lambda_store=lambda_store,gamma_store=gamma_store,omega_store=omega_store,
                lambda2_store=lambda2_store,tau_store=tau_store,
                lambda_A_endo_store=lambda_A_endo_store,lambda_A_exo_store=lambda_A_exo_store,lambda_L_store=lambda_L_store,
                nu_A_endo_store=nu_A_endo_store,nu_A_exo_store=nu_A_exo_store,nu_L_store=nu_L_store,
                tau_A_endo_store=tau_A_endo_store,tau_A_exo_store=tau_A_exo_store,tau_L_store=tau_L_store,
                zeta_A_endo_store=zeta_A_endo_store,zeta_A_exo_store=zeta_A_exo_store,zeta_L_store=zeta_L_store)
  out  <- list(Y=Y,X=X,W=W,store=store,post=post)
  return(out)
}

#' @name .BVAR_linear_R
#' @importFrom stochvol svsample_fast_cpp specify_priors get_default_fast_sv
#' @importFrom MASS ginv mvrnorm
#' @importFrom methods is
#' @importFrom stats rnorm rgamma runif dnorm
#' @noRd
.BVAR_linear_R <- function(Yraw,Wraw,Exraw,lags,draws,burnin,thin,cons,trend,sv,prior,hyperpara,verbose,setting_store){
  #----------------------------------------INPUTS----------------------------------------------------#
  plag     <- lags[1]
  plagstar <- lags[2]
  pmax     <- max(lags)
  
  Traw  <- nrow(Yraw)
  M     <- ncol(Yraw)
  K     <- M*plag
  Ylag  <- .mlag(Yraw,plag)
  nameslags <- NULL
  for (ii in 1:plag) nameslags <- c(nameslags,rep(paste("Ylag",ii,sep=""),M))
  colnames(Ylag) <- nameslags
  
  Mstar        <- hyperpara$Mstar
  Kstar        <- Mstar*(plagstar+1)
  wexnames     <- rep("Wex",Mstar)
  wexnameslags <- NULL
  for (ii in 1:plagstar) wexnameslags <- c(wexnameslags,rep(paste("Wexlag",ii,sep=""),Mstar))
  if(!is.null(Wraw)){
    wexo             <- TRUE
    Wexlag           <- .mlag(Wraw,plagstar)
    colnames(Wraw)   <- wexnames
    colnames(Wexlag) <- wexnameslags
  }else{
    wexo             <- FALSE
    Wexlag           <- NULL
  }
  
  texo <- FALSE; Mex <- 0; exnames <- NULL
  if(nrow(Exraw) != 1){
    Mex             <- ncol(Exraw)
    texo            <- TRUE
    exnames         <- rep("Tex",Mex)
    colnames(Exraw) <- exnames
  }
  nameslags_end <- c(nameslags,wexnames,wexnameslags,exnames)
  
  Xraw          <- cbind(Ylag,Wraw,Wexlag)
  if(texo) Xraw <- cbind(Xraw,Exraw)
  X             <- Xraw[(pmax+1):nrow(Xraw),,drop=FALSE]
  Y             <- Yraw[(pmax+1):Traw,,drop=FALSE]
  bigT          <- nrow(X)
  
  if(cons){
    X                    <- cbind(X,1)
    colnames(X)[ncol(X)] <- "cons"
    nameslags_end        <- c(nameslags_end,"cons")
  }
  if(trend){
    X                    <- cbind(X,seq(1,bigT))
    colnames(X)[ncol(X)] <- "trend"
    nameslags_end        <- c(nameslags_end,"trend")
  }
  
  k     <- ncol(X)
  k_end <- ncol(Ylag) + Kstar + ifelse(cons,1,0) + ifelse(trend,1,0)
  v     <- (M*(M-1))/2
  n     <- K*M
  nstar <- Kstar*M
  #---------------------------------------------------------------------------------------------------------
  # HYPERPARAMETERS
  #---------------------------------------------------------------------------------------------------------
  prmean      = hyperpara$prmean
  a_1         = hyperpara$a_1
  b_1         = hyperpara$b_1
  crit_eig    = hyperpara$crit_eig
  Bsigma      = hyperpara$Bsigma
  a0          = hyperpara$a0
  b0          = hyperpara$b0
  bmu         = hyperpara$bmu
  Bmu         = hyperpara$Bmu
  # prior == 1: MN
  lambda1     = hyperpara$lambda1
  lambda2     = hyperpara$lambda2
  lambda3     = hyperpara$lambda3
  lambda4     = hyperpara$lambda4
  # prior == 2: SSVS
  tau00       = hyperpara$tau0
  tau11       = hyperpara$tau1
  p_i         = hyperpara$p_i
  kappa0      = hyperpara$kappa0
  kappa1      = hyperpara$kappa1
  q_ij        = hyperpara$q_ij
  # prior == 3: NG
  d_lambda    = hyperpara$d_lambda
  e_lambda    = hyperpara$e_lambda
  tau_theta   = hyperpara$tau_theta
  sample_tau  = hyperpara$sample_tau
  #---------------------------------------------------------------------------------------------------------
  # STORE SETTINGS
  #---------------------------------------------------------------------------------------------------------
  save_shrink_MN   = setting_store$shrink_MN
  save_shrink_SSVS = setting_store$shrink_SSVS
  save_shrink_NG   = setting_store$shrink_NG
  save_shrink_HS   = setting_store$shrink_HS
  save_vola_pars   = setting_store$vola_pars
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
  A_draw    <- A_OLS
  SIGMA     <- array(SIGMA_OLS, c(M,M,bigT))
  Em        <- Em_str <- E_OLS
  L_draw    <- diag(M)
  L_drawinv <- diag(M)
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
  accept3 <- 0
  scale1  <- .43
  scale2  <- .43
  scale3  <- .43
  sigma_sq  <- matrix(0,M,1) #vector which stores the residual variance
  for (i in 1:M){
    Ylag_i        <- .mlag(Yraw[,i],plag)
    Ylag_i        <- Ylag_i[(plag+1):nrow(Ylag_i),,drop=FALSE]
    Y_i           <- Yraw[(plag+1):nrow(Yraw),i,drop=FALSE]
    Ylag_i        <- cbind(Ylag_i,seq(1,nrow(Y_i)))
    alpha_i       <- solve(crossprod(Ylag_i))%*%crossprod(Ylag_i,Y_i)
    sigma_sq[i,1] <- (1/(nrow(Y_i)-plag-1))*t(Y_i-Ylag_i%*%alpha_i)%*%(Y_i-Ylag_i%*%alpha_i)
  }
  if(wexo){
    sigma_wex <- matrix(0,Mstar,1)
    for (j in 1:Mstar){
      Ywex_i <- .mlag(Wraw[,j],plagstar)
      Ywex_i <- Ywex_i[(plag+1):Traw,]
      Yw_i   <- Wraw[(plag+1):Traw,j,drop=FALSE]
      Ywex_i <- cbind(Ywex_i,seq(1,nrow(Yw_i)))
      alpha_w <- solve(crossprod(Ywex_i))%*%t(Ywex_i)%*%Yw_i
      sigma_wex[j,1] <- (1/(nrow(Yw_i)-plag-1))*t(Yw_i-Ywex_i%*%alpha_w)%*%(Yw_i-Ywex_i%*%alpha_w)
    }
  }else{
    sigma_wex <- NULL
  }
  
  # MN prior
  if(prior == 1){
    theta <- .get_V(k=k,M=M,Mstar,plag,plagstar,lambda1,lambda2,lambda3,lambda4,sigma_sq,sigma_wex,trend,wexo)
    post1 <- sum(dnorm(as.vector(A_draw),a_prior,sqrt(as.vector(theta)),log=TRUE))+dgamma(lambda1,0.01,0.01,log=TRUE)+log(lambda1) # correction term
    post2 <- sum(dnorm(as.vector(A_draw),a_prior,sqrt(as.vector(theta)),log=TRUE))+dgamma(lambda2,0.01,0.01,log=TRUE)+log(lambda2) # correction term
    post3 <- sum(dnorm(as.vector(A_draw),a_prior,sqrt(as.vector(theta)),log=TRUE))+dgamma(lambda3,0.01,0.01,log=TRUE)+log(lambda3) # correction term
  }
  
  # SSVS prior
  if(prior == 2){
    gamma  <-  matrix(1,k,M)
    sigma_alpha  <-  sqrt(diag(kronecker(SIGMA_OLS,XtXinv)))
    tau0 <- matrix(NA_real_, k, M); tau1 <- matrix(NA_real_, k, M)
    ii <- 1
    for(mm in 1:M){
      for(kk in 1:k){
        tau0[kk,mm] <- tau00*sigma_alpha[ii]
        tau1[kk,mm] <- tau11*sigma_alpha[ii]
        ii <- ii+1
      }
    }
  }
  
  # NG stuff
  if(prior == 3){
    lambda2_A       <- matrix(0.01,pmax+1,2)
    A_tau           <- matrix(tau_theta,pmax+1,2)
    colnames(A_tau) <- colnames(lambda2_A) <- c("endo","exo")
    rownames(A_tau) <- rownames(lambda2_A) <- paste("lag.",seq(0,pmax),sep="")
    A_tuning        <- matrix(.43,pmax+1,2)
    A_accept        <- matrix(0,pmax+1,2)
    lambda2_A[1,1]  <- A_tau[1,1] <- A_tuning[1,1] <- A_accept[1,1] <- NA
  }

  # HS stuff
  if(prior == 4){
    lambda_A_endo <- nu_A_endo <- rep(1, n)
    lambda_A_exo  <- nu_A_exo  <- rep(1, nstar)
    lambda_L      <- nu_L      <- rep(1, v)
    tau_A_endo  <- tau_A_exo  <- tau_L <- 1
    zeta_A_endo <- zeta_A_exo <- zeta_L <- 1
  }
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
  L_tau     <- tau_theta
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
  ntot  <- draws+burnin
  
  # thinning
  count        <- 0
  thindraws    <- draws/thin
  thin.draws   <- seq(burnin+1,ntot,by=thin)
  #---------------------------------------------------------------------------------------------------------
  # STORAGES
  #---------------------------------------------------------------------------------------------------------
  A_store      <- array(0,        c(k_end,M,thindraws))
  L_store      <- array(NA_real_, c(M,M,thindraws))
  res_store    <- array(NA_real_, c(bigT,M,thindraws))
  # SV
  Sv_store     <- array(NA_real_, c(bigT,M,thindraws))
  if(save_vola_pars){
    pars_store <- array(NA_real_, c(4,M,thindraws))
  }else{
    pars_store <- NULL
  }
  # MN
  if(save_shrink_MN){
    lambda_store <- array(NA_real_, c(3,1,thindraws))
  }else{
    lambda_store <- NULL
  }
  # SSVS
  if(save_shrink_SSVS){
    gamma_store  <- array(0,        c(k_end,M,thindraws))
    omega_store  <- array(NA_real_, c(M,M,thindraws))
  }else{
    gamma_store <- omega_store <- NULL
  }
  # NG
  if(save_shrink_NG){
    theta_store  <- array(0,        c(k_end,M,thindraws))
    lambda2_store<- array(NA_real_, c(pmax+1,3,thindraws))
    tau_store    <- array(NA_real_, c(pmax+1,3,thindraws))
  }else{
    theta_store <- lambda2_store <- tau_store <- NULL
  }
  # HS
  if(save_shrink_HS){
    lambda_A_endo_store <- array(0, c(n, thindraws))
    lambda_A_exo_store  <- array(0, c(nstar, thindraws))
    lambda_L_store      <- array(0, c(v, thindraws))
    nu_A_endo_store     <- array(0, c(n, thindraws))
    nu_A_exo_store      <- array(0, c(nstar, thindraws))
    nu_L_store          <- array(0, c(v, thindraws))
    tau_A_endo_store    <- array(0, c(1, thindraws))
    tau_A_exo_store     <- array(0, c(1, thindraws))
    tau_L_store         <- array(0, c(1, thindraws))
    zeta_A_endo_store   <- array(0, c(1, thindraws))
    zeta_A_exo_store    <- array(0, c(1, thindraws))
    zeta_L_store        <- array(0, c(1, thindraws))
  }else{
    lambda_A_endo_store <- lambda_A_exo_store <- lambda_L_store <- NULL
    nu_A_endo_store     <- nu_A_exo_store     <- nu_L_store     <- NULL
    tau_A_endo_store    <- tau_A_exo_store    <- tau_L_store    <- NULL
    zeta_A_endo_store   <- zeta_A_exo_store   <- zeta_L_store   <- NULL
  }
  #---------------------------------------------------------------------------------------------------------
  # MCMC LOOP
  #---------------------------------------------------------------------------------------------------------
  for(irep in 1:ntot){
    #----------------------------------------------------------------------------
    # Step 1: Sample coefficients
    for(mm in 1:M){
      A0_draw = A_draw
      A0_draw[,mm] <- 0
      
      ztilde <- as.vector((Y - X%*%A0_draw)%*%t(L_drawinv[mm:M,,drop=FALSE])) * exp(-0.5*as.vector(Sv_draw[,mm:M,drop=FALSE]))
      xtilde <- (L_drawinv[mm:M,mm,drop=FALSE] %x% X) * exp(-0.5*as.vector(Sv_draw[,mm:M,drop=FALSE]))
      
      V_post <- try(chol2inv(chol(crossprod(xtilde)+diag(1/theta[,mm]))),silent=TRUE)
      if(is(V_post,"try-error")) V_post <- try(solve(crossprod(xtilde)+diag(1/theta[,mm])),silent=TRUE)
      if(is(V_post,"try-error")) V_post <- ginv(crossprod(xtilde)+diag(1/theta[,mm]))
      A_post <- V_post%*%(crossprod(xtilde,ztilde)+diag(1/theta[,mm])%*%A_prior[,mm])
      
      A.draw.i <- try(A_post+t(chol(V_post))%*%rnorm(ncol(X)),silent=TRUE)
      if(is(A.draw.i,"try-error")) A.draw.i <- t(mvrnorm(1,A_post,V_post))
      A_draw[,mm] <- A.draw.i
      Em[,mm] <- Y[,mm]-X%*%A.draw.i
    }
    rownames(A_draw) <- colnames(X)
    # Step 1b: Sample coefficients in L matrix
    if(M > 1){
      for(mm in 2:M){
        eps.m <- Em[,mm]*exp(-0.5*Sv_draw[,mm,drop=TRUE])
        eps.x <- Em[,1:(mm-1),drop=FALSE]*exp(-0.5*Sv_draw[,mm,drop=TRUE])
        
        L_post <- try(chol2inv(chol(crossprod(eps.x)+diag(1/L_prior[mm,1:(mm-1)],mm-1,mm-1))),silent=TRUE)
        if(is(L_post,"try-error")) L_post <- try(solve(crossprod(eps.x)+diag(1/L_prior[mm,1:(mm-1)],mm-1,mm-1)),silent=TRUE)
        if(is(L_post,"try-error")) L_post <- ginv(crossprod(eps.x)+diag(1/L_prior[mm,1:(mm-1)],mm-1,mm-1))
        l_post <- L_post%*%(crossprod(eps.x,eps.m)+diag(1/L_prior[mm,1:(mm-1)],mm-1,mm-1)%*%l_prior[mm,1:(mm-1)])
        
        L.draw.i <- try(l_post+t(chol(L_post))%*%rnorm(length(1:(mm-1))),silent=TRUE)
        if(is(L.draw.i,"try-error")) L.draw.i <- t(mvrnorm(1,l_post,L_post))
        L_draw[mm,1:(mm-1)] <- L.draw.i
      }
    }
    # Step 1c: Compute Em_str
    L_drawinv = solve(L_draw)
    Em_str    = Y%*%t(L_drawinv) - X%*%A_draw%*%t(L_drawinv)
    # for (mm in 1:M){
    #   if (mm==1){
    #     Y.i <- Y[,mm]*exp(-0.5*Sv_draw[,mm])
    #     X.i <- X*exp(-0.5*Sv_draw[,mm])
    #     
    #     V_post <- try(chol2inv(chol(crossprod(X.i)+diag(1/theta[,mm]))),silent=TRUE)
    #     if (is(V_post,"try-error")) V_post <- ginv(crossprod(X.i)+diag(1/theta[,mm]))
    #     A_post <- V_post%*%(crossprod(X.i,Y.i)+diag(1/theta[,mm])%*%A_prior[,mm])
    #     
    #     A.draw.i <- try(A_post+t(chol(V_post))%*%rnorm(ncol(X.i)),silent=TRUE)
    #     if (is(A.draw.i,"try-error")) A.draw.i <- mvrnorm(1,A_post,V_post)
    #     A_draw[,mm] <- A.draw.i
    #     Em[,mm] <-  Em_str[,mm] <- Y[,mm]-X%*%A.draw.i
    #   }else{
    #     Y.i <- Y[,mm]*exp(-0.5*Sv_draw[,mm])
    #     X.i <- cbind(X,Em[,1:(mm-1)])*exp(-0.5*Sv_draw[,mm])
    #     
    #     V_post <- try(chol2inv(chol((crossprod(X.i)+diag(1/c(theta[,mm],L_prior[mm,1:(mm-1)]))))),silent=TRUE)
    #     if (is(V_post,"try-error")) V_post <- ginv((crossprod(X.i)+diag(1/c(theta[,mm],L_prior[mm,1:(mm-1)]))))
    #     A_post <- V_post%*%(crossprod(X.i,Y.i)+diag(1/c(theta[,mm],L_prior[mm,1:(mm-1)]))%*%c(A_prior[,mm],l_prior[mm,1:(mm-1)]))
    #     
    #     A.draw.i <- try(A_post+t(chol(V_post))%*%rnorm(ncol(X.i)),silent=TRUE)
    #     if (is(A.draw.i,"try-error")) A.draw.i <- mvrnorm(1,A_post,V_post)
    #     
    #     A_draw[,mm] <- A.draw.i[1:ncol(X)]
    #     Em[,mm] <- Y[,mm]-X%*%A.draw.i[1:ncol(X)]
    #     Em_str[,mm] <- Y[,mm]-X%*%A.draw.i[1:ncol(X)]-Em[,1:(mm-1),drop=FALSE]%*%A.draw.i[(ncol(X)+1):ncol(X.i),drop=FALSE]
    #     L_draw[mm,1:(mm-1)] <- A.draw.i[(ncol(X)+1):ncol(X.i)]
    #   }
    # }
    # rownames(A_draw) <- colnames(X)
    #----------------------------------------------------------------------------
    # Step 2: different shrinkage prior setups
    # MN
    if(prior==1){
      #Step for the first shrinkage parameter (own lags)
      lambda1.prop = exp(rnorm(1,0,scale1))*lambda1
      lambda1.prop = ifelse(lambda1.prop<1e-16,1e-16,lambda1.prop)
      lambda1.prop = ifelse(lambda1.prop>1e+16,1e+16,lambda1.prop)
      theta1.prop  = .get_V(k,M,Mstar,plag,plagstar,lambda1.prop,lambda2,lambda3,lambda4,sigma_sq,sigma_wex,trend,wexo)
      post1.prop   = sum(dnorm(as.vector(A_draw),a_prior,sqrt(as.vector(theta1.prop)),log=TRUE))+dgamma(lambda1.prop,0.01,0.01,log=TRUE)+log(lambda1.prop) # correction term
      if((post1.prop-post1)>log(runif(1,0,1))){
        lambda1    = lambda1.prop
        theta      = theta1.prop
        post1      = post1.prop
        accept1    = accept1+1
      }
      
      #Step for the second shrinkage parameter (cross equation)
      lambda2.prop = exp(rnorm(1,0,scale2))*lambda2
      lambda2.prop = ifelse(lambda2.prop<1e-16,1e-16,lambda2.prop)
      lambda2.prop = ifelse(lambda2.prop>1e+16,1e+16,lambda2.prop)
      theta2.prop  = .get_V(k,M,Mstar,plag,plagstar,lambda1,lambda2.prop,lambda3,lambda4,sigma_sq,sigma_wex,trend,wexo)
      post2.prop   = sum(dnorm(as.vector(A_draw),a_prior,sqrt(as.vector(theta2.prop)),log=TRUE))+dgamma(lambda2.prop,0.01,0.01,log=TRUE)+log(lambda2.prop) # correction term
      if((post2.prop-post2)>log(runif(1,0,1))){
        lambda2    = lambda2.prop
        theta      = theta2.prop
        post2      = post2.prop
        accept2    = accept2+1
      }
      
      #Step for the final shrinkage parameter (weakly exogenous)
      if(wexo){ # do only update if weakly exogenous are present
        lambda3.prop = exp(rnorm(1,0,scale3))*lambda3
        lambda3.prop = ifelse(lambda3.prop<1e-16,1e-16,lambda3.prop)
        lambda3.prop = ifelse(lambda3.prop>1e+16,1e+16,lambda3.prop)
        theta3.prop  = .get_V(k,M,Mstar,plag,plagstar,lambda1,lambda2,lambda3,lambda3.prop,sigma_sq,sigma_wex,trend,wexo)
        post3.prop   = sum(dnorm(as.vector(A_draw),a_prior,sqrt(as.vector(theta3.prop)),log=TRUE))+dgamma(lambda3.prop,0.01,0.01,log=TRUE)+log(lambda3.prop)
        if((post3.prop-post3)>log(runif(1,0,1))){
          lambda3    = lambda3.prop
          theta      = theta3.prop
          post3      = post3.prop
          accept3    = accept3+1
        }
      }
      
      if (irep<(0.5*burnin)){
        if((accept1/irep)<0.15) scale1 = 0.99*scale1
        if((accept1/irep)>0.3)  scale1 = 1.01*scale1
        if((accept2/irep)<0.15) scale2 = 0.99*scale2
        if((accept2/irep)>0.3)  scale2 = 1.01*scale2
        if((accept3/irep)<0.15) scale3 = 0.99*scale3
        if((accept3/irep)>0.3)  scale3 = 1.01*scale3
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
      if(M>1){
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
      } # END-if M>1
    }
    # NG
    if(prior==3){
      # Normal-Gamma for Covariances
      if(M>1){
        lambda2_L    <- rgamma(1,d_lambda+L_tau*v,e_lambda+L_tau/2*sum(L_prior[lower.tri(L_prior)]))
        #Step VI: Sample the prior scaling factors for covariances from GIG
        for(mm in 2:M){
          for(ii in 1:(mm-1)){
            temp <- do_rgig1(lambda = L_tau-0.5, 
                             chi    = (L_draw[mm,ii] - l_prior[mm,ii])^2, 
                             psi    = L_tau*lambda2_L)
            temp <- ifelse(temp<1e-7,1e-7,ifelse(temp>1e+7,1e+7,temp))
            L_prior[mm,ii] <- temp
          }
        }
        if(sample_tau){
          #Sample L_tau through a simple RWMH step
          L_tau_prop       = exp(rnorm(1,0,L_tuning))*L_tau
          L_tau_prop       = ifelse(L_tau_prop<1e-16,1e-16,L_tau_prop)
          L_tau_prop       = ifelse(L_tau_prop>1e+16,1e+16,L_tau_prop)
          post_L_tau_prop  = .atau_post(atau=L_tau_prop, thetas=L_prior[lower.tri(L_prior)], k=v, lambda2=lambda2_L)
          post_L_tau_old   = .atau_post(atau=L_tau,      thetas=L_prior[lower.tri(L_prior)], k=v, lambda2=lambda2_L)
          post.diff        = post_L_tau_prop-post_L_tau_old+log(L_tau_prop)-log(L_tau)
          post.diff        = ifelse(is.nan(post.diff),-Inf,post.diff)
          if (post.diff > log(runif(1,0,1))){
            L_tau      <- L_tau_prop
            L_accept   <- L_accept+1
          }
          if (irep<(0.5*burnin)){
            if ((L_accept/irep)>0.3)  L_tuning <- 1.01*L_tuning
            if ((L_accept/irep)<0.15) L_tuning <- 0.99*L_tuning
          }
        }
      } # END-if M>1
      # Norml-Gamma for weakly exogenous
      if(wexo){
        for(ss in 0:plagstar){
          if(ss==0) slct.i <- which(rownames(A_draw)=="Wex") else slct.i <- which(rownames(A_draw)==paste("Wexlag",ss,sep=""))
          A.lag.star  <- A_draw[slct.i,,drop=FALSE]
          A.lag.prior <- A_prior[slct.i,,drop=FALSE]
          theta.lag   <- theta[slct.i,,drop=FALSE]
          
          if (ss==0){
            lambda2_A[ss+1,2] <- rgamma(n     = 1,
                                        shape = d_lambda + A_tau[ss+1,2]*Mstar*M,
                                        rate  = e_lambda + A_tau[ss+1,2]/2*sum(theta.lag))
          }else{
            lambda2_A[ss+1,2] <- rgamma(n     = 1,
                                        shape = d_lambda + A_tau[ss+1,2]*Mstar^2,
                                        rate  = e_lambda + A_tau[ss+1,2]*0.5*prod(lambda2_A[1:ss,2])*sum(theta.lag))
          }
          for(ii in 1:Mstar){
            for(mm in 1:M){
              temp <- do_rgig1(lambda = A_tau[ss+1,2]-0.5,
                               chi    = (A.lag.star[ii,mm] - A.lag.prior[ii,mm])^2,
                               psi    = A_tau[ss+1,2]*prod(lambda2_A[1:(ss+1),2]))
              temp <- ifelse(temp<1e-7,1e-7,ifelse(temp>1e+7,1e+7,temp))
              theta.lag[ii,mm] <- temp
            }
          }
          theta[slct.i,] <- theta.lag
          
          if(sample_tau){
            #Sample a_tau through a simple RWMH step (on-line tuning of the MH scaling within the first 50% of the burn-in phase)
            A_tau_prop       = exp(rnorm(1,0,A_tuning[ss+1,2]))*A_tau[ss+1,2]
            A_tau_prop       = ifelse(A_tau_prop<1e-16,1e-16,A_tau_prop)
            A_tau_prop       = ifelse(A_tau_prop>1e+16,1e+16,A_tau_prop)
            post_A_tau_prop  = .atau_post(atau=A_tau_prop,    thetas=as.vector(theta.lag),lambda2 = prod(lambda2_A[1:(ss+1),2]))
            post_A_tau_old   = .atau_post(atau=A_tau[ss+1,2], thetas=as.vector(theta.lag),lambda2 = prod(lambda2_A[1:(ss+1),2]))
            post.diff        = post_A_tau_prop - post_A_tau_old + log(A_tau_prop) - log(A_tau[ss+1,2])
            post.diff        = ifelse(is.nan(post.diff),-Inf,post.diff)
            if (post.diff > log(runif(1,0,1))){
              A_tau[ss+1,2]    <- A_tau_prop
              A_accept[ss+1,2] <- A_accept[ss+1,2]+1
            }
            if (irep<(0.5*burnin)){
              if ((A_accept[ss+1,2]/irep)>0.3)  A_tuning[ss+1,2] <- 1.01*A_tuning[ss+1,2]
              if ((A_accept[ss+1,2]/irep)<0.15) A_tuning[ss+1,2] <- 0.99*A_tuning[ss+1,2]
            }
          }
        }
      }
      # Normal-Gamma for endogenous variables
      for(ss in 1:plag){
        slct.i    <- which(rownames(A_draw)==paste("Ylag",ss,sep=""))
        A.lag     <- A_draw[slct.i,,drop=FALSE]
        A.prior   <- A_prior[slct.i,,drop=FALSE]
        theta.lag <- theta[slct.i,,drop=FALSE]
        
        if(ss==1){
          lambda2_A[ss+1,1] <- rgamma(n     = 1,
                                      shape = d_lambda + A_tau[ss+1,1]*M^2,
                                      rate  = e_lambda + A_tau[ss+1,1]/2*sum(theta.lag))
        }else{
          lambda2_A[ss+1,1] <- rgamma(n     = 1,
                                      shape = d_lambda + A_tau[ss+1,1]*M^2,
                                      rate  = e_lambda + A_tau[ss+1,1]/2*prod(lambda2_A[2:(ss+1),1])*sum(theta.lag))
        }
        for(ii in 1:M){
          for(mm in 1:M){
            temp <- do_rgig1(lambda = A_tau[ss+1,1] - 0.5,
                             chi    = (A.lag[ii,mm] - A.prior[ii,mm])^2,
                             psi    = A_tau[ss+1,1]*prod(lambda2_A[2:(ss+1),1]))
            temp <- ifelse(temp<1e-7,1e-7,ifelse(temp>1e+7,1e+7,temp))
            theta.lag[ii,mm] <- temp
          }
        }
        theta[slct.i,] <- theta.lag
        if(sample_tau){
          #Sample a_tau through a simple RWMH step (on-line tuning of the MH scaling within the first 50% of the burn-in phase)
          A_tau_prop      = exp(rnorm(1,0,A_tuning[ss+1,1]))*A_tau[ss+1,1]
          A_tau_prop      = ifelse(A_tau_prop<1e-16,1e-16,A_tau_prop)
          A_tau_prop      = ifelse(A_tau_prop>1e+16,1e+16,A_tau_prop)
          post_A_tau_prop = .atau_post(atau=A_tau_prop,    thetas=as.vector(theta.lag),lambda2=prod(lambda2_A[2:(ss+1),1]))
          post_A_tau_old  = .atau_post(atau=A_tau[ss+1,1], thetas=as.vector(theta.lag),lambda2=prod(lambda2_A[2:(ss+1),1]))
          post.diff       = post_A_tau_prop-post_A_tau_old+log(A_tau_prop)-log(A_tau[ss+1,1])
          post.diff       = ifelse(is.nan(post.diff),-Inf,post.diff)
          if (post.diff > log(runif(1,0,1))){
            A_tau[ss+1,1]    = A_tau_prop
            A_accept[ss+1,1] = A_accept[ss+1,1]+1
          }
          if (irep<(0.5*burnin)){
            if ((A_accept[ss+1,1]/irep)>0.3)  A_tuning[ss+1,1] <- 1.01*A_tuning[ss+1,1]
            if ((A_accept[ss+1,1]/irep)<0.15) A_tuning[ss+1,1] <- 0.99*A_tuning[ss+1,1]
          }
        }
      }
    }
    if(prior == 4){
      # local shrinkage parameters - L
      lambda_L = 1 / rgamma(n     = v, 
                            shape = 1, 
                            rate  = 1 / nu_L + 0.5 * as.vector(L_draw[lower.tri(L_draw)])^2 / tau_L)
      nu_L = 1 /rgamma(n     = v, 
                       shape = 1,
                       rate  = 1 + 1 / lambda_L)
      # global shrinkage parameter - L
      WSSR_L = sum(as.vector(L_draw[lower.tri(L_draw)])^2 / lambda_L)
      tau_L  = 1 / rgamma(n     = 1, 
                          shape = (v + 1)/2, 
                          rate  = 1 / zeta_L + 0.5*WSSR_L)
      zeta_L = 1 / rgamma(n     = 1, 
                          shape = 1, 
                          rate = 1 + 1 / tau_L)
      # update prior VCV
      L_prior[lower.tri(L_prior)] = tau_L * lambda_L
      
      ############# - A endo
      slct.i    <- grep("Ylag",rownames(A_draw))
      # local shrinkage parameter - A endo
      lambda_A_endo = 1 / rgamma(n     = n, 
                                 shape = 1,
                                 rate  = 1 / nu_A_endo + 0.5 * as.vector(A_draw[slct.i,])^2 / tau_A_endo)
      nu_A_endo     = 1 / rgamma(n     = n, 
                                 shape = 1,
                                 rate  = 1 + 1 / lambda_A_endo)
      # global shrinkage parameter - A endo
      WSSR_A_endo = sum(as.vector(A_draw[slct.i,])^2 / lambda_A_endo)
      tau_A_endo  = 1 / rgamma(n     = 1,
                               shape = (n + 1)/2,
                               rate  = 1 / zeta_A_endo + 0.5*WSSR_A_endo)
      zeta_A_endo = 1 / rgamma(n     = 1,
                               shape = 1,
                               rate  = 1 + 1 / tau_A_endo)
      # update prior VCV
      theta[slct.i,] <- tau_A_endo * lambda_A_endo
      
      ############# - A endo
      if(wexo){
        slct.w <- grep("Wex", rownames(A_draw))
        # local shrinkage parameter - A exo
        lambda_A_exo = 1 / rgamma(n     = nstar, 
                                  shape = 1,
                                  rate  = 1 / nu_A_exo + 0.5 * as.vector(A_draw[slct.w,])^2 / tau_A_exo)
        nu_A_exo     = 1 / rgamma(n     = nstar, 
                                  shape = 1,
                                  rate  = 1 + 1 / lambda_A_exo)
        # global shrinkage parameter - A exo
        WSSR_A_exo = sum(as.vector(A_draw[slct.w,])^2 / lambda_A_exo)
        tau_A_exo  = 1 / rgamma(n     = 1,
                                shape = (nstar + 1)/2,
                                rate  = 1 / zeta_A_endo + 0.5*WSSR_A_exo)
        zeta_A_exo = 1 / rgamma(n     = 1,
                                shape = 1,
                                rate  = 1 + 1 / tau_A_exo)
        # update prior VCV
        theta[slct.w,] <- tau_A_exo * lambda_A_exo
      }
    } 
    #----------------------------------------------------------------------------
    # Step 3: Sample variances
    if(sv){
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
      A_tmp <- matrix(0, k_end, M)
      if(wexo){
        A_tmp <- A_draw
      }else{
        A_tmp[1:K,] = A_draw[1:K,]
        if(cons) A_tmp[K+Kstar+1,] = A_draw[K+1,]
        if(trend) A_tmp[K+Kstar+ifelse(cons,1,0)+1,] = A_draw[K+ifelse(cons,1,0)+1,]
      }
      A_store[,,count]   = A_tmp
      L_store[,,count]   = L_draw
      res_store[,,count] = Y-X%*%A_draw
      # SV
      Sv_store[,,count] = Sv_draw
      if(save_vola_pars){
        pars_store[,,count] = pars_var
      }
      # MN
      if(save_shrink_MN){
        lambda_store[,,count] = c(lambda1,lambda2,lambda3)
      }
      # SSVS
      if(save_shrink_SSVS){
        gamma_tmp = matrix(0, k_end, 1)
        if(wexo){
          gamma_tmp = gamma
        }else{
          gamma_tmp[1:K,] = gamma[1:K,]
          if(cons) gamma_tmp[K+Kstar+1,] = gamma[K+1,]
          if(trend) gamma_tmp[K+Kstar+ifelse(cons,1,0)+1,] = gamma[K+ifelse(cons,1,0)+1,]
        }
        gamma_store[,,count] = gamma_tmp
        omega_store[,,count] = omega
      }
      # NG
      if(save_shrink_NG){
        theta_store[,,count]                = theta
        lambda2_store[1,3,count]            = lambda2_L
        lambda2_store[1:(plag+1),1:2,count] = lambda2_A
        tau_store[1,3,count]                = L_tau
        tau_store[1:(plag+1),1:2,count]     = A_tau
      }
      # HS
      if(save_shrink_HS){
        lambda_A_endo_store[,count] = lambda_A_endo
        lambda_A_exo_store[,count]  = lambda_A_exo
        lambda_L_store[,count]      = lambda_L
        nu_A_endo_store[,count]     = nu_A_endo
        nu_A_exo_store[,count]      = nu_A_exo
        nu_L_store[,count]          = nu_L
        tau_A_endo_store[,count]    = tau_A_endo
        tau_A_exo_store[,count]     = tau_A_exo
        tau_L_store[,count]         = tau_L
        zeta_A_endo_store[,count]   = zeta_A_endo
        zeta_A_exo_store[,count]    = zeta_A_exo
        zeta_L_store[,count]        = zeta_L
      }
    }
  }
  #---------------------------------------------------------------------------------------------------------
  # END ESTIMATION
  #---------------------------------------------------------------------------------------------------------
  dimnames(A_store)=list(nameslags_end,colnames(A_OLS),NULL)
  
  ret <- list(Y=Y,X=X,A_store=A_store,L_store=L_store,Sv_store=Sv_store,pars_store=pars_store,res_store=res_store,
              MN=list(
                lambda_store=lambda_store
                ),
              SSVS=list(
                gamma_store=gamma_store,
                omega_store=omega_store),
              NG=list(
                theta_store=theta_store,
                lambda2_store=lambda2_store,
                tau_store=tau_store),
              HS=list(
                lambda_A_endo_store=lambda_A_endo_store,
                lambda_A_exo_store=lambda_A_exo_store,
                lambda_L_store=lambda_L_store,
                nu_A_endo_store=nu_A_endo_store,
                nu_A_exo_store=nu_A_exo_store,
                nu_L_store=nu_L_store,
                tau_A_endo_store=tau_A_endo_store,
                tau_A_exo_store=tau_A_exo_store,
                tau_L_store=tau_L_store,
                zeta_A_endo_store=zeta_A_endo_store,
                zeta_A_exo_store=zeta_A_exo_store,
                zeta_L_store=zeta_L_store
              ))
  
  return(ret)
}

#' @name .gvar.stacking.wrapper
#' @importFrom stats median
#' @importFrom utils memory.limit
#' @noRd
.gvar.stacking.wrapper<-function(xglobal,plag,globalpost,draws,thin,trend,eigen,trim,verbose){
  results <- tryCatch(
    {
      bigT      <- nrow(xglobal)
      bigK      <- ncol(xglobal)
      cN        <- names(globalpost)
      thindraws <- draws/thin
      F_large   <- array(NA, dim=c(bigK,bigK,plag,thindraws))
      trim.info <- "No trimming"
      
      ## call Rcpp
      # Rcpp::sourceCpp("./src/gvar_stacking.cpp")
      out <- gvar_stacking(xglobal    = xglobal,
                           plag       = as.integer(plag), 
                           globalpost = globalpost, 
                           draws      = as.integer(draws),
                           thin       = as.integer(thin), 
                           trend      = trend, 
                           eigen      = TRUE, 
                           verbose    = verbose)
      A_large    <- out$A_large
      for(pp in 1:plag){
        F_large[,,pp,] <- out$F_large[,((bigK*(pp-1))+1):(bigK*pp),,drop=FALSE]
      }
      S_large    <- out$S_large
      Ginv_large <- out$Ginv_large
      F.eigen    <- out$F_eigen
      dimnames(S_large)[[1]]<-dimnames(S_large)[[2]]<-dimnames(Ginv_large)[[1]]<-dimnames(Ginv_large)[[2]]<-dimnames(A_large)[[1]]<-colnames(xglobal)
      names <- c(paste(rep(colnames(xglobal),plag),".",rep(seq(1,plag),each=bigK),sep=""),"cons")
      if(trend) names <- c(names,"trend")
      dimnames(A_large)[[2]]<-names
      
      # kick out in-stable draws
      if(eigen){
        idx<-which(F.eigen<trim)
        
        F_large     <- F_large[,,,idx,drop=FALSE]
        S_large     <- S_large[,,idx,drop=FALSE]
        Ginv_large  <- Ginv_large[,,idx,drop=FALSE]
        A_large     <- A_large[,,idx,drop=FALSE]
        F.eigen     <- F.eigen[idx]
        
        if(length(idx)<10){
          stop("Less than 10 stable draws have been found. Please re-estimate the model.")
        }
        
        trim.info <- round((length(idx)/thindraws)*100,2)
        trim.info <- paste("Trimming leads to ",length(idx) ," (",trim.info,"%) stable draws out of ",thindraws," total draws.",sep="")
      }
      
      results<-list(S_large=S_large,F_large=F_large,Ginv_large=Ginv_large,A_large=A_large,F.eigen=F.eigen,trim.info=trim.info)
    },
    error=function(cond){
      if(cond$message == "vector memory exhausted (limit reached?)"){
        message("BGVAR incurred an error due to memory exhaustiveness.\n See original error message:")
        message(cond)
        if(.get_os() == "osx"){
          message("\nWe advise you to do the following on Mac OS: The enviornment variable R_MAX_VSIZE can be adjusted to to specify the maximal vector heap size.")
          message("\nCurrent value of R_MAX_VSIZE: ", Sys.getenv("R_MAX_VSIZE"))
          message("\nAdjust this parameter as follows: Sys.setenv(R_MAX_VSIZE='100GB')")
        }else if(.get_os() == "windows"){
          message("\nWe advise you to do the following on Windows OS: The memory limit can be adjusted to specify the maximal vector heap size.")
          message("\nCurrent value of memory.limit(): ", memory.limit())
          message("\nAdjust this parameter as follows: memory.limit(size=100000).")
        }else if(.get_os() == "linux"){
          message("\nWe advise you to do the following on Linux: The memory limit can be adjusted to specify the maximal vector heap size. The maximal vector heap size of physical and virtual memory is set to the maximum of 16Gb. \n Adjust this parameter as follows: unix::rlimit_as(Inf).")
        }
        message("\n In any case, increasin the thinning factor (argument 'thin' of 'bgvar') reduces memory requirements.")
      }
      return(NULL)
    },
    warning=function(cond){},
    finally={}
  )
  return(results)
}

#' @name .gvar.stacking
#' @importFrom abind adrop
#' @importFrom Matrix bdiag
#' @importFrom stats median
#' @importFrom utils txtProgressBar setTxtProgressBar
#' @noRd
.gvar.stacking<-function(xglobal,plag,globalpost,draws,thin,trend,eigen=FALSE,trim=NULL,verbose=TRUE){
  # initialize objects here
  bigT <- nrow(xglobal) 
  bigK <- ncol(xglobal)
  cN   <- names(globalpost)
  
  thindraws <- draws/thin
  F.eigen   <- numeric(thindraws)
  trim.info <- "No trimming"
  
  A_large     <- array(NA_real_, dim=c(bigK,bigK*plag+1+ifelse(trend,1,0),thindraws))
  S_large     <- array(NA_real_, dim=c(bigK,bigK,thindraws))
  Ginv_large  <- array(NA_real_, dim=c(bigK,bigK,thindraws))
  F_large     <- array(NA_real_, dim=c(bigK,bigK,plag,thindraws))
  dimnames(S_large)[[1]]<-dimnames(S_large)[[2]]<-dimnames(Ginv_large)[[1]]<-dimnames(Ginv_large)[[2]]<-dimnames(A_large)[[1]]<-colnames(xglobal)
  
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
      A   <- cbind(diag(ncol(VAR$Y)),-t(adrop(VAR$store$Lambda0store[,,irep,drop=FALSE],drop=3)))
      
      for(pp in 1:plag){
        assign(paste("B",pp,sep=""),cbind(t(adrop(VAR$store$Phistore[[pp]][,,irep,drop=FALSE],drop=3)),
                                    t(adrop(VAR$store$Lambdastore[[pp]][,,irep,drop=FALSE],drop=3))))
        if(cc==1) assign(paste("H",pp,sep=""), get(paste("B",pp,sep=""))%*%W)
        if(cc>1)  assign(paste("H",pp,sep=""), rbind(get(paste("H",pp,sep="")),get(paste("B",pp,sep=""))%*%W))
      }
      G            <- rbind(G,A%*%W)
      a0           <- rbind(a0,VAR$store$a0store[,irep,drop=FALSE])
      if(trend) a1 <- rbind(a1,VAR$store$a1store[,irep,drop=FALSE])
      S_post[[cc]] <- adrop(VAR$store$SIGMAmed_store[,,irep,drop=FALSE],drop=3)
    }
    G.inv  <- solve(G)
    S_large[,,irep] <- as.matrix(bdiag(S_post))
    b0     <- G.inv%*%a0
    if(trend) b1 <- G.inv%*%a1 else b1 <- NULL
    Ginv_large[,,irep] <- G.inv
    
    ALPHA <- NULL
    for (kk in 1:plag){
      assign(paste("F",kk,sep=""),G.inv%*%get(paste("H",kk,sep="")))
      F_large[,,kk,irep] <- get(paste("F",kk,sep=""))
      ALPHA <- cbind(ALPHA,F_large[,,kk,irep])
    }
    
    ALPHA <- cbind(ALPHA,b0,b1)
    A_large[,,irep]<-ALPHA
    
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
    if(verbose) setTxtProgressBar(pb, irep)
  }
  
  # kick out in-stable draws
  if(!is.null(trim)){
    if(trim==TRUE) trim <- 1.05
    idx<-which(F.eigen<trim)
    
    F_large     <- F_large[,,,idx,drop=FALSE]
    S_large     <- S_large[,,idx,drop=FALSE]
    Ginv_large  <- Ginv_large[,,idx,drop=FALSE]
    A_large     <- A_large[,,idx,drop=FALSE]
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

#' @name .get_os
#' @noRd
.get_os <- function(){
  sysinf <- Sys.info()
  if (!is.null(sysinf)){
    os <- sysinf['sysname']
    if (os == 'Darwin')
      os <- "osx"
  } else { ## mystery machine
    os <- .Platform$OS.type
    if (grepl("^darwin", R.version$os))
      os <- "osx"
    if (grepl("linux-gnu", R.version$os))
      os <- "linux"
  }
  return(tolower(os))
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
    colNames[[cc]] <- gsub(paste(cN[cc],"\\.",sep=""),"",colNames[[cc]])
    varNames[[cc]] <- gsub(paste(cN[cc],"\\.",sep=""),"",varNames[[cc]])
  }
  colNames <- unique(unlist(colNames))
  varNames <- unique(unlist(varNames))
  shrink <- array(NA,dim=c(length(varNames),length(colNames),N))
  dimnames(shrink)[[1]] <- varNames; dimnames(shrink)[[2]] <- colNames; dimnames(shrink)[[3]] <- cN
  for(cc in 1:N){
    aux <- country.shrink[[cc]]
    rownames(aux)<-gsub(paste(cN[cc],"\\.",sep=""),"",rownames(aux))
    colnames(aux)<-gsub(paste(cN[cc],"\\.",sep=""),"",colnames(aux))
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
.irf.sign.zero <- function(xdat,plag,n.ahead,Ginv,Fmat,Smat,shocklist,...){
  bigT          <- nrow(xdat)
  bigK          <- ncol(xdat)
  varNames      <- colnames(xdat)
  shock.idx     <- shocklist$shock.idx
  shock.cidx    <- shocklist$shock.cidx
  MaxTries      <- shocklist$MaxTries
  S.cube        <- shocklist$S.cube
  P.cube        <- shocklist$P.cube
  Z.cube        <- shocklist$Z.cube
  shock.horz    <- shocklist$shock.horz
  shock.order   <- shocklist$shock.order
  H.restr       <- length(shock.horz)
  N.restr       <- bigK*H.restr
  N             <- length(shock.idx)
  no.zero.restr <- shocklist$no.zero.restr
  #-----------------------------------------------------------------------------
  # create P0G
  P0G <- diag(bigK); colnames(P0G) <- rownames(P0G) <- varNames
  for(cc in 1:N){
    idx <- shock.idx[[cc]]
    if(shock.cidx[cc]){
      temp <- try(t(chol(Smat[idx,idx,drop=FALSE])),silent=TRUE) 
      if(is(temp,"try-error")){
        return(list(impl=NA,rot=NA,icounter=NA))
      }
      P0G[idx,idx] <- temp
    }else{
      P0G[idx,idx] <- Smat[idx,idx,drop=FALSE]
    }
  }
  # create dynamic multiplier
  PHIx <- array(0,c(bigK,bigK,plag+n.ahead+1)); dimnames(PHIx)[[1]] <- dimnames(PHIx)[[2]] <- varNames
  PHIx[,,plag+1]  <-  diag(bigK)
  for (ihor in (plag+2):(plag+n.ahead+1)){
    acc = matrix(0,bigK,bigK)
    for (pp in 1:plag){
      acc  <-  acc + Fmat[,,pp]%*%PHIx[,,ihor-pp]
    }
    PHIx[,,ihor]  <-  acc
  }
  PHI  <-  PHIx[,,(plag+1):(plag+n.ahead+1)]
  #-----------------------------------------------------------------------------
  irf.restr         <- matrix(NA, N.restr, bigK)
  invGSigma_u       <- Ginv%*%P0G
  for(hh in 1:H.restr){
    # ARRW: Definition 1
    if(shock.horz[hh]!=Inf) irf.hh<-PHI[,,shock.horz[hh]]%*%invGSigma_u
    # ARRW: Definition 2
    #if(sign.horizon[hh]==Inf) irf.hh <- solve(A0-A0%*%Cm[1:M,]%*%do.call("rbind",rep(list(diag(M)),p)))
    irf.restr[((hh-1)*bigK+1):(bigK*hh),1:bigK] <- irf.hh
  }
  colnames(irf.restr) <- varNames
  rownames(irf.restr) <- paste(rep(varNames,H.restr),".",
                               rep(shock.horz,each=bigK),sep="")
  #-----------------------------------------------------------------------------
  # reorder - important!!
  Z.cube <- Z.cube[,,shock.order]
  # draw rotation matrix here
  icounter <- 0
  condall <- 0
  impresp<-Q_bar<-NA
  while(condall == 0 && icounter < MaxTries){
    Q <- diag(bigK)
    for(cc in 1:N){
      idx <- shock.idx[[cc]]
      Kidx <- length(idx)
      if(shock.cidx[cc]){
        randMat <- matrix(rnorm(Kidx^2),Kidx,Kidx)
        Qc <- matrix(0, Kidx, Kidx)
        if(no.zero.restr[cc]){
          QR <- qr(randMat)
          Qc <- qr.Q(QR)
        }else{
          for(kk in 1:Kidx){
            Z.temp <- Z.cube[,,idx[kk]]
            Z.temp <- Z.temp[rowSums(abs(Z.temp))!=0,,drop=F]
            if(nrow(Z.temp)==0){
              Z.temp <- matrix(0, 1, N.restr)
            }
            if(all(Z.temp==0) && kk>1){
              R <- c()
            }else{
              R <- Z.temp%*%irf.restr[,idx]
            }
            if(kk > 1){R <- rbind(R, t(Qc[,(1:(kk-1)), drop=FALSE]))}
            
            NU  <- Null(t(R))
            x_j <- randMat[,kk,drop=FALSE]
            
            q_j <- NU%*%(t(NU)%*%x_j/sqrt(as.numeric(crossprod(t(NU)%*%x_j))))
            Qc[,kk] <- q_j
          }
        }
        Q[idx,idx] <- Qc
      }
    }
    colnames(Q) <- varNames[shock.order]; rownames(Q) <- varNames
    Q_bar <- Q[,varNames]
    # Q_bar <- Q%*%diag(((diag(Q)>0)-(diag(Q)<0)))
    # check irf
    irf.check <- irf.restr%*%Q_bar
    colnames(irf.check) <- varNames
    rownames(irf.check) <- paste(rep(varNames,H.restr),".",rep(shock.horz,each=bigK),sep="")
    signCheck <- matrix(NA,bigK,1)
    for(ss in 1:bigK){
      STemp <- S.cube[,ss,drop=FALSE]
      if(sum(abs(STemp))==0){
        signCheck[ss,] <- TRUE
        next
      }
      PDiag <- diag(N.restr); diag(PDiag) <- sign(P.cube[,ss,drop=TRUE]>runif(N.restr))
      IrfCheckTemp <- sign(irf.check[,ss,drop = FALSE])
      signCheck[ss,] <- t(IrfCheckTemp)%*%PDiag%*%STemp==sum(abs(PDiag%*%STemp))
    }
    condall <- prod(signCheck)
    icounter <- icounter + 1
  }
  # compute impulses
  st_impulses <- array(NA,c(bigK,bigK,n.ahead+1));dimnames(st_impulses)[[1]] <- dimnames(st_impulses)[[2]] <- varNames
  Cmhat <- invGSigma_u%*%Q_bar
  for(ihor in 1:(n.ahead+1)){
    st_impulses[,,ihor] <- as.matrix((PHI[,,ihor]%*%Cmhat))
  }
  
  if(icounter==MaxTries){
    st_impulses <- Q_bar <- NA
  }
  # end rotation matrix loop ----------------------------------------------------------------------------
  return(list(impl=st_impulses,rot=Q_bar,icounter=icounter))
}

#' @name .irf.chol
#' @noRd
.irf.chol <- function(xdat,plag,n.ahead,Ginv,Fmat,Smat,shocklist,...){
  bigT       <- nrow(xdat)
  bigK       <- ncol(xdat)
  varNames   <- colnames(xdat)
  shock.idx  <- shocklist$shock.idx
  shock.cidx <- shocklist$shock.cidx
  N          <- length(shock.idx)
  # create P0G
  P0G <- diag(bigK); colnames(P0G) <- rownames(P0G) <- varNames
  for(cc in 1:N){
    idx <- shock.idx[[cc]]
    if(shock.cidx[cc]){
      P0G[idx,idx] <- t(chol(Smat[idx,idx,drop=FALSE])) # calculate local cholesky factor of gcov
    }else{
      P0G[idx,idx] <- Smat[idx,idx,drop=FALSE]
    }
  }
  # create dynamic multiplier
  PHIx <- array(0,c(bigK,bigK,plag+n.ahead+1)); dimnames(PHIx)[[1]] <- dimnames(PHIx)[[2]] <- varNames
  PHIx[,,plag+1] <- diag(bigK)
  for (ihor in (plag+2):(plag+n.ahead+1)){
    acc = matrix(0,bigK,bigK)
    for (pp in 1:plag){
      acc  <-  acc + Fmat[,,pp]%*%PHIx[,,ihor-pp]
    }
    PHIx[,,ihor]  <-  acc
  }
  PHI  <-  PHIx[,,(plag+1):(plag+n.ahead+1)]
  # compute shock
  invGSigma_u  <- Ginv%*%P0G
  # computing impulse response function
  irfa  <- array(0,c(bigK,bigK,n.ahead+1)); dimnames(irfa)[[2]] <- varNames
  for (ihor in 1:(n.ahead+1)){
    irfa[,,ihor] <- PHI[,,ihor]%*%invGSigma_u
  }
  # define output
  out <- list(impl=irfa,rot=NULL,icounter=1)
  return(out)
}

#' @name .irf.girf
#' @noRd
.irf.girf <- function(xdat,plag,n.ahead,Ginv,Fmat,Smat, ...){
  bigT     <- nrow(xdat)
  bigK     <- ncol(xdat)
  varNames <- colnames(xdat)
  # create dynamic multiplier
  PHIx <- array(0,c(bigK,bigK,plag+n.ahead+1)); dimnames(PHIx)[[1]] <- dimnames(PHIx)[[2]] <- varNames
  PHIx[,,plag+1]  <-  diag(bigK)
  for (ihor in (plag+2):(plag+n.ahead+1)){
    acc = matrix(0,bigK,bigK)
    for (pp in 1:plag){
      acc  <-  acc + Fmat[,,pp]%*%PHIx[,,ihor-pp]
    }
    PHIx[,,ihor]  <-  acc
  }
  PHI  <-  PHIx[,,(plag+1):(plag+n.ahead+1)]
  # create shock
  invGSigma_u <- Ginv%*%Smat
  # computing impulse response function
  irfa  <- array(0,c(bigK,bigK,n.ahead+1)); dimnames(irfa)[[1]]  <- varNames
  for (ihor in 1:(n.ahead+1)){
    irfa[,,ihor] <- PHI[,,ihor]%*%invGSigma_u
  }
  # define output
  out <- list(impl=irfa,rot=NULL,icounter=1)
  return(out)
}

#' @name .irf.girf.sims
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
.impulsdtrf <- function(B,smat,nstep){
  
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
