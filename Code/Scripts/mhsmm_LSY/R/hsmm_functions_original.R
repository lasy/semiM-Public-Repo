
predict.hsmm_original <- function(object,newdata,method="viterbi",...) {
  if(missing(newdata)) stop("newdata missing!")
  else x=newdata
  J = object$J
  m=-1e300
  if(mode(x)=="numeric" | mode(x)=="integer") {
    warning('x is a primitive vector.  Assuming single sequence.')
    x0 = x
    N = NROW(x)
    NN = c(0,N)
    if(N<1) stop("N less than one")
  }
  else{
    N = x$N
    NN = cumsum(c(0,x$N))
    x0 = x$x
  }
  #statehat = integer(NROW(x))
  statehat = NA
  if(method=="viterbi") {
    M = nrow(object$model$d)    
    loga = as.double(log(object$model$transition))
    loga[loga==-Inf]=m
    logstart = as.double(log(object$model$init))
    logstart[logstart==-Inf|is.nan(logstart)] = m
    d = apply(object$model$d,2,function(x) x/sum(x))
    D = apply(d,2,function(x) rev(cumsum(rev(x))))
    d = log(d)
    d[d==-Inf]=m
    D = log(D)
    D[D==-Inf]=m
    loglik=0
    for(i in 1:length(N)) {
      if(NCOL(x0)==1)    b = log(unlist(sapply(1:J,function(state) object$f(x0[(NN[i]+1):NN[i+1]],state,object$model))))
      else    b = log(unlist(sapply(1:J,function(state) object$f(x0[(NN[i]+1):NN[i+1],],state,object$model))))   
      b[b==-Inf]=m
      tmp = .C("viterbi",
               a=loga,
               pi=logstart,
               p=as.double(b),
               d=as.double(d),
               D=as.double(D),
               timelength=as.integer(N[i]),
               J=as.integer(J), 
               M=as.integer(rep(M,J)),
               alpha = double(N[i]*J),
               statehat=integer(N[i]),
               psi_state0=integer(N[i]*J),
               psi_time0=integer(N[i]*J)          
               ,PACKAGE='mhsmm')
      print(range(b))
      print(head(    tmp$alpha ))
      print(tail(    tmp$alpha ))
      loglik=loglik+max(tmp$alpha[N[i]*(1:J)])
      statehat[(NN[i]+1):NN[i+1]] = tmp$statehat+1
    }
    ans <- list(x=x,s=statehat,N=N,loglik=loglik)
  }
  else if(method=="smoothed") {
    M = nrow(object$model$d)    
    m <- object$model
    m$dens.emission <- object$f
    tmp <- hsmmfit(x,m,object$mstep,maxit=1,M=M)
    ans <- list(x=x$x,s=tmp$yhat,N=x$N,p=matrix(tmp$estep_variables$gamma,ncol=object$J))
  }
  else stop(paste("Unavailable prediction method",method))
  
  class(ans) <- 'hsmm.data'  
  ans
  # tmp
}






hsmmfit.original <- function(x,model,mstep=NULL,M=NA,maxit=100,lock.transition=FALSE,lock.d=FALSE,graphical=FALSE, debug = FALSE) {
  sojourn.distribution=model$sojourn$type
  tol=1e-4
  ksmooth.thresh = 1e-20 #this is a threshold for which d(u) values to use - if we throw too many weights in the default density() seems to work quite poorly
  shiftthresh = 1e-20 #threshold for effective "0" when considering d(u)
  J = nrow(model$transition)
  model$J = J
  
  if(is.null(mstep)) 
    if(is.null(model$mstep)) stop("mstep not specified")
  else  mstep=model$mstep      
  
  .check.hsmmspec(model)
  
  f=model$dens.emission
  
  if(mode(x)=="numeric" | mode(x)=="integer") {
    warning('x is a primitive vector.  Assuming single sequence.')
    NN = NROW(x)    
  }
  else{
    NN = x$N
    x = x$x
  }
  if(is.na(M)) M = max(NN)      
  if(length(model$init)!=J) stop("length(model$init)!=J")
  if(NROW(x)!=sum(NN)) stop("NROW(x)!=sum(NN)")
  model <- .build_d(model,M)
  
  new.model = model
  ll = rep(NA,maxit)
  rm(model)
  
  
  for(it in 1:maxit) {
    if(debug){cat("it: ",it,"\n")}
    if(graphical)   plot.hsmm(list(model=new.model,J=J))
    p = sapply(1:J,function(state) f(x,state,new.model))
    print(range(p))
    if(any(is.na(p)|p==Inf)) stop("NAs detected in b(x), check your supplied density function")
    if(any(apply(p,1,max)==0)) stop("Some values have 0 pdf for all states!  Check your model parameters")
    if(graphical) matplot(p, type = "l", col = hsmm$states$colors)
    
    #    print(paste("Iteration",it))
    # E-STEP
    estep_variables  = .C("backward",
                          transition=as.double(new.model$transition),
                          init=as.double(new.model$init),
                          p=as.double(p),
                          d=as.double(new.model$d),
                          D=as.double(new.model$D),
                          timelength=as.integer(NN),
                          J=as.integer(J),
                          M=as.integer(rep(M,J)),
                          L1 = double(NROW(x)*J),N = double(NROW(x)),
                          eta = double(M*J),
                          F1=double(J*NROW(x)),
                          si=double(J*NROW(x)),
                          gamma=double(J*NROW(x)),
                          nsequences=as.integer(length(NN)),
                          totallength=NROW(x),
                          G=double(J*NROW(x)),
                          PACKAGE='mhsmm')
    
    #M-Step
    
    #    N.debug[[it]] = estep_variables$N
    if(any(is.nan(estep_variables$gamma))) {
      warning("NaNs detected in gamma.  Exiting...")
      return(estep_variables)
    }
    if(any(estep_variables$gamma<0)) estep_variables$gamma = zapsmall(estep_variables$gamma)      
    if(any(estep_variables$eta<0)) estep_variables$eta = zapsmall(estep_variables$eta)      
    if(any(estep_variables$N<0))  estep_variables$N = zapsmall(estep_variables$N)
    
    
    old.model = new.model
    state_wt <- matrix(estep_variables$gamma,ncol=J)
    #if(any(colSums(state_wt)==0)) stop("Error: at least one state has an expected number of occurences equal to 0.\n This may be caused by bad starting parameters are insufficent sample size")
    #new.model$parms.emission = mstep(x,state_wt_mod)
    # I replaced the two lines above with the following paragraph:
    if(debug){cat("colSums(state_wt) = ",colSums(state_wt),"\n")}
    if(any(colSums(state_wt)==0) & debug) cat("!!! at least one state has an expected number of occurences equal to 0.\n")
    state_wt_mod = state_wt 
    unvisited_states = which(colSums(state_wt) == 0)
    if(length(unvisited_states)>0){ state_wt_mod[,unvisited_states] = 1}
    old.parms.emission = new.model$parms.emission
    new.parms.emission = mstep(x,state_wt_mod)
    for(state in unvisited_states){for(par.n in length(new.parms.emission)){ new.parms.emission[[par.n]][[state]] = old.parms.emission[[par.n]][[state]]}}
    new.model$parms.emission = new.parms.emission
    
    if(lock.d) {
      cat("lock.d\n")
      new.model$d = old.model$d
      new.model$D = old.model$D
    }
    else {
      if(sojourn.distribution=="nonparametric") {
        new.model$d = apply(matrix(estep_variables$eta,ncol=J),2,function(x) x/sum(x))
        new.model$sojourn$d <- new.model$d
      }
      else if(sojourn.distribution=="ksmoothed-nonparametric") {
        new.model$d = apply(matrix(estep_variables$eta+1e-100,ncol=J),2,function(x) x/sum(x))    
        for(i in 1:J) {
          new.model$d[,i] = approx(density(which(new.model$d[,i]>ksmooth.thresh),weights=new.model$d[which(new.model$d[,i]>ksmooth.thresh),i],from=1,n=M),xout=1:M)$y
          new.model$d[is.na(new.model$d[,i]),i] = 0
          new.model$d[,i] = (new.model$d[,i]+1e-300)/sum(new.model$d[,i])
        }
        new.model$sojourn$d <- new.model$d                
      }
      
      else if(sojourn.distribution=="poisson") {
        new.model$d = apply(matrix(estep_variables$eta,ncol=J),2,function(x) x/sum(x))
        new.model$sojourn$lambda = numeric(J)
        new.model$sojourn$shift = numeric(J)            
        for(i in 1:J) {
          eta = new.model$d[,i]
          maxshift =  match(TRUE,eta>shiftthresh)
          Mtmp = tail(which(eta>shiftthresh),1)
          new.model$sojourn$shift[i] = which.max(sapply(1:maxshift, function(shift) .dpois.hsmm.sojourn(x = maxshift:Mtmp,lambda=((maxshift:Mtmp)-shift)%*%eta[maxshift:Mtmp],shift=shift,log=TRUE)%*%eta[maxshift:Mtmp]))
          new.model$sojourn$lambda[i] = ((new.model$sojourn$shift[i]:Mtmp)-new.model$sojourn$shift[i])%*%eta[new.model$sojourn$shift[i]:Mtmp]         
          new.model$d[,i] = .dpois.hsmm.sojourn(1:M,new.model$sojourn$lambda[i],new.model$sojourn$shift[i])
        }
      }
      else if(sojourn.distribution=="nbinom") {       
        new.model$d = matrix(nrow=M,ncol=J)
        new.model$sojourn$size = numeric(J)
        new.model$sojourn$shift = integer(J)
        new.model$sojourn$mu = numeric(J)                            
        new.model$sojourn$prob = numeric(J)                            
        eta = matrix(estep_variables$eta,ncol=J)
        for(i in 1:J) { 
          tmp = .fitnbinom(eta[,i])          
          new.model$sojourn$shift[i] = tmp[1]
          new.model$sojourn$size[i] =  tmp[2]
          new.model$sojourn$mu[i] =  tmp[3]
          new.model$sojourn$prob[i] =  tmp[4]
          new.model$d[,i] =  .dnbinom.hsmm.sojourn(1:M,new.model$sojourn$size[i],new.model$sojourn$prob[i],new.model$sojourn$shift[i])
        }
      }
      else if(sojourn.distribution=="gamma") {
        new.model$d = matrix(estep_variables$eta,ncol=J)
        new.model$sojourn$shape = numeric(J)
        new.model$sojourn$scale = numeric(J)
        for(i in 1:J) {           
          tmp = gammafit(1:M,wt=new.model$d[,i])
          new.model$sojourn$shape[i] = tmp$shape
          new.model$sojourn$scale[i] = tmp$scale
          new.model$d[,i] = dgamma(1:M,shape=tmp$shape,scale=tmp$scale)              
        }          
      }
      else if(sojourn.distribution=="logarithmic") {
        new.model$d = apply(matrix(estep_variables$eta+1e-100,ncol=J),2,function(x) x/sum(x))
        new.model$sojourn$shape = numeric(J)
        for(i in 1:J) {           
          new.model$sojourn$shape[i] = .logdistrfit(wt=new.model$d[,i])
          new.model$d[,i] = .dlog(1:M,new.model$sojourn$shape[i])
        }
      }
      else if(sojourn.distribution=="lnorm") {
        eta = matrix(estep_variables$eta,ncol=J)
        new.model$d = matrix(nrow=M,ncol=J)
        new.model$sojourn$meanlog = numeric(J)
        new.model$sojourn$s.dlog = numeric(J)
        for(i in 1:J) {           
          new.model$sojourn$meanlog[i] = weighted.mean(log(1:M),eta[,i])
          new.model$sojourn$s.dlog[i] = sqrt(cov.wt(data.frame(log(1:M)),eta[,i])$cov)
          new.model$d[,i] = dlnorm(1:M,new.model$sojourn$meanlog[i],new.model$sojourn$s.dlog[i])
          new.model$d[,i] = new.model$d[,i]/sum(new.model$d[,i])              
        }
      }
      else stop("Invalid sojourn distribution")         
      new.model$D = apply(new.model$d,2,function(x) rev(cumsum(rev(x))))
    }
    
    #estimate transition matrix
    if(lock.transition) {
      cat("lock.transition\n")
      new.model$init=old.model$init
      new.model$transition = old.model$transition
    }
    else {
      new.model$init=estep_variables$init
      new.model$init[new.model$init<0]=0
      new.model$transition = matrix(estep_variables$transition,ncol=J)
      new.model$transition[new.model$transition<0]=0
    }
    
    ll[it]=sum(log(estep_variables$N))
    new.model$J = J
    if(it>2) if(abs(ll[it]-ll[it-1])<tol) break()
    
  } #end iterations of EM algorithm
  
  class(new.model) <- "hsmmspec"
  ret = list(loglik=ll[!is.na(ll)],
             model=new.model,
             estep_variables=estep_variables,
             M=M,
             J=J,
             NN=NN,
             f=f,
             mstep=mstep,
             yhat=apply(matrix(estep_variables$gamma,ncol=J),1,which.max))
  
  class(ret) <- "hsmm"
  ret
}


