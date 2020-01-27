
# POISSON

mstep.pois <- function (x, wt) 
{
  k = ncol(wt)
  lambda = numeric(k)
  for (i in 1:k) lambda[i]=weighted.mean(x,wt[,i])
  list(lambda=lambda)
}

dpois.hsmm <- function (x, j, model) dpois(x,model$parms.emission$lambda[j])

rpois.hsmm <- function (j, model)  rpois(1, model$parms.emission$lambda[j])


# NORMAL (UNIVARIATE)

mstep.norm <- function(x,wt) {
  k = ncol(wt)
  mu = numeric(k)
  sigma = numeric(k)
  for(i in 1:k) {
    tmp = cov.wt(data.frame(x[!is.na(x)]),wt[!is.na(x),i])
    mu[i] = tmp$center
    sigma[i] = tmp$cov
  }
  list(mu=mu,sigma=sigma)
}

dnorm.hsmm <- function(x,j,model) {
  ret = dnorm(x,model$parms.emission$mu[j],sqrt(model$parms.emission$sigma[j]))
  ret[is.na(ret)] = 1
  ret           
}

rnorm.hsmm <- function(j,model)  rnorm(1,model$parms.emission$mu[j],sqrt(model$parms.emission$sigma[j]))


# MULTIVARIATE NORMAL (NOT CENSORED)

rmvnorm.hsmm <- function(j,model) 
  rmvnorm(1,mean=model$parms.emission$mu[[j]],sigma=model$parms.emission$sigma[[j]])

mstep.mvnorm <- function(x, wt) {
  idx <-  apply(is.na(x),1,any) # Find rows with NA's (cov.wt does not like them)
  x  <- x[!idx,,drop=FALSE]
  wt <- wt[!idx,,drop=FALSE]
  emission <- list(mu = list(), sigma = list())
  for (i in 1:ncol(wt)) {  
    tmp <- cov.wt(x, wt[, i])
    emission$mu[[i]] <- tmp$center
    emission$sigma[[i]] <- tmp$cov
  }
  emission
} 

dmvnorm.hsmm <- function(x, j, model) {
  ans <- dmvnorm(x, mean = model$parms.emission$mu[[j]],
                 sigma = model$parms.emission$sigma[[j]])
  ans[is.na(ans)] <- 1
  ans 
}


# MULTIVARIATE NORMAL WITH CENSORING 

rmvnorm.hsmm.censored <- function(j,model) # This one is the same
  rmvnorm(1,mean=model$parms.emission$mu[[j]],sigma=model$parms.emission$sigma[[j]])

mstep.mvnorm.censored <- function(x, wt) {
  emission <- list(mu = list(), sigma = list())
  for (i in 1:ncol(wt)) {  
    tmp <- cov.wt.censored(x, wt[, i])
    emission$mu[[i]] <- tmp$center
    emission$sigma[[i]] <- tmp$cov
  }
  emission
} 

dmvnorm.hsmm.censored <- function(x, c, j, model, w = NULL) {
  N = ifelse(is.matrix(x),ncol(x),length(x))
  if(is.null(w)){w = rep(1/N,N)}
  if(length(w)!=N){error("w must be the same length as the number of variables in x")}
  # we replace the missing value by the most likely observations in that state
  if(sum(c == 0)>0){
    mu_obs = matrix(model$parms.emission$mu[[j]], ncol = ncol(x), nrow = nrow(x), byrow = TRUE)
    x[c == 0] = mu_obs[c == 0]
  }
  ans = dmvnorm(x, mean = model$parms.emission$mu[[j]],
                sigma = model$parms.emission$sigma[[j]])
  ans[is.na(ans)] <- 1
  ans = ans* rowSums(c*w) # we multiply by a factor proportional to the weighted number of missing variables
  ans 
}


##############

cov.wt.censored = function(x, wt = rep(1/nrow(x), nrow(x)), cor = FALSE, center = TRUE){
  if (is.data.frame(x)) x <- as.matrix(x)
  else if (!is.matrix(x)) stop("'x' must be a matrix or a data frame")
  #if (!all(is.finite(x))) warning("'x' contains not-finite values")
  n <- nrow(x)
  if (with.wt <- !missing(wt)) {
    if (length(wt) != n) 
      stop("length of 'wt' must equal the number of rows in 'x'")
    if (any(wt < 0) || (s <- sum(wt)) == 0) 
      stop("weights must be non-negative and not all zero")
    wt <- wt/s
  }
  wt_mat = matrix(wt,ncol = ncol(x), nrow = n, byrow = FALSE) * (!is.na(x))
  wt_mat = sweep(wt_mat, 2, colSums(wt_mat), FUN = "/", check.margin = FALSE)
  
  if (is.logical(center)) {
    center <- if (center) colSums(wt_mat * x, na.rm = TRUE) else 0 #colSums(wt * x, na.rm = TRUE)
  }
  else {
    if (length(center) != ncol(x)) 
      stop("length of 'center' must equal the number of columns in 'x'")
  }
  x <- sqrt(wt_mat) * sweep(x, 2, center, check.margin = FALSE) #x <- sqrt(wt) * sweep(x, 2, center, check.margin = FALSE)
  
  #cov = var(x, na.rm = TRUE, use = "pairwise.complete.obs") * nrow(x) /(1 - sum(wt^2))
  cov = var(x, na.rm = TRUE, use = "pairwise.complete.obs") * nrow(x) / (1 -  t(wt_mat) %*% wt_mat)
  
  y <- list(cov = cov, center = center, n.obs = n)
  if (with.wt) 
    y$wt <- wt; y$wt_mat = wt_mat
  if (cor) {
    Is <- 1/sqrt(diag(cov))
    R <- cov
    R[] <- Is * cov * rep(Is, each = nrow(cov))
    y$cor <- R
  }
  y
}




