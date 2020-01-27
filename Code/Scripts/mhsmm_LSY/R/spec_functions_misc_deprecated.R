
# works but only with diagonal sigmas
dmvnorm.hsmm.censored.diag <- function(x,c, j, model) {
  ans <- dmvnorm_censored_diag(x, c, mean = model$parms.emission$mu[[j]],
                               sigma = model$parms.emission$sigma[[j]])
  ans[is.na(ans)] <- 1
  ans 
}


# works but is very slow
dmvnorm.hsmm.censored_v2 <- function(x, c, j, model) {
  c = (c==1)
  ans = sapply(1:nrow(x), function(i) dmvnorm.hsmm.censored.row(i,x,c,j,model))
  ans[is.na(ans)] <- 1
  ans 
}

# does not do what I want
dmvnorm.hsmm.censored_v3 <- function(x, c, j, model) {
  c = (c==1)
  mus_obs = matrix(model$parms.emission$mu[[j]], ncol = ncol(x), nrow = nrow(x), byrow = TRUE)
  xx = x
  xx[!c] = mus_obs[!c]
  ans = dmvnorm(xx, mean = model$parms.emission$mu[[j]],
                sigma = model$parms.emission$sigma[[j]])
  ans[is.na(ans)] <- 1
  ans 
}




dmvnorm_censored_diag = function (x, c, mean = rep(0, p), sigma = diag(p), log = FALSE) 
{
  if (is.vector(x)) 
    x <- matrix(x, ncol = length(x))
  p <- ncol(x)
  if (!missing(mean)) {
    if (!is.null(dim(mean))) 
      dim(mean) <- NULL
    if (length(mean) != p) 
      stop("mean and sigma have non-conforming size")
  }
  if (!missing(sigma)) {
    if (p != ncol(sigma)) 
      stop("x and sigma have non-conforming size")
    if (!isSymmetric(sigma, tol = sqrt(.Machine$double.eps), 
                     check.attributes = FALSE)) 
      stop("sigma must be a symmetric matrix")
    # we further check that the co-variance are NULL
    if(!(all(sigma[upper.tri(sigma)] == 0)&all(sigma[lower.tri(sigma)] == 0))){
      stop("sigma must be a diagonal matrix")
    }
  }
  dec <- tryCatch(chol(sigma), error = function(e) e)
  if (inherits(dec, "error")) {
    x.is.mu <- colSums(t(x) != mean) == 0
    logretval <- rep.int(-Inf, nrow(x))
    logretval[x.is.mu] <- Inf
  }
  else {
    DEC = matrix(diag(dec),ncol(x), nrow(x))
    DEC[t(c) == 0] = NA
    #tmp <- backsolve(dec, t(x) - mean, transpose = TRUE)
    tmp = (t(x) - mean) / DEC
    #rss <- colSums(tmp^2)
    rss <- colSums(tmp^2, na.rm = TRUE)
    #logretval <- -sum(log(diag(dec))) - 0.5 * p * log(2 * 
    #                                                    pi) - 0.5 * rss
    logretval <- -colSums(log(DEC),na.rm = TRUE) - 0.5 * rowSums(c) * log(2 * pi) - 0.5 * rss
  }
  names(logretval) <- rownames(x)
  if (log) 
    logretval
  else exp(logretval)
}



dmvnorm.hsmm.censored.row = function(i, x, c,j, model){
  l = which(c[i,])
  ans = dmvnorm(x[i,l], mean = model$parms.emission$mu[[j]][l],
                sigma = as.matrix(model$parms.emission$sigma[[j]][l,l])
  )
  return(ans)
}


