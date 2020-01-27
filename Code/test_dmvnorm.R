

##########

my_dmvnorm = function (x, c, mean = rep(0, p), sigma = diag(p), log = FALSE) 
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


dmvnorm_row = function(i, x, c, mu, sigma){
  j = which(c[i,] == 1)
  d = dmvnorm(x[i,j], mean = mu[j], sigma = as.matrix(sigma[j,j]))
}


###############


n = 10000
mu = c(1:3)
sigma = diag(rep(1,3),3,3)

x = rmvnorm(n, mean = mu, sigma = sigma)
c = matrix(sample(c(0,1), n*3, replace = TRUE), n, 3)
c[which(apply(c,1,sum) == 0),1] = 1
I = matrix(1, n , 3)
xNA = x
xNA[c==0] = NA


tic()
ans1 = dmvnorm(x, mean = mu, sigma = sigma)
toc()


tic()
ans2 = sapply(1:n, function(i) dmvnorm_row(i,x,I,mu, sigma))
toc()

all(ans1 == ans2)


tic()
ans3 = sapply(1:n, function(i) dmvnorm_row(i,x,c,mu, sigma))
toc()


tic()
ans4 = dmvnorm(xNA, mean = mu, sigma = sigma)
toc()

sum(is.na(ans4))
all(ans3 == ans4)
all(ans3 == ans4, na.rm = TRUE)


tic()
ans5 = my_dmvnorm(x, c, mean = mu, sigma = sigma, log = FALSE)
toc()

cat(all(ans3 == ans5), "\n")
plot(ans3, ans5)

