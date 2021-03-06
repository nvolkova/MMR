#########################################################################
### Functions for NMF: extracting signatures, contributions,
### extracting new signatures simultaneously, Qoasi-Poisson models
### M. Gerstung, N. Volkova, EMBL-EBI 2016-2017
#########################################################################

# Extraction of signatures for given exposures
# D - matrix of mutation counts, samples x mutation types
# P - exposures, samples x signatures
# maxIter - maximal number of iterations in NMF algorithm
# tol - lower bound for the change of signature matrix between iterations
# div.err - lower bound for the change of the divergence between D and its factorization between iterations
nmSolve <- function(D, P, maxIter = 10000, tol=1e-5, div.err=1e-7) {
  n <- nrow(D)
  mask <- !is.na(D)
  m <- ncol(D)
  s <- ncol(P)
  rP <- rep(colSums(P), m)
  tP <- t(P)
  D <- as.matrix(D)
  P <- as.matrix(P)
  E1 <- E2 <- matrix(runif(s * m, 1e-3, 1), ncol = m)
  err <- 2*tol
  D[is.na(D)] <- 0
  
  iter <- 1
  divergence.old <- mean(D*log(D/(P %*% (E2 + .Machine$double.eps))) - D + P%*%E2, na.rm=T)
  div.change <- 2 * div.err
  
  while (iter < maxIter & err > tol & abs(div.change) > div.err) {
    E1 <- E2
    E2 <- E1 * (tP %*% ((mask*D)/(mask*(P %*% (E1 + .Machine$double.eps)) + .Machine$double.eps)))/rP
    iter <- iter + 1
    err <- mean(abs(E2 - E1)/(E1+.Machine$double.eps), na.rm=TRUE)
    divergence <- mean(D*log(D/(P %*% (E2 + .Machine$double.eps))) - D + P%*%E2, na.rm=T) # KL distance from D to P%*%E2
    div.change <- divergence.old - divergence
    divergence.old = divergence
    if(iter %% 100 == 0) cat(round(-log10(err)))
  }
  cat("\n")
  if(iter == maxIter) warning(paste("No convergence after",iter, "iterations."))
  E2
}

# Fitting the signatures (extracting exposures)
# D - matrix of mutation counts, samples x mutation types
# E - signatures, signatures x mutation types
# maxIter - maximal number of iterations in NMF algorithm
# tol - lower bound for the change of exposure matrix between iterations
# div.err - lower bound for the change of the divergence between D and its factorization between iterations
nmFit <- function(D, E, maxIter = 10000, tol=1e-5, div.err=1e-7) {
  n <- nrow(D)
  m <- ncol(D)
  s <- nrow(E)
  tE <- t(E)
  rE <- rep(rowSums(E),each=n)
  D <- as.matrix(D)
  E <- as.matrix(E)
  P1 <- P2 <- matrix(runif(s * n, 1e-3, 1), nrow = n)
  err <- 2*tol
  
  iter <- 1
  divergence.old <- mean(D * log(D/(P2%*%(E + .Machine$double.eps))) - D + P2%*%E, na.rm=T)
  div.change <- 2 * div.err
  while (iter < maxIter & err > tol & abs(div.change) > div.err) {
    P1 <- P2
    P2 <- P1 * ((D/((P1 +.Machine$double.eps) %*% (E +.Machine$double.eps))) %*% tE) / (rE+.Machine$double.eps)
    iter <- iter + 1
    err <- mean(abs(P2 - P1)/(P1+.Machine$double.eps), na.rm=TRUE)
    divergence <- mean(D*log(D/((P2 + .Machine$double.eps)%*%(E+.Machine$double.eps))) - D + P2%*%E, na.rm=T) # KL distance from D to P%*%E
    div.change <- divergence.old - divergence
    divergence.old = divergence
    if(iter %% 100 == 0) cat(round(-log10(err)))
    # add likelihood convergence
  }
  cat("\n")
  if(iter == maxIter) warning(paste("No convergence after",iter, "iterations."))
  P2
}