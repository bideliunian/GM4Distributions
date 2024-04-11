#########################
## generate gaussian data with
## sample size: n
## mean: mu 
## precision matrix: Omega
######################
rMVNormP = function(n, mu, Omega){
  p = length(mu)
  Z = matrix(rnorm(p*n), p, n)
  U = chol(Omega) # By default R's chol fxn returns upper cholesky factor
  X = backsolve(U, Z) # more efficient and stable than acctually inverting
  X = sweep(X, 1, mu, FUN=`+`)
  
  return(t(X))
}

