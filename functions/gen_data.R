#####################
# Generates multidimensional data for a graphical model.
#
# @param n sample size.
# @param p number of nodes.
# @param m dimension of vectors on each nodes.
# @param model.
#   - "AR2": AR2 model with Theta_{j, j-1} = 0.4*I_m and Theta_{j, j-2} = 0.2*I_m.
#   - "Random": Random block-sparse model.
#   - "TwoHubs": Two hubs connected model.
#   - "Theta1": Model with Theta_{1, j} = 0.2*I_m.
# @param run.ind Random seed.
# @return An n by p*d matrix.
#############################
require(MASS)

gen_data = function(n, p, m, model, run.ind){
  
  if (m < 2) {stop('m must be greater or equal to 2')}
  if (p < 2) {stop('p must be greater or equal to 2')}
  
  Omega = matrix(0, nrow = m*p, ncol = m*p)
  Omega0 = matrix(0, nrow = p, ncol = p)
  
  if(model == 'model1'){
    A = diag(m)
    for (j in 1:(p-1)) {
      index.j = (((j-1)*m+1):(j*m))
      Omega[index.j,(index.j+m)] = 0.4*A 
      Omega0[j, j+1] = 1
      if(j != p-1) {
        Omega[index.j,index.j+2*m] = 0.2*A
        Omega0[j, j+2] = 1
      } 
    }
    Omega = Omega + t(Omega) + diag(m*p)
    Omega0 = Omega0 + t(Omega0)
  }
  
  if(model == 'model2'){
    delta = 0.3
    set.seed(run.ind)
    for(i in 1:(p/2-1)){
      index.i = (((i-1)*m+1):(i*m))
      for(j in (i+1):(p/2)){
        index.j = (((j-1)*m+1):(j*m))
        c = rbinom(n=1, size = 1, prob = 0.3)
        if(c){
          Omega[index.i, index.j] = 0.3*diag(m)
          Omega0[i, j] = 1
        }
      }
    }
    for(i in (p/2+1):(p-1)){
      index.i = (((i-1)*m+1):(i*m))
      for(j in (i+1):p){
        index.j = (((j-1)*m+1):(j*m))
        c = rbinom(n=1, size = 1, prob = 0.1)
        if(c){
          Omega[index.i, index.j] = 0.3*diag(m)
          Omega0[i, j] = 1
        }
      }
    }
    Omega = Omega + t(Omega)
    diag(Omega) = 0
    ee = min(eigen(Omega, only.values = T)$values)
    diag(Omega) = ifelse(ee < 0, -ee + 1.0, 1.0)
    Omega0 = Omega0 + t(Omega0) + diag(p)
  }
  
  if(model == 'model3'){
    set.seed(run.ind)
    seed = run.ind
    G.true = edge.gen(p=p, sparsity=0.99, hubnumber=2, hubsparsity=0.5)$G
    Omega = Omega.gen(G=G.true, M=m)
    Omega = hub.trans(Omega, p)
    diag(Omega) = 0
    ee = min(eigen(Omega, only.values = T)$values)
    diag(Omega) = ifelse(ee < 0, -ee + 1.0, 1.0)
    Omega0 = G.true
  }
  
  if(model == 'model4'){
    for (j in 1:(p-1)) {
      index.j = (j*m+1):((j+1)*m)
      if(j <= p-1){
        Omega[(1:m),(index.j)] = 0.2*diag(m)
      }
      Omega0[1,(2:p)] = 1
    }
    Omega = Omega + t(Omega)
    diag(Omega) = 0
    ee = min(eigen(Omega, only.values = T)$values)
    diag(Omega) = ifelse(ee < 0, -ee + 0.1, 0.1)
    Omega0 = Omega0 + t(Omega0) + diag(p)
  }
  
  if(model == 'model5'){
    set.seed(run.ind)
    seed = run.ind
    
    G1 = edge.gen(p=p/2, sparsity=0.99, hubnumber=1, hubsparsity=0.5)$G
    G2 = edge.gen(p=p/2, sparsity=0.99, hubnumber=1, hubsparsity=0.5)$G
    G = matrix(0, p, p)
    G[1:(p/2), 1:(p/2)] = G1
    G[(p/2+1):p, (p/2+1):p] = G2
    Omega = Omega.gen.new(G=G, m=m)
    Omega0 = G
  }
  
  Sigma = solve(Omega)
  diag.sigma = matrix(0, m*p, m*p)
  for (j in 1:p){
    index.j = (((j-1)*m+1):(j*m))
    diag.sigma[index.j, index.j] = matpower(Sigma[index.j, index.j], -1/2)
  }
  
  cov.X = diag.sigma %*% Sigma %*% diag.sigma
  X.pre = rmvnorm(n, sigma = cov.X)
  X = X.pre
  
  return(list(X = X, Omega = Omega, Omega0 = Omega0, X.pre = X.pre))
}