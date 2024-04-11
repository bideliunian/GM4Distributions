# library(mvtnorm)
# library(fda)
# library(matrixcalc)

############ pathes
function.path = "functions"
working.path = "simulations"
save.path = "simulations/results"


# source all function scipts from the function path
function.sources = list.files(function.path, pattern="\\.R$", full.names=TRUE)
sapply(function.sources, source)

# Global Parameter Settings
seed = 2024
set.seed(seed)
n = 300
p = 10
m = 5 # number of basis used to generate data
tau = 100 # number of observations
model = 'model3'
time.start = proc.time()

####################################
#     PART 1: DATA GENERATION      #
####################################

#    Generate Random Functions and Observation h_ijk       
# h_ijk = Fourier basis func(t_k) %*% delta_ij + iid error

# 0. Generate delta
delta.list = gen_data(n=n, p=p, m=m, model=model, run.ind=seed)
delta = delta.list$X
Omega = delta.list$Omega
Omega0 = delta.list$Omega0

# 1. Observation time
obs.time = seq(1/tau, 1, 1/tau) # vector of observation time points of delta

# 2. Fourier basis function for data generation
b.mat = obs.fourier.bases(obs.time, m)

# 3. Observations h_ijk
h = array(0, c(n, p, tau))
for(i in 1:n){
  for(j in 1:p){
    h[i,j,] = b.mat %*% matrix(delta[i, ((j-1)*m+1) : (j*m)], ncol=1) + rnorm(tau, 0, 0.5)
  }
}
data <- aperm(h, c(1, 3, 2))

####################################
#######     PART 2: FAPO      ######
####################################
ridge = 10
result = round(fapo(data, ridge) / max(fapo(data, ridge)), 2)
