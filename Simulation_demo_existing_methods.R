library(MASS)
library(glasso)
library(clime)
library(GUEST)
library(huge)
library(mgm)
library(gRim)

set.seed(723)
n = 1000; p = 10; m = 20
mu = rep(0,p); mu_Y = rep(0,m)
sigma_eta = diag(0.1,p); sigma_delta = diag(0.1,m)
Omega = diag(2,m)
Omega0 = matrix(c(2,0,1,0,1,0,2,0,0,0,1,0,2,0,0,0,0,0,
         2,1,1,0,0,1,2),5,5)
Omega[1:5,1:5] = Omega0
S = diag(1,p); X = mvrnorm(n,mu,S)
### Scenario I ###
B = matrix(0,p,m)
B1 = matrix(c(0.2,0,0,0,0.6,0.5,0.1,0,0.9,0.3,0,0.4,0.3,
     0,0,0.7,0,0,0.3,0.8),4,5)
B[1:4,1:5] = B1
S_Y = solve(Omega)
Y = X %*% B + mvrnorm(n,mu_Y,S_Y)

### Generation of error-prone variables ###
W = Y + mvrnorm(n,mu_Y,sigma_delta)
Z = X + mvrnorm(n,mu,sigma_eta)

evaluation = function(U,V) {
d1 = dim(U)[1]; d2 = dim(U)[2]
spe = 0; sen = 0
for(i in 1:d1) {
for(j in 1:d2) {
spe = spe + (U[i,j] ==0 & V[i,j] ==0)*1
sen = sen + (U[i,j] !=0 & V[i,j] !=0)*1
  }
    }
spe = spe / length(which(U==0))
sen = sen / length(which(U!=0))
return(round(c(spe,sen),3))
}

S = cor(W); SW = S - sigma_delta
Method_glasso = glasso(S,rho=0.05)$wi
Method_clime = clime(W)$Omega[[100]]
Method_GUEST = boost.graph(W,thre = 0, ite1 = 30, ite2 = 0, ite3 = 0, sigma_e=0.1, rep = 1)$w
Method_huge = huge.glasso(W)$icov[[10]]
Method_mgm = mgm(data=W, type=rep("g",dim(W)[2]))$pairwise$wadj
Method_gRim = fit_ggm_grips(S,nobs = nrow(W))$K
Method_Wainwright = glasso(SW,rho=0.05)$wi
### Make evaluation

evaluation(Omega, Method_glasso)
evaluation(Omega, (Method_clime>0.1)*1)
evaluation(Omega, Method_huge)
evaluation(Omega, Method_GUEST)
evaluation(Omega, Method_mgm)
evaluation(Omega, Method_gRim)
evaluation(Omega, Method_Wainwright)

