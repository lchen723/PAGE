library(MASS)
library(PAGE)
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
B0 = matrix(c(0.3,0,0,0,0.6,0.5,0.3,0,0.9,0.3,0,0.4,0.3,
     0,0,0.7,0,0,0.3,0.7),4,5)
B[1:4,1:5] = B0
S_Y = solve(Omega)
Y = X %*% B + mvrnorm(n,mu_Y,S_Y)



### Generation of error-prone variables ###
W = Y + mvrnorm(n,mu_Y,sigma_delta)
Z = X + mvrnorm(n,mu,sigma_eta)

#######################################################
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
#######################################################

### Run the conditional likelihood estimation and its 
### model averaging in Sections 3.2 and 3.3
CLE = Cond_Gaussian(W, Z, sigma_eta= diag(0.1,p), sigma_delta=diag(0.1,m), 
       alpha_1 = 0.5, alpha_2 = 0.5, max_iter = 1, tol = 1e-6, 
       label_name = TRUE)

### Run the model averaging of the conditional likelihood estimation
CLE_MA = Cond_Gaussian(W, Z, sigma_eta= diag(0.1,p), sigma_delta=diag(0.1,m), 
          alpha_1 = 0.5, alpha_2 = 0.5, alpha_1_list = seq(0.1,0.6,0.1),
          alpha_2_list = seq(0.1,0.6,0.1),max_iter = 1, tol = 1e-6, 
          label_name = TRUE)



### Make evaluation
Bhat_CLE = (CLE$B>0.1)*1
gammahat_CLE = (abs(CLE$gamma)>0.1)*1
evaluation(B, Bhat_CLE)
evaluation(Omega, gammahat_CLE)

Bhat_CLE_MA = (CLE_MA$B>0.1)*1
gammahat_CLE_MA = (abs(CLE_MA$gamma)>0.1)*1
evaluation(B, Bhat_CLE_MA)
evaluation(Omega, gammahat_CLE_MA)



