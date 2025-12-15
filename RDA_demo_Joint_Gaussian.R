### install required R packages and data setup
library(spls); library(PAGE); library(MASS); library(lattice)
library(randomForest)
data(yeast)
W = yeast$y; Z = yeast$x
mRNA = colnames(W); gene = colnames(Z)

Z = matrix(as.vector(unlist(Z)),542,106)
W = matrix(as.vector(unlist(W)),542,18)
p = dim(Z)[2]; m = dim(W)[2]; n = dim(W)[1]

Sigma_eta = 0.003 # it can be replaced by 0.005 or 0.007
Sigma_delta = 0.1 # it can be replaced by 0.15 or 0.2

DA1 = Joint_Gaussian(W, Z, sigma_eta= diag(Sigma_eta,p), 
      sigma_delta= diag(Sigma_delta,m), alpha_1 = .1, alpha_2 = .1, 
      label_name = mRNA)


DA3 = Joint_Gaussian(W, Z, sigma_eta= diag(Sigma_eta,p), 
      sigma_delta= diag(Sigma_delta,m), alpha_1 = .1, alpha_2 = .1, 
      alpha_1_list = seq(0.1,0.6,0.1), alpha_2_list = seq(0.1,0.6,0.1), 
      label_name = mRNA)

B1 = DA1$B
B3 = DA3$B
sum((W - Z%*%B1)^2)/ 18
sum((W - Z%*%B3)^2)/ 18

