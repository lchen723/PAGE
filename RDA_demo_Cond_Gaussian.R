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

DA2 = Cond_Gaussian(W, Z, sigma_eta= diag(Sigma_eta,p), 
      sigma_delta= diag(Sigma_delta,m), alpha_1 = 0.5, alpha_2 = 0.5,
      max_iter = 1, tol = 1e-6, label_name = mRNA)

DA4 = Cond_Gaussian(W, Z, sigma_eta= diag(Sigma_eta,p), 
      sigma_delta= diag(Sigma_delta,m), alpha_1 = 0.5, alpha_2 = 0.5, 
      alpha_1_list = seq(0.1,0.6,0.1), alpha_2_list = seq(0.1,0.6,0.1), 
      max_iter = 1, tol = 1e-6, label_name = mRNA) 


B2 = DA2$B
B4 = DA4$B

sum((W - Z%*%B2)^2)/ 18
sum((W - Z%*%B4)^2)/ 18

