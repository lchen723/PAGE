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

### analysis of the nonparametric method:
### run the function "NP_Graph"
DA0 = NP_Graph(W, Z, sigma_eta = diag(Sigma_eta,p), rho=0.2, 
              sigma_delta = Sigma_delta, label_name = mRNA)

### summarize important variables
est0 = DA0$importance_matrix
id0 = NULL; B0 = NULL
for(i in 1:18){
id0 = c(id0,which(est0[,i] >=10))
B0 = c(B0,gene[which(est0[,i] >=10)])
}
barchart(as.vector(B0))

### print PSE
DA0$PSE

### draw estimated curves based on variables corrected 
### by regression calibration
id = which(gene == "SWI6_YPD")  ## can be replaced by STE12_YPD, NDD1_YPD, 
                                ## SWI6_YPD or MBP1_YPD
W_Hat = DA0$W_hat[,1]  #1 or 9
Z_Hat = DA0$Z_hat[,id]

DATA = data.frame(W_Hat, Z_Hat)
model = randomForest(W_Hat ~ Z_Hat, data = DATA)
pred = predict(model)

ord_1 = order(Z_Hat)
smooth_data_1 = lowess(Z_Hat[ord_1], pred[ord_1], f = 0.5) 
plot(smooth_data_1$x, smooth_data_1$y,
     type = "l",     
     xlab = "SWI6_YPD",
     ylab = "alpha0",  #0 or 56
     lwd = 2)
