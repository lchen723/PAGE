### install required R packages and data setup
library(spls); library(PAGE); library(MASS); library(lattice)
library(MASS)
library(glasso)
library(clime)
library(GUEST)
library(huge)
library(mgm)
library(gRim)
library(randomForest)


data(yeast)
W = yeast$y; Z = yeast$x
mRNA = colnames(W); gene = colnames(Z)

Z = matrix(as.vector(unlist(Z)),542,106)
W = matrix(as.vector(unlist(W)),542,18)
p = dim(Z)[2]; m = dim(W)[2]; n = dim(W)[1]


S = cor(W); SW1 = S - diag(0.1,m); SW15 = S - diag(0.15,m); SW2 = S - diag(0.2,m)
Method_glasso = glasso(S,rho=0.05)$wi
Method_clime = clime(W)$Omega[[100]]
Method_GUEST0 = boost.graph(W,thre = 0, ite1 = 30, ite2 = 0, ite3 = 0, sigma_e=0, rep = 1)$w
Method_GUEST1 = boost.graph(W,thre = 0, ite1 = 30, ite2 = 0, ite3 = 0, sigma_e=0.1, rep = 1)$w
Method_GUEST15 = boost.graph(W,thre = 0, ite1 = 30, ite2 = 0, ite3 = 0, sigma_e=0.15, rep = 1)$w
Method_GUEST2 = boost.graph(W,thre = 0, ite1 = 30, ite2 = 0, ite3 = 0, sigma_e=0.2, rep = 1)$w
Method_huge = huge.glasso(W)$icov[[10]]
Method_mgm = mgm(data=W, type=rep("g",dim(W)[2]))$pairwise$wadj
Method_gRim = fit_ggm_grips(S,nobs = nrow(W))$K
Method_Wainwright1 = glasso(SW1,rho=0.05)$wi
Method_Wainwright15 = glasso(SW15,rho=0.1)$wi
Method_Wainwright2 = glasso(SW2,rho=0.1)$wi

length(which(Method_glasso!=0))/2
length(which(Method_clime!=0))/2
length(which(Method_GUEST0!=0))/2
length(which(Method_GUEST1!=0))/2
length(which(Method_GUEST15!=0))/2
length(which(Method_GUEST2!=0))/2
length(which(Method_huge!=0))/2
length(which(Method_mgm!=0))/2
length(which(Method_gRim!=0))/2
length(which(Method_Wainwright1!=0))/2
length(which(Method_Wainwright15!=0))/2
length(which(Method_Wainwright2!=0))/2


graph_mRNA = function(M) {
net = M
    net = network::network(net, directed = FALSE)
    network::network.vertex.names(net) = paste0("Y", network::network.vertex.names(net))
    graph = GGally::ggnet2(net, size = 10, node.color = "lightgray", 
        label = mRNA, label.size = 3, mode = "circle")

return(graph)

}

graph_mRNA(Method_glasso)
graph_mRNA(Method_clime)
graph_mRNA(Method_GUEST0)
graph_mRNA(Method_GUEST1)
graph_mRNA(Method_GUEST15)
graph_mRNA(Method_GUEST2)
graph_mRNA(Method_huge)
graph_mRNA(Method_mgm)
graph_mRNA(Method_gRim)
graph_mRNA(Method_Wainwright1)
graph_mRNA(Method_Wainwright15)
graph_mRNA(Method_Wainwright2)


