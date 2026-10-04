# Project - Principal Component Analysis (PCA)
# PART 2

library(stats)     
library(factoextra)
library(ggplot2)    
library(ggpubr)   

#data loading and preparation

danePCA <- read.csv("danePCA.csv", sep = ";")
head(danePCA)
dim(danePCA)
str(danePCA)
X <- danePCA[ ,-c(1, 2)]
X <- as.matrix(X)     

#install.packages("FactoMineR")
library("FactoMineR")
dane.pca <- PCA(X, scale.unit = TRUE, ncp = 5, graph = TRUE)

#descriptive statistics and PCA computation
n <- nrow(X)
wektor_srednich <- colMeans(X)
macierz_kowariancji <- cov(X)

pca_result <- prcomp(X, center = TRUE, scale. = TRUE) # PCA z centrowaniem i skalowaniem 
Y <- pca_result$x
Y <- prcomp(X, center = TRUE, scale. = TRUE)
#wektor_srednich_std <- colMeans(Y)
Sy <- cov(Y)

#eigenvalues and eigenvectors extraction
#eigenvalues and eigenvectors decomposition
eigen_decomp <- eigen(Sy)

#eigenvalues
eval <- eigen_decomp$values

#eigenvectors
evec <- eigen_decomp$vectors


#extracting eigenvalues
#install.packages("factoextra")
library("factoextra")
pca.eigen <- get_eigenvalue(dane.pca)

#scree Plot visualization
fviz_eig(pca_result, choice = "variance", geom = "bar", ncp = 10, main = "Scree Plot")



#variables correlation plot
viz_pca_var(pca_result, geom = c("arrow", "text"), repel = TRUE, col.var = "black", arrowsize = 0.5, labelsize = 4, main = "Korelacja zmiennych z PC1 i PC2") 
pca_var <- get_pca_var(pca_result) 


#scree Plot with Kaiser criterion line (variance = 1)
fviz_eig(pca_result, choice = "variance", geom = "bar", ncp = 10, main = "Scree Plot z linią wariancji = 1", addlabels = TRUE) + geom_hline(yintercept = 1, linetype = "dashed", color = "red") 
n_components <- sum(pca_eigen$variance > 1) 
cat("Liczba składowych z wariancją > 1:", n_components, "\n") 

















