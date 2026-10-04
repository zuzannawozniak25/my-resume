#Project - Principal Component Analysis (PCA)
#part 1 - Manual Implementation and Diagnostics

#load data and extract feature matrix

danePCA <- read.csv("danePCA.csv", sep = ";")
head(danePCA)
dim(danePCA)
str(danePCA)
X <- danePCA[ ,-c(1, 2)]

#compute mean vector, covariance matrix, and standardize data (Y)
n <- nrow(X)

#mean vector
m1 <- rep(1,n)
wektor_srednich <- (1/n)*t(X)%*%(m1)
twa <- t(wektor_srednich)
XX <- as.matrix(X)           

#covariance matrix
macierz_kowariancji <- (1/(n-1))*t(XX)%*%(diag(n)-(1/n)*m1%*%t(m1))%*%XX
macierz_kowariancji

#data standarization

library(matlib)

d <- sqrt(diag(macierz_kowariancji))
m <- diag(1/d)

Y <- (XX-m1%*%t(wektor_srednich))%*%m

#standardized data mean vector
wektor_srednich_std <- (1/n)*t(Y)%*%(m1)
tstd <- t(wektor_srednich_std)
#standardized data covariance matrix
Sy <- (1/(n-1))*t(Y)%*%(diag(n)-(1/n)*m1%*%t(m1))%*%Y

#compute eigenvalues and eigenvectors
eigen_decomp <- eigen(Sy)

#eigenvalues
eval <- eigen_decomp$values

#eigenvectors
evec <- eigen_decomp$vectors

#variance explained by principal components

wariancja <- eval
proporcja_wariancji <- wariancja / sum(wariancja)
skumulowana_proporcja <- cumsum(proporcja_wariancji)
PCA <- data.frame(
  PC = paste0("PC", 1:length(eval)),
  wariancja = round(wariancja, 4),
  proporcja_wariancji = round(proporcja_wariancji, 4),
  skumulowana_proporcja = round(skumulowana_proporcja, 4)
)
tPCA <- t(PCA[, -1])  
colnames(tPCA) <- PCA$PC

#linear combinations for PC1 and PC2 & scree Plots

z1 <- paste0("z1 = ", paste0(round(evec[,1], 4), " * y", 1:30, collapse = " + "))
z2 <- paste0("z2 = ", paste0(round(evec[,2], 4), " * y", 1:30, collapse = " + "))

nazwy_zmiennych <- paste0("y", 1:30)

z_1 <- data.frame(zmienna = nazwy_zmiennych, wspolczynnik = z1)
z_2 <- data.frame(zmienna = nazwy_zmiennych, wspolczynnik = z2)

#scree plot (absolute variance)
plot(eval, type = "b", main = "Wariancja (bezwzględna)", xlab = "skladowe glówne", ylab = "wariancja", 
     pch = 19, col = "black")
abline(h = 1, col = "red", lty = 2)

#scree plot (relative variance)
plot(proporcja_wariancji, type = "b", main = "Wyjaśniana wariancja (względna)", 
     xlab = "skladowe glówne", ylab = "Wyjasniana wariancja", pch = 19, col = "black")


#correlation between original variables and principal components

Z <- Y%*%evec
korelacje <- matrix(0, nrow = ncol(Y), ncol = ncol(Z))
for (j in 1:ncol(Z)) {
  korelacje[, j] <- evec[, j] * sqrt(eval[j])
}

r2 <- korelacje^2
suma_r2 <- rowSums(r2)

x <- korelacje[,1]  #PC1
y <- korelacje[,2]  #PC2

nazwy_zmiennych <- as.character(1:30)
plot(x, y, xlim = c(-2, 2), ylim = c(-1, 1), asp = 1,
     xlab = "PC1", ylab = "PC2",
     main = "Korelacje zmiennych PC1 i PC2",
     pch = 19, cex = 0.4, col = "black")
text(x, y, labels = nazwy_zmiennych, pos = 3, cex = 0.5, col = "black")
theta <- seq(0, 2 * pi, length.out = 100)
lines(cos(theta), sin(theta), col = "red")

#component selection (Kaiser criterion: variance > 1) and cumulative variance

skladowe_1 <- which(eval > 1)
n_skladowych <- length(skladowe_1)

plot(skumulowana_proporcja, 
     type = "b", 
     main = "Skumulowana proporcja wariancji",
     xlab = "składowe główne", 
     ylab = "wyjaśniana wariancja",
     pch = 19, cex = 0.5)
abline(v = n_skladowych, col = "red", lty = 2)
abline(h = skumulowana_proporcja[n_skladowych], col = "red", lty = 2)

##
wariancja_wyjasniona <- skumulowana_proporcja[skladowe_1[n_skladowych]]
cat("Wartość wyjaśnianej wariancji przez wybrane składowe:", round(wariancja_wyjasniona, 4), "\n")

#data projection onto PC1, PC2, and PC3
PC1 <- Z[,1]
PC2 <- Z[,2]
PC3 <- Z[,3]

stripchart(PC1,
           type = "p", 
           main = "Projekcja na PC1",
           xlab = "PC1",
           pch = 19, col = "blue", cex = 0.5)

plot(PC1, PC2, 
     main = "Projekcja na PC1 i PC2",
     xlab = "PC1", 
     ylab = "PC2",
     pch = 19, col = "black", cex = 0.5)

#install.packages("scatterplot3d")
library(scatterplot3d)
scatterplot3d(PC1, PC2, PC3,
              main = "Projekcja na PC1, PC2, PC3",
              xlab = "PC1", ylab = "PC2", zlab = "PC3",
              pch = 19, color = "black", cex.symbols = 0.5,
              xlim = c(-15,5), ylim = c(-10,5))

#2D projection with target variable (Malignant vs Benign) and 95% confidence ellipses
diagnosis <- danePCA[, 2]
kolory <- c("M" = "red", "B" = "blue")  #kolory dla klas

#scatter plot colored by class
plot(PC1, PC2,
     main = "Projekcja na PC1 i PC2 (złośliwe vs łagodne)",
     xlab = "PC1", ylab = "PC2",
     col = kolory[diagnosis], pch = 19, cex = 0.6)

legend("bottomleft", legend = c("Złośliwe", "Łagodne"),
       col = c("red", "blue"), pch = 19)

#function to draw 95% confidence ellipse for a given group
rysuj_elipse <- function(grupa, kolor) {
  idx <- diagnosis == grupa
  srodek <- colMeans(cbind(PC1[idx], PC2[idx]))
  cov_mat <- cov(cbind(PC1[idx], PC2[idx]))
  eig <- eigen(cov_mat)
  skale <- sqrt(qchisq(0.95, df = 2) * eig$values)
  wektory <- eig$vectors
  
  theta <- seq(0, 2 * pi, length.out = 100)
  elipsa_punkty <- t(srodek + wektory %*% diag(skale) %*% rbind(cos(theta), sin(theta)))
  
  lines(elipsa_punkty, col = kolor, lwd = 2)
  points(srodek[1], srodek[2], pch = 4, col = kolor, cex = 2, lwd = 2)  # Krzyżyk w środku
}

#draw ellipses for both groups
for (g in c("M", "B")) {
  rysuj_elipse(g, kolory[g])
}


par(mar = c(5, 5, 4, 2) + 0.1)
