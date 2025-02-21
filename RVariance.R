# Provides the expectation and variance-covariance of LRS
# Can be used for any strutured population projeciton model
# We use the model of the illustration.

rm(list = ls())


## Inputs

s <- 2  # Number of states

# Survival probabilities (expectation of survival)
S <- matrix(c(0.3, 0, 0.2, 0.7), nrow = 2, byrow = TRUE)

# Reproduction probabilities (expectation of reproduction)
F <- matrix(c(0.5, 0.5, 0, 0.5), nrow = 2, byrow = TRUE)

# Variance of reproduction
VF <- matrix(c(0.25, 0.25, 0, -0.25, 0, -0.25,0 , 0.25), nrow = 4, byrow = TRUE)



# Intermediary outputs
# Vector of expectation of reproduction
fT <- rep(1, s) %*% F

# Variance survival
J <- matrix(0, s^2, s)
for (i in 1:s) {
  ei <- rep(0, s)
  ei[i] <- 1
  J <- J + kronecker(ei,ei)%*%t(ei) }
VS <- J %*% S - kronecker(S, S) %*% J  # (eq.4)

# Projection of abundance vector over time
n0 <- c(1,1)
Vn1 <- (VS + VF) %*% n0
Vn1
matrix(Vn1, s, s)

# Expectation of LRS
R <- F %*% solve(diag(s) - S)  # (eq 8, next generation matrix)
R

rT <- rep(1, s) %*% F %*% solve(diag(s) - S)  # Expectation of Total LRS

# Variance-Covariance of LRS
N <- solve(diag(s) - S)  # Fundamental matrix: expected time spent in stages

VR <- (VF + kronecker(R, R) %*% VS) %*% N  # (eq.9)
VR

matrix(VR[, 1], s, s)  # Variance-covariance of LRS for a type 2 (born medium) individual

# Variance of total LRS
VrT <- (kronecker(rT, rT) %*% VS + matrix(1, 1, s^2) %*% VF) %*% solve(diag(s) - S)  # (eq.11)
