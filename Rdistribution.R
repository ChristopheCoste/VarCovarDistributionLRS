#install.packages("expm")   
library(expm)

# Clear workspace
rm(list = ls())

# Inputs
epsilon1 <- 0.00001  # Tolerance for approximation of maximum age
epsilon1 <-0.000000000001
s <- 2  # Number of classes
S <- matrix(c(0.3, 0, 0.2, 0.7), nrow = s, byrow = TRUE)  # Expectation of survival

# Probability distribution of reproduction
# We store the probability to have k newborns at a time step into a vector
Fcal <- matrix(c(0.5, 0, 0.5, 1), nrow = 2, byrow = TRUE)  
alpha <- nrow(Fcal) - 1  # Maximum number of offspring that can be produced by an individual in a time step


# Computing the maximum age
omega <- 0
smax <- 1
while (smax > epsilon1) {
  omega <- omega + 1
  smax <- max(rowSums(S %^% omega))
}
omega  # Age at which less than a proportion epsilon of individuals are expected to still be alive

# Probability distribution of the total LRS
maxr <- alpha * omega  # Maximum coefficient of polynomial of r = maximum LRS allowed by the approximation + 1

Rcal <- matrix(0, maxr + 1, s)  
Rcal[1:(alpha + 1), ] <- Fcal  # Initialization (first line of equation 18)

for (age in (omega - 1):1) {  # Iteration of equation 18
  Qcal <- Rcal %*% S  
  Qcal[1, ] <- Qcal[1, ] + 1 - colSums(S)
  
  for (j in 1:s) {  # Convolution of equation 18
    vect <- rev(convolve(rev(Qcal[, j]), Fcal[, j], type = "open"))  
    Rcal[, j] <- vect[1:(maxr + 1)]
  }
}
Rcal


# Plot results
kmax <- 15  # Max number for figure
plot(0:kmax, Rcal[1:(kmax + 1), 1], type = "l", main = "Polyrt")
plot(0:kmax, Rcal[1:(kmax + 1), 2], type = "l", main = "Polyrt")

#install.packages("ggplot2")  
#install.packages("gridExtra")  
library(ggplot2)
library(gridExtra)  # For arranging plots

kmax <- 15  # Max number for figure

# Create data frames for plotting
df1 <- data.frame(LRS = 0:kmax, Probability = Rcal[1:(kmax + 1), 1])
df2 <- data.frame(LRS = 0:kmax, Probability = Rcal[1:(kmax + 1), 2])

# First plot
p1 <- ggplot(df1, aes(x = LRS, y = Probability)) +
  geom_bar(stat = "identity", fill = "blue") +
  labs(x = "Lifetime reproductive success", y = "Probability", title = "Total LRS for type 'small'") +
  theme_minimal()

# Second plot
p2 <- ggplot(df2, aes(x = LRS, y = Probability)) +
  geom_bar(stat = "identity", fill = "red") +
  labs(x = "Lifetime reproductive success", y = "Probability", title = "Total LRS for type 'large'") +
  theme_minimal()

# Arrange the plots side by side
grid.arrange(p1, p2, ncol = 2)

