set.seed(2025)

# g(α) function as per theory
g = function(a) {
  return(gamma((1 + a)/2) * gamma((1 - a)/2) / pi)
}

# Our proposed scale estimator using median as pilot location
f1 = function(x, a) {
  y = abs(x - median(x))
  return((mean(y^a) / g(a))^(1 / a))
}
# Geometric mean-based scale estimator (α = 0 limit)
library(psych)
f3 = function(x) {
  y = abs(x - median(x))
  return(geometric.mean(y))
}

# MLE scale estimator
library(cauchypca)
m1 = function(x) {
  mle = cauchy.mle(x)
  return(as.vector(mle$param)[2])
}

# Simulation parameters
n = 14
runs = 1000
alpha = c(seq(-0.99, -0.01, 0.01), seq(0.01, 0.99, 0.01))
q = length(alpha)
mu = 5
sigma = 2

# Storage
x = matrix(0, runs, n)
sig_hat1 = matrix(0, runs, q)
sig_est1 = rep(0, q)
med = rep(0, runs)
qd = rep(0, runs)
mle_sc = rep(0, runs)
est0_1 = rep(0, runs)
bias1 = rep(0, q)
MSE1 = rep(0, q)
v1 = rep(0, q)

# Main simulation loop
for (k in 1:runs) {
  x[k, ] = rcauchy(n, mu, sigma)
  
  med[k] = median(x[k, ])
  qd[k] = (quantile(x[k, ], 0.75) - quantile(x[k, ], 0.25)) / 2
  mle_sc[k] = m1(x[k, ])
  est0_1[k] = f3(x[k, ])
  
  for (i in 1:q) {
    sig_hat1[k, i] = f1(x[k, ], alpha[i])
  }
}

# Bias, variance, MSE calculations
for (i in 1:q) {
  sig_est1[i] = mean(sig_hat1[, i])
  bias1[i] = sig_est1[i] - sigma
  v1[i] = var(sig_hat1[, i])
  MSE1[i] = v1[i] + bias1[i]^2
}

# Special case for α = 0 (geometric mean)
bias0_1 = mean(est0_1 - sigma)
MSE0_1 = mean((est0_1 - sigma)^2)

# Absolute Bias and MSE of sample QD and MLE
abs(mean(qd - sigma))
abs(mean(mle_sc - sigma))
mean((qd - sigma)^2)
mean((mle_sc - sigma)^2)

# Construct α vector including α = 0, and extended bias/MSE
alpha_new = seq(-0.25, 0.25, 0.01)
bias1_new = c(bias1[75:99], bias0_1, bias1[100:124])
MSE1_new = c(MSE1[75:99], MSE0_1, MSE1[100:124])

# Identify α with minimum absolute bias and MSE
alpha_new[which.min(abs(bias1_new))]
abs(bias1_new)[which.min(abs(bias1_new))]
MSE1_new[which.min(abs(bias1_new))]

alpha_new[which.min(MSE1_new)]
abs(bias1_new)[which.min(MSE1_new)]
MSE1_new[which.min(MSE1_new)]

# Compare to QD and MLE scale
alpha_new[which.min(abs(abs(bias1_new) - abs(mean(qd - sigma))))]
abs(bias1_new)[which.min(abs(abs(bias1_new) - abs(mean(qd - sigma))))]
MSE1_new[which.min(abs(abs(bias1_new) - abs(mean(qd - sigma))))]

alpha_new[which.min(abs(MSE1_new - mean((qd - sigma)^2)))]
abs(bias1_new)[which.min(abs(MSE1_new - mean((qd - sigma)^2)))]
MSE1_new[which.min(abs(MSE1_new - mean((qd - sigma)^2)))]

alpha_new[which.min(abs(MSE1_new - mean((mle_sc - sigma)^2)))]
abs(bias1_new)[which.min(abs(MSE1_new - mean((mle_sc - sigma)^2)))]
MSE1_new[which.min(abs(MSE1_new - mean((mle_sc - sigma)^2)))]

# Plot: Absolute Bias vs alpha
alpha_0 = seq(-0.5, 0.5, 0.01)
bias1_0 = c(bias1[50:99], bias0_1, bias1[100:149])
plot(alpha_0, abs(bias1_0), col = 2, type = "l", lwd = 3,
     xlim = c(-0.5, 0.5), ylim = c(0, 0.20),
     xlab = expression(alpha), ylab = "Absolute Bias",
     main = expression(paste("Abs. Bias: ", hat(sigma[alpha]), " vs QD & MLE")))
abline(h = c(0, abs(mean(qd - sigma)), abs(mean(mle_sc - sigma))), col = c(1, 4, 6),
       lwd = 3)
legend(-0.1, 0.08, legend = c(expression(hat(sigma[alpha])), "QD", "MLE"),       
       col = c(2, 4, 6),
       lty = 1, lwd = 3, cex = 0.6)

# Plot: MSE vs alpha
MSE1_0 = c(MSE1[50:99], MSE0_1, MSE1[100:149])
plot(alpha_0, MSE1_0, col = 2, type = "l",  lwd = 3, xlim = c(-0.5, 0.5),
     ylim = c(0, 5),
     xlab = expression(alpha), ylab = "MSE",
     main = expression(paste("MSE: ", hat(sigma[alpha]), " vs QD & MLE")))
abline(h = c(0, mean((qd - sigma)^2), mean((mle_sc - sigma)^2)), col = c(1, 4, 6),
       lwd = 3)
legend(-0.15, 1.25, legend = c(expression(hat(sigma[alpha])), "QD", "MLE"),
       col = c(2, 4, 6), lty = 1, lwd = 3, cex = 0.6)

# Histograms
hist(qd, prob = TRUE, xlab = "QD", main = "Histogram of QD", 
     xlim = c(0, 12), ylim = c(0,0.6))
hist(mle_sc, prob = TRUE, xlab = "MLE", main = "Histogram of MLE",
     xlim = c(0, 12), ylim = c(0,0.6))

hist(sig_hat1[, which.min(abs(bias1_new))], prob = TRUE,
     xlab = expression(hat(sigma[alpha])), main = "Estimator (α min Abs. Bias)",
     xlim = c(0, 12), ylim = c(0,0.6))
hist(sig_hat1[, which.min(MSE1_new)], prob = TRUE,
     xlab = expression(hat(sigma[alpha])), main = "Estimator (α min MSE)",
     xlim = c(0, 12), ylim = c(0,0.6))
