rm(list = ls())

#######################################################################
## functions for mle

rbzinb <- function(n, phi1, phi2, beta10, beta11, beta20, beta21, tau1, tau2, w) {
  
  y <- expand.grid(0:30, 0:30)
  names(y) <- c("y1", "y2")
  y11 <- y$y1
  y22 <- y$y2
  y = cbind(y11, y22)
  
  d <- 1 - exp(-1)
  
  x1 = cbind(1, runif(n))
  x2 = cbind(1, runif(n))
  
  mu1 = exp(x1 %*% c(beta10, beta11))
  mu2 = exp(x2 %*% c(beta20, beta21))
  
  
  c1 <- phi1 + (1 - phi1)*(1 + d * tau1 * mu1)^(-tau1^{-1}); c2 <- phi2 + (1 - phi2)*(1 + d * tau2 * mu2)^(-tau2^{-1})
  
  ry = matrix(0, n, 2)
  
  for(i in 1:n){
    temp <- c()
    
    for(j in 1:nrow(y)){
      y1 <- y11[j]
      y2 <- y22[j]
      if (y1 == 0 & y2 == 0) {
        p <- (1 + w * (exp(-y1) - c1[i]) * (exp(-y2) - c2[i])) *
          (phi1 + (1 - phi1) * dnbinom(0, mu = mu1[i], size = 1/tau1)) *
          (phi2 + (1 - phi2) * dnbinom(0, mu = mu2[i], size = 1/tau2))
      } else if (y1 == 0 & y2 > 0) {
        p <- (1 + w * (exp(-y1) - c1[i]) * (exp(-y2) - c2[i])) *
          (phi1 + (1 - phi1) * dnbinom(0, mu = mu1[i], size = 1/tau1)) *
          (1 - phi2) * dnbinom(y2, mu = mu2[i], size = 1/tau2)
      } else if (y1 > 0 & y2 == 0) {
        p <- (1 + w * (exp(-y1) - c1[i]) * (exp(-y2) - c2[i])) *
          (1 - phi1) * dnbinom(y1, mu = mu1[i], size = 1/tau1) *
          (phi2 + (1 - phi2) * dnbinom(0, mu = mu2[i], size = 1/tau2))
      } else {
        p <- (1 + w * (exp(-y1) - c1[i]) * (exp(-y2) - c2[i])) *
          (1 - phi1) * dnbinom(y1, mu = mu1[i], size = 1/tau1) *
          (1 - phi2) * dnbinom(y2, mu = mu2[i], size = 1/tau2)
      }
      temp[j] <- p
    }
    
    id = min(which(runif(1) <= cumsum(temp)))
    id = ifelse(id == Inf, nrow(y), id)
    ry[i,] = y[id,]
    
  }
  return(list(y = ry, x = cbind(x1[, 2], x2[, 2])))
}

logL <- function(param){
  
  beta_inter <- c(param[1], param[3])
  beta_coff <- c(param[2], param[4])
  tau <- c(param[5], param[6])
  w <- param[7]
  
  y <- sample_data$y
  x <- sample_data$x
  
  d <- 1 - exp(-1)
  
  v1 <- matrix(0, nrow(y), ncol(y))   # A, B
  
  for(k in 1:2){
    v1[, k] <- ifelse(y[, k] == 0,
                      log((1 + tau[k] * exp(beta_inter[k] + x[, k] * beta_coff[k]))^(-tau[k]^{-1})),
                      log(dnbinom(y[, k], mu = exp(beta_inter[k] + x[, k] * beta_coff[k]), size = 1/tau[k])))
  }
  
  D1 <- exp(-y[, 1]) - ((1 + d * param[5] * exp(param[1] + param[2] *  x[, 1]))^{-param[5]^{-1}}) 
  D2 <- exp(-y[, 2]) - ((1 + d * param[6] * exp(param[3]  + param[4] * x[, 2]))^{-param[6]^{-1}})
  vv2 = ifelse(1 + param[7] * D1 * D2 < 0, 1e-10, 1 + param[7] * D1 * D2)
  
  v2 = log(vv2)
  
  logbp <- -(sum(v1 + v2))
  
  return(logbp)
}

#######################################################################
# score test for phi1 and phi2
d <- 1 
k <- 1
result_score <- c()
iter <- 1000


while(1) {
  set.seed(k)
  sample_data <- rbzinb(n = 50, phi1 = 0, phi2 = 0, beta10 = .5, beta11 = 1, beta20 = .5, beta21 = 1, tau1 = .5, tau2 = 1, w = 0.2)
  
  if (class(try(optim(par = rep(.5, 7), logL, method = "L-BFGS-B",
                      lower = c(-Inf, -Inf, -Inf, -Inf, 1e-10, 1e-10, -Inf), upper = c(Inf, Inf, Inf, Inf, Inf, Inf, Inf)), silent = TRUE)) == "try-error") {
    next
    
  } else {
    par <- optim(par = rep(.5, 7), logL, method = "L-BFGS-B",
                 lower = c(-Inf, -Inf, -Inf, -Inf, 1e-10, 1e-10, -Inf), upper = c(Inf, Inf, Inf, Inf, Inf, Inf, Inf))$par
  }
  
  d <- 1 - exp(-1)
  beta_inter <- c(par[1], par[3]); beta_coff <- c(par[2], par[4])
  tau <- c(par[5], par[6]); w <- par[7]
  
  x <- sample_data$x; y <- sample_data$y
  
  mu <- cbind(exp(beta_inter[1]  + beta_coff[1] * x[, 1]), exp(beta_inter[2]  + beta_coff[2] * x[, 2]))
  temp1 <- cbind((1 + d * tau[1] * mu[, 1])^(-tau[1])^{-1}, (1 + d * tau[2] * mu[, 2])^(-tau[2])^{-1})
  temp2 <- cbind((1 + tau[1] * mu[, 1])^(-tau[1])^{-1}, (1 + tau[2] * mu[, 2])^(-tau[2])^{-1})
  D <- cbind(exp(-y[ ,1]) - temp1[, 1], exp(-y[, 2]) - temp1[, 2])
  
  ## garbage
  
  tau1_term <- ifelse(y[, 1] == 0, tau[1]^{-2} * log(1 + tau[1] * mu[, 1]) - ((tau[1]^{-1} * mu[, 1]) / (1 + mu[, 1] * tau[1])),
                      -sapply(1:nrow(y), function(i){
                        if(y[i, 1] == 0){
                          return(0)
                        } else{
                          sv = sapply(0:(y[i, 1] - 1), function(v){ ((1 + tau[1] * v) * tau[1])^{-1} })
                          return(sum(sv))
                        }
                      }) + (y[, 1] / tau[1]) - ((y[, 1] * mu[, 1]) / (1 + tau[1] * mu[, 1])) + tau[1]^{-2} * log(1 + tau[1] * mu[, 1]) - tau[1]^{-1} * mu[, 1] * (1 + tau[1] * mu[, 1])^{-1})
  
  tau2_term <- ifelse(y[, 2] == 0, tau[2]^{-2} * log(1 + tau[2] * mu[, 2]) - ((tau[2]^{-1} * mu[, 2]) / (1 + mu[, 2] * tau[2])),
                      -sapply(1:nrow(y), function(i){
                        if(y[i, 2] == 0){
                          return(0)
                        } else{
                          sv = sapply(0:(y[i, 2] - 1), function(v){ ((1 + tau[2] * v) * tau[2])^{-1} })
                          return(sum(sv))
                        }
                      }) + (y[, 2] / tau[2]) - ((y[, 2] * mu[, 2]) / (1 + tau[2] * mu[, 2])) + tau[2]^{-2} * log(1 + tau[2] * mu[, 2]) - tau[2]^{-1} * mu[, 2] * (1 + tau[2] * mu[, 2])^{-1})
  
  ## partial vector
  
  
  d.phi1 <- - ((w * D[, 2] * (1 - temp1[, 1])) / (1 + w * D[, 1] * D[, 2])) + (((1 - (1 + tau[1] * mu[, 1])^(-tau[1]^{-1})) / (temp2[, 1])) * (y[, 1] == 0) - (y[, 1] != 0))
  d.phi2 <- - ((w * D[, 1] * (1 - temp1[, 2])) / (1 + w * D[, 1] * D[, 2])) + (((1 - (1 + tau[2] * mu[, 2])^(-tau[2]^{-1})) / (temp2[, 2])) * (y[, 2] == 0) - (y[, 2] != 0))
  d.b10 <- (w * D[, 2] * (1 + d * tau[1] * mu[, 1])^(-tau[1]^{-1}-1) * d * mu[, 1]) /(1 + w * D[, 1] * D[, 2]) + ifelse(y[, 1] == 0, -(1 + tau[1] * mu[, 1])^{-1} * mu[, 1], (y[, 1] - mu[, 1]) / (1 + tau[1] * mu[, 1]))
  d.b11 <- (w * D[, 2] * (1 + d * tau[1] * mu[, 1])^(-tau[1]^{-1}-1) * d * mu[, 1] * x[, 1]) /(1 + w * D[, 1] * D[, 2]) + ifelse(y[, 1] == 0, -(1 + tau[1] * mu[, 1])^{-1} * mu[, 1] * x[, 1], (y[, 1] * x[, 1] - mu[, 1] * x[, 1]) / (1 + tau[1] * mu[, 1]))
  d.b20 <- (w * D[, 1] * (1 + d * tau[2] * mu[, 2])^(-tau[2]^{-1}-1) * d * mu[, 2]) /(1 + w * D[, 1] * D[, 2]) + ifelse(y[, 2] == 0, -(1 + tau[2] * mu[, 2])^{-1} * mu[, 2], (y[, 2] - mu[, 2]) / (1 + tau[2] * mu[, 2]))
  d.b21 <- (w * D[, 1] * (1 + d * tau[2] * mu[, 2])^(-tau[2]^{-1}-1) * d * mu[, 2] * x[, 2]) /(1 + w * D[, 1] * D[, 2]) + ifelse(y[, 2] == 0, -(1 + tau[2] * mu[, 2])^{-1} * mu[, 2] * x[, 2], (y[, 2] * x[, 2] - mu[, 2] * x[, 2]) / (1 + tau[2] * mu[, 2]))
  d.tau1 <- - (w * D[, 2] * temp1[, 1] * (tau[1]^{-2} * log(1 + d * tau[1] * mu[, 1]) - ((tau[1]^{-1} * d * mu[, 1]) / (1 + d * mu[, 1] * tau[1])))) / (1 + w * D[, 1] * D[, 2]) + tau1_term
  d.tau2 <- - (w * D[, 1] * temp1[, 2] * (tau[2]^{-2} * log(1 + d * tau[2] * mu[, 2]) - ((tau[2]^{-1} * d * mu[, 2]) / (1 + d * mu[, 2] * tau[2])))) / (1 + w * D[, 1] * D[, 2]) + tau2_term
  d.w <- (D[, 1] * D[, 2]) / (1 + w * D[, 1] * D[, 2])
  
  
  ## Information matrix
  
  gr_df <- data.frame(d.phi1, d.phi2, d.b10, d.b11, d.b20, d.b21, d.tau1, d.tau2, d.w)
  
  I <- matrix(c(colSums(gr_df * d.phi1),
                colSums(gr_df * d.phi2),
                colSums(gr_df * d.b10),
                colSums(gr_df * d.b11),
                colSums(gr_df * d.b20),
                colSums(gr_df * d.b21),
                colSums(gr_df * d.tau1),
                colSums(gr_df * d.tau2),
                colSums(gr_df * d.w)), 9, 9)
  
  score <- matrix(c(sum(d.phi1), sum(d.phi2), rep(0, 7)), 1, 9)
  result_score[k] <- score %*% solve(I) %*% t(score)
  
  k <- k + 1
  
  if (k %% 100 == 0) {
    print(k)
    print(Sys.time())
  }
  
  if (k == iter + 1) {
    break
  }
}

qq = c(0.99, 0.95, 0.9, 0.8)
qchi = qchisq(qq, df = 2)
for(i in 1:4){
  cat("true:", 1-qq[i], " 예측:", mean(result_score > qchi[i], na.rm = TRUE), "\n")
}

