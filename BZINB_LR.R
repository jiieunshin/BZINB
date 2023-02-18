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

logL_H0 <- function(param){
  
  beta_inter <- c(param[1], param[3])
  beta_coff <- c(param[2], param[4])
  tau <- c(param[5], param[6])
  w <- param[7]
  
  y <- sample_data$y
  x <- sample_data$x
  
  d <- 1 - exp(-1)
  
  v1 <- matrix(0, nrow(y), ncol(y))   # A, B
  
  for(k in 1:2){
    v11 = dnbinom(y[, k], mu = exp(beta_inter[k] + x[, k] * beta_coff[k]), size = 1/tau[k])
    v11 = ifelse(v11>0, v11, 1e-7)
    v1[, k] <- ifelse(y[, k] == 0,
                      log((1 + tau[k] * exp(beta_inter[k] + x[, k] * beta_coff[k]))^(-tau[k]^{-1})),
                      log(v11)
                      )
  }
  
  D1 <- exp(-y[, 1]) - ((1 + d * param[5] * exp(param[1] + param[2] *  x[, 1]))^{-param[5]^{-1}}) 
  D2 <- exp(-y[, 2]) - ((1 + d * param[6] * exp(param[3]  + param[4] * x[, 2]))^{-param[6]^{-1}})
  vv2 = ifelse(1 + param[7] * D1 * D2 < 0, 1e-7, 1 + param[7] * D1 * D2)
  
  v2 = log(vv2)
  
  logbp <- -(sum(v1 + v2))
  
  return(logbp)
}

logL_H1 <- function(param){
  
  beta_inter <- c(param[1], param[3])
  beta_coff <- c(param[2], param[4])
  tau <- c(param[5], param[6])
  w <- param[7]
  phi <- c(param[8], param[9])

  y <- sample_data$y
  x <- sample_data$x
  
  d <- 1 - exp(-1)
  
  v1 <- matrix(0, nrow(y), ncol(y))   # A, B
  
  for(k in 1:2){
    v112 = (1-phi[k]) * dnbinom(y[, k], mu = exp(beta_inter[k] + x[, k] * beta_coff[k]), size = 1/tau[k])
    v112 = ifelse(v112>0, v112, 1e-7)
    v1[, k] <- ifelse(y[, k] == 0,
                      log(phi[k] + (1-phi[k]) * (1 + tau[k] * exp(beta_inter[k] + x[, k] * beta_coff[k]))^(-tau[k]^{-1})),
                      log(v112)
                      )
  }
  
  D1 <- exp(-y[, 1]) - (phi[1] + (1-phi[1]) * (1 + d * param[5] * exp(param[1] + param[2] *  x[, 1]))^{-param[5]^{-1}})
  D2 <- exp(-y[, 2]) - (phi[2] + (1-phi[2]) * (1 + d * param[6] * exp(param[3]  + param[4] * x[, 2]))^{-param[6]^{-1}})
  vv2 = ifelse(1 + param[7] * D1 * D2 < 0, 1e-7, 1 + param[7] * D1 * D2)
  
  v2 = log(vv2)
  
  logbp <- -(sum(v1 + v2))
  
  return(logbp)
}


#######################################################################
# score test for phi1 and phi2
d = 1-exp(-1)
result_score <- lapply(1:25, function(j) 0)
w_seq = c(-2, -1, 0, 1, 2)
tau_seq = c(0.1, 0.5, 1, 1.5, 2)
t = 0
for(ta in tau_seq){
  cat("start tau=", ta, '\n')
  print(Sys.time())
  
  for(w in w_seq){
    cat("start w=", w, '\n')
    print(Sys.time())
    t = t + 1
  for(k in 1:1000) {
    set.seed(k)
    sample_data <- rbzinb(n = 100, phi1 = 0, phi2 = 0, beta10 = .5, beta11 = 1, beta20 = .5, beta21 = 1, tau1 = 1, tau2 = tau, w = w)
    
    if (class(try(optim(par = rep(.5, 9), logL_H1, method = "L-BFGS-B",
                        lower = c(-Inf, -Inf, -Inf, -Inf, 1e-10, 1e-10, -Inf, 1e-10, 1e-10),
                        upper = c(Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf)), silent = TRUE)) == "try-error" |
        class(try(optim(par = rep(.5, 7), logL_H0, method = "L-BFGS-B",
                        lower = c(-Inf, -Inf, -Inf, -Inf, 1e-10, 1e-10, -Inf),
                        upper = c(Inf, Inf, Inf, Inf, Inf, Inf, Inf)), silent = TRUE)) == "try-error" ) {
      
      result_score[[t]][k] = NA
      next
    } else{
      
      logL_H0val = optim(par = rep(.5, 7), logL_H0, method = "L-BFGS-B",
                         lower = c(-Inf, -Inf, -Inf, -Inf, 1e-10, 1e-10, -Inf),
                         upper = c(Inf, Inf, Inf, Inf, Inf, Inf, Inf))$value
      
      logL_H1val = optim(par = rep(.5, 9), logL_H1, method = "L-BFGS-B",
                         lower = c(-Inf, -Inf, -Inf, -Inf, 1e-10, 1e-10, -Inf, 1e-10, 1e-10),
                         upper = c(Inf, Inf, Inf, Inf, Inf, Inf, Inf, 1-(1e-10), 1-(1e-10)))$value
      
      result_score[[t]][k] = 2 * (logL_H0val - logL_H1val)
      
      }
    }
  }
}
qq = c(0.99, 0.95, 0.9, 0.8)
qchi = qchisq(qq, df = 2)
for(i in 1:25){
  cat("tau =", tau_seq[i], " w =", w_seq[i], '\n')
  for(j in 1:4)
  cat("true: prob", 1-qq[j], " 예측:", mean(result_score[[i]] > qchi[j], na.rm = TRUE), "\n")
}

