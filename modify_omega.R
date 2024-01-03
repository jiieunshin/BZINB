rm(list = ls())

library(parallel)
library(dplyr)
library(MASS)

# numCores <- detectCores() - 1
cl <- makeCluster(45)
result_list = list()
boot_list = list()

please <- function(j) {
  library(dplyr)
  library(MASS)
  w = 0
  RBNBZI <- function(n, b10, b11, b20, b21, gam10, gam11, w, tau1, tau2, grid = 20) {
    
    find_idx <- function(mat) {
      if(sum(mat) == 0){
        return(c(grid-1, grid-1))
      } else{
        for (row in 1:grid) {
          for (col in 1:grid) {
            if (mat[row, col] == TRUE) {
              return(c(row - 1, col - 1))
            }
          }
        }
      }
    }
    
    d <- 1 - exp(-1)
    
    x1 <- cbind(1, runif(n))
    x2 <- cbind(1, runif(n))
    z <- cbind(1, runif(n))
    
    mu1 <- exp(x1 %*% c(b10, b11));
    mu2 <- exp(x2 %*% c(b20, b21));
    phi <- exp(z %*% c(gam10, gam11)) / (1 + exp(z %*% c(gam10, gam11)))
    
    c1 <- (1 + d * mu1 * tau1)^{-1/tau1}
    c2 <- (1 + d * mu2 * tau2)^{-1/tau2}
    
    y <- matrix(0, nrow = n, ncol = 2)
    
    for (i in 1:n) {
      PX <- cbind(dnbinom(0:(grid - 1), mu = mu1[i], size = 1/tau1), dnbinom(0:(grid - 1), mu = mu2[i], size = 1/tau2))
      PTMPY <- cbind(exp(-(0:(grid - 1))) - (1 + d * mu1[i] * tau1)^{-1/tau1}, exp(-(0:(grid - 1))) - (1 + d * mu2[i] * tau2)^{-1/tau2})
      
      TXY = matrix(PX[, 1], nrow = grid, ncol = 1) %*% matrix(PX[, 2], nrow = 1, ncol = grid)  * 
        (1 + w * matrix(PTMPY[, 1], nrow = grid, ncol = 1) %*% matrix(PTMPY[, 2], nrow = 1, ncol = grid))
      
      SXY = matrix(cumsum(TXY %>% t() %>% as.vector()), nrow = 20, ncol = 20, byrow = TRUE)
      
      RAN_N <- runif(1)
      
      if (RAN_N < phi[i]) {
        y[i, ] <- c(0, 0)
      } else {
        RAN_A <- runif(1)
        mat = RAN_A < SXY
        find_idx(mat)
        y[i, ] <- find_idx(mat)
      }
    }
    return(list(y = y, x1 = x1, x2 = x2, z = z, c1 = c1, c2 = c2))
  }
  
  logH0 <- function(params) {
    beta1 <- params[1:2]
    beta2 <- params[3:4]
    gamma <- params[5:6]
    tau <- params[7:8]
    
    x1 <- sample_data$x1
    x2 <- sample_data$x2
    z <- sample_data$z
    
    y <- sample_data$y
    
    mu1 <- exp(x1 %*% beta1)
    mu2 <- exp(x2 %*% beta2)
    phi <- exp(z %*% gamma) / (1 + exp(z %*% gamma))
    
    ind <- rowSums(y) == 0
    
    
    v1 <- log(phi + (1 - phi) * (1 + tau[1] * mu1)^{-1/tau[1]} * (1 + tau[2] * mu2)^{-1/tau[2]}) * ind
    v2 <- log((1 - phi) * dnbinom(x = y[, 1], mu = mu1, size = 1/tau[1]) * dnbinom(x = y[, 2], mu = mu2, size = 1/tau[2])) * !ind
    
    logbp <- sum(v1) + sum(v2)
    
    return(-logbp)
  }
  
  
  n_grid <- c(100, 200, 500)
  tau2_grid <- c(0.4, 0.6, 0.8, 1.0, 1.2)
  g2_grid <- c(-0.997, -0.186, 0.794)
  
  param_grid <- t(expand.grid(g2_grid, tau2_grid, n_grid))
  
  # for(j in 6:30){
  # for(j in 1:45){
  iteration = 1000
  result = result2 =  c()
  param_result <- matrix(0, iteration, 8)
  gr_result <- matrix(0, iteration, 9)
  for(i in 1:iteration) {
    set.seed(i)
    
    sample_data <- RBNBZI(n = param_grid[3, j], b10 = .2, b11 = .4, b20 = .4, b21 = .8, gam10 = -1.2, gam11 = param_grid[1, j],
                          w = 0, tau1 = .8, tau2 = param_grid[2, j])
    one <- try(optim(par = rep(.5, 8), fn = logH0, lower = c(rep(-Inf, 6), 1e-10, 1e-10), method = "L")$par, silent = TRUE)
    
    if (class(one) == "try-error") {
      next
    } else {
      param <- one
      param_result[i,] <- param
    }
    
    x1 <- sample_data$x1
    x2 <- sample_data$x2
    z <- sample_data$z
    y <- sample_data$y
    d <- 1 - exp(-1)
    ind <- (y[, 1] == 0 & y[, 2] == 0)
    
    beta1 <- param[1:2]
    beta2 <- param[3:4]
    gamma <- param[5:6]
    tau <- param[7:8]
    
    mu1 <- exp(x1 %*% beta1)
    mu2 <- exp(x2 %*% beta2)
    phi <- exp(z %*% gamma) / (1 + exp(z %*% gamma)) + 1e-10
    # psi <- phi/(1-phi)
    
    c <- cbind((1 + d * mu1 * tau[1])^{-1/tau[1]}, (1 + d * mu2 * tau[2])^{-1/tau[2]})
    D <- (1 + w * (exp(-y[, 1]) - c[, 1]) * (exp(-y[, 2]) - c[, 2]))
    
    dbeta10 <- (- (1 - phi) * (1 + tau[1] * mu1)^{-1/tau[1]-1} * (1 + tau[2] * mu2)^{-1/tau[2]} * mu1 * x1[, 1] / (phi + (1 - phi) * (1 + tau[1] * mu1)^{-1/tau[1]} * (1 + tau[2] * mu2)^{-1/tau[2]})) * ind +
      ((y[, 1] - mu1) * x1[, 1] / (1 + tau[1] * mu1)) * !ind
    dbeta11 <- (- (1 - phi) * (1 + tau[1] * mu1)^{-1/tau[1]-1} * (1 + tau[2] * mu2)^{-1/tau[2]} * mu1 * x1[, 2] / (phi + (1 - phi) * (1 + tau[1] * mu1)^{-1/tau[1]} * (1 + tau[2] * mu2)^{-1/tau[2]})) * ind +
      ((y[, 1] - mu1) * x1[, 2] / (1 + tau[1] * mu1)) * !ind
    dbeta20 <- (- (1 - phi) * (1 + tau[2] * mu2)^{-1/tau[2]-1} * (1 + tau[1] * mu1)^{-1/tau[1]} * mu2 * x2[, 1] / (phi + (1 - phi) * (1 + tau[1] * mu1)^{-1/tau[1]} * (1 + tau[2] * mu2)^{-1/tau[2]})) * ind +
      ((y[, 2] - mu2) * x2[, 1] / (1 + tau[2] * mu2)) * !ind
    dbeta21 <- (- (1 - phi) * (1 + tau[2] * mu2)^{-1/tau[2]-1} * (1 + tau[1] * mu1)^{-1/tau[1]} * mu2 * x2[, 2] / (phi + (1 - phi) * (1 + tau[1] * mu1)^{-1/tau[1]} * (1 + tau[2] * mu2)^{-1/tau[2]})) * ind +
      ((y[, 2] - mu2) * x2[, 2] / (1 + tau[2] * mu2)) * !ind
    dgam10 <- (z[, 1] * exp(z %*% gamma) / (1 + exp(z %*% gamma))^2 * (1 - (1 + tau[1] * mu1)^{-1/tau[1]} * (1 + tau[2] * mu2)^{-1/tau[2]}) / (phi + (1 - phi) * (1 + tau[1] * mu1)^{-1/tau[1]} * (1 + tau[2] * mu2)^{-1/tau[2]})) * ind -
      (z[, 1] * exp(z %*% gamma) / (1 + exp(z %*% gamma))^2 / (1 - phi)) * !ind
    dgam11 <- (z[, 2] * exp(z %*% gamma) / (1 + exp(z %*% gamma))^2 * (1 - (1 + tau[1] * mu1)^{-1/tau[1]} * (1 + tau[2] * mu2)^{-1/tau[2]}) / (phi + (1 - phi) * (1 + tau[1] * mu1)^{-1/tau[1]} * (1 + tau[2] * mu2)^{-1/tau[2]})) * ind -
      (z[, 2] * exp(z %*% gamma) / (1 + exp(z %*% gamma))^2 / (1 - phi)) * !ind
    
    
    
    aa1 <- (tau[1]^{-2}*log(1+tau[1]*mu1) - tau[1]^{-1}*mu1 / (1+tau[1]*mu1))
    
    aa2 <- ((1-phi) * (1 + tau[2] * mu2)^{-1/tau[2]} * (1 + tau[1] * mu1)^{-1/tau[1]} * aa1 ) / (phi + (1-phi) * (1 + tau[2] * mu2)^{-1/tau[2]} * (1 + tau[1] * mu1)^{-1/tau[1]}) * ind
    
    aa3 = sapply(1:nrow(y), function(m){
      if(y[m, 1] == 0){
        val = 0
      } else{
        sv = sapply(1:y[m, 1], function(v){ (y[m, 1] - v)/(1 - v * tau[1] + y[m, 1] * tau[1]) })
        val = sum(sv)
      }
      return(val) 
    })
    
    aa4 <- ( aa3 - y[, 1] * mu1 / (1 + tau[1] * mu1) - mu1 * tau[1]^{-1} / (1 + tau[1] * mu1) + tau[1]^{-2} * log(1 + tau[1] * mu1) ) * !ind
    
    # 
    ab1 <- (tau[2]^{-2}*log(1+tau[2]*mu2) - tau[2]^{-1}*mu2 / (1+tau[2]*mu2))
    
    ab2 <- ((1-phi) * (1 + tau[2] * mu2)^{-1/tau[2]} * (1 + tau[1] * mu1)^{-1/tau[1]} * ab1 ) / (phi + (1-phi) * (1 + tau[2] * mu2)^{-1/tau[2]} * (1 + tau[1] * mu1)^{-1/tau[1]}) * ind
    
    ab3 = sapply(1:nrow(y), function(m){
      if(y[m, 2] == 0){
        val = 0
      } else{
        sv = sapply(1:y[m, 2], function(v){ (y[m, 2] - v)/(1 - v * tau[2] + y[m, 2] * tau[2]) })
        val = sum(sv)
      }
      return(val) 
    })
    
    ab4 <- ( ab3 - y[, 2] * mu2 / (1 + tau[2] * mu2) - mu2 * tau[2]^{-1} / (1 + tau[2] * mu2) + tau[2]^{-2} * log(1 + tau[2] * mu2) ) * !ind
    
    g_tau1 <- (aa2 + aa4)
    g_tau2 <- (ab2 + ab4)
    
    g_w <- ( (1-phi) * (1 + tau[2] * mu2)^{-1/tau[2]} * (1 + tau[1] * mu1)^{-1/tau[1]} * (1 - c[, 1]) *  (1 - c[, 2]) 
             / (phi + (1-phi) * (1 + tau[2] * mu2)^{-1/tau[2]} * (1 + tau[1] * mu1)^{-1/tau[1]} * (1 + w * (1 - c[, 1]) *  (1 - c[, 2]))) ) * ind + 
      ((exp(-y[, 1]) - c[, 1]) *  (exp(-y[, 2]) - c[, 2]) / (1 + w * (exp(-y[, 1]) - c[, 1]) *  (exp(-y[, 2]) - c[, 2]))) * !ind
    
    
    gr_df <- data.frame(g_w, dbeta10, dbeta11, dbeta20, dbeta21, dgam10, dgam11, g_tau1, g_tau2)
    gr_result[i,] <- colSums(gr_df)
    
    I <- matrix(c(colSums(gr_df * g_w),
                  colSums(gr_df * dbeta10),
                  colSums(gr_df * dbeta11),
                  colSums(gr_df * dbeta20),
                  colSums(gr_df * dbeta21),
                  colSums(gr_df * dgam10),
                  colSums(gr_df * dgam11),
                  colSums(gr_df * g_tau1),
                  colSums(gr_df * g_tau2)), 9, 9)
    
    score <- matrix(c(sum(g_w), rep(0, 8)), 1, 9)
    
    gr_df2 <- data.frame(dbeta10, dbeta11, dbeta20, dbeta21, dgam10, dgam11, g_tau1, g_tau2, g_w)
    I2 <- matrix(c(colSums(gr_df * dbeta10),
                   colSums(gr_df * dbeta11),
                   colSums(gr_df * dbeta20),
                   colSums(gr_df * dbeta21),
                   colSums(gr_df * dgam10),
                   colSums(gr_df * dgam11),
                   colSums(gr_df * g_tau1),
                   colSums(gr_df * g_tau2),
                   colSums(gr_df * g_w)), 9, 9)
    
    
    sing <- ginv(I)
    
    if (class(sing)[[1]] == "try-error") {
      next
    } else {
      
    }
    result[i] <- score %*% ginv(I) %*% t(score)
    V = I2[9,9] - I2[9,1:8] %*% ginv(I2[1:8,1:8]) %*% as.matrix(I2[1:8,9])
    
    result2[i] <- sum(g_w) / sqrt(V)
  }
  # }
  # return(result)
  return(list(result, result2))
}

result_omega <- parLapply(cl, X = 1:45, please)

## 결과
n_grid <- c(100, 200, 500)
tau2_grid <- c(0.4, 0.6, 0.8, 1.0, 1.2)
g2_grid <- c(-0.997, -0.186, 0.794)

param_grid <- expand.grid(g2_grid, tau2_grid, n_grid)

qval = c(0.05, 0.1)
upp = qchisq(qval, 1, lower.tail = F)
low = qchisq(qval, 1)

qval = c(0.05, 0.1)/2
upp = qnorm(qval, 1, lower.tail = F)
low = qnorm(qval, 1)

# 스코어 결과
for(ii in 1:45){
  cat("n =", param_grid[ii, 3], "tau2 =", param_grid[ii, 2], "phi =", param_grid[ii, 1], '\n')
  for(qq in 1:2){
    cat("true q = ", qval[qq], "estimated prob =", round(mean(result_omega[[ii]][[1]] >= upp[qq], na.rm = T), 3), '\n')
    # cat("true q = ", qval[qq]*2, "estimated prob =", round(mean(result_omega[[ii]] >= upp[qq] | result_omega[[ii]] <= low[qq], na.rm = T), 3), '\n')
  }
}



# 스코어 결과
for(ii in 1:45){
  cat("n =", param_grid[ii, 3], "tau2 =", param_grid[ii, 2], "phi =", param_grid[ii, 1], '\n')
  for(qq in 1:2){
    # cat("true q = ", qval[qq], "estimated prob =", round(mean(result_omega[[ii]][[1]] >= upp[qq], na.rm = T), 3), '\n')
    cat("true q = ", qval[qq]*2, "estimated prob =", round(mean(result_omega[[ii]][[1]] >= upp[qq] | result_omega[[ii]][[1]] <= low[qq], na.rm = T), 3), '\n')
  }
}
















