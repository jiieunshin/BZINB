rm(list = ls())

library(parallel)
library(dplyr)
library(MASS)

# numCores <- detectCores() - 1
cl <- makeCluster(45)
result_list = list()
boot_list = list()

please <- function(num) {
  # for(num in 2:15){
  library(dplyr)
  library(MASS)
  RBZI <- function(n, phi, b10, b11, b20, b21, tau1, tau2, w, grid = 20) {
    
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
    
    mu1 <- exp(x1 %*% c(b10, b11));
    mu2 <- exp(x2 %*% c(b20, b21));
    
    c1 <- (1 + d * mu1 * tau1)^{-1/tau1};
    c2 <- (1 + d * mu2 * tau2)^{-1/tau2};
    
    y <- matrix(0, nrow = n, ncol = 2)
    
    for (i in 1:n) {
      PX <- cbind(dnbinom(0:(grid - 1), mu = mu1[i], size = 1/tau1), dnbinom(0:(grid - 1), mu = mu2[i], size = 1/tau2))
      PTMPY <- cbind(exp(-(0:(grid - 1))) - (1 + d * mu1[i] * tau1)^{-1/tau1}, exp(-(0:(grid - 1))) - (1 + d * mu2[i] * tau2)^{-1/tau2})
      
      TXY = matrix(PX[, 1], nrow = grid, ncol = 1) %*% matrix(PX[, 2], nrow = 1, ncol = grid)  * 
        (1 + w * matrix(PTMPY[, 1], nrow = grid, ncol = 1) %*% matrix(PTMPY[, 2], nrow = 1, ncol = grid))
      
      SXY = matrix(cumsum(TXY %>% t() %>% as.vector()), nrow = 20, ncol = 20, byrow = TRUE)
      
      # RAN_N <- runif(1)
      
      # if (RAN_N < phi) {
      #   y[i, ] <- c(0, 0)
      # } else {
      RAN_A <- runif(1)
      mat = RAN_A < SXY
      find_idx(mat)
      y[i, ] <- find_idx(mat)
      # }
    }
    return(list(y = y, x1 = x1, x2 = x2, c1 = c1, c2 = c2))
  }
  
  logL <- function(param) {
    beta1 <- param[c(1, 2)]
    beta2 <- param[c(3, 4)]
    tau <- param[c(5, 6)] + 1e-10
    w <- param[7]
    
    d <- 1 - exp(-1)
    y <- sample_data$y; x1 <- sample_data$x1; x2 <- sample_data$x2;
    
    mu1 <- c(exp(x1 %*% beta1)); mu2 <- c(exp(x2 %*% beta2))
    c <- cbind((1 + d * mu1 * tau[1])^{-1/tau[1]}, (1 + d * mu2 * tau[2])^{-1/tau[2]})
    
    v1 <- dnbinom(x = y[, 1], mu = mu1, size = 1/tau[1])
    v2 <- dnbinom(x = y[, 2], mu = mu2, size = 1/tau[2])
    v3 <- (1 + w * (exp(-y[, 1]) - c[, 1]) * (exp(-y[, 2]) - c[, 2]))
    
    v1 <- ifelse(v1 <= 0, log(1e-15), log(v1))
    v2 <- ifelse(v2 <= 0, log(1e-15), log(v2))
    v3 <- ifelse(v3 <= 0, log(1e-15), log(v3))
    
    logbp <- v1 + v2 + v3
    return( -sum(logbp))
  }
  
  
  ## Iteration
  # i = 1
  iteration = 1000
  result <- result2 <- c()
  modify_result <- c()
  
  n_grid <- c(100, 200, 500)
  tau_grid <- c(0.4, 0.6, 0.8, 1.0, 1.2)
  w_grid <- c(-1, 0, 1)
  
  param_grid <- expand.grid(w_grid, tau_grid, n_grid)
  # param_result <- matrix(0, 100, 7)
  
  for(i in 1:iteration) {
    set.seed(i)
    sample_data <- RBZI(n = param_grid[num, 3], 0, .2, .4, .4, .8, 0.8, param_grid[num, 2], param_grid[num, 1])
    
    
    one <- try(optim(par = rep(0, 7), lower = c(rep(-Inf, 4), 1e-10, 1e-10, -Inf), logL, method = "L-BFGS-B"), silent = TRUE)
    
    if (class(one) == "try-error") {
      next
    } else {
      param <- one$par
      # param_result[i,] <- param
    }
    
    beta_i <- param[c(1, 3)]
    beta_c <- param[c(2, 4)]
    tau <- param[(c(5, 6))]
    w <- param[7]
    
    d <- 1 - exp(-1)
    
    y <- sample_data$y; x <- cbind(sample_data$x1[, 2], sample_data$x2[, 2])
    mu <- cbind(exp(sample_data$x1 %*% c(beta_i[1], beta_c[1])), exp(sample_data$x2 %*% c(beta_i[2], beta_c[2])))
    c <- cbind((1 + d * mu[, 1] * tau[1])^{-1/tau[1]}, (1 + d * mu[, 2] * tau[2])^{-1/tau[2]})
    
    Ind <- (y[, 1] == 0 & y[, 2] == 0)
    
    
    # a11 <- (tau[1]^{-2} * log(1 + tau[1] * mu[, 1]) - tau[1]^{-1} * mu[, 1] / (1 + tau[1] * mu[, 1]) - w * (1 - c[, 2]) * (1 + d * tau[1] * mu[, 1])^{-1/tau[1]} * (tau[1]^{-2} * log(1 + d * tau[1] * mu[, 1]) - d * tau[1]^{-1} * mu[, 1] / (1 + d * tau[1] * mu[, 1])) / (1 + w * (1 - c[, 1]) * (1 - c[, 2]))) * (y[, 1] == 0 & y[, 2] == 0)
    # a12 <- (- w * (exp(-y[, 2]) - c[, 2]) * (1 + d * tau[1] * mu[, 1])^{-1/tau[1]} * (tau[1]^{-2} * log(1 + d * tau[1] * mu[, 1]) - d * tau[1]^{-1} * mu[, 1] / (1 + d * tau[1] * mu[, 1])) / (1 + w * (exp(-y[, 1]) - c[, 1]) * (exp(-y[, 2]) - c[, 2]))) * !(y[, 1] == 0 & y[, 2] == 0)
    # a13 <- (tau[1]^{-2} * log(1 + tau[1] * mu[, 1]) - sapply(1:nrow(y), function(i){
    #   sv = sapply(0:(y[i, 1] - 1), function(v){tau[1]^{-2} / (v + tau[1]^{-1})})
    #   return(sum(sv))
    # }) + (y[, 1] - mu[, 1]) * tau[1]^{-1} / (1 + tau[1] * mu[, 1])) * (y[, 1] > 0 & y[, 2] == 0)
    # a14 <- (tau[1]^{-2} * log(1 + tau[1] * mu[, 1]) - tau[1]^{-1} * mu[, 1] / (1 + tau[1] * mu[, 1])) * (y[, 1] == 0 & y[, 2] > 0)
    # 
    # a21 <- (tau[2]^{-2} * log(1 + tau[2] * mu[, 2]) - tau[2]^{-1} * mu[, 2] / (1 + tau[2] * mu[, 2]) - w * (1 - c[, 1]) * (1 + d * tau[2] * mu[, 2])^{-1/tau[2]} * (tau[2]^{-2} * log(1 + d * tau[2] * mu[, 2]) - d * tau[2]^{-1} * mu[, 2] / (1 + d * tau[2] * mu[, 2])) / (1 + w * (1 - c[, 1]) * (1 - c[, 2]))) * (y[, 1] == 0 & y[, 2] == 0)
    # a22 <- (- w * (exp(-y[, 1]) - c[, 1]) * (1 + d * tau[2] * mu[, 2])^{-1/tau[2]} * (tau[2]^{-2} * log(1 + d * tau[2] * mu[, 2]) - d * tau[2]^{-1} * mu[, 2] / (1 + d * tau[2] * mu[, 2])) / (1 + w * (exp(-y[, 1]) - c[, 1]) * (exp(-y[, 2]) - c[, 2]))) * !(y[, 1] == 0 & y[, 2] == 0)
    # a23 <- (tau[2]^{-2} * log(1 + tau[2] * mu[, 2]) - sapply(1:nrow(y), function(i){
    #   sv = sapply(0:(y[i, 2] - 1), function(v){tau[2]^{-2} / (v + tau[2]^{-1})})
    #   return(sum(sv))
    # }) + (y[, 2] - mu[, 2]) * tau[2]^{-1} / (1 + tau[2] * mu[, 2])) * (y[, 1] == 0 & y[, 2] > 0)
    # a24 <- (tau[2]^{-2} * log(1 + tau[2] * mu[, 2]) - tau[2]^{-1} * mu[, 2] / (1 + tau[2] * mu[, 2])) * (y[, 1] > 0 & y[, 2] == 0)
    
    # g_tau1 <- (a11 + a12 + a13 + a14)
    # g_tau2 <- (a21 + a22 + a23 + a24)
    
    g_b10 <- (y[, 1] - mu[, 1]) / (1 + tau[1] * mu[, 1]) + w * d * mu[, 1] * (exp(-y[, 2]) - c[, 2]) * (1 + d * tau[1] * mu[, 1])^{-1/tau[1] - 1} / (1 + w * (exp(-y[, 1]) - c[, 1]) * (exp(-y[, 2]) - c[, 2]))
    g_b11 <- (y[, 1] - mu[, 1]) * x[, 1] / (1 + tau[1] * mu[, 1]) + w * d * mu[, 1] * x[, 1] * (exp(-y[, 2]) - c[, 2]) * (1 + d * tau[1] * mu[, 1])^{-1/tau[1] - 1} / (1 + w * (exp(-y[, 1]) - c[, 1]) * (exp(-y[, 2]) - c[, 2]))
    g_b20 <- (y[, 2] - mu[, 2]) / (1 + tau[2] * mu[, 2]) + w * d * mu[, 2] * (exp(-y[, 1]) - c[, 1]) * (1 + d * tau[2] * mu[, 2])^{-1/tau[2] - 1} / (1 + w * (exp(-y[, 1]) - c[, 1]) * (exp(-y[, 2]) - c[, 2]))
    g_b21 <- (y[, 2] - mu[, 2]) * x[, 2] / (1 + tau[2] * mu[, 2]) + w * d * mu[, 2] * x[, 2] * (exp(-y[, 1]) - c[, 1]) * (1 + d * tau[2] * mu[, 2])^{-1/tau[2] - 1} / (1 + w * (exp(-y[, 1]) - c[, 1]) * (exp(-y[, 2]) - c[, 2]))

    g_w <- (exp(-y[, 1]) - c[, 1]) *  (exp(-y[, 2]) - c[, 2]) / (1 + w * (exp(-y[, 1]) - c[, 1]) *  (exp(-y[, 2]) - c[, 2]))
    
    g_psi <- (1 - (1 + tau[2] * mu[, 2])^{-1/tau[2]} * (1 + tau[1] * mu[, 1])^{-1/tau[1]} * (1 + w * (1 - c[, 1]) * (1 - c[, 2]))) / ((1 + tau[2] * mu[, 2])^{-1/tau[2]} * (1 + tau[1] * mu[, 1])^{-1/tau[1]} * (1 + w * (1 - c[, 1]) * (1 - c[, 2]))) * Ind + (-1) * !Ind

    aa1 <- (tau[1]^{-2} * log(1 + tau[1] * mu[, 1]) - tau[1]^{-1} * mu[, 1] / (1 + tau[1] * mu[, 1]))
    
    aa2 <- ((1 + tau[2] * mu[, 2])^{-1/tau[2]} * (1 + tau[1] * mu[, 1])^{-1/tau[1]} * aa1 ) / ((1 + tau[2] * mu[, 2])^{-1/tau[2]} * (1 + tau[1] * mu[, 1])^{-1/tau[1]}) * Ind
    
    aa3 = sapply(1:nrow(y), function(m){
      if(y[m, 1] == 0){
        val = 0
      } else{
        sv = sapply(1:y[m, 1], function(v){ (y[m, 1] - v)/(1 - v * tau[1] + y[m, 1] * tau[1]) })
        val = sum(sv)
      }
      return(val) 
    })
    
    aa4 <- ( aa3 - y[, 1] * mu[, 1] / (1 + tau[1] * mu[, 1]) - mu[, 1] * tau[1]^{-1} / (1 + tau[1] * mu[, 1]) + tau[1]^{-2} * log(1 + tau[1] * mu[, 1]) ) * !Ind
    
    # 
    ab1 <- (tau[2]^{-2} * log(1 + tau[2] * mu[, 2]) - tau[2]^{-1} * mu[, 2] / (1 + tau[2] * mu[, 2]))
    
    ab2 <- ((1 + tau[2] * mu[, 2])^{-1/tau[2]} * (1 + tau[1] * mu[, 1])^{-1/tau[1]} * ab1 ) / ((1 + tau[2] * mu[, 2])^{-1/tau[2]} * (1 + tau[1] * mu[, 1])^{-1/tau[1]}) * Ind
    
    ab3 = sapply(1:nrow(y), function(m){
      if(y[m, 2] == 0){
        val = 0
      } else{
        sv = sapply(1:y[m, 2], function(v){ (y[m, 2] - v)/(1 - v * tau[2] + y[m, 2] * tau[2]) })
        val = sum(sv)
      }
      return(val) 
    })
    
    ab4 <- ( ab3 - y[, 2] * mu[, 2] / (1 + tau[2] * mu[, 2]) - mu[, 2] * tau[2]^{-1} / (1 + tau[2] * mu[, 2]) + tau[2]^{-2} * log(1 + tau[2] * mu[, 2]) ) * !Ind
    
    g_tau1 <- (aa2 + aa4)
    g_tau2 <- (ab2 + ab4)
    
    gr_df <- data.frame(g_psi, g_b10, g_b11, g_b20, g_b21, g_tau1, g_tau2, g_w)
    
    I <- matrix(c(colSums(gr_df * g_psi),
                  colSums(gr_df * g_b10),
                  colSums(gr_df * g_b11),
                  colSums(gr_df * g_b20),
                  colSums(gr_df * g_b21),
                  colSums(gr_df * g_tau1),
                  colSums(gr_df * g_tau2),
                  colSums(gr_df * g_w)), 8, 8)
    
    score <- matrix(c(sum(g_psi), rep(0, 7)), 1, 8)

    gr_df2 <- data.frame(g_b10, g_b11, g_b20, g_b21, g_tau1, g_tau2, g_w, g_psi)
    I2 <- matrix(c(colSums(gr_df2 * g_b10),
                  colSums(gr_df2 * g_b11),
                  colSums(gr_df2 * g_b20),
                  colSums(gr_df2 * g_b21),
                  colSums(gr_df2 * g_tau1),
                  colSums(gr_df2 * g_tau2),
                  colSums(gr_df2 * g_w),
                  colSums(gr_df2 * g_psi)), 8, 8)
    
    
    sing <- ginv(I)
    
    if (class(sing)[[1]] == "try-error") {
      next
    } else {
      
    }
    result[i] <- score %*% ginv(I) %*% t(score)
    V = I2[8,8] - I2[8,1:7] %*% ginv(I2[1:7,1:7]) %*% as.matrix(I2[1:7,8])
    
    result2[i] <- sum(g_psi) / sqrt(V)
  }
  return(list(result, result2))
  # return(result)
}



result_list1000 <- parLapply(cl, X = 1:45, please)

stopCluster(cl)

## 결과
n_grid <- c(100, 200, 500)
tau_grid <- c(0.4, 0.6, 0.8, 1.0, 1.2)
w_grid <- c(-1, 0, 1)

param_grid <- expand.grid(w_grid, tau_grid, n_grid)

qval = c(0.05, 0.1)
upp = qchisq(qval, 1, lower.tail = F)
low = qchisq(qval, 1)

qval = c(0.05, 0.1)/2
upp = qnorm(qval, lower.tail = F)
low = qnorm(qval)


# 스코어 결과
for(ii in 1:45){
  cat("n =", param_grid[ii, 3], "tau2 =", param_grid[ii, 2], "w =", param_grid[ii, 1], '\n')
  for(qq in 1:2){
    cat("true q = ", qval[qq], "estimated prob =", round(mean(result_list1000[[ii]][[1]]>= upp[qq], na.rm = T), 3), '\n')
      # cat("true q = ", qval[qq]*2, "estimated prob =", round(mean(result_list1000[[ii]][[1]]>= upp[qq] | result_list1000[[ii]][[1]]<= low[qq], na.rm = T), 3), '\n')
  }
}

# 스코어 결과
for(ii in 1:45){
  cat("n =", param_grid[ii, 3], "tau2 =", param_grid[ii, 2], "w =", param_grid[ii, 1], '\n')
  for(qq in 1:2){
    # cat("true q = ", qval[qq], "estimated prob =", round(mean(result_list1000[[ii]][[1]]>= upp[qq], na.rm = T), 3), '\n')
    cat("true q = ", qval[qq]*2, "estimated prob =", round(mean(result_list1000[[ii]][[2]]>= upp[qq] | result_list1000[[ii]][[2]]<= low[qq], na.rm = T), 3), '\n')
  }
}

#histogram
id = c(16, 19, 22, 25, 28, 31, 34, 37, 40, 43)
par(mfrow=c(2,5))
for(ii in id){
  y = result_list10000[[ii]][[1]]
  qqplot(qchisq(ppoints(500), df=1), y, main = paste0("n=", param_grid[ii,3], ", tau2=", param_grid[ii,2], ", w=", param_grid[ii,1]), ylab = 'Simulated score tests', xlab = 'chi-square(1) quantiles')
  qqline(y, distribution = function(p) qchisq(p, df = 1), probs = c(0.001, 0.999), col = 2)
  
}

id = c(16, 19, 22, 25, 28, 31, 34, 37, 40, 43)
par(mfrow=c(2,5))
for(ii in id){
  y = result_list10000[[ii]][[3]]
  qqplot(qnorm(ppoints(500)), y, main = paste0("n=", param_grid[ii,3], ", tau2=", param_grid[ii,2], ", w=", param_grid[ii,1]), ylab = 'Simulated score tests', xlab = 'N(0,1) quantiles')
  qqline(y, distribution = function(p) qnorm(p), probs = c(0.001, 0.999), col = 2)
  
}



save(result_list1000, file = "~/phi_score1000.RData")
load("~/phi_score.RData")
