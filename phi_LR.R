rm(list = ls())

library(parallel)
library(dplyr)

numCores <- detectCores() - 1
cl <- makeCluster(numCores)

result = list()

phi_LR1000 <-  parLapply(cl = cl, 1:45, function(j){
  
  library(dplyr)
  library(MASS)
  n_grid <- c(100, 200, 500)
  tau_grid <- c(0.2, 0.4, 0.8, 1.2, 1.6)
  w_grid <- c(-1, 0, 1)
  
  param_grid <- expand.grid(w_grid, tau_grid, n_grid)
  
  ### random number generator
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
    
    c1 <- (1 + d * tau1 * mu1)^{-1/tau1};
    c2 <- (1 + d * tau2 * mu2)^{-1/tau2};
    
    y <- matrix(0, nrow = n, ncol = 2)
    
    for (i in 1:n) {
      PX <- cbind(dnbinom(0:(grid - 1), mu = mu1[i], size = 1/tau1), dnbinom(0:(grid - 1), mu = mu2[i], size = 1/tau2))
      PTMPY <- cbind(exp(-(0:(grid - 1))) - c1[i], exp(-(0:(grid - 1))) - c2[i])
      
      TXY = matrix(PX[, 1], nrow = grid, ncol = 1) %*% matrix(PX[, 2], nrow = 1, ncol = grid)  * 
        (1 + w * matrix(PTMPY[, 1], nrow = grid, ncol = 1) %*% matrix(PTMPY[, 2], nrow = 1, ncol = grid))
      
      SXY = matrix(cumsum(TXY %>% t() %>% as.vector()), nrow = 20, ncol = 20, byrow = TRUE)
      
      RAN_N <- runif(1)
      
      if (RAN_N < phi) {
        y[i, ] <- c(0, 0)
      } else {
        RAN_A <- runif(1)
        mat = RAN_A < SXY
        find_idx(mat)
        y[i, ] <- find_idx(mat)
      }
    }
    return(list(y = y, x1 = x1, x2 = x2, c1 = c1, c2 = c2))
  }
  
  
  ### log likelihood
  logL0 <- function(param) {
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
  
  ### log likelihood
  logL1 <- function(param) {
    phi <- param[1]
    beta1 <- param[c(2, 3)]
    beta2 <- param[c(4, 5)]
    tau <- param[c(6, 7)]
    w <- param[8]
    
    d <- 1 - exp(-1)
    
    y <- sample_data$y
    x1 <- sample_data$x1
    x2 <- sample_data$x2
    
    mu1 <- exp(x1 %*% beta1)
    mu2 <- exp(x2 %*% beta2)
    
    c1 <- (1 + d * mu1 * tau[1])^{-1/tau[1]}
    c2 <- (1 + d * mu2 * tau[2])^{-1/tau[2]}
    
    ind <- (y[, 1] == 0 & y[, 2] == 0)
    
    v1 <- phi + (1 - phi) * (1 + tau[1] * mu1)^{-1/tau[1]} * (1 + tau[2] * mu2)^{-1/tau[2]} * (1 + w * (1 - c1) * (1 - c2))
    v2 <- (1 - phi) * dnbinom(x = y[, 1], mu = mu1, size = 1/tau[1]) * dnbinom(x = y[, 2], mu = mu2, size = 1/tau[2]) * (1 + w * (exp(-y[, 1]) - c1) * (exp(-y[, 2]) - c2))    
    
    v1 <- log(v1[ind])
    v2 <- log(v2[!ind])
    
    return(-(sum(v1) + sum(v2)))
  }
  
  # iteration
  result_score <- c()
  
  i <- 1
  k <- 1
  oo <- list()
  tt <- list()
  ### parallel loop
  while(1) {
    set.seed(k)
    
    k <- k + 1
    sample_data <- RBZI(n = param_grid[j, 3], 0, .2, .4, .4, .8, param_grid[j, 2], 1, param_grid[j, 1])
    
    one <- try(optim(par = rep(.1, 8), fn = logL1, method = "L", lower = c(rep(-Inf, 5),  1e-10, 1e-10, -Inf)), silent = TRUE)
    two <- try(optim(par = rep(.1, 7), fn = logL0, method = "L", lower = c(rep(-Inf, 4),  1e-10, 1e-10, -Inf)), silent = TRUE)
    
    if (class(one) == "try-error" | class(two) == "try-error") {
      
      next
    }
    
    print(i)
    
    result_score[i] <- 2 * (two$value - one$value)
    tt[[i]] <- two$par
    oo[[i]] <- one$par
    i <- i + 1
    
    if (i == 1001) {
      break
    }
    
  }
  
  return(result_score)
  
})


### 결과 정리 ###
# 파라미터들
n_grid <- c(100, 200, 500)
tau_grid <- c(1.2, 1.4, 1.6, 1.8, 2)
w_grid <- c(-1, 0, 1)

param_grid <- expand.grid(w_grid, tau_grid, n_grid)

qval = c(0.1, 0.05)/2

upp = qchisq(qval, 1, lower.tail = F)
low = qchisq(1-qval, 1, lower.tail = F)

# 양측
for(ii in 1:45){
  cat('-------- \n')
  cat("n =", param_grid[ii, 3], "tau =", param_grid[ii, 2], "w =", param_grid[ii, 1], '\n')
  for(qq in 1:2){
    # cat("true q = ", qval[qq], "estimated prob =", round((mean(phi_LR1000[[ii]] >= upp[qq], na.rm = T)), 3), '\n')
    cat("true q = ", qval[qq]*2, "estimated prob =", round((mean(phi_LR1000[[ii]] >= upp[qq] | phi_LR1000[[ii]] <= low[qq], na.rm = T)), 3), '\n')
  }
}


#histogram
id = c(16, 19, 22, 25, 28, 31, 34, 37, 40, 43)
par(mfrow=c(2,5))
for(ii in id){
  y = phi_LR10000[[ii]]
  qqplot(qchisq(ppoints(500), df=1), y, main = paste0("n=", param_grid[ii,3], ", tau2=", param_grid[ii,2], ", w=", param_grid[ii,1]), ylab = 'Simulated score tests', xlab = 'chi-square(1) quantiles')
  qqline(y, distribution = function(p) qchisq(p, df = 1), probs = c(0.001, 0.999), col = 2)
  
}



###################