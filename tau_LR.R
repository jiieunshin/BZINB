rm(list = ls())
library(parallel)

n_cores = detectCores() - 1
cl <- makeCluster(45)
clusterEvalQ(cl, c("MASS", "dplyr"))
j=1
tau_LR10000 = parLapply(cl = cl, 1:45, function(j){
  library(dplyr)
  library(MASS)
  n_grid <- c(100, 200, 500)
  g2_grid <- c(-0.997, -0.186, 0.352, 0.794, 1.20)
  w_grid <- c(-1, 0, 1)
  
  param_grid <- t(expand.grid(w_grid, g2_grid, n_grid))
  
  RBZI <- function(n, b10, b11, b20, b21, gam10, gam11, w, tau1, tau2, grid = 20) {
    
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
  
  logL0 <- function(param) {
    beta1 <- param[c(1, 2)]
    beta2 <- param[c(3, 4)]
    gamma <- param[c(5, 6)]
    w <- param[7]
    
    d <- 1 - exp(-1)
    y <- sample_data$y; x1 <- sample_data$x1; x2 <- sample_data$x2; z <- sample_data$z
    
    mu1 <- c(exp(x1 %*% beta1)); mu2 <- c(exp(x2 %*% beta2))
    c1 <- exp(mu1 * exp(-1) - mu1); c2 <- exp(mu2 * exp(-1) - mu2)
    
    phi <- c(exp(z %*% gamma) / (1 + exp(z %*% gamma)))
    zg <- c(exp(z %*% gamma))
    ind <- (y[, 1] == 0 & y[, 2] == 0)
    
    v1 <- sum(log(exp(z %*% gamma) + dpois(x = 0, lambda = mu1) * dpois(x = 0, lambda = mu2) * (1 + w * (1 - c1) * (1 - c2))) * ind)
    
    v21 <- (dpois(x = y[, 1], lambda = mu1, log = T) + dpois(x = y[, 2], lambda = mu2, log = T)) * !ind
    v22 <- c((1 + w * (exp(-y[, 1]) - c1) * (exp(-y[, 2]) - c2)) * !ind)
    v22 <- ifelse(v22 <= 0, log(1e-15), log(v22)) * !ind
    v2 <- (v21 + v22)*!ind
    logbp <- sum(v1) + sum(v2) - sum(log(1+exp(z %*% gamma)))
    
    return(-logbp)
  }
  
  logL1 <- function(param) {
    beta1 <- param[c(1, 2)]
    beta2 <- param[c(3, 4)]
    gamma <- param[c(5, 6)]
    tau <- param[c(7, 8)] + 1e-10
    w <- param[9]
    # print(tau)
    d <- 1 - exp(-1)
    y <- sample_data$y; x1 <- sample_data$x1; x2 <- sample_data$x2; z <- sample_data$z
    # tau = ifelse(tau < 0, 1e-8, tau)
    
    mu1 <- c(exp(x1 %*% beta1)); mu2 <- c(exp(x2 %*% beta2))
    c1 <- (1 + d*tau[1]*mu1)^{-1/tau[1]}; c2 <- (1 + d*tau[2]*mu2)^{-1/tau[2]}
    
    phi <- c(exp(z %*% gamma) / (1 + exp(z %*% gamma)))
    zg <- c(exp(z %*% gamma))
    ind <- (y[, 1] == 0 & y[, 2] == 0)
    
    # tau = ifelse(tau < 0, 1e-8, tau)
    
    v1 <- sum(log(exp(z %*% gamma) + dnbinom(x = 0, mu = mu1, size = 1/tau[1]) * dnbinom(x = 0, mu = mu2, size = 1/tau[2]) * (1 + w * (1 - c1) * (1 - c2))) * ind)
    # print(tau)
    v21 <- (dnbinom(x = y[, 1], mu = mu1, size = 1/tau[1], log = T) + dnbinom(x = y[, 2], mu = mu2, size = 1/tau[2], log = T)) * !ind
    v22 <- c((1 + w * (exp(-y[, 1]) - c1) * (exp(-y[, 2]) - c2)) * !ind)
    v22 <- ifelse(v22 <= 0, log(1e-15), log(v22)) * !ind
    v2 <- (v21 + v22)*!ind
    logbp <- sum(v1) + sum(v2) - sum(log(1+exp(z %*% gamma)))
    
    
    return(-logbp)
  }
  
  # 돌리기
  iteration = 10000

  result <- c()
  param_est <- matrix(0, iteration, 7)
  for(i in 1:iteration) {
    set.seed(i)
    print(i)
    sample_data <- RBZI(n = param_grid[3, j], .2, .4, .4, .8, -1.2, param_grid[2, j], param_grid[1, j], tau1 = 1e-10, tau2 = 1e-10)

    c1 = sample_data$c1
    c2 = sample_data$c2
    
    w1 = -1/max(c1*c2, 1-c1-c2+c1*c2)  #### 추가함
    w2 = 1/max(c1*(1-c1), c2*(1-c2))
    
    if ((class(try(optim(par = rep(0, 7), logL0, upper = rep(10, 7), method = "L-BFGS-B"), silent = TRUE)) == "try-error") |
        (class(try(optim(par = rep(0, 9), logL1, upper = rep(10, 9), lower = c(rep(-Inf, 6),rep(1e-10, 2), -10), method = "L-BFGS-B"), silent = TRUE)) == "try-error")) {
      
      L0 = NA
      L1 = NA

    } else{
      L0 <- optim(par = rep(0, 7), logL0, upper = rep(10, 7), method = "L-BFGS-B")$value
      L1 <- optim(par = rep(0, 9), logL1, upper = rep(10, 9), lower = c(rep(-Inf, 6),rep(1e-10, 2), -10), method = "L-BFGS-B")$value
      param_est[i,] <- optim(par = rep(0, 7), logL0, upper = rep(10, 7), method = "L-BFGS-B")$par
    }

    result[i] <- 2*(L0-L1)
    
  } # end while 
  return(list('score' = result))
  # return(result)
})  # end func

for(qq in 1:4){
  cat("true q = ", qval[qq], "estimated prob =", round(sum(result[1:i]>= upp[qq], na.rm = T)/length(result[1:i]), 4), '\n')
}


############################

## n=500일때 w가 마이너스로 쏠리는지 확인할 것
n_grid <- c(100, 200, 500)
phi_grid <- c(0.1, 0.2, 0.3, 0.4, 0.5)
w_grid <- c(-1.5, 0,  1.5)

param_grid <- t(expand.grid(w_grid, phi_grid, n_grid))

qval = c(0.1, 0.05)
upp = qchisq(qval, 2, lower.tail = F) / 4 + qchisq(qval, 1, lower.tail = F)/2
upp = qchisq(qval, 2, lower.tail = F)

hist(result, breaks = 200)
abline(v = upp, col = "blue")
head(param_grid)

# 단측
for(ii in 1:45){
    cat("n =", param_grid[3, ii], "gam2 =", param_grid[2, ii], "w =", param_grid[1, ii], '\n')
    for(qq in 1:2){
      cat("true q = ", qval[qq], "estimated prob =", round(mean(tau_LR10000[[ii]]$score>= upp[qq], na.rm = T), 3), '\n')
    }
}

id = c(16, 19, 22, 25, 28, 31, 34, 37, 40, 43)
par(mfrow=c(2,5))
for(ii in id){
  y = tau_LR10000[[ii]]$score
  qqplot(qchisq(ppoints(500), df=1), y, main = paste0("n=", param_grid[3,ii], ", tau2=", param_grid[2,ii], ", w=", param_grid[1,ii]), ylab = 'Simulated score tests', xlab = 'chi-square(1) quantiles')
  qqline(y, distribution = function(p) qchisq(p, df = 1), probs = c(0.001, 0.999), col = 2)
}

##
save(out, file = '~/tau_LR.RData')
load('C:/Users/jieun/Dropbox/BZINB/zero-zero결과/dispersion.RData')
