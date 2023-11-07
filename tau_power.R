rm(list = ls())

library(parallel)
library(dplyr)
library(MASS)

# numCores <- detectCores() - 1
cl <- makeCluster(30)
result_list = list()
boot_list = list()

please <- function(num) {
  # for(num in 2:15){
  library(dplyr)
  library(MASS)
  
  
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
  
  
  logL_nb <- function(param) {
    beta1 <- param[c(1, 2)]
    beta2 <- param[c(3, 4)]
    gamma <- param[c(5, 6)]
    tau <- param[c(7, 8)]
    w <- param[9]
    
    d <- 1 - exp(-1)
    y <- sample_data$y; x1 <- sample_data$x1; x2 <- sample_data$x2; z <- sample_data$z
    
    mu1 <- c(exp(x1 %*% beta1)); mu2 <- c(exp(x2 %*% beta2))
    c1 <- (1 + d * mu1 * tau[1])^{-1/tau[1]}; c2 <- (1 + d * mu2 * tau[2])^{-1/tau[2]}
    
    phi <- c(exp(z %*% gamma) / (1 + exp(z %*% gamma)))
    zg <- c(exp(z %*% gamma))
    ind <- (y[, 1] == 0 & y[, 2] == 0)
    
    v1 <- sum(log(exp(z %*% gamma) + dnbinom(0, mu = mu1, size = 1/tau[1]) * dnbinom(0, mu = mu2, size = 1/tau[2]) * (1 + w * (1 - c1) * (1 - c2))) * ind)
    
    v21 <- (dnbinom(x = y[, 1], mu = mu1, size = 1/tau[1], log = T) + dnbinom(x = y[, 2], mu = mu2, size = 1/tau[2], log = T)) * !ind
    v22 <- c((1 + w * (exp(-y[, 1]) - c1) * (exp(-y[, 2]) - c2)) * !ind)
    v22 <- ifelse(v22 <= 0, log(1e-15), log(v22)) * !ind
    v2 <- (v21 + v22)*!ind
    logbp <- sum(v1) + sum(v2) - sum(log(1+exp(z %*% gamma)))
    
    return(-logbp)
  }
  
  ## Iteration
  # i = 1
  iteration = 1000
  result <- c()
  boot_cric <- matrix(0, 1000, 4)
  
  n_grid <- c(100, 200)
  # g2_grid <- c(0.2, 0.8, 1.2)
  tau_grid <- c(1e-10, 0.1, 0.2, 0.4, 0.8)
  w_grid <- c(-1.5, 0, 1.5)
  
  param_grid <- expand.grid(w_grid, tau_grid, n_grid)
  
  for(i in 1:iteration) {
    set.seed(i)
    # sample_data <- RBPZI(n = param_grid[num, 3], .2, .4, .4, .8, -1.2, param_grid[num, 2], param_grid[num, 1])
    sample_data <- RBNBZI(n = param_grid[num, 3], .2, .4, .4, .8, -1.2, 1.2, 
                          w = param_grid[num, 1], tau1 = 0.4, tau2 = param_grid[num, 2])
    
    one <- try(optim(par = rep(0, 9), fn = logL_nb, lower = c(rep(-5, 6), 1e-10, 1e-10, -10), method = "L-BFGS-B"), silent = TRUE)
    
    if (class(one) == "try-error") {
      next
    } else {
      param <- one$par
    }
    
    beta_i <- param[c(1, 3)]
    beta_c <- param[c(2, 4)]
    gam <- param[(c(5, 6))]
    tau <- param[(c(7, 8))]
    w <- param[9]
    
    d <- 1 - exp(-1)
    
    ## Data
    y <- sample_data$y;
    x <- cbind(sample_data$x1[, 2], sample_data$x2[, 2]);
    z <- sample_data$z[,2]
    
    ## Parameter
    mu <- cbind(exp(sample_data$x1 %*% c(beta_i[1], beta_c[1])), exp(sample_data$x2 %*% c(beta_i[2], beta_c[2])))
    phi <- exp(sample_data$z %*% gam)/(1+exp(sample_data$z %*% gam))
    psi <- phi / (1 - phi)
    c <- cbind((1 + d * mu[, 1] * tau[1])^{-1/tau[1]}, (1 + d * mu[, 2] * tau[2])^{-1/tau[2]})
    D <- (1 + w * (exp(-y[, 1]) - c[, 1]) * (exp(-y[, 2]) - c[, 2]))
    
    ## Indicator 
    Ind <- (y[, 1] == 0 & y[, 2] == 0)
    
    ## Gradient
    
    ## GAMMA
    dgam0 <- (exp(cbind(1, z) %*% gam) / (exp(cbind(1, z) %*% gam) + exp(-mu[, 1] - mu[, 2]) * D)) * Ind -
      (exp(cbind(1, z) %*% gam) / (1 + exp(cbind(1, z) %*% gam)))
    
    dgam1 <- (z * exp(cbind(1, z) %*% gam) / (exp(cbind(1, z) %*% gam) + exp(-mu[, 1] - mu[, 2]) * D)) * Ind -
      (z * exp(cbind(1, z) %*% gam) / (1 + exp(cbind(1, z) %*% gam)))
    
    
    ## BETA1
    dbeta10 <- (exp(-mu[, 1] - mu[, 2]) * mu[, 1] * (w * d * (1 - exp(-d * mu[, 2])) * exp(-d * mu[, 1]) - D) / (psi + exp(-mu[, 1] - mu[, 2]) * D)) * Ind +
      ((y[, 1] - mu[, 1]) + w * d * mu[, 1] * (exp(-y[, 2]) - exp(-d * mu[, 2])) * exp(-d * mu[, 1]) / D) * !Ind
    
    dbeta11 <- (x[, 1] * exp(-mu[, 1] - mu[, 2]) * mu[, 1] * (w * d * (1 - exp(-d * mu[, 2])) * exp(-d * mu[, 1]) - D) / (psi + exp(-mu[, 1] - mu[, 2]) * D)) * Ind +
      ((y[, 1] - mu[, 1]) * x[, 1] + w * d * mu[, 1] * x[, 1] * (exp(-y[, 2]) - exp(-d * mu[, 2])) * exp(-d * mu[, 1]) / D) * !Ind
    
    
    ## BETAA2
    dbeta20 <- (exp(-mu[, 1] - mu[, 2]) * mu[, 2] * (w * d * (1 - exp(-d * mu[, 1])) * exp(-d * mu[, 2]) - D) / (psi + exp(-mu[, 1] - mu[, 2]) * D)) * Ind +
      ((y[, 2] - mu[, 2]) + w * d * mu[, 2] * (exp(-y[, 1]) - exp(-d * mu[, 1])) * exp(-d * mu[, 2]) / D) * !Ind
    
    dbeta21 <- (x[, 2] * exp(-mu[, 1] - mu[, 2]) * mu[, 2] * (w * d * (1 - exp(-d * mu[, 1])) * exp(-d * mu[, 2]) - D) / (psi + exp(-mu[, 1] - mu[, 2]) * D)) * Ind +
      ((y[, 2] - mu[, 2]) * x[, 2] + w * d * mu[, 2] * x[, 2] * (exp(-y[, 1]) - exp(-d * mu[, 1])) * exp(-d * mu[, 2]) / D) * !Ind
    
    
    ## OMEGA
    dw = (exp(-mu[, 1] - mu[, 2]) * (1 - exp(-d * mu[, 1])) * (1 - exp(-d * mu[, 2])) / (psi + exp(-mu[, 1] - mu[, 2]) * D)) * Ind +
      ((exp(-y[, 1]) - exp(-d * mu[, 1])) * (exp(-y[, 2]) - exp(-d * mu[, 2])) / D) * !Ind
    
    
    ## TAU
    a1 = sapply(1:nrow(y), function(i){ 
      sv = 0
      if(y[i, 1] > 0){
        for(v in 1:y[i, 1]){
          sv = sv + (y[i,1] - v)
        }
      }
      return (sv)
    })
    
    a2 = sapply(1:nrow(y), function(i){ 
      sv = 0
      if(y[i, 2] > 0){
        for(v in 1:y[i, 2]){
          sv = sv + (y[i,2] - v)
        }
      }
      return (sv)
    })
    
    dtau1 <- (exp(-mu[, 1] -mu[ ,2]) * D * (mu[, 1]^2 / 2) / (exp(cbind(1, z) %*% gam) + exp(-mu[, 1] - mu[, 2]) * D)) * Ind -
      (w * exp(-mu[, 1] - mu[, 2]) * exp(-d * mu[, 1]) * (1 - exp(-d * mu[, 2])) * (mu[, 1]^2 * d / 2)  / (exp(cbind(1, z) %*% gam) + exp(-mu[, 1] - mu[, 2]) * D)) * Ind +
      (a1 - mu[, 1] * y[, 1] + mu[, 1]^2/ 2 - w * (exp(-y[, 2]) - c[, 2]) * exp(-d * mu[, 1]) * (d * mu[, 1]^2 / 2) / (1 + D)) * !Ind
    
    dtau2 <- (exp(-mu[, 1] -mu[ ,2]) * D * (mu[, 2]^2 / 2) / (exp(cbind(1, z) %*% gam) + exp(-mu[, 1] - mu[, 2]) * D)) * Ind -
      (w * exp(-mu[, 1] - mu[, 2]) * exp(-d * mu[, 2]) * (1 - exp(-d * mu[, 1])) * (mu[, 2]^2 * d / 2)  / (exp(cbind(1, z) %*% gam) + exp(-mu[, 1] - mu[, 2]) * D)) * Ind +
      (a2 - mu[, 2] * y[, 2] + mu[, 2]^2/ 2 - w * (exp(-y[, 1]) - c[, 1]) * exp(-d * mu[, 2]) * (d * mu[, 2]^2 / 2) / (1 + D)) * !Ind
    
    
    
    gr_df <- data.frame(dtau1, dtau2, dgam0, dgam1, dbeta10, dbeta11, dbeta20, dbeta21, dw)
    
    I <- matrix(c(colSums(gr_df * dtau1),
                  colSums(gr_df * dtau2),
                  colSums(gr_df * dgam0),
                  colSums(gr_df * dgam1),
                  colSums(gr_df * dbeta10),
                  colSums(gr_df * dbeta11),
                  colSums(gr_df * dbeta20),
                  colSums(gr_df * dbeta21),
                  colSums(gr_df * dw)), 9, 9)
    
    
    
    score <- colSums(gr_df)
    
    sing <- try(ginv(I), silent = TRUE)
    
    if (class(sing)[[1]] == "try-error") {
      next
    } else {
      
    }
    result[i] <- t(score) %*% ginv(I) %*% score
    
    logLB <- function(paramB) {
      beta1B <- paramB[c(1, 2)]
      beta2B <- paramB[c(3, 4)]
      gammaB <- paramB[c(5, 6)]
      tauB <- paramB[c(7, 8)]
      wB <- paramB[9]
      
      d <- 1 - exp(-1)
      yB <- yB; x1 <- sample_data$x1; x2 <- sample_data$x2; z <- sample_data$z
      
      mu1B <- c(exp(x1 %*% beta1B)); mu2B <- c(exp(x2 %*% beta2B))
      c1B <- (1 + d * mu1B * tauB[1])^{-1/tauB[1]}; c2B <- (1 + d *mu2B * tauB[2])^{-1/tauB[2]}
      
      phiB <- c(exp(z %*% gammaB) / (1 + exp(z %*% gammaB)))
      zB <- c(exp(z %*% gammaB))
      indB <- (yB[, 1] == 0 & yB[, 2] == 0)
      
      v1 <- sum(log(exp(z %*% gammaB) + dnbinom(x = 0, mu = mu1B, size = 1/tauB[1]) * dnbinom(x = 0, mu = mu2B, size = 1/tauB[2]) * (1 + wB * (1 - c1B) * (1 - c2B))) * indB)
      
      v21 <- (dnbinom(x = yB[, 1], mu = mu1B, size = 1/tauB[1], log = T) + dnbinom(x = yB[, 2], mu = mu2B, size = 1/tauB[2], log = T)) * !indB
      v22 <- c((1 + wB * (exp(-yB[, 1]) - c1B) * (exp(-yB[, 2]) - c2B)) * !indB)
      v22 <- ifelse(v22 <= 0, log(1e-15), log(v22)) * !indB
      v2 <- (v21 + v22)*!indB
      logbp <- sum(v1) + sum(v2) - sum(log(1+exp(z %*% gammaB)))
      
      return(-logbp)
    }
    
    ## Bootstrap step
    B <- 1000
    # b <- 1
    K = nrow(y)
    grid = 20
    score_b <- c()
    
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
    
    for(b in 1:B) {
      yB <- matrix(0, nrow = K, ncol = 2)
      
      for (k in 1:K) {
        PX <- cbind(dnbinom(0:(grid - 1), mu = mu[, 1][k], size = 1/tau[1]), dnbinom(0:(grid - 1), mu = mu[, 2][k], size = 1/tau[2]))
        PTMPY <- cbind(exp(-(0:(grid - 1))) -(1 + d * mu[, 1][k] * tau[1])^{-1/tau[1]}, exp(-(0:(grid - 1))) - (1 + d * mu[, 2][k] * tau[1])^{-1/tau[1]})
        
        TXY = matrix(PX[, 1], nrow = grid, ncol = 1) %*% matrix(PX[, 2], nrow = 1, ncol = grid)  * 
          (1 + w * matrix(PTMPY[, 1], nrow = grid, ncol = 1) %*% matrix(PTMPY[, 2], nrow = 1, ncol = grid))
        
        SXY = matrix(cumsum(TXY %>% t() %>% as.vector()), nrow = 20, ncol = 20, byrow = TRUE)
        
        RAN_N <- runif(1)
        
        if (RAN_N < phi[k]) {
          yB[k, ] <- c(0, 0)
        } else {
          RAN_A <- runif(1)
          mat = RAN_A < SXY
          find_idx(mat)
          yB[k, ] <- find_idx(mat)
        }
      }
      paramB <- try(optim(par = rep(0, 9), logLB, lower = c(rep(-5, 6), 1e-10, 1e-10, -10), method = "L-BFGS-B")$par, silent = TRUE)
      
      if (class(paramB) == "try-error") {
        next
      } else {
        
      }
      
      beta_iB <- paramB[c(1, 3)]
      beta_cB <- paramB[c(2, 4)]
      gamB <- paramB[(c(5, 6))]
      tauB <- paramB[(c(7, 8))]
      wB <- paramB[9]
      
      ## Parameter
      muB <- cbind(exp(sample_data$x1 %*% c(beta_iB[1], beta_cB[1])), exp(sample_data$x2 %*% c(beta_iB[2], beta_cB[2])))
      phiB <- exp(sample_data$z %*% gamB)/(1+exp(sample_data$z %*% gamB))
      psiB <- phiB / (1 - phiB)
      cB <- cbind((1 + d * muB[, 1] * tauB[1])^{-1/tauB[1]}, (1 + d * muB[, 2] * tauB[2])^{-1/tauB[2]})
      DB <- (1 + wB * (exp(-yB[, 1]) - cB[, 1]) * (exp(-yB[, 2]) - cB[, 2]))
      
      ## Indicator 
      IndB <- (yB[, 1] == 0 & yB[, 2] == 0)
      
      ## Gradient
      
      ## GAMMA
      dgam0B <- (exp(cbind(1, z) %*% gamB) / (exp(cbind(1, z) %*% gamB) + exp(-muB[, 1] - muB[, 2]) * DB)) * IndB -
        (exp(cbind(1, z) %*% gamB) / (1 + exp(cbind(1, z) %*% gamB)))
      
      dgam1B <- (z * exp(cbind(1, z) %*% gamB) / (exp(cbind(1, z) %*% gamB) + exp(-muB[, 1] - muB[, 2]) * DB)) * IndB -
        (z * exp(cbind(1, z) %*% gamB) / (1 + exp(cbind(1, z) %*% gamB)))
      
      ## BETA1
      dbeta10B <- (exp(-muB[, 1] - muB[, 2]) * muB[, 1] * (wB * d * (1 - exp(-d * muB[, 2])) * exp(-d * muB[, 1]) - DB) / (psiB + exp(-muB[, 1] - muB[, 2]) * DB)) * IndB +
        ((yB[, 1] - muB[, 1]) + wB * d * muB[, 1] * (exp(-yB[, 2]) - exp(-d * muB[, 2])) * exp(-d * muB[, 1]) / DB) * !IndB
      
      dbeta11B <- (x[, 1] * exp(-muB[, 1] - muB[, 2]) * muB[, 1] * (wB * d * (1 - exp(-d * muB[, 2])) * exp(-d * muB[, 1]) - DB) / (psiB + exp(-muB[, 1] - muB[, 2]) * DB)) * IndB +
        ((yB[, 1] - muB[, 1]) * x[, 1] + wB * d * muB[, 1] * x[, 1] * (exp(-yB[, 2]) - exp(-d * muB[, 2])) * exp(-d * muB[, 1]) / DB) * !IndB
      
      
      ## BETAA2
      dbeta20B <- (exp(-muB[, 1] - muB[, 2]) * muB[, 2] * (wB * d * (1 - exp(-d * muB[, 1])) * exp(-d * muB[, 2]) - DB) / (psiB + exp(-muB[, 1] - muB[, 2]) * DB)) * IndB +
        ((yB[, 2] - muB[, 2]) + wB * d * muB[, 2] * (exp(-yB[, 1]) - exp(-d * muB[, 1])) * exp(-d * muB[, 2]) / DB) * !IndB
      
      dbeta21B <- (x[, 2] * exp(-muB[, 1] - muB[, 2]) * muB[, 2] * (wB * d * (1 - exp(-d * muB[, 1])) * exp(-d * muB[, 2]) - DB) / (psiB + exp(-muB[, 1] - muB[, 2]) * DB)) * IndB +
        ((yB[, 2] - muB[, 2]) * x[, 2] + wB * d * muB[, 2] * x[, 2] * (exp(-yB[, 1]) - exp(-d * muB[, 1])) * exp(-d * muB[, 2]) / DB) * !IndB
      
      
      ## OMEGA
      dwB = (exp(-muB[, 1] - muB[, 2]) * (1 - exp(-d * muB[, 1])) * (1 - exp(-d * muB[, 2])) / (psiB + exp(-muB[, 1] - muB[, 2]) * DB)) * IndB +
        ((exp(-yB[, 1]) - exp(-d * muB[, 1])) * (exp(-yB[, 2]) - exp(-d * muB[, 2])) / DB) * !IndB
      
      
      ## TAU
      a1B = sapply(1:nrow(yB), function(i){ 
        sv = 0
        if(yB[i, 1] > 0){
          for(v in 1:yB[i, 1]){
            sv = sv + (yB[i, 1] - v)
          }
        }
        return (sv)
      })
      
      a2B = sapply(1:nrow(yB), function(i){ 
        sv = 0
        if(yB[i, 2] > 0){
          for(v in 1:yB[i, 2]){
            sv = sv + (yB[i, 2] - v)
          }
        }
        return (sv)
      })
      
      dtau1B <- (exp(-muB[, 1] -muB[ ,2]) * DB * (muB[, 1]^2 / 2) / (exp(cbind(1, z) %*% gamB) + exp(-muB[, 1] - muB[, 2]) * DB)) * IndB -
        (wB * exp(-muB[, 1] - muB[, 2]) * exp(-d * muB[, 1]) * (1 - exp(-d * muB[, 2])) * (muB[, 1]^2 * d / 2)  / (exp(cbind(1, z) %*% gamB) + exp(-muB[, 1] - muB[, 2]) * DB)) * IndB +
        (a1B - muB[, 1] * yB[, 1] + muB[, 1]^2/ 2 - wB * (exp(-yB[, 2]) - cB[, 2]) * exp(-d * muB[, 1]) * (d * muB[, 1]^2 / 2) / (1 + DB)) * !IndB
      
      dtau2B <- (exp(-muB[, 1] -muB[ ,2]) * DB * (muB[, 2]^2 / 2) / (exp(cbind(1, z) %*% gamB) + exp(-muB[, 1] - muB[, 2]) * DB)) * IndB -
        (wB * exp(-muB[, 1] - muB[, 2]) * exp(-d * muB[, 2]) * (1 - exp(-d * muB[, 1])) * (muB[, 2]^2 * d / 2)  / (exp(cbind(1, z) %*% gamB) + exp(-muB[, 1] - muB[, 2]) * DB)) * IndB +
        (a2B - muB[, 2] * yB[, 2] + muB[, 2]^2/ 2 - wB * (exp(-yB[, 1]) - cB[, 1]) * exp(-d * muB[, 2]) * (d * muB[, 2]^2 / 2) / (1 + DB)) * !IndB
      
      
      
      gr_dfB <- data.frame(dtau1B, dtau2B, dgam0B, dgam1B, dbeta10B, dbeta11B, dbeta20B, dbeta21B, dwB)
      
      colSums(gr_dfB)
      
      IB <- matrix(c(colSums(gr_dfB * dtau1B),
                     colSums(gr_dfB * dtau2B),
                     colSums(gr_dfB * dgam0B),
                     colSums(gr_dfB * dgam1B),
                     colSums(gr_dfB * dbeta10B),
                     colSums(gr_dfB * dbeta11B),
                     colSums(gr_dfB * dbeta20B),
                     colSums(gr_dfB * dbeta21B),
                     colSums(gr_dfB * dwB)), 9, 9)
      
      scoreB <- colSums(gr_dfB)
      song <- try(ginv(IB), silent = TRUE)
      
      if (class(song)[[1]] == "try-error") {
        next
      } else {
        
      }
      
      score_b[b] <- t(scoreB) %*% ginv(IB) %*% scoreB
      
    }
    
    boot_cric[i, ] <- quantile(score_b, c(.99, .95, .9, .8), na.rm = TRUE)
    
    if(i %% 10 == 0){
      print(i)
    }
    
  }
  # result_list[[num]] = result
  # boot_list[[num]] = boot_cric
  return(list(result, boot_cric))
}



result_list <- parLapply(cl, X = 1:30, please)

stopCluster(cl)


## 결과
# parameters for type1 error
n_grid <- c(50,100, 200)
g2_grid <- c(0.2, 0.8, 1.2)
w_grid <- c(-1.5, -1, 0, 1, 1.5)
param_grid <- expand.grid(w_grid, g2_grid, n_grid)

# parameters for power
n_grid <- c(100, 200)
# g2_grid <- c(0.2, 0.8, 1.2)
tau_grid <- c(1e-10, 0.1, 0.2, 0.4, 0.8)
w_grid <- c(-1.5, 0, 1.5)

param_grid <- expand.grid(w_grid, tau_grid, n_grid)


qval = c(0.01, 0.05, 0.1, 0.2)
upp = qchisq(qval, 2, lower.tail = F)
for(qq in 1:4) print(mean(result >= upp[qq], na.rm = T))

result_list[[1]]$result


result_df = matrix(0, 30, 4)
# 스코어 결과
for(ii in 1:45){
  cat("n =", param_grid[ii, 3], "gam2 =", param_grid[ii, 2], "w =", param_grid[ii, 1], '\n')
  for(qq in 1:4){
    cat("true q = ", qval[qq], "estimated prob =", round(sum(result_list_alpha[[ii]][[1]]>= upp[qq], na.rm = T)/length(result_list_alpha[[ii]][[1]]), 4), '\n')
    # result_df[ii, qq] <- round(sum(result_list_alpha[[ii]][[1]]>= upp[qq], na.rm = T)/length(result_list_alpha[[ii]][[1]]), 4)
  }
}


# 붓스트랩 결과
result_boot_df = matrix(0, 30, 1)
for(ii in 16:30){
  cat("n =", param_grid[ii, 3], "gam2 =", param_grid[ii, 2], "w =", param_grid[ii, 1], '\n')
  # aa = rowSums(sapply(1:1000, function(j) result_list_alpha[[ii]][[2]][,3] < result_list_power[[ii]][[1]][j]), na.rm = T)/1000
  # print(aa)
  result_boot_df[ii,] <- sum(sapply(1:1000, function(j){ result_list_alpha[[45]][[2]][j,2] < result_list_power[[ii]][[1]][j] }), na.rm = T)/1000
  cat("\n")
}

save(result_list, file = "~/bootstrap.RData")

##########################
# plot(0, type = "n", xlim = range(tau_grid), ylim = range(c(result_df, result_boot_df)), ylab = "power", xlab = "tau")
plot(0, type = "n", xlim = range(tau_grid), ylim = range(0, 1), ylab = "power", xlab = "tau", main = expression(w == 1.5))
sel100 = (seq(1, 15, 3) + 2)  # 0: w=-1.5, 1: w=0, 2: w=1.5
sel200 = (seq(16, 30, 3) + 2)

abline(h = 0.05, col = "gray")

for(ii in 1:4){
  if(ii == 1) points(tau_grid, result_df[sel100,2], type = 'b', pch = 1, col = "red", lty = 1, lwd =2)
  if(ii == 2) points(tau_grid, result_df[sel200,2], type = 'b', pch = 19, col = "pink", lty = 2, lwd =2)
  if(ii == 3) points(tau_grid, result_boot_df[sel100,1], type = 'b', pch = 0, col = "blue", lty = 3, lwd =2)
  if(ii == 4) points(tau_grid, result_boot_df[sel200,1], type = 'b', pch = 15, col = "skyblue", lty = 4, lwd =2)
}

legend("bottomright", 
       c(expression(Chi^2~(2)~", "~n~"=100"),
         "Bootstrap, n=100", 
         expression(Chi^2~(2)~", "~n~"=200"),
         "Bootstrap, n=200"),
       lwd = 2, lty = 1:4, col = c("red", "pink", "blue", "skyblue"), pch = c(1, 19, 0, 15))

par(mfrow=c(1,3))

##########################
load(file = "~/bootstrap.RData")

# expression(Chi^2~(1)/2 + Chi^2~(2)/4~", "~n~"=100")










