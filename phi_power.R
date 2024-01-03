rm(list = ls())

library(parallel)
library(dplyr)
library(MASS)

numCores <- detectCores() - 1
cl <- makeCluster(numCores)

phi_power_fun <- function(num) {
  
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
  
  logL_H0 <- function(param) {
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
  
  logL_H1 <- function(param) {
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
    
    v1 <- ifelse(v1[ind] <= 0, log(1e-15), log(v1[ind]))
    v2 <- ifelse(v2[!ind] <= 0, log(1e-15), log(v2[!ind]))
    
    return(-(sum(v1) + sum(v2)))
  }
  
  ## Iteration
  iteration = 1000
  
  result <- matrix(0, nrow = iteration, ncol = 2)
  
  n_grid <- c(200, 500)
  tau_grid <- c(0.6, 1.2)
  phi_grid <- c(0, 0.025, 0.05, 0.1, 0.2, 0.4)
  w_grid <- c(-1, 0, 1)
  
  param_grid <- expand.grid(w_grid, phi_grid, tau_grid, n_grid)
  
  for(i in 1:iteration) {
    set.seed(i)
    sample_data <- RBZI(n = param_grid[num, 4], param_grid[num, 2], .2, .4, .4, .8, .8, param_grid[num, 3], param_grid[num, 1])
    
    one <- try(optim(par = rep(0, 7), fn = logL_H0, method = "L", lower = c(rep(-Inf, 4),  1e-10, 1e-10, -Inf)), silent = TRUE)
    two <- try(optim(par = rep(0, 8), fn = logL_H1, method = "L", lower = c(rep(-Inf, 5),  1e-10, 1e-10, -Inf)), silent = TRUE)
    
    if (class(one) == "try-error" | class(two) == "try-error") {
      next
    } 
    
    param <- one$par
    
    ## LR-test
    L0val <- one$value
    L1val <- two$value
    LR <- 2 * (L0val - L1val)
    
    beta_i <- param[c(1, 3)]
    beta_c <- param[c(2, 4)]
    tau <- param[(c(5, 6))]
    w <- param[7]
    
    d <- 1 - exp(-1)
    
    y <- sample_data$y; x <- cbind(sample_data$x1[, 2], sample_data$x2[, 2])
    mu <- cbind(exp(sample_data$x1 %*% c(beta_i[1], beta_c[1])), exp(sample_data$x2 %*% c(beta_i[2], beta_c[2])))
    c <- cbind((1 + d * mu[, 1] * tau[1])^{-1/tau[1]}, (1 + d * mu[, 2] * tau[2])^{-1/tau[2]})
    
    Ind <- (y[, 1] == 0 & y[, 2] == 0)
    
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
    
    g_b10 <- (y[, 1] - mu[, 1]) / (1 + tau[1] * mu[, 1]) + w * d * mu[, 1] * (exp(-y[, 2]) - c[, 2]) * (1 + d * tau[1] * mu[, 1])^{-1/tau[1] - 1} / (1 + w * (exp(-y[, 1]) - c[, 1]) * (exp(-y[, 2]) - c[, 2]))
    g_b11 <- (y[, 1] - mu[, 1]) * x[, 1] / (1 + tau[1] * mu[, 1]) + w * d * mu[, 1] * x[, 1] * (exp(-y[, 2]) - c[, 2]) * (1 + d * tau[1] * mu[, 1])^{-1/tau[1] - 1} / (1 + w * (exp(-y[, 1]) - c[, 1]) * (exp(-y[, 2]) - c[, 2]))
    g_b20 <- (y[, 2] - mu[, 2]) / (1 + tau[2] * mu[, 2]) + w * d * mu[, 2] * (exp(-y[, 1]) - c[, 1]) * (1 + d * tau[2] * mu[, 2])^{-1/tau[2] - 1} / (1 + w * (exp(-y[, 1]) - c[, 1]) * (exp(-y[, 2]) - c[, 2]))
    g_b21 <- (y[, 2] - mu[, 2]) * x[, 2] / (1 + tau[2] * mu[, 2]) + w * d * mu[, 2] * x[, 2] * (exp(-y[, 1]) - c[, 1]) * (1 + d * tau[2] * mu[, 2])^{-1/tau[2] - 1} / (1 + w * (exp(-y[, 1]) - c[, 1]) * (exp(-y[, 2]) - c[, 2]))
    
    g_w <- (exp(-y[, 1]) - c[, 1]) *  (exp(-y[, 2]) - c[, 2]) / (1 + w * (exp(-y[, 1]) - c[, 1]) *  (exp(-y[, 2]) - c[, 2]))
    g_psi <- (1 - (1 + tau[2] * mu[, 2])^{-1/tau[2]} * (1 + tau[1] * mu[, 1])^{-1/tau[1]} * (1 + w * (1 - c[, 1]) * (1 - c[, 2]))) / ((1 + tau[2] * mu[, 2])^{-1/tau[2]} * (1 + tau[1] * mu[, 1])^{-1/tau[1]} * (1 + w * (1 - c[, 1]) * (1 - c[, 2]))) * Ind + (-1) * !Ind
    
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
    
    sing <- ginv(I)
    
    if (class(sing)[[1]] == "try-error") {
      next
    } else {
      
    }
    score <- score %*% ginv(I) %*% t(score)
    
    ## Result-Store
    
    result[i, 1] <- score
    result[i, 2] <- LR
  }
  return(result)
}

phi_power <- parLapply(cl, X = 1:72, phi_power_fun)
stopCluster(cl)


#####################################################

load("C:/Users/jieun/Dropbox/BZINB/zero-zero/graph/tau_power_score.RData")

n_grid <- c(200, 500)
tau_grid <- c(0.6, 1.2)
phi_grid <- c(0, 0.025, 0.05, 0.1, 0.2, 0.4)
w_grid <- c(-1, 0, 1)

param_grid <- expand.grid(w_grid, phi_grid, tau_grid, n_grid)

qval = c(0.05)
upp = qchisq(qval, 1, lower.tail = F)
low = qchisq(qval, 1)


# Open PNG file for writing
tiff('phi_power.tif', units = "px", res = 300, width = 2600, height = 1500)

# Set the outer margin
par(oma = c(0, 0, 0, 0))

# Set the inner margin
par(mar = c(5, 4, 4, 2) + 0.1)


par(mfrow = c(2, 3))

kk = 0
pn = 0
for(tt in tau_grid){
  pn = pn + 1
  for(ww in w_grid){
    kk = kk + 1
    
    out <- data.frame(n = rep(0, 12), phi = rep(0, 12), LR = rep(0, 12), score = rep(0, 12))
    
    # Score_mixed 결과
    k = 0
    for(ii in 1:72){
      if(param_grid[ii, 3] == tt & param_grid[ii, 1] == ww){  ## tau1 조정
        cat("n =", param_grid[ii, 4], "tau1 =", param_grid[ii, 3], "phi =", param_grid[ii, 2], "w =", param_grid[ii, 1], '\n')
        k = k + 1
        print(ii)
        out$n[k] <- param_grid[ii, 4]
        out$phi[k] <- param_grid[ii, 2]
        out$LR[k] <- round(mean(phi_power[[ii]][, 2] >= upp , na.rm = T), 3)
        out$score[k] <- round(mean(phi_power[[ii]][, 1] >= upp , na.rm = T), 3)
      }
    }
    
    # if(ww == -1.5){
    #   out_full <- out  
    # } else{
    #   out_full <- rbind(out_full, out)
    # }
    print(out$score[1])
    
    # 그림 그리기
    # n=200, LR
    plot(out$phi[out$n == 500], out$LR[out$n == 500], lty = 1, type = 'l', col = 'darkgray',
         xlab = expression(phi), ylab = "estimated power", main = bquote(tau[2] == .(tau_grid[pn])  ~ ", " ~ w == .(ifelse(ww == 0, 0, ww))), ylim = c(0, 1))
    points(out$phi[out$n == 500], out$LR[out$n == 500], lty = 1, lwd = 1, pch = 15, col = 'darkgray')
    
    # n=100, LR
    points(out$phi[out$n == 200], out$LR[out$n == 200], lty = 2, type = 'l', col = 'darkgray',
           xlab = expression(phi), ylab = "estimated power")
    points(out$phi[out$n == 200], out$LR[out$n == 200], lty = 2, lwd = 1, pch = 15, col = 'darkgray')
    
    # n=200, score
    points(out$phi[out$n == 500], out$score[out$n == 500], lty = 1, type = 'l', col = 'black',
           xlab = expression(phi), ylab = "estimated power")
    points(out$phi[out$n == 500], out$score[out$n == 500], lwd = 1, pch = 16, col = 'black')
    
    # n=100, score
    points(out$phi[out$n == 200], out$score[out$n == 200], lty = 2, type = 'l', col = 'black',
           xlab = expression(phi), ylab = "estimated power")
    points(out$phi[out$n == 200], out$score[out$n == 200], lty = 2, lwd = 1, pch = 16, col = 'black')
    
    abline(h = 0.05)
    
    if(kk == 6)
      legend("bottomright", legend = c('n=500, score', 'n=200, score', 'n=500, LR', 'n=200, LR'),
             col = c("black", "black", "darkgray","darkgray"), lty = c(1, 3, 1, 3), lwd = rep(1, 4),
             pch = c(16, 16, 15, 15))
    
  }
}
dev.off()
