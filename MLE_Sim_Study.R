library(dplyr)
library(future.apply)

log_likelihood <- function(r, p, y) {
  
  if (r <= 0 || p <= 0 || p >= 1) {
    return(-Inf)
  }
  loglik <- sum(dnbinom(y, size = r, prob = p, log = TRUE))
  if (is.nan(loglik)) {
    return(-Inf)
  } else {
    return(loglik)
  }

}

MLE_sim <- function(s, true_r, true_p, starting_vals) {
  
  sample <- rnbinom(s, true_r, true_p)
  
  mle_est <- optim(starting_vals, \(pars) -log_likelihood(pars[1], pars[2], sample))$par
  
  c(mle_est[1], mle_est[1] - true_r, (mle_est[1] - true_r)^2, mle_est[2], mle_est[2] - true_p, (mle_est[2] - true_p)^2)
  
}

MLE_uncertainty <- function(s, true_r, true_p, nBootstrapSamples, starting_vals) {
  
  dat <- rnbinom(s, true_r, true_p)
  
  bootstrap_MLE <- data.frame(t(replicate(nBootstrapSamples, {
    
    bootstrap_sample <- sample(dat, replace = TRUE)
    
    optim(starting_vals, \(pars) -log_likelihood(pars[1], pars[2], bootstrap_sample))$par # better with more data points
    
  })))
  
  names(bootstrap_MLE) <- c("r", "p")
  
  c(
    between(true_r, quantile(bootstrap_MLE$r, .025), quantile(bootstrap_MLE$r, .975)), # Coverage
    quantile(bootstrap_MLE$r, .975) - quantile(bootstrap_MLE$r, .025), # Width
    between(true_p, quantile(bootstrap_MLE$p, .025), quantile(bootstrap_MLE$p, .975)), # Coverage
    quantile(bootstrap_MLE$p, .975) - quantile(bootstrap_MLE$p, .025) # Width
  )
  
}

plan(multisession, workers = 8)

MLE_effectiveness <- cbind(
  t(Reduce("+", future_replicate(10000, sapply(c(10, 100, 1000), \(s) MLE_sim(s, 5, .1, c(1, .5))), simplify = FALSE)) / 10000),
  t(sapply(c(10, 100, 1000), \(s) rowMeans(future_replicate(100, MLE_uncertainty(s, 5, .1, 1000, c(1, .5))))))
)

MLE_effectiveness <- MLE_effectiveness[,c(1:3,7:8,4:6,9:10)]

MLE_effectiveness <- data.frame(rbind(MLE_effectiveness[,c(1:5)], MLE_effectiveness[,c(6:10)]))

names(MLE_effectiveness) <- c("est", "bias", "mse", "coverage", "width")

MLE_effectiveness <- MLE_effectiveness %>% mutate(n = rep(c("n = 10", "n = 100", "n = 1000"), 2))

MLE_effectiveness2 <- cbind(
  t(Reduce("+", future_replicate(10000, sapply(c(10, 100, 1000), \(s) MLE_sim(s, 10, .5, c(1, .5))), simplify = FALSE)) / 10000),
  t(sapply(c(10, 100, 1000), \(s) rowMeans(future_replicate(100, MLE_uncertainty(s, 10, .5, 1000, c(1, .5))))))
)

MLE_effectiveness2 <- MLE_effectiveness2[,c(1:3,7:8,4:6,9:10)]

MLE_effectiveness2 <- data.frame(rbind(MLE_effectiveness2[,c(1:5)], MLE_effectiveness2[,c(6:10)]))

names(MLE_effectiveness2) <- c("est", "bias", "mse", "coverage", "width")

MLE_effectiveness2 <- MLE_effectiveness2 %>% mutate(n = rep(c("n = 10", "n = 100", "n = 1000"), 2))

MLE_effectiveness3 <- cbind(
  t(Reduce("+", future_replicate(10000, sapply(c(10, 100, 1000), \(s) MLE_sim(s, 25, .9, c(1, .5))), simplify = FALSE)) / 10000),
  t(sapply(c(10, 100, 1000), \(s) rowMeans(future_replicate(100, MLE_uncertainty(s, 5, .1, 1000, c(25, .9))))))
)

MLE_effectiveness3 <- MLE_effectiveness3[,c(1:3,7:8,4:6,9:10)]

MLE_effectiveness3 <- data.frame(rbind(MLE_effectiveness3[,c(1:5)], MLE_effectiveness3[,c(6:10)]))

names(MLE_effectiveness3) <- c("est", "bias", "mse", "coverage", "width")

MLE_effectiveness3 <- MLE_effectiveness3 %>% mutate(n = rep(c("n = 10", "n = 100", "n = 1000"), 2))

save(
  list = c("MLE_effectiveness", "MLE_effectiveness2", "MLE_effectiveness3"),
  file = "Project_624_Simulation_Study_MLE.RData"
)