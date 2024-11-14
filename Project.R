library(rvest)

# Goals Messi scored (count data)

lionel_messi_url <- "https://fbref.com/en/players/d70ce98e/dom_lg/Lionel-Messi-Domestic-League-Stats"

lionel_messi_page_html <- lionel_messi_url %>%
  rvest::read_html()  

lionel_messi <- lionel_messi_url %>%
  read_html %>%
  html_nodes('table') %>%
  html_table() %>%
  .[[1]]  %>%
  setNames(make.unique(unlist(.[1,]))) %>% 
  slice(-1L) %>% 
  filter(Matches == "Matches" & Season != "" & Season != "2022-2023") %>% # 
  dplyr::select(1:18)

lionel_messi_gls_df <- lionel_messi %>% dplyr::select(Gls)

lionel_messi_gls <- as.numeric(lionel_messi_gls_df$Gls)

# 3 methods on real data

# MLE

# log_likelihood <- function(r, p, y) {
#   ifelse(
#     p < 0 | p > 1 | r < 0, 
#     -Inf, 
#     sum(lfactorial(y + r - 1) - (lfactorial(r - 1) + lfactorial(y)) + y*log(1 - p) + r*log(p))
#   )
# }

log_likelihood <- function(r, p, y) {
  ifelse(p < 0 | p > 1 | r < 0, -Inf, sum(dnbinom(y, r, p, log = TRUE)))
}

mle_est <- optim(c(10, .3), \(pars) -log_likelihood(pars[1], pars[2], lionel_messi_gls))$par

optim(c(10, .3), \(pars) -ifelse(pars[2] < 0 | pars[2] > 1 | pars[1] < 0,  -Inf, sum(dnbinom(lionel_messi_gls, pars[1], pars[2], log = TRUE))))$par

# Method of Moments

method_of_moments <- function(y) {
  
  r <- (mean(y)^2) / (var(y) - mean(y))
  
  p <- 1 - (mean(y) / var(y))
  
  return(c(r, p))
  
}

mom_est <- method_of_moments(lionel_messi_gls)

# ????

# Bayes estimate under squared error loss

y <- lionel_messi_gls

log_prior <- function(r, p) {
  dgamma(r, .01, .01, log = TRUE) + dbeta(p, 1, 1, log = TRUE)
}

log_likelihood <- function(r, p) {
  sum(dnbinom(y, r, p, log = TRUE))
}

log_posterior <- function(r, p) {
  log_likelihood(r, p) + log_prior(r, p)
}

nSamples <- 100000

draws <- data.frame(matrix(0, nrow = nSamples, ncol = 2))

# Initial state
draws[1,] <- mle_est

for (i in 2:nSamples) {
  
  r <- draws[i - 1, 1]
  p <- draws[i - 1, 2]
  
  # Propose rstar from a normal distribution centered at r
  proposal <- rnorm(1, r, 1)
  
  # Calculate the log Metropolis ratio
  log_m_ratio <- ifelse(proposal < 0, -Inf, log_posterior(proposal, p) - log_posterior(r, p))
  
  # Accept the proposal with the appropriate probability
  r <- ifelse(log_m_ratio > log(runif(1)), proposal, r)
  
  # Gibbs update of p
  # Full conditional is conjugate
  p <- rbeta(1, 1 + r*length(y), 1 + sum(y))
  
  draws[i,] <- c(r, p)
  
}

plot(draws$X1, type = 'l')
plot(draws$X2, type = 'l')

coda::effectiveSize(draws)

# Bayes estimate under SEL
colMeans(draws)

# Negative binomial

# 3 different settings of r (size) and p (prob)

# r <- 5, 10, 25
# p <- .1, .5, .9


# 3 different n 10, 100, 1000

dat <- rnbinom(100, 5, .1)

# bootstrap_dat <- sample(dat, replace = TRUE)

# optim(c(1, .1), \(pars) -log_likelihood(pars[1], pars[2], bootstrap_dat))$par

bootstrap_MLE <- data.frame(t(replicate(1000, {
  
  bootstrap_sample <- sample(dat, replace = TRUE)
  
  optim(c(1, .1), \(pars) -log_likelihood(pars[1], pars[2], bootstrap_sample))$par
  
})))

names(bootstrap_MLE) <- c("r", "p")

colMeans(bootstrap_MLE)
t(apply(bootstrap_MLE, 2, quantile, c(.025, .975)))

# Bias
mean(sapply(1:1000, \(i) bootstrap_MLE[i, 1] - 5))
mean(sapply(1:1000, \(i) bootstrap_MLE[i, 2] - .1))

colMeans(bootstrap_MLE)[1] - 5

# Just average bias or uncertainty of bias??

# MLE is point estimate??

# MSE
mean(sapply(1:1000, \(i) (bootstrap_MLE[i, 1] - 5)^2))
mean(sapply(1:1000, \(i) (bootstrap_MLE[i, 2] - .1)^2))

mean((bootstrap_MLE$r - 5)^2)

mean((colMeans(bootstrap_MLE)[1] - 5)^2)

# Coverage
as.numeric(between(5, quantile(bootstrap_MLE$r, .025), quantile(bootstrap_MLE$r, .975)))
between(.1, quantile(bootstrap_MLE$p, .025), quantile(bootstrap_MLE$p, .975))

c(
  mean(bootstrap_MLE$r), # Point estimate
  mean(bootstrap_MLE$r - 5), # Bias
  mean((bootstrap_MLE$r - 5)^2), # MSE
  between(5, quantile(bootstrap_MLE$r, .025), quantile(bootstrap_MLE$r, .975)), # Coverage
  mean(bootstrap_MLE$p), # Point estimate
  mean(bootstrap_MLE$p - .1), # Bias
  mean((bootstrap_MLE$p - .1)^2), # MSE
  between(.1, quantile(bootstrap_MLE$p, .025), quantile(bootstrap_MLE$p, .975)) # Coverage
)

# Run the bootstrap procedure multiple times?
# How many times ??

MLE_simulation_study <- data.frame(t(replicate(100, {
  
  bootstrap_MLE <- data.frame(t(replicate(1000, {
    
    bootstrap_sample <- sample(dat, replace = TRUE)
    
    optim(c(1, .1), \(pars) -log_likelihood(pars[1], pars[2], bootstrap_sample))$par
    
  })))
  
  names(bootstrap_MLE) <- c("r", "p")
  
  c(
    mean(bootstrap_MLE$r), # Point estimate
    mean(bootstrap_MLE$r - 5), # Bias
    mean((bootstrap_MLE$r - 5)^2), # MSE
    between(5, quantile(bootstrap_MLE$r, .025), quantile(bootstrap_MLE$r, .975)), # Coverage
    mean(bootstrap_MLE$p), # Point estimate
    mean(bootstrap_MLE$p - .1), # Bias
    mean((bootstrap_MLE$p - .1)^2), # MSE
    between(.1, quantile(bootstrap_MLE$p, .025), quantile(bootstrap_MLE$p, .975)) # Coverage
  )
  
})))

colMeans(MLE_simulation_study)

bootstrap_MoM <- replicate(1000, {
  
  bootstrap_sample <- sample(dat, replace = TRUE)
  
  method_of_moments(bootstrap_sample)
  
})

rowMeans(bootstrap_MoM)
cbind(t(apply(bootstrap_MoM, 1, quantile, c(.025, .975))))

MoM_simulation_study <- data.frame(t(replicate(100, {
  
  bootstrap_MoM <- data.frame(t(replicate(1000, {
    
    bootstrap_sample <- sample(dat, replace = TRUE)
    
    method_of_moments(bootstrap_sample)
    
  })))
  
  names(bootstrap_MoM) <- c("r", "p")
  
  c(
    mean(bootstrap_MoM$r), # Point estimate
    mean(bootstrap_MoM$r - 5), # Bias
    mean((bootstrap_MoM$r - 5)^2), # MSE
    between(5, quantile(bootstrap_MoM$r, .025), quantile(bootstrap_MoM$r, .975)), # Coverage
    mean(bootstrap_MoM$p), # Point estimate
    mean(bootstrap_MoM$p - .1), # Bias
    mean((bootstrap_MoM$p - .1)^2), # MSE
    between(.1, quantile(bootstrap_MoM$p, .025), quantile(bootstrap_MoM$p, .975)) # Coverage
  )
  
})))

colMeans(MoM_simulation_study)

nbinom_mcmc <- function(dat, iter = 1000) {
  
  log_prior <- function(r, p) {
    dgamma(r, .01, .01, log = TRUE) + dbeta(p, 1, 1, log = TRUE)
  }
  
  log_likelihood <- function(r, p) {
    ifelse(p < 0 | p > 1 | r < 0, -Inf, sum(dnbinom(dat, r, p, log = TRUE)))
  }
  
  log_posterior <- function(r, p) {
    log_likelihood(r, p) + log_prior(r, p)
  }
  
  draws <- data.frame(matrix(0, nrow = iter, ncol = 2))
  
  # Initial state
  draws[1,] <- c(1, .1)
  
  for (i in 2:iter) {
    
    r <- draws[i - 1, 1]
    p <- draws[i - 1, 2]
    
    # Propose rstar from a normal distribution centered at r
    proposal <- rnorm(1, r, 1)
    
    # Calculate the log Metropolis ratio
    log_m_ratio <- ifelse(proposal < 0, -Inf, log_posterior(proposal, p) - log_posterior(r, p))
    
    # Accept the proposal with the appropriate probability
    r <- ifelse(log_m_ratio > log(runif(1)), proposal, r)
    
    # Gibbs update of p
    # Full conditional is conjugate
    p <- rbeta(1, 1 + r*length(dat), 1 + sum(dat))
    
    draws[i,] <- c(r, p)
    
  }
  
  return(draws)
  
}

Bayes_simulation_study <- data.frame(t(replicate(100, {
  
  markov_chain <- nbinom_mcmc(dat)
  
  names(markov_chain) <- c("r", "p")
  
  c(
    mean(markov_chain$r), # Point estimate
    mean(markov_chain$r - 5), # Bias
    mean((markov_chain$r - 5)^2), # MSE
    between(5, quantile(markov_chain$r, .025), quantile(markov_chain$r, .975)), # Coverage
    mean(markov_chain$p), # Point estimate
    mean(markov_chain$p - .1), # Bias
    mean((markov_chain$p - .1)^2), # MSE
    between(.1, quantile(markov_chain$p, .025), quantile(markov_chain$p, .975)) # Coverage
  )
  
})))

colMeans(Bayes_simulation_study)

simulation_study_df <- data.frame(rbind(
  colMeans(MLE_simulation_study),
  colMeans(MoM_simulation_study),
  colMeans(Bayes_simulation_study)
))

simulation_study <- function(n, true_r, true_p, nSamples, nReplicates) {
  
  dat <- rnbinom(n, true_r, true_p)
  
  MLE_simulation_study <- data.frame(t(replicate(nReplicates, {
    
    bootstrap_MLE <- data.frame(t(replicate(nSamples, {
      
      bootstrap_sample <- sample(dat, replace = TRUE)
      
      optim(c(1, .1), \(pars) -log_likelihood(pars[1], pars[2], bootstrap_sample))$par
      
    })))
    
    names(bootstrap_MLE) <- c("r", "p")
    
    c(
      mean(bootstrap_MLE$r), # Point estimate
      mean(bootstrap_MLE$r - true_r), # Bias
      mean((bootstrap_MLE$r - true_r)^2), # MSE
      between(true_r, quantile(bootstrap_MLE$r, .025), quantile(bootstrap_MLE$r, .975)), # Coverage
      mean(bootstrap_MLE$p), # Point estimate
      mean(bootstrap_MLE$p - true_p), # Bias
      mean((bootstrap_MLE$p - true_p)^2), # MSE
      between(true_p, quantile(bootstrap_MLE$p, .025), quantile(bootstrap_MLE$p, .975)) # Coverage
    )
    
  })))
  
  MoM_simulation_study <- data.frame(t(replicate(nReplicates, {
    
    bootstrap_MoM <- data.frame(t(replicate(nSamples, {
      
      bootstrap_sample <- sample(dat, replace = TRUE)
      
      method_of_moments(bootstrap_sample)
      
    })))
    
    names(bootstrap_MoM) <- c("r", "p")
    
    c(
      mean(bootstrap_MoM$r), # Point estimate
      mean(bootstrap_MoM$r - true_r), # Bias
      mean((bootstrap_MoM$r - true_r)^2), # MSE
      between(true_r, quantile(bootstrap_MoM$r, .025), quantile(bootstrap_MoM$r, .975)), # Coverage
      mean(bootstrap_MoM$p), # Point estimate
      mean(bootstrap_MoM$p - true_p), # Bias
      mean((bootstrap_MoM$p - true_p)^2), # MSE
      between(true_p, quantile(bootstrap_MoM$p, .025), quantile(bootstrap_MoM$p, .975)) # Coverage
    )
    
  })))
  
  Bayes_simulation_study <- data.frame(t(replicate(nReplicates, {
    
    markov_chain <- nbinom_mcmc(dat, nSamples)
    
    names(markov_chain) <- c("r", "p")
    
    c(
      mean(markov_chain$r), # Point estimate
      mean(markov_chain$r - true_r), # Bias
      mean((markov_chain$r - true_r)^2), # MSE
      between(true_r, quantile(markov_chain$r, .025), quantile(markov_chain$r, .975)), # Coverage
      mean(markov_chain$p), # Point estimate
      mean(markov_chain$p - true_p), # Bias
      mean((markov_chain$p - true_p)^2), # MSE
      between(true_p, quantile(markov_chain$p, .025), quantile(markov_chain$p, .975)) # Coverage
    )
    
  })))
  
  simulation_study_df <- data.frame(rbind(
    colMeans(MLE_simulation_study),
    colMeans(MoM_simulation_study),
    colMeans(Bayes_simulation_study)
  ))
  
  names(simulation_study_df) <- c("r_est", "r_bias", "r_mse", "r_coverage", "p_est", "p_bias", "p_mse", "p_coverage")
  
  return(simulation_study_df)
  
}

setting1_sim_study <- simulation_study(10, 5, .1, 1000, 1000)
setting2_sim_study <- simulation_study(100, 5, .1, 1000, 1000)
setting3_sim_study <- simulation_study(1000, 5, .1, 1000, 1000)

setting4_sim_study <- simulation_study(10, 10, .5, 1000, 1000)
setting5_sim_study <- simulation_study(100, 10, .5, 1000, 1000)
setting6_sim_study <- simulation_study(1000, 10, .5, 1000, 1000)

setting7_sim_study <- simulation_study(10, 25, .9, 1000, 1000)
setting8_sim_study <- simulation_study(100, 25, .9, 1000, 1000)
setting9_sim_study <- simulation_study(1000, 25, .9, 1000, 1000)

markov_chain <- nbinom_mcmc(dat)

names(markov_chain) <- c("r", "p")

c(
  mean(markov_chain$r), # Point estimate
  mean(markov_chain$r - 5), # Bias
  mean((markov_chain$r - 5)^2), # MSE
  between(5, quantile(markov_chain$r, .025), quantile(markov_chain$r, .975)), # Coverage
  mean(markov_chain$p), # Point estimate
  mean(markov_chain$p - .1), # Bias
  mean((markov_chain$p - .1)^2), # MSE
  between(.1, quantile(markov_chain$p, .025), quantile(markov_chain$p, .975)) # Coverage
)

plot(thing$X1, type = 'l')
plot(thing$X2, type = 'l')

coda::effectiveSize(draws)

colMeans(thing)
t(apply(thing, 2, quantile, c(.025, .975)))

# MC error or credible interval??

# Same number of bootstrap samples & MCMC iterations ??

# Explore the data and apply the three estimators and uncertainty intervals using the chosen distribution to the data.

quantile(rnbinom(1000, mle_est[1], mle_est[2]), c(.025, .975))

quantile(rnbinom(1000, mom_est[1], mom_est[2]), c(.025, .975))



# missclassification rate => they are distinct

# all low => distinct clusters

# How are they different?

# 4 groups => nor more than 3 discriminant functions

# 1st one = most important way

# easier to predict for lower k

# missclassification rate with go down as k goes down... diminishing returns

# elbow method