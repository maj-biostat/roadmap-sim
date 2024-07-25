source("./R/init.R")



g_silo = c("e", "l" , "c") 
g_pr_silo <- c(0.3, 0.5, 0.2)
names(g_pr_silo) <- g_silo

g_pr_jnt <- matrix(
  c(0.4, 0.6, 0.7, 0.3, 0.5, 0.5), 3, 2, byrow = T
)
colnames(g_pr_jnt) <- c("knee", "hip")
rownames(g_pr_jnt) <- g_silo



g_pr_e_surg <- c(0.9, 0.1, 0)

# preference for one/two-stage under revision
g_pr_l_surg <- c(0.3, 0.7)

# g_pr_l_surg <- rbind(
#   dair = c(0.2, 0.25, 0.56),
#   rev = c(0, 0.3, 0.7)
# )
g_pr_c_surg <- c(0.2, 0.2, 0.6)



N = 1e5
mu = 0
# all effects are sum to zero
# silo (e, l, c)
b_silo = c(0.4, 0.1, -0.5)
b_silo = c(qlogis(0.7), qlogis(0.6), qlogis(0.6))
# jnt (knee, hip)
b_jnt = c(-0.4, 0.4)
# domains:
# surgery 

b_pref = c(0, -0.1, 0.1)
# non-rand: dair, one, two
# rand: dair, rev (one), rev (two)
b_d1 = c(-0.4, -0.4, -0.4, 0, 1, 0.2)
# ab backbone duration 
# non-rand (received dair), non-rand (received one), non-rand (received two),
# 6wk, 12wk
b_d2 = c(-0.2, -0.4, 0.6)
# ext proph duration (non-rand, 12wk, 7day)
b_d3 = c(-0.4, 0.1, 0.3)
# ab choice (non-rand, no-rif, rif)
b_d4 = c(0.3, -0.1, -0.2)


get_sim_data <- function(
    N = 1000,
    
    mu = 0,
    
    # all effects are sum to zero
    
    # silo (e, l, c)
    b_silo = c(0.4, 0.1, -0.5),
    
    # b_silo is ref log-odds for each silo
    # b_silo = c(qlogis(0.7), qlogis(0.6), qlogis(0.6)),
    
    
    # jnt (knee, hip)
    b_jnt = c(-0.4, 0.4),
    
    # domains:
    # surgery 
    # preference effects under randomised revision
    # non-rand, pref one, pref two
    b_pref = c(0, -0.1, 0.1),
    # non-rand: dair, one, two
    # rand: dair, rev (one), rev (two)
    b_d1 = c(-0.4, -0.4, -0.4, 0, 1, 0.2),
    # ab backbone duration 
    # non-rand (received dair), non-rand (received one), non-rand (received two),
    # 6wk, 12wk
    b_d2 = c(-0.2, -0.4, 0.6),
    # ext proph duration (non-rand, 12wk, 7day)
    b_d3 = c(-0.4, 0.1, 0.3),
    # ab choice (non-rand, no-rif, rif)
    b_d4 = c(0.3, -0.1, -0.2)
    ){
  
  # stopifnot(all.equal(0, sum(b_silo), tol = sqrt(.Machine$double.eps)))
  stopifnot(all.equal(0, sum(b_jnt), tol = sqrt(.Machine$double.eps)))
  stopifnot(all.equal(0, sum(b_d1), tol = sqrt(.Machine$double.eps)))
  stopifnot(all.equal(0, sum(b_d2), tol = sqrt(.Machine$double.eps)))
  stopifnot(all.equal(0, sum(b_d3), tol = sqrt(.Machine$double.eps)))
  stopifnot(all.equal(0, sum(b_d4), tol = sqrt(.Machine$double.eps)))
  
  
  K_s <- length(b_silo)
  K_j <- length(b_jnt)
  K_p <- length(b_pref)
  K_d1 <- length(b_d1)
  K_d2 <- length(b_d2)
  K_d3 <- length(b_d3)
  K_d4 <- length(b_d4)
  

  d <- data.table()
  d[, pt := 1:N]
  
  # Stratification
  d[, silo := sample(1:K_s, size = N, replace = T, prob = g_pr_silo)]
  d[silo == 1, jnt := sample(1:K_j, size = .N, replace = T, prob = g_pr_jnt[1, ])]
  d[silo == 2, jnt := sample(1:K_j, size = .N, replace = T, prob = g_pr_jnt[2, ])]
  d[silo == 3, jnt := sample(1:K_j, size = .N, replace = T, prob = g_pr_jnt[3, ])]
  
  # Domains
  
  # Surgery
  # early - non-randomised
  # preference isn't applicable as there are no randomised options
  d[silo == 1, pref := 1]
  d[silo == 1, d1 := sample(1:3, size = .N, replace = T, prob = g_pr_e_surg)]
  
  # late (two-stage rand) - assumes full compliance
  d[silo == 2, pref := sample(2:3, size = .N, replace = T, prob = g_pr_l_surg)]
  d[silo == 2, d1 := rbinom(.N, 1, 0.5)]
  # for those that received revision (indexes 4, 5 and 6 into the b_d1 pars
  # correspond to dair, rev( one-stage), rev(two-stage))
  d[silo == 2 & d1 == 1, d1 := pref + 3]
  # for those that received dair
  d[silo == 2 & d1 == 0, d1 := 4]
  
  # late (two-stage rand) - assumes full compliance
  d[silo == 3, pref := 1]
  d[silo == 3, d1 := sample(1:3, size = .N, replace = T, prob = g_pr_c_surg)]
  
  # Outcome
  # d[, eta := mu + b_silo[silo] + b_jnt[jnt] + b_pref[pref] + b_d1[d1]]
  d[, eta := mu + b_silo[silo] + b_jnt[jnt] + b_d1[d1] ]
  d[, y := rbinom(N, 1, plogis(eta))]
 
  
  # Convert design to lower dim
  
  X_s <- diag(K_s)
  X_j <- diag(K_j)
  X_p <- diag(K_p)
  X_d1 <- diag(K_d1)
  X_d2 <- diag(K_d2)
  X_d3 <- diag(K_d3)
  X_d4 <- diag(K_d4)
  
  # correlation matrix to enforce sum to zero
  S_s <- X_s - (1 / K_s )
  S_j <- X_j - (1 / K_j )
  S_p <- X_p - (1 / K_p )
  S_d1 <- X_1 - (1 / K_d1 )
  S_d2 <- X_2 - (1 / K_d2 )
  S_d3 <- X_3 - (1 / K_d3 )
  S_d4 <- X_4 - (1 / K_d4 )
  
  # decomposition eigen vectors
  Q_s <- eigen(S_s)$vector[, 1:(K_s - 1)]
  Q_j <- eigen(S_j)$vector[, 1:(K_j - 1)]
  Q_p <- eigen(S_p)$vector[, 1:(K_p - 1)]
  Q_d1 <- eigen(S_d1)$vector[, 1:(K_d1 - 1)]
  Q_d2 <- eigen(S_d2)$vector[, 1:(K_d2 - 1)]
  Q_d3 <- eigen(S_d3)$vector[, 1:(K_d3 - 1)]
  Q_d4 <- eigen(S_d4)$vector[, 1:(K_d4 - 1)]
  
  # transformed pars
  b_s_s <- t(Q_s) %*% b_silo
  b_j_s <- t(Q_j) %*% b_jnt
  b_p_s <- t(Q_p) %*% b_pref
  b_d1_s <- t(Q_d1) %*% b_d1
  b_d2_s <- t(Q_d2) %*% b_d2
  b_d3_s <- t(Q_d3) %*% b_d3
  b_d4_s <- t(Q_d4) %*% b_d4
  
  # full rank design components
  X_s_s <- X_s %*% Q_s
  X_j_s <- X_j %*% Q_j
  X_p_s <- X_p %*% Q_p
  X_d1_s <- X_d1 %*% Q_d1
  X_d2_s <- X_d2 %*% Q_d2
  X_d3_s <- X_d3 %*% Q_d3
  X_d4_s <- X_d4 %*% Q_d4
  
  
  d_grid <- CJ(silo = 1:K_s, jnt = 1:K_j, p = 1:K_p, d1 = 1:K_d1 )
  # irrelevant combinations
  d_grid <- d_grid[! (silo == 1 & d1 %in% 4:6)]
  d_grid <- d_grid[! (silo == 1 & p %in% 2:3)]
  d_grid <- d_grid[! (silo == 3 & d1 %in% 4:6)]
  d_grid <- d_grid[! (silo == 3 & p %in% 2:3)]
  
  X_full_s <- matrix(NA, nrow(d_grid), 10)
  for(i in 1:nrow(d_grid)){
    X_full_s[i, ] <- c(
      X_s_s[d_grid[i, silo], ],
      X_j_s[d_grid[i, jnt], ],
      X_p_s[d_grid[i, p], ],
      X_d1_s[d_grid[i, d1], ]
      )
  }
  
  
  list(
    d = d, 
    
    mu = mu, 
    b_silo = b_silo, b_jnt = b_jnt, b_pref = b_pref, 
    b_d1 = b_d1, b_d2 = b_d2, b_d3 = b_d3, b_d4 = b_d4, 
    
    K_s = K_s, K_j = K_j, K_p  = K_p,
    K_d1 = K_d1, K_d2 = K_d2 , K_d3 = K_d3, K_d4 = K_d4,
    
    X_s_s  = X_s_s, X_j_s  = X_j_s, X_p_s  = X_p_s, 
    X_d1_s = X_d1_s, X_d2_s = X_d2_s,
    X_d3_s = X_d3_s, X_d4_s = X_d4_s, 
    X_full_s = X_full_s
  )
  
}


main <- function(){
  
  m1 <- cmdstanr::cmdstan_model("stan/model09.stan")
  ll <- get_sim_data(N = 1e3)
  d_mod <- ll$d[, .(y = sum(y), n = .N), keyby = .(silo, jnt, d1)]

  ld <- list(
    N = nrow(d_mod), y = d_mod$y, n = d_mod$n, 
    silo = d_mod$silo, 
    jnt = d_mod$jnt,
    # pref = d_mod$pref,
    d1 = d_mod$d1,
    
    nrXs = nrow(ll$X_s_s), ncXs = ncol(ll$X_s_s),
    Xsdes = ll$X_s_s, ss = rep(1, ncol(ll$X_s_s)),
    
    nrXj = nrow(ll$X_j_s), ncXj = ncol(ll$X_j_s), 
    Xjdes = ll$X_j_s, sj = rep(1, ncol(ll$X_j_s)),
    
    nrXp = nrow(ll$X_p_s), ncXp = ncol(ll$X_p_s), 
    Xpdes = ll$X_p_s, sp = rep(1, ncol(ll$X_p_s)),
    
    nrXd1 = nrow(ll$X_d1_s), ncXd1 = ncol(ll$X_d1_s), 
    Xd1des = ll$X_d1_s, sd1 = rep(1, ncol(ll$X_d1_s)),
    
    prior_only = 0
  )
  
  f1 <- m1$sample(
    ld, iter_warmup = 1000, iter_sampling = 1000,
    parallel_chains = 2, chains = 2, refresh = 0, show_exceptions = F,
    max_treedepth = 13)
  
  
  #
  m_post_s <- f1$draws(variables = c("mu"), format = "matrix")
  rbind(
    par = ll$mu,
    post_mu = colMeans(m_post_s)
  )
  
  #
  m_post_s <- f1$draws(variables = c("bs"), format = "matrix")
  post_mu <- m_post_s %*% t(cbind(ll$X_s_s))
  rbind(
    par = ll$b_silo,
    post_mu = colMeans(post_mu)
  )

  # #
  m_post_s <- f1$draws(variables = c("bj"), format = "matrix")
  post_mu <- m_post_s %*% t(cbind(ll$X_j_s))
  rbind(
    par = ll$b_jnt,
    post_mu = colMeans(post_mu)
  )
  # 
  # #
  # m_post_s <- f1$draws(variables = c("bp"), format = "matrix")
  # post_mu <- m_post_s %*% t(cbind(ll$X_p_s))
  # rbind(
  #   par = ll$b_pref,
  #   post_mu = colMeans(post_mu)  
  # )
  
  #
  # m_post_s <- f1$draws(variables = c("bd1"), format = "matrix")
  # post_mu <- m_post_s %*% t(cbind(ll$X_d1_s))
  # rbind(
  #   par = ll$b_d1,
  #   post_mu = colMeans(post_mu)  
  # )
  
}
