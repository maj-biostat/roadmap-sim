source("./R/init.R")



g_silo = c("e", "c", "l") 
g_pr_silo <- c(0.3, 0.2, 0.5)
names(g_pr_silo) <- g_silo

g_pr_jnt <- matrix(
  c(0.4, 0.6, 0.5, 0.5, 0.7, 0.3), 3, 2, byrow = T
)
colnames(g_pr_jnt) <- c("knee", "hip")
rownames(g_pr_jnt) <- g_silo


g_pr_e_surg <- c(0.85, 0.15)
g_pr_e_pref <- rbind(
  dair = c(0.85, 0.1, 0.05),
  rev = c(0, 2/3, 1/3)
)

g_pr_l_surg <- c(0.2, 0.8)
g_pr_l_pref <- rbind(
  dair = c(0.2, 0.24, 0.56),
  rev = c(0, 0.3, 0.7)
)

g_pr_c_surg <- c(0.2, 0.8)
g_pr_c_pref <- rbind(
  dair = c(0.2, 0.2, 0.6),
  rev = c(0, 0.25, 0.75)
)


get_sim_01_data <- function(
    N = 1000,
    
    mu = 0,
    b_silo = c(0.5, -0.3, -0.2),
    b_jnt = c(-0.4, 0.4),
    b_pref = c(0.2, -0.2),
    b_d1 = c(0.2, -0.1, -0.1), 
    b_g1 = rbind(
      c(-0.80,  0.00,  0.80),
      c( 0.40, -0.20, -0.20),
      c( 0.40,  0.20, -0.60)
    ),
    b_d2 = c(-0.2, -0.4, 0.6),
    b_d3 = c(-0.4, 0.1, 0.3),
    b_d4 = c(0.3, -0.1, -0.2)
    ){
  
  stopifnot(all.equal(0, sum(b_silo), tol = sqrt(.Machine$double.eps)))
  stopifnot(all.equal(0, sum(b_jnt), tol = sqrt(.Machine$double.eps)))
  stopifnot(all.equal(0, sum(b_g1), tol = sqrt(.Machine$double.eps)))
  stopifnot(all.equal(0, sum(b_d1), tol = sqrt(.Machine$double.eps)))

  stopifnot(all.equal(colSums(b_g1), rep(0, ncol(b_g1)),
                      tolerance = sqrt(.Machine$double.eps)))
  stopifnot(all.equal(rowSums(b_g1), rep(0, nrow(b_g1)),
                      tolerance = sqrt(.Machine$double.eps)))
  
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
  
  K_g1 <- K_s * K_d1
  
  d <- data.table()
  
  # Stratification
  d[, silo := sample(1:K_s, size = N, replace = T, prob = g_pr_silo)]
  d[silo == 1, jnt := sample(1:K_j, size = .N, replace = T, prob = g_pr_jnt[1, ])]
  d[silo == 2, jnt := sample(1:K_j, size = .N, replace = T, prob = g_pr_jnt[2, ])]
  # silo 3 is the late acute
  d[silo == 3, jnt := sample(1:K_j, size = .N, replace = T, prob = g_pr_jnt[3, ])]
  
  # Domains
  
  # Surgery - assignment to surgery dictates pretty much everything else.
  
  # early - non-randomised
  # preference is collinear with treatment selection
  d[silo == 1, pref_rev := sample(1:2, size = .N, replace = T, prob = g_pr_e_pref[2, 2:3])]
  d[silo == 1, d1 := sample(c(-98, -99), size = .N, replace = T, prob = g_pr_e_surg)]
  d[silo == 1 & d1 == -98, d1 := 1]
  # under revision, d1 has high probability of being the preference
  d[silo == 1 & d1 == -99 & pref_rev == 1, d1 := sample(2:3, .N, replace = T, prob = c(0.95, 0.05))]
  d[silo == 1 & d1 == -99 & pref_rev == 2, d1 := sample(2:3, .N, replace = T, prob = c(0.05, 0.95))]

  d[silo == 2, pref_rev := sample(1:2, size = .N, replace = T, prob = g_pr_c_pref[2, 2:3])]
  d[silo == 2, d1 := sample(c(-98, -99), size = .N, replace = T, prob = g_pr_c_surg)]
  d[silo == 2 & d1 == -98, d1 := 1]
  # under revision, d1 has high probability of being the preference
  d[silo == 2 & d1 == -99 & pref_rev == 1, d1 := sample(2:3, .N, replace = T, prob = c(0.95, 0.05))]
  d[silo == 2 & d1 == -99 & pref_rev == 2, d1 := sample(2:3, .N, replace = T, prob = c(0.05, 0.95))]
  
  # late acute is a randomised assignment for surgical domain
  d[silo == 3, pref_rev := sample(1:2, size = .N, replace = T, prob = g_pr_l_pref[2, 2:3])]
  # 1:1 rand assignment to dair (3) or rev (4) 
  d[silo == 3, d1 := sample(c(-98, -99), size = .N, replace = T)]
  # selection surgical type dair (1), one(2), two(3)
  # dair is fixed
  d[silo == 3 & d1 == -98, d1 := 1]
  # under revision, d1 has high probability of being the preference
  d[silo == 3 & d1 == -99 & pref_rev == 1, d1 := sample(2:3, .N, replace = T, prob = c(0.95, 0.05))]
  d[silo == 3 & d1 == -99 & pref_rev == 2, d1 := sample(2:3, .N, replace = T, prob = c(0.05, 0.95))]
  
  
  # Antibiotic backbone duration - only applicable to units recv one-stage
  
  d[d1 == 1, d2 := 1]
  d[d1 == 3, d2 := 1]
  d[d1 == 2, d2 := sample(2:3, size = .N, replace = T)]
  
  # d[, ix := rbinom(.N, 1, 0.6)]
  # d[ix == 0, d2 := 1]
  # d[ix == 1, d2 := sample(2:3, size = .N, replace = T)]
  
  # Extended prophylaxis - only applicable to units recv two-stage
  
  d[d1 == 1, d3 := 1]
  d[d1 == 2, d3 := 1]
  d[d1 == 3, d3 := sample(2:3, size = .N, replace = T)]

  # Outcome
  d[, mu := mu]
  d[, b_silo := b_silo[silo]]
  d[, b_jnt := b_jnt[jnt]]
  d[, b_pref := b_pref[pref_rev]]
  d[, b_d1 := b_d1[d1]]
  d[, b_g1 := b_g1[cbind(d$silo, d$d1)]]
  d[, b_d2 := b_d2[d2]]
  d[, eta := mu + b_silo + b_jnt + b_pref + b_d1 + b_d2     ]

  d[, y := rbinom(.N, 1, plogis(eta))]

  # Convert design to lower dim
  
  X_s <- diag(K_s)
  X_j <- diag(K_j)
  X_p <- diag(K_p)
  X_d1 <- diag(K_d1)
  X_d2 <- diag(K_d2)
  X_d3 <- diag(K_d3)
  X_d4 <- diag(K_d4)
  
  X_g1 <- diag(K_g1)
  
  # correlation matrix to enforce sum to zero
  S_s <- X_s - (1 / K_s )
  S_j <- X_j - (1 / K_j )
  S_p <- X_p - (1 / K_p )
  S_d1 <- X_d1 - (1 / K_d1 )
  S_d2 <- X_d2 - (1 / K_d2 )
  S_d3 <- X_d3 - (1 / K_d3 )
  S_d4 <- X_d4 - (1 / K_d4 )
  S_g1 <- kronecker(S_s, S_d1)
  
  # decomposition eigen vectors
  Q_s <- eigen(S_s)$vector[, 1:(K_s - 1)]
  Q_j <- eigen(S_j)$vector[, 1:(K_j - 1)]
  Q_p <- eigen(S_p)$vector[, 1:(K_p - 1)]
  Q_d1 <- eigen(S_d1)$vector[, 1:(K_d1 - 1)]
  Q_d2 <- eigen(S_d2)$vector[, 1:(K_d2 - 1)]
  Q_d3 <- eigen(S_d3)$vector[, 1:(K_d3 - 1)]
  Q_d4 <- eigen(S_d4)$vector[, 1:(K_d4 - 1)]
  Q_g1 <- kronecker(Q_s, Q_d1)
  
  # transformed pars
  b_s_s <- t(Q_s) %*% b_silo
  b_j_s <- t(Q_j) %*% b_jnt
  b_p_s <- t(Q_p) %*% b_pref
  b_d1_s <- t(Q_d1) %*% b_d1
  b_d2_s <- t(Q_d2) %*% b_d2
  b_d3_s <- t(Q_d3) %*% b_d3
  b_d4_s <- t(Q_d4) %*% b_d4
  
  b_g1_s <- t(Q_g1) %*% as.numeric(b_g1)
  
  # full rank design components
  X_s_s <- X_s %*% Q_s
  X_j_s <- X_j %*% Q_j
  X_p_s <- X_p %*% Q_p
  X_d1_s <- X_d1 %*% Q_d1
  X_d2_s <- X_d2 %*% Q_d2
  X_d3_s <- X_d3 %*% Q_d3
  X_d4_s <- X_d4 %*% Q_d4
  X_g1_s <- X_g1 %*% Q_g1
  
  # round(as.numeric(X_g1_s %*% b_g1_s), 3)
  # as.numeric(b_g1)

  list(
    d = d, 
    
    mu = mu, 
    b_silo = b_silo, b_jnt = b_jnt, b_pref = b_pref, 
    b_g1 = b_g1, 
    b_d1 = b_d1, b_d2 = b_d2, b_d3 = b_d3, b_d4 = b_d4, 
    
    K_s = K_s, K_j = K_j, K_p  = K_p,
    K_g1 = K_g1,
    K_d1 = K_d1, K_d2 = K_d2 , K_d3 = K_d3, K_d4 = K_d4,
    
    X_s_s  = X_s_s, X_j_s  = X_j_s, X_p_s  = X_p_s, 
    X_g1_s = X_g1_s,
    X_d1_s = X_d1_s, X_d2_s = X_d2_s,
    X_d3_s = X_d3_s, X_d4_s = X_d4_s
    # , 
    # X_full_s = X_full_s
  )
  
}






# basically the same as get_sim_01 but removes the sum to zero 
# manipulations and associated checks
get_sim_02_data <- function(
    N = 1000,
    
    mu = 0,
    b_silo = c(0, -0.3, -0.2),
    b_jnt = c(0, 0.4),
    b_pref = c(0, -0.2),
    b_d1 = c(0, -0.1, -0.1),
    b_d2 = c(0, -0.4, 0.6),
    b_d3 = c(0, 0.1, 0.3),
    b_d4 = c(0, -0.1, -0.2)
){

  K_s <- length(b_silo)
  K_j <- length(b_jnt)
  K_p <- length(b_pref)
  K_d1 <- length(b_d1)
  K_d2 <- length(b_d2)
  K_d3 <- length(b_d3)
  K_d4 <- length(b_d4)
  
  d <- data.table()
  
  # Stratification
  d[, silo := sample(1:K_s, size = N, replace = T, prob = g_pr_silo)]
  d[silo == 1, jnt := sample(1:K_j, size = .N, replace = T, prob = g_pr_jnt[1, ])]
  d[silo == 2, jnt := sample(1:K_j, size = .N, replace = T, prob = g_pr_jnt[2, ])]
  # silo 3 is the late acute
  d[silo == 3, jnt := sample(1:K_j, size = .N, replace = T, prob = g_pr_jnt[3, ])]
  
  # Domains
  
  # Surgery - assignment to surgery dictates pretty much everything else.
  
  # early - non-randomised
  # preference is collinear with treatment selection
  d[silo == 1, pref_rev := sample(1:2, size = .N, replace = T, prob = g_pr_e_pref[2, 2:3])]
  d[silo == 1, d1 := sample(c(-98, -99), size = .N, replace = T, prob = g_pr_e_surg)]
  d[silo == 1 & d1 == -98, d1 := 1]
  # under revision, d1 has high probability of being the preference
  d[silo == 1 & d1 == -99 & pref_rev == 1, d1 := sample(2:3, .N, replace = T, prob = c(0.95, 0.05))]
  d[silo == 1 & d1 == -99 & pref_rev == 2, d1 := sample(2:3, .N, replace = T, prob = c(0.05, 0.95))]
  
  d[silo == 2, pref_rev := sample(1:2, size = .N, replace = T, prob = g_pr_c_pref[2, 2:3])]
  d[silo == 2, d1 := sample(c(-98, -99), size = .N, replace = T, prob = g_pr_c_surg)]
  d[silo == 2 & d1 == -98, d1 := 1]
  # under revision, d1 has high probability of being the preference
  d[silo == 2 & d1 == -99 & pref_rev == 1, d1 := sample(2:3, .N, replace = T, prob = c(0.95, 0.05))]
  d[silo == 2 & d1 == -99 & pref_rev == 2, d1 := sample(2:3, .N, replace = T, prob = c(0.05, 0.95))]
  
  # late acute is a randomised assignment for surgical domain
  d[silo == 3, pref_rev := sample(1:2, size = .N, replace = T, prob = g_pr_l_pref[2, 2:3])]
  # 1:1 rand assignment to dair (3) or rev (4) 
  d[silo == 3, d1 := sample(c(-98, -99), size = .N, replace = T)]
  # selection surgical type dair (1), one(2), two(3)
  # dair is fixed
  d[silo == 3 & d1 == -98, d1 := 1]
  # under revision, d1 has high probability of being the preference
  d[silo == 3 & d1 == -99 & pref_rev == 1, d1 := sample(2:3, .N, replace = T, prob = c(0.95, 0.05))]
  d[silo == 3 & d1 == -99 & pref_rev == 2, d1 := sample(2:3, .N, replace = T, prob = c(0.05, 0.95))]
  
  
  # Antibiotic backbone duration - only applicable to units recv one-stage
  
  d[d1 == 1, d2 := 1]
  d[d1 == 3, d2 := 1]
  d[d1 == 2, d2 := sample(2:3, size = .N, replace = T)]
  
  # d[, ix := rbinom(.N, 1, 0.6)]
  # d[ix == 0, d2 := 1]
  # d[ix == 1, d2 := sample(2:3, size = .N, replace = T)]
  
  # Extended prophylaxis - only applicable to units recv two-stage
  
  d[d1 == 1, d3 := 1]
  d[d1 == 2, d3 := 1]
  d[d1 == 3, d3 := sample(2:3, size = .N, replace = T)]
  
  # Outcome
  d[, mu := mu]
  d[, b_silo := b_silo[silo]]
  d[, b_jnt := b_jnt[jnt]]
  d[, b_pref := b_pref[pref_rev]]
  d[, b_d1 := b_d1[d1]]
  d[, b_d2 := b_d2[d2]]
  d[, eta := mu + b_silo + b_jnt + b_pref + b_d1 + b_d2      ]
  
  d[, y := rbinom(.N, 1, plogis(eta))]
  
  list(
    d = d, 
    
    mu = mu, 
    b_silo = b_silo, b_jnt = b_jnt, b_pref = b_pref,  
    b_d1 = b_d1, b_d2 = b_d2, b_d3 = b_d3, b_d4 = b_d4, 
    
    K_s = K_s, K_j = K_j, K_p  = K_p,
    K_d1 = K_d1, K_d2 = K_d2 , K_d3 = K_d3, K_d4 = K_d4
  )
  
}



main <- function(){
  mu = 0
  b_silo = c(0.5, -0.3, -0.2)
  b_jnt = c(-0.4, 0.4)
  b_pref = c(0.2, -0.2)
  b_d1 = c(0.3, -0.2, -0.1)
  b_g1 = rbind(
    c(-0.60,  0.20,  0.40),
    c( 0.20,  0.00, -0.20),
    c( 0.40, -0.20, -0.20)
  )
  b_d2 = c(-0.2, -0.4, 0.6)
  b_d3 = c(-0.4, 0.1, 0.3)
  b_d4 = c(0.3, -0.1, -0.2)
  
  m1 <- cmdstanr::cmdstan_model("stan/model10.stan")
  ll <- get_sim_01_data(
    N = 1e6, 
    b_silo = b_silo, b_jnt = b_jnt, b_pref = b_pref,
    b_g1 = b_g1, b_d1 = b_d1,
    # 6wk, 12wk
    b_d2 = b_d2,
    # ext proph duration (non-rand, 12wk, 7day)
    b_d3 = b_d3,
    # ab choice (non-rand, no-rif, rif)
    b_d4 = b_d4
    )
  
  d_mod1 <- ll$d[, .(y = sum(y), n = .N, eta = round(unique(eta), 3)), 
                 keyby = .(silo, jnt, pref_rev, d1, d2)]
  d_mod1[, eta_obs := round(qlogis(y / n), 3)]
  d_mod1[, p_obs := y / n]
  d_mod1
  
  
  
  ld <- list(
    N = nrow(d_mod1), 
    
    y = d_mod1[, y], 
    n = d_mod1[, n], 
    
    silo = d_mod1[, silo], 
    jnt = d_mod1[, jnt],
    pref = d_mod1[, pref_rev],
    d1 = d_mod1[, d1],
    d2 = d_mod1[, d2],
    
    g1 = ((d_mod1$silo - 1) * ll$K_d1) + d_mod1$d1,
    
    nrXs = nrow(ll$X_s_s), ncXs = ncol(ll$X_s_s), Xsdes = ll$X_s_s, ss = rep(1, ncol(ll$X_s_s)),
    nrXj = nrow(ll$X_j_s), ncXj = ncol(ll$X_j_s), Xjdes = ll$X_j_s, sj = rep(1, ncol(ll$X_j_s)),
    nrXp = nrow(ll$X_p_s), ncXp = ncol(ll$X_p_s), Xpdes = ll$X_p_s, sp = rep(1, ncol(ll$X_p_s)),
    
    nrXd1 = nrow(ll$X_d1_s), ncXd1 = ncol(ll$X_d1_s), Xd1des = ll$X_d1_s, sd1 = rep(1, ncol(ll$X_d1_s)),
    nrXd2 = nrow(ll$X_d2_s), ncXd2 = ncol(ll$X_d2_s), Xd2des = ll$X_d2_s, sd2 = rep(1, ncol(ll$X_d2_s)),
    
    nrXg1 = nrow(ll$X_g1_s), ncXg1 = ncol(ll$X_g1_s), Xg1des = ll$X_g1_s, sg1 = rep(1, ncol(ll$X_g1_s)),
    
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
  m_post_s <- f1$draws(variables = c("bp"), format = "matrix")
  post_mu <- m_post_s %*% t(cbind(ll$X_p_s))
  rbind(
    par = ll$b_pref,
    post_mu = colMeans(post_mu)
  )
  
  # 
  m_post_s <- f1$draws(variables = c("bg1"), format = "matrix")
  post_mu <- m_post_s %*% t(cbind(ll$X_g1_s))
  rbind(
    par = as.numeric(t(ll$b_g1)),
    post_mu = colMeans(post_mu)
  )
  
  #
  m_post_s <- f1$draws(variables = c("bd1"), format = "matrix")
  post_mu <- m_post_s %*% t(cbind(ll$X_d1_s))
  rbind(
    par = ll$b_d1,
    post_mu = colMeans(post_mu)
  )
  
  
  
}


# based on the sum to zero contrained parameterisation
sim_01 <- function(){
  
  
  m1 <- cmdstanr::cmdstan_model("stan/model10.stan")
  
  mu = 0
  b_silo = c(0.5, -0.3, -0.2)
  b_jnt = c(-0.4, 0.4)
  b_pref = c(0.2, -0.2)
  b_d1 = c(0.3, -0.2, -0.1)
  b_g1 = rbind(
    c(-0.60,  0.20,  0.40),
    c( 0.20,  0.00, -0.20),
    c( 0.40, -0.20, -0.20)
  )
  b_d2 = c(0, 0, 0)
  b_d3 = c(-0.4, 0.1, 0.3)
  b_d4 = c(0.3, -0.1, -0.2)
  
  labs <- c("mu", 
            paste0("bs",seq_along(b_silo)),
            paste0("bj",seq_along(b_jnt)),
            paste0("bp",seq_along(b_pref)),
            paste0("bd1",seq_along(b_d1)),
            # ,
            "bd1_1", "bd1_2", "bd1_3",
            
            # paste0("bg1",seq_along(b_g1)),
            paste0("bd2",seq_along(b_d2)),
            "bd2_2", "bd2_3"
            # paste0("bd3",seq_along(b_d3))
  )
  
  n_sim <- 500
  d_res <- data.table(do.call(rbind, mclapply(1:n_sim, FUN = function(ii){
    
    ll <- get_sim_01_data(
      N = 2.5e3, 
      mu = mu,
      b_silo = b_silo,
      b_jnt = b_jnt,
      b_pref = b_pref,
      b_d1 = b_d1,
      b_g1 = b_g1,
      b_d2 = b_d2,  # 6wk, 12wk
      b_d3 = b_d3,   # ext proph duration (non-rand, 12wk, 7day)
      b_d4 = b_d4   # ab choice (non-rand, no-rif, rif)
    )
    
    d_mod1 <- ll$d[, .(y = sum(y), n = .N, eta = round(unique(eta), 3)), 
                   keyby = .(silo, jnt, pref_rev, d1, d2)]
    d_mod1[, eta_obs := round(qlogis(y / n), 3)]
    d_mod1[, p_obs := y / n]
    # d_mod1
    
    # only late silo (silo = 3) get rand surgery
    d_tmp <- d_mod1[silo == 3]
    X_pred_d1_1 <- matrix(NA, nrow(d_tmp), 2 + 1 + 1 + 2 + 2)
    X_pred_d1_2 <- matrix(NA, nrow(d_tmp), 2 + 1 + 1 + 2 + 2)
    X_pred_d1_3 <- matrix(NA, nrow(d_tmp), 2 + 1 + 1 + 2 + 2)
    for(i in 1:nrow(d_tmp)){

      X_pred_d1_1[i, ] <- c(
        ll$X_s_s[d_tmp[i, silo], ],      # 2
        ll$X_j_s[d_tmp[i, jnt], ],       # 1
        ll$X_p_s[d_tmp[i, pref_rev], ],  # 1
        # when I set this to 1 (dair)
        ll$X_d1_s[1, ],                  # 2
        # then anything other than d2 = 1 is not permitted.
        ll$X_d2_s[1, ]        # 2
      )

      X_pred_d1_2[i, ] <- c(
        ll$X_s_s[d_tmp[i, silo], ],      # 2
        ll$X_j_s[d_tmp[i, jnt], ],       # 1
        ll$X_p_s[d_tmp[i, pref_rev], ],  # 1
        # when I set this to 2 (one-stage)
        ll$X_d1_s[2, ],                  # 2
        # then d2 is permitted
        ll$X_d2_s[1, ]        # 2
      )

      X_pred_d1_3[i, ] <- c(
        ll$X_s_s[d_tmp[i, silo], ],      # 2
        ll$X_j_s[d_tmp[i, jnt], ],       # 1
        ll$X_p_s[d_tmp[i, pref_rev], ],  # 1
        # when I set this to 3 (two-stage)
        ll$X_d1_s[3, ],                  # 2
        # then anything other than d2 = 1 is not permitted.
        ll$X_d2_s[1, ]        # 2
      )
    }
    nd1p <- d_tmp$n
     
    # 
    # only one-stage (d1 = 2) get rand ab backbone duration (d2)
    d_tmp <- d_mod1[d1 == 2]
    X_pred_d2_2 <- matrix(NA, nrow(d_tmp), 2 + 1 + 1 + 2 + 2)
    X_pred_d2_3 <- matrix(NA, nrow(d_tmp), 2 + 1 + 1 + 2 + 2)
    for(i in 1:nrow(d_tmp)){

      X_pred_d2_2[i, ] <- c(
        ll$X_s_s[d_tmp[i, silo], ],      # 2
        ll$X_j_s[d_tmp[i, jnt], ],       # 1
        ll$X_p_s[d_tmp[i, pref_rev], ],  # 1
        ll$X_d1_s[d_tmp[i, d1], ],       # 2

        ll$X_d2_s[2, ]                   # 2
      )

      X_pred_d2_3[i, ] <- c(
        ll$X_s_s[d_tmp[i, silo], ],      # 2
        ll$X_j_s[d_tmp[i, jnt], ],       # 1
        ll$X_p_s[d_tmp[i, pref_rev], ],  # 1
        ll$X_d1_s[d_tmp[i, d1], ],       # 2
        ll$X_d2_s[3, ]                   # 2
      )
    }
    nd2p <- d_tmp$n
    
    
    ld <- list(
      N = nrow(d_mod1), 
      
      y = d_mod1[, y], 
      n = d_mod1[, n], 
      
      silo = d_mod1[, silo], 
      jnt = d_mod1[, jnt],
      pref = d_mod1[, pref_rev],
      d1 = d_mod1[, d1],
      d2 = d_mod1[, d2],
      
      # g1 = ((d_mod1$silo - 1) * ll$K_d1) + d_mod1$d1,
      
      nrXs = nrow(ll$X_s_s), ncXs = ncol(ll$X_s_s), Xsdes = ll$X_s_s, ss = rep(1, ncol(ll$X_s_s)),
      nrXj = nrow(ll$X_j_s), ncXj = ncol(ll$X_j_s), Xjdes = ll$X_j_s, sj = rep(1, ncol(ll$X_j_s)),
      nrXp = nrow(ll$X_p_s), ncXp = ncol(ll$X_p_s), Xpdes = ll$X_p_s, sp = rep(1, ncol(ll$X_p_s)),
      
      nrXd1 = nrow(ll$X_d1_s), ncXd1 = ncol(ll$X_d1_s), Xd1des = ll$X_d1_s, sd1 = rep(1, ncol(ll$X_d1_s)),
      nrXd2 = nrow(ll$X_d2_s), ncXd2 = ncol(ll$X_d2_s), Xd2des = ll$X_d2_s, sd2 = rep(1, ncol(ll$X_d2_s)),
      
      nrd1p = nrow(X_pred_d1_1), ncd1p = ncol(X_pred_d1_1),
      Xd1p1 = X_pred_d1_1, Xd1p2 = X_pred_d1_2, Xd1p3 = X_pred_d1_3,
      nd1p = nd1p,
      # 
      nrd2p = nrow(X_pred_d2_2), ncd2p = ncol(X_pred_d2_2),
      Xd2p1 = X_pred_d2_2, Xd2p2 = X_pred_d2_3,
      nd2p = nd2p,
      
      prior_only = 0
    )
    
    snk <- capture.output(
      f1 <- m1$pathfinder(ld, num_paths=20, single_path_draws=200,
                          history_size=50, max_lbfgs_iters=100,
                          refresh = 0, draws = 2000)
      )
    
    m_post_s <- f1$draws(variables = c("mu_pop"), format = "matrix")
    v_out <- c(colMeans(m_post_s))
    
    #
    m_post_s <- f1$draws(variables = c("bs"), format = "matrix")
    post_mu <- m_post_s %*% t(cbind(ll$X_s_s))
    v_out <- c(v_out,  colMeans(post_mu))
    
    # 
    m_post_s <- f1$draws(variables = c("bj"), format = "matrix")
    post_mu <- m_post_s %*% t(cbind(ll$X_j_s))
    v_out <- c(v_out,  colMeans(post_mu))
    
    #
    m_post_s <- f1$draws(variables = c("bp"), format = "matrix")
    post_mu <- m_post_s %*% t(cbind(ll$X_p_s))
    v_out <- c(v_out,  colMeans(post_mu))
    
    #
    m_post_s <- f1$draws(variables = c("bd1"), format = "matrix")
    post_mu <- m_post_s %*% t(cbind(ll$X_d1_s))
    v_out <- c(v_out,  colMeans(post_mu))
    
    
    m_post_s <- f1$draws(variables = c("bd1_1"), format = "matrix")
    v_out <- c(v_out,  colMeans(m_post_s))
    m_post_s <- f1$draws(variables = c("bd1_2"), format = "matrix")
    v_out <- c(v_out,  colMeans(m_post_s))
    m_post_s <- f1$draws(variables = c("bd1_3"), format = "matrix")
    v_out <- c(v_out,  colMeans(m_post_s))
    # 
    # #
    # # m_post_s <- f1$draws(variables = c("bg1"), format = "matrix")
    # # post_mu <- m_post_s %*% t(cbind(ll$X_g1_s))
    # # v_out <- c(v_out,  colMeans(post_mu))
    # 
    # #
    m_post_s <- f1$draws(variables = c("bd2"), format = "matrix")
    post_mu <- m_post_s %*% t(cbind(ll$X_d2_s))
    v_out <- c(v_out,  colMeans(post_mu))
    # 
    m_post_s <- f1$draws(variables = c("bd2_2"), format = "matrix")
    v_out <- c(v_out,  colMeans(m_post_s))
    m_post_s <- f1$draws(variables = c("bd2_3"), format = "matrix")
    v_out <- c(v_out,  colMeans(m_post_s))
    
    names(v_out) <- labs
    
    v_out
      
    }, mc.cores = 6)))
  
  d_res <- melt(d_res, measure.vars = names(d_res))
  d_res[, variable := factor(variable, levels = labs)]
  
  d_tru <- data.table(
    variable = labs,
    value = c(
      mu,
      b_silo,
      b_jnt,
      b_pref,
      b_d1,
      # ,
      b_d1,
      # as.numeric(t(b_g1)),
      b_d2,
      b_d2[2:3]
      # ,
      # b_d3
    )
  )
  d_tru[, variable := factor(variable, levels = labs)]
  
  ggplot(d_res, aes(x = value, group = variable)) +
    geom_density() +
    geom_vline(data = d_res[, .(value = mean(value)), keyby = variable],
               aes(xintercept = value), col = 1, lwd = 0.3) +
    geom_vline(data = d_tru,
               aes(xintercept = value), col = 2, lwd = 0.3, lty = 2) +
    facet_wrap(~variable, scales = "free_x")
  
  d_res[, .(mu = mean(value)), keyby = variable]
  
}





# based on the fixed reference level parameterisation
sim_02 <- function(){
  
  m1 <- cmdstanr::cmdstan_model("stan/model11.stan")
  
  mu = 0
  b_silo = c(0, -0.3, -0.2)
  b_jnt = c(0, 0.4)
  b_pref = c(0, -0.2)
  b_d1 = c(0, -0.5, 0.2)
  b_d2 = c(0, 0, 0)
  b_d3 = c(0, 0.1, 0.3)
  b_d4 = c(0, -0.1, -0.2)
  
  labs <- c("mu", 
            paste0("bs",seq_along(b_silo)),
            paste0("bj",seq_along(b_jnt)),
            paste0("bp",seq_along(b_pref)),
            paste0("bd1",seq_along(b_d1)),
            # ,
            "bd1_1", "bd1_2", "bd1_3",
            
            # paste0("bg1",seq_along(b_g1)),
            paste0("bd2",seq_along(b_d2)),
            "bd2_2", "bd2_3"
            # paste0("bd3",seq_along(b_d3))
  )
  
  n_sim <- 500
  d_res <- data.table(do.call(rbind, mclapply(1:n_sim, FUN = function(ii){
    
    ll <- get_sim_02_data(
      N = 2.5e3, 
      mu = mu,
      b_silo = b_silo,
      b_jnt = b_jnt,
      b_pref = b_pref,
      b_d1 = b_d1,
      b_d2 = b_d2,  # 6wk, 12wk
      b_d3 = b_d3,   # ext proph duration (non-rand, 12wk, 7day)
      b_d4 = b_d4   # ab choice (non-rand, no-rif, rif)
    )
    
    d_mod1 <- ll$d[, .(y = sum(y), n = .N, eta = round(unique(eta), 3)), 
                   keyby = .(silo, jnt, pref_rev, d1, d2)]
    d_mod1[, eta_obs := round(qlogis(y / n), 3)]
    d_mod1[, p_obs := y / n]
    # d_mod1
    
    # only late silo (silo = 3) get rand surgery
    d_d1_pred <- d_mod1[silo == 3]
     
    # only one-stage (d1 = 2) get rand ab backbone duration (d2)
    d_d2_pred <- d_mod1[d1 == 2]
    
    ld <- list(
      N = nrow(d_mod1), 
      
      y = d_mod1[, y], 
      n = d_mod1[, n], 
      
      K_silo = ll$K_s, K_jnt = ll$K_j, K_pref = ll$K_p, 
      K_d1 = ll$K_d1, K_d2 = ll$K_d2, 
      silo = d_mod1[, silo], 
      jnt = d_mod1[, jnt],
      pref = d_mod1[, pref_rev],
      d1 = d_mod1[, d1],
      d2 = d_mod1[, d2],
   
      N_ix_d1_1 = length(d_mod1[d1 == 1, which = T]),
      N_ix_d1_2 = length(d_mod1[d1 == 2, which = T]),
      N_ix_d1_3 = length(d_mod1[d1 == 3, which = T]),
      ix_d1_1 = d_mod1[d1 == 1, which = T],
      ix_d1_2 = d_mod1[d1 == 2, which = T],
      ix_d1_3 = d_mod1[d1 == 3, which = T],
      
      # nrd2p = nrow(X_pred_d2_2), ncd2p = ncol(X_pred_d2_2),
      # Xd2p1 = X_pred_d2_2, Xd2p2 = X_pred_d2_3,
      # nd2p = d_d2_pred$n,
      
      nrd1p = nrow(d_d1_pred),
      d1_s = d_d1_pred$silo,
      d1_j = d_d1_pred$jnt,
      d1_p = d_d1_pred$pref_rev,
      # d1 is intervened being set to index par 2 or 3
      # should d2 be logically constrained by d1 and the need for a fair comparison???
      d1_d2 = d_d1_pred$d2,
      # of like with like
      nd1p = d_d1_pred$n,
      
      nrd2p = nrow(d_d2_pred),
      d2_s = d_d2_pred$silo,
      d2_j = d_d2_pred$jnt,
      d2_p = d_d2_pred$pref_rev,
      d2_d1 = d_d2_pred$d1,
      # d2 is intervened being set to index par 2 or 3
      nd2p = d_d2_pred$n,
      
      prior_only = 0
    )
    
    snk <- capture.output(
      f1 <- m1$pathfinder(ld, num_paths=20, single_path_draws=200,
                          history_size=50, max_lbfgs_iters=100,
                          refresh = 0, draws = 2000)
    )
    
    m_post_s <- f1$draws(variables = c("mu"), format = "matrix")
    v_out <- c(colMeans(m_post_s))
    
    #
    m_post_s <- f1$draws(variables = c("bs"), format = "matrix")
    v_out <- c(v_out,  colMeans(m_post_s))
    
    # 
    m_post_s <- f1$draws(variables = c("bj"), format = "matrix")
    v_out <- c(v_out,  colMeans(m_post_s))
    
    #
    m_post_s <- f1$draws(variables = c("bp"), format = "matrix")
    v_out <- c(v_out,  colMeans(m_post_s))
    
    #
    m_post_s <- f1$draws(variables = c("bd1"), format = "matrix")
    v_out <- c(v_out,  colMeans(m_post_s))
   
    m_post_s <- f1$draws(variables = c("bd1_1"), format = "matrix")
    v_out <- c(v_out,  colMeans(m_post_s))
    m_post_s <- f1$draws(variables = c("bd1_2"), format = "matrix")
    v_out <- c(v_out,  colMeans(m_post_s))
    m_post_s <- f1$draws(variables = c("bd1_3"), format = "matrix")
    v_out <- c(v_out,  colMeans(m_post_s))
    # 
    # #
    # # m_post_s <- f1$draws(variables = c("bg1"), format = "matrix")
    # # post_mu <- m_post_s %*% t(cbind(ll$X_g1_s))
    # # v_out <- c(v_out,  colMeans(post_mu))
    # 
    # #
    m_post_s <- f1$draws(variables = c("bd2"), format = "matrix")
    v_out <- c(v_out,  colMeans(m_post_s))
    # 
    m_post_s <- f1$draws(variables = c("bd2_2"), format = "matrix")
    v_out <- c(v_out,  colMeans(m_post_s))
    m_post_s <- f1$draws(variables = c("bd2_3"), format = "matrix")
    v_out <- c(v_out,  colMeans(m_post_s))
    
    names(v_out) <- labs
    
    v_out
    
  }, mc.cores = 6)))
  
  d_res <- melt(d_res, measure.vars = names(d_res))
  d_res[, variable := factor(variable, levels = labs)]
  
  d_tru <- data.table(
    variable = labs,
    value = c(
      mu,
      b_silo,
      b_jnt,
      b_pref,
      b_d1,
      b_d1,
      # as.numeric(t(b_g1)),
      b_d2,
      b_d2[2:3]
      # ,
      # b_d3
    )
  )
  d_tru[, variable := factor(variable, levels = labs)]
  
  ggplot(d_res, aes(x = value, group = variable)) +
    geom_density() +
    geom_vline(data = d_res[, .(value = mean(value)), keyby = variable],
               aes(xintercept = value), col = 1, lwd = 0.3) +
    geom_vline(data = d_tru,
               aes(xintercept = value), col = 2, lwd = 0.3, lty = 2) +
    facet_wrap(~variable, scales = "free_x")
  
  d_res[, .(mu = mean(value)), keyby = variable]
  
  b_d1
  mean(d_res[variable == "bd1_2", value] - d_res[variable == "bd1_1", value])
  mean(d_res[variable == "bd1_3", value] - d_res[variable == "bd1_1", value])
  mean(d_res[variable == "bd1_3", value] - d_res[variable == "bd1_2", value])
  
  b_d2
  mean(d_res[variable == "bd2_3", value] - d_res[variable == "bd2_2", value])
  
}