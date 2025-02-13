# date 4/9/2024 - the focus is on sim_02.


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
    # dair, one, two-stage
    b_d1 = c(0, -0.1, -0.1),
    # b_d1_ex = c(0, -0.1, -0.1, -0.1),
    # 
    b_d2 = c(0, -0.4, 0.6),
    b_d3 = c(0, 0.1, 0.3),
    b_d4 = c(0, -0.1, -0.2)
){

  # number of silos, joints, preference
  K_s <- length(b_silo)
  K_j <- length(b_jnt)
  K_p <- length(b_pref)
  # number of parameters in each domain
  K_d1 <- length(b_d1)
  K_d2 <- length(b_d2)
  K_d3 <- length(b_d3)
  K_d4 <- length(b_d4)
  
  d <- data.table()
  
  # Stratification
  d[, silo := sample(1:K_s, size = N, replace = T, prob = g_pr_silo)]
  d[silo == 1, jnt := sample(1:K_j, size = .N, replace = T, prob = g_pr_jnt[1, ])]
  d[silo == 2, jnt := sample(1:K_j, size = .N, replace = T, prob = g_pr_jnt[2, ])]
  # NOTE - I have set silo 3 as the one for late acute
  d[silo == 3, jnt := sample(1:K_j, size = .N, replace = T, prob = g_pr_jnt[3, ])]
  
  # Domains
  
  # Surgery - assignment to surgery dictates pretty much everything else.
  
  # early - non-randomised
  # clinician decides on what surgery type occurs for silo 1 and 2
  # pref_rev refers to the preferred revision type, either 
  # one-stage or two-stage.
  d[silo == 1, pref_rev := sample(1:2, size = .N, replace = T, prob = g_pr_e_pref[2, 2:3])]
  d[silo == 1, d1 := sample(c(-98, -99), size = .N, replace = T, prob = g_pr_e_surg)]
  # dair
  d[silo == 1 & d1 == -98, d1 := 1]
  # or revision, here d1 is set has very high probability of being whatever the 
  # preference was
  d[silo == 1 & d1 == -99 & pref_rev == 1, d1 := sample(2:3, .N, replace = T, prob = c(0.95, 0.05))]
  d[silo == 1 & d1 == -99 & pref_rev == 2, d1 := sample(2:3, .N, replace = T, prob = c(0.05, 0.95))]
  
  # don't think this kind of approach (trying to adjust for non-randomised)
  # would work as the silos are entirely collinear
  # d[silo == 1 & d1 == 1, d1_ex := 2]
  # d[silo == 1 & d1 == 2, d1_ex := 3]
  # d[silo == 1 & d1 == 3, d1_ex := 4]
  
  # exactly same approach for chronic silo but using different probabilities
  d[silo == 2, pref_rev := sample(1:2, size = .N, replace = T, prob = g_pr_c_pref[2, 2:3])]
  d[silo == 2, d1 := sample(c(-98, -99), size = .N, replace = T, prob = g_pr_c_surg)]
  d[silo == 2 & d1 == -98, d1 := 1]
  # under revision, d1 has high probability of being the preference
  d[silo == 2 & d1 == -99 & pref_rev == 1, d1 := sample(2:3, .N, replace = T, prob = c(0.95, 0.05))]
  d[silo == 2 & d1 == -99 & pref_rev == 2, d1 := sample(2:3, .N, replace = T, prob = c(0.05, 0.95))]
  
  # for late acute we have a randomised assignment process for surgical domain
  # rather than selection
  # for the purposes here, the simulation approach is similar to that above
  # but the initial dair vs rev assignment is 1:1 see below
  d[silo == 3, pref_rev := sample(1:2, size = .N, replace = T, prob = g_pr_l_pref[2, 2:3])]
  # 1:1 rand assignment to dair (3) or rev (4) 
  d[silo == 3, d1 := sample(c(-98, -99), size = .N, replace = T)]
  # selection surgical type dair (1), one(2), two(3)
  # dair is fixed
  d[silo == 3 & d1 == -98, d1 := 1]
  # under revision, d1 has high probability of being the preference
  d[silo == 3 & d1 == -99 & pref_rev == 1, d1 := sample(2:3, .N, replace = T, prob = c(0.95, 0.05))]
  d[silo == 3 & d1 == -99 & pref_rev == 2, d1 := sample(2:3, .N, replace = T, prob = c(0.05, 0.95))]
  
  # entry indictor for d1
  # d[silo != 3 , g1 := 1]
  # d[silo == 3 , g1 := 2]
  
  # Antibiotic backbone duration - only applicable to units recv one-stage
  # all others just get set to the reference index.
  # Comparisons are made directly via the full linear predictor rather than
  # relying on specific parameters.
  d[d1 %in% c(1, 3), d2 := 1]
  
  # We have to have some small group with d1 = 2 and d2 = 1, i.e. just receiving 
  # the hospital/surgeon select duration and therefore not receiving either
  # of the specific randomised treatment arms (6wk vs 12wk). If we do not make this
  # assumption we will hit an identification problem and that would mean that we 
  # could not recover the true parameters and our results without fairly major bias.
  # We assume that of those receiving one-stage revision,
  # 70% enter into the duration domain (the remainder being assigned to the ref index)
  d[d1 %in% c(2), d2 := sample(1:3, size = .N, replace = T, prob = c(0.3, 0.35, 0.35))]
  
  # indicator of entry into d2
  # d[d1 %in% c(1, 3), g2 := 1]
  # d[d1 %in% c(2), g2 := 2]
  
  # d[, ix := rbinom(.N, 1, 0.6)]
  # d[ix == 0, d2 := 1]
  # d[ix == 1, d2 := sample(2:3, size = .N, replace = T)]
  
  # Extended prophylaxis - only applicable to units recv two-stage
  
  d[d1 %in% c(1, 2), d3 := 1]
  
  # Again, we have to have some small group with d1 = 3 and d2 = 1 otherwise we have
  # an identification problem that would mean that we cannot recover the true
  # parameters and our results will be biased.
  # Here we assume that of those receiving two-stage revision,
  # 90% enter into the Extended prophylaxis domain
  d[d1 %in% c(3), d3 := sample(1:3, size = .N, replace = T, prob = c(0.1, 0.45, 0.45))]
  
  # indicator of entry into d3
  # d[d1 %in% c(1, 2), g3 := 1]
  # d[d1 %in% c(3), g3 := 2]
  
  # Choice domain is independent to others but only 60% of the population
  # enter it.
  
  d[, g4 := sample(1:2, size = .N, replace = T, prob = c(0.4, 0.6))]
  d[g4 == 1, d4 := 1]
  d[g4 == 2, d4 := sample(2:3, size = .N, replace = T)]
  d[, g4 := NULL]
  
  # Outcome
  d[, mu := mu]
  d[, b_silo := b_silo[silo]]
  d[, b_jnt := b_jnt[jnt]]
  d[, b_pref := b_pref[pref_rev]]
  d[, b_d1 := b_d1[d1]]
  d[, b_d2 := b_d2[d2]]
  d[, b_d3 := b_d3[d3]]
  d[, b_d4 := b_d4[d4]]
  d[, eta := mu + b_silo + b_jnt + b_pref + b_d1 + b_d2 + b_d3 + b_d4      ]
  
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




# 2024-09-16 - current approach being tested.
# based on the fixed reference level parameterisation
sim_02 <- function(){
  
  m11 <- cmdstanr::cmdstan_model("stan/model11.stan")
  
  mc_cores <- 1
  if(unname(Sys.info()[1]) == "Darwin"){
    log_info("On mac, reset cores to 5")
    mc_cores <- 5
  } else if (unname(Sys.info()[1]) == "Linux"){
    log_info("On linux, reset cores to 40")
    mc_cores <- 40
  }
  
  mu = 1
  b_silo = c(0, -0.3, 0.2)
  b_jnt = c(0, -0.4)
  # What is the interpretation of this? That those for whom the preference 
  # for revision type is two stage have a lower log odds of trt success.
  b_pref = c(0, -0.5)
  
  # For the intervention parameters the way that they enter the model is 
  # somewhat irrelevant in that all the comparisons of interest are derived
  # via g-computation where we are making a prediction on each unit and 
  # then averaging over these for the comparisons we need to report.
  
  # Assume that everyone gets either dair, one or two (no-one gets no surgery). 
  # Within d1 we do not differentiate between randomised and non-randomised 
  # assignment on the assumption that revision preference and additional 
  # covariates in the model would give us exchangeable groups.
  # It isn't entirely clear (to me) how we will adjust for non-randomised 
  # intervention in the present setup given that non-randomised surgical treatment
  # is collinear with silo membership, i.e. after accounting for silo, there is 
  # no further variance to account for (not under our current data generation 
  # assumptions anyway).
  # In the final model we may attempt this by including another vector of 
  # parameters for d1 to 
  # adjust for all of non-rand dair, one, two.
  # In the g-comp for other domains, e.g. d2, you would incl this parameter set 
  # on a unit by unit basis.
  # This would be setup as a vector with four entries, the first being set to zero
  # corresponding to randomised treatment and the subsequent being adjustment 
  # for non-randomised dair, one and two respectively.
  # However, I am skeptical as to whether this would work and am not 
  # investigating it further.
  
  # Note the bias we are seeing at the moment is purely a function of the sample
  # size and the small number of simulations. If you increase the sample size 
  # you get consistency and unbiasedness. The problem with small sample sizes
  # is that you get a lot of extremes where the observed/empircal log-odds 
  # would be infinite in a number of the covariate groups.
  
  b_d1 = c(0, 0.8, 0.8)
  # non-randomised surgical intervention adj
  # b_d1_ex = c(0, 0.1, -0.3, 0.25)
  # non-randomised entries correspond to reference level
  b_d2 = c(0, -0.5, 2)
  b_d3 = c(0, 0.1, -0.3)
  b_d4 = c(0, 0.1, -0.2)
  #
  
  labs <- c(
    "gbd1",
    "gbd2",
    "gbd3",
    "gbd4"
  )
  
  n_sim <- 5000
  
  d_res <- data.table(do.call(rbind, mclapply(1:n_sim, FUN = function(ii){
  
    ll <- get_sim_02_data(
      N = 2.5e3,
      # 
      # N = 1e6,
      mu = mu,
      b_silo = b_silo,
      b_jnt = b_jnt,
      b_pref = b_pref,
      b_d1 = b_d1,
      b_d2 = b_d2,  # 6wk, 12wk
      b_d3 = b_d3,   # ext proph duration (non-rand, 12wk, 7day)
      b_d4 = b_d4   # ab choice (non-rand, no-rif, rif)
    )
    
    # convert from binary representation to binomial (successes/trials)
    d_mod <- ll$d[, .(y = sum(y), n = .N, eta = round(unique(eta), 3)), 
                   keyby = .(silo, jnt, pref_rev, d1, d2, d3, d4)]
    d_mod[, eta_obs := round(qlogis(y / n), 3)]
    d_mod[, p_obs := y / n]
    
    # really don't need to do this for all of the domains but it provides
    # a uniform approach for generating the parameters of interest
    
    # restrict to silo 3 (which I have set to be the late acute silo) for gcomp
    # for d1 randomised comparisons
    # here I just use the silo assignment to imply the right subset of units
    # but in practice there may need to be an indicator variable to for this.
    d_mod_d1 <- d_mod[silo == 3]
    
    # restrict to one-stage for gcomp
    # here I just use the treatment assignment to imply the right subset of units
    # but in practice there may need to be an indicator variable to for this.
    d_mod_d2 <- d_mod[d1 == 2]
    
    # restrict to two-stage for gcomp
    # here I just use the treatment assignment to imply the right subset of units
    # but in practice there may need to be an indicator variable to for this.
    d_mod_d3 <- d_mod[d1 == 3]
    
    # restrict to d4 randomised group 
    # here I just use the treatment assignment to imply the right subset of units
    # but in practice there may need to be an indicator variable to for this.
    d_mod_d4 <- d_mod[d4 %in% 2:3]
    
    ld <- list(
      N = nrow(d_mod), 
      
      y = d_mod[, y], 
      n = d_mod[, n], 
      
      K_silo = ll$K_s, K_jnt = ll$K_j, K_pref = ll$K_p, 
      K_d1 = ll$K_d1, K_d2 = ll$K_d2, K_d3 = ll$K_d3, K_d4 = ll$K_d4, 
      silo = d_mod[, silo], 
      jnt = d_mod[, jnt],
      pref = d_mod[, pref_rev],
      d1 = d_mod[, d1],
      d2 = d_mod[, d2],
      d3 = d_mod[, d3],
      d4 = d_mod[, d4],
      
      # cohort for surgical domain g-comp complicated due to one-stage/two-stage
      # considerations
      N_d1 = nrow(d_mod_d1),
      d1_s = d_mod_d1[, silo], d1_j = d_mod_d1[, jnt], d1_p = d_mod_d1[, pref_rev],
      # all d1 assignments
      d1_d1 = d_mod_d1[, d1], 
      d1_d2 = d_mod_d1[, d2],
      d1_d3 = d_mod_d1[, d3],
      d1_d4 = d_mod_d1[, d4],
      # number of trials within each strata
      n_d1 = d_mod_d1[, n],
      # sample size of those where preference is for one-stage
      N_d1_p1 = d_mod_d1[pref_rev == 1, .N],
      # sample size of those where preference is for two-stage
      N_d1_p2 = d_mod_d1[pref_rev == 2, .N],
      # indexes for those with preference for one-stage
      ix_d1_p1 = d_mod_d1[pref_rev == 1, which = T],
      # indexes for those with preference for two-stage
      ix_d1_p2 = d_mod_d1[pref_rev == 2, which = T],
      # number of trials within each of these subsets
      n_d1_p1 = d_mod_d1[pref_rev == 1, n],
      n_d1_p2 = d_mod_d1[pref_rev == 2, n],
      
      
      N_d2 = nrow(d_mod_d2),
      d2_s = d_mod_d2[, silo], d2_j = d_mod_d2[, jnt], d2_p = d_mod_d2[, pref_rev],
      d2_d1 = d_mod_d2[, d1], 
      d2_d2 = d_mod_d2[, d2],
      d2_d3 = d_mod_d2[, d3],
      d2_d4 = d_mod_d2[, d4],
      n_d2 = d_mod_d2[, n],
      
      N_d3 = nrow(d_mod_d3),
      d3_s = d_mod_d3[, silo], d3_j = d_mod_d3[, jnt], d3_p = d_mod_d3[, pref_rev],
      d3_d1 = d_mod_d3[, d1], 
      d3_d2 = d_mod_d3[, d2],
      d3_d3 = d_mod_d3[, d3],
      d3_d4 = d_mod_d3[, d4],
      n_d3 = d_mod_d3[, n],
      
      N_d4 = nrow(d_mod_d4),
      d4_s = d_mod_d4[, silo], d4_j = d_mod_d4[, jnt], d4_p = d_mod_d4[, pref_rev],
      d4_d1 = d_mod_d4[, d1], 
      d4_d2 = d_mod_d4[, d2],
      d4_d3 = d_mod_d4[, d3],
      d4_d4 = d_mod_d4[, d4],
      n_d4 = d_mod_d4[, n],
      
      prior_only = 0
    )
    
    # pathfinder is actually a bit biased for logistic but for quick and 
    # dirty ballpark it is ok.
    # snk <- capture.output(
    #   f1 <- m1$pathfinder(ld, num_paths=20, single_path_draws=200,
    #                       history_size=50, max_lbfgs_iters=100,
    #                       refresh = 0, draws = 2000)
    # )
    
    f1 <- m11$sample(
      ld, iter_warmup = 1000, iter_sampling = 2000,
      parallel_chains = 1, chains = 1, refresh = 0, show_exceptions = F, 
      max_treedepth = 13)
    
    d_post <- data.table(f1$draws(variables = c("lor_d1", "lor_d2", "lor_d3", "lor_d4"), format = "matrix"))
    colMeans(d_post)

    
  }, mc.cores = mc_cores)))
  
  
  message("Done - up to here...", Sys.time())
  # 
  
  d_fig <- melt(d_res, measure.vars = names(d_res))
  d_tru <- data.table(
    variable = paste0("lor_d", 1:4),
    value = c(
      b_d1[3],
      b_d2[3] - b_d2[2],
      b_d3[3] - b_d3[2],
      b_d4[3] - b_d4[2]
    )
  )
  
  ggplot(d_fig, aes(x = value, group = variable)) +
    geom_density() +
    geom_vline(data = d_fig[, .(value = mean(value)), keyby = variable],
               aes(xintercept = value, group = variable)) +
    geom_vline(data = d_tru, aes(xintercept = value, group = variable), col = 2) +
    facet_wrap(~variable, scales = "free_x")
  
  
  
  
  
}




# identification
demo_identification <- function(){
  
  m1 <- cmdstanr::cmdstan_model("stan/model12.stan")
  N <- 2e7
  
  set.seed(1)
  dd <- data.table(
    s = sample(1:3, size = N, replace = T)
  )
  dd[s == 1, w := 1]
  # We need to have some of the s = 2 group be allocated to the first level of w
  # otherwise we cannot recover the true parameters.
  dd[s == 2, w := sample(1:3, size = .N, replace = T, prob = c(0.1,0.45,0.45))]
  dd[s == 3, w := 1]
  
  b_s = c(0, -2, 2)
  b_w = c(0, 0.5, -0.3)
  dd[, eta := 0 + b_s[s] + b_w[w]]
  dd[, y := rbinom(.N, 1, prob = plogis(eta))]
  # unique(dd[s == 1])
  # unique(dd[s == 2])
  # unique(dd[s == 3])
  # unique(dd[, .(s, w, eta)])
  dd[, s := factor(s)]
  dd[, w := factor(w)]
  d_mod <- dd[, .(y = sum(y), n = .N), keyby = .(s, w)]
  
  #
  ld <- list(
    N = nrow(d_mod), y = d_mod$y, n = d_mod$n, s = d_mod$s, w = d_mod$w, 
    K_s = length(unique(d_mod$s)), K_w = length(unique(d_mod$w)), 
    N_ix_s_1 = nrow(d_mod[s == 1]), 
    N_ix_s_2 = nrow(d_mod[s == 2]), 
    N_ix_s_3 = nrow(d_mod[s == 3]),
    ix_s_1 = d_mod[s == 1, which = T],
    ix_s_2 = d_mod[s == 2, which = T],
    ix_s_3 = d_mod[s == 3, which = T],
    prior_only = 0
  )
  snk <- capture.output(
    f1 <- m1$pathfinder(ld, num_paths=20, single_path_draws=200,
                        history_size=50, max_lbfgs_iters=100,
                        refresh = 0, draws = 2000)
  )
  m_post_s <- f1$draws(variables = c("eta"), format = "matrix")
  # can get back to the probability of response for each combination 
  # present in the data
  v_out <- c(colMeans(m_post_s))
  cbind(unique(dd[, .(s, w, eta)][order(s, w)]),
        round(v_out, 3))
  
  # but the true parameters will never be recovered because we simply 
  # don't have the information necessary to do so.
  m_post_s <- f1$draws(variables = c("bs", "bw"), format = "matrix")
  v_out <- c(colMeans(m_post_s))
  rbind(round(v_out, 3), 
        c(b_s, b_w))
  
  m_post_s <- f1$draws(variables = c("bs_delta", "bw_delta"), format = "matrix")
  message("bs_delta gives one-stage rel to dair and two-stage rel to dair")
  message("bw_delta is randomised 12wk vs 6wk or vice versa")
  v_out <- c(colMeans(m_post_s))
  rbind(est = round(v_out, 3),
        tru = c(b_s[2:3], b_w[3] - b_w[2]))

}





# NOTES FROM JT

#   I've had a brief look, just up to the null case example.
# One thing that jumps out at me is  the calculations of cond_p1 and cond_p2.
# You have:
# 
#   // one
#   cond_p1 = inv_logit(a0 + b1[1] + X2*b2 +         X4*b4); 
#   // two
#   cond_p2 = inv_logit(a0 + b1[2] + X2*b2 + X3*b3 + X4*b4);      
# 
# If we are just "intervening" on surgery type, i.e. assigning patients to either DAIR or revision,
# then cond_p1 would only apply to patients who would receive a one-stage if assigned to revision.
# Similarly, cond_p2 would only apply to patients who would receive a two-stage if assigned to revision.
# 
# It seems like you give everyone a one-stage and then give everyone a two-stage irrespective of preferences.
# Instead, if "revision" is the intervention, then only individuals for whom a one-stage was preferred would have gotten a one-stage 
# if assigned to "revision" instead of DAIR similarly for two-stage.
# 
# For mine, this would be something like (just doing it in R), for individual risks
# 
# p_dair_i <- plogis(a0 + X2 %*% b2)
# id_one <- which(X2[, 1] == 0) # X2[, 1] is preferred surgery type, 0 for one-stage
# p_one_i <- plogis(a0 + b1[1] + cbind(X2, X4)[id_one, ] %**% c(b2, b4))
# id_two <-  which(X2[, 1] == 1)# X2[, 1] is preferred surgery type, 1 for two-stage
# p_two_i <- plogis(a0 + b1[2] + cbind(X2, X3, X4)[id_two, ] %**% c(b2, b3, b4))
# 
# Weighted-risks then being
# 
# p_dair <- mean(p_dair_i)
# p_one <- sum(d_s[id_one, n] * p_one_i / sum(d_s[id_one, n]))
# p_two <- sum(d_s[id_two, n] * p_two_i / sum(d_s[id_two, n]))
# 
# And the risk under revision (the intervention)
# 
# p_rev <- sum(d_s[id_one, n]) / sum(d_s[, n]) * p_one + (1 - sum(d_s[id_one, n]) / sum(d_s[, n])) * p_two
# 
# Does this make sense?


rand1 <- function(){
  
  factorial_design <- CJ(
    d1 = c("dair", "rev"),
    d2 = c("6wk", "12wk"),
    d3 = c("0", "12wk"),
    d4 = c("no-rif", "rif")
  )
  factorial_design
  N <- 2000
  ix <- sample(1:nrow(factorial_design), size = N, replace = T)
  N_site <- 50
  site <- sample(1:N_site, N, replace = T, prob = as.numeric(rdirichlet(1, rep(1, N_site))))
  rev_pref <- sample(c("one-stage", "two-stage"), N, replace = T, prob = c(2/3, 1/3))
  
  # observed data
  d <- data.table(factorial_design[ix])  
  d[, `:=`(site = site, rev_pref = rev_pref)]
  head(d)
  
  # at the margins these are fine  
  d[, .N, keyby = .(site, rev_pref, d1)]
  d[, .N, keyby = .(site, rev_pref, d2)]
  d[, .N, keyby = .(site, rev_pref, d3)]
  d[, .N, keyby = .(site, rev_pref, d4)]
  
  # now consider based on what would be revealed
  
  merge(CJ(site = 1:N_site, d2 = unique(factorial_design$d2)), 
        # d2 is only relevant for those receiving one-stage
        d[d1 != "dair" & rev_pref == "one-stage", .N, keyby = .(site, d2)],
        all = T)
  
  merge(CJ(site = 1:N_site, d3 = unique(factorial_design$d3)), 
        # d3 is only relevant for those receiving two-stage
        d[d1 != "dair" & rev_pref == "two-stage", .N, keyby = .(site, d3)],
        all = T)
  
  
  
}




tmp <- function(){

  qlogis(0.6)
    
  lo <- rlogis(1e5, 0.4, 1)
  plot(density(lo))
  p <- plogis(lo)
  plot(density(p))
  
  
  
}
