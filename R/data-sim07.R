library(data.table)

# Data generation script takes over from the roadmap.data pkg as of 2024-12-11.


# late silo has been put last here
g_silo = c("e", "l", "c") 
g_pr_silo <- c(0.3, 0.5, 0.2)
names(g_pr_silo) <- g_silo

# matrix of the distribution of implicated joint by silo
g_pr_jnt <- matrix(
  c(0.4, 0.6, 0.7, 0.3, 0.5, 0.5), 3, 2, byrow = T
)
colnames(g_pr_jnt) <- c("knee", "hip")
rownames(g_pr_jnt) <- g_silo

# early silo randomisation probabilities to dair/rev
g_pr_e_surg <- c(0.85, 0.15)
names(g_pr_e_surg) <- c("dair", "rev")
# early silo preference to surgical type
g_pr_e_pref <- rbind(
  # for those that receive dair, the preferences are
  dair = c(0.85, 0.1, 0.05),
  # for those that receive rev, the preferences are spread across the 
  # revision options (one-stage or two-stage)
  rev = c(0, 2/3, 1/3)
)

# late silo randomisation probabilities to dair/rev
g_pr_l_surg <- c(0.5, 0.5)
names(g_pr_l_surg) <- c("dair", "rev")
# late silo preference to surgical type
g_pr_l_pref <- rbind(
  # for those that receive dair, the preferences are
  dair = c(0.2, 0.24, 0.56),
  # for those that receive rev, the preferences are spread across the 
  # revision options (one-stage or two-stage)
  rev = c(0, 0.3, 0.7)
)

# as above but for chronic silo
g_pr_c_surg <- c(0.2, 0.8)
names(g_pr_c_surg) <- c("dair", "rev")
g_pr_c_pref <- rbind(
  dair = c(0.2, 0.2, 0.6),
  rev = c(0, 0.25, 0.75)
)



get_sim07_trial_data <- function(
    l_spec
    ){
  
  if(is.null(l_spec$ia)){
    ia <- 1
  } else {
    ia <- l_spec$ia
  }
  if(is.null(l_spec$is)){
    is <- 1
  } else {
    is <- l_spec$is
  }
  if(is.null(l_spec$ie)){
    ie <- 1
  } else {
    ie <- l_spec$ie
  }
  if(is.null(l_spec$t0)){
    t0 <- 1
  } else {
    t0 <- l_spec$t0
  }
  
  d <- data.table(
    ia = ia,
    id = is:ie,
    t0 = t0,
    s = sample(1:3, l_spec$N[ia], replace = T, prob = l_spec$p_s_alloc)
  )
  d[s == 1, `:=`(
    # ctl/trt allocations ignoring dependencies
    d1_alloc = rbinom(.N, 1, l_spec$l_e$p_d1_alloc),
    # 70% entery d2
    d2_entry = rbinom(.N, 1, l_spec$l_e$p_d2_entry),
    d2_alloc = rbinom(.N, 1, l_spec$l_e$p_d2_alloc),
    # 90% enter d3
    d3_entry = rbinom(.N, 1, l_spec$l_e$p_d3_entry),
    d3_alloc = rbinom(.N, 1, l_spec$l_e$p_d3_alloc),
    # 60% enter d4
    d4_entry = rbinom(.N, 1, l_spec$l_e$p_d4_entry),
    d4_alloc = rbinom(.N, 1, l_spec$l_e$p_d4_alloc),
    # preference directs type of revision (0 rev(1), 1 rev(2))
    pref = rbinom(.N, 1, l_spec$l_e$p_pref)  + 1 
  )]
  d[s == 2, `:=`(
    # ctl/trt allocations ignoring dependencies
    d1_alloc = rbinom(.N, 1, l_spec$l_l$p_d1_alloc),
    # 70% entery d2
    d2_entry = rbinom(.N, 1, l_spec$l_l$p_d2_entry),
    d2_alloc = rbinom(.N, 1, l_spec$l_l$p_d2_alloc),
    # 90% enter d3
    d3_entry = rbinom(.N, 1, l_spec$l_l$p_d3_entry),
    d3_alloc = rbinom(.N, 1, l_spec$l_l$p_d3_alloc),
    # 60% enter d4
    d4_entry = rbinom(.N, 1, l_spec$l_l$p_d4_entry),
    d4_alloc = rbinom(.N, 1, l_spec$l_l$p_d4_alloc),
    # preference directs type of revision (0 rev(1), 1 rev(2))
    pref = rbinom(.N, 1, l_spec$l_l$p_pref)  + 1 
  )]
  d[s == 3, `:=`(
    # ctl/trt allocations ignoring dependencies
    d1_alloc = rbinom(.N, 1, l_spec$l_c$p_d1_alloc),
    # 70% entery d2
    d2_entry = rbinom(.N, 1, l_spec$l_c$p_d2_entry),
    d2_alloc = rbinom(.N, 1, l_spec$l_c$p_d2_alloc),
    # 90% enter d3
    d3_entry = rbinom(.N, 1, l_spec$l_c$p_d3_entry),
    d3_alloc = rbinom(.N, 1, l_spec$l_c$p_d3_alloc),
    # 60% enter d4
    d4_entry = rbinom(.N, 1, l_spec$l_c$p_d4_entry),
    d4_alloc = rbinom(.N, 1, l_spec$l_c$p_d4_alloc),
    # preference directs type of revision (0 rev(1), 1 rev(2))
    pref = rbinom(.N, 1, l_spec$l_c$p_pref)  + 1 
  )]
  
  # dair gets dair, revision gets split
  d[d1_alloc == 0, d1 := 1]
  d[d1_alloc == 1 & pref == 1, d1 := 2]
  d[d1_alloc == 1 & pref == 2, d1 := 3]
  
  d[d1 == 2 & d2_entry == 0, d2 := 1]
  d[d1 == 2 & d2_entry == 1, d2 := 2 + d2_alloc]
  
  d[d1 == 3 & d3_entry == 0, d3 := 1]
  d[d1 == 3 & d3_entry == 1, d3 := 2 + d3_alloc]
  
  d[d4_entry == 0, d4 := 1]
  d[d4_entry == 1, d4 := 2 + d4_alloc]
  
  # d[, .N, keyby = .(s, pref, d1, d2, d3, d4)]
  
  bd1 <- c(l_spec$l_e$bd1, l_spec$l_l$bd1, l_spec$l_c$bd1)
  # index for d1 is a function of 
  # d1 (3 levels - dair, rev(1), rev(2)) and silo membership (1:3)
  d[, d1_ix := d1 + (3 * (s - 1))]

  # bd1 is irrelevant since d1 == 1 is the ref group, fixed at zero
  d[d1 == 1, eta := l_spec$mu + l_spec$bs[s] + l_spec$bp[pref] + l_spec$bd4[d4]]
  # pref is irrelevant as d1 = 2 only occurs if pref = 0
  d[d1 == 2, eta := l_spec$mu + l_spec$bs[s] + bd1[d1_ix] + l_spec$bd2[d2] + l_spec$bd4[d4]]
  # but here pref is relevant as d1 = 3 only if pref = 1
  d[d1 == 3, eta := l_spec$mu + l_spec$bs[s] + l_spec$bp[pref] + bd1[d1_ix] + l_spec$bd3[d3] + l_spec$bd4[d4]]
  
  d[, `:=`(s = factor(s), 
           d1 = factor(d1), 
           d2 = factor(d2, levels = 1:4), 
           d3 = factor(d3, levels = 1:4), 
           d4 = factor(d4))]
  
  d[, p := plogis(eta)]
  d[, y := rbinom(.N, 1, p)]
  
  d
}


get_sim07_stan_data <- function(d_all){
  
  # convert from binary representation to binomial (successes/trials)
  d_mod <- d_all[, .(y = sum(y), n = .N, eta = round(unique(eta), 3)), 
                 keyby = .(s, pref, d1, d2, d3, d4)]
  
  
  d_mod[, `:=`(
    s = as.integer(s),
    d1 = as.integer(d1),
    d2 = as.integer(d2),
    d3 = as.integer(d3),
    d4 = as.integer(d4)
  )]
  
  d_mod[is.na(d2), d2 := 999]
  d_mod[is.na(d3), d3 := 999]
  
  K_d1 <- length(unique(d_all$d1))
  d_mod[, d1_ix := d1 + (K_d1 * (s - 1))]
  
  d_mod[, eta_obs := qlogis(y / n)]
  d_mod[, p_obs := y / n]
  
  # g-comp setup for all of the domains to provide
  # a uniform approach for generating the parameters of interest.
  
  # Surgical -----
  
  # restrict to silo 2 (late acute silo) for gcomp for d1 randomised comparisons
  # here I just use the silo assignment to imply the right subset of units
  # but in practice there may need to be an indicator variable to for this.
  d_mod_d1 <- d_mod[s == 2]
  
  
  # AB duration ----
  
  # restrict to one-stage for gcomp
  # surgery assignment being one-stage (d1 at level 2)
  # permits entry
  d_mod_d2 <- d_mod[d1 == 2]
  
  # Ext proph duration ----
  
  # restrict to two-stage for gcomp - 
  # surgery assignment being two-stage (d1 at level 3)
  # permits entry
  d_mod_d3 <- d_mod[d1 == 3]
  
  # AB choice ----
  
  # restrict to d4 randomised group - levels 2 and 3 are the randomised groups.
  # The units having d4 set to 1 were not included in ab choice
  d_mod_d4 <- d_mod[d4 %in% 2:3]
  
  ld <- list(
    # full dataset
    N = nrow(d_mod), 
    y = d_mod[, y], 
    n = d_mod[, n], 
    s = d_mod[, s], 
    pref = d_mod[, pref],
    d1 = d_mod[, d1],
    d2 = d_mod[, d2],
    d3 = d_mod[, d3],
    d4 = d_mod[, d4],
    
    # Number of levels for silos, joints, pref and each trt.
    K_s = length(unique(d_mod$s)), 
    K_p = length(unique(d_mod$pref)), 
    K_d1 = d_mod[, length(unique(d1))], 
    K_d2 = d_mod[d2 != 999, length(unique(d2))], 
    K_d3 = d_mod[d3 != 999, length(unique(d3))], 
    K_d4 = d_mod[, length(unique(d4))], 
    
    
    # cohort for surgical domain g-comp complicated due to one-stage/two-stage
    # considerations
    N_d1 = nrow(d_mod_d1),
    d1_s = d_mod_d1[, s], 
    d1_p = d_mod_d1[, pref],
    # all d1 assignments
    d1_d1 = d_mod_d1[, d1], 
    d1_d2 = d_mod_d1[, d2],
    d1_d3 = d_mod_d1[, d3],
    d1_d4 = d_mod_d1[, d4],
    # number of trials within each strata
    n_d1 = d_mod_d1[, n],
    # sample size of those where preference is for one-stage
    N_d1_p1 = d_mod_d1[pref == 1, .N],
    # sample size of those where preference is for two-stage
    N_d1_p2 = d_mod_d1[pref == 2, .N],
    # indexes for those with preference for one-stage
    ix_d1_p1 = d_mod_d1[pref == 1, which = T],
    # indexes for those with preference for two-stage
    ix_d1_p2 = d_mod_d1[pref == 2, which = T],
    # number of trials within each of these subsets
    n_d1_p1 = d_mod_d1[pref == 1, n],
    n_d1_p2 = d_mod_d1[pref == 2, n],
    
    prop_p1 = d_all[, .(wgt = .N/nrow(d_all)), keyby = pref][pref == 1, wgt],
    prop_p2 = d_all[, .(wgt = .N/nrow(d_all)), keyby = pref][pref == 2, wgt],
    
    N_d2 = nrow(d_mod_d2),
    d2_s = d_mod_d2[, s], 
    d2_p = d_mod_d2[, pref],
    d2_d1 = d_mod_d2[, d1], 
    d2_d2 = d_mod_d2[, d2],
    d2_d3 = d_mod_d2[, d3],
    d2_d4 = d_mod_d2[, d4],
    n_d2 = d_mod_d2[, n],
    
    N_d3 = nrow(d_mod_d3),
    d3_s = d_mod_d3[, s], 
    d3_p = d_mod_d3[, pref],
    d3_d1 = d_mod_d3[, d1], 
    d3_d2 = d_mod_d3[, d2],
    d3_d3 = d_mod_d3[, d3],
    d3_d4 = d_mod_d3[, d4],
    n_d3 = d_mod_d3[, n],
    
    N_d4 = nrow(d_mod_d4),
    d4_s = d_mod_d4[, s], 
    d4_p = d_mod_d4[, pref],
    d4_d1 = d_mod_d4[, d1], 
    d4_d2 = d_mod_d4[, d2],
    d4_d3 = d_mod_d4[, d3],
    d4_d4 = d_mod_d4[, d4],
    n_d4 = d_mod_d4[, n],
    
    prior_only = 0
  )
  
  list(
    d_mod = d_mod,
    d_mod_d1 = d_mod_d1,
    d_mod_d2 = d_mod_d2,
    d_mod_d3 = d_mod_d3,
    d_mod_d4 = d_mod_d4,
    
    ld = ld
  )
  
}






get_sim07_enrol_time_int <- function(N = 2500, lambda = 1.52,
                           rho = function(t) pmin(t/360, 1)){
  
  c(0, poisson::nhpp.event.times(lambda, N - 1, rho))
}




