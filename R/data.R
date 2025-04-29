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
g_pr_c_pref <- rbind(
  dair = c(0.2, 0.2, 0.6),
  rev = c(0, 0.25, 0.75)
)



# Probabilities of strata membership for data generation
get_pop_spec_int <- function(){
  
  # proportion allocated to each silo
  r_silo <- c("early" = 0.3, "late" = 0.5, "chronic" = 0.2)
  # proportion with each site of infection by silo
  r_joint <- array(
    c(0.4,0.7,0.5,0.6,0.3,0.5), dim = c(3, 2),
    dimnames = list(c("early", "late", "chronic"), c("knee", "hip")))
  
  # randomisation probabilities for domain a by silo
  r_a <- list(
    early = NA,
    late = c("dair" = 0.5, "rev" = 0.5),
    chronic = NA
  )
  
  # preference proportions by silo
  r_a_q <- list(
    early = c("dair" = 0.9, "one" = 0.1, "two" = 0.0),
    late = c("dair" = 0.2, "one" = 0.24, "two" = 0.56),
    chronic = c("dair" = 0.2, "one" = 0.2, "two" = 0.6)
  )
  
  # randomisation within the antibiotic duration domain
  r_b <- list(
    dair = NA,
    one = c("long" = 0.5, "short" = 0.5),
    two = c("long" = 0.5, "short" = 0.5)
  )
  
  # randomisation within the antibiotic duration domain
  r_c <- list(
    dair = NA,
    one = c("long" = 0.5, "short" = 0.5),
    two = c("long" = 0.5, "short" = 0.5)
  )
  
  # randomisation within the rifampacin domain
  r_d <- c("norif" = 0.5, "rif" = 0.5)
  
  pop_spec <- list(
    r_silo = r_silo,
    r_joint = r_joint,
    r_a = r_a,
    r_a_q = r_a_q,
    r_b = r_b,
    r_c = r_c
  )
  
  pop_spec
}



# Really just defines the model parameter structure, the values will generally
# be overwritten in code.
get_sim_spec_int <- function(){
  
  
  a0 <- qlogis(0.65)
  
  # trying to simplify life:
  # just incorporate a shift for silo, don't worry about joint site
  
  # shift for late, knee
  # m1 <- qlogis(0.55) - a0
  # shift for chronic, knee
  # m2 <- qlogis(0.60) - a0
  # shift for hip
  # m3 <- qlogis(0.75) - a0
  # additional shift for late, hip
  # m4 <- qlogis(0.60) - a0 - m1 - m3
  # additional shift for chronic, hip
  # m5 <- qlogis(0.65) - a0 - m2 - m3
  # sanity check
  # X <- model.matrix(~l*j, data = CJ(j = factor(1:2), l = factor(1:3))[, .(l,j)])
  
  # a_l_j <- qlogis(
  #   array(
  #     c(plogis(X %*% c(a0, m1, m2, m3, m4, m5))),
  #     dim = c(3,2),
  #     dimnames = list(
  #       c("early", "late", "chronic"),
  #       c("knee", "hip")
  #     )
  #   )
  # )
  # p_l_j <- plogis(a_l_j)
  
  # the estimates for late and chronic are just weighted combination
  # where the weights correspond to the proportion of joints for each silo
  
  # shift for late
  m1 <- qlogis(0.57) - a0
  # shift for chronic
  m2 <- qlogis(0.63) - a0
  
  b <- c(
    "erx" = -0.1, "erx-r1" = -0.05, "erx-r2" = 0.05,
    "r1" = 0.2, "r2" = 0.09,
    # "edx" = -0.07,
    "r1d" = 0.3, "r2d" = 0.1,
    "efx" = 0.25, "f" = 0.15
  )
  
  sim_spec <- list(
    # a_l_j = a_l_j,
    # p_l_j = p_l_j,
    a0 = a0, m = c("l1" = m1, "l2" = m2),
    b = b
  )
  
  sim_spec
  
}

# Fields in data output:
# id: pt id
# l: strata 0 early, 1 late, 2 chronic
# l1: strata indicator for late
# l2: strata indicator for chronic
# j: joint indicator knee/hip
# er: revealed indicator surgery
# ed: revealed indicator duration
# ef: revealed indicator choice
# erx: 1 - revealed = non-revealed indicator for surgery
# edx: 1 - revealed = non-revealed indicator for duration
# efx: 1 - revealed = non-revealed indicator for choice
# r: randomisation for surgery domain, dair vs rev hard-coded restriction to late
# sr: preferred surgery at elicited at baseline (0 dair, 1 one-stage, 2 two-stage)
# sra: indicator derived from sr for preference for two-stage
# ra: allocated surgical approach accounting for whether rand in surg or not
# ic: indicator for treatment switch (what was planned was not received).
# rp: indicator of performed surgical approach (0 dair, 1 rev)
# srp: performed surgical approach (0 dair, 1 one-stage, 2 two-stage)
# srp1: indicator for one-stage performed
# srp2: indicator for two-stage performed
# Duration domain depends on what surgery was received, NOT what was planned.
# Because of the questions of interest:
# For one-stage long (0), short (1)
# For two-stage short (0), long (1)
# d: indicator for short/long duration
# f: indicator for choice (0 no-rif, 1 rif)
# t0: enrolment time
# eta_y: log-odds treatment success
# p_y: pr treatment success
# y: observed outcome (0 fail, 1 success)

# N = 1e6
# idx_s = 1
# set.seed(1)
# pop_spec <- get_pop_spec()
get_design_int <- function(N = 2500, pop_spec = NULL, idx_s = 1){
  
  stopifnot(!is.null(pop_spec))
  
  # silo (l)
  l <- sample(0:2, N, replace = T, prob = pop_spec$r_silo )
  # convenience indicators for late and chronic
  l1 = as.numeric(l == 1)
  l2 = as.numeric(l == 2)
  
  # joint is matrix, rows are silo, cols are knee and hip. 
  # j = 1 indicates hip
  j <- rbinom(N, 1, pop_spec$r_joint[l+1, 2])
  # reveal for late only with negligible numbers not entering in for randomised
  # treatment
  
  # if rand has been shut off then do not enrol any more into surgical domain
  if(all(is.na(unlist(pop_spec$r_a)))){
    # no surgical allocation revealed.
    er <- rep(0, N)
  } else{
    # only late silo can enter surgical domain
    er <- as.numeric(l == 1)
    # but say 2% are never revealed in this silo (for whatever reason)
    i_rec <- as.logical(rbinom(er[l==1], 1, 0.02))
    er[l==1][i_rec] <- 0
  }
  
  # randomised assignment for surgery - with silo specific allocation probs.
  r <- rep(NA, N)
  for(i in 1:length(pop_spec$r_a)){
    # r_a has allocation probs for early, late and chronic silo
    # pick the current silo
    z <- pop_spec$r_a[[i]]
    # if the elements within the silo are all non-na then randomise to the
    # set of interventions 
    if(all(!is.na(z))){
      # allocate the surgical trt based on the the number of trt for the strata
      r[l==(i-1)] <- sample(0:(length(z)-1), sum(l==(i-1)), TRUE, z)
    } else {
      # otherwise just set to zero. er will be set to zero so these will be
      # ignored anyway
      r[l==(i-1)] <- 0
    }
  }
  # fix those that are in the late silo but were not revealed, see above
  r[er == 0] <- 0
  
  
  # Just take r_a_q to indicate the surgery type that actually happens
  # as preference is not really required any more in the sims.
  
  # For each silo r_a_q indicates the prob of each surg type (dair, one, two)
  # We still will record preference for the trial.
  dtmp <- data.table(cbind(l, er, r, srp = rep(NA, N)))
  # For each silo.
  for(i in 1:length(pop_spec$r_a_q)){
    z <- pop_spec$r_a_q[[i]]
    
    # non-randomised trt just gets whatever r_a_q implies
    dtmp[l == i - 1 & er == 0, srp := sample(0:(length(z)-1), .N, TRUE, z)]
    # randomised trt to dair gets dair
    dtmp[l == i - 1 & er == 1 & r == 0, srp := 0]
    # randomised trt to rev gets one/two based on normalised probs
    dtmp[l == i - 1 & er == 1 & r == 1,
         srp := sample(1:(length(z)-1), .N, TRUE, z[2:length(z)]/sum(z[2:length(z)]))]
    
  }
  
  # for kicks, set 2% of the late silo to not revealed for randomised surgery
  # and set them to receive dair
  # ic <- rbinom(dtmp[l == 1, .N], 1, 0.02)
  # if(all(ic == 0)){ ic[1] <- 1}
  # dtmp[l == 1, ic := ic]
  # dtmp[l != 1, ic := 0]
  # dtmp[ic == 1, `:=`(er = 0, r = 0, srp = 0)]
  
  srp <- dtmp$srp
  
  # rp is now redundant.....?
  
  # was the procedure that was performed dair or revision?
  # rp <- as.numeric(srp %in% 1:2)
  
  # reveal for duration
  
  # Reveal only occurs for the period over which the effects are being evaluated.
  # Once the quest has been answered we stop revealling the randomisation.
  # Presumably we also stop randomising.
  dtmp <- data.table(srp, ed = rep(0, N))
  if(all(!is.na(pop_spec$r_b$one))){
    dtmp[srp == 1, ed := 1]
  }
  if(all(!is.na(pop_spec$r_b$two))){
    dtmp[srp == 2, ed := 1]
  }
  
  ed <- dtmp$ed
  
  # rand to long (0), short (1) for one-stage
  # rand to short (0), long (1) for two-stage
  
  dtmp <- data.table(cbind(srp, ed, d = rep(NA, N)))
  # d <- rep(NA, N)
  # For each surgery type
  for(i in 1:length(pop_spec$r_b)){
    z <- pop_spec$r_b[[i]]
    if(all(!is.na(z))){
      # d[srp==(i-1)] <- sample(0:(length(z)-1), sum(srp==(i-1)), TRUE, z)
      dtmp[srp == (i - 1), d := sample(0:(length(z)-1), .N, TRUE, z)]
    } else {
      dtmp[srp == (i - 1), d := 0]
    }
  }
  
  d <- dtmp$d
  
  # reveal for ab choice
  
  # Again, reveal only occurs for the period that effects are being
  # evaluated. Once the question is answered we reveal no more pts to rand trt.
  
  if(all(is.na(pop_spec$r_c))){
    ef <- rep(0, N)
    f <- rep(0, N)
  } else {
    # 60% reveal ab choice
    ef <- rbinom(N, 1, 0.6)
    f <- as.numeric((ef == 1) * rbinom(N, 1, pop_spec$r_c))
  }
  
  D <- data.table(
    l, l1, l2, j,
    er, ed, ef,
    erx = 1-er, edx = 1-ed, efx=1-ef,
    r,
    # sr, sra, ra,
    # ic,
    # rp,
    srp,
    srp0 = as.numeric(srp == 0),
    srp1 = as.numeric(srp == 1),
    srp2 = as.numeric(srp == 2),
    d,
    f
  )
  
  D[, id := idx_s:(N+idx_s - 1)]
  setcolorder(D, "id")
  
  D
}


get_mvn_test_data <- function(N = 30, n_time = 2, 
                    b_trt = 2, rho = 0.5, s_id = c(1, 1), s_e = 1){
  total_obs <- N * n_time
  
  cov_mat <- matrix(c(s_id[1]^2, 
                      s_id[1] * s_id[2] * rho, 
                      s_id[1] * s_id[2] * rho, 
                      s_id[2]^2), nrow = 2)
  
  d <- data.table(
    rbind(
      rmvnorm(N, mean = c(0, 0), sigma = cov_mat),
      rmvnorm(N, mean = c(0, b_trt), sigma = cov_mat)
    )
  )
  d[, id := 1:(2*N)]
  d[, trt := rep(c(0, 1), each = N)]
  d <- melt(d, measure.vars = c("V1", "V2"))
  d[variable == "V1", time := 0]
  d[variable == "V2", time := 1]
  d[, y := value + s_e]
  d[, variable := NULL]
  setkey(d, id, time)
  d
}

# Generates trial data for cohort size specified by N along with pop_spec and
# sim_spec (if provided otherwise defaults used) according to linear predictor
# specified in g.
# Fields in data output:
# id: pt id
# l: strata 0 early, 1 late, 2 chronic
# l1: strata indicator for late
# l2: strata indicator for chronic
# j: joint indicator knee/hip
# er: revealed indicator surgery
# ed: revealed indicator duration
# ef: revealed indicator choice
# erx: 1 - revealed = non-revealed indicator for surgery
# edx: 1 - revealed = non-revealed indicator for duration
# efx: 1 - revealed = non-revealed indicator for choice
# r: randomisation for surgery domain, dair vs rev hard-coded restriction to late
# sr: preferred surgery at elicited at baseline (0 dair, 1 one-stage, 2 two-stage)
# sra: indicator derived from sr for preference for two-stage
# ra: allocated surgical approach accounting for whether rand in surg or not
# ic: indicator for treatment switch (what was planned was not received).
# rp: indicator of performed surgical approach (0 dair, 1 rev)
# srp: performed surgical approach (0 dair, 1 one-stage, 2 two-stage)
# srp1: indicator for one-stage performed
# srp2: indicator for two-stage performed
# Duration domain depends on what surgery was received, NOT what was planned.
# Because of the questions of interest:
# For one-stage long (0), short (1)
# For two-stage short (0), long (1)
# d: indicator for short/long duration
# f: indicator for choice (0 no-rif, 1 rif)
# t0: enrolment time
# eta_y: log-odds treatment success
# p_y: pr treatment success
# y: observed outcome (0 fail, 1 success)
get_trial_data_int <- function(
    N = 1000,
    
    # reference level log odds of response
    mu = 0,
    # silo effects
    # silo 1 early,  silo 2 LATE ACUTE, silo 3 chronic
    b_silo = c(0, -0.2, -0.3),
    b_jnt = c(0, 0.4),
    b_pref = c(0, -0.2),
    # dair, one, two-stage
    b_d1 = c(0, -0.1, -0.1),
    # b_d1_ex = c(0, -0.1, -0.1, -0.1),
    # is short duration (6 wk) is non-inferior to long duration (12 wk) ab?
    # misc, 12wk, 6wk
    b_d2 = c(0, -0.4, 0.6),
    # is a 12 wk duration is superior to no extended prophylaxis.
    # misc, none, 12wk
    b_d3 = c(0, 0.1, 0.3),
    # is rifampicin is superior to no rifampicin.
    # misc, none, rif
    b_d4 = c(0, -0.1, -0.2),
    
    dec_sup = list(
      surg = NA,
      ext_proph = NA,
      ab_choice = NA
    ),
    dec_ni = list(
      ab_dur = NA
    ),
    # dec_eq = list(
    #   ab_dur = NA
    # ),
    dec_sup_fut = list(
      surg = NA,
      ext_proph = NA,
      ab_choice = NA
    ),
    dec_ni_fut = list(
      ab_dur = NA
    ),
    # dec_inf = list(
    #   ab_dur = NA
    # ),
    
    idx_s = NULL,
    t0 = NULL,
    id_analys = NULL
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
  # location of infection knee/hip
  d[silo == 1, jnt := sample(1:K_j, size = .N, replace = T, prob = g_pr_jnt[1, ])]
  d[silo == 2, jnt := sample(1:K_j, size = .N, replace = T, prob = g_pr_jnt[2, ])]
  d[silo == 3, jnt := sample(1:K_j, size = .N, replace = T, prob = g_pr_jnt[3, ])]
  
  # Domain treatment effects
  
  ## Surgical domain--------
  
  # Early disease receive non-randomised intervention in that they are not 
  # part of the experimental considerations so the allocation here is really 
  # just arbitrary based on preferences for this cohort.
  
  # clinician decides on what surgery type occurs for silo 1 and 2
  # pref_rev refers to the preferred revision type, either 
  # one-stage or two-stage.
  d[silo == 1, pref_rev := sample(1:2, size = .N, replace = T, prob = g_pr_e_pref[2, 2:3])]
  d[silo == 1, d1 := sample(c(-98, -99), size = .N, replace = T, prob = g_pr_e_surg)]
  # dair - replace the -98 with 1
  d[silo == 1 & d1 == -98, d1 := 1]
  # revision - split the -99 into 2 or 3 corresponding to one-stage or two-stage
  # the probability of one-stage/two-stage is very highly weighted towards the
  # preference, but sometimes the preferred revision type will not occur
  d[silo == 1 & d1 == -99 & pref_rev == 1, d1 := sample(2:3, .N, replace = T, prob = c(0.95, 0.05))]
  d[silo == 1 & d1 == -99 & pref_rev == 2, d1 := sample(2:3, .N, replace = T, prob = c(0.05, 0.95))]
  
  # don't think this kind of approach (trying to adjust for non-randomised)
  # would work as the silos are entirely collinear
  # d[silo == 1 & d1 == 1, d1_ex := 2]
  # d[silo == 1 & d1 == 2, d1_ex := 3]
  # d[silo == 1 & d1 == 3, d1_ex := 4]
  
  
  # Here we are addressing the late silo data who are the only cohort to enter
  # into the surgical domain randomised treatment.
  
  # If a decision is yet to be made on superiority or futility for superiority
  # in the surgical domain then we randomise the surgical domain.
  # However, if revision has been deemed superior then everyone gets a one-stage
  # or two-stage revision dependent on the clinician preference for the unit.
  # Alternatively, if superiority assessment has been deemed futile then all units
  # get DAIR.
  if(is.na(dec_sup$surg) & is.na(dec_sup_fut$surg)){
    # Default situation, we are allocating randomised trt
    
    # For late acute, we have a randomised assignment process for surgical domain
    # rather than selection. For the purposes here, the simulation approach is 
    # similar to that above but the initial dair vs rev assignment is 1:1. 
    # All the silo are assume to be revealed to randomised trt. In practice this
    # probably won't be the case and so we would need some indicator of the reveal
    # to select the appropriate subset of data for the analysis.
    d[silo == 2, pref_rev := sample(1:2, size = .N, replace = T, prob = g_pr_l_pref[2, 2:3])]
    # 1:1 rand assignment to dair (3) or rev (4) 
    d[silo == 2, d1 := sample(c(-98, -99), size = .N, replace = T)]
    # selection surgical type dair (1), one(2), two(3)
    # dair is fixed, i.e. preference does not come into play
    d[silo == 2 & d1 == -98, d1 := 1]
    # under revision, d1 (where d1 indicates one or two-stage) preference will
    # direct which is received so if the preference is for one stage then with
    # high probability d1 will be one-stage (similarly for two-stage)
    d[silo == 2 & d1 == -99 & pref_rev == 1, d1 := sample(2:3, .N, replace = T, prob = c(0.95, 0.05))]
    d[silo == 2 & d1 == -99 & pref_rev == 2, d1 := sample(2:3, .N, replace = T, prob = c(0.05, 0.95))]
    
  } else if (!is.na(dec_sup$surg)) {
    
    # Revision has been deemed superior - allocation is now either one-stage or
    # two-stage dependent on preference.
    
    d[silo == 2, pref_rev := sample(1:2, size = .N, replace = T, prob = g_pr_l_pref[2, 2:3])]
    
    # under revision, d1 has high probability of being the preference
    d[silo == 2 & pref_rev == 1, d1 := sample(2:3, .N, replace = T, prob = c(0.95, 0.05))]
    d[silo == 2 & pref_rev == 2, d1 := sample(2:3, .N, replace = T, prob = c(0.05, 0.95))]
    
  } else if (!is.na(dec_sup_fut$surg)) {
    
    d[silo == 2, pref_rev := sample(1:2, size = .N, replace = T, prob = g_pr_l_pref[2, 2:3])]
    
    # revision is futile, everyone now gets dair irrespective of preference
    d[silo == 2 , d1 := 1]
    
  }
  
  
  
  # exactly same approach as early silo for the chronic silo but using different probabilities
  d[silo == 3, pref_rev := sample(1:2, size = .N, replace = T, prob = g_pr_c_pref[2, 2:3])]
  d[silo == 3, d1 := sample(c(-98, -99), size = .N, replace = T, prob = g_pr_c_surg)]
  d[silo == 3 & d1 == -98, d1 := 1]
  # under revision, d1 has high probability of being the preference
  d[silo == 3 & d1 == -99 & pref_rev == 1, d1 := sample(2:3, .N, replace = T, prob = c(0.95, 0.05))]
  d[silo == 3 & d1 == -99 & pref_rev == 2, d1 := sample(2:3, .N, replace = T, prob = c(0.05, 0.95))]
  
  
  ## Antibiotic duration domain--------
  
  # Antibiotic backbone duration - only applicable to units that recv 
  # one-stage revision.
  # All others just get set to the reference index.
  
  # If the surgical intervention was dair or two-stage, set d2 to 1, which is a 
  # soc (whatever that happens to be)
  d[d1 %in% c(1, 3), d2 := 1]
  
  # We have to have a small set of pt with d1 = 2 and d2 = 1, i.e. just receiving 
  # the hospital/surgeon select duration. If we do not make this
  # assumption we will hit an identification problem and that would mean that we 
  # could not recover the true parameters without fairly major bias.
  if(is.na(dec_ni$ab_dur) & is.na(dec_ni_fut$ab_dur) ){
  # if(is.na(dec_ni$ab_dur) & is.na(dec_ni_fut$ab_dur) &
  #    is.na(dec_eq$ab_dur) & is.na(dec_inf$ab_dur)){
    # Default situation, we are allocating randomised trt
    
    # We assume that of those receiving one-stage revision,
    # 70% enter into the duration domain (the remainder being assigned to the ref index)
    d[d1 %in% c(2), d2 := sample(1:3, size = .N, replace = T, prob = c(0.3, 0.35, 0.35))]
  
  } else if (!is.na(dec_ni$ab_dur) ) {    
  # } else if (!is.na(dec_ni$ab_dur) | !is.na(dec_eq$ab_dur)) {
    # 6wks is NI to 12wks (or equivalent) so we assume that those receiving one-stage revision
    # all receive the 6wk trt
    d[d1 %in% c(2), d2 := 3]
    
  }  else if (!is.na(dec_ni_fut$ab_dur) ) {
  # }  else if (!is.na(dec_ni_fut$ab_dur) | !is.na(dec_inf$ab_dur)) {
    
    # 6wk ni assessment is futile or inferior, everyone now gets 12wk
    d[d1 %in% c(2), d2 := 2]
  } 
  
  ## Extended prophylaxis ----
  
  # Extended prophylaxis - only applicable to units recv two-stage, the remainder
  # get the soc (whatever that happens to be), i.e. d3 set to 1
  d[d1 %in% c(1, 2), d3 := 1]
  
  # Again, we have to have some small group with d1 = 3 and d2 = 1 otherwise we have
  # an identification problem that would mean that we cannot recover the true
  # parameters and our results will be biased.
  
  if(is.na(dec_sup$ext_proph) & is.na(dec_sup_fut$ext_proph)){
    # Default situation, we are allocating randomised trt
    
    # Here we assume that of those receiving two-stage revision,
    # 90% enter into the Extended prophylaxis domain
    d[d1 %in% c(3), d3 := sample(1:3, size = .N, replace = T, prob = c(0.1, 0.45, 0.45))]
    
  } else if (!is.na(dec_sup$ext_proph)){
    
    # 12wks is superior to none so we assume that those receiving two-stage revision
    # all receive the 12wk trt
    d[d1 %in% c(3), d3 := 3]
  } else if (!is.na(dec_sup_fut$ext_proph)){
    
    # 12wks is futile, everyone now gets none
    d[d1 %in% c(3), d3 := 2]
  }
  
  ## Antibiotic choice ----
  
  # Choice domain is independent to others but only 60% of the population
  # enter it.
  
  d[, g4 := sample(1:2, size = .N, replace = T, prob = c(0.4, 0.6))]
  d[g4 == 1, d4 := 1]
  
  if(is.na(dec_sup$ab_choice) & is.na(dec_sup_fut$ab_choice)){
    
    # Default situation
    d[g4 == 2, d4 := sample(2:3, size = .N, replace = T)]
    d[, g4 := NULL]
    
  } else if (!is.na(dec_sup$ab_choice)){
    
    # Superiority decision, all get allocated to rif
    d[g4 == 2, d4 := 3]
    d[, g4 := NULL]
  } else if (!is.na(dec_sup_fut$ab_choice)){
    
    # rif is futile, everyone now gets none
    d[g4 == 2, d4 := 2]
    d[, g4 := NULL]
  }
  
  
  
  # Outcome
  d[, mu := mu]
  d[, b_silo := b_silo[silo]]
  d[, b_jnt := b_jnt[jnt]]
  d[, b_pref := b_pref[pref_rev]]
  d[, b_d1 := b_d1[d1]]
  d[, b_d2 := b_d2[d2]]
  d[, b_d3 := b_d3[d3]]
  d[, b_d4 := b_d4[d4]]
  
  # Think about preference as something that is induced in the clinician by
  # the status of the patient, i.e. it is a kind of baseline adjustment.
  d[, eta := mu + b_silo + b_jnt + b_pref + b_d1 + b_d2 + b_d3 + b_d4      ]
  
  d[, y := rbinom(.N, 1, plogis(eta))]
  
  
  
  # add a pt id
  if(!is.null(idx_s)){
    d[, id := idx_s:(N+idx_s - 1)]
  } 
  
  if(!is.null(t0)){
    d[, t0 := t0]
  } else{
    d[, t0 := rep(0.0, .N)]
  }
  
  if(!is.null(id_analys)){
    d[, id_analys := id_analys]
  } else{
    d[, id_analys := NA]
  }
  
  setcolorder(d, c("id", "t0", "id_analys"))
  
  
  list(
    d = d, 
    
    mu = mu, 
    b_silo = b_silo, b_jnt = b_jnt, b_pref = b_pref,  
    b_d1 = b_d1, b_d2 = b_d2, b_d3 = b_d3, b_d4 = b_d4, 
    
    K_s = K_s, K_j = K_j, K_p  = K_p,
    K_d1 = K_d1, K_d2 = K_d2 , K_d3 = K_d3, K_d4 = K_d4
  )
  
  
}


# is rev superior to dair
get_stan_data_surgery_int <- function(d_all){
  
  # Only applies to late silo and presently assuming that all pt are revealed.
  # In practice only a subset will be revealed so we would need some idea
  # of the proportion and indicators of such.
  
  # Here we are solely interested in conditioning on intervention and 
  # joint effects. Again, in practice we would need to condition on other 
  # elements such as site and time of analysis (trend adjustment) and 
  # baseline characteristics.
  
  d_mod <- d_all[, .(y = sum(y), n = .N),
                 keyby = .(silo, jnt, pref_rev, d1)]

  d_mod_d1 <- d_mod[silo == 2]
  
  
  ld <- list(
    
    N = nrow(d_mod_d1), 
    
    y = d_mod_d1[, y], 
    n = d_mod_d1[, n], 
    
    K_jnt = 2, K_pref = 2, K_d1 = 3, 
    
    jnt = d_mod_d1[, jnt],
    pref = d_mod_d1[, pref_rev],
    d1 = d_mod_d1[, d1],
    
    # sample size of those where preference is for one-stage
    N_p1 = d_mod_d1[pref_rev == 1, .N],
    # sample size of those where preference is for two-stage
    N_p2 = d_mod_d1[pref_rev == 2, .N],
    
    # indexes for those with preference for one-stage
    ix_p1 = d_mod_d1[pref_rev == 1, which = T],
    # indexes for those with preference for two-stage
    ix_p2 = d_mod_d1[pref_rev == 2, which = T],
    # number of trials within each of these subsets
    n_p1 = d_mod_d1[pref_rev == 1, n],
    n_p2 = d_mod_d1[pref_rev == 2, n],
    
    prior_only = 0
  )
  
  list(
    d_s = d_mod_d1,
    ld = ld
  )
  
}

# is short duration (6 wk) non-inferior to long duration (12 wk) antibiotic?
get_stan_data_ab_dur_int <- function(d_all){
  
  # only applicable to one-stage revision noting that this may occur
  # in all silos
  d_s <- d_all[d1 == 2 & d2 %in% 2:3, ]
  
  # subset based on ed 
  # - indicates revealed to randomised antibiotic duration intervention
  # 12wk
  d_s[d2 == 2, x := 0]
  # 6wk
  d_s[d2 == 3, x := 1]
  
  # all silos contribute so need to account for that in model
  d_s <- d_s[, .(y = sum(y), n = .N),
             keyby = .(silo, x)]
  
  d_s[silo == 1, `:=`(l2 = 0, l3 = 0)]
  d_s[silo == 2, `:=`(l2 = 1, l3 = 0)]
  d_s[silo == 3, `:=`(l2 = 0, l3 = 1)]
  
  ld <- list(
    N = nrow(d_s), y = d_s$y, n = d_s$n,
    x = d_s$x, 
    
    l2 = d_s$l2, 
    l3 = d_s$l3, 
    pri_b_sd = rep(1, 1),
    pri_m_sd = rep(1, 2),
    prior_only = 0
  )
  
  list(
    d_s = d_s,
    ld = ld
  )
  
}

# is 12 wk duration superior to no extended prophylaxis.
get_stan_data_ext_proph_int <- function(d_all){
  
  # only applicable to two-stage revision noting that this may occur
  # in all silos
  d_s <- d_all[d1 == 3 & d3 %in% 2:3, ]
  
  # subset based on ed 
  # - indicates revealed to randomised antibiotic duration intervention
  # none
  d_s[d3 == 2, x := 0]
  # 12wk
  d_s[d3 == 3, x := 1]
  
  # all silos contribute so need to account for that in model
  d_s <- d_s[, .(y = sum(y), n = .N),
             keyby = .(silo, x)]
  
  d_s[silo == 1, `:=`(l2 = 0, l3 = 0)]
  d_s[silo == 2, `:=`(l2 = 1, l3 = 0)]
  d_s[silo == 3, `:=`(l2 = 0, l3 = 1)]
  
  ld <- list(
    N = nrow(d_s), y = d_s$y, n = d_s$n,
    x = d_s$x, 
    
    l2 = d_s$l2, 
    l3 = d_s$l3, 
    pri_b_sd = rep(1, 1),
    pri_m_sd = rep(1, 2),
    prior_only = 0
  )
  
  list(
    d_s = d_s,
    ld = ld
  )
  
}



# is rifampicin superior to no rifampicin.
get_stan_data_ab_choice_int <- function(d_all){
  
  # only applicable to two-stage revision
  d_s <- d_all[d4 %in% 2:3, ]
  
  # subset based on ed 
  # - indicates revealed to randomised antibiotic duration intervention
  # none
  d_s[d4 == 2, x := 0]
  # rif
  d_s[d4 == 3, x := 1]
  
  # all silos contribute so need to account for that in model
  # d_s <- d_s[, .(y = sum(y), n = .N),
  #            keyby = .(silo, jnt, x)]
  
  d_s <- d_s[, .(y = sum(y), n = .N),
             keyby = .(silo, x)]
  
  d_s[silo == 1, `:=`(l2 = 0, l3 = 0)]
  d_s[silo == 2, `:=`(l2 = 1, l3 = 0)]
  d_s[silo == 3, `:=`(l2 = 0, l3 = 1)]
  
  ld <- list(
    N = nrow(d_s), y = d_s$y, n = d_s$n,
    x = d_s$x, 
    
    l2 = d_s$l2, 
    l3 = d_s$l3, 
    pri_b_sd = rep(1, 1),
    pri_m_sd = rep(1, 2),
    prior_only = 0
  )
  
  list(
    d_s = d_s,
    ld = ld
  )
  
}

get_stan_data_all_int <- function(d_all){
  
  # convert from binary representation to binomial (successes/trials)
  d_mod <- d_all[, .(y = sum(y), n = .N, eta = round(unique(eta), 3)), 
                keyby = .(silo, jnt, pref_rev, d1, d2, d3, d4)]
  d_mod[, eta_obs := round(qlogis(y / n), 3)]
  d_mod[, p_obs := y / n]
  
  # really don't need to do this for all of the domains but it provides
  # a uniform approach for generating the parameters of interest
  
  # Surgical -----
  
  # restrict to silo 2 (late acute silo) for gcomp for d1 randomised comparisons
  # here I just use the silo assignment to imply the right subset of units
  # but in practice there may need to be an indicator variable to for this.
  d_mod_d1 <- d_mod[silo == 2]
  
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
    silo = d_mod[, silo], 
    jnt = d_mod[, jnt],
    pref = d_mod[, pref_rev],
    d1 = d_mod[, d1],
    d2 = d_mod[, d2],
    d3 = d_mod[, d3],
    d4 = d_mod[, d4],
    
    # Number of levels for silos, joints, pref and each trt.
    K_silo = length(unique(d_all$silo)), 
    K_jnt = length(unique(d_all$jnt)), 
    K_pref = length(unique(d_all$pref_rev)), 
    K_d1 = length(unique(d_all$d1)), 
    K_d2 = length(unique(d_all$d2)), 
    K_d3 = length(unique(d_all$d3)), 
    K_d4 = length(unique(d_all$d4)), 
    
    
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
    
    prop_p1 = d_all[, .(wgt = .N/nrow(d_all)), keyby = pref_rev][pref_rev == 1, wgt],
    prop_p2 = d_all[, .(wgt = .N/nrow(d_all)), keyby = pref_rev][pref_rev == 2, wgt],
    
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
  
  list(
    d_mod = d_mod,
    d_mod_d1 = d_mod_d1,
    d_mod_d2 = d_mod_d2,
    d_mod_d3 = d_mod_d3,
    d_mod_d4 = d_mod_d4,
    
    ld = ld
  )
  
}

get_empirical_risk <- function(
    nsim = 1000,
    N = 1000,
    mu = 0,
    b_silo = c(0, -0.3, -0.2),
    b_jnt = c(0, 0.4),
    b_pref = c(0, -0.2),
    # dair, one, two-stage
    b_d1 = c(0, -0.1, -0.2),
    b_d2 = c(0, -0.4, 0.6),
    b_d3 = c(0, 0.1, 0.3),
    b_d4 = c(0, -0.1, -0.2)
){
  
  # m1 <- cmdstanr::cmdstan_model("stan/model11.stan")
  # mu <- 0
  # # silo effects
  # b_silo <- c(0, -0.3, -0.2)
  # b_jnt <- c(0, 0.4)
  # b_pref <- c(0, -0.2)
  # # dair, one, two-stage
  # b_d1 <- c(0, -0.6, -0.6)
  # b_d2 <- c(0, -0.4, 0.6)
  # b_d3 <- c(0, 0.1, 0.3)
  # b_d4 <- c(0, -0.1, -0.2)
  
  nsim <- 5000
  N = 2.5e3
  i <- 1
  
  tic()
  
  d_simres <- rbindlist(parallel::mclapply(
    # X=1:g_cfgsc$nsim, mc.cores = g_cfgsc$mc_cores, FUN=function(ix) {
    X=1:nsim, mc.cores = g_cfgsc$mc_cores, FUN=function(ix) {
      
      log_info("Risk estimate sim ", ix);
      
      l_new <- get_trial_data_int(
        N = N, 
        
        # reference level log odds of response
        mu = mu,
        # silo effects
        b_silo = b_silo,
        b_jnt = b_jnt, # knee/hip
        b_pref = b_pref,
        # dair, one, two-stage
        b_d1 = b_d1,
        b_d2 = b_d2,
        b_d3 = b_d3,
        b_d4 = b_d4,
        
        # enrolment to all groups
        dec_sup = list(
          surg = NA,
          ext_proph = NA,
          ab_choice = NA
        ),
        dec_ni = list(
          ab_dur = NA
        ),
        dec_sup_fut = list(
          surg = NA,
          ext_proph = NA,
          ab_choice = NA
        ),
        dec_ni_fut = list(
          ab_dur = NA
        ),
        
        # utility items needed by default data gen
        idx_s = 1,
        t0 = 1:N,
        id_analys = 1
      )
      
      d_all <- copy(l_new$d)
      
      lsd <- get_stan_data_all_int(d_all)
      
      # fit model - does it matter that I continue to fit the model after the
      # decision is made...?
      # snk <- capture.output(
      # f_1 <- m1$sample(
      #   lsd$ld, iter_warmup = 1000, iter_sampling = 2000,
      #   parallel_chains = 1, chains = 1, refresh = 0, show_exceptions = T,
      #   max_treedepth = 11)
      # )
      
      # snk <- capture.output(
      f_1 <- m1$pathfinder(lsd$ld, num_paths=20, single_path_draws=200,
                          history_size=50, max_lbfgs_iters=100,
                          refresh = 0, draws = 2000)
      # )
      
      d_post <- data.table(f_1$draws(
        variables = c(
          "prf_1", "prf_2", # observed proportion of preferring each type
          "nu_d1_1", "nu_d1_23", "lor_d1", # log odds and lor
          "p_d1_1", "p_d1_23", "rd_d1",
          
          "nu_d2_2", "nu_d2_3", "lor_d2",
          "p_d2_2", "p_d2_3", "rd_d2",
          
          "nu_d3_2", "nu_d3_3", "lor_d3",
          "p_d3_2", "p_d3_3", "rd_d3",
          
          "nu_d4_2", "nu_d4_3", "lor_d4",
          "p_d4_2", "p_d4_3", "rd_d4"
          
          ),   # risk scale
        format = "matrix"))
      
      hist(d_post$rd_d1)
      
      data.table(t(colMeans(d_post)))
      
      
      
    }), idcol = "sim")
  
  toc()
  
  
  d_simres <- melt(d_simres, id.vars = "sim")
  
  d_simres_smry <- d_simres[, .(mu = mean(value)), keyby = .(variable)]
  
  d_simres_smry
  
  # Sanity check - are we at least close to the log-ors implied by the params?
  d_fig <- d_simres[variable %in% c("lor_d1", "lor_d2", "lor_d3", "lor_d4")]
  ggplot(d_fig, aes(x = value, group = variable)) +
    geom_histogram() +
    geom_vline(
      data = d_fig[, .(mu = mean(value), 
                       min = min(value),
                       max = max(value)), keyby = .(variable)],
      aes(xintercept = mu, group = variable), col = "red"
    ) +
    facet_wrap(~variable, scales = "free_x") 
  
  
  # preference based on sample estimates.
  d_fig <- d_simres[variable %in% c("prf_1", "prf_2")]
  ggplot(d_fig, aes(x = value, group = variable)) +
    geom_histogram() +
    geom_vline(
      data = d_fig[, .(mu = mean(value), 
                       min = min(value),
                       max = max(value)), keyby = .(variable)],
      aes(xintercept = mu, group = variable), col = "red"
    ) +
    facet_wrap(~variable, scales = "free_x") 
  
  # primary risk differences of interest
  d_fig <- d_simres[variable %in% c("rd_d1", "rd_d2", "rd_d3", "rd_d4")]
  ggplot(d_fig, aes(x = value, group = variable)) +
    geom_histogram() +
    geom_vline(
      data = d_fig[, .(mu = mean(value), 
                       min = min(value),
                       max = max(value)), keyby = .(variable)],
      aes(xintercept = mu, group = variable), col = "red"
    ) +
    facet_wrap(~variable, scales = "free_x") 
  
  
  # risk scale estimates 
  d_fig <- d_simres[variable %in% c("p_d1_1", "p_d1_23",
                                    "p_d2_2", "p_d2_3",
                                    "p_d3_2", "p_d3_3",
                                    "p_d4_2", "p_d4_3"
                                    )]
  ggplot(d_fig, aes(x = value, group = variable)) +
    geom_histogram() +
    geom_vline(
      data = d_fig[, .(mu = mean(value), 
                       min = min(value),
                       max = max(value)), keyby = .(variable)],
      aes(xintercept = mu, group = variable), col = "red"
    ) +
    facet_wrap(~variable, scales = "free_x") 
  
  
  # 
  
}


get_enrol_time_int <- function(N = 2500, lambda = 1.52,
                           rho = function(t) pmin(t/360, 1)){
  
  c(0, poisson::nhpp.event.times(lambda, N - 1, rho))
}

get_design_opts_int <- function(){
  
  d_x <- CJ(
    l = c("e", "l", "c"),
    er = c("rand", "nonrand"),
    r = c("dair", "rev"),
    sr = c("dair", "one", "two"),
    ed = c("rand", "nonrand"),
    d = c("short", "long"),
    ef = c("rand", "nonrand"),
    f = c("norif", "rif")
  )
  
  d_x <- d_x[!(er == "nonrand" & ed == "nonrand" & ef == "nonrand")]
  d_x <- d_x[!(l == "e" & er == "rand")]
  d_x <- d_x[!(l == "c" & er == "rand")]
  d_x <- d_x[!(l == "l" & er == "nonrand")]
  
  d_x <- d_x[!(l == "e" & sr == "two")]
  
  d_x <- d_x[!(r == "dair" & sr == "one")]
  d_x <- d_x[!(r == "dair" & sr == "two")]
  d_x <- d_x[!(r == "rev" & sr == "dair")]
  d_x <- d_x[!(r == "dair" & ed == "rand")]
  
  d_x <- d_x[!(ed == "nonrand" & r == "rev")]
  d_x <- d_x[!(ef == "nonrand" & f == "rif")]
  
  fwrite(d_x, "x.csv")
  
  
  
}




