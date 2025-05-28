library(data.table)
library(pbapply)
library(ggplot2)
library(patchwork)




risk_pars_surg <- function(
  p_d1_alloc = 0.5,
  p_d2_entry = 0.7,
  p_d2_alloc = 0.5,
  p_d3_entry = 0.9,
  p_d3_alloc = 0.5,
  p_d4_entry = 0.6,
  p_d4_alloc = 0.5,
  # 40% pref for two-stage
  p_pref = 0.4,
  # log-odds
  mu = 0.7,
  # log-or
  # different baseline risk for rev
  bp = -0.4,
  bd1 = c(0, 0, 0),
  bd2 = c(0, 0.1, 0),
  bd3 = c(0, 0, -0.2),
  bd4 = c(0, 0.1, 0.2)
  ){
  
  # dair
  # averages over combinations from f1_1:  ~ 1 + pref + d4
  # mu + bd4[nonrand]
  # mu + bd4[norif]
  # mu + bd4[rif]
  # mu + bp + bd4[nonrand]
  # mu + bp + bd4[norif]
  # mu + bp + bd4[rif]
  # if 60% prefer rev(1) (if they had recv rev) and 60% enter ab choice
  # then 60% weight is given to the first three and 40% to the last
  # in the first three 60% enter ab choice so we have 
  # 0.6 * 0.6 = 0.36 allocated to norif and rif and the rest left to nonrand
  # in the second three, using the same logic
  # 0.6 * 0.4 = 0.24 allocated to no rif and rif and the rest left to nonrand
  # so we want:
  # 0.24 * (mu + bd4[nonrand])
  # 0.18 * (mu + bd4[norif])
  # 0.18 * (mu + bd4[rif])
  # 0.16 * (mu + bp + bd4[nonrand])
  # 0.12 * (mu + bp + bd4[norif])
  # 0.12 * (mu + bp + bd4[rif])
  
  prop_tru <- numeric(6)
  prop_tru[1] <- (1-p_pref) - (p_d4_entry * (1-p_pref))
  prop_tru[2] <- p_d4_entry * (1-p_pref)/2
  prop_tru[3] <- p_d4_entry * (1-p_pref)/2
  prop_tru[4] <- (p_pref) - (p_d4_entry * (p_pref))
  prop_tru[5] <- p_d4_entry * (p_pref)/2
  prop_tru[6] <- p_d4_entry * (p_pref)/2
  
  p_dair <- 
    prop_tru[1] * plogis(mu + bd4[1]) + 
    prop_tru[2] * plogis(mu + bd4[2]) + 
    prop_tru[3] * plogis(mu + bd4[3]) + 
    prop_tru[4] * plogis(mu + bp + bd4[1]) + 
    prop_tru[5] * plogis(mu + bp + bd4[2]) + 
    prop_tru[6] * plogis(mu + bp + bd4[3])
  
  
  # rev 1
  # averages over combinations from f1_2:  ~ 1 + d2 + d4
  # mu + bd1[2] + bd2[1] + bd4[nonrand]
  # mu + bd1[2] + bd2[1] + bd4[norif]
  # mu + bd1[2] + bd2[1] + bd4[rif]    
  # mu + bd1[2] + bd2[2] + bd4[nonrand]
  # mu + bd1[2] + bd2[2] + bd4[norif]
  # mu + bd1[2] + bd2[2] + bd4[rif]
  # mu + bd1[2] + bd2[3] + bd4[nonrand]
  # mu + bd1[2] + bd2[3] + bd4[norif]
  # mu + bd1[2] + bd2[3] + bd4[rif]
  # 70% of rev(1) enter into d2, 60% enter d4
  # 0.6 * 0.3 = 0.18 allocated to norif and rif (rest 0.3 - 0.18 for nonrand)
  # 0.6 * 0.35 = 0.21 allocated to norif and rif (rest for nonrand)
  # 0.6 * 0.35 = 0.21 allocated to norif and rif (rest for nonrand)
  
  prop_tru <- numeric(9)
  prop_tru[1] <- (1-p_d2_entry) - (p_d4_entry * (1-p_d2_entry))
  prop_tru[2] <- p_d4_entry * (1-p_d2_entry) / 2
  prop_tru[3] <- p_d4_entry * (1-p_d2_entry) / 2
  prop_tru[4] <- (p_d2_entry/2) - (p_d4_entry * (p_d2_entry/2))
  prop_tru[5] <- p_d4_entry * (p_d2_entry/2) / 2
  prop_tru[6] <- p_d4_entry * (p_d2_entry/2) / 2
  prop_tru[7] <- (p_d2_entry/2) - (p_d4_entry * (p_d2_entry/2))
  prop_tru[8] <- p_d4_entry * (p_d2_entry/2) / 2
  prop_tru[9] <- p_d4_entry * (p_d2_entry/2) / 2
  
  p_rev_1 <- 
    prop_tru[1] * plogis(mu + bd1[2] + bd2[1] + bd4[1]) + 
    prop_tru[2] * plogis(mu + bd1[2] + bd2[1] + bd4[2]) + 
    prop_tru[3] * plogis(mu + bd1[2] + bd2[1] + bd4[3]) + 
    prop_tru[4] * plogis(mu + bd1[2] + bd2[2] + bd4[1]) + 
    prop_tru[5] * plogis(mu + bd1[2] + bd2[2] + bd4[2]) + 
    prop_tru[6] * plogis(mu + bd1[2] + bd2[2] + bd4[3]) + 
    prop_tru[7] * plogis(mu + bd1[2] + bd2[3] + bd4[1]) + 
    prop_tru[8] * plogis(mu + bd1[2] + bd2[3] + bd4[2]) + 
    prop_tru[9] * plogis(mu + bd1[2] + bd2[3] + bd4[3])  
  
  # rev 2
  # averages over combinations from f1_2:   ~ 1 + d3 + d4
  # pref is included as intercept = mu + bp for this cohort
  # mu + bp + bd1[3] + bd3[1] + bd4[nonrand]
  # mu + bp + bd1[3] + bd3[1] + bd4[norif]
  # mu + bp + bd1[3] + bd3[1] + bd4[rif]    
  # mu + bp + bd1[3] + bd3[2] + bd4[nonrand]
  # mu + bp + bd1[3] + bd3[2] + bd4[norif]
  # mu + bp + bd1[3] + bd3[2] + bd4[rif]
  # mu + bp + bd1[3] + bd3[3] + bd4[nonrand]
  # mu + bp + bd1[3] + bd3[3] + bd4[norif]
  # mu + bp + bd1[3] + bd3[3] + bd4[rif]
  # 90% of rev(2) enter into d3, 60% enter d4
  # 0.6 * 0.1 = 0.06 allocated to norif and rif (rest 0.1 - 0.6 for nonrand)
  # 0.6 * 0.45 = 0.27 allocated to norif and rif (rest for nonrand)
  # 0.6 * 0.45 = 0.27 allocated to norif and rif (rest for nonrand)
  
  prop_tru <- numeric(9)
  prop_tru[1] <- (1-p_d3_entry) - (p_d4_entry * (1-p_d3_entry))
  prop_tru[2] <- p_d4_entry * (1-p_d3_entry) / 2
  prop_tru[3] <- p_d4_entry * (1-p_d3_entry) / 2
  prop_tru[4] <- (p_d3_entry/2) - (p_d4_entry * (p_d3_entry/2))
  prop_tru[5] <- p_d4_entry * (p_d3_entry/2) / 2
  prop_tru[6] <- p_d4_entry * (p_d3_entry/2) / 2
  prop_tru[7] <- (p_d3_entry/2) - (p_d4_entry * (p_d3_entry/2))
  prop_tru[8] <- p_d4_entry * (p_d3_entry/2) / 2
  prop_tru[9] <- p_d4_entry * (p_d3_entry/2) / 2
  
  # contribution of non-rand, bp needs to be included here
  p_rev_2 <- 
    prop_tru[1] * plogis(mu + bp + bd1[3] + bd3[1] + bd4[1]) + 
    prop_tru[2] * plogis(mu + bp + bd1[3] + bd3[1] + bd4[2]) + 
    prop_tru[3] * plogis(mu + bp + bd1[3] + bd3[1] + bd4[3]) + 
    prop_tru[4] * plogis(mu + bp + bd1[3] + bd3[2] + bd4[1]) + 
    prop_tru[5] * plogis(mu + bp + bd1[3] + bd3[2] + bd4[2]) + 
    prop_tru[6] * plogis(mu + bp + bd1[3] + bd3[2] + bd4[3]) + 
    prop_tru[7] * plogis(mu + bp + bd1[3] + bd3[3] + bd4[1]) + 
    prop_tru[8] * plogis(mu + bp + bd1[3] + bd3[3] + bd4[2]) + 
    prop_tru[9] * plogis(mu + bp + bd1[3] + bd3[3] + bd4[3])  
  
  # rev
  prop_tru <- c(1 - p_pref, p_pref)
  p_rev <- prop_tru[1] * p_rev_1 + prop_tru[2] * p_rev_2
  
  
  c(rd = p_rev - p_dair,
    p_dair = p_dair, p_rev = p_rev, 
    p_rev_1 = p_rev_1, p_rev_2 = p_rev_2)
  
}


risk_pars_dur <- function(
    p_d2_entry = 0.7,
    p_d2_alloc = 0.5,
    p_d4_entry = 0.6,
    p_d4_alloc = 0.5,
    # log-odds
    mu = 0.7,
    # log-or
    # different baseline risk for rev
    bp = -0.4,
    bd1 = c(0, 0, 0),
    bd2 = c(0, 0.1, 0),
    bd3 = c(0, 0, -0.2),
    bd4 = c(0, 0.1, 0.2)
){
  
  # only model that is relevant is f1_2 ~ 1 + d2 + d4 which incorporates
  # the shift for rev(1) surgery into the intercept, pref has zero contribution
  
  # 12 weeks
  # averages over combinations
  # mu + d2[2] + bd4[nonrand]
  # mu + d2[2] + bd4[norif]
  # mu + d2[2] + bd4[rif]
  
  prop_tru <- numeric(3)
  prop_tru[1] <- (1-p_d4_entry) 
  prop_tru[2] <- p_d4_entry*p_d4_alloc
  prop_tru[3] <- p_d4_entry*p_d4_alloc
  
  p_dur_12wk <- 
    prop_tru[1] * plogis(mu + bd1[2] + bd2[2] + bd4[1]) +
    prop_tru[2] * plogis(mu + bd1[2] + bd2[2] + bd4[2]) +
    prop_tru[3] * plogis(mu + bd1[2] + bd2[2] + bd4[3]) 
  
  # 6 weeks
  
  prop_tru <- numeric(3)
  prop_tru[1] <- (1-p_d4_entry) 
  prop_tru[2] <- p_d4_entry*p_d4_alloc
  prop_tru[3] <- p_d4_entry*p_d4_alloc
  
  p_dur_06wk <- 
    prop_tru[1] * plogis(mu + bd1[2] + bd2[3] + bd4[1]) +
    prop_tru[2] * plogis(mu + bd1[2] + bd2[3] + bd4[2]) +
    prop_tru[3] * plogis(mu + bd1[2] + bd2[3] + bd4[3]) 
  
  c(rd_dur = p_dur_06wk - p_dur_12wk,
    p_dur_12wk = p_dur_12wk, p_dur_06wk = p_dur_06wk)
  
}

risk_pars_extp <- function(
    p_d1_alloc = 0.5,
    p_d2_entry = 0.7,
    p_d2_alloc = 0.5,
    p_d3_entry = 0.9,
    p_d3_alloc = 0.5,
    p_d4_entry = 0.6,
    p_d4_alloc = 0.5,
    # 40% pref for two-stage
    p_pref = 0.4,
    # log-odds
    mu = 0.7,
    # log-or
    # different baseline risk for rev
    bp = -0.4,
    bd1 = c(0, 0, 0),
    bd2 = c(0, 0.1, 0),
    bd3 = c(0, 0, -0.2),
    bd4 = c(0, 0.1, 0.2)
){
  
  # only model that is relevant is     f1_3 ~ 1 + d3 + d4 which incorporates
  # the shifts for both rev(2) surgery and pref into the intercept
  
  # 0 ext proph
  # averages over combinations
  # mu + bp + bd1[3] + bd2[2] + bd4[nonrand]
  # mu + bp + bd1[3] + bd2[2] + bd4[norif]
  # mu + bp + bd1[3] + bd2[2] + bd4[rif]
  
  prop_tru <- numeric(3)
  prop_tru[1] <- (1-p_d4_entry) 
  prop_tru[2] <- p_d4_entry*p_d4_alloc
  prop_tru[3] <- p_d4_entry*p_d4_alloc
  
  p_extp_0wk <- 
    prop_tru[1] * plogis(mu + bp + bd1[3] + bd3[2] + bd4[1]) +
    prop_tru[2] * plogis(mu + bp + bd1[3] + bd3[2] + bd4[2]) +
    prop_tru[3] * plogis(mu + bp + bd1[3] + bd3[2] + bd4[3]) 
  
  # 12 weeks
  
  prop_tru <- numeric(3)
  prop_tru[1] <- (1-p_d4_entry) 
  prop_tru[2] <- p_d4_entry*p_d4_alloc
  prop_tru[3] <- p_d4_entry*p_d4_alloc
  
  p_extp_12wk <-  
    prop_tru[1] * plogis(mu + bp + bd1[3] + bd3[3] + bd4[1]) +
    prop_tru[2] * plogis(mu + bp + bd1[3] + bd3[3] + bd4[2]) +
    prop_tru[3] * plogis(mu + bp + bd1[3] + bd3[3] + bd4[3]) 
  
  c(rd_extp = p_extp_12wk - p_extp_0wk,
    p_extp_0wk = p_extp_0wk, p_extp_12wk = p_extp_12wk)
  
}

risk_pars_choice <- function(
    p_d1_alloc = 0.5,
    p_d2_entry = 0.7,
    p_d2_alloc = 0.5,
    p_d3_entry = 0.9,
    p_d3_alloc = 0.5,
    p_d4_entry = 0.6,
    p_d4_alloc = 0.5,
    # 40% pref for two-stage
    p_pref = 0.4,
    # log-odds
    mu = 0.7,
    # log-or
    # different baseline risk for rev
    bp = -0.4,
    bd1 = c(0, 0, 0),
    bd2 = c(0, 0.1, 0),
    bd3 = c(0, 0, -0.2),
    bd4 = c(0, 0.1, 0.2)
){
  # y ~ 1 + pref + d4
  # y ~ 1 + d2 + d4
  # y ~ 1 + d3 + d4
  
  # ab choice included everywhere so risk is weighted over all models

  # we average over pref for the first model
  prop_tru <- numeric(2)
  prop_tru[1] <- (1-p_pref) 
  prop_tru[2] <- p_pref
  p_choice_norif_1 <-  
    prop_tru[1] * plogis(mu + bd4[2]) +
    prop_tru[2] * plogis(mu + bp + bd4[2]) 
  p_choice_rif_1 <-  
    prop_tru[1] * plogis(mu + bd4[3]) +
    prop_tru[2] * plogis(mu + bp + bd4[3]) 
  
  # average over bd2
  prop_tru <- numeric(3)
  prop_tru[1] <- (1-p_d2_entry) 
  prop_tru[2] <- p_d2_entry * p_d2_alloc
  prop_tru[3] <- p_d2_entry * p_d2_alloc
  
  p_choice_norif_2 <-  
    prop_tru[1] * plogis(mu + bd1[2] + bd2[1] + bd4[2]) +
    prop_tru[2] * plogis(mu + bd1[2] + bd2[2] + bd4[2]) +
    prop_tru[3] * plogis(mu + bd1[2] + bd2[3] + bd4[2]) 
  p_choice_rif_2 <-  
    prop_tru[1] * plogis(mu + bd1[2] + bd2[1] + bd4[3]) +
    prop_tru[2] * plogis(mu + bd1[2] + bd2[2] + bd4[3]) +
    prop_tru[3] * plogis(mu + bd1[2] + bd2[3] + bd4[3]) 
  
  # average over bd2
  prop_tru <- numeric(3)
  prop_tru[1] <- (1-p_d3_entry) 
  prop_tru[2] <- p_d3_entry * p_d3_alloc
  prop_tru[3] <- p_d3_entry * p_d3_alloc
  
  p_choice_norif_3 <-  
    prop_tru[1] * plogis(mu + bd1[2] + bd3[1] + bd4[2]) +
    prop_tru[2] * plogis(mu + bd1[2] + bd3[2] + bd4[2]) +
    prop_tru[3] * plogis(mu + bd1[2] + bd3[3] + bd4[2]) 
  p_choice_rif_3 <-  
    prop_tru[1] * plogis(mu + bd1[2] + bd3[1] + bd4[3]) +
    prop_tru[2] * plogis(mu + bd1[2] + bd3[2] + bd4[3]) +
    prop_tru[3] * plogis(mu + bd1[2] + bd3[3] + bd4[3]) 
  
  # model estimates need to be weighted by the relevant cohort contributions
  prop_tru <- numeric(3)
  prop_tru[1] <- p_d1_alloc
  prop_tru[2] <- p_d1_alloc*(1-p_pref)
  prop_tru[3] <- p_d1_alloc*p_pref
  
  p_choice_norif <- 
    prop_tru[1] * p_choice_norif_1 +
    prop_tru[2] * p_choice_norif_2 +
    prop_tru[3] * p_choice_norif_3
  
  p_choice_rif <- 
    prop_tru[1] * p_choice_rif_1 +
    prop_tru[2] * p_choice_rif_2 +
    prop_tru[3] * p_choice_rif_3 
  
  c(rd_choice = p_choice_rif - p_choice_norif,
    p_choice_norif_1 = p_choice_norif_1, 
    p_choice_norif_2 = p_choice_norif_2,
    p_choice_norif_3 = p_choice_norif_3, 
    p_choice_rif_1 = p_choice_rif_1,
    p_choice_rif_2 = p_choice_rif_2,
    p_choice_rif_3 = p_choice_rif_3)

}


test_pars <- function(){
  
  p_d1_alloc = 0.5
  p_d2_entry = 0.7
  p_d2_alloc = 0.5
  p_d3_entry = 0.9
  p_d3_alloc = 0.5
  p_d4_entry = 0.6
  p_d4_alloc = 0.5
  # 40% pref for two-stage
  p_pref = 0.4
  # log-odds
  mu = 0.7
  # log-or
  # different baseline risk for rev
  bp = -0.4
  bd1 = c(0, 0, 0)
  bd2 = c(0, 0.1, 0)
  bd3 = c(0, 0, -0.2)
  bd4 = c(0, 0.1, 0.2)
  
  risk_pars_surg(bd1 = c(0, 0, 0))
  
  risk_pars_dur(bd1 = c(0, 0.5, 0), 
                bd2 = c(0, 0.1, 0), 
                bd4 = c(0, 0.1, 0.2))
  
  risk_pars_extp(bd1 = c(0, 0.5, 0.1),
                 bd2 = c(0, 0.1, 0.5), 
                 bd3 = c(0, 0.4, -0.1),
                 bd4 = c(0, 0.1, 0.2))
  
  risk_pars_choice(bd4 = c(0, 0, 0))
  
  obj_f1 <- function(x){
    res <- risk_pars_surg(bd1 = c(0, x, x),
                          bd2 = bd2,
                          bd3 = bd3,
                          bd4 = bd4)
    (res["rd"] - 0)^2
  }
  obj_f2 <- function(x){
    res <- risk_pars_surg(bd1 = c(0, 0, x))
    (res["rd"] - 0)^2
  }
  obj_f3 <- function(x){
    res <- risk_pars_surg(bd1 = c(0, x, 0))
    (res["rd"] - 0)^2
  }
  
  # use to work out the null scenario in terms of risk when using 
  # log odds ratios to set up the simulation
  
  optimize(obj_f1, c(0, 1))$minimum
  
  # to get a null scenario on the risk scale when both rev(1) and rev(2) 
  # are non-zero, you need to set the logor for both to:
  (lor <- optimize(obj_f1, c(0, 1))$minimum)
  risk_pars_surg(bd1 = c(0, lor, lor))
  # to get a null scenario on the risk scale when rev(1) is zero, you need 
  # to set the rev(2) logor to:
  optimize(obj_f2, c(0, 1))$minimum
  # to get a null scenario on the risk scale when rev(2) is zero, you need 
  # to set rev(1) logor to:
  optimize(obj_f3, c(0, 1))$minimum
  
}



multi_model_approach <- function(){
  
  N_sim <- 1000
  mc_cores <- 4
  
  p_d1_alloc = 0.5
  p_d2_entry = 0.7
  p_d2_alloc = 0.5
  p_d3_entry = 0.9
  p_d3_alloc = 0.5
  p_d4_entry = 0.6
  p_d4_alloc = 0.5
  
  p_pref = 0.4
  # log-odds
  mu = 0.7
  # log-or
  # different baseline risk for rev
  bp = -0.4
  bd1 = c(0, 0, 0)
  bd2 = c(0, 0.1, 0)
  bd3 = c(0, 0, -0.2)
  bd4 = c(0, 0.1, 0.2)

  r <- pbapply::pblapply(X=1:N_sim, cl = mc_cores, FUN = function(ix){
    
    N <- 1e3
    d <- data.table(
      # ctl/trt allocations ignoring dependencies
      d1_alloc = rbinom(N, 1, p_d1_alloc),
      # 70% entery d2
      d2_entry = rbinom(N, 1, p_d2_entry),
      d2_alloc = rbinom(N, 1, p_d2_alloc),
      # 90% enter d3
      d3_entry = rbinom(N, 1, p_d3_entry),
      d3_alloc = rbinom(N, 1, p_d3_alloc),
      # 60% enter d4
      d4_entry = rbinom(N, 1, p_d4_entry),
      d4_alloc = rbinom(N, 1, p_d4_alloc),
      # preference directs type of revision
      pref = rbinom(N, 1, p_pref)
    )
    
    # dair gets dair, revision gets split
    d[d1_alloc == 0, d1 := 1]
    d[d1_alloc == 1 & pref == 0, d1 := 2]
    d[d1_alloc == 1 & pref == 1, d1 := 3]
    
    d[d1 == 2 & d2_entry == 0, d2 := 1]
    d[d1 == 2 & d2_entry == 1, d2 := 2 + d2_alloc]
    
    d[d1 == 3 & d3_entry == 0, d3 := 1]
    d[d1 == 3 & d3_entry == 1, d3 := 2 + d3_alloc]
    
    d[d4_entry == 0, d4 := 1]
    d[d4_entry == 1, d4 := 2 + d4_alloc]
    
    d[, .N, keyby = .(pref, d1, d2, d3, d4)]

    # bd1 is irrelevant since d1 == 1 is the ref group, fixed at zero
    d[d1 == 1, eta := mu + bp*pref + bd4[d4]]
    # pref is irrelevant as d1 = 2 only occurs if pref = 0
    d[d1 == 2, eta := mu + bd1[2] + bd2[d2] + bd4[d4]]
    # but here pref is relevant as d1 = 3 only if pref = 1
    d[d1 == 3, eta := mu + bp*pref + bd1[3] + bd3[d3] + bd4[d4]]
    
    d[, `:=`(d1 = factor(d1), d2 = factor(d2), d3 = factor(d3), d4 = factor(d4))]
    
    d[, p := plogis(eta)]
    d[, y := rbinom(N, 1, p)]
    
    f1_1 <- glm(
      y ~ 1 + pref + d4 , data = d, subset= d1==1, family = binomial,
      control = glm.control(epsilon = 1e-4, maxit = 100, trace = FALSE))
    # pref is zero for everyone
    f1_2 <- glm(
      y ~ 1 + d2 + d4, data = d, subset= d1==2, family = binomial,
      control = glm.control(epsilon = 1e-4, maxit = 100, trace = FALSE))
    # effect of pref will get rolled into the intercept
    f1_3 <- glm(
      y ~ 1 + d3 + d4, data = d, subset= d1==3, family = binomial,
      control = glm.control(epsilon = 1e-4, maxit = 100, trace = FALSE))
    
    # surgical domain
    
    # dair
    d_dair <- copy(d)
    d_dair[, d1 := factor(1)]
    # d4 stays at whatever it was
    d_dair[, p_hat := predict(f1_1, newdata = d_dair, type = "response")]
    p_dair_hat <- mean(d_dair$p_hat)
    
    # rev(1)
    d_rev_1 <- copy(d)
    d_rev_1[, d1 := factor(2)]
    d_rev_1[, pref := 0]
    d_rev_1[d2_entry == 1, d2 := factor(d[d2_entry == 1, d2_alloc + 2], levels = 1:3)]
    d_rev_1[d2_entry == 0, d2 := factor(1, levels = 1:3)]
    # and d4 stays at whatever it was
    d_rev_1[, p_hat := predict(f1_2, newdata = d_rev_1, type = "response")]
    p_rev_1_hat <- mean(d_rev_1$p_hat)
    
    # rev(2)
    d_rev_2 <- copy(d)
    d_rev_2[, d1 := factor(3)]
    # pref doesn't matter because the model had to roll it up into the intercept
    # d_rev_2[, pref := 1]
    d_rev_2[d3_entry == 1, d3 := factor(d[d3_entry == 1, d3_alloc + 2], levels = 1:3)]
    d_rev_2[d3_entry == 0, d3 := factor(1, levels = 1:3)]
    # and d4 stays at whatever it was
    d_rev_2[, p_hat := predict(f1_3, newdata = d_rev_2, type = "response")]
    p_rev_2_hat <- mean(d_rev_2$p_hat)
    
    # now aggregate the one and two-stage revision
    d_prop <- d[, .(prop = .N/nrow(d)), keyby = pref]
    p_rev_hat <- p_rev_1_hat * d_prop[pref == 0, prop] + p_rev_2_hat * d_prop[pref == 1, prop]
    
    rd_hat <- p_rev_hat - p_dair_hat

    
    d_uniq <- unique(d[, .(pref, d1, d2, d3, d4, p)])
    setkey(d_uniq, pref, d1, d2, d3, d4)
    d_uniq[d1 == 1, p_hat := predict(f1_1, newdata = d_uniq[d1 == 1], type = "response")]
    d_uniq[d1 == 2, p_hat := predict(f1_2, newdata = d_uniq[d1 == 2], type = "response")]
    d_uniq[d1 == 3, p_hat := predict(f1_3, newdata = d_uniq[d1 == 3], type = "response")]
    
    # don't really need to do this for all, but:
    risk_surg <- risk_pars_surg(
      p_d1_alloc, 
      p_d2_entry, p_d2_alloc,
      p_d3_entry, p_d3_alloc,
      p_d4_entry, p_d4_alloc, 
      p_pref, 
      mu, bp, bd1, bd2, bd3, bd4)
    
    list(
      d_uniq = d_uniq,
      d_g = data.table(
        p_dair_ref = risk_surg["p_dair"],
        p_dair_hat = p_dair_hat,
        p_rev_1_ref = risk_surg["p_rev_1"],
        p_rev_1_hat = p_rev_1_hat,
        p_rev_2_ref = risk_surg["p_rev_2"],
        p_rev_2_hat = p_rev_2_hat,
        p_rev_ref = risk_surg["p_rev"],
        p_rev_hat = p_rev_hat,
        rd_ref = risk_surg["rd"],
        rd_hat = rd_hat
      )
    )
    
  })
  
  
  d_g <- rbindlist(lapply(r, function(z) z$d_g))
  
  d_g[, bias_rd := rd_hat - rd_ref]
  d_g[, bias_p_dair := p_dair_hat - p_dair_ref]
  d_g[, bias_p_rev := p_rev_hat - p_rev_ref]
  d_g[, bias_p_rev_1 := p_rev_1_hat - p_rev_1_ref]
  d_g[, bias_p_rev_2 := p_rev_2_hat - p_rev_2_ref]
  
  d_fig <- d_g[, .(bias_rd, bias_p_dair, bias_p_rev)]
  d_fig <- melt(d_fig, measure.vars = names(d_fig))
  
  p1 <- ggplot(d_fig, aes(x = value)) +
    geom_density() +
    geom_vline(data = d_fig[, .(mu = mean(value)), keyby = variable],
               aes(xintercept = mu), lwd = 0.2) + 
    facet_wrap(~variable, nrow = 1)
  #
  
  d_fig <- d_g[, .(p_dair_hat, p_rev_hat)]
  d_fig <- melt(d_fig, measure.vars = names(d_fig))
  
  p2 <- ggplot(d_fig, aes(x = value)) +
    geom_density() +
    geom_vline(data = d_fig[, .(mu = mean(value)), keyby = variable],
               aes(xintercept = mu), lwd = 0.2) +
    facet_wrap(~variable,  nrow = 1)
  
  d_fig <- d_g[, .(rd_hat)]
  p3 <- ggplot(d_fig, aes(x = rd_hat)) +
    geom_density() +
    geom_vline(data = d_fig[, .(mu = mean(rd_hat))],
               aes(xintercept = mu), lwd = 0.2) 
  
  p1 / (p2 + p3)
  
}




single_model_approach <- function(){
  
  N_sim <- 1000
  mc_cores <- 4
  
  r <- pbapply::pblapply(X=1:N_sim, cl = mc_cores, FUN = function(ix){
    
    N <- 1e4
    d <- data.table(
      # ctl/trt allocations ignoring dependencies
      d1_alloc = rbinom(N, 1, 0.5),
      # 70% entery d2
      d2_entry = rbinom(N, 1, 0.7),
      d2_alloc = rbinom(N, 1, 0.5),
      # 90% enter d3
      d3_entry = rbinom(N, 1, 0.9),
      d3_alloc = rbinom(N, 1, 0.5),
      # 60% enter d4
      d4_entry = rbinom(N, 1, 0.6),
      d4_alloc = rbinom(N, 1, 0.5),
      # preference directs type of revision
      pref = rbinom(N, 1, 0.4)
    )
    
    # dair gets dair, revision gets split
    d[d1_alloc == 0, d1 := 1]
    d[d1_alloc == 1 & pref == 0, d1 := 2]
    d[d1_alloc == 1 & pref == 1, d1 := 3]
    
    d[d1 == 2 & d2_entry == 0, d2 := 1]
    d[d1 == 2 & d2_entry == 1, d2 := 2 + d2_alloc]
    
    d[d1 == 3 & d3_entry == 0, d3 := 1]
    d[d1 == 3 & d3_entry == 1, d3 := 2 + d3_alloc]
    
    d[d4_entry == 0, d4 := 1]
    d[d4_entry == 1, d4 := 2 + d4_alloc]
    
    d[, .N, keyby = .(pref, d1, d2, d3, d4)]
    
    mu <- 0.7
    # different baseline risk for rev
    bp <- -0.4
    # actually the null case when you think about it on the risk scale
    bd1 <- c(0, 0.01956613, 0.01956613)
    bd2 <- c(0, 0.1, 0)
    bd3 <- c(0, 0, -0.2)
    bd4 <- c(0, 0.1, 0.2)
    
    # bp is irrelevant since it is multiplied by 0
    d[d1 == 1, eta := mu + bd1[1] + bd4[d4]]
    d[d1 == 2, eta := mu + bd1[2] + bd2[d2] + bd4[d4]]
    # bp is relevant since it is multiplied by 1
    d[d1 == 3, eta := mu + bp*pref + bd1[3] + bd3[d3] + bd4[d4]]
    
    # these are the undefined records, those for whom there is no d2 nor d3
    d[is.na(d2), d2 := 4]
    d[is.na(d3), d3 := 4]
    
    d[, `:=`(d1 = factor(d1), d2 = factor(d2), d3 = factor(d3),
             d4 = factor(d4))]
    
    d[, p := plogis(eta)]
    d[, y := rbinom(N, 1, p)]
    
    f2 <- glm(
      y ~ 1 + pref + d1 + d2 + d3 + d4, data = d, family = binomial,
      control = glm.control(epsilon = 1e-4, maxit = 100, trace = FALSE))
    
    # summary(f2)
    # lapply(list(f1_1, f1_2, f1_3), function(z) summary(z)$coef)
    # dair
    d_dair <- copy(d)
    d_dair[, d1 := factor(1)]
    d_dair[, d2 := factor(4)]
    d_dair[, d3 := factor(4)]
    # d4 as it was
    d_dair[, p_hat := predict(f2, newdata = d_dair, type = "response")]
    
    # see earlier, noting that bd2[4] and bd3[4] do not exist
    prop_tru <- c(0.24, 0.18, 0.18, 0.16, 0.12, 0.12)
    p_dair_ref <- 
      prop_tru[1] * plogis(mu + bd4[1]) + 
      prop_tru[2] * plogis(mu + bd4[2]) + 
      prop_tru[3] * plogis(mu + bd4[3]) + 
      prop_tru[4] * plogis(mu + bp + bd4[1]) + 
      prop_tru[5] * plogis(mu + bp + bd4[2]) + 
      prop_tru[6] * plogis(mu + bp + bd4[3])
    p_dair_hat <- mean(d_dair$p_hat)
    
    # rev(1)
    d_rev_1 <- copy(d)
    d_rev_1[, d1 := factor(2)]
    d_rev_1[, pref := 0]
    d_rev_1[d2_entry == 1, d2 := factor(d[d2_entry == 1, d2_alloc + 2], levels = 1:4)]
    d_rev_1[d2_entry == 0, d2 := factor(1)]
    d_rev_1[, d3 := factor(4)]
    # d4 as is
    d_rev_1[, p_hat := predict(f2, newdata = d_rev_1, type = "response")]
    # see earlier
    prop_tru <- c(0.12, 0.09, 0.09,
                  0.14, 0.105, 0.105,
                  0.14, 0.105, 0.105)
    p_rev_1_ref <- 
      prop_tru[1] * plogis(mu + bd1[2] + bd2[1] + bd4[1]) + 
      prop_tru[2] * plogis(mu + bd1[2] + bd2[1] + bd4[2]) + 
      prop_tru[3] * plogis(mu + bd1[2] + bd2[1] + bd4[3]) + 
      prop_tru[4] * plogis(mu + bd1[2] + bd2[2] + bd4[1]) + 
      prop_tru[5] * plogis(mu + bd1[2] + bd2[2] + bd4[2]) + 
      prop_tru[6] * plogis(mu + bd1[2] + bd2[2] + bd4[3]) + 
      prop_tru[7] * plogis(mu + bd1[2] + bd2[3] + bd4[1]) + 
      prop_tru[8] * plogis(mu + bd1[2] + bd2[3] + bd4[2]) + 
      prop_tru[9] * plogis(mu + bd1[2] + bd2[3] + bd4[3])  
    p_rev_1_hat <- mean(d_rev_1$p_hat)
    
    # rev(2)
    d_rev_2 <- copy(d)
    d_rev_2[, d1 := factor(3)]
    d_rev_2[, pref := 1]
    d_rev_2[, d2 := factor(4)]
    d_rev_2[d3_entry == 1, d3 := factor(d[d3_entry == 1, d3_alloc + 2], levels = 1:4)]
    d_rev_2[d3_entry == 0, d3 := factor(1, levels = 1:4)]
    # d4 as is
    d_rev_2[, p_hat := predict(f2, newdata = d_rev_2, type = "response")]
    
    prop_tru <- c(0.04, 0.03, 0.03,
                  0.18, 0.135, 0.135,
                  0.18, 0.135, 0.135)
    # contribution of non-rand, bp needs to be included here
    p_rev_2_ref <- 
      prop_tru[1] * plogis(mu + bp + bd1[3] + bd3[1] + bd4[1]) + 
      prop_tru[2] * plogis(mu + bp + bd1[3] + bd3[1] + bd4[2]) + 
      prop_tru[3] * plogis(mu + bp + bd1[3] + bd3[1] + bd4[3]) + 
      prop_tru[4] * plogis(mu + bp + bd1[3] + bd3[2] + bd4[1]) + 
      prop_tru[5] * plogis(mu + bp + bd1[3] + bd3[2] + bd4[2]) + 
      prop_tru[6] * plogis(mu + bp + bd1[3] + bd3[2] + bd4[3]) + 
      prop_tru[7] * plogis(mu + bp + bd1[3] + bd3[3] + bd4[1]) + 
      prop_tru[8] * plogis(mu + bp + bd1[3] + bd3[3] + bd4[2]) + 
      prop_tru[9] * plogis(mu + bp + bd1[3] + bd3[3] + bd4[3])  
    p_rev_2_hat <- mean(d_rev_2$p_hat)
    
    prop_tru <- c(0.6, 0.4)
    p_rev_ref <- p_rev_1_ref * prop_tru[1] + p_rev_2_ref * prop_tru[2]
    
    d_prop <- d[, .(prop = .N/nrow(d)), keyby = pref]
    p_rev_hat <- p_rev_1_hat * d_prop[pref == 0, prop] + p_rev_2_hat * d_prop[pref == 1, prop]
    
    rd_ref <- p_rev_ref - p_dair_ref
    rd_hat <- p_rev_hat - p_dair_hat
    
    d_uniq <- unique(d[, .(pref, d1, d2, d3, d4, p)])
    setkey(d_uniq, pref, d1, d2, d3, d4)
    d_uniq[d1 == 1, p_hat := predict(f2, newdata = d_uniq[d1 == 1], type = "response")]
    d_uniq[d1 == 2, p_hat := predict(f2, newdata = d_uniq[d1 == 2], type = "response")]
    d_uniq[d1 == 3, p_hat := predict(f2, newdata = d_uniq[d1 == 3], type = "response")]
    
    
    list(
      d_uniq = d_uniq,
      d_g = data.table(
        p_dair_ref = p_dair_ref,
        p_dair_hat = p_dair_hat,
        p_rev_1_ref = p_rev_1_ref,
        p_rev_1_hat = p_rev_1_hat,
        p_rev_2_ref = p_rev_2_ref,
        p_rev_2_hat = p_rev_2_hat,
        p_rev_ref = p_rev_ref,
        p_rev_hat = p_rev_hat,
        rd_ref = rd_ref,
        rd_hat = rd_hat
      )
    )
  })
  
  d_g <- rbindlist(lapply(r, function(z) z$d_g))
  
  d_g[, bias_rd := rd_hat - rd_ref]
  d_g[, bias_p_dair := p_dair_hat - p_dair_ref]
  d_g[, bias_p_rev := p_rev_hat - p_rev_ref]
  d_g[, bias_p_rev_1 := p_rev_1_hat - p_rev_1_ref]
  d_g[, bias_p_rev_2 := p_rev_2_hat - p_rev_2_ref]
  
  d_fig <- d_g[, .(p_dair_hat, p_rev_1_hat, p_rev_2_hat, p_rev_hat)]
  d_fig <- melt(d_fig, measure.vars = names(d_fig))
  
  p1 <- ggplot(d_fig, aes(x = value)) +
    geom_density() +
    geom_vline(data = d_fig[, .(mu = mean(value)), keyby = variable],
               aes(xintercept = mu), lwd = 0.2) + 
    facet_wrap(~variable, nrow = 1)
  #
  
  d_fig <- d_g[, .(p_dair_hat, p_rev_hat)]
  d_fig <- melt(d_fig, measure.vars = names(d_fig))
  
  p2 <- ggplot(d_fig, aes(x = value)) +
    geom_density() +
    geom_vline(data = d_fig[, .(mu = mean(value)), keyby = variable],
               aes(xintercept = mu), lwd = 0.2) +
    facet_wrap(~variable,  nrow = 1)
  
  d_fig <- d_g[, .(rd_hat)]
  p3 <- ggplot(d_fig, aes(x = rd_hat)) +
    geom_density() +
    geom_vline(data = d_fig[, .(mu = mean(rd_hat))],
               aes(xintercept = mu), lwd = 0.2) 
  
  p1 / (p2 + p3)
  
}

