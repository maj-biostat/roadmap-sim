library(data.table)
library(pbapply)
library(ggplot2)
library(patchwork)


risk_pars_surg <- function(
    # parameters would be those from the late acute list
    p_d1_alloc = 0.5,
    p_d2_entry = 0.7,
    p_d2_alloc = 0.5,
    p_d3_entry = 0.9,
    p_d3_alloc = 0.5,
    p_d4_entry = 0.6,
    p_d4_alloc = 0.5,
    # pref towards two-stage if assigned to rev
    p_pref = 0.4,
    # log-odds
    bs = c(0.9, 0.8, 0.7),
    # log-or
    # different baseline risk for rev
    bp = -0.4,
    bd1 = c(0, 0, 0),
    bd2 = c(0, 0.1, 0),
    bd3 = c(0, 0, -0.2),
    bd4 = c(0, 0.1, 0.2)
  ){
  
  # dair
  # averages over combinations from f1_1:  ~ -1 + s + pref + d4 (in late acute)
  # bs[2] + bd4[nonrand]
  # bs[2] + bd4[norif]
  # bs[2] + bd4[rif]
  # bs[2] + bp + bd4[nonrand]
  # bs[2] + bp + bd4[norif]
  # bs[2] + bp + bd4[rif]
  # For example:
  # if 60% prefer rev(1) (if they had recv rev) and 60% enter ab choice
  # then 60% weight is given to the first three and 40% to the last.
  # 60% enter ab choice unconditionally, so we have 
  # 0.6 * 0.6 = 0.36 allocated across norif and rif and the rest go into nonrand
  # in the second three combinations, using the same logic
  # 0.6 * 0.4 = 0.24 allocated to no rif and rif and the rest go into nonrand
  # so we would want:
  # 0.24 * invlogit(bs[2] + bd4[nonrand])
  # 0.18 * invlogit(bs[2] + bd4[norif])
  # 0.18 * invlogit(bs[2] + bd4[rif])
  # 0.16 * invlogit(bs[2] + bp + bd4[nonrand])
  # 0.12 * invlogit(bs[2] + bp + bd4[norif])
  # 0.12 * invlogit(bs[2] + bp + bd4[rif])
  
  prop_tru <- numeric(6)
  prop_tru[1] <- (1-p_pref) - (p_d4_entry * (1-p_pref))
  prop_tru[2] <- p_d4_entry * (1-p_pref)/2
  prop_tru[3] <- p_d4_entry * (1-p_pref)/2
  prop_tru[4] <- (p_pref) - (p_d4_entry * (p_pref))
  prop_tru[5] <- p_d4_entry * (p_pref)/2
  prop_tru[6] <- p_d4_entry * (p_pref)/2
  
  p_dair <- 
    prop_tru[1] * plogis(bs[2] + bd4[1]) + 
    prop_tru[2] * plogis(bs[2] + bd4[2]) + 
    prop_tru[3] * plogis(bs[2] + bd4[3]) + 
    prop_tru[4] * plogis(bs[2] + bp + bd4[1]) + 
    prop_tru[5] * plogis(bs[2] + bp + bd4[2]) + 
    prop_tru[6] * plogis(bs[2] + bp + bd4[3])
  
  # rev 1
  # averages over combinations from f1_2:  ~ -1 + s + d2 + d4
  # bs[2] + bd1[2] + bd2[1] + bd4[nonrand]
  # bs[2] + bd1[2] + bd2[1] + bd4[norif]
  # bs[2] + bd1[2] + bd2[1] + bd4[rif]    
  # bs[2] + bd1[2] + bd2[2] + bd4[nonrand]
  # bs[2] + bd1[2] + bd2[2] + bd4[norif]
  # bs[2] + bd1[2] + bd2[2] + bd4[rif]
  # bs[2] + bd1[2] + bd2[3] + bd4[nonrand]
  # bs[2] + bd1[2] + bd2[3] + bd4[norif]
  # bs[2] + bd1[2] + bd2[3] + bd4[rif]
  # For example, 
  # 70% of rev(1) enter into d2, 60% enter d4
  # 0.3 enter the nonrand treatment in d2
  # 0.6 * 0.3 = 0.18 allocated for norif and rif (rest 0.3 - 0.18 for nonrand)
  # Below the 0.35 are the 0.7*0.5 allocations
  # 0.6 * 0.35 = 0.21 allocated to norif and rif (rest 0.35 - 0.21 for nonrand)
  # 0.6 * 0.35 = 0.21 allocated to norif and rif (rest 0.35 - 0.21 for nonrand)
  
  prop_tru <- numeric(9)
  # non-rand d2
  prop_tru[1] <- (1-p_d2_entry) - (p_d4_entry * (1-p_d2_entry))
  prop_tru[2] <- p_d4_entry * (1-p_d2_entry) / 2
  prop_tru[3] <- p_d4_entry * (1-p_d2_entry) / 2
  # rand d2 level 2
  prop_tru[4] <- (p_d2_entry/2) - (p_d4_entry * (p_d2_entry/2))
  prop_tru[5] <- p_d4_entry * (p_d2_entry/2) / 2
  prop_tru[6] <- p_d4_entry * (p_d2_entry/2) / 2
  # rand d2 level 3
  prop_tru[7] <- (p_d2_entry/2) - (p_d4_entry * (p_d2_entry/2))
  prop_tru[8] <- p_d4_entry * (p_d2_entry/2) / 2
  prop_tru[9] <- p_d4_entry * (p_d2_entry/2) / 2
  
  p_rev_1 <- 
    prop_tru[1] * plogis(bs[2] + bd1[2] + bd2[1] + bd4[1]) + 
    prop_tru[2] * plogis(bs[2] + bd1[2] + bd2[1] + bd4[2]) + 
    prop_tru[3] * plogis(bs[2] + bd1[2] + bd2[1] + bd4[3]) + 
    prop_tru[4] * plogis(bs[2] + bd1[2] + bd2[2] + bd4[1]) + 
    prop_tru[5] * plogis(bs[2] + bd1[2] + bd2[2] + bd4[2]) + 
    prop_tru[6] * plogis(bs[2] + bd1[2] + bd2[2] + bd4[3]) + 
    prop_tru[7] * plogis(bs[2] + bd1[2] + bd2[3] + bd4[1]) + 
    prop_tru[8] * plogis(bs[2] + bd1[2] + bd2[3] + bd4[2]) + 
    prop_tru[9] * plogis(bs[2] + bd1[2] + bd2[3] + bd4[3])  
  
  # rev 2
  # averages over combinations from f1_2:   ~ -1 + bs + d3 + d4
  # pref is included as intercept = mu + bp for this cohort
  # bs[2] + bp + bd1[3] + bd3[1] + bd4[nonrand]
  # bs[2] + bp + bd1[3] + bd3[1] + bd4[norif]
  # bs[2] + bp + bd1[3] + bd3[1] + bd4[rif]    
  # bs[2] + bp + bd1[3] + bd3[2] + bd4[nonrand]
  # bs[2] + bp + bd1[3] + bd3[2] + bd4[norif]
  # bs[2] + bp + bd1[3] + bd3[2] + bd4[rif]
  # bs[2] + bp + bd1[3] + bd3[3] + bd4[nonrand]
  # bs[2] + bp + bd1[3] + bd3[3] + bd4[norif]
  # bs[2] + bp + bd1[3] + bd3[3] + bd4[rif]
  # 90% of rev(2) enter into d3, 60% enter d4
  # 0.6 * 0.1 = 0.06 allocated to norif and rif (rest 0.1 - 0.06 for nonrand)
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
    prop_tru[1] * plogis(bs[2] + bp + bd1[3] + bd3[1] + bd4[1]) + 
    prop_tru[2] * plogis(bs[2] + bp + bd1[3] + bd3[1] + bd4[2]) + 
    prop_tru[3] * plogis(bs[2] + bp + bd1[3] + bd3[1] + bd4[3]) + 
    prop_tru[4] * plogis(bs[2] + bp + bd1[3] + bd3[2] + bd4[1]) + 
    prop_tru[5] * plogis(bs[2] + bp + bd1[3] + bd3[2] + bd4[2]) + 
    prop_tru[6] * plogis(bs[2] + bp + bd1[3] + bd3[2] + bd4[3]) + 
    prop_tru[7] * plogis(bs[2] + bp + bd1[3] + bd3[3] + bd4[1]) + 
    prop_tru[8] * plogis(bs[2] + bp + bd1[3] + bd3[3] + bd4[2]) + 
    prop_tru[9] * plogis(bs[2] + bp + bd1[3] + bd3[3] + bd4[3])  
  
  # rev
  prop_tru <- c(1 - p_pref, p_pref)
  p_rev <- prop_tru[1] * p_rev_1 + prop_tru[2] * p_rev_2
  
  
  c(rd = p_rev - p_dair,
    p_dair = p_dair, p_rev = p_rev, 
    p_rev_1 = p_rev_1, p_rev_2 = p_rev_2)
  
}


risk_pars_dur <- function(
    p_s_alloc = c(0.3, 0.5, 0.2),
    l_e, l_l, l_c,
    # log-odds
    bs = c(0.9, 0.8, 0.7),
    # log-or
    # different baseline risk for rev
    bp = -0.4,
    bd1 = c(0, 0, 0),
    bd2 = c(0, 0.1, 0),
    bd3 = c(0, 0, -0.2),
    bd4 = c(0, 0.1, 0.2)
){
  
  # only model that is relevant is f1_2 ~ -1 + s + d2 + d4 which incorporates
  # the shift for rev(1) surgery into the intercept, pref has zero contribution
  
  # 12 weeks
  # averages over combinations
  # bs[1] + d2[2] + bd4[nonrand]
  # bs[1] + d2[2] + bd4[norif]
  # bs[1] + d2[2] + bd4[rif]
  # bs[2] + d2[2] + bd4[nonrand]
  # bs[2] + d2[2] + bd4[norif]
  # bs[2] + d2[2] + bd4[rif]
  # bs[3] + d2[2] + bd4[nonrand]
  # bs[3] + d2[2] + bd4[norif]
  # bs[3] + d2[2] + bd4[rif]
  
  prop_tru <- numeric(9)
  prop_tru[1] <- (p_s_alloc[1]) * (1-l_e$p_d4_entry)
  prop_tru[2] <- (p_s_alloc[1]) * l_e$p_d4_entry*l_e$p_d4_alloc
  prop_tru[3] <- (p_s_alloc[1]) * l_e$p_d4_entry*l_e$p_d4_alloc
  prop_tru[4] <- (p_s_alloc[2]) * (1-l_l$p_d4_entry)
  prop_tru[5] <- (p_s_alloc[2]) * l_l$p_d4_entry*l_l$p_d4_alloc
  prop_tru[6] <- (p_s_alloc[2]) * l_l$p_d4_entry*l_l$p_d4_alloc
  prop_tru[7] <- (p_s_alloc[3]) * (1-l_e$p_d4_entry)
  prop_tru[8] <- (p_s_alloc[3]) * l_c$p_d4_entry*l_c$p_d4_alloc
  prop_tru[9] <- (p_s_alloc[3]) * l_c$p_d4_entry*l_c$p_d4_alloc
  
  p_dur_12wk <- 
    prop_tru[1] * plogis(bs[1] + bd1[2] + bd2[2] + bd4[1]) +
    prop_tru[2] * plogis(bs[1] + bd1[2] + bd2[2] + bd4[2]) +
    prop_tru[3] * plogis(bs[1] + bd1[2] + bd2[2] + bd4[3]) +
    prop_tru[4] * plogis(bs[2] + bd1[2] + bd2[2] + bd4[3]) +
    prop_tru[5] * plogis(bs[2] + bd1[2] + bd2[2] + bd4[3]) +
    prop_tru[6] * plogis(bs[2] + bd1[2] + bd2[2] + bd4[3]) +
    prop_tru[7] * plogis(bs[3] + bd1[2] + bd2[2] + bd4[3]) +
    prop_tru[8] * plogis(bs[3] + bd1[2] + bd2[2] + bd4[3]) +
    prop_tru[9] * plogis(bs[3] + bd1[2] + bd2[2] + bd4[3])   
  
  # same weights used for the 06 under 1:1 allocation
  
  p_dur_6wk <- 
    prop_tru[1] * plogis(bs[1] + bd1[2] + bd2[3] + bd4[1]) +
    prop_tru[2] * plogis(bs[1] + bd1[2] + bd2[3] + bd4[2]) +
    prop_tru[3] * plogis(bs[1] + bd1[2] + bd2[3] + bd4[3]) +
    prop_tru[4] * plogis(bs[2] + bd1[2] + bd2[3] + bd4[3]) +
    prop_tru[5] * plogis(bs[2] + bd1[2] + bd2[3] + bd4[3]) +
    prop_tru[6] * plogis(bs[2] + bd1[2] + bd2[3] + bd4[3]) +
    prop_tru[7] * plogis(bs[3] + bd1[2] + bd2[3] + bd4[3]) +
    prop_tru[8] * plogis(bs[3] + bd1[2] + bd2[3] + bd4[3]) +
    prop_tru[9] * plogis(bs[3] + bd1[2] + bd2[3] + bd4[3])   
  
  c(rd_dur = p_dur_6wk - p_dur_12wk,
    p_dur_12wk = p_dur_12wk, p_dur_6wk = p_dur_6wk)
  
}

risk_pars_extp <- function(
    p_s_alloc = c(0.3, 0.5, 0.2),
    l_e, l_l, l_c,
    # log-odds
    bs = c(0.9, 0.8, 0.7),
    # log-or
    # different baseline risk for rev
    bp = -0.4,
    bd1 = c(0, 0, 0),
    bd2 = c(0, 0.1, 0),
    bd3 = c(0, 0, -0.2),
    bd4 = c(0, 0.1, 0.2)
){
  
  # only model that is relevant is     f1_3 ~ -1 + s + d3 + d4 which incorporates
  # the shifts for both rev(2) surgery and pref into the intercept
  
  prop_tru <- numeric(9)
  prop_tru[1] <- (p_s_alloc[1]) * (1-l_e$p_d4_entry)
  prop_tru[2] <- (p_s_alloc[1]) * l_e$p_d4_entry*l_e$p_d4_alloc
  prop_tru[3] <- (p_s_alloc[1]) * l_e$p_d4_entry*l_e$p_d4_alloc
  prop_tru[4] <- (p_s_alloc[2]) * (1-l_l$p_d4_entry)
  prop_tru[5] <- (p_s_alloc[2]) * l_l$p_d4_entry*l_l$p_d4_alloc
  prop_tru[6] <- (p_s_alloc[2]) * l_l$p_d4_entry*l_l$p_d4_alloc
  prop_tru[7] <- (p_s_alloc[3]) * (1-l_e$p_d4_entry)
  prop_tru[8] <- (p_s_alloc[3]) * l_c$p_d4_entry*l_c$p_d4_alloc
  prop_tru[9] <- (p_s_alloc[3]) * l_c$p_d4_entry*l_c$p_d4_alloc
  
  p_extp_0wk <- 
    prop_tru[1] * plogis(bs[1] + bd1[2] + bd3[2] + bd4[1]) +
    prop_tru[2] * plogis(bs[1] + bd1[2] + bd3[2] + bd4[2]) +
    prop_tru[3] * plogis(bs[1] + bd1[2] + bd3[2] + bd4[3]) +
    prop_tru[4] * plogis(bs[2] + bd1[2] + bd3[2] + bd4[3]) +
    prop_tru[5] * plogis(bs[2] + bd1[2] + bd3[2] + bd4[3]) +
    prop_tru[6] * plogis(bs[2] + bd1[2] + bd3[2] + bd4[3]) +
    prop_tru[7] * plogis(bs[3] + bd1[2] + bd3[2] + bd4[3]) +
    prop_tru[8] * plogis(bs[3] + bd1[2] + bd3[2] + bd4[3]) +
    prop_tru[9] * plogis(bs[3] + bd1[2] + bd3[2] + bd4[3])   
  
  # same weights used under 1:1 allocation
  
  p_extp_12wk <- 
    prop_tru[1] * plogis(bs[1] + bd1[2] + bd3[3] + bd4[1]) +
    prop_tru[2] * plogis(bs[1] + bd1[2] + bd3[3] + bd4[2]) +
    prop_tru[3] * plogis(bs[1] + bd1[2] + bd3[3] + bd4[3]) +
    prop_tru[4] * plogis(bs[2] + bd1[2] + bd3[3] + bd4[3]) +
    prop_tru[5] * plogis(bs[2] + bd1[2] + bd3[3] + bd4[3]) +
    prop_tru[6] * plogis(bs[2] + bd1[2] + bd3[3] + bd4[3]) +
    prop_tru[7] * plogis(bs[3] + bd1[2] + bd3[3] + bd4[3]) +
    prop_tru[8] * plogis(bs[3] + bd1[2] + bd3[3] + bd4[3]) +
    prop_tru[9] * plogis(bs[3] + bd1[2] + bd3[3] + bd4[3])   
  
  c(rd_extp = p_extp_12wk - p_extp_0wk,
    p_extp_0wk = p_extp_0wk, p_extp_12wk = p_extp_12wk)
  
}

risk_pars_choice <- function(
    p_s_alloc = c(0.3, 0.5, 0.2),
    l_e, l_l, l_c,
    # log-odds
    bs = c(0.9, 0.8, 0.7),
    # log-or
    # different baseline risk for rev
    bp = -0.4,
    bd1 = c(0, 0, 0),
    bd2 = c(0, 0.1, 0),
    bd3 = c(0, 0, -0.2),
    bd4 = c(0, 0.1, 0.2)
){
  
  
  # ab choice included everywhere, so risk is weighted over all models

  # first model
  # y ~ -1 + s + pref + d4  for d1 = 1 (all silo)
  # average over silo and pref for the first model
  # combinations
  # bs[1] + bd4[norif]
  # bs[2] + bd4[norif]
  # bs[3] + bd4[norif]
  # bs[1] + bp + bd4[norif]
  # bs[2] + bp + bd4[norif]
  # bs[3] + bp + bd4[norif]
  
  prop_tru <- numeric(6)
  prop_tru[1] <- p_s_alloc[1] * (1-l_e$p_pref)
  prop_tru[2] <- p_s_alloc[2] * (1-l_l$p_pref)
  prop_tru[3] <- p_s_alloc[3] * (1-l_c$p_pref)
  prop_tru[4] <- p_s_alloc[1] * l_e$p_pref
  prop_tru[5] <- p_s_alloc[2] * l_l$p_pref
  prop_tru[6] <- p_s_alloc[3] * l_c$p_pref
  
  p_choice_norif_1 <-  
    prop_tru[1] * plogis(bs[1] + bd4[2]) +
    prop_tru[2] * plogis(bs[2] + bd4[2]) +
    prop_tru[3] * plogis(bs[3] + bd4[2]) +
    prop_tru[4] * plogis(bs[1] + bp + bd4[2]) +
    prop_tru[5] * plogis(bs[2] + bp + bd4[2]) +
    prop_tru[6] * plogis(bs[3] + bp + bd4[2])
  
  p_choice_rif_1 <-  
    prop_tru[1] * plogis(bs[1] + bd4[3]) +
    prop_tru[2] * plogis(bs[2] + bd4[3]) +
    prop_tru[3] * plogis(bs[3] + bd4[3]) +
    prop_tru[4] * plogis(bs[1] + bp + bd4[3]) +
    prop_tru[5] * plogis(bs[2] + bp + bd4[3]) +
    prop_tru[6] * plogis(bs[3] + bp + bd4[3])
  
  # second model
  # y ~ -1 + s + d2 + d4  for d1 = 2 (all silo)
  # average over bd2 all pref at zero since is relevant only for rev(1)
  # combinations
  # bs[1] + d2[1] + bd4[norif]
  # bs[2] + d2[1] + bd4[norif]
  # bs[3] + d2[1] + bd4[norif]
  # bs[1] + d2[2] + bd4[norif]
  # bs[2] + d2[2] + bd4[norif]
  # bs[3] + d2[2] + bd4[norif]
  # bs[1] + d2[3] + bd4[norif]
  # bs[2] + d2[3] + bd4[norif]
  # bs[3] + d2[3] + bd4[norif]
  
  prop_tru <- numeric(9)
  prop_tru[1] <- p_s_alloc[1] * (1-l_e$p_d2_entry) 
  prop_tru[2] <- p_s_alloc[2] * (1-l_l$p_d2_entry)
  prop_tru[3] <- p_s_alloc[3] * (1-l_c$p_d2_entry)
  # technically the l_e$p_d2_alloc should be 1 - l_e$p_d2_alloc but we have
  # 1:1 randomisation so it doesn't matter.
  prop_tru[4] <- p_s_alloc[1] * (l_e$p_d2_entry * l_e$p_d2_alloc)
  prop_tru[5] <- p_s_alloc[2] * (l_l$p_d2_entry * l_l$p_d2_alloc)
  prop_tru[6] <- p_s_alloc[3] * (l_c$p_d2_entry * l_c$p_d2_alloc)
  prop_tru[7] <- p_s_alloc[1] * (l_e$p_d2_entry * l_e$p_d2_alloc)
  prop_tru[8] <- p_s_alloc[2] * (l_l$p_d2_entry * l_l$p_d2_alloc)
  prop_tru[9] <- p_s_alloc[3] * (l_c$p_d2_entry * l_c$p_d2_alloc)
  
  p_choice_norif_2 <-  
    prop_tru[1] * plogis(bs[1] + bd1[2] + bd2[1] + bd4[2]) +
    prop_tru[2] * plogis(bs[2] + bd1[2] + bd2[1] + bd4[2]) +
    prop_tru[3] * plogis(bs[3] + bd1[2] + bd2[1] + bd4[2]) +
    prop_tru[4] * plogis(bs[1] + bd1[2] + bd2[2] + bd4[2]) +
    prop_tru[5] * plogis(bs[2] + bd1[2] + bd2[2] + bd4[2]) +
    prop_tru[6] * plogis(bs[3] + bd1[2] + bd2[2] + bd4[2]) +
    prop_tru[7] * plogis(bs[1] + bd1[2] + bd2[3] + bd4[2]) +
    prop_tru[8] * plogis(bs[2] + bd1[2] + bd2[3] + bd4[2]) +
    prop_tru[9] * plogis(bs[3] + bd1[2] + bd2[3] + bd4[2]) 
  
  p_choice_rif_2 <-  
    prop_tru[1] * plogis(bs[1] + bd1[2] + bd2[1] + bd4[3]) +
    prop_tru[2] * plogis(bs[2] + bd1[2] + bd2[1] + bd4[3]) +
    prop_tru[3] * plogis(bs[3] + bd1[2] + bd2[1] + bd4[3]) +
    prop_tru[4] * plogis(bs[1] + bd1[2] + bd2[2] + bd4[3]) +
    prop_tru[5] * plogis(bs[2] + bd1[2] + bd2[2] + bd4[3]) +
    prop_tru[6] * plogis(bs[3] + bd1[2] + bd2[2] + bd4[3]) +
    prop_tru[7] * plogis(bs[1] + bd1[2] + bd2[3] + bd4[3]) +
    prop_tru[8] * plogis(bs[2] + bd1[2] + bd2[3] + bd4[3]) +
    prop_tru[9] * plogis(bs[3] + bd1[2] + bd2[3] + bd4[3]) 
  
  
  # third model
  
  # y ~ -1 + s + d3 + d4  for d1 = 3 (all silo)
  # combinations
  # bs[1] + bp + d3[1] + bd4[norif]
  # bs[2] + bp + d3[1] + bd4[norif]
  # bs[3] + bp + d3[1] + bd4[norif]
  # bs[1] + bp + d3[2] + bd4[norif]
  # bs[2] + bp + d3[2] + bd4[norif]
  # bs[3] + bp + d3[2] + bd4[norif]
  # bs[1] + bp + d3[3] + bd4[norif]
  # bs[2] + bp + d3[3] + bd4[norif]
  # bs[3] + bp + d3[3] + bd4[norif]
  
  prop_tru <- numeric(9)
  # pr of being in silo * pr of being in non-rand set
  prop_tru[1] <- p_s_alloc[1] * (1-l_e$p_d3_entry) 
  prop_tru[2] <- p_s_alloc[2] * (1-l_l$p_d3_entry)
  prop_tru[3] <- p_s_alloc[3] * (1-l_c$p_d3_entry)
  # technically the l_e$p_d2_alloc should be 1 - l_e$p_d2_alloc but we have
  # 1:1 randomisation so it doesn't matter.
  prop_tru[4] <- p_s_alloc[1] * (l_e$p_d3_entry * l_e$p_d3_alloc)
  prop_tru[5] <- p_s_alloc[2] * (l_l$p_d3_entry * l_l$p_d3_alloc)
  prop_tru[6] <- p_s_alloc[3] * (l_c$p_d3_entry * l_c$p_d3_alloc)
  prop_tru[7] <- p_s_alloc[1] * (l_e$p_d3_entry * l_e$p_d3_alloc)
  prop_tru[8] <- p_s_alloc[2] * (l_l$p_d3_entry * l_l$p_d3_alloc)
  prop_tru[9] <- p_s_alloc[3] * (l_c$p_d3_entry * l_c$p_d3_alloc)
  
  p_choice_norif_3 <-  
    prop_tru[1] * plogis(bs[1] + bp + bd1[3] + bd3[1] + bd4[2]) +
    prop_tru[2] * plogis(bs[2] + bp + bd1[3] + bd3[1] + bd4[2]) +
    prop_tru[3] * plogis(bs[3] + bp + bd1[3] + bd3[1] + bd4[2]) +
    prop_tru[4] * plogis(bs[1] + bp + bd1[3] + bd3[2] + bd4[2]) +
    prop_tru[5] * plogis(bs[2] + bp + bd1[3] + bd3[2] + bd4[2]) +
    prop_tru[6] * plogis(bs[3] + bp + bd1[3] + bd3[2] + bd4[2]) +
    prop_tru[7] * plogis(bs[1] + bp + bd1[3] + bd3[3] + bd4[2]) +
    prop_tru[8] * plogis(bs[2] + bp + bd1[3] + bd3[3] + bd4[2]) +
    prop_tru[9] * plogis(bs[3] + bp + bd1[3] + bd3[3] + bd4[2]) 
  
  p_choice_rif_3 <-  
    prop_tru[1] * plogis(bs[1] + bp + bd1[3] + bd3[1] + bd4[3]) +
    prop_tru[2] * plogis(bs[2] + bp + bd1[3] + bd3[1] + bd4[3]) +
    prop_tru[3] * plogis(bs[3] + bp + bd1[3] + bd3[1] + bd4[3]) +
    prop_tru[4] * plogis(bs[1] + bp + bd1[3] + bd3[2] + bd4[3]) +
    prop_tru[5] * plogis(bs[2] + bp + bd1[3] + bd3[2] + bd4[3]) +
    prop_tru[6] * plogis(bs[3] + bp + bd1[3] + bd3[2] + bd4[3]) +
    prop_tru[7] * plogis(bs[1] + bp + bd1[3] + bd3[3] + bd4[3]) +
    prop_tru[8] * plogis(bs[2] + bp + bd1[3] + bd3[3] + bd4[3]) +
    prop_tru[9] * plogis(bs[3] + bp + bd1[3] + bd3[3] + bd4[3]) 
  
  # model estimates need to be weighted by the relevant cohort contributions
  
  p_choice_norif <- 
    # those allocated to dair (contribute to model 1)
    ((p_s_alloc[1] * (1-l_e$p_d1_alloc)) + 
       (p_s_alloc[2] * (1-l_l$p_d1_alloc)) + 
       (p_s_alloc[3] * (1-l_c$p_d1_alloc))) * p_choice_norif_1 +
    # those allocated to rev(1) (contribute to model 2)
    ((p_s_alloc[1] * (l_e$p_d1_alloc * (1-l_e$p_pref))) + 
       (p_s_alloc[2] * (l_l$p_d1_alloc * (1-l_l$p_pref))) + 
       (p_s_alloc[3] * (l_c$p_d1_alloc * (1-l_c$p_pref)))) * p_choice_norif_2 +
    # those allocated to rev(2) (contribute to model 3)
    ((p_s_alloc[1] * (l_e$p_d1_alloc * l_e$p_pref)) + 
       (p_s_alloc[2] * (l_l$p_d1_alloc * l_l$p_pref)) + 
       (p_s_alloc[3] * (l_c$p_d1_alloc * l_c$p_pref))) * p_choice_norif_3 
    
  p_choice_rif <- 
    # those allocated to dair (contribute to model 1)
    ((p_s_alloc[1] * (1-l_e$p_d1_alloc)) + 
       (p_s_alloc[2] * (1-l_l$p_d1_alloc)) + 
       (p_s_alloc[3] * (1-l_c$p_d1_alloc))) * p_choice_rif_1 +
    # those allocated to rev(1) (contribute to model 2)
    ((p_s_alloc[1] * (l_e$p_d1_alloc * (1-l_e$p_pref))) + 
       (p_s_alloc[2] * (l_l$p_d1_alloc * (1-l_l$p_pref))) + 
       (p_s_alloc[3] * (l_c$p_d1_alloc * (1-l_c$p_pref)))) * p_choice_rif_2 +
    # those allocated to rev(2) (contribute to model 3)
    ((p_s_alloc[1] * (l_e$p_d1_alloc * l_e$p_pref)) + 
       (p_s_alloc[2] * (l_l$p_d1_alloc * l_l$p_pref)) + 
       (p_s_alloc[3] * (l_c$p_d1_alloc * l_c$p_pref))) * p_choice_rif_3 

  
  c(rd_choice = p_choice_rif - p_choice_norif,
    p_choice_rif = p_choice_rif, 
    p_choice_norif = p_choice_norif,
    p_choice_norif_1 = p_choice_norif_1, 
    p_choice_norif_2 = p_choice_norif_2,
    p_choice_norif_3 = p_choice_norif_3, 
    p_choice_rif_1 = p_choice_rif_1,
    p_choice_rif_2 = p_choice_rif_2,
    p_choice_rif_3 = p_choice_rif_3)

}


trial_data <- function(l_spec){
  
  d <- data.table(
    s = sample(1:3, l_spec$N, replace = T, prob = l_spec$p_s_alloc)
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
    pref = rbinom(.N, 1, l_spec$l_e$p_pref)
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
    pref = rbinom(.N, 1, l_spec$l_l$p_pref)
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
    pref = rbinom(.N, 1, l_spec$l_c$p_pref)
  )]
  
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
  
  # d[, .N, keyby = .(s, pref, d1, d2, d3, d4)]
  
  # bd1 is irrelevant since d1 == 1 is the ref group, fixed at zero
  d[d1 == 1, eta := l_spec$bs[s] + l_spec$bp*pref + l_spec$bd4[d4]]
  # pref is irrelevant as d1 = 2 only occurs if pref = 0
  d[d1 == 2, eta := l_spec$bs[s] + l_spec$bd1[2] + l_spec$bd2[d2] + l_spec$bd4[d4]]
  # but here pref is relevant as d1 = 3 only if pref = 1
  d[d1 == 3, eta := l_spec$bs[s] + l_spec$bp*pref + l_spec$bd1[3] + l_spec$bd3[d3] + l_spec$bd4[d4]]
  
  d[, `:=`(s = factor(s), d1 = factor(d1), d2 = factor(d2), d3 = factor(d3), d4 = factor(d4))]
  
  d[, p := plogis(eta)]
  d[, y := rbinom(.N, 1, p)]
  
  d
}


multi_model_approach <- function(l_spec){
  
  r <- pbapply::pblapply(X=1:l_spec$N_sim, cl = l_spec$mc_cores, FUN = function(ix){
    
    d <- trial_data(l_spec)
    
    f1_1 <- glm(
      y ~ -1 + s + pref + d4 , data = d, subset= d1==1, family = binomial,
      control = glm.control(epsilon = 1e-4, maxit = 100, trace = FALSE))
    # pref is zero for everyone
    f1_2 <- glm(
      y ~ -1 + s + d2 + d4, data = d, subset= d1==2, family = binomial,
      control = glm.control(epsilon = 1e-4, maxit = 100, trace = FALSE))
    # effect of pref will get rolled into the intercept
    f1_3 <- glm(
      y ~ -1 + s + d3 + d4, data = d, subset= d1==3, family = binomial,
      control = glm.control(epsilon = 1e-4, maxit = 100, trace = FALSE))
    
    # coef(f1_1)
    # coef(f1_2)
    # coef(f1_3)
    
    # surgical domain -----
    
    # dair
    d_surg_dair <- copy(d[s == 2])
    d_surg_dair[, d1 := factor(1)]
    # d4 stays at whatever it was
    d_surg_dair[, p_hat := predict(f1_1, newdata = d_surg_dair, type = "response")]
    p_surg_dair_hat <- mean(d_surg_dair$p_hat)
    
    # rev(1)
    d_surg_rev_1 <- copy(d[s == 2 & pref == 0])
    d_surg_rev_1[, d1 := factor(2)]
    d_surg_rev_1[, pref := 0]
    d_surg_rev_1[d2_entry == 1, d2 := factor(d_surg_rev_1[d2_entry == 1, d2_alloc + 2], levels = 1:3)]
    d_surg_rev_1[d2_entry == 0, d2 := factor(1, levels = 1:3)]
    # and d4 stays at whatever it was
    d_surg_rev_1[, p_hat := predict(f1_2, newdata = d_surg_rev_1, type = "response")]
    p_surg_rev_1_hat <- mean(d_surg_rev_1$p_hat)
    
    # rev(2)
    d_surg_rev_2 <- copy(d[s == 2 & pref == 1])
    d_surg_rev_2[, d1 := factor(3)]
    # pref doesn't matter because the model had to roll it up into the intercept
    # d_rev_2[, pref := 1]
    d_surg_rev_2[d3_entry == 1, d3 := factor(d_surg_rev_2[d3_entry == 1, d3_alloc + 2], levels = 1:3)]
    d_surg_rev_2[d3_entry == 0, d3 := factor(1, levels = 1:3)]
    # and d4 stays at whatever it was
    d_surg_rev_2[, p_hat := predict(f1_3, newdata = d_surg_rev_2, type = "response")]
    p_surg_rev_2_hat <- mean(d_surg_rev_2$p_hat)
    
    # now aggregate the one and two-stage revision
    d_prop <- d[s == 2, .(prop = .N/nrow(d[s == 2])), keyby = pref]
    p_surg_rev_hat <- p_surg_rev_1_hat * d_prop[pref == 0, prop] + p_surg_rev_2_hat * d_prop[pref == 1, prop]
    
    rd_surg_hat <- p_surg_rev_hat - p_surg_dair_hat
    
    # don't really need to do this for all, but:
    risk_surg <- risk_pars_surg(
      l_spec$l_l$p_d1_alloc, 
      l_spec$l_l$p_d2_entry, l_spec$l_l$p_d2_alloc,
      l_spec$l_l$p_d3_entry, l_spec$l_l$p_d3_alloc,
      l_spec$l_l$p_d4_entry, l_spec$l_l$p_d4_alloc, 
      l_spec$l_l$p_pref, 
      l_spec$bs, l_spec$bp, l_spec$bd1, l_spec$bd2, l_spec$bd3, l_spec$bd4)
    
    
    # duration domain -----
    
    d_dur_12wk <- copy(d[d1 == 2])
    d_dur_12wk[, d2 := factor(2)]
    d_dur_12wk[, p_hat := predict(f1_2, newdata = d_dur_12wk, type = "response")]
    p_dur_12wk_hat <- mean(d_dur_12wk$p_hat)
    
    d_dur_6wk <- copy(d[d1 == 2])
    d_dur_6wk[, d2 := factor(3)]
    d_dur_6wk[, p_hat := predict(f1_2, newdata = d_dur_6wk, type = "response")]
    p_dur_6wk_hat <- mean(d_dur_6wk$p_hat)
    
    rd_dur_hat <- p_dur_6wk_hat - p_dur_12wk_hat
    
    risk_dur <- risk_pars_dur(
      l_spec$p_s_alloc, l_spec$l_e, l_spec$l_l, l_spec$l_c, 
      l_spec$bs, l_spec$bp, l_spec$bd1, l_spec$bd2, l_spec$bd3, l_spec$bd4)
    
    # extp domain -----
    
    d_extp_0wk <- copy(d[d1 == 3])
    d_extp_0wk[, d3 := factor(2)]
    d_extp_0wk[, p_hat := predict(f1_3, newdata = d_extp_0wk, type = "response")]
    p_extp_0wk_hat <- mean(d_extp_0wk$p_hat)
    
    d_extp_12wk <- copy(d[d1 == 3])
    d_extp_12wk[, d3 := factor(3)]
    d_extp_12wk[, p_hat := predict(f1_3, newdata = d_extp_12wk, type = "response")]
    p_extp_12wk_hat <- mean(d_extp_12wk$p_hat)
    
    rd_extp_hat <- p_extp_12wk_hat - p_extp_0wk_hat
    
    risk_extp <- risk_pars_extp(
      l_spec$p_s_alloc, l_spec$l_e, l_spec$l_l, l_spec$l_c, 
      l_spec$bs, l_spec$bp, l_spec$bd1, l_spec$bd2, l_spec$bd3, l_spec$bd4)
    
    # choice domain -----
    
    d_choice_norif <- copy(d[d4 != 1 & d1 == 1])
    d_choice_norif[, d4 := factor(2)]
    d_choice_norif[, p_hat := predict(f1_1, newdata = d_choice_norif, type = "response")]
    p_choice_norif_1_hat <- mean(d_choice_norif$p_hat)
    
    d_choice_norif <- copy(d[d4 != 1 & d1 == 2])
    d_choice_norif[, d4 := factor(2)]
    d_choice_norif[, p_hat := predict(f1_2, newdata = d_choice_norif, type = "response")]
    p_choice_norif_2_hat <- mean(d_choice_norif$p_hat)
    
    d_choice_norif <- copy(d[d4 != 1 & d1 == 3])
    d_choice_norif[, d4 := factor(2)]
    d_choice_norif[, p_hat := predict(f1_3, newdata = d_choice_norif, type = "response")]
    p_choice_norif_3_hat <- mean(d_choice_norif$p_hat)
    
    d_choice_rif <- copy(d[d4 != 1 & d1 == 1])
    d_choice_rif[, d4 := factor(3)]
    d_choice_rif[, p_hat := predict(f1_1, newdata = d_choice_rif, type = "response")]
    p_choice_rif_1_hat <- mean(d_choice_rif$p_hat)
    
    d_choice_rif <- copy(d[d4 != 1 & d1 == 2])
    d_choice_rif[, d4 := factor(3)]
    d_choice_rif[, p_hat := predict(f1_2, newdata = d_choice_rif, type = "response")]
    p_choice_rif_2_hat <- mean(d_choice_rif$p_hat)
    
    d_choice_rif <- copy(d[d4 != 1 & d1 == 3])
    d_choice_rif[, d4 := factor(3)]
    d_choice_rif[, p_hat := predict(f1_3, newdata = d_choice_rif, type = "response")]
    p_choice_rif_3_hat <- mean(d_choice_rif$p_hat)
    
    
    d_prop_s <- d[, .(prop = .N/nrow(d)), keyby = s]
    d_prop_e_d1 <- d[s == 1, .(prop = .N/nrow(d[s == 1])), keyby = d1]
    d_prop_l_d1 <- d[s == 2, .(prop = .N/nrow(d[s == 2])), keyby = d1]
    d_prop_c_d1 <- d[s == 3, .(prop = .N/nrow(d[s == 3])), keyby = d1]
    
    d_prop_e_pref <- d[s == 1 , .(prop = .N/nrow(d[s == 1])), keyby = pref]
    d_prop_l_pref <- d[s == 2 , .(prop = .N/nrow(d[s == 2])), keyby = pref]
    d_prop_c_pref <- d[s == 3 , .(prop = .N/nrow(d[s == 3])), keyby = pref]
    
    p_choice_norif_hat <- 
      # those allocated to dair (contribute to model 1)
      ((d_prop_s[s==1, prop] * (d_prop_e_d1[d1 == 1, prop])) + 
         (d_prop_s[s==2, prop] * (d_prop_l_d1[d1 == 1, prop])) + 
         (d_prop_s[s==3, prop] * (d_prop_c_d1[d1 == 1, prop]))) * p_choice_norif_1_hat +
      # those allocated to rev(1) (contribute to model 2)
      ((d_prop_s[s==1, prop] * (d_prop_e_d1[d1 == 2, prop])) + 
         (d_prop_s[s==2, prop] * (d_prop_l_d1[d1 == 2, prop])) + 
         (d_prop_s[s==3, prop] * (d_prop_c_d1[d1 == 2, prop]))) * p_choice_norif_2_hat +
      # those allocated to rev(2) (contribute to model 3)
      ((d_prop_s[s==1, prop] * (d_prop_e_d1[d1 == 3, prop])) + 
         (d_prop_s[s==2, prop] * (d_prop_l_d1[d1 == 3, prop])) + 
         (d_prop_s[s==3, prop] * (d_prop_c_d1[d1 == 3, prop]))) * p_choice_norif_3_hat 
    
    p_choice_rif_hat <- 
      # those allocated to dair (contribute to model 1)
      ((d_prop_s[s==1, prop] * (d_prop_e_d1[d1 == 1, prop])) + 
         (d_prop_s[s==2, prop] * (d_prop_l_d1[d1 == 1, prop])) + 
         (d_prop_s[s==3, prop] * (d_prop_c_d1[d1 == 1, prop]))) * p_choice_rif_1_hat +
      # those allocated to rev(1) (contribute to model 2)
      ((d_prop_s[s==1, prop] * (d_prop_e_d1[d1 == 2, prop])) + 
         (d_prop_s[s==2, prop] * (d_prop_l_d1[d1 == 2, prop])) + 
         (d_prop_s[s==3, prop] * (d_prop_c_d1[d1 == 2, prop]))) * p_choice_rif_2_hat +
      # those allocated to rev(2) (contribute to model 3)
      ((d_prop_s[s==1, prop] * (d_prop_e_d1[d1 == 3, prop])) + 
         (d_prop_s[s==2, prop] * (d_prop_l_d1[d1 == 3, prop])) + 
         (d_prop_s[s==3, prop] * (d_prop_c_d1[d1 == 3, prop]))) * p_choice_rif_3_hat 
    
    rd_choice_hat <- p_choice_rif_hat - p_choice_norif_hat
    
    risk_choice_ref <- risk_pars_choice(
      l_spec$p_s_alloc, l_spec$l_e, l_spec$l_l, l_spec$l_c, 
      l_spec$bs, l_spec$bp, l_spec$bd1, l_spec$bd2, l_spec$bd3, l_spec$bd4)
    

    data.table(
      p_surg_dair_ref = risk_surg["p_dair"],
      p_surg_dair_hat = p_surg_dair_hat,
      p_surg_rev_1_ref = risk_surg["p_rev_1"],
      p_surg_rev_1_hat = p_surg_rev_1_hat,
      p_surg_rev_2_ref = risk_surg["p_rev_2"],
      p_surg_rev_2_hat = p_surg_rev_2_hat,
      p_surg_rev_ref = risk_surg["p_rev"],
      p_surg_rev_hat = p_surg_rev_hat,
      rd_surg_ref = risk_surg["rd"],
      rd_surg_hat = rd_surg_hat,
      
      p_dur_12wk_ref = risk_dur["p_dur_12wk"],
      p_dur_12wk_hat = p_dur_12wk_hat,
      p_dur_6wk_ref = risk_dur["p_dur_6wk"],
      p_dur_6wk_hat = p_dur_6wk_hat,
      rd_dur_ref = risk_dur["rd_dur"],
      rd_dur_hat = rd_dur_hat,
      
      p_extp_0wk_ref = risk_extp["p_extp_0wk"],
      p_extp_0wk_hat = p_extp_0wk_hat,
      p_extp_12wk_ref = risk_extp["p_extp_12wk"],
      p_extp_12wk_hat = p_extp_12wk_hat,
      rd_extp_ref = risk_extp["rd_extp"],
      rd_extp_hat = rd_extp_hat,
      
      p_choice_norif_ref = risk_choice_ref["p_choice_norif"],
      p_choice_norif_hat = p_choice_norif_hat,
      p_choice_rif_ref = risk_choice_ref["p_choice_rif"],
      p_choice_rif_hat = p_choice_rif_hat,
      rd_choice_ref = risk_choice_ref["rd_choice"],
      rd_choice_hat = rd_choice_hat
    )
    
    
  })
  
  r
  
  
  
}




single_model_approach <- function(l_spec){

  r <- pbapply::pblapply(X=1:l_spec$N_sim, cl = l_spec$mc_cores, FUN = function(ix){
    
    d <- trial_data(l_spec)
    
    # arbitrary, levels don't get used but have to exist for the sake of the
    # initial setup.
    d[is.na(d2), d2 := factor(4)]
    d[is.na(d3), d3 := factor(4)]
    
    f2 <- glm(
      y ~ -1 + s + pref + d1 + d2:(d1==2) + d3:(d1==3) + d4, 
      data = d, 
      family = binomial,
      control = glm.control(epsilon = 1e-4, maxit = 100, trace = FALSE))
    # will give you a lot of NA in summary since over parameterised but it gets us to the 
    # predictions we are targeting
    # summary(f2)
    
    
    # surgical domain ------
    
    # dair
    d_surg_dair <- copy(d[s == 2])
    d_surg_dair[, d1 := factor(1)]
    d_surg_dair[, d2 := factor(4)]
    d_surg_dair[, d3 := factor(4)]
    # d4 stays at whatever it was
    d_surg_dair[, p_hat := predict(f2, newdata = d_surg_dair, type = "response")]
    p_surg_dair_hat <- mean(d_surg_dair$p_hat)
    
    # rev(1)
    d_surg_rev_1 <- copy(d[s == 2 & pref == 0])
    d_surg_rev_1[, d1 := factor(2)]
    d_surg_rev_1[, pref := 0]
    d_surg_dair[, d3 := factor(4)]
    d_surg_rev_1[d2_entry == 1, d2 := factor(d_surg_rev_1[d2_entry == 1, d2_alloc + 2], levels = 1:3)]
    d_surg_rev_1[d2_entry == 0, d2 := factor(1, levels = 1:3)]
    # and d4 stays at whatever it was
    d_surg_rev_1[, p_hat := predict(f2, newdata = d_surg_rev_1, type = "response")]
    p_surg_rev_1_hat <- mean(d_surg_rev_1$p_hat)
    
    # rev(2)
    d_surg_rev_2 <- copy(d[s == 2 & pref == 1])
    d_surg_rev_2[, d1 := factor(3)]
    d_surg_dair[, d2 := factor(4)]
    # pref doesn't matter because the model had to roll it up into the intercept
    # d_rev_2[, pref := 1]
    d_surg_rev_2[d3_entry == 1, d3 := factor(d_surg_rev_2[d3_entry == 1, d3_alloc + 2], levels = 1:3)]
    d_surg_rev_2[d3_entry == 0, d3 := factor(1, levels = 1:3)]
    # and d4 stays at whatever it was
    d_surg_rev_2[, p_hat := predict(f2, newdata = d_surg_rev_2, type = "response")]
    p_surg_rev_2_hat <- mean(d_surg_rev_2$p_hat)
    
    # now aggregate the one and two-stage revision
    d_prop <- d[s == 2, .(prop = .N/nrow(d[s == 2])), keyby = pref]
    p_surg_rev_hat <- p_surg_rev_1_hat * d_prop[pref == 0, prop] + p_surg_rev_2_hat * d_prop[pref == 1, prop]
    
    rd_surg_hat <- p_surg_rev_hat - p_surg_dair_hat
    
    # don't really need to do this for all the sims, but:
    risk_surg <- risk_pars_surg(
      l_spec$l_l$p_d1_alloc, 
      l_spec$l_l$p_d2_entry, l_spec$l_l$p_d2_alloc,
      l_spec$l_l$p_d3_entry, l_spec$l_l$p_d3_alloc,
      l_spec$l_l$p_d4_entry, l_spec$l_l$p_d4_alloc, 
      l_spec$l_l$p_pref, 
      l_spec$bs, l_spec$bp, l_spec$bd1, l_spec$bd2, l_spec$bd3, l_spec$bd4)
    
    
    # duration domain -----
    
    d_dur_12wk <- copy(d[d1 == 2])
    d_dur_12wk[, d2 := factor(2)]
    # has to be set as the hold out set
    d_dur_12wk[, d3 := factor(4)]
    d_dur_12wk[, p_hat := predict(f2, newdata = d_dur_12wk, type = "response")]
    p_dur_12wk_hat <- mean(d_dur_12wk$p_hat)
    
    d_dur_6wk <- copy(d[d1 == 2])
    d_dur_6wk[, d2 := factor(3)]
    # has to be set as the hold out set
    d_dur_12wk[, d3 := factor(4)]
    d_dur_6wk[, p_hat := predict(f2, newdata = d_dur_6wk, type = "response")]
    p_dur_6wk_hat <- mean(d_dur_6wk$p_hat)
    
    rd_dur_hat <- p_dur_6wk_hat - p_dur_12wk_hat
    
    risk_dur <- risk_pars_dur(
      l_spec$p_s_alloc, l_spec$l_e, l_spec$l_l, l_spec$l_c, 
      l_spec$bs, l_spec$bp, l_spec$bd1, l_spec$bd2, l_spec$bd3, l_spec$bd4)
    
    # extp domain -----
    
    d_extp_0wk <- copy(d[d1 == 3])
    d_extp_0wk[, d3 := factor(2)]
    # has to be set as the hold out set
    d_dur_12wk[, d2 := factor(4)]
    d_extp_0wk[, p_hat := predict(f2, newdata = d_extp_0wk, type = "response")]
    p_extp_0wk_hat <- mean(d_extp_0wk$p_hat)
    
    d_extp_12wk <- copy(d[d1 == 3])
    d_extp_12wk[, d3 := factor(3)]
    # has to be set as the hold out set
    d_dur_12wk[, d2 := factor(4)]
    d_extp_12wk[, p_hat := predict(f2, newdata = d_extp_12wk, type = "response")]
    p_extp_12wk_hat <- mean(d_extp_12wk$p_hat)
    
    rd_extp_hat <- p_extp_12wk_hat - p_extp_0wk_hat
    
    risk_extp <- risk_pars_extp(
      l_spec$p_s_alloc, l_spec$l_e, l_spec$l_l, l_spec$l_c, 
      l_spec$bs, l_spec$bp, l_spec$bd1, l_spec$bd2, l_spec$bd3, l_spec$bd4)
    
    
    # choice domain -----
    
    # not sure how to do this without doing it within d1 levels
    # if you 
    d_choice_norif <- copy(d[d4 != 1 & d1 == 1])
    d_choice_norif[, d4 := factor(2)]
    # hold out conditional on dair
    d_choice_norif[, d2 := factor(4)]
    d_choice_norif[, d3 := factor(4)]
    d_choice_norif[, p_hat := predict(f2, newdata = d_choice_norif, type = "response")]
    p_choice_norif_1_hat <- mean(d_choice_norif$p_hat)
    
    d_choice_norif <- copy(d[d4 != 1 & d1 == 2])
    d_choice_norif[, d4 := factor(2)]
    # hold out conditional on rev(1)
    d_choice_norif[, d3 := factor(4)]
    d_choice_norif[, p_hat := predict(f2, newdata = d_choice_norif, type = "response")]
    p_choice_norif_2_hat <- mean(d_choice_norif$p_hat)
    
    d_choice_norif <- copy(d[d4 != 1 & d1 == 3])
    d_choice_norif[, d4 := factor(2)]
    # hold out conditional on rev(1)
    d_choice_norif[, d2 := factor(4)]
    d_choice_norif[, p_hat := predict(f2, newdata = d_choice_norif, type = "response")]
    p_choice_norif_3_hat <- mean(d_choice_norif$p_hat)
    
    d_choice_rif <- copy(d[d4 != 1 & d1 == 1])
    d_choice_rif[, d4 := factor(3)]
    d_choice_norif[, d2 := factor(4)]
    d_choice_norif[, d3 := factor(4)]
    d_choice_rif[, p_hat := predict(f2, newdata = d_choice_rif, type = "response")]
    p_choice_rif_1_hat <- mean(d_choice_rif$p_hat)
    
    d_choice_rif <- copy(d[d4 != 1 & d1 == 2])
    d_choice_rif[, d4 := factor(3)]
    d_choice_norif[, d3 := factor(4)]
    d_choice_rif[, p_hat := predict(f2, newdata = d_choice_rif, type = "response")]
    p_choice_rif_2_hat <- mean(d_choice_rif$p_hat)
    
    d_choice_rif <- copy(d[d4 != 1 & d1 == 3])
    d_choice_rif[, d4 := factor(3)]
    d_choice_norif[, d2 := factor(4)]
    d_choice_rif[, p_hat := predict(f2, newdata = d_choice_rif, type = "response")]
    p_choice_rif_3_hat <- mean(d_choice_rif$p_hat)
    
    
    d_prop_s <- d[, .(prop = .N/nrow(d)), keyby = s]
    d_prop_e_d1 <- d[s == 1, .(prop = .N/nrow(d[s == 1])), keyby = d1]
    d_prop_l_d1 <- d[s == 2, .(prop = .N/nrow(d[s == 2])), keyby = d1]
    d_prop_c_d1 <- d[s == 3, .(prop = .N/nrow(d[s == 3])), keyby = d1]
    
    d_prop_e_pref <- d[s == 1 , .(prop = .N/nrow(d[s == 1])), keyby = pref]
    d_prop_l_pref <- d[s == 2 , .(prop = .N/nrow(d[s == 2])), keyby = pref]
    d_prop_c_pref <- d[s == 3 , .(prop = .N/nrow(d[s == 3])), keyby = pref]
    
    p_choice_norif_hat <- 
      # those allocated to dair (contribute to model 1)
      ((d_prop_s[s==1, prop] * (d_prop_e_d1[d1 == 1, prop])) + 
         (d_prop_s[s==2, prop] * (d_prop_l_d1[d1 == 1, prop])) + 
         (d_prop_s[s==3, prop] * (d_prop_c_d1[d1 == 1, prop]))) * p_choice_norif_1_hat +
      # those allocated to rev(1) (contribute to model 2)
      ((d_prop_s[s==1, prop] * (d_prop_e_d1[d1 == 2, prop])) + 
         (d_prop_s[s==2, prop] * (d_prop_l_d1[d1 == 2, prop])) + 
         (d_prop_s[s==3, prop] * (d_prop_c_d1[d1 == 2, prop]))) * p_choice_norif_2_hat +
      # those allocated to rev(2) (contribute to model 3)
      ((d_prop_s[s==1, prop] * (d_prop_e_d1[d1 == 3, prop])) + 
         (d_prop_s[s==2, prop] * (d_prop_l_d1[d1 == 3, prop])) + 
         (d_prop_s[s==3, prop] * (d_prop_c_d1[d1 == 3, prop]))) * p_choice_norif_3_hat 
    
    p_choice_rif_hat <- 
      # those allocated to dair (contribute to model 1)
      ((d_prop_s[s==1, prop] * (d_prop_e_d1[d1 == 1, prop])) + 
         (d_prop_s[s==2, prop] * (d_prop_l_d1[d1 == 1, prop])) + 
         (d_prop_s[s==3, prop] * (d_prop_c_d1[d1 == 1, prop]))) * p_choice_rif_1_hat +
      # those allocated to rev(1) (contribute to model 2)
      ((d_prop_s[s==1, prop] * (d_prop_e_d1[d1 == 2, prop])) + 
         (d_prop_s[s==2, prop] * (d_prop_l_d1[d1 == 2, prop])) + 
         (d_prop_s[s==3, prop] * (d_prop_c_d1[d1 == 2, prop]))) * p_choice_rif_2_hat +
      # those allocated to rev(2) (contribute to model 3)
      ((d_prop_s[s==1, prop] * (d_prop_e_d1[d1 == 3, prop])) + 
         (d_prop_s[s==2, prop] * (d_prop_l_d1[d1 == 3, prop])) + 
         (d_prop_s[s==3, prop] * (d_prop_c_d1[d1 == 3, prop]))) * p_choice_rif_3_hat 
    
    rd_choice_hat <- p_choice_rif_hat - p_choice_norif_hat
    
    risk_choice_ref <- risk_pars_choice(
      l_spec$p_s_alloc, l_spec$l_e, l_spec$l_l, l_spec$l_c, 
      l_spec$bs, l_spec$bp, l_spec$bd1, l_spec$bd2, l_spec$bd3, l_spec$bd4)
    
    
    data.table(
      p_surg_dair_ref = risk_surg["p_dair"],
      p_surg_dair_hat = p_surg_dair_hat,
      p_surg_rev_1_ref = risk_surg["p_rev_1"],
      p_surg_rev_1_hat = p_surg_rev_1_hat,
      p_surg_rev_2_ref = risk_surg["p_rev_2"],
      p_surg_rev_2_hat = p_surg_rev_2_hat,
      p_surg_rev_ref = risk_surg["p_rev"],
      p_surg_rev_hat = p_surg_rev_hat,
      rd_surg_ref = risk_surg["rd"],
      rd_surg_hat = rd_surg_hat,
      
      p_dur_12wk_ref = risk_dur["p_dur_12wk"],
      p_dur_12wk_hat = p_dur_12wk_hat,
      p_dur_6wk_ref = risk_dur["p_dur_6wk"],
      p_dur_6wk_hat = p_dur_6wk_hat,
      rd_dur_ref = risk_dur["rd_dur"],
      rd_dur_hat = rd_dur_hat,
      
      p_extp_0wk_ref = risk_extp["p_extp_0wk"],
      p_extp_0wk_hat = p_extp_0wk_hat,
      p_extp_12wk_ref = risk_extp["p_extp_12wk"],
      p_extp_12wk_hat = p_extp_12wk_hat,
      rd_extp_ref = risk_extp["rd_extp"],
      rd_extp_hat = rd_extp_hat,
      
      p_choice_norif_ref = risk_choice_ref["p_choice_norif"],
      p_choice_norif_hat = p_choice_norif_hat,
      p_choice_rif_ref = risk_choice_ref["p_choice_rif"],
      p_choice_rif_hat = p_choice_rif_hat,
      rd_choice_ref = risk_choice_ref["rd_choice"],
      rd_choice_hat = rd_choice_hat
    )
    
    
  })
  
  r
  
  
  
}


parse_results <- function(r){
  
  # results -----
  d_g <- rbindlist(r)
  
  d_g[, bias_surg_rd := rd_surg_hat - rd_surg_ref]
  d_g[, bias_dur_rd := rd_dur_hat - rd_dur_ref]
  d_g[, bias_extp_rd := rd_extp_hat - rd_extp_ref]
  d_g[, bias_choice_rd := rd_choice_hat - rd_choice_ref]
  
  d_fig <- d_g[, .(bias_surg_rd, bias_dur_rd, bias_extp_rd, bias_choice_rd)]
  d_fig_1 <- melt(d_fig, measure.vars = names(d_fig))
  
  
  d_fig <- d_g[, .(p_surg_dair_hat, p_surg_rev_hat,
                   p_dur_12wk_hat, p_dur_6wk_hat,
                   p_extp_0wk_hat, p_extp_12wk_hat,
                   p_choice_norif_hat, p_choice_rif_hat)]
  d_fig_2 <- melt(d_fig, measure.vars = names(d_fig))
  
  d_fig <- d_g[, .(rd_surg_hat, rd_dur_hat, rd_extp_hat, rd_choice_hat)]
  d_fig_3 <- melt(d_fig, measure.vars = names(d_fig))
  
  
  list(
    # bias
    d_fig_1 = d_fig_1,
    # risk
    d_fig_2 = d_fig_2, 
    # risk difference
    d_fig_3 = d_fig_3
  )
}

plot_results <- function(r){
  
  l_f <- parse_results(r)
  
  p1 <- ggplot(l_f$d_fig_1, aes(x = value)) +
    geom_density() +
    geom_vline(data = l_f$d_fig_1[, .(mu = mean(value)), keyby = variable],
               aes(xintercept = mu), lwd = 0.2) +
    facet_wrap(~variable, nrow = 1) +
    theme(
      axis.title = element_blank()
    )
  
  p2 <- ggplot(l_f$d_fig_2, aes(x = value)) +
    geom_density() +
    geom_vline(data = l_f$d_fig_2[, .(mu = mean(value)), keyby = variable],
               aes(xintercept = mu), lwd = 0.2) +
    facet_wrap(~variable,  nrow = 4)+
    theme(
      axis.title.x = element_blank()
    )
  
  p3 <- ggplot(l_f$d_fig_3, aes(x = value)) +
    geom_density() +
    geom_vline(data = l_f$d_fig_3[, .(mu = mean(value)), keyby = variable],
               aes(xintercept = mu), lwd = 0.2)  +
    facet_wrap(~variable,  nrow = 2)+
    theme(
      axis.title.x = element_blank()
    )
  
  layout <- "
    AAAA
    BBCC
    BBCC
  "
  print(p1 + p2 + p3 +
    plot_layout(design = layout))
  
}


examples <- function(){
  
  # setup for simulations, data generation and model parameters
  l_spec <- list(
    
    # number of sims and cores to use
    N_sim = 1000,
    mc_cores = 4,
    
    
    # sample size
    N = 3e3,
    # silo allocation
    p_s_alloc = c(0.3, 0.5, 0.2),
    
    # early silo
    l_e = list(
      # prob of rev in surg domain (assuming all enter surg)
      p_d1_alloc = 0.15,
      # probability of entering d2 assuming you are in an eligible set
      # note that entry into duration same across all silos
      p_d2_entry = 0.7,
      # 1:1 randomisation (note that if you deviate from this then some of the
      # risk calcs with need updating so that the right weight is applied to 
      # each of the randomised arm contributions)
      p_d2_alloc = 0.5,
      # probability of entering d3 assuming you are in an eligible set
      # entry into extp same across all silos
      p_d3_entry = 0.9,
      p_d3_alloc = 0.5,
      # probability of entering d4
      # entry into choice same across all silos
      p_d4_entry = 0.6,
      p_d4_alloc = 0.5,
      # preference for two-stage
      p_pref = 0.35
    ),
    l_l = list(
      # prob of rev in surg domain (assuming all enter surg)
      p_d1_alloc = 0.5,
      p_d2_entry = 0.7,
      p_d2_alloc = 0.5,
      p_d3_entry = 0.9,
      p_d3_alloc = 0.5,
      p_d4_entry = 0.6,
      p_d4_alloc = 0.5,
      # preference for two-stage
      p_pref = 0.7
    ),
    l_c = list(
      # prob of rev in surg domain (assuming all enter surg)
      p_d1_alloc = 0.8,
      p_d2_entry = 0.7,
      p_d2_alloc = 0.5,
      p_d3_entry = 0.9,
      p_d3_alloc = 0.5,
      p_d4_entry = 0.6,
      p_d4_alloc = 0.5,
      # preference for two-stage
      p_pref = 0.75
    ),
    
    # model parameters
    
    # log-odds (early, late, acute)
    bs = c(0.9, 0.8, 0.7),
    
    # log-or
    # different baseline risk for rev
    bp = -0.4,
    bd1 = c(0, 0, 0) ,
    
    # index 4 is never referenced but needs to be there so that the 
    # calcs in the single model approach don't end up with NA
    bd2 = c(0, 0, 0, 999),
    bd3 = c(0, 0, 0, 999),
    bd4 = c(0, 0, 0)
    
  )
  
  
  
  
  # objective function to target effect on relative log scale corresponding
  # to a zero effect on the absolute risk scale.
  # tweak as necessary to target arbitrary effect size, or whatever variation 
  # is of interest.
  obj_surg_f1 <- function(x, l_spec){
    
    l_spec$bd1 = c(0, x, x)
    
    res <- risk_pars_surg(
      l_spec$l_l$p_d1_alloc, 
      l_spec$l_l$p_d2_entry, l_spec$l_l$p_d2_alloc,
      l_spec$l_l$p_d3_entry, l_spec$l_l$p_d3_alloc,
      l_spec$l_l$p_d4_entry, l_spec$l_l$p_d4_alloc, 
      l_spec$l_l$p_pref, 
      l_spec$bs, l_spec$bp, 
      l_spec$bd1, 
      l_spec$bd2, 
      l_spec$bd3, 
      l_spec$bd4
    )
    # square error
    (res["rd"] - 0)^2
  }
  
  # use to work out the null scenario in terms of risk when using 
  # log odds ratios to set up the simulation
  
  # to get a null scenario on the risk scale when both rev(1) and rev(2) 
  # are non-zero, you need to set the logor for both to:
  optimize(obj_surg_f1, c(-1, 1), l_spec)$minimum
  
  
  
  # examples -------
  
  # null effects in all domain on log odds scale translates to null on abs 
  # pr scale in all domain
  single_model_approach(l_spec) |>
    plot_results()
  multi_model_approach(l_spec) |>
    plot_results()
  
  
  # introduce abx duration effect 
  l_spec$bd1 <- c(0, 0, 0)
  l_spec$bd2 <- c(0, 0, 0.9, 999)
  l_spec$bd3 <- c(0, 0, 0, 999)
  l_spec$bd4 <- c(0, 0, 0)
  
  # unbiased but reveals small +ve risk difference due to nonlinear transform
  # see bottom right set of plots
  single_model_approach(l_spec) |>
    plot_results()
  # should be same (seems to run much quicker)
  multi_model_approach(l_spec) |>
    plot_results()
  
  # update so that we resolve to zero on risk scale:
  # adjust surg domain in order to represent null effect on 
  # risk scale
  (lor_d1_rev <- optimize(obj_surg_f1, c(-1, 1), l_spec)$minimum)
  l_spec$bd1 <- c(0, lor_d1_rev, lor_d1_rev)
  # rd_surg_hat should now be centred on zero
  single_model_approach(l_spec) |>
    plot_results()
  # same
  multi_model_approach(l_spec) |>
    plot_results()
  
  # occurs if effects in d3 alone
  l_spec$bd1 <- c(0, 0, 0)
  l_spec$bd2 <- c(0, 0, 0, 999)
  l_spec$bd3 <- c(0, 0, -0.5, 999)
  l_spec$bd4 <- c(0, 0, 0)
  single_model_approach(l_spec) |>
    plot_results()
  (lor_d1_rev <- optimize(obj_surg_f1, c(-1, 1), l_spec)$minimum)
  l_spec$bd1 <- c(0, lor_d1_rev, lor_d1_rev)
  # notice a bias here for ext proph but haven;t determined why 
  single_model_approach(l_spec) |>
    plot_results()
  # same here 
  multi_model_approach(l_spec) |>
    plot_results()
  # bias not due to insufficient sims 
  # will need to check on this...
  l_spec$N_sim <- 5000
  multi_model_approach(l_spec) |>
    plot_results()
  l_spec$N_sim <- 1000
  
  
  # but not if effects in d4 alone, i.e. no adjustment required to resolve to 
  # the null in the surgical domain on the risk scale.
  l_spec$bd1 <- c(0, 0, 0)
  l_spec$bd2 <- c(0, 0, 0, 999)
  l_spec$bd3 <- c(0, 0, 0, 999)
  l_spec$bd4 <- c(0, 0, 1.2)
  single_model_approach(l_spec) |>
    plot_results()
  # surg transform is already zero
  (lor_d1_rev <- optimize(obj_surg_f1, c(-1, 1), l_spec)$minimum)
  
}



