# Platform trial in Prosthetic joint infection management
#
# Primary outcome is treatment successs (binary yes/no) defined as being alive
# with clinical cure, no use of antibiotics and prosthetic in place.
#
# The primary outcome is always 12 months post "Platform entry" irrespective of
# the combination of domains a pt enters into. Platform entry is the timepoint
# when the patient has met core eligibility criteria, given informed consent for
# the platform, and been randomised across all domains. Reveal to allocation
# within specific domains occurs later, aligned with the eligibility timing for
# the domain.
#
# In all of the following, the following definitions apply:
#
# D1 = surgical management domain. D2 = antibiotic duration domain following the
# first stage procedure. D3 = extended prophylaxis domain following the second
# stage procedure. D4 = adjunctive Antibiotic type domain.
#
# DAIR=Debridement, Antibiotics and Implant Retention
#
# 1s = One stage revision is so called because there is a single surgery
# procedure where the infected prosthetic is removed and replaced. 
# 2s = A two stage procedure that implants an antibiotic infused temporary
# prosthetic at the first surgery and then replaces this with the final
# prosthetic in the second surgery.
#
# D1 randomised treatment options are dair or revision. If assigned to revision,
# the clinician still gets to decide whether a 1s or 2s occurs.
#
# D2 randomised treatment options are 12 weeks and 6 weeks antibiotic following
# the first surgery. Receipt of D2 trts is to occur after the surgery but no
# more than 7 days after surgery and no later than 12 weeks after platform
# entry. Only patients having had 1s revision can enter into d2. For 2s pts, the
# antibiotic duration received after the first surgery will not be randomised.
#
# D3 randomised treatment options are 12 weeks Extended Prophylaxis and no
# extended prophylaxis (none). Receipt is after second surgery but no later than
# 6 months after platform entry date. Only patients having had 2s revision can
# enter into d2. The antibiotic duration that was received after the first
# surgery will not have been randomised.
#
# D4 randomised treatment options are Backbone regimen with or without
# adjunctive rifampicin. Reveal (if applicable) occurs bewteen 48h and 7 days
# following surgery. Pts need to have a relevant pathogen profile but there is 
# no dependency on the type of surgery.
# 
# Silo E = Early stage Prosthetic joint infection (Early pt <=30 days post
# implant)
#
# Silo L = Late acute PJI (Late acute pt = >30 days post implant and <= 21 days
# of symptoms at diagnosis.)
#
# Silo C = Chronic PJI (Chronic pt =>30 days post implant AND a sinus or >21
# days of symptoms)
#
# Setup for each silo:
#
# Silo E patients: 
#
# D1 no randomisation option. The clinician decides whether the pt gets dair, 1s
# or 2s.
#
# D2 randomisation may be revealed for pts having had 1s.
#
# D3 randomisation may be revealed for pts having had 2s.
#
# D4 randomisation may be revealed for any pt in this silo having the relevant
# pathogen profile
#
# Silo L:
#
# D1 randomisation as above to dair or revision. The clinician decides whether
# the pt gets dair, 1s or 2s.
#
# D2 randomisation may be revealed for pts having had 1s.
#
# D3 randomisation may be revealed for pts having had 2s.
#
# D4 randomisation may be revealed for any pt in this silo having the relevant
# pathogen profile
#
# Silo L:
#
# D1 no randomisation option. The clinician decides whether the pt gets dair, 1s
# or 2s.
#
# D2 randomisation may be revealed for pts having had 1s.
#
# D3 randomisation may be revealed for pts having had 2s.
#
# D4 randomisation may be revealed for any pt in this silo having the relevant
# pathogen profile
#
# NOTE:
#
# It is conceivable that a pt enters into one or more domain for rand trt but
# not the others. However, They may still receive some treatment option that was
# not randomised but is relevant to one of the domains. For example, a pt might
# get a two stage revision, does not enter d3 but does get 12 wks of ext
# prophylaxis. This might arise because the pt had the trt outside of the time
# window defined by the protocol, or perhaps a clinical decision or some other
# reason. The pt still remains part of the analysis population for the other
# domains where they did enter into randomised treatment, e.g. the antibiotic
# type domain.
#



library(data.table)
library(ggplot2)
library(fastglm)
library(parallel)


sim09_cohort <- function(l_spec){
  
  # need to incorporate joint. bar the obvious, not sure which way should do this...
  
  id_cohort <- l_spec$is:l_spec$ie
  N_cohort <- length(id_cohort)
  
  d <- data.table(id = id_cohort, t_0 = l_spec$t_0[id_cohort])
  
  # baseline / design covariates
  d[, silo     := sample(l_spec$lab_silo, .N, replace = TRUE, prob = l_spec$pr_silo)]
  # standardised prognostic score (frailty/infection severity)
  d[, sev := rnorm(.N)]  
  # start building linear predictor
  d[, lp_base := l_spec$b_0 + l_spec$b_sev * sev]
  
  # D1: surgical
  # late silo: randomised dair vs revision; early/chronic: clinician-directed
  # assume all late silo are randomised (best case scenario for this trt effect)
  d[, rand_d1 := as.integer(silo == "l")]
  d[, d1_trt := NA_character_]
  # assign trt to those in late acute
  i_ix <- which(d$silo == "l")
  n_ix <- length(i_ix)
  if (n_ix > 0L) {
    d1_rev <- runif(n_ix) < l_spec$pr_d1_rand
    d$d1_trt[i_ix] <- fifelse(
      d1_rev,
      "revision",
      "dair"
    )
  }
  
  # compute probability of two stage based on sev
  # b_sev_stage is a parameter related to decision process for 1s or 2s rev
  # higher sev leads to increased chance of deciding on 2s
  d[, d1_arm := NA_character_]
  # at this stage all those with trt = rev are in the late silo
  i_ix <- which(d$d1_trt == "revision")
  n_ix <- length(i_ix)
  if (n_ix > 0L) {
    stage_p <- plogis(l_spec$b_sev_stage * d$sev[i_ix])
    d$d1_arm[i_ix] <- fifelse(
      runif(n_ix) < stage_p,
      "2s",
      "1s"
    )
  }
  # just set this for completeness
  i_ix <- which(d$d1_trt == "dair")
  d$d1_arm[i_ix] <- "dair"
  # still got to fill in dair
  
  # non-randomised assignment is assumed to be based on severity (sev is centred
  # on zero so for avg sev, we have multinom with p_k = 0.32, 0.42, 0.26.
  # positive sev increases prob mass for 2 stage. could configure these parameters
  # but they are pretty arbitrary anyway so not really much point. it is more
  # to just have a formal mechanism for assignment in place.
  i_ix <- which(d$rand_d1 == 0L)
  n_ix <- length(i_ix)
  if (n_ix > 0L) {
    sev_nr <- d$sev[i_ix]
    u_dair = 0
    u_1s   = 0.3 - 0.4 * sev_nr
    u_2s   = -0.2 + 0.9 * sev_nr
    
    e1 <- exp(u_1s)
    e2 <- exp(u_2s)
    denom <- 1 + e1 + e2
    
    p_dair <- 1 / denom
    p_1s   <- e1 / denom
    
    r <- runif(n_ix)
    d$d1_arm[i_ix] <- fifelse(r < p_dair, "dair",
                              fifelse(r < p_dair + p_1s, "1s", "2s")
                              )
     
  }
  # domain 1 linear predictor increment
  d[, lp_d1_inc := l_spec$b_d1[d1_arm]]
  
  # D2: antibiotic durati on (one-stage patients only, any silo and both rand/non-rand)
  # assume all pt for which this applicable will reach this reveal coz only a 
  # few wks after surgery
  d[, reveal_d2 := NA_integer_]
  d[, d2_trt := NA_character_]
  d[, lp_d2_inc := 0.0]
  
  # almost certain reveal in absence of non-zero b_reveal_d2 parameter
  # in practice, I think it is probably a lot more complicated that this as we
  # have a situation where the reveal process is not easily represented.
  # this is just an approximation there are also aspects like a pt or clin 
  # preference/suitability wrt d2
  i_ix <- which(d$d1_arm == "1s")
  n_ix <- length(i_ix)
  
  if (n_ix > 0L) {
    # reveal decision
    p_reveal_d2 <- plogis(8 - l_spec$b_reveal_d2 * d$sev[i_ix])
    d$reveal_d2[i_ix] <- as.integer(runif(n_ix) < p_reveal_d2)
  }
  
  # Randomise those reaching D2
  i_ix <- which(d$reveal_d2 == 1L); n_ix <- length(i_ix)
  if (n_ix > 0L) {
    d$d2_trt[i_ix] <- fifelse(runif(n_ix) < l_spec$pr_d2_rand, "12wk", "6wk" )
    d$lp_d2_inc[i_ix] <- l_spec$b_d2[d$d2_trt[i_ix]]
  }
  
  # D3: ext proh (two-stage patients, conditioned on surviving to the second-stage 
  # reveal) higher sev means less chance of reveal
  d[, reveal_d3 := NA_integer_]
  d[, d3_trt := NA_character_]
  d[, lp_d3_inc := 0.0]
  
  i_ix <- which(d$d1_arm == "2s"); n_ix <- length(i_ix)
  if (n_ix > 0L) {
    # Survival/reveal
    p_reveal_d3 <- plogis(8 - l_spec$b_reveal_d3 * d$sev[i_ix])
    d$reveal_d3[i_ix] <- as.integer(runif(n_ix) < p_reveal_d3)
  }
  
  i_ix <- which(d$reveal_d3 == 1L); n_ix <- length(i_ix)
  if (n_ix > 0L) {
    d$d3_trt[i_ix] <- fifelse(runif(n_ix) < l_spec$pr_d3_rand, "12wk", "none")
    d$lp_d3_inc[i_ix] <- l_spec$b_d3[d$d3_trt[i_ix]]
  }
  
  
  # D4: rif (condit on pathogen-profile eligibility) assume reveal is basically
  # immediate and there is no survival consideration
  d[, d4_elig := as.integer(runif(.N) < l_spec$pr_d4_elig)]
  
  d[, d4_trt := NA_character_]
  d[, lp_d4_inc := 0.0]
  
  i_ix <- which(d$d4_elig == 1L); n_ix <- length(i_ix)
  if (n_ix > 0L) {
    d$d4_trt[i_ix] <- fifelse(runif(n_ix) < l_spec$pr_d4_rand, "rif", "norif")
    d$lp_d4_inc[i_ix] <- l_spec$b_d4[d$d4_trt[i_ix]]
  }
  
  d[, lp := lp_base + lp_d1_inc + lp_d2_inc + lp_d3_inc + lp_d4_inc
      # fcoalesce(lp_d2_inc, 0) +
      # fcoalesce(lp_d3_inc, 0) +
      # fcoalesce(lp_d4_inc, 0) 
    ]
  
  d[, y := rbinom(.N, 1, plogis(lp))]
  
  
  
  d[]
  
}

sim09_cfg <- function(N_tot = 2000){
  
  l_spec <- structure(
    list(
      is = 1,
      ie = N_tot,
      b_d2 = c(`12wk` = 0, `6wk` = 0), 
      b_d3 = c(none = 0, `12wk` = 0), 
      b_d4 = c(norif = 0, rif = 0), 
      desc = "RD = 0 in all domains", 
      n_sim = 500L, mc_cores = 40L, nex = 3L, 
      N_pt = N_tot, 
      pr_silo = c(e = 0.3, l = 0.5, c = 0.2), 
      lab_silo = c("e", "l", "c"), 
      pr_knee = 0.55, 
      pr_d1_rand = 0.5, 
      pr_d4_elig = 0.65, 
      pr_d2_rand = 0.5, 
      pr_d3_rand = 0.5, 
      pr_d4_rand = 0.5, 
      b_sev_stage = 0.2, 
      b_reveal_d2 = 0, 
      b_reveal_d3 = 0, 
      b_0 = 0.2, 
      b_sev = 0, 
      b_knee = c(0, 0.8), 
      b_d1 = c(dair = 0, `1s` = 0.2, `2s` = 0.3), 
      lab_b_d1_a = c("dair", "revision"), 
      lab_b_d1_b = c("dair", "1s", "2s"), 
      lab_b_d2 = c("12wk", "6wk"), 
      lab_b_d3 = c("none", "12wk"), 
      lab_b_d4 = c("norif", "rif"), 
      b_d2_x_d4 = 0, 
      b_d3_x_d4 = 0, 
      ex_trial_ix = c(1, 4, 9))
  )
  
  l_spec
  
}

sim09_ex_fit_d1 <- function(l_spec){
  
  d_sim <- rbindlist(parallel::mclapply(
    X=1:l_spec$n_sim, mc.cores = 4, FUN=function(ii) {
      
      d <- sim09_cohort(l_spec)
      d_mod <- d[rand_d1 == 1L, ]
      d_mod[, silo := factor(silo, levels = l_spec$lab_silo)]
      
      d_mod[, d1_trt := factor(d1_trt, levels = l_spec$lab_b_d1_a)]
      d_mod[, d1_arm := factor(d1_arm, levels = l_spec$lab_b_d1_b)]
      
      d_mod[, d2_trt := factor(d2_trt, levels = c(l_spec$lab_b_d2))]
      d_mod[, d3_trt := factor(d3_trt, levels = c(l_spec$lab_b_d3))]
      
      d_mod[, regimen := NA_character_]
      d_mod[d1_arm == "dair", regimen := "dair"]
      d_mod[d1_arm == "2s", regimen := paste0("2s:", d3_trt)]
      d_mod[d1_arm == "1s", regimen := paste0("1s:", d2_trt)]
      d_mod[, regimen := factor(regimen, levels = c(
        "dair", 
        "1s:12wk",
        "1s:6wk",
        "1s:NA",
        "2s:none",
        "2s:12wk",
        "2s:NA"
      ))]
      
      d_reg <- d_mod[regimen != "dair", .(N_reg = .N), keyby = regimen]
      d_reg[, w := N_reg / sum(d_reg$N_reg)]
      
      X <- model.matrix(~ sev + regimen, data = d_mod)
      f_d1 <- fastglm(x = X, y = d_mod$y, family = binomial)
      
      d_mod <- d[d1_arm == "1s" & reveal_d2 == 1, ]
      d_mod[, silo := factor(silo, levels = l_spec$lab_silo)]
      d_mod[, trt := factor(d2_trt, levels = l_spec$lab_b_d2)]

      X <- model.matrix(~ trt + sev + silo, data = d_mod)
      f_d2 <- fastglm(x = X, d_mod$y, family = binomial)

      d_mod <- d[d1_arm == "2s" & reveal_d3 == 1, ]
      d_mod[, silo := factor(silo, levels = l_spec$lab_silo)]
      d_mod[, trt := factor(d3_trt, levels = l_spec$lab_b_d3)]

      X <- model.matrix(~ trt + sev + silo, data = d_mod)
      f_d3 <- fastglm(x = X, d_mod$y, family = binomial)

      d_mod <- d[d4_elig == 1, ]
      d_mod[, silo := factor(silo, levels = l_spec$lab_silo)]
      d_mod[, trt := factor(d4_trt, levels = l_spec$lab_b_d4)]

      X <- model.matrix(~ trt + sev + silo, data = d_mod)
      f_d4 <- fastglm(x = X, d_mod$y, family = binomial)

      cf_d1 <- summary(f_d1)$coef[, "Estimate"]
      cf_d2 <- summary(f_d2)$coef[, "Estimate"]
      cf_d3 <- summary(f_d3)$coef[, "Estimate"]
      cf_d4 <- summary(f_d4)$coef[, "Estimate"]

      d_out <- data.table(
        domain = c(
          rep("d1", length(cf_d1)),
          rep("d2", length(cf_d2)),
          rep("d3", length(cf_d3)),
          rep("d4", length(cf_d4))
        ),
        par = names(c(cf_d1, cf_d2, cf_d3, cf_d4)),
        est = c(cf_d1, cf_d2, cf_d3, cf_d4)
      )
      d_out[, ix := 1:.N]

      d_reg[, par := paste0("regimen", regimen)]
      d_out <- base::merge(
        d_out, d_reg[, .(par, w)], by = "par", all.x = T
      )
      d_out[is.na(w), w := 1]
      setorder(d_out, "ix")

      d_out
    }
  ), idcol = "id_sim")
  
  
  
}

sim09_ex_d1 <- function(){
  
  set.seed(1)
  l_spec <- sim09_cfg()
  # override b_reveal_d2 param
  l_spec$b_reveal_d2 <- 0
  d_res_1 <- sim09_ex_fit_d1(l_spec)
  
  # still ok
  l_spec$b_reveal_d2 <- 6
  d <- sim09_cohort(l_spec)
  d_res_2 <- sim09_ex_fit_d1(l_spec)
  
  # cannot recover d1 effect for arm1s as the effect of d2 is 
  # incorporated into the arm1s estimate
  # the arm1s is less likely to achive trt success due to the possibility
  # of receiving a detrimental treatment in d2
  l_spec$b_d2[2] <- -0.4
  d <- sim09_cohort(l_spec)
  d_res_3 <- sim09_ex_fit_d1(l_spec)
  
  d_res <- base::merge(d_res_1, d_res_2, by = c("domain", "par"))
  d_res <- base::merge(d_res, d_res_3, by = c("domain", "par"))
  colnames(d_res) <- c("domain", "par", "result 1", "result 2", "result 3")
  
  kableExtra::kable(
    d_res, digits = 2, format = "simple"
  )
  
  # The domain specific model for d1 gives the effect of the surgical strategy 
  # together with whatever downstream treatment opportunities the surgical  
  # intervention creates.
  
  
}


