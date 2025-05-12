# Experiment with independent set of models with reduced linear predictor.


source("./R/init.R")
source("./R/data.R")

# Command line arguments list scenario (true dose response),
# the number of simulations to run, the number of cores
# to use and the simulator to use.
args = commandArgs(trailingOnly=TRUE)

# Load cfg based on cmd line args.
if (length(args)<1) {
  log_info("Setting default run method (does nothing)")
  args[1] = "run_none_sim_05"
  args[2] = "./sim05/cfg-sim05-sc01-v01.yml"
} else {
  log_info("Run method ", args[1])
  log_info("Scenario config ", args[2])
}

# Log setup
f_log_sim <- file.path("./logs", "log-sim.txt")
log_appender(appender_file(f_log_sim))
log_info("*** START UP ***")
# log_threshold(TRACE)

f_cfgsc <- file.path("./etc", args[2])
g_cfgsc <- config::get(file = f_cfgsc)
stopifnot("Config is null" = !is.null(g_cfgsc))

ix <- 1
m1 <- cmdstanr::cmdstan_model("stan/model-sim-05-a.stan")
m2 <- cmdstanr::cmdstan_model("stan/model-sim-05-b.stan")

output_dir_mcmc <- paste0(getwd(), "/tmp")

# Sanity check:
# The target enrolment rate is 1.5 participants per day, i.e. 
# about 10 per 7 day week.
# n_per_wk = numeric(100)
# for(i in 1:100){
#   loc_t0 <- get_enrol_time_int(max(g_cfgsc$N_pt), lambda, rho)
#   # ignore the first 400 as these are the ramp up phase.
#   d_fig <- data.table(
#     i = seq_along(loc_t0),
#     t0 = loc_t0
#   )
#   # d_fig[i %in% 400:420]
#   wks_total <- max(loc_t0) / 7
#   n_per_wk[i] <- max(g_cfgsc$N_pt) / wks_total
# }
# n_per_wk


# Main trial loop.
run_trial <- function(
    ix,
    mu = 0,
    # silo effects
    # silo 1 as the one for late acute
    # silo 2 as the one for late acute
    # silo 3 is LATE ACUTE
    b_silo = c(0, -0.3, -0.2),
    b_jnt = c(0, 0.4),
    b_pref = c(0, -0.2),
    # dair, one, two-stage
    # b_d1 = c(0, -0.1, -0.2),
    b_d1 = c(0, 0.693, 0.693, 0, 0, 0, 0, 0.693, 0.693),
    # 
    b_d2 = c(0, -0.4, 0.6),
    b_d3 = c(0, 0.1, 0.3),
    b_d4 = c(0, -0.1, -0.2),
    l_prior = list(
      pri_mu = c(0, 0.47),
      pri_b_silo = c(0, 1),
      pri_b_jnt = c(0, 1),
      pri_b_prf = c(0, 1),
      pri_b_trt = c(0, 1)
    ),
    return_posterior = F
    ){
  
  log_info("Entered  run_trial for trial ", ix)
  
  # Get enrolment times for arbitrary large set of patients
  # Simpler to produce enrol times all in one hit rather than try to do them 
  # incrementally with a non-hom pp.
  
  # events per day
  lambda = 1.52
  # ramp up over 12 months 
  rho = function(t) pmin(t/360, 1)
  
  loc_t0 <- get_enrol_time_int(max(g_cfgsc$N_pt), lambda, rho)
  
  # loop controls
  stop_enrol <- FALSE
  ii <- 1 # interim number
  N_analys <- length(g_cfgsc$N_pt)

  # posterior summaries
  d_post_smry_1 <- CJ(
    id_analys = 1:N_analys,
    domain = 1:4,
    par = c("p0", "p1", "rd", "lor")
  )
  d_post_smry_1[, mu := NA_real_]
  d_post_smry_1[, med := NA_real_]
  d_post_smry_1[, se := NA_real_]
  d_post_smry_1[, q_025 := NA_real_]
  d_post_smry_1[, q_975 := NA_real_]
  
  d_post_smry_2 <- data.table()

  # superiority probs
  pr_sup <- array(
    NA, 
    # num analysis x num domains
    dim = c(N_analys, 4),
    dimnames = list(1:N_analys, paste0("d", 1:4))
  )
  # probability of futility wrt superiority decision
  pr_sup_fut  <- array(
    NA, 
    # num analysis x num domains
    dim = c(N_analys, 4),
    dimnames = list(1:N_analys, paste0("d", 1:4))
  )
  # non-inferiority probs - trt is ni to reference (soc)
  pr_ni <- array(
    NA, 
    # num analysis x num domains
    dim = c(N_analys, 4),
    dimnames = list(1:N_analys, paste0("d", 1:4))
  )
  # non-inferiority probs
  pr_ni_fut <- array(
    NA, 
    # num analysis x num domains
    dim = c(N_analys, 4),
    dimnames = list(1:N_analys, paste0("d", 1:4))
  )
  
  # units informing estimates
  n_units <- array(
    NA, 
    # num analysis x num domains
    dim = c(N_analys, 4),
    dimnames = list(1:N_analys, paste0("d", 1:4))
  )
 
  # decisions 
  # superior, ni, 
  # futile (for superiority - idiotic)
  # futile (for ni - idiotic)
  g_dec_type <- c("sup", 
                  "ni", 
                  "fut_sup", 
                  "fut_ni"
                  )

  decision <- array(
    NA,
    # num analysis x num domains x num decision types
    dim = c(N_analys, 4, length(g_dec_type)),
    dimnames = list(
      1:N_analys, paste0("d", 1:4), g_dec_type)
  )

  # store all simulated trial pt data
  d_all <- data.table()
  
  # decisions - only include domains for which they are evaluated
  dec_sup = list(
    surg = NA,
    ext_proph = NA,
    ab_choice = NA
  )
  dec_ni = list(
    ab_dur = NA
  )
  dec_sup_fut = list(
    surg = NA,
    ext_proph = NA,
    ab_choice = NA
  )
  dec_ni_fut = list(
    ab_dur = NA
  )
  
  if(return_posterior){
    l_post <- list()
  }
  
  
  ## LOOP -------
  while(!stop_enrol){
    
    log_info("Trial ", ix, " analysis ", ii)
  
    # next chunk of data on pts.
    if(ii == 1){
      N_c <- g_cfgsc$N_pt[ii]
      # starting pt index in data
      is <- 1
      ie <- is + N_c - 1
    } else {
      N_c <- g_cfgsc$N_pt[ii] - g_cfgsc$N_pt[ii-1]
      is <- g_cfgsc$N_pt[ii-1] + 1
      ie <- is + N_c - 1
    }
    
    N = N_c
    
    # id and time
    idx_s = is
    t0 = loc_t0[is:ie]
    id_analys = NULL
    
    # Our analyses only occur on those that have reached 12 months post 
    # randomisation. As such, we are assuming that the analysis takes place
    # 12 months following the last person to be enrolled in the current 
    # analysis set.
    l_new <- get_trial_data_int(
      N = N_c, 
      
      # true parameters
      mu = mu,
      b_silo = b_silo,
      b_jnt = b_jnt,
      b_pref = b_pref,
      b_d1 = b_d1,
      b_d2 = b_d2,
      b_d3 = b_d3,
      b_d4 = b_d4,
      
      dec_sup = dec_sup,
      dec_ni = dec_ni,
      dec_sup_fut = dec_sup_fut,
      dec_ni_fut = dec_ni_fut,
      
      idx_s = is,
      t0 = loc_t0[is:ie],
      id_analys = ii
    )
    
    log_info("Trial ", ix, " new data generated ", ii)

    # combine the existing and new data
    d_all <- rbind(d_all, l_new$d)
    
    # create stan data format based on the relevant subsets of pt
    lsd <- get_stan_data_all_int(d_all)
    
    lsd$ld$pri_mu <- l_prior$pri_mu
    lsd$ld$pri_b_silo <- l_prior$pri_b_silo
    lsd$ld$pri_b_jnt <- l_prior$pri_b_jnt
    lsd$ld$pri_b_prf <- l_prior$pri_b_prf
    lsd$ld$pri_b_trt <- l_prior$pri_b_trt
    
    foutname <- paste0(
      format(Sys.time(), format = "%Y%m%d%H%M%S"),
      "-sim-", ix, "-intrm-", ii)
    
    # fit model - does it matter that I continue to fit the model after the
    # decision is made...?
    snk <- capture.output(
      f_1 <- m1$sample(
        lsd$ld, iter_warmup = 1000, iter_sampling = 2000,
        parallel_chains = 1, chains = 1, refresh = 0, show_exceptions = F,
        max_treedepth = 11,
        output_dir = output_dir_mcmc,
        output_basename = foutname
        )
    )
    
    log_info("Trial ", ix, " fitted models ", ii)

    # extract posterior - marginal probability of outcome by group
    # dair vs rev
    d_post <- data.table(f_1$draws(
      variables = c(
        
        "lor_d1", # log odds and lor
        "p_d1_1", "p_d1_23", "rd_d1",
        
        "lor_d2",
        "p_d2_2", "p_d2_3", "rd_d2",
        
        "lor_d3",
        "p_d3_2", "p_d3_3", "rd_d3",
        
        "lor_d4",
        "p_d4_2", "p_d4_3", "rd_d4"
        
      ),   # risk scale
      format = "matrix"))
    
    
    if(return_posterior){
      l_post[[ii]] <- copy(d_post)
    }
    
    d_post_long <- melt(d_post, measure.vars = names(d_post))
    d_post_long[variable %like% "d1", domain := 1]
    d_post_long[variable %like% "d2", domain := 2]
    d_post_long[variable %like% "d3", domain := 3]
    d_post_long[variable %like% "d4", domain := 4]
    
    # d_post_long[variable %like% "lor", .(mu = mean(value), sd = sd(value)), keyby = variable]
  
    d_post_smry_1[id_analys == ii, 
                  mu := d_post_long[, mean(value), 
                                    keyby = .(domain, variable)]$V1] 
    d_post_smry_1[id_analys == ii, 
                  med := d_post_long[, median(value), 
                                    keyby = .(domain, variable)]$V1]
    d_post_smry_1[id_analys == ii, 
                  se := d_post_long[, sd(value), 
                                     keyby = .(domain, variable)]$V1]
    d_post_smry_1[id_analys == ii, 
                  q_025 := d_post_long[, quantile(value, prob = 0.025), 
                                    keyby = .(domain, variable)]$V1]
    d_post_smry_1[id_analys == ii, 
                  q_975 := d_post_long[, quantile(value, prob = 0.975), 
                                    keyby = .(domain, variable)]$V1]

    # These should produce the same answers as the generated quantities
    # outputs.
    d_post_chk <- data.table(f_1$draws(
      variables = c(
        
        "bd1", "bd2", "bd3", "bd4",
        
        "nu_d1_1", "nu_d1_2", "nu_d1_3",
        "nu_d2_2", "nu_d2_3",
        "nu_d3_2", "nu_d3_3",
        "nu_d4_2", "nu_d4_3"
        
      ),   # risk scale
      format = "matrix"))
    d_post_chk <- melt(d_post_chk, measure.vars = names(d_post_chk))
    d_post_chk[variable %like% "d1", domain := 1]
    d_post_chk[variable %like% "d2", domain := 2]
    d_post_chk[variable %like% "d3", domain := 3]
    d_post_chk[variable %like% "d4", domain := 4]
    
    d_post_smry_2 <- rbind(
      d_post_smry_2,
      d_post_chk[, .(analys = ii,
                     mu = mean(value),
                     med = median(value),
                     se = sd(value),
                     q_025 = quantile(value, prob = 0.025),
                     q_975 = quantile(value, prob = 0.975)),
                 keyby = .(variable, domain)]
    )
    
    log_info("Trial ", ix, " extracted posterior ", ii)

    # this is how many units are used in the g-comp up to this analysis
    n_units[ii, ] <- c(
      sum(lsd$ld$n_d1),
      sum(lsd$ld$n_d2),
      sum(lsd$ld$n_d3),
      sum(lsd$ld$n_d4)
    )
    
    # superiority is implied by a high probability that the risk diff 
    # greater than zero
    pr_sup[ii, ]  <-   c(
      d_post[, mean(p_d1_23 - p_d1_1 > g_cfgsc$dec_ref$delta_sup)],
      d_post[, mean(p_d2_3 - p_d2_2 > g_cfgsc$dec_ref$delta_sup)],
      d_post[, mean(p_d3_3 - p_d3_2 > g_cfgsc$dec_ref$delta_sup)],
      d_post[, mean(p_d4_3 - p_d4_2 > g_cfgsc$dec_ref$delta_sup)]
    )
    # futility for the superiority decision is implied by a low probability 
    # that the risk diff is greater than some small value (2% difference)
    pr_sup_fut[ii, ]  <-   c(
      d_post[, mean(p_d1_23 - p_d1_1 > g_cfgsc$dec_ref$delta_sup_fut)],
      d_post[, mean(p_d2_3 - p_d2_2 > g_cfgsc$dec_ref$delta_sup_fut)],
      d_post[, mean(p_d3_3 - p_d3_2 > g_cfgsc$dec_ref$delta_sup_fut)],
      d_post[, mean(p_d4_3 - p_d4_2 > g_cfgsc$dec_ref$delta_sup_fut)]
    )
    
    # ni is implied by a high probability that the risk diff reduces the 
    # efficacy by no more than some small amount (here 2%)
    pr_ni[ii, ]  <-   c(
      d_post[, mean(p_d1_23 - p_d1_1 > g_cfgsc$dec_ref$delta_ni)],
      d_post[, mean(p_d2_3 - p_d2_2 > g_cfgsc$dec_ref$delta_ni)],
      d_post[, mean(p_d3_3 - p_d3_2 > g_cfgsc$dec_ref$delta_ni)],
      d_post[, mean(p_d4_3 - p_d4_2 > g_cfgsc$dec_ref$delta_ni)]
    )
    # futility for the ni decision is implied by a low probability that 
    # the risk diff suggests any efficacy   
    pr_ni_fut[ii, ]  <-   c(
      d_post[, mean(p_d1_23 - p_d1_1 > g_cfgsc$dec_ref$delta_ni_fut)],
      d_post[, mean(p_d2_3 - p_d2_2 > g_cfgsc$dec_ref$delta_ni_fut)],
      d_post[, mean(p_d3_3 - p_d3_2 > g_cfgsc$dec_ref$delta_ni_fut)],
      d_post[, mean(p_d4_3 - p_d4_2 > g_cfgsc$dec_ref$delta_ni_fut)]
    )

    log_info("Trial ", ix, " calculated decision quantities ", ii)
    
    decision[ii, , "sup"] <- pr_sup[ii, ] > g_cfgsc$dec_probs$thresh_sup
    decision[ii, , "ni"] <- pr_ni[ii, ] > g_cfgsc$dec_probs$thresh_ni
    
    # futility rules
    # taken to imply negligible chance of being superior or ni
    
    # negligible chance of being superior
    decision[ii, , "fut_sup"] <- pr_sup_fut[ii, ] < g_cfgsc$dec_probs$thresh_fut_sup
    # taken to imply negligible chance of being ni 
    decision[ii, , "fut_ni"] <- pr_ni_fut[ii, ] < g_cfgsc$dec_probs$thresh_fut_ni
    
    log_info("Trial ", ix, " compared to thresholds ", ii)
    
    # Earlier decisions are retained - once superiority has been decided, 
    # we retain this conclusion irrespective of subsequent post probabilities.
    # The following simply overwrites any decision reversal.
    # This means that there could be inconsistency with a silo pr_sup and the 
    # decision reported.
    decision[1:ii, , "sup"] <- apply(decision[1:ii, , "sup", drop = F], 2, function(z){ cumsum(z) > 0 })
    decision[1:ii, , "ni"] <- apply(decision[1:ii, , "ni", drop = F], 2, function(z){ cumsum(z) > 0 })
    decision[1:ii, , "fut_sup"] <- apply(decision[1:ii, , "fut_sup", drop = F], 2, function(z){ cumsum(z) > 0 })
    decision[1:ii, , "fut_ni"] <- apply(decision[1:ii, , "fut_ni", drop = F], 2, function(z){ cumsum(z) > 0 })
  
    # superiority decisions apply to domains 1, 3 and 4
    if(any(decision[ii, , "sup"])){
      # since there are only two treatments per cell, if a superiority decision 
      # is made then we have answered all the questions and we can stop 
      # enrolling into that cell. if there were more than two treatments then
      # we would need to take a different approach.
      
      # if we have found it to be superior to dair then all subsequent get 
      # revision. the value of the dec_sup$surg doesn't matter at this time,
      # it just needs to be something other than NA. the logic for selecting
      # the data is in the get_data function.
      # this decision would only impact late acute
      if(decision[ii, "d1", "sup"]){
        dec_sup$surg <- 3
      }
      # we are not assessing superiority for domain 2 (antibiotic duration)
      
      # extproph domain relates to all silo so decision on b_r2d impacts all cohorts
      if(decision[ii, "d3", "sup"]){
        dec_sup$ext_proph <- 3
      }
      # choice domain relates to all silo so decision on b_f impacts all cohorts
      if(decision[ii, "d4", "sup"]){
        dec_sup$ab_choice <- 3
      }
    }
    # stop enrolling if futile wrt superiority decision
    if(any(decision[ii, , "fut_sup"])){
      if(decision[ii, "d1", "fut_sup"]){
        dec_sup_fut$surg <- 3
      }
      if(decision[ii, "d3", "fut_sup"]){
        dec_sup_fut$ext_proph <- 3
      }
      if(decision[ii, "d4", "fut_sup"]){
        dec_sup_fut$ab_choice <- 3
      }
    }
    
    
    # ni decisions apply to domains 2
    # if short is ni to long then stop enrolment
    # for this randomised comparison
    if(any(decision[ii, , "ni"])){
      if(decision[ii, "d2", "ni"]){
        dec_ni$ab_dur <- 3
      }
    }
    
    if(any(decision[ii, , "fut_ni"])){
      if(decision[ii, "d2", "fut_ni"]){
        dec_ni_fut$ab_dur <- 3
      }
    }
    
    # have we answered all questions of interest?
    if(
      # if rev (in late acute silo) is superior to dair (or superiority decision
      # decision is futile to pursue)
      (decision[ii, "d1", "sup"] | decision[ii, "d1", "fut_sup"]) &
      # if 6wk backbone duration (in one-stage units) is non-inferior to 12wk
      # (or decision for non-inferiority is futile to pursue)
      (decision[ii, "d2", "ni"] | decision[ii, "d2", "fut_ni"] ) &
      # if 12wk ext-proph duration (in two-stage units) is sup to none
      # (or decision for sup is futile to pursue)
      (decision[ii, "d3", "sup"] | decision[ii, "d3", "fut_sup"]) &
      # if rif (in all relevant units) is sup to none 
      # (or decision for sup is futile to pursue)
      (decision[ii, "d4", "sup"] | decision[ii, "d4", "fut_sup"])
    ){
      log_info("Stop trial all questions addressed ", ix)
      stop_enrol <- T  
    } 
    
    log_info("Trial ", ix, " updated allocation control ", ii)
    
    # next interim
    ii <- ii + 1
    
    if(ii > N_analys){
      stop_enrol <- T  
    }
  }
  
  # did we stop (for any reason) prior to the final interim?
  stop_at <- N_analys
  
  # any na's in __any__ of the decision array rows means we stopped early:
  # we can just use sup to test this:
  if(any(is.na(decision[, 1, "sup"]))){
    # interim where the stopping rule was met
    stop_at <- min(which(is.na(decision[, 1, "sup"]))) - 1
    
    if(stop_at < N_analys){
      log_info("Stopped at analysis ", stop_at, " filling all subsequent entries")
      decision[(stop_at+1):N_analys, , "sup"] <- decision[rep(stop_at, N_analys-stop_at), , "sup"]
      decision[(stop_at+1):N_analys, , "ni"] <- decision[rep(stop_at, N_analys-stop_at), , "ni"]
      decision[(stop_at+1):N_analys, , "fut_sup"] <- decision[rep(stop_at, N_analys-stop_at), , "fut_sup"]
      decision[(stop_at+1):N_analys, , "fut_ni"] <- decision[rep(stop_at, N_analys-stop_at), , "fut_ni"]
    }
  }
  
  l_ret <- list(
    # data collected in the trial
    d_all = d_all[, .(y = sum(y), .N), keyby = .(id_analys, silo, jnt, pref_rev, d1, d2, d3, d4)],
    
    d_post_smry_1 = d_post_smry_1,
    d_post_smry_2 = d_post_smry_2,
    
    # number of units used in g-computation by domain
    n_units = n_units,
    
    decision = decision,
    
    pr_sup = pr_sup,
    pr_ni = pr_ni,
    pr_sup_fut = pr_sup_fut,
    pr_ni_fut = pr_ni_fut,
    
    stop_at = stop_at
  )
  
  if(return_posterior){
    l_ret$l_post <- l_post
  }
  # 
  
  
  return(l_ret)
}



model_consistency_check <- function(){
  
  mc_cores <- 1
  if(unname(Sys.info()[1]) == "Darwin"){
    log_info("On mac, reset cores to 5")
    mc_cores <- 5
    N_sim <- 250
  } else if (unname(Sys.info()[1]) == "Linux"){
    log_info("On linux, reset cores to 40")
    mc_cores <- 40
    N_sim <- 5000
  }
  
  dec_sup = list(
    surg = NA,
    ext_proph = NA,
    ab_choice = NA
  )
  dec_ni = list(
    ab_dur = NA
  )
  dec_sup_fut = list(
    surg = NA,
    ext_proph = NA,
    ab_choice = NA
  )
  dec_ni_fut = list(
    ab_dur = NA
  )
  
  mu <- 0
  b_silo <- c(0.0, -0.3, -0.2)
  b_jnt <- c(0, 0.4)
  b_pref <- c(0.0, -0.2)
  
  # dair, one, two-stage, we compare avg of one and two stage rev to dair
  b_d1 <- c(0, 0.693, 0.693, 0, 0, 0, 0, 0.693, 0.693)
  # always ref, 12wk, 6wk as we are assessing if 6wk ni to 12wk
  b_d2 <- c(0, -0.25, 0.5)
  # always ref, 0, 12wk as we are assessing if 12wk sup to none
  b_d3 <- c(0, -0.25, 0.5)
  # always ref, none, rif as we are assessing if rif is sup to none
  b_d4 <- c(0, -0.25, 0.5)
  N <- 2500
  
  r <- pbapply::pblapply(X=1:N_sim, cl = mc_cores, FUN = function(ix){
    
    l_new <- get_trial_data_int(
      N = N, 
      
      # reference level log odds of response
      mu = mu,
      # silo effects
      # silo 1 as the one for late acute
      # silo 2 as the one for late acute
      # silo 3 is LATE ACUTE
      b_silo = b_silo,
      b_jnt = b_jnt,
      b_pref = b_pref,
      # dair, one, two-stage
      b_d1 = b_d1,
      b_d2 = b_d2,
      b_d3 = b_d3,
      b_d4 = b_d4,
      
      dec_sup = dec_sup,
      dec_ni = dec_ni,
      dec_sup_fut = dec_sup_fut,
      dec_ni_fut = dec_ni_fut,
      
      idx_s = 1,
      t0 = rep(1, N),
      id_analys = 1
    )
    
    # create stan data format based on the relevant subsets of pt
    lsd <- get_stan_data_all_int(l_new$d)
    lsd$ld$pri_mu <- c(0, 0.47)
    lsd$ld$pri_b_silo <- c(0, 1)
    lsd$ld$pri_b_jnt <- c(0, 1)
    lsd$ld$pri_b_prf <- c(0, 1)
    lsd$ld$pri_b_trt <- c(0, 1)
    
    # dtmp <- lsd$d_mod_d2[d2 %in% c(2, 3), .(y = sum(y), n = sum(n)), keyby = d2]
    # dtmp[, p_obs := y/n]
    # dtmp[, lo_obs := log(p_obs / (1-p_obs))]
    # dtmp
    # dtmp[d2 == 3, lo_obs] - dtmp[d2 == 2, lo_obs]
    
    f_1 <- m1$sample(
      lsd$ld, iter_warmup = 1000, iter_sampling = 3000,
      parallel_chains = 1, chains = 1, refresh = 0, show_exceptions = F,
      max_treedepth = 11)
    
    # f_1 <- m1$pathfinder(lsd$ld, num_paths=20, single_path_draws=200,
    #                      history_size=50, max_lbfgs_iters=100,
    #                      refresh = 0, draws = 2000)
    
    d_post_1 <- data.table(f_1$draws(
      variables = c(
        
        "lor_d1", 
        "lor_d2",
        "lor_d3",
        "lor_d4"
        
      ),   
      format = "matrix"))
    
    c(colMeans(d_post_1))
    
  })
  
  
  d_r <- data.table(do.call(rbind, r))
  names(d_r) <- paste(rep("m1", len = ncol(d_r)), names(d_r), sep = "_")
  colMeans(d_r)
  
  d_fig <- melt(d_r, measure.vars = names(d_r))
  
  ggplot(d_fig, aes(x = value, group = variable)) +
    geom_density() +
    geom_vline(data = d_fig[, .(value = mean(value)), keyby = variable],
               aes(xintercept = value), col = 1, lwd = 0.3) +
    facet_wrap(~variable, ncol = ncol(d_r)/2)
  
  d_tbl <- copy(d_fig)
  d_tbl[, model := substr(variable, 1, 2)]
  d_tbl[, variable := gsub("m1_", "", variable, fixed = T)]
  d_tbl <- dcast(
    d_tbl[, .(mu = mean(value)), keyby = .(model, variable)],
    variable ~ model, value.var = "mu")
  d_tbl[]
  
}

model_prior_predictive <- function(){
  
  dec_sup = list(
    surg = NA,
    ext_proph = NA,
    ab_choice = NA
  )
  dec_ni = list(
    ab_dur = NA
  )
  dec_sup_fut = list(
    surg = NA,
    ext_proph = NA,
    ab_choice = NA
  )
  dec_ni_fut = list(
    ab_dur = NA
  )
  
  mu <- 0.65
  b_silo <- c(0.0, -0.3, -0.2)
  b_jnt <- c(0, 0.4)
  b_pref <- c(0.0, -0.2)
  
  # dair, one, two-stage, we compare avg of one and two stage rev to dair
  
  # null (late acute silo)
  # b_d1 <- c(0, 0.693, 0.693, 0, 0, 0, 0, 0.693, 0.693)
  # "Moderate (OR 1.75) surgical revision effect (one-stage only)"
  b_d1 <- c(0, 0.693, 0.693, -0.1, 0.4596, -0.1, 0.1, 0.4, 0.2)
  
  # always ref, 12wk, 6wk as we are assessing if 6wk ni to 12wk
  b_d2 <- c(0, -0.25, 0.5)
  # always ref, 0, 12wk as we are assessing if 12wk sup to none
  b_d3 <- c(0, -0.25, 0.5)
  # always ref, none, rif as we are assessing if rif is sup to none
  b_d4 <- c(0, -0.25, 0.5)
  N <- 2500
  
  l_new <- get_trial_data_int(
    N = N, 
    
    # reference level log odds of response
    mu = mu,
    # silo effects
    # silo 1 as the one for late acute
    # silo 2 as the one for late acute
    # silo 3 is LATE ACUTE
    b_silo = b_silo,
    b_jnt = b_jnt,
    b_pref = b_pref,
    # dair, one, two-stage
    b_d1 = b_d1,
    b_d2 = b_d2,
    b_d3 = b_d3,
    b_d4 = b_d4,
    
    dec_sup = dec_sup,
    dec_ni = dec_ni,
    dec_sup_fut = dec_sup_fut,
    dec_ni_fut = dec_ni_fut,
    
    idx_s = 1,
    t0 = rep(1, N),
    id_analys = 1
  )
  
  # create stan data format based on the relevant subsets of pt
  lsd <- get_stan_data_all_int(l_new$d)
  lsd$ld$pri_mu <- c(0.8, 0.47)
  lsd$ld$pri_b_silo <- c(0, 1)
  lsd$ld$pri_b_jnt <- c(0, 1)
  lsd$ld$pri_b_prf <- c(0, 1)
  lsd$ld$pri_b_trt <- c(0, 1)
  
  lsd$ld$prior_only <- 1
  
  # dtmp <- lsd$d_mod_d2[d2 %in% c(2, 3), .(y = sum(y), n = sum(n)), keyby = d2]
  # dtmp[, p_obs := y/n]
  # dtmp[, lo_obs := log(p_obs / (1-p_obs))]
  # dtmp
  # dtmp[d2 == 3, lo_obs] - dtmp[d2 == 2, lo_obs]
  
  f_1 <- m1$sample(
    lsd$ld, iter_warmup = 1000, iter_sampling = 3000,
    parallel_chains = 1, chains = 1, refresh = 0, show_exceptions = F,
    max_treedepth = 11)
  
  
  
  
  # extract posterior - marginal probability of outcome by group
  # dair vs rev
  d_post <- data.table(f_1$draws(
    variables = c(
      "p_hat_pred"
    ),   
    format = "matrix"))
  
  hist(d_post$p_hat_pred)
  summary(d_post$p_hat_pred)
  
  ggplot(d_fig, aes(x = p_hat_pred)) + 
    geom_density() +
    facet_grid(~d1)
  
}


run_sim_05 <- function(){
  
  log_info(paste0(match.call()[[1]]))
  
  if(unname(Sys.info()[1]) == "Darwin"){
    log_info("On mac, reset cores to 5")
    g_cfgsc$mc_cores <- 5
  }
  
  N = max(g_cfgsc$N_pt)
  
  mu <- g_cfgsc$cov$mu
  b_silo <- unlist(g_cfgsc$cov$silo)
  b_jnt <- unlist(g_cfgsc$cov$jnt)
  b_pref <- unlist(g_cfgsc$cov$pref)
  # dair, one, two-stage, we compare avg of one and two stage rev to dair
  b_d1 <- unlist(g_cfgsc$d1)
  # always ref, 12wk, 6wk as we are assessing if 6wk ni to 12wk
  b_d2 <- unlist(g_cfgsc$d2)
  # always ref, 0, 12wk as we are assessing if 12wk sup to none
  b_d3 <- unlist(g_cfgsc$d3)
  # always ref, none, rif as we are assessing if rif is sup to none
  b_d4 <- unlist(g_cfgsc$d4)
  
  list(
    mu = mu,
    b_silo = b_silo,
    b_jnt = b_jnt,
    b_pref = b_pref,
    # dair, one, two-stage, we compare avg of one and two stage rev to dair
    b_d1 = b_d1,
    # always ref, 12wk, 6wk as we are assessing if 6wk ni to 12wk
    b_d2 = b_d2,
    # always ref, 0, 12wk as we are assessing if 12wk sup to none
    b_d3 = b_d3,
    # always ref, none, rif as we are assessing if rif is sup to none
    b_d4 = b_d4
  )
  
  
  
  l_prior <- list(
    pri_mu = unlist(g_cfgsc$pri$mu),
    pri_b_silo = unlist(g_cfgsc$pri$b_silo) ,
    pri_b_jnt = unlist(g_cfgsc$pri$b_jnt) ,
    pri_b_prf = unlist(g_cfgsc$pri$b_prf) ,
    pri_b_trt = unlist(g_cfgsc$pri$b_trt) 
  )
  l_prior
  
  log_info("Starting simulation with following parameters:");
  log_info("N: ", paste0(g_cfgsc$N_pt, collapse = ", "));
  log_info("b_silo: ", paste0(b_silo, collapse = ", "));
  log_info("b_jnt: ", paste0(b_jnt, collapse = ", "));
  log_info("b_pref: ", paste0(b_pref, collapse = ", "));
  log_info("b_d1: ", paste0(b_d1, collapse = ", "));
  log_info("b_d2: ", paste0(b_d2, collapse = ", "));
  log_info("b_d3: ", paste0(b_d3, collapse = ", "));
  log_info("b_d4: ", paste0(b_d4, collapse = ", "));
  
  return_posterior = F
  
  e = NULL
  log_info("Starting simulation")
  r <- parallel::mclapply(
    X=1:g_cfgsc$nsim, mc.cores = g_cfgsc$mc_cores, FUN=function(ix) {
    # X=1:5, mc.cores = g_cfgsc$mc_cores, FUN=function(ix) {
      log_info("Simulation ", ix);
      ll <- tryCatch({
        run_trial(
          ix,
          mu = mu,
          b_silo = b_silo, b_jnt = b_jnt, b_pref = b_pref,
          # dair, one, two-stage
          b_d1 = b_d1,
          b_d2 = b_d2,
          b_d3 = b_d3,
          b_d4 = b_d4,
          l_prior = l_prior,
          return_posterior = return_posterior
          )
      },
      error=function(e) {
        log_info("ERROR in MCLAPPLY LOOP (see terminal output):")
        message(" ERROR in MCLAPPLY LOOP " , e);
        log_info("Traceback (see terminal output):")
        message(traceback())
        stop(paste0("Stopping with error ", e))
      })

      ll
    })
  
  
  
  log_info("Length of result set ", length(r))
  log_info("Sleep for 5 before processing")
  Sys.sleep(5)
  
  for(i in 1:length(r)){
    log_info("Element at index ",i, " is class ", class(r[[i]]))
    if(any(class(r[[i]]) %like% "try-error")){
      log_info("Element at index ",i, " has content ", r[[i]])  
    }
    log_info("Element at index ",i, " has names ", 
             paste0(names(r[[i]]), collapse = ", "))
  }
  
  
  d_pr_sup <- data.table()
  for(i in 1:length(r)){
    
    log_info("Appending pr_sup for result ", i)
    
    if(is.recursive(r[[i]])){
      d_pr_sup <- rbind(
        d_pr_sup,
        cbind(
          sim = i, analys = as.integer(rownames(r[[i]]$pr_sup)), r[[i]]$pr_sup
        ) 
      )  
    } else {
      log_info("Value for r at this index is not recursive ", i)
      message("r[[i]] ", r[[i]])
      message(traceback())
      stop(paste0("Stopping due to non-recursive element "))
    }
    
  }

  d_pr_ni <- data.table(
    do.call(rbind, lapply(1:length(r), function(i){ 
      cbind(
        sim = i, analys = as.integer(rownames(r[[i]]$pr_ni)), r[[i]]$pr_ni) 
    } )))
  
  d_pr_sup_fut <- data.table(
    do.call(rbind, lapply(1:length(r), function(i){ 
      cbind(
        sim = i, analys = as.integer(rownames(r[[i]]$pr_sup_fut)), r[[i]]$pr_sup_fut) 
    } )))
  
  d_pr_ni_fut <- data.table(
    do.call(rbind, lapply(1:length(r), function(i){ 
      cbind(
        sim = i, analys = as.integer(rownames(r[[i]]$pr_ni_fut)), r[[i]]$pr_ni_fut) 
    } )))

  # decisions based on probability thresholds
  d_decision <- rbind(
    
    data.table(do.call(rbind, lapply(1:length(r), function(i){ 
      m <- r[[i]]$decision[, , "sup"]
      cbind(sim = i, analys = as.integer(rownames(m)), quant = "sup", m) 
    } ))), 
    
    data.table(do.call(rbind, lapply(1:length(r), function(i){ 
      m <- r[[i]]$decision[, , "ni"]
      cbind(sim = i, analys = as.integer(rownames(m)), quant = "ni", m) 
    } ))),
    
    data.table(do.call(rbind, lapply(1:length(r), function(i){ 
      m <- r[[i]]$decision[, , "fut_sup"]
      cbind(sim = i, analys = as.integer(rownames(m)), quant = "fut_sup", m) 
    } ))),
    
    data.table(do.call(rbind, lapply(1:length(r), function(i){ 
      m <- r[[i]]$decision[, , "fut_ni"]
      cbind(sim = i, analys = as.integer(rownames(m)), quant = "fut_ni", m) 
    } )))
    
  )
  d_decision[, `:=`(sim = as.integer(sim), analys = as.integer(analys),
                d1 = as.logical(d1), 
                d2 = as.logical(d2), 
                d3 = as.logical(d3),  
                d4 = as.logical(d4) 
                )]
  
  d_post_smry_1 <- rbindlist(lapply(1:length(r), function(i){ 
    r[[i]]$d_post_smry_1
  } ), idcol = "sim")
  
  d_post_smry_2 <- rbindlist(lapply(1:length(r), function(i){ 
    r[[i]]$d_post_smry_2
  } ), idcol = "sim")
  
  # data from each simulated trial
  d_all <- rbindlist(lapply(1:length(r), function(i){ 
    r[[i]]$d_all
  } ), idcol = "sim")
  
  # number of units informing estimates using g-comp by analys
  d_n_units <- data.table(
    do.call(rbind, lapply(1:length(r), function(i){ 
      cbind(
        sim = i, analys = as.integer(rownames(r[[i]]$n_units)), r[[i]]$n_units) 
    } )))
  
  # number of units assigned to each trt level
  d_n_assign <- data.table(
    do.call(rbind, lapply(1:length(r), function(i){ 
      
      rbind(
        cbind(sim = i, silo = 1, domain = 1, analys = as.integer(rownames(r[[i]]$n_assign[,,"d1","s1"])), r[[i]]$n_assign[,,"d1","s1"]),
        cbind(sim = i, silo = 1, domain = 2, analys = as.integer(rownames(r[[i]]$n_assign[,,"d2","s1"])), r[[i]]$n_assign[,,"d2","s1"]),
        cbind(sim = i, silo = 1, domain = 3, analys = as.integer(rownames(r[[i]]$n_assign[,,"d2","s1"])), r[[i]]$n_assign[,,"d2","s1"]),
        cbind(sim = i, silo = 1, domain = 4, analys = as.integer(rownames(r[[i]]$n_assign[,,"d2","s1"])), r[[i]]$n_assign[,,"d2","s1"]),
        
        cbind(sim = i, silo = 2, domain = 1, analys = as.integer(rownames(r[[i]]$n_assign[,,"d1","s2"])), r[[i]]$n_assign[,,"d1","s2"]),
        cbind(sim = i, silo = 2, domain = 2, analys = as.integer(rownames(r[[i]]$n_assign[,,"d2","s2"])), r[[i]]$n_assign[,,"d2","s2"]),
        cbind(sim = i, silo = 2, domain = 3, analys = as.integer(rownames(r[[i]]$n_assign[,,"d2","s2"])), r[[i]]$n_assign[,,"d2","s2"]),
        cbind(sim = i, silo = 2, domain = 4, analys = as.integer(rownames(r[[i]]$n_assign[,,"d2","s2"])), r[[i]]$n_assign[,,"d2","s2"]),
        
        cbind(sim = i, silo = 3, domain = 1, analys = as.integer(rownames(r[[i]]$n_assign[,,"d1","s3"])), r[[i]]$n_assign[,,"d1","s3"]),
        cbind(sim = i, silo = 3, domain = 2, analys = as.integer(rownames(r[[i]]$n_assign[,,"d2","s3"])), r[[i]]$n_assign[,,"d2","s3"]),
        cbind(sim = i, silo = 3, domain = 3, analys = as.integer(rownames(r[[i]]$n_assign[,,"d2","s3"])), r[[i]]$n_assign[,,"d2","s3"]),
        cbind(sim = i, silo = 3, domain = 4, analys = as.integer(rownames(r[[i]]$n_assign[,,"d2","s3"])), r[[i]]$n_assign[,,"d2","s3"])
      )
      
    } )))


  l <- list(
    cfg = g_cfgsc,
    
    # reference risk
    # d_risk_smry = d_risk_smry,
    
    d_pr_sup = d_pr_sup, 
    d_pr_ni = d_pr_ni,
    
    d_pr_sup_fut = d_pr_sup_fut, 
    d_pr_ni_fut = d_pr_ni_fut,
    
    d_decision = d_decision,
    d_post_smry_1 = d_post_smry_1,
    d_all = d_all,
    
    d_n_units = d_n_units,
    d_n_assign = d_n_assign,
    
    d_post_smry_2 = d_post_smry_2
    # d_grp = d_grp
    )
  
  toks <- unlist(tstrsplit(args[2], "[-.]"))
  fname <- paste0("data/sim05/sim05-", toks[4], "-", toks[5], "-", format(Sys.time(), "%Y%m%d-%H%M%S"), ".qs")
  
  log_info("Saving results file ", fname)
  
  qs::qsave(l, file = fname)
}

run_none_sim_05 <- function(){
  log_info("run_none_sim_05: Nothing doing here bud.")
}

main_sim_05 <- function(){
  funcname <- paste0(args[1], "()")
  log_info("Main, invoking ", funcname)
  eval(parse(text=funcname))
}

main_sim_05()


