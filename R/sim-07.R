# Experiment with independent set of models with reduced linear predictor.

source("./R/init.R")
source("./R/data-sim07.R")

# Command line arguments list scenario (true dose response),
# the number of simulations to run, the number of cores
# to use and the simulator to use.
args = commandArgs(trailingOnly=TRUE)

# Load cfg based on cmd line args.
if (length(args)<1) {
  log_info("Setting default run method (does nothing)")
  args[1] = "run_none_sim_07"
  args[2] = "./sim07/cfg-sim07-sc01-v02.yml"
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
m1 <- cmdstanr::cmdstan_model("stan/model-sim-07-a.stan")

output_dir_mcmc <- paste0(getwd(), "/tmp")



# Main trial loop.
run_trial <- function(
    ix,
    l_spec,
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
  
  loc_t0 <- get_sim07_enrol_time_int(sum(l_spec$N), lambda, rho)
  
  # loop controls
  stop_enrol <- FALSE
  l_spec$ia <- 1 # interim number
  N_analys <- length(l_spec$N)

  # posterior summaries
  d_post_smry_1 <- CJ(
    ia = 1:N_analys,
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
    
    log_info("Trial ", ix, " analysis ", l_spec$ia)
    
    # next chunk of data on pts.
    if(l_spec$ia == 1){
      # starting pt index in data
      l_spec$is <- 1
      l_spec$ie <- l_spec$is + l_spec$N[l_spec$ia] - 1
    } else {
      l_spec$is <- l_spec$N[l_spec$ia-1] + 1
      l_spec$ie <- l_spec$is + l_spec$N[l_spec$ia] - 1
    }
    
    # id and time
    l_spec$t0 = loc_t0[l_spec$is:l_spec$ie]
    
    # Our analyses only occur on those that have reached 12 months post 
    # randomisation. As such, we are assuming that the analysis takes place
    # 12 months following the last person to be enrolled in the current 
    # analysis set.
    
    d <- get_sim07_trial_data(
      l_spec,
      dec_sup, dec_ni, dec_sup_fut, dec_ni_fut)
    
    log_info("Trial ", ix, " new data generated ", l_spec$ia)

    # combine the existing and new data
    d_all <- rbind(d_all, d)
    
    # create stan data format based on the relevant subsets of pt
    lsd <- get_sim07_stan_data(d_all)
    
    lsd$ld$pri_mu <- l_spec$prior$mu
    lsd$ld$pri_bs <- l_spec$prior$bs
    lsd$ld$pri_bp <- l_spec$prior$bp
    lsd$ld$pri_b1 <- l_spec$prior$bd1
    lsd$ld$pri_b2 <- l_spec$prior$bd2
    lsd$ld$pri_b3 <- l_spec$prior$bd3
    lsd$ld$pri_b4 <- l_spec$prior$bd4
    
    foutname <- paste0(
      format(Sys.time(), format = "%Y%m%d%H%M%S"),
      "-sim-", ix, "-intrm-", l_spec$ia)
    
    
    
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
    
    log_info("Trial ", ix, " fitted models ", l_spec$ia)

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
      l_post[[l_spec$ia]] <- copy(d_post)
    }
    
    d_post_long <- melt(d_post, measure.vars = names(d_post))
    d_post_long[variable %like% "d1", domain := 1]
    d_post_long[variable %like% "d2", domain := 2]
    d_post_long[variable %like% "d3", domain := 3]
    d_post_long[variable %like% "d4", domain := 4]
    
    # d_post_long[variable %like% "lor", .(mu = mean(value), sd = sd(value)), keyby = variable]
  
    d_post_smry_1[ia == l_spec$ia, 
                  mu := d_post_long[, mean(value), 
                                    keyby = .(domain, variable)]$V1] 
    d_post_smry_1[ia == l_spec$ia, 
                  med := d_post_long[, median(value), 
                                    keyby = .(domain, variable)]$V1]
    d_post_smry_1[ia == l_spec$ia, 
                  se := d_post_long[, sd(value), 
                                     keyby = .(domain, variable)]$V1]
    d_post_smry_1[ia == l_spec$ia, 
                  q_025 := d_post_long[, quantile(value, prob = 0.025), 
                                    keyby = .(domain, variable)]$V1]
    d_post_smry_1[ia == l_spec$ia, 
                  q_975 := d_post_long[, quantile(value, prob = 0.975), 
                                    keyby = .(domain, variable)]$V1]

    # These should produce the same answers as the generated quantities
    # outputs.
    d_post_chk <- data.table(f_1$draws(
      variables = c(
        
        "bd1", "bd2", "bd3", "bd4"
        
      ),  
      format = "matrix"))
    d_post_chk <- melt(d_post_chk, measure.vars = names(d_post_chk))
    d_post_chk[variable %like% "d1", domain := 1]
    d_post_chk[variable %like% "d2", domain := 2]
    d_post_chk[variable %like% "d3", domain := 3]
    d_post_chk[variable %like% "d4", domain := 4]
    
    d_post_smry_2 <- rbind(
      d_post_smry_2,
      d_post_chk[, .(ia = l_spec$ia,
                     mu = mean(value),
                     med = median(value),
                     se = sd(value),
                     q_025 = quantile(value, prob = 0.025),
                     q_975 = quantile(value, prob = 0.975)),
                 keyby = .(variable, domain)]
    )
    
    log_info("Trial ", ix, " extracted posterior ", l_spec$ia)

    # this is how many units are used in the g-comp up to this analysis
    n_units[l_spec$ia, ] <- c(
      sum(lsd$ld$n_d1),
      sum(lsd$ld$n_d2),
      sum(lsd$ld$n_d3),
      sum(lsd$ld$n_d4)
    )
    
    # superiority is implied by a high probability that the risk diff 
    # greater than zero
    pr_sup[l_spec$ia, ]  <-   c(
      d_post[, mean(p_d1_23 - p_d1_1 > l_spec$delta$sup)],
      d_post[, mean(p_d2_3 - p_d2_2 > l_spec$delta$sup)],
      d_post[, mean(p_d3_3 - p_d3_2 > l_spec$delta$sup)],
      d_post[, mean(p_d4_3 - p_d4_2 > l_spec$delta$sup)]
    )
    # futility for the superiority decision is implied by a low probability 
    # that the risk diff is greater than some small value (2% difference)
    pr_sup_fut[l_spec$ia, ]  <-   c(
      d_post[, mean(p_d1_23 - p_d1_1 > l_spec$delta$sup_fut)],
      d_post[, mean(p_d2_3 - p_d2_2 > l_spec$delta$sup_fut)],
      d_post[, mean(p_d3_3 - p_d3_2 > l_spec$delta$sup_fut)],
      d_post[, mean(p_d4_3 - p_d4_2 > l_spec$delta$sup_fut)]
    )
    
    # ni is implied by a high probability that the risk diff reduces the 
    # efficacy by no more than some small amount (here 2%)
    pr_ni[l_spec$ia, ]  <-   c(
      d_post[, mean(p_d1_23 - p_d1_1 > l_spec$delta$ni)],
      d_post[, mean(p_d2_3 - p_d2_2 > l_spec$delta$ni)],
      d_post[, mean(p_d3_3 - p_d3_2 > l_spec$delta$ni)],
      d_post[, mean(p_d4_3 - p_d4_2 > l_spec$delta$ni)]
    )
    # futility for the ni decision is implied by a low probability that 
    # the risk diff suggests any efficacy   
    pr_ni_fut[l_spec$ia, ]  <-   c(
      d_post[, mean(p_d1_23 - p_d1_1 > l_spec$delta$ni_fut)],
      d_post[, mean(p_d2_3 - p_d2_2 > l_spec$delta$ni_fut)],
      d_post[, mean(p_d3_3 - p_d3_2 > l_spec$delta$ni_fut)],
      d_post[, mean(p_d4_3 - p_d4_2 > l_spec$delta$ni_fut)]
    )

    log_info("Trial ", ix, " calculated decision quantities ", l_spec$ia)
    
    decision[l_spec$ia, , "sup"] <- pr_sup[l_spec$ia, ] > l_spec$thresh$sup
    decision[l_spec$ia, , "ni"] <- pr_ni[l_spec$ia, ] > l_spec$thresh$ni
    
    # futility rules
    # taken to imply negligible chance of being superior or ni
    
    # negligible chance of being superior
    decision[l_spec$ia, , "fut_sup"] <- pr_sup_fut[l_spec$ia, ] < l_spec$thresh$fut_sup
    # taken to imply negligible chance of being ni 
    decision[l_spec$ia, , "fut_ni"] <- pr_ni_fut[l_spec$ia, ] < l_spec$thresh$fut_ni
    
    log_info("Trial ", ix, " compared to thresholds ", l_spec$ia)
    
    # Earlier decisions are retained - once superiority has been decided, 
    # we retain this conclusion irrespective of subsequent post probabilities.
    # The following simply overwrites any decision reversal.
    # This means that there could be inconsistency with a silo pr_sup and the 
    # decision reported.
    decision[1:l_spec$ia, , "sup"] <- apply(decision[1:l_spec$ia, , "sup", drop = F], 2, function(z){ cumsum(z) > 0 })
    decision[1:l_spec$ia, , "ni"] <- apply(decision[1:l_spec$ia, , "ni", drop = F], 2, function(z){ cumsum(z) > 0 })
    decision[1:l_spec$ia, , "fut_sup"] <- apply(decision[1:l_spec$ia, , "fut_sup", drop = F], 2, function(z){ cumsum(z) > 0 })
    decision[1:l_spec$ia, , "fut_ni"] <- apply(decision[1:l_spec$ia, , "fut_ni", drop = F], 2, function(z){ cumsum(z) > 0 })
  
    # superiority decisions apply to domains 1, 3 and 4
    if(any(decision[l_spec$ia, , "sup"])){
      # since there are only two treatments per cell, if a superiority decision 
      # is made then we have answered all the questions and we can stop 
      # enrolling into that cell. if there were more than two treatments then
      # we would need to take a different approach.
      
      # if we have found it to be superior to dair then all subsequent get 
      # revision. the value of the dec_sup$surg doesn't matter at this time,
      # it just needs to be something other than NA. the logic for selecting
      # the data is in the get_data function.
      # this decision would only impact late acute
      if(decision[l_spec$ia, "d1", "sup"]){
        dec_sup$surg <- 3
      }
      # we are not assessing superiority for domain 2 (antibiotic duration)
      
      # extproph domain relates to all silo so decision on b_r2d impacts all cohorts
      if(decision[l_spec$ia, "d3", "sup"]){
        dec_sup$ext_proph <- 3
      }
      # choice domain relates to all silo so decision on b_f impacts all cohorts
      if(decision[l_spec$ia, "d4", "sup"]){
        dec_sup$ab_choice <- 3
      }
    }
    # stop enrolling if futile wrt superiority decision
    if(any(decision[l_spec$ia, , "fut_sup"])){
      if(decision[l_spec$ia, "d1", "fut_sup"]){
        dec_sup_fut$surg <- 3
      }
      if(decision[l_spec$ia, "d3", "fut_sup"]){
        dec_sup_fut$ext_proph <- 3
      }
      if(decision[l_spec$ia, "d4", "fut_sup"]){
        dec_sup_fut$ab_choice <- 3
      }
    }
    
    
    # ni decisions apply to domains 2
    # if short is ni to long then stop enrolment
    # for this randomised comparison
    if(any(decision[l_spec$ia, , "ni"])){
      if(decision[l_spec$ia, "d2", "ni"]){
        dec_ni$ab_dur <- 3
      }
    }
    
    if(any(decision[l_spec$ia, , "fut_ni"])){
      if(decision[l_spec$ia, "d2", "fut_ni"]){
        dec_ni_fut$ab_dur <- 3
      }
    }
    
    # have we answered all questions of interest?
    if(
      # if rev (in late acute silo) is superior to dair (or superiority decision
      # decision is futile to pursue)
      (decision[l_spec$ia, "d1", "sup"] | decision[l_spec$ia, "d1", "fut_sup"]) &
      # if 6wk backbone duration (in one-stage units) is non-inferior to 12wk
      # (or decision for non-inferiority is futile to pursue)
      (decision[l_spec$ia, "d2", "ni"] | decision[l_spec$ia, "d2", "fut_ni"] ) &
      # if 12wk ext-proph duration (in two-stage units) is sup to none
      # (or decision for sup is futile to pursue)
      (decision[l_spec$ia, "d3", "sup"] | decision[l_spec$ia, "d3", "fut_sup"]) &
      # if rif (in all relevant units) is sup to none 
      # (or decision for sup is futile to pursue)
      (decision[l_spec$ia, "d4", "sup"] | decision[l_spec$ia, "d4", "fut_sup"])
    ){
      log_info("Stop trial all questions addressed ", ix)
      stop_enrol <- T  
    } 
    
    log_info("Trial ", ix, " updated allocation control ", l_spec$ia)
    
    # next interim
    l_spec$ia <- l_spec$ia + 1
    
    if(l_spec$ia > N_analys){
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
    d_all = d_all[, .(y = sum(y), .N), keyby = .(ia, s, pref, d1, d2, d3, d4)],
    
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


# run_sim_07_tmp <- function(){
#   
#   
#   for(i in 1:lsd$ld$N){
#     # d1_ix[i] = 
#     cat(lsd$ld$d1[i] + (lsd$ld$K_d1 * (lsd$ld$s[i] - 1)), "\n")
#   } 
#   # // for g-comp to pick up the correct surgical domain parameter
#   for(i in 1:lsd$ld$N_d2){
#     # d2_d1_ix[i] = 
#     cat(lsd$ld$d2_d1[i] + (lsd$ld$K_d1 * (lsd$ld$d2_s[i] - 1)), "\n")
#   }
#   for(i in 1:lsd$ld$N_d3){
#     # d3_d1_ix[i] = 
#     cat(lsd$ld$d3_d1[i] + (lsd$ld$K_d1 * (lsd$ld$d3_s[i] - 1)), "\n")
#   }
#   for(i in 1:lsd$ld$N_d4){
#     # d4_d1_ix[i]
#     cat(lsd$ld$d4_d1[i] + (lsd$ld$K_d1 * (lsd$ld$d4_s[i] - 1)), "\n")
#   }
#   
#   
#   
#   
# }

run_sim_07_model_check <- function(){
  
  
  log_info(paste0(match.call()[[1]]))
  
  if(unname(Sys.info()[1]) == "Darwin"){
    log_info("On mac, reset cores to 5")
    mc_cores <- 5
  }
  
  l_spec <- list(
    N = 3e3,
    # silo allocation
    p_s_alloc = c(0.3, 0.5, 0.2),
    l_e = list(),
    l_l = list(),
    l_c = list(),
    # model parameters
    # intercept is early silo
    mu = NA, # 0.9,
    bs = c(0, -0.1, -0.2),
    # different baseline risk for rev
    bp = -0.4,
    # index 4 is never referenced but needs to be there so that the 
    # calcs in the single model approach don't end up with NA
    bd2 = c(0, 0, 0, 999),
    bd3 = c(0, 0, 0, 999),
    bd4 = c(0, 0, 0)
  )
  # N by analysis
  l_spec$N <- g_cfgsc$N_pt
  l_spec$l_e$p_d1_alloc <- g_cfgsc$e_p_d1_alloc
  l_spec$l_e$p_d2_entry <- g_cfgsc$e_p_d2_entry
  l_spec$l_e$p_d2_alloc <- g_cfgsc$e_p_d2_alloc
  l_spec$l_e$p_d3_entry <- g_cfgsc$e_p_d3_entry
  l_spec$l_e$p_d3_alloc <- g_cfgsc$e_p_d3_alloc
  l_spec$l_e$p_d4_entry <- g_cfgsc$e_p_d4_entry
  l_spec$l_e$p_d4_alloc <- g_cfgsc$e_p_d4_alloc
  # preference for two-stage
  l_spec$l_e$p_pref <- g_cfgsc$e_p_pref
  
  l_spec$l_l$p_d1_alloc <- g_cfgsc$l_p_d1_alloc
  l_spec$l_l$p_d2_entry <- g_cfgsc$l_p_d2_entry
  l_spec$l_l$p_d2_alloc <- g_cfgsc$l_p_d2_alloc
  l_spec$l_l$p_d3_entry <- g_cfgsc$l_p_d3_entry
  l_spec$l_l$p_d3_alloc <- g_cfgsc$l_p_d3_alloc
  l_spec$l_l$p_d4_entry <- g_cfgsc$l_p_d4_entry
  l_spec$l_l$p_d4_alloc <- g_cfgsc$l_p_d4_alloc
  # preference for two-stage
  l_spec$l_l$p_pref <- g_cfgsc$l_p_pref
  
  l_spec$l_c$p_d1_alloc <- g_cfgsc$c_p_d1_alloc
  l_spec$l_c$p_d2_entry <- g_cfgsc$c_p_d2_entry
  l_spec$l_c$p_d2_alloc <- g_cfgsc$c_p_d2_alloc
  l_spec$l_c$p_d3_entry <- g_cfgsc$c_p_d3_entry
  l_spec$l_c$p_d3_alloc <- g_cfgsc$c_p_d3_alloc
  l_spec$l_c$p_d4_entry <- g_cfgsc$c_p_d4_entry
  l_spec$l_c$p_d4_alloc <- g_cfgsc$c_p_d4_alloc
  # preference for two-stage
  l_spec$l_c$p_pref <- g_cfgsc$c_p_pref
  
  # model params
  l_spec$mu <- g_cfgsc$bmu
  l_spec$bs <- unlist(g_cfgsc$bs)
  l_spec$bp <- unlist(g_cfgsc$bp)
  # dair, one, two-stage, we compare avg of one and two stage rev to dair
  l_spec$l_e$bd1 <- unlist(g_cfgsc$bed1)
  l_spec$l_l$bd1 <- unlist(g_cfgsc$bld1)
  l_spec$l_c$bd1 <- unlist(g_cfgsc$bcd1)
  # always ref, 12wk, 6wk as we are assessing if 6wk ni to 12wk
  l_spec$bd2 <- unlist(g_cfgsc$bd2)
  # always ref, 0, 12wk as we are assessing if 12wk sup to none
  l_spec$bd3 <- unlist(g_cfgsc$bd3)
  # always ref, none, rif as we are assessing if rif is sup to none
  l_spec$bd4 <- unlist(g_cfgsc$bd4)
  
  l_spec$prior <- list()
  # location, scale
  l_spec$prior$mu <- unlist(g_cfgsc$pri_bmu)
  l_spec$prior$bs <- unlist(g_cfgsc$pri_bs)
  l_spec$prior$bp <- unlist(g_cfgsc$pri_bp)
  l_spec$prior$bd1 <- unlist(g_cfgsc$pri_bd1)
  l_spec$prior$bd2 <- unlist(g_cfgsc$pri_bd2)
  l_spec$prior$bd3 <- unlist(g_cfgsc$pri_bd3)
  l_spec$prior$bd4 <- unlist(g_cfgsc$pri_bd4)
  
  l_spec$delta <- list()
  l_spec$delta$sup <- g_cfgsc$dec_delta_sup
  l_spec$delta$sup_fut <- g_cfgsc$dec_delta_sup_fut
  l_spec$delta$ni <- g_cfgsc$dec_delta_ni
  l_spec$delta$ni_fut <- g_cfgsc$dec_delta_ni_fut
  
  # domain specific
  l_spec$thresh <- list()
  l_spec$thresh$sup <- unlist(g_cfgsc$dec_thresh_sup)
  l_spec$thresh$ni <- unlist(g_cfgsc$dec_thresh_fut_sup)
  l_spec$thresh$sup_fut <- unlist(g_cfgsc$dec_thresh_ni)
  l_spec$thresh$ni_fut <- unlist(g_cfgsc$dec_thresh_fut_ni)
  
  
  l_spec$N <- 3e3
  
  str(l_spec)
  
  
  log_info("Starting simulation with following parameters:");
  log_info("N: ", paste0(l_spec$N, collapse = ", "));
  log_info("b_silo: ", paste0(l_spec$bs, collapse = ", "));
  log_info("b_pref: ", paste0(l_spec$bp, collapse = ", "));
  log_info("b_d1: ", paste0(c(l_spec$l_e$bd1, l_spec$l_l$bd1, l_spec$l_c$bd1), collapse = ", "));
  log_info("b_d2: ", paste0(l_spec$bd2, collapse = ", "));
  log_info("b_d3: ", paste0(l_spec$bd3, collapse = ", "));
  log_info("b_d4: ", paste0(l_spec$bd4, collapse = ", "));
  
  
  N_sim <- 500
  r <- pbapply::pblapply(X=1:N_sim, cl = mc_cores, FUN = function(ix){
    
    d <- get_sim07_trial_data(l_spec)
    # d[]
    
    # combine the existing and new data
    d_all <- copy(d)
    
    # create stan data format based on the relevant subsets of pt
    lsd <- get_sim07_stan_data(d_all)
    
    lsd$ld$pri_mu <- l_spec$prior$mu
    lsd$ld$pri_bs <- l_spec$prior$bs
    lsd$ld$pri_bp <- l_spec$prior$bp
    lsd$ld$pri_b1 <- l_spec$prior$bd1
    lsd$ld$pri_b2 <- l_spec$prior$bd2
    lsd$ld$pri_b3 <- l_spec$prior$bd3
    lsd$ld$pri_b4 <- l_spec$prior$bd4
    
    foutname <- paste0(
      format(Sys.time(), format = "%Y%m%d%H%M%S"),
      "-sim-", 1, "-intrm-", ix)
    
    # fit model - does it matter that I continue to fit the model after the
    # decision is made...?
    f_1 <- m1$sample(
      lsd$ld, iter_warmup = 1000, iter_sampling = 2000,
      parallel_chains = 1, chains = 1, refresh = 0, show_exceptions = F,
      max_treedepth = 11,
      output_dir = output_dir_mcmc,
      output_basename = foutname
    )
    
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
    colMeans(d_post)
    
    
  })
  d_fig <- data.table(do.call(rbind, r))
  
  d_fig <- melt(d_fig, measure.vars = names(d_fig))
  
  ggplot(d_fig, aes(x = value, group = variable)) +
    geom_density() +
    geom_vline(data = d_fig[, .(mu = mean(value)), keyby = variable],
               aes(xintercept = mu)) +
    facet_wrap2(~variable, scales = "free_x")
  #

  
}



run_sim_07 <- function(){
  
  log_info(paste0(match.call()[[1]]))
  
  if(unname(Sys.info()[1]) == "Darwin"){
    log_info("On mac, reset cores to 5")
    g_cfgsc$mc_cores <- 5
  }
  
  l_spec <- list(
    N = 3e3,
    # silo allocation
    p_s_alloc = c(0.3, 0.5, 0.2),
    l_e = list(),
    l_l = list(),
    l_c = list(),
    # model parameters
    # intercept is early silo
    mu = NA, # 0.9,
    bs = c(0, -0.1, -0.2),
    # different baseline risk for rev
    bp = -0.4,
    # index 4 is never referenced but needs to be there so that the 
    # calcs in the single model approach don't end up with NA
    bd2 = c(0, 0, 0, 999),
    bd3 = c(0, 0, 0, 999),
    bd4 = c(0, 0, 0)
  )
  # N by analysis
  l_spec$N <- g_cfgsc$N_pt
  l_spec$l_e$p_d1_alloc <- g_cfgsc$e_p_d1_alloc
  l_spec$l_e$p_d2_entry <- g_cfgsc$e_p_d2_entry
  l_spec$l_e$p_d2_alloc <- g_cfgsc$e_p_d2_alloc
  l_spec$l_e$p_d3_entry <- g_cfgsc$e_p_d3_entry
  l_spec$l_e$p_d3_alloc <- g_cfgsc$e_p_d3_alloc
  l_spec$l_e$p_d4_entry <- g_cfgsc$e_p_d4_entry
  l_spec$l_e$p_d4_alloc <- g_cfgsc$e_p_d4_alloc
  # preference for two-stage
  l_spec$l_e$p_pref <- g_cfgsc$e_p_pref
  
  l_spec$l_l$p_d1_alloc <- g_cfgsc$l_p_d1_alloc
  l_spec$l_l$p_d2_entry <- g_cfgsc$l_p_d2_entry
  l_spec$l_l$p_d2_alloc <- g_cfgsc$l_p_d2_alloc
  l_spec$l_l$p_d3_entry <- g_cfgsc$l_p_d3_entry
  l_spec$l_l$p_d3_alloc <- g_cfgsc$l_p_d3_alloc
  l_spec$l_l$p_d4_entry <- g_cfgsc$l_p_d4_entry
  l_spec$l_l$p_d4_alloc <- g_cfgsc$l_p_d4_alloc
  # preference for two-stage
  l_spec$l_l$p_pref <- g_cfgsc$l_p_pref
  
  l_spec$l_c$p_d1_alloc <- g_cfgsc$c_p_d1_alloc
  l_spec$l_c$p_d2_entry <- g_cfgsc$c_p_d2_entry
  l_spec$l_c$p_d2_alloc <- g_cfgsc$c_p_d2_alloc
  l_spec$l_c$p_d3_entry <- g_cfgsc$c_p_d3_entry
  l_spec$l_c$p_d3_alloc <- g_cfgsc$c_p_d3_alloc
  l_spec$l_c$p_d4_entry <- g_cfgsc$c_p_d4_entry
  l_spec$l_c$p_d4_alloc <- g_cfgsc$c_p_d4_alloc
  # preference for two-stage
  l_spec$l_c$p_pref <- g_cfgsc$c_p_pref
  
  # model params
  l_spec$mu <- g_cfgsc$bmu
  l_spec$bs <- unlist(g_cfgsc$bs)
  l_spec$bp <- unlist(g_cfgsc$bp)
  # dair, one, two-stage, we compare avg of one and two stage rev to dair
  l_spec$l_e$bd1 <- unlist(g_cfgsc$bed1)
  l_spec$l_l$bd1 <- unlist(g_cfgsc$bld1)
  l_spec$l_c$bd1 <- unlist(g_cfgsc$bcd1)
  # always ref, 12wk, 6wk as we are assessing if 6wk ni to 12wk
  l_spec$bd2 <- unlist(g_cfgsc$bd2)
  # always ref, 0, 12wk as we are assessing if 12wk sup to none
  l_spec$bd3 <- unlist(g_cfgsc$bd3)
  # always ref, none, rif as we are assessing if rif is sup to none
  l_spec$bd4 <- unlist(g_cfgsc$bd4)
  
  l_spec$prior <- list()
  # location, scale
  l_spec$prior$mu <- unlist(g_cfgsc$pri_bmu)
  l_spec$prior$bs <- unlist(g_cfgsc$pri_bs)
  l_spec$prior$bp <- unlist(g_cfgsc$pri_bp)
  l_spec$prior$bd1 <- unlist(g_cfgsc$pri_bd1)
  l_spec$prior$bd2 <- unlist(g_cfgsc$pri_bd2)
  l_spec$prior$bd3 <- unlist(g_cfgsc$pri_bd3)
  l_spec$prior$bd4 <- unlist(g_cfgsc$pri_bd4)
  
  l_spec$delta <- list()
  l_spec$delta$sup <- g_cfgsc$dec_delta_sup
  l_spec$delta$sup_fut <- g_cfgsc$dec_delta_sup_fut
  l_spec$delta$ni <- g_cfgsc$dec_delta_ni
  l_spec$delta$ni_fut <- g_cfgsc$dec_delta_ni_fut
  
  # domain specific
  l_spec$thresh <- list()
  l_spec$thresh$sup <- unlist(g_cfgsc$dec_thresh_sup)
  l_spec$thresh$fut_sup <- unlist(g_cfgsc$dec_thresh_fut_sup)
  l_spec$thresh$ni <- unlist(g_cfgsc$dec_thresh_ni)
  l_spec$thresh$fut_ni <- unlist(g_cfgsc$dec_thresh_fut_ni)
  
  log_info("Starting simulation with following parameters:");
  log_info("N: ", paste0(l_spec$N, collapse = ", "));
  log_info("b_silo: ", paste0(l_spec$bs, collapse = ", "));
  log_info("b_pref: ", paste0(l_spec$bp, collapse = ", "));
  log_info("b_d1: ", paste0(c(l_spec$l_e$bd1, l_spec$l_l$bd1, l_spec$l_c$bd1), collapse = ", "));
  log_info("b_d2: ", paste0(l_spec$bd2, collapse = ", "));
  log_info("b_d3: ", paste0(l_spec$bd3, collapse = ", "));
  log_info("b_d4: ", paste0(l_spec$bd4, collapse = ", "));
  
  # str(l_spec)
  
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
          l_spec,
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
          sim = i, ia = as.integer(rownames(r[[i]]$pr_sup)), r[[i]]$pr_sup
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
        sim = i, ia = as.integer(rownames(r[[i]]$pr_ni)), r[[i]]$pr_ni) 
    } )))
  
  d_pr_sup_fut <- data.table(
    do.call(rbind, lapply(1:length(r), function(i){ 
      cbind(
        sim = i, ia = as.integer(rownames(r[[i]]$pr_sup_fut)), r[[i]]$pr_sup_fut) 
    } )))
  
  d_pr_ni_fut <- data.table(
    do.call(rbind, lapply(1:length(r), function(i){ 
      cbind(
        sim = i, ia = as.integer(rownames(r[[i]]$pr_ni_fut)), r[[i]]$pr_ni_fut) 
    } )))

  # decisions based on probability thresholds
  d_decision <- rbind(
    
    data.table(do.call(rbind, lapply(1:length(r), function(i){ 
      m <- r[[i]]$decision[, , "sup"]
      cbind(sim = i, ia = as.integer(rownames(m)), quant = "sup", m) 
    } ))), 
    
    data.table(do.call(rbind, lapply(1:length(r), function(i){ 
      m <- r[[i]]$decision[, , "ni"]
      cbind(sim = i, ia = as.integer(rownames(m)), quant = "ni", m) 
    } ))),
    
    data.table(do.call(rbind, lapply(1:length(r), function(i){ 
      m <- r[[i]]$decision[, , "fut_sup"]
      cbind(sim = i, ia = as.integer(rownames(m)), quant = "fut_sup", m) 
    } ))),
    
    data.table(do.call(rbind, lapply(1:length(r), function(i){ 
      m <- r[[i]]$decision[, , "fut_ni"]
      cbind(sim = i, ia = as.integer(rownames(m)), quant = "fut_ni", m) 
    } )))
    
  )
  d_decision[, `:=`(sim = as.integer(sim), ia = as.integer(ia),
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
        sim = i, ia = as.integer(rownames(r[[i]]$n_units)), r[[i]]$n_units) 
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
    d_post_smry_2 = d_post_smry_2,
    
    d_all = d_all
    
    # d_n_units = d_n_units,
    # d_n_assign = d_n_assign,
    # d_grp = d_grp
    )
  
  toks <- unlist(tstrsplit(args[2], "[-.]"))
  fname <- paste0("data/sim07/sim07-", toks[4], "-", toks[5], "-", format(Sys.time(), "%Y%m%d-%H%M%S"), ".qs")
  
  log_info("Saving results file ", fname)
  
  qs::qsave(l, file = fname)
}

run_none_sim_07 <- function(){
  log_info("run_none_sim_07: Nothing doing here bud.")
}

main_sim_07 <- function(){
  funcname <- paste0(args[1], "()")
  log_info("Main, invoking ", funcname)
  eval(parse(text=funcname))
}

main_sim_07()


