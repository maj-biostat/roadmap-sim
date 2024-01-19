source("./R/init.R")

# Command line arguments list scenario (true dose response),
# the number of simulations to run, the number of cores
# to use and the simulator to use.
args = commandArgs(trailingOnly=TRUE)

if (length(args)<1) {
  log_info("Setting default run method (does nothing)")
  args[1] = "run_none_sim_01"
  args[2] = "config-sc06.yml"
} else {
  log_info("Run method ", args[1])
  log_info("Scenario config ", args[2])
}

# Logs
f_log_sim <- file.path("./logs", "log-sim.txt")
log_appender(appender_file(f_log_sim))
# message(Sys.time(), " Log file initialised ", f_log)
log_info("*** START UP ***")
# log_threshold(TRACE)

f_cfgsc <- file.path("./etc", args[2])
g_cfgsc <- config::get(file = f_cfgsc)
stopifnot("Config is null" = !is.null(g_cfgsc))


ix <- 1
m1 <- cmdstanr::cmdstan_model("stan/model-sim-02.stan")

run_trial <- function(ix){
  
  log_info("Entered  run_trial for trial ", ix)
  
  # generate data

  log_trace("Create sim_spec", ix)
  sim_spec <- get_sim_spec()
  sim_spec$b_a_late["rev"] <- g_cfgsc$b_a_l_2
  sim_spec$b_a_chronic["two"] <- g_cfgsc$b_a_c_2
  sim_spec$b_b1_late_one["w12p1"] <- g_cfgsc$b_b1_l_2
  sim_spec$b_b2_late_two["w12p2"] <- g_cfgsc$b_b2_l_2
  sim_spec$b_b1_chronic_one["w12p1"] <- g_cfgsc$b_b1_c_2
  sim_spec$b_b2_chronic_two["w12p2"] <- g_cfgsc$b_b2_c_2
  sim_spec$b_c["rif"] <- g_cfgsc$b_c_2
  
  # just create data, nothing else
  log_debug("Create data", ix)
  ll <- get_trial_data(N = 2500, pop_spec = NULL, sim_spec = sim_spec)
  
  # analyses
  lsd <- get_stan_data(ll$d_i)

  f1 <- m1$sample(
    lsd$ld, iter_warmup = 1000, iter_sampling = 2000,
    parallel_chains = 2, chains = 2, refresh = 0, show_exceptions = F, 
    max_treedepth = 13)
  
  post_1 <- data.table(f1$draws(variables = c(
    "alpha", "gamma_c", 
    "b_a_l", 
    "b_b1_l", 
    "b_b2_l",
    "b_a_c",
    "b_b1_c", 
    "b_b2_c",
    "b_c"
  ), format = "matrix"))
  

  # test outcome
  cols <- names(post_1)
  cols <- gsub("[","_",cols,fixed = T)
  cols <- gsub("]","",cols,fixed = T)
  names(post_1) <- cols
  
  # effects
  effs <- c(
    # domain a, late (revision) and chronic (two-stage)
    "b_a_l_2", "b_a_c_2", 
    # domain b, (late/revision one stage pts)
    # wk12p1 (ref is wk6p1)
    "b_b1_l_2", 
    # wk12p2 (ref is day7p2)
    "b_b2_l_2", 
    # wk12p1 (ref is wk6p1)
    "b_b1_c_2", 
    # wk12p2 (ref is day7p2)
    "b_b2_c_2",
    # rif (ref is no-rif)
    "b_c_2")
  
  pr_sup <- post_1[, sapply(.SD, function(z){mean(z>0)}), .SDcols = effs]
  win <- pr_sup > g_cfgsc$d_sup
  
  # return results
  list(
    pr_sup = pr_sup, 
    win = win
  )
  
  
  
}

run_sim_01 <- function(){
  log_info(paste0(match.call()[[1]]))

  e = NULL
  log_info("Starting simulation")
  r <- parallel::mclapply(
    X=1:g_cfgsc$nsim, mc.cores = g_cfgsc$mc_cores, FUN=function(ix) {
      log_info("Simulation ", ix);
      ll <- tryCatch({
        run_trial(ix)
      },
      error=function(e) {
        log_info(" ERROR " , e);
        stop(paste0("Stopping with error ", e))
      })
      
      ll
    })
  
  d_pr_sup <- data.table(do.call(rbind, lapply(1:length(r), function(i){ r[[i]]$pr_sup } )))
  d_win <- data.table(do.call(rbind, lapply(1:length(r), function(i){ r[[i]]$win } )))
  
  l <- list(
    cfg = g_cfgsc,
    d_pr_sup = d_pr_sup, 
    d_win = d_win
    )
  
  fname <- paste0("data/sim01-", g_cfgsc$sc, "-", format(Sys.time(), "%Y%m%d-%H%M%S"), ".qs")
  qs::qsave(l, file = fname)
}

run_none_sim_01 <- function(){
  log_info("run_none_sim_01: Nothing doing here bud.")
}

main_sim_01 <- function(){
  funcname <- paste0(args[1], "()")
  log_info("Main, invoking ", funcname)
  eval(parse(text=funcname))
}

main_sim_01()