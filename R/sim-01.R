source("./R/init.R")

# Command line arguments list scenario (true dose response),
# the number of simulations to run, the number of cores
# to use and the simulator to use.
args = commandArgs(trailingOnly=TRUE)

if (length(args)<1) {
  log_info("Setting default run method (does nothing)")
  args[1] = "run_none_sim_01"
  args[2] = "cfg-sim01-sc01-v01.yml"
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
m1 <- cmdstanr::cmdstan_model("stan/model-sim-04.stan")

run_trial <- function(ix){
  
  log_info("Entered  run_trial for trial ", ix)
  
  # generate data

  log_trace("Create sim_spec", ix)
  sim_spec <- get_sim_spec()
  
  sim_spec$a0 <- qlogis(g_cfgsc$p_a0)
  sim_spec$m['l1'] <- qlogis(g_cfgsc$m_l1) - sim_spec$a0
  sim_spec$m['l2'] <- qlogis(g_cfgsc$m_l2) - sim_spec$a0
  
  sim_spec$b['erx-r0'] <- g_cfgsc$b_erx_r0
  sim_spec$b['erx-r1'] <- g_cfgsc$b_erx_r1
  sim_spec$b['erx-r2'] <- g_cfgsc$b_erx_r2
  
  sim_spec$b['r1'] <- g_cfgsc$b_r1
  sim_spec$b['r2'] <- g_cfgsc$b_r2
  # no longer required for the data generation process
  # sim_spec$b['edx'] <- g_cfgsc$b_edx
  sim_spec$b['r1d'] <- g_cfgsc$b_r1d
  sim_spec$b['r2d'] <- g_cfgsc$b_r2d
  sim_spec$b['efx'] <- g_cfgsc$b_efx
  sim_spec$b['f'] <- g_cfgsc$b_f

  
  # just create data, nothing else
  log_debug("Create data", ix)
  ll <- get_trial_data(N = g_cfgsc$N_pt, pop_spec = NULL, sim_spec = sim_spec)
  
  # analyses
  lsd <- get_stan_data(ll$d)

  f1 <- m1$sample(
    lsd$ld, iter_warmup = 1000, iter_sampling = 2000,
    parallel_chains = 2, chains = 2, refresh = 0, show_exceptions = F, 
    max_treedepth = 13)
  
  post_1 <- data.table(f1$draws(variables = c("a0", "m", "b"), format = "matrix"))
  post_1_cols <- names(post_1)
  post_1_cols <- gsub("[","_",post_1_cols,fixed = T)
  post_1_cols <- gsub("]","",post_1_cols,fixed = T)
  names(post_1) <- post_1_cols
  
  m <- post_1[, sapply(.SD, function(z){
    c(mu = mean(z), med = median(z),
      q_025 = unname(quantile(z, prob = 0.025)),
      q_975 = unname(quantile(z, prob = 0.975)))}),
    .SDcols = g_mod4_pars]
  post_smry_1 <- t(m)
  
  post_fx <- data.table(cbind(
    b_r = ll$d[er == 1 & r == 1, mean(srp1)] * post_1$b_4 + 
      ll$d[er == 1 & r == 1, mean(srp2)] * post_1$b_5,
    # short vs long (one-stage)
    b_r1d = post_1$b_6,
    # short vs long (two-stage)
    # now based on a single parameter
    b_r2d = post_1$b_7,
    # rif vs no-rif
    b_f = post_1$b_9
  ))
  
  m <- post_fx[, sapply(.SD, function(z){
    c(mu = mean(z), med = median(z),
      q_025 = unname(quantile(z, prob = 0.025)),
      q_975 = unname(quantile(z, prob = 0.975)))}),
    .SDcols = g_fx]
  post_smry_2 <- t(m)
  
  pr_sup <- post_fx[, sapply(.SD, function(z){mean(z > log(g_cfgsc$delta_sup))}), .SDcols = g_fx]
  pr_sup_fut <- post_fx[, sapply(.SD, function(z){mean(z > log(g_cfgsc$delta_sup_fut))}), .SDcols = g_fx]
  
  pr_trt_ni_ref <- post_fx[, sapply(.SD, function(z){mean(z > log(1/g_cfgsc$delta_ni))}), .SDcols = g_fx]
  pr_trt_ni_ref_fut <- post_fx[, sapply(.SD, function(z){mean(z > log(1/g_cfgsc$delta_ni_fut))}), .SDcols = g_fx]
 
  sup <- pr_sup > g_cfgsc$thresh_sup
  trt_ni_ref <- pr_trt_ni_ref > g_cfgsc$thresh_non_inf
  fut_sup <- pr_sup_fut < g_cfgsc$thresh_fut_sup
  fut_trt_ni_ref  <- pr_trt_ni_ref < g_cfgsc$thresh_fut_ni
 
  # return results
  list(
    d_grp = ll$d[, .(y = sum(y), .N), keyby = .(l, er, r, srp, ed, d, ef, f)],
    pr_sup = pr_sup,
    pr_sup_fut = pr_sup_fut,
    pr_trt_ni_ref = pr_trt_ni_ref,
    pr_trt_ni_ref_fut = pr_trt_ni_ref_fut,
    sup = sup, 
    trt_ni_ref = trt_ni_ref,
    fut_sup = fut_sup,
    fut_trt_ni_ref = fut_trt_ni_ref,
    post_smry_1 = post_smry_1,
    post_smry_2 = post_smry_2
  )
  
  
  
}

run_sim_01 <- function(){
  log_info(paste0(match.call()[[1]]))

  e = NULL
  log_info("Starting simulation")   #
  r <- parallel::mclapply(
    X=1:g_cfgsc$nsim, mc.cores = g_cfgsc$mc_cores, FUN=function(ix) {
    # X=1:5, mc.cores = g_cfgsc$mc_cores, FUN=function(ix) {
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
  d_pr_sup_fut <- data.table(do.call(rbind, lapply(1:length(r), function(i){ r[[i]]$fut_sup } )))
  d_pr_trt_ni_ref <- data.table(do.call(rbind, lapply(1:length(r), function(i){ r[[i]]$pr_trt_ni_ref } )))
  d_pr_trt_ni_ref_fut <- data.table(do.call(rbind, lapply(1:length(r), function(i){ r[[i]]$pr_trt_ni_ref_fut } )))
  # no pr_fut

  d_sup <- data.table(do.call(rbind, lapply(1:length(r), function(i){ r[[i]]$sup } )))
  d_trt_ni_ref <- data.table(do.call(rbind, lapply(1:length(r), function(i){ r[[i]]$trt_ni_ref } )))
  d_fut_sup <- data.table(do.call(rbind, lapply(1:length(r), function(i){ r[[i]]$fut_sup } )))
  d_fut_trt_ni_ref <- data.table(do.call(rbind, lapply(1:length(r), function(i){ r[[i]]$fut_trt_ni_ref } )))
  
  # posterior parameter summary
  d_post_smry_2 <- data.table(do.call(rbind, lapply(1:length(r), function(i){ 
    m <- data.table(
      sim = i, par = rownames(r[[i]]$post_smry_2), r[[i]]$post_smry_2)
    m
  } )))
  
  # data from each simulated trial
  d_grp <- rbindlist(lapply(1:length(r), function(i){ 
    r[[i]]$d_grp
  } ), idcol = "sim")
  
  l <- list(
    cfg = g_cfgsc,
    d_pr_sup = d_pr_sup, 
    d_pr_sup_fut = d_pr_sup_fut, 
    d_pr_trt_ni_ref = d_pr_trt_ni_ref, 
    
    d_sup = d_sup,
    d_trt_ni_ref = d_trt_ni_ref,
    d_fut_sup = d_fut_sup,
    d_fut_trt_ni_ref = d_fut_trt_ni_ref,
    d_post_smry_2 = d_post_smry_2,
    d_grp = d_grp
    )
  
  toks <- unlist(tstrsplit(args[2], "[-.]"))
  fname <- paste0("data/sim01-", toks[3], "-", toks[4], "-", format(Sys.time(), "%Y%m%d-%H%M%S"), ".qs")
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


