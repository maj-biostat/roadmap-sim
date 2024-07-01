# Uses posterior predictive distribution to evaluate the likelihood of a 
# trial conclusion in the next X pts. Does not evaluate at the max sample
# size so probably of limited value. Nevertheless, implementation is basically
# in place and so have left as is.

source("./R/init.R")

# Command line arguments list scenario (true dose response),
# the number of simulations to run, the number of cores
# to use and the simulator to use.
args = commandArgs(trailingOnly=TRUE)

# Load cfg based on cmd line args.
if (length(args)<1) {
  log_info("Setting default run method (does nothing)")
  args[1] = "run_none_sim_04"
  args[2] = "cfg-sim04-sc01-v05.yml"
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
m1 <- cmdstanr::cmdstan_model("stan/model-sim-04.stan")



# Main trial loop.
run_trial <- function(ix){
  
  log_info("Entered  run_trial for trial ", ix)
  
  # initialise simulation parameters
  # get_pop_spec and get_sim_spec are defined in the roadmap.data package
  # These functions allow us to simulate data according to our spec and also 
  # to switch treatments on and off. 
  pop_spec <- roadmap.data::get_pop_spec()
  sim_spec <- roadmap.data::get_sim_spec()
  # Intercept
  sim_spec$a0 <- qlogis(g_cfgsc$p_a0)
  # Strata
  sim_spec$m['l1'] <- qlogis(g_cfgsc$m_l1) - sim_spec$a0
  sim_spec$m['l2'] <- qlogis(g_cfgsc$m_l2) - sim_spec$a0
  # Surgery 
  sim_spec$b['erx'] <- g_cfgsc$b_erx
  sim_spec$b['r1'] <- g_cfgsc$b_r1
  sim_spec$b['r2'] <- g_cfgsc$b_r2
  # Duration
  sim_spec$b['edx'] <- g_cfgsc$b_edx
  sim_spec$b['r1d'] <- g_cfgsc$b_r1d
  sim_spec$b['r2d'] <- g_cfgsc$b_r2d
  # Choice
  sim_spec$b['efx'] <- g_cfgsc$b_efx
  sim_spec$b['f'] <- g_cfgsc$b_f

  # enrol times
  # easier to produce enrol times all in one hit rather than try to do them 
  # incrementally with a non-hom pp.
  loc_t0 <- roadmap.data::get_enrol_time(max(g_cfgsc$N_pt))

  # loop controls
  stop_enrol <- FALSE
  ii <- 1 # interim number
  N_analys <- length(g_cfgsc$N_pt)

  # tracking/output quantities
  smry_lab <- c("mu", "med", "q_025", "q_975")
  # all parameters
  post_smry_1 <- array(
    NA, 
    dim = c(N_analys, length(g_mod4_pars), length(smry_lab)),
    dimnames = list(1:N_analys, g_mod4_pars, smry_lab)
  )
  # primary quantities of interest
  post_smry_2 <- array(
    NA, 
    dim = c(N_analys, length(g_fx), length(smry_lab)),
    dimnames = list(1:N_analys, g_fx, smry_lab)
  )
  # superiority probs
  pr_sup <- array(
    NA, 
    dim = c(N_analys, length(g_fx)),
    dimnames = list(1:N_analys, g_fx)
  )
  
  # non-inferiority probs
  pr_trt_ni_ref <- array(
    NA, 
    dim = c(N_analys, length(g_fx)),
    dimnames = list(1:N_analys, g_fx)
  )
  
  # decisions
  g_dec_type <- c("sup", 
                  "trt_ni_ref", 
                  "fut_sup", 
                  "fut_trt_ni_ref"
                  )

  decision <- array(
    NA,
    dim = c(N_analys, length(g_fx), length(g_dec_type)),
    dimnames = list(
      1:N_analys, g_fx, g_dec_type)
  )

  # store all simulated trial pt data
  d_all <- data.table()
  
  while(!stop_enrol){
  
    # next chunk of data on pts.
    if(ii == 1){
      N_c <- g_cfgsc$N_pt[ii]
      # starting pt index in data
      is <- 1
      ie <- is + N_c - 1
    } else {
      N_c <- g_cfgsc$N_pt[ii] - g_cfgsc$N_pt[ii-1]
      is <- g_cfgsc$N_pt[ii-1] + 1
      ie <- is + N_c
    }
    
    # next block of pts that have reached 12 months post rand
    l_new <- get_trial_data(
      N = N_c, 
      pop_spec = pop_spec, 
      sim_spec = sim_spec, 
      idx_s = is, entry_times = FALSE)
  
    # add in a counter for tracking analysis
    l_new$d[, analys := ii]

    # combine the existing and new data
    d_all <- rbind(d_all, l_new$d)
    d_all[, t0 := loc_t0[1:.N]]
    
    # create stand data based on enrolled pt
    lsd <- get_stan_data(d_all)

    # fit model
    # f1 <- m1$sample(
    #   lsd$ld, iter_warmup = 1000, iter_sampling = 2000,
    #   parallel_chains = 2, chains = 2, refresh = 0, show_exceptions = F, 
    #   max_treedepth = 13)
    
    snk <- capture.output(
      f1 <- m1$pathfinder(lsd$ld, num_paths=10, single_path_draws=40,
                             history_size=50, max_lbfgs_iters=100,
                             refresh = 0)
    )

    # extract post
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
    post_smry_1[ii, ,] <- t(m)
   
    # primary quantities of interest - recheck these
    post_fx <- data.table(cbind(
      # revision vs dair
      # b_r = post_1$b_2 + d_all[, mean(srp2)] * post_1$b_2,
      # need to think through whether this weighting should be from the whole
      # pop or just those revealed to the surgery domain. think it is the 
      # latter as implemented below.  
      b_r = post_1$b_2 + d_all[er == 1, mean(srp2)] * post_1$b_2,
      # short vs long (one-stage)
      b_r1d = post_1$b_5,
      # short vs long (two-stage)
      b_r2d = post_1$b_6,
      # rif vs no-rif
      b_f = post_1$b_8
    ))
    
    m <- post_fx[, sapply(.SD, function(z){
      c(mu = mean(z), med = median(z),
        q_025 = unname(quantile(z, prob = 0.025)),
        q_975 = unname(quantile(z, prob = 0.975)))}),
      .SDcols = g_fx]
    post_smry_2[ii, ,] <- t(m)
    
    # 9.	Hypothesis testing – the Core protocol, Stats appendix and DSA all 
    # describe superiority (OR<1), non-inf (OR>1.2) and futility – see below. 
    # The hypotheses tested, as outlines in the DSAs will be: 
    
    # a.	A – Surg Domain – Early Silo – Revision is superior to DAIR
    # b.	B – AB Choice – All silos combined – Rifampicin is superior to no-rifampicin
    
    # are these to be silo specific
    # c.	AB Duration – One stage revision – 6 weeks is non-inferior to 12 weeks (assumed to be soc)
    # d.	AB Duration – Two-stage revision – 12 weeks is superior to 7 days (soc)
    
    pr_sup[ii, ]  <- post_fx[, sapply(.SD, function(z){mean(z > log(g_cfgsc$delta_sup))}), .SDcols = g_fx]
    pr_trt_ni_ref[ii, ] <- post_fx[, sapply(.SD, function(z){mean(z > log(1/g_cfgsc$delta_ni))}), .SDcols = g_fx]

    # evaluate decisions
    decision[ii, , "sup"] <- pr_sup[ii, ] > g_cfgsc$thresh_sup
    decision[ii, , "trt_ni_ref"] <- pr_trt_ni_ref[ii, ] > g_cfgsc$thresh_non_inf

    # posterior predictive

    jj <- 1
    idx_draw <- sort(sample(1:nrow(post_1), g_cfgsc$nsimpred, replace = F))
    fut_sup <- array(
      NA, dim = c(g_cfgsc$nsimpred, length(g_fx)), dimnames = list(1:g_cfgsc$nsimpred, g_fx))
    fut_trt_ni_ref <- array(
      NA, dim = c(g_cfgsc$nsimpred, length(g_fx)), dimnames = list(1:g_cfgsc$nsimpred, g_fx))
    
    for(jj in 1:g_cfgsc$nsimpred){
      
      log_debug("Predictive run ", jj, " for trial ", ix)
      
      d_all_pp <- copy(d_all)
      d_all_pp[, analys := NULL]
      
      sim_spec_pp <- list()
      sim_spec_pp$a0 <- post_1$a0[idx_draw[jj]]
      sim_spec_pp$m["l1"] <- post_1$m_1[idx_draw[jj]]
      sim_spec_pp$m["l2"] <- post_1$m_2[idx_draw[jj]]
      sim_spec_pp$b["erx"] <- post_1$b_1[idx_draw[jj]]
      sim_spec_pp$b["r1"] <- post_1$b_2[idx_draw[jj]]
      sim_spec_pp$b["r2"] <- post_1$b_3[idx_draw[jj]]
      sim_spec_pp$b["edx"] <- post_1$b_4[idx_draw[jj]]
      sim_spec_pp$b["r1d"] <- post_1$b_5[idx_draw[jj]]
      sim_spec_pp$b["r2d"] <- post_1$b_6[idx_draw[jj]]
      sim_spec_pp$b["efx"] <- post_1$b_7[idx_draw[jj]]
      sim_spec_pp$b["f"] <- post_1$b_8[idx_draw[jj]]
      
      l_pp <- get_trial_data(
        N = 500, 
        pop_spec = pop_spec, 
        sim_spec = sim_spec_pp, 
        idx_s = max(d_all_pp$id) + 1, entry_times = FALSE)

      # combine the existing and new data
      d_all_pp <- rbind(d_all_pp, l_pp$d)
      
      # create stand data based on enrolled pt
      lsd_pp <- get_stan_data(d_all_pp)
      
      # fit model
      snk <- capture.output(
        f1_pp <- m1$pathfinder(lsd_pp$ld, num_paths=10, single_path_draws=40,
                             history_size=50, max_lbfgs_iters=100,
                             refresh = 0)
      )

      post_1_pp <- data.table(f1_pp$draws(variables = c("a0", "m", "b"), format = "matrix"))
      post_1_cols <- names(post_1_pp)
      post_1_cols <- gsub("[","_",post_1_cols,fixed = T)
      post_1_cols <- gsub("]","",post_1_cols,fixed = T)
      names(post_1_pp) <- post_1_cols
      
      post_fx_pp <- data.table(cbind(  
        b_r = post_1_pp$b_2 + d_all_pp[er == 1, mean(srp2)] * post_1_pp$b_2,
        # short vs long (one-stage)
        b_r1d = post_1_pp$b_5,
        # short vs long (two-stage)
        b_r2d = post_1_pp$b_6,
        # rif vs no-rif
        b_f = post_1_pp$b_8
      ))
      
      pr_sup_pp  <- post_fx_pp[, sapply(.SD, function(z){mean(z > log(1/g_cfgsc$delta_sup_fut))}), .SDcols = g_fx]
      pr_trt_ni_ref_pp <- post_fx[, sapply(.SD, function(z){mean(z > log(1/g_cfgsc$delta_ni_fut))}), .SDcols = g_fx]
      
      # if the probability of achieving the criteria is small then we 
      # deem the study futile
      fut_sup[jj, ] <- pr_sup_pp < g_cfgsc$thresh_fut_sup
      fut_trt_ni_ref[jj, ] <- pr_trt_ni_ref_pp < g_cfgsc$thresh_fut_ni
      
    }
    
    decision[1:ii, , "fut_sup"] <- colMeans(fut_sup) > g_cfgsc$thresh_fut_sup
    decision[1:ii, , "fut_trt_ni_ref"] <- colMeans(fut_trt_ni_ref) > g_cfgsc$thresh_fut_ni
    
    
    # Earlier decisions are retained - once superiority has been decided, 
    # we retain this conclusion irrespective of subsequent post probabilities.
    # This means that there could be inconsistency with a silo pr_sup and the 
    # decision reported.
    decision[1:ii, , "sup"] <- apply(decision[1:ii, , "sup", drop = F], 2, function(z){ cumsum(z) > 0 })
    decision[1:ii, , "trt_ni_ref"] <- apply(decision[1:ii, , "trt_ni_ref", drop = F], 2, function(z){ cumsum(z) > 0 })
    decision[1:ii, , "fut_sup"] <- apply(decision[1:ii, , "fut_sup", drop = F], 2, function(z){ cumsum(z) > 0 })
    decision[1:ii, , "fut_trt_ni_ref"] <- apply(decision[1:ii, , "fut_trt_ni_ref", drop = F], 2, function(z){ cumsum(z) > 0 })

    # update pop_spec
    # if we set the alloc probs to NA, the data generation does not enter pt
    # into the randomised comparison for that domain
    
    if(any(decision[ii, , "sup"])){
      
      # since there are only two treatments per cell, if a superiority decision 
      # is made then we have answered all the questions and we can stop 
      # enrolling into that cell. if there were more than two treatments then
      # we would probably need to take a different approach.
      
      # the only surgical intervention is for revision versus dair in the late cohort
      # so b_r only relates to that group, so it is a silo specific decision
      # note that r_a$early and r_a$chronic are NA from the start
      if(decision[ii, "b_r", "sup"]){
        pop_spec$r_a$late['dair'] <- NA
        pop_spec$r_a$late['rev'] <- NA
      }
      # choice domain relates to all silo so decision on b_f impacts all cohorts
      if(decision[ii, "b_f", "sup"]){
        pop_spec$r_c['norif'] <- NA
        pop_spec$r_c['rif'] <- NA
      }
      # choice domain relates to all silo so decision on b_f impacts all cohorts
      if(decision[ii, "b_r2d", "sup"]){
        pop_spec$r_b$two['long'] <- NA
        pop_spec$r_b$two['short'] <- NA
      }
    }
    
    # ni only applies to duration. if short is ni to long then stop enrolment
    # for this randomised comparison
    if(any(decision[ii, , "trt_ni_ref"])){
      
      if(decision[ii, "b_r1d", "trt_ni_ref"]){
        pop_spec$r_b$one['long'] <- NA
        pop_spec$r_b$one['short'] <- NA
      }
    }
    
    # stop enrolling if futility decision wrt superiority
    if(any(decision[ii, , "fut_sup"])){
      
      if(decision[ii, "b_r", "fut_sup"]){
        # for late stage silo only
        pop_spec$r_a$late['dair'] <- NA
        pop_spec$r_a$late['rev'] <- NA
      }
      if(decision[ii, "b_f", "fut_sup"]){
        pop_spec$r_c['norif'] <- NA
        pop_spec$r_c['rif'] <- NA
      }
      if(decision[ii, "b_r2d", "fut_sup"]){
        pop_spec$r_b$two['long'] <- NA
        pop_spec$r_b$two['short'] <- NA
      }
    }
    
    if(any(decision[ii, , "fut_trt_ni_ref"])){
      if(decision[ii, "b_r1d", "fut_trt_ni_ref"]){
        pop_spec$r_b$one['long'] <- NA
        pop_spec$r_b$one['short'] <- NA
      }
    }
    
    # stop sim if all questions answered - for the above setup this is 
    # redundant since these can no longer happen in all domains
    
    
    if(
      # if b_r is superior or futile
      (decision[ii, "b_r", "sup"] | decision[ii, "b_r", "fut_sup"]) &
      # if b_f is superior or futile
      (decision[ii, "b_f", "sup"] | decision[ii, "b_f", "fut_sup"]) &
      # if one stage short is non-inferior or futile for non-inferiority
      (decision[ii, "b_r1d", "trt_ni_ref"] | decision[ii, "b_r1d", "fut_trt_ni_ref"]) &
      # if two stage short is inferior or futile for inferiority
      (decision[ii, "b_r2d", "sup"] | decision[ii, "b_r2d", "fut_sup"])
    ){
      log_info("Stop trial ", ix, " as all questions addressed.")
      stop_enrol <- T  
    } 
    
    # next interim
    ii <- ii + 1
    
    if(ii > N_analys){
      stop_enrol <- T  
    }
  }
  
  # did we stop (for any reason) prior to the final interim?
  stop_at <- N_analys
  # any na's in any of the decision array rows means we stopped early:
  if(any(is.na(decision[, 1, "sup"]))){
    # interim where the stopping rule was met
    stop_at <- min(which(is.na(decision[, 1, "sup"]))) - 1
    
    if(stop_at < N_analys){
      message("stopped at analysis ", stop_at)

      log_info("Stopped at analysis ", stop_at, " filling all subsequent entries")
      decision[(stop_at+1):N_analys, , "sup"] <- decision[rep(stop_at, N_analys-stop_at), , "sup"]
      decision[(stop_at+1):N_analys, , "trt_ni_ref"] <- decision[rep(stop_at, N_analys-stop_at), , "trt_ni_ref"]
      decision[(stop_at+1):N_analys, , "fut_sup"] <- decision[rep(stop_at, N_analys-stop_at), , "fut_sup"]
      decision[(stop_at+1):N_analys, , "fut_trt_ni_ref"] <- decision[rep(stop_at, N_analys-stop_at), , "fut_trt_ni_ref"]

    }
  }
  
  list(
    d_grp = d_all[, .(y = sum(y), .N), keyby = .(analys, l, er, r, srp, ed, d, ef, f)],
    post_smry_1 = post_smry_1,
    post_smry_2 = post_smry_2,
    decision = decision,
    pr_sup = pr_sup,
    pr_trt_ni_ref = pr_trt_ni_ref,
    # decision made on analysing the data from this analysis
    # e.g. if all treatments were superior after analysing the data at the 
    # third analysis then stop_at = 3
    stop_at = stop_at
  )
}

run_sim_04 <- function(){
  log_info(paste0(match.call()[[1]]))

  e = NULL
  log_info("Starting simulation")
  r <- parallel::mclapply(
    X=1:g_cfgsc$nsim, mc.cores = g_cfgsc$mc_cores, FUN=function(ix) {
    # X=1:10, mc.cores = g_cfgsc$mc_cores, FUN=function(ix) {
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
  
  # put results into data.tables
  d_pr_sup <- data.table(
    do.call(rbind, lapply(1:length(r), function(i){ 
      cbind(
        sim = i, analys = as.integer(rownames(r[[i]]$pr_sup)), r[[i]]$pr_sup) 
      } )))

  d_pr_trt_ni_ref <- data.table(
    do.call(rbind, lapply(1:length(r), function(i){ 
      cbind(
        sim = i, analys = as.integer(rownames(r[[i]]$pr_trt_ni_ref)),r[[i]]$pr_trt_ni_ref) 
    } )))

  # decisions
  d_decision <- rbind(
    
    data.table(do.call(rbind, lapply(1:length(r), function(i){ 
      m <- r[[i]]$decision[, , "sup"]
      cbind(sim = i, analys = as.integer(rownames(m)), quant = "sup", m) 
    } ))), 
    
    data.table(do.call(rbind, lapply(1:length(r), function(i){ 
      m <- r[[i]]$decision[, , "trt_ni_ref"]
      cbind(sim = i, analys = as.integer(rownames(m)), quant = "trt_ni_ref", m) 
    } ))),
    
    data.table(do.call(rbind, lapply(1:length(r), function(i){ 
      m <- r[[i]]$decision[, , "fut_sup"]
      cbind(sim = i, analys = as.integer(rownames(m)), quant = "fut_sup", m) 
    } ))),
    
    data.table(do.call(rbind, lapply(1:length(r), function(i){ 
      m <- r[[i]]$decision[, , "fut_trt_ni_ref"]
      cbind(sim = i, analys = as.integer(rownames(m)), quant = "fut_trt_ni_ref", m) 
    } )))
    
  )
  d_decision[, `:=`(sim = as.integer(sim), analys = as.integer(analys),
                b_r = as.logical(b_r), 
                b_r1d = as.logical(b_r1d), 
                b_r2d = as.logical(b_r2d),  
                b_f = as.logical(b_f) 
                )]
  
  # posterior parameter summary
  d_post_smry_1 <- data.table(do.call(rbind, lapply(1:length(r), function(i){ 
    # build data.table(analys, smry, a0, m_1, ..., b_8)
    m <- data.table(
      analys = rep(
        # dim gives numer of interim, number of params, number of smry stats
        1:(dim(r[[i]]$post_smry_1)[1]), 
        len = (dim(r[[i]]$post_smry_1)[1]) * dim(r[[i]]$post_smry_1)[3]),
      # summary quantity (mean, med etc)
      smry = rep(dimnames(r[[i]]$post_smry_1)[[3]], 
                 each = dim(r[[i]]$post_smry_1)[1]),
      # summary value
      apply(r[[i]]$post_smry_1,2,c)
    )
    # convert to long
    m <- melt(m, id.vars = c("analys", "smry"), variable.name = "par")
    m <- dcast(m, analys + par ~ smry, value.var = "value")
    m <- cbind(sim = i, m)
    m
  } )))
  d_post_smry_1[, analys := as.integer(as.character(analys))]
  
  # posterior parameter summary
  d_post_smry_2 <- data.table(do.call(rbind, lapply(1:length(r), function(i){ 
    m <- data.table(
      analys = rep(
        # dim gives numer of interim, number of params, number of smry stats
        1:(dim(r[[i]]$post_smry_2)[1]), 
        len = (dim(r[[i]]$post_smry_2)[1]) * dim(r[[i]]$post_smry_2)[3]),
      # summary quantity (mean, med etc)
      smry = rep(dimnames(r[[i]]$post_smry_2)[[3]], 
                 each = dim(r[[i]]$post_smry_2)[1]),
      # summary value
      apply(r[[i]]$post_smry_2,2,c)
    )
    # convert to long
    m <- melt(m, id.vars = c("analys", "smry"), variable.name = "par")
    m <- dcast(m, analys + par ~ smry, value.var = "value")
    m <- cbind(sim = i, m)
    m
  } )))
  d_post_smry_2[, analys := as.integer(as.character(analys))]
  
  # data from each simulated trial
  d_grp <- rbindlist(lapply(1:length(r), function(i){ 
    r[[i]]$d_grp
  } ), idcol = "sim")


  l <- list(
    cfg = g_cfgsc,
    d_pr_sup = d_pr_sup, 
    d_pr_trt_ni_ref = d_pr_trt_ni_ref,
    d_decision = d_decision,
    d_post_smry_1 = d_post_smry_1,
    d_post_smry_2 = d_post_smry_2,
    d_grp = d_grp
    )
  
  toks <- unlist(tstrsplit(args[2], "[-.]"))
  fname <- paste0("data/sim04-", toks[3], "-", toks[4], "-", format(Sys.time(), "%Y%m%d-%H%M%S"), ".qs")
  qs::qsave(l, file = fname)
}

run_none_sim_04 <- function(){
  log_info("run_none_sim_04: Nothing doing here bud.")
}

main_sim_04 <- function(){
  funcname <- paste0(args[1], "()")
  log_info("Main, invoking ", funcname)
  eval(parse(text=funcname))
}

main_sim_04()


