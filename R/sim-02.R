source("./R/init.R")

# Command line arguments list scenario (true dose response),
# the number of simulations to run, the number of cores
# to use and the simulator to use.
args = commandArgs(trailingOnly=TRUE)

if (length(args)<1) {
  log_info("Setting default run method (does nothing)")
  args[1] = "run_none_sim_02"
  args[2] = "cfg-sim02-sc01-v01.yml"
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
  
  # init
  pop_spec <- get_pop_spec()
  sim_spec <- get_sim_spec()
  sim_spec$a0 <- qlogis(g_cfgsc$p_a0)
  sim_spec$m['l1'] <- qlogis(g_cfgsc$m_l1) - sim_spec$a0
  sim_spec$m['l2'] <- qlogis(g_cfgsc$m_l2) - sim_spec$a0
  sim_spec$b['erx'] <- g_cfgsc$b_erx
  sim_spec$b['r1'] <- g_cfgsc$b_r1
  sim_spec$b['r2'] <- g_cfgsc$b_r2
  sim_spec$b['edx'] <- g_cfgsc$b_edx
  sim_spec$b['r1d'] <- g_cfgsc$b_r1d
  sim_spec$b['r2d'] <- g_cfgsc$b_r2d
  sim_spec$b['efx'] <- g_cfgsc$b_efx
  sim_spec$b['f'] <- g_cfgsc$b_f

  # enrol times
  # easier to produce enrol times all in one hit rather than try to do them 
  # incrementally with a non-hom pp.
  loc_t0 <- get_enrol_time(max(g_cfgsc$N_pt))

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
  # inferiority probs
  pr_inf <- array(
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
  pr_ref_ni_trt <- array(
    NA, 
    dim = c(N_analys, length(g_fx)),
    dimnames = list(1:N_analys, g_fx)
  )
  # no need to track pr_fut as just taken as a treshold for pr_sup
  # i.e. fut iff pr_sup < 0.05 (or some other threshold)
  
  # decisions
  g_dec_type <- c("sup", "inf", "trt_ni_ref", "ref_ni_trt", "fut")
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
    f1 <- m1$sample(
      lsd$ld, iter_warmup = 1000, iter_sampling = 2000,
      parallel_chains = 2, chains = 2, refresh = 0, show_exceptions = F, 
      max_treedepth = 13)

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
      b_r2d = post_1$b_5 + post_1$b_6,
      # rif vs no-rif
      b_f = post_1$b_8
    ))
    
    m <- post_fx[, sapply(.SD, function(z){
      c(mu = mean(z), med = median(z),
        q_025 = unname(quantile(z, prob = 0.025)),
        q_975 = unname(quantile(z, prob = 0.975)))}),
      .SDcols = g_fx]
    post_smry_2[ii, ,] <- t(m)
    
  
    pr_sup[ii, ]  <- post_fx[, sapply(.SD, function(z){mean(z > log(g_cfgsc$delta_sup))}), .SDcols = g_fx]
    
    # assume inferiority is just defined as 1 - pr_sup[ii, ], i.e. Pr(z < log(g_cfgsc$delta_sup))
    # 12 weeks is coded as reference arm but want to determine if
    # 12 weeks is superior to 7 days (soc) in the two stage group and so we 
    # answer via the equivalent question namely the probability that 
    # 7 days is inferior to 12 weeks (assuming that we can define superiority 
    # and inferiority as symmetrical - I see no reason why we cannot)
    pr_inf[ii, ]  <- 1 - pr_sup[ii, ]
    
    # non-inferiority
    # The reference groups are dair, long and no-rif respectively for each of the
    # domains but there is an interest in doing ni comparisons such as 
    # "for the one-stage surgery cohort is 6-wks non-inferior to 12-wks"
    # which necessitates taking the mirror image of the log-or as it is 
    # originally parameterised.
    # Anyway, that's why I produce both of these summaries. There is an 
    # empirical example in Untitled1.R (see non_inferiority_2 func) if 
    # you need further detail.
    pr_trt_ni_ref[ii, ] <- post_fx[, sapply(.SD, function(z){mean(z > log(1/g_cfgsc$delta_ni))}), .SDcols = g_fx]
    pr_ref_ni_trt[ii, ] <- post_fx[, sapply(.SD, function(z){mean(-z > log(1/g_cfgsc$delta_ni))}), .SDcols = g_fx]

    # evaluate decisions
    decision[ii, , "sup"] <- pr_sup[ii, ] > g_cfgsc$thresh_sup
    decision[ii, , "inf"] <- pr_inf[ii, ] > g_cfgsc$thresh_sup
    decision[ii, , "trt_ni_ref"] <- pr_trt_ni_ref[ii, ] > g_cfgsc$thresh_non_inf
    decision[ii, , "ref_ni_trt"] <- pr_ref_ni_trt[ii, ] > g_cfgsc$thresh_non_inf
    # taken to imply negligible chance of ever being superior
    decision[ii, , "fut"] <- pr_sup[ii, ] < g_cfgsc$thresh_fut
    
    # Earlier decisions are retained - once superiority has been decided, 
    # we retain this conclusion irrespective of subsequent post probabilities.
    # This means that there could be inconsistency with a silo pr_sup and the 
    # decision reported.
    decision[1:ii, , "sup"] <- apply(decision[1:ii, , "sup", drop = F], 2, function(z){ cumsum(z) > 0 })
    decision[1:ii, , "inf"] <- apply(decision[1:ii, , "inf", drop = F], 2, function(z){ cumsum(z) > 0 })
    decision[1:ii, , "trt_ni_ref"] <- apply(decision[1:ii, , "trt_ni_ref", drop = F], 2, function(z){ cumsum(z) > 0 })
    decision[1:ii, , "ref_ni_trt"] <- apply(decision[1:ii, , "ref_ni_trt", drop = F], 2, function(z){ cumsum(z) > 0 })
    decision[1:ii, , "fut"] <- apply(decision[1:ii, , "fut", drop = F], 2, function(z){ cumsum(z) > 0 })

    # update pop_spec
    if(any(decision[ii, , "sup"])){
      # since there are only two treatments per cell, if a superiority decision 
      # is made then we have answered all the questions and we can stop 
      # enrolling into that cell. if there were more than two treatments then
      # we may take a different approach.
      
      # stop enrolling into the relevant domain, no further info 
      # accrues to the parameter estimates that led to the stopping decision
      if(decision[ii, "b_r", "sup"]){
        pop_spec$r_a$late['dair'] <- NA
        pop_spec$r_a$late['rev'] <- NA
      }
      if(decision[ii, "b_r1d", "sup"]){
        pop_spec$r_b$one['long'] <- NA
        pop_spec$r_b$one['short'] <- NA
      }
      if(decision[ii, "b_r2d", "sup"]){
        pop_spec$r_b$two['long'] <- NA
        pop_spec$r_b$two['short'] <- NA
      }
      if(decision[ii, "b_f", "sup"]){
        pop_spec$r_c['norif'] <- NA
        pop_spec$r_c['rif'] <- NA
      }
    }
    # same approach under inferiority. we stop enrolment into the 
    # respective domain. 
    if(any(decision[ii, , "inf"])){
      if(decision[ii, "b_r", "inf"]){
        pop_spec$r_a$late['dair'] <- NA
        pop_spec$r_a$late['rev'] <- NA
      }
      if(decision[ii, "b_r1d", "inf"]){
        pop_spec$r_b$one['long'] <- NA
        pop_spec$r_b$one['short'] <- NA
      }
      if(decision[ii, "b_r2d", "inf"]){
        pop_spec$r_b$two['long'] <- NA
        pop_spec$r_b$two['short'] <- NA
      }
      if(decision[ii, "b_f", "inf"]){
        pop_spec$r_c['norif'] <- NA
        pop_spec$r_c['rif'] <- NA
      }
    }
    
    # dont take any action for ni decisions. just keep going as is:
    # if(any(decision[ii, , "trt_ni_ref"])){
    # }
    # if(any(decision[ii, , "ref_ni_trt"])){
    # }
    
    # stop enrolling if futility decision made
    if(any(decision[ii, , "fut"])){
      if(decision[ii, "b_r", "fut"]){
        pop_spec$r_a$late['dair'] <- NA
        pop_spec$r_a$late['rev'] <- NA
      }
      if(decision[ii, "b_r1d", "fut"]){
        pop_spec$r_b$one['long'] <- NA
        pop_spec$r_b$one['short'] <- NA
      }
      if(decision[ii, "b_r2d", "fut"]){
        pop_spec$r_b$two['long'] <- NA
        pop_spec$r_b$two['short'] <- NA
      }
      if(decision[ii, "b_f", "fut"]){
        pop_spec$r_c['norif'] <- NA
        pop_spec$r_c['rif'] <- NA
      }
    }
    
    # stop sim if all questions answered
    
    if(all(decision[ii, , "sup"])){
      log_info("Stop trial all sup", ix)
      stop_enrol <- T  
    } 
    if(all(decision[ii, , "inf"])){
      log_info("Stop trial all inf", ix)
      stop_enrol <- T  
    } 
    if(all(decision[ii, , "fut"])){
      log_info("Stop trial all fut", ix)
      stop_enrol <- T  
    } 
    
    # next interim
    ii <- ii + 1
    
    if(ii > N_analys){
      stop_enrol <- T  
    }
  }
  
  # did we stop (for any reason) prior to the final interim?
  stop_at <- 5
  # review the first column of the superiority decision
  if(any(is.na(decision[, 1, "sup"]))){
    # interim where the stopping rule was met
    stop_at <- min(which(is.na(decision[, 1, "sup"]))) - 1
    
    if(stop_at < N_analys){
      message("stopped at analysis ", stop_at)

      log_info("Stopped at analysis ", stop_at, " filling all subsequent entries")
      decision[(stop_at+1):N_analys, , "sup"] <- decision[(stop_at), , "sup"]
      decision[(stop_at+1):N_analys, , "inf"] <- decision[(stop_at), , "inf"]
      decision[(stop_at+1):N_analys, , "trt_ni_ref"] <- decision[(stop_at), , "trt_ni_ref"]
      decision[(stop_at+1):N_analys, , "ref_ni_trt"] <- decision[(stop_at), , "ref_ni_trt"]
      decision[(stop_at+1):N_analys, , "fut"] <- decision[(stop_at), , "fut"]

    }
  }
  
  list(
    d_grp = d_all[, .(y = sum(y), .N), keyby = .(analys, l, er, r, srp, ed, d, ef, f)],
    post_smry_1 = post_smry_1,
    post_smry_2 = post_smry_2,
    decision = decision,
    pr_sup = pr_sup,
    pr_inf = pr_inf,
    pr_ref_ni_trt = pr_ref_ni_trt,
    pr_trt_ni_ref = pr_trt_ni_ref,
    # decision made on analysing the data from this analysis
    # e.g. if all treatments were superior after analysing the data at the 
    # third analysis then stop_at = 3
    stop_at = stop_at
  )
}

run_sim_02 <- function(){
  log_info(paste0(match.call()[[1]]))

  e = NULL
  log_info("Starting simulation")
  r <- parallel::mclapply(
    X=1:g_cfgsc$nsim, mc.cores = g_cfgsc$mc_cores, FUN=function(ix) {
    # X=1:20, mc.cores = g_cfgsc$mc_cores, FUN=function(ix) {
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
        sim = i, analys = as.integer(rownames(r[[i]]$pr_sup)),r[[i]]$pr_sup) 
      } )))
  
  d_pr_inf <- data.table(
    do.call(rbind, lapply(1:length(r), function(i){ 
      cbind(
        sim = i, analys = as.integer(rownames(r[[i]]$pr_inf)),r[[i]]$pr_inf) 
    } )))
  
  d_pr_trt_ni_ref <- data.table(
    do.call(rbind, lapply(1:length(r), function(i){ 
      cbind(
        sim = i, analys = as.integer(rownames(r[[i]]$pr_trt_ni_ref)),r[[i]]$pr_trt_ni_ref) 
    } )))
  
  d_pr_ref_ni_trt <- data.table(
    do.call(rbind, lapply(1:length(r), function(i){ 
      cbind(
        sim = i, analys = as.integer(rownames(r[[i]]$pr_ref_ni_trt)),r[[i]]$pr_ref_ni_trt) 
    } )))
  
  # decisions
  d_decision <- rbind(
    data.table(do.call(rbind, lapply(1:length(r), function(i){ 
      m <- r[[i]]$decision[, , "sup"]
      cbind(sim = i, analys = as.integer(rownames(m)), quant = "sup", m) 
    } ))), 
    
    data.table(do.call(rbind, lapply(1:length(r), function(i){ 
      m <- r[[i]]$decision[, , "inf"]
      cbind(sim = i, analys = as.integer(rownames(m)), quant = "inf", m) 
    } ))), 
    
    data.table(do.call(rbind, lapply(1:length(r), function(i){ 
      m <- r[[i]]$decision[, , "trt_ni_ref"]
      cbind(sim = i, analys = as.integer(rownames(m)), quant = "trt_ni_ref", m) 
    } ))),
    
    data.table(do.call(rbind, lapply(1:length(r), function(i){ 
      m <- r[[i]]$decision[, , "ref_ni_trt"]
      cbind(sim = i, analys = as.integer(rownames(m)), quant = "ref_ni_trt", m) 
    } ))),
    
    data.table(do.call(rbind, lapply(1:length(r), function(i){ 
      m <- r[[i]]$decision[, , "fut"]
      cbind(sim = i, analys = as.integer(rownames(m)), quant = "fut", m) 
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
    d_pr_inf = d_pr_inf,
    d_pr_trt_ni_ref = d_pr_trt_ni_ref,
    d_pr_ref_ni_trt = d_pr_ref_ni_trt,
    d_decision = d_decision,
    d_post_smry_1 = d_post_smry_1,
    d_post_smry_2 = d_post_smry_2,
    d_grp = d_grp
    )
  
  toks <- unlist(tstrsplit(args[2], "[-.]"))
  fname <- paste0("data/sim02-", toks[3], "-", toks[4], "-", format(Sys.time(), "%Y%m%d-%H%M%S"), ".qs")
  qs::qsave(l, file = fname)
}

run_none_sim_02 <- function(){
  log_info("run_none_sim_02: Nothing doing here bud.")
}

main_sim_02 <- function(){
  funcname <- paste0(args[1], "()")
  log_info("Main, invoking ", funcname)
  eval(parse(text=funcname))
}

main_sim_02()


