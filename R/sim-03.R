# Sequential design based on truncated set of decisions that are of primary
# interest to the researchers.
# Uses posterior at current interim to coordinate decisions based on 
# pre-specified thresholds. 

source("./R/init.R")

# Command line arguments list scenario (true dose response),
# the number of simulations to run, the number of cores
# to use and the simulator to use.
args = commandArgs(trailingOnly=TRUE)

# Load cfg based on cmd line args.
if (length(args)<1) {
  log_info("Setting default run method (does nothing)")
  args[1] = "run_none_sim_03"
  args[2] = "cfg-sim03-sc01-v01.yml"
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

run_trial <- function(ix){
  
  log_info("Entered  run_trial for trial ", ix)
  
  # initialise simulation parameters
  # get_pop_spec and get_sim_spec are defined in the roadmap.data package
  # These functions allow us to simulate data according to our spec and also 
  # to switch treatments on and off. 
  pop_spec <- roadmap.data::get_pop_spec()
  sim_spec <- roadmap.data::get_sim_spec()
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

  # enrol times
  # easier to produce big list of enrol times all in one hit rather than do this
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
  # probability of futility wrt superiority decision
  pr_sup_fut  <- array(
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
  # non-inferiority probs
  pr_trt_ni_ref_fut <- array(
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
    
    # next block of pts that are assumed to have reached 12 months post rand
    l_new <- roadmap.data::get_trial_data(
      N = N_c, 
      pop_spec = pop_spec, 
      sim_spec = sim_spec, 
      idx_s = is, entry_times = FALSE)
  
    # add in a counter for tracking analysis
    l_new$d[, analys := ii]

    # combine the existing and new data and add enrolment time
    d_all <- rbind(d_all, l_new$d)
    d_all[, t0 := loc_t0[1:.N]]
    
    # convert all data into form suitable for stan
    lsd <- roadmap.data::get_stan_data(d_all)

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
      # Query whether this weighting based on observed srp2 (indicator for 
      # two-stage performed) should be from the whole pop or just those 
      # revealed for the surgery domain. think it is the latter as currently
      # implemented below.  The weight calculation is restricted to those 
      # that receive revision.
      b_r = d_all[er == 1 & r == 1, mean(srp1)] * post_1$b_4 + 
        d_all[er == 1 & r == 1, mean(srp2)] * post_1$b_5,
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
    post_smry_2[ii, ,] <- t(m)

    # The hypotheses tested, as outlines in the DSAs will be: 
    
    # a.	A – Surg Domain – Late Silo – Revision is superior to DAIR
    # b.	B – AB Choice – All silos combined – Rifampicin is superior to no-rifampicin
    
    # are these to be silo specific
    # c.	AB Duration – One stage revision – 6 weeks is non-inferior to 12 weeks (assumed to be soc)
    # d.	AB Duration – Two-stage revision – 12 weeks is superior to 7 days (soc)
    
    pr_sup[ii, ]  <- post_fx[, sapply(.SD, function(z){mean(z > log(g_cfgsc$delta_sup))}), .SDcols = g_fx]
    
    # just allows me to put a different increment in relative to the superiority assessment
    pr_sup_fut[ii, ]  <- post_fx[, sapply(.SD, function(z){mean(z > log(g_cfgsc$delta_sup_fut))}), .SDcols = g_fx]
    
    # non-inferiority
    # all treatment terms coded such that NI relates directly to comparison of interest
    # e.g. for one stage rev, is short non-inferior to long (soc)
    pr_trt_ni_ref[ii, ] <- post_fx[, sapply(.SD, function(z){mean(z > log(1/g_cfgsc$delta_ni))}), .SDcols = g_fx]
    
    # allows for a different increment in relative to the ni assessment
    pr_trt_ni_ref_fut[ii, ] <- post_fx[, sapply(.SD, function(z){mean(z > log(1/g_cfgsc$delta_ni_fut))}), .SDcols = g_fx]

    # evaluate decisions
    decision[ii, , "sup"] <- pr_sup[ii, ] > g_cfgsc$thresh_sup
    decision[ii, , "trt_ni_ref"] <- pr_trt_ni_ref[ii, ] > g_cfgsc$thresh_non_inf
    # taken to imply negligible chance of being superior
    decision[ii, , "fut_sup"] <- pr_sup_fut[ii, ] < g_cfgsc$thresh_fut_sup
    # taken to imply negligible chance of being ni 
    # this is not symmetrical with the approach taken for superiority.
    decision[ii, , "fut_trt_ni_ref"] <- pr_trt_ni_ref_fut[ii, ] < g_cfgsc$thresh_fut_ni
    
    # Earlier decisions are retained - once superiority has been decided, 
    # we retain this conclusion irrespective of subsequent post probabilities.
    # The following simply overwrites any decision reversal.
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
      # duration domain relates to all silo so decision on b_r2d impacts all cohorts
      if(decision[ii, "b_r2d", "sup"]){
        pop_spec$r_b$two['long'] <- NA
        pop_spec$r_b$two['short'] <- NA
      }
    }
    
    # not considering inferiority
    
    # ni only applies to duration. if short is ni to long then stop enrolment
    # for this randomised comparison
    if(any(decision[ii, , "trt_ni_ref"])){
      
      if(decision[ii, "b_r1d", "trt_ni_ref"]){
        pop_spec$r_b$one['long'] <- NA
        pop_spec$r_b$one['short'] <- NA
      }
    }
    
    # stop enrolling if futile wrt superiority decision
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
    
    # stop sim if all questions answered - 
    
    if(
      # if b_r is superior or futile
      (decision[ii, "b_r", "sup"] | decision[ii, "b_r", "fut_sup"]) &
      # if b_f is superior or futile
      (decision[ii, "b_f", "sup"] | decision[ii, "b_f", "fut_sup"]) &
      # if one stage short is non-inferior or futile for non-inferiority
      (decision[ii, "b_r1d", "trt_ni_ref"] | decision[ii, "b_r1d", "fut_trt_ni_ref"]) &
      # if two stage short is superior or futile for superior
      (decision[ii, "b_r2d", "sup"] | decision[ii, "b_r2d", "fut_sup"])
    ){
      log_info("Stop trial all questions addressed ", ix)
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
    pr_sup_fut = pr_sup_fut, 
    pr_trt_ni_ref = pr_trt_ni_ref,
    pr_trt_ni_ref_fut = pr_trt_ni_ref_fut,
    # decision made on analysing the data from this analysis
    # e.g. if all treatments were superior after analysing the data at the 
    # third analysis then stop_at = 3
    stop_at = stop_at
  )
}

run_sim_03 <- function(){
  log_info(paste0(match.call()[[1]]))

  e = NULL
  log_info("Starting simulation")
  r <- parallel::mclapply(
    X=1:g_cfgsc$nsim, mc.cores = g_cfgsc$mc_cores, FUN=function(ix) {
    # X=1:100, mc.cores = g_cfgsc$mc_cores, FUN=function(ix) {
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
  
  d_pr_sup_fut <- data.table(
    do.call(rbind, lapply(1:length(r), function(i){ 
      cbind(
        sim = i, analys = as.integer(rownames(r[[i]]$pr_sup_fut)), r[[i]]$pr_sup_fut) 
    } )))

  d_pr_trt_ni_ref <- data.table(
    do.call(rbind, lapply(1:length(r), function(i){ 
      cbind(
        sim = i, analys = as.integer(rownames(r[[i]]$pr_trt_ni_ref)),r[[i]]$pr_trt_ni_ref) 
    } )))
  
  d_pr_trt_ni_ref_fut <- data.table(
    do.call(rbind, lapply(1:length(r), function(i){ 
      cbind(
        sim = i, analys = as.integer(rownames(r[[i]]$pr_trt_ni_ref_fut)), r[[i]]$pr_trt_ni_ref_fut) 
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
    d_pr_sup_fut = d_pr_sup_fut,
    d_pr_trt_ni_ref = d_pr_trt_ni_ref,
    d_pr_trt_ni_ref_fut = d_pr_trt_ni_ref_fut,
    d_decision = d_decision,
    d_post_smry_1 = d_post_smry_1,
    d_post_smry_2 = d_post_smry_2,
    d_grp = d_grp
    )
  
  toks <- unlist(tstrsplit(args[2], "[-.]"))
  fname <- paste0("data/sim03-", toks[3], "-", toks[4], "-", format(Sys.time(), "%Y%m%d-%H%M%S"), ".qs")
  qs::qsave(l, file = fname)
}

run_none_sim_03 <- function(){
  log_info("run_none_sim_03: Nothing doing here bud.")
}

main_sim_03 <- function(){
  funcname <- paste0(args[1], "()")
  log_info("Main, invoking ", funcname)
  eval(parse(text=funcname))
}

main_sim_03()


