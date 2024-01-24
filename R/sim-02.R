source("./R/init.R")

# Command line arguments list scenario (true dose response),
# the number of simulations to run, the number of cores
# to use and the simulator to use.
args = commandArgs(trailingOnly=TRUE)

if (length(args)<1) {
  log_info("Setting default run method (does nothing)")
  args[1] = "run_none_sim_02"
  args[2] = "cfg-sim02-sc01-v02.yml"
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
  
  # init
  log_trace("Create sim_spec", ix)
  
  pop_spec <- get_pop_spec()
  setkey(pop_spec$r_a, silo, a)
  setkey(pop_spec$r_b, silo, qa, b)
  setkey(pop_spec$r_c, c)
  
  sim_spec <- get_sim_spec()
  sim_spec$b_a_late["rev"] <- g_cfgsc$b_a_l_2
  sim_spec$b_a_chronic["two"] <- g_cfgsc$b_a_c_2
  sim_spec$b_b1_late_one["w12p1"] <- g_cfgsc$b_b1_l_2
  sim_spec$b_b2_late_two["w12p2"] <- g_cfgsc$b_b2_l_2
  sim_spec$b_b1_chronic_one["w12p1"] <- g_cfgsc$b_b1_c_2
  sim_spec$b_b2_chronic_two["w12p2"] <- g_cfgsc$b_b2_c_2
  sim_spec$b_c["rif"] <- g_cfgsc$b_c_2

  
  # enrol times
  # easier to produce enrol times all in one hit rather than try to do them 
  # incrementally with a non-hom pp.
  loc_t0 <- get_enrol_time(max(g_cfgsc$N_pt))

  # loop controls
  stop_enrol <- FALSE
  ii <- 1 # interim number
  N_analys <- length(g_cfgsc$N_pt)

  # tracking/output quantities
  par_out <- array(
    NA, 
    dim = c(N_analys, length(g_effs), 4),
    dimnames = list(1:N_analys, g_effs, c("mu", "med", "q_025", "q_975"))
  )
  pr_sup <- array(
    NA, 
    dim = c(N_analys, length(g_effs)),
    dimnames = list(1:N_analys, g_effs)
  )
  trig <- array(
    NA,
    dim = c(N_analys, length(g_effs), 3),
    dimnames = list(1:N_analys, g_effs, c("sup", "inf", "fut"))
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
    
    l_new <- get_trial_data(N = N_c, pop_spec, sim_spec, is, entry_times = FALSE)

    # combine the existing and new data
    d_all <- rbind(d_all, l_new$d)
    d_all[, t0 := loc_t0[1:.N]]
    
    # create stand data based on enrolled pt
    d_i <- get_indexes(d_all, sim_spec)
    lsd <- get_stan_data(d_i)
    
    
    # fit model
    f1 <- m1$sample(
      lsd$ld, iter_warmup = 1000, iter_sampling = 2000,
      parallel_chains = 2, chains = 2, refresh = 0, show_exceptions = F, 
      max_treedepth = 13)
    
    # extract post
    post_1 <- data.table(f1$draws(variables = c(g_pars), format = "matrix"))
    cols <- names(post_1)
    cols <- gsub("[","_",cols,fixed = T)
    cols <- gsub("]","",cols,fixed = T)
    names(post_1) <- cols
   
    # par out
    m <- post_1[, sapply(.SD, function(z){
      c(mu = mean(z), med = median(z), 
        q_025 = unname(quantile(z, prob = 0.025)), 
        q_975 = unname(quantile(z, prob = 0.975)))}), 
      .SDcols = g_effs]
    par_out[ii, ,] <- t(m)
     
    # probability superiority
    pr_sup[ii, ] <- post_1[, sapply(.SD, function(z){mean(z>0)}), .SDcols = g_effs]
    
    # evaluate triggers
    trig[ii, , "sup"] <- pr_sup[ii, ] > g_cfgsc$d_sup
    trig[ii, , "inf"] <- (1-pr_sup[ii, ]) > g_cfgsc$d_inf
    trig[ii, , "fut"] <- pr_sup[ii, ] < g_cfgsc$d_fut
    
    # Earlier decisions are retained - once superiority has been decided, 
    # we retain this conclusion irrespective of subsequent post probabilities.
    # This means that there could be inconsistency with a silo pr_sup and the 
    # decision reported.
    trig[1:ii, , "sup"] <- apply(trig[1:ii, , "sup", drop = F], 2, function(z){ cumsum(z) > 0 })
    trig[1:ii, , "inf"] <- apply(trig[1:ii, , "inf", drop = F], 2, function(z){ cumsum(z) > 0 })
    trig[1:ii, , "fut"] <- apply(trig[1:ii, , "fut", drop = F], 2, function(z){ cumsum(z) > 0 })
    

    # update pop_spec
    # keys were initialised earlier - mandatory.
    setkey(pop_spec$r_a, silo, a)
    setkey(pop_spec$r_b, silo, qa, b)
    setkey(pop_spec$r_c, c)
    if(any(trig[ii, , "sup"])){
      if(trig[ii, "b_a_l_2", "sup"]){
        pop_spec$r_a[.("late", "dair"), p := 0]
        pop_spec$r_a[.("late", "rev"), p := 1]
      }
      if(trig[ii, "b_a_c_2", "sup"]){
        pop_spec$r_a[.("chronic", "one"), p := 0]
        pop_spec$r_a[.("chronic", "two"), p := 1]
      }
      if(trig[ii, "b_b1_l_2", "sup"]){
        pop_spec$r_b[.("late", "one", "w06p1"), p := 0]
        pop_spec$r_b[.("late", "one", "w12p1"), p := 1]
      }
      if(trig[ii, "b_b2_l_2", "sup"]){
        pop_spec$r_b[.("late", "two", "d07p2"), p := 0]
        pop_spec$r_b[.("late", "two", "w12p2"), p := 1]
      }
      if(trig[ii, "b_b1_c_2", "sup"]){
        pop_spec$r_b[.("chronic", "one", "w06p1"), p := 0]
        pop_spec$r_b[.("chronic", "one", "w12p1"), p := 1]
      }
      if(trig[ii, "b_b2_c_2", "sup"]){
        pop_spec$r_b[.("chronic", "two", "d07p2"), p := 0]
        pop_spec$r_b[.("chronic", "two", "w12p2"), p := 1]
      }
      # not a silo-specific decision - stops in all
      if(trig[ii, "b_c_2", "sup"]){
        pop_spec$r_c[.("norif"), p := 0]
        pop_spec$r_c[.("rif"), p := 1]
      }
    } 
    if(any(trig[ii, , "inf"])){
      if(trig[ii, "b_a_l_2", "inf"]){
        pop_spec$r_a[.("late", "dair"), p := 1]
        pop_spec$r_a[.("late", "rev"), p := 0]
      }
      if(trig[ii, "b_a_c_2", "inf"]){
        pop_spec$r_a[.("chronic", "one"), p := 1]
        pop_spec$r_a[.("chronic", "two"), p := 0]
      }
      if(trig[ii, "b_b1_l_2", "inf"]){
        pop_spec$r_b[.("late", "one", "w06p1"), p := 1]
        pop_spec$r_b[.("late", "one", "w12p1"), p := 0]
      }
      if(trig[ii, "b_b2_l_2", "inf"]){
        pop_spec$r_b[.("late", "two", "d07p2"), p := 1]
        pop_spec$r_b[.("late", "two", "w12p2"), p := 0]
      }
      if(trig[ii, "b_b1_c_2", "inf"]){
        pop_spec$r_b[.("chronic", "one", "w06p1"), p := 1]
        pop_spec$r_b[.("chronic", "one", "w12p1"), p := 0]
      }
      if(trig[ii, "b_b2_c_2", "inf"]){
        pop_spec$r_b[.("chronic", "two", "d07p2"), p := 1]
        pop_spec$r_b[.("chronic", "two", "w12p2"), p := 0]
      }
      # not a silo-specific decision - stops in all
      if(trig[ii, "b_c_2", "inf"]){
        pop_spec$r_c[.("norif"), p := 1]
        pop_spec$r_c[.("rif"), p := 0]
      }
    }
    if(any(trig[ii, , "fut"])){
      if(trig[ii, "b_a_l_2", "fut"]){
        pop_spec$r_a[.("late", "dair"), p := 1]
        pop_spec$r_a[.("late", "rev"), p := 0]
      }
      if(trig[ii, "b_a_c_2", "fut"]){
        pop_spec$r_a[.("chronic", "one"), p := 1]
        pop_spec$r_a[.("chronic", "two"), p := 0]
      }
      if(trig[ii, "b_b1_l_2", "fut"]){
        pop_spec$r_b[.("late", "one", "w06p1"), p := 1]
        pop_spec$r_b[.("late", "one", "w12p1"), p := 0]
      }
      if(trig[ii, "b_b2_l_2", "fut"]){
        pop_spec$r_b[.("late", "two", "d07p2"), p := 1]
        pop_spec$r_b[.("late", "two", "w12p2"), p := 0]
      }
      if(trig[ii, "b_b1_c_2", "fut"]){
        pop_spec$r_b[.("chronic", "one", "w06p1"), p := 1]
        pop_spec$r_b[.("chronic", "one", "w12p1"), p := 0]
      }
      if(trig[ii, "b_b2_c_2", "fut"]){
        pop_spec$r_b[.("chronic", "two", "d07p2"), p := 1]
        pop_spec$r_b[.("chronic", "two", "w12p2"), p := 0]
      }
      # not a silo-specific decision - stops in all
      if(trig[ii, "b_c_2", "fut"]){
        pop_spec$r_c[.("norif"), p := 1]
        pop_spec$r_c[.("rif"), p := 0]
      }
    }
    
    # trig
    # list(pop_spec$r_a, pop_spec$r_b, pop_spec$r_c)
    
    
    # test if all stopped...
    if(all(trig[ii, , "sup"]) | all(trig[ii, , "inf"]) | all(trig[ii, , "fut"])){
      stop_enrol <- T  
    }
    
    # next interim
    ii <- ii + 1
    
    if(ii > N_analys){
      stop_enrol <- T  
    }
  }
  
  if(ii < N_analys){
    log_info("Stopped at analysis ", ii - 1, " filling in missing entries")
    
    trig[ii:N_analys, , "sup"] <- trig[(ii-1), , "sup"]
    trig[ii:N_analys, , "inf"] <- trig[(ii-1), , "inf"]
    trig[ii:N_analys, , "fut"] <- trig[(ii-1), , "fut"]
  }

  
  list(
    d_grp = d_all[, .(y = sum(y), .N), keyby = .(silo, joint, ea, a, qa, eb, b, ec, c)],
    par_out = par_out,
    trig = trig,
    pr_sup = pr_sup  
  )
}

run_sim_02 <- function(){
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
  
  # put results into data.tables
  d_pr_sup <- data.table(
    do.call(rbind, lapply(1:length(r), function(i){ 
      cbind(
        sim = i, analys = as.integer(rownames(r[[i]]$pr_sup)),r[[i]]$pr_sup) 
      } )))
  
  d_trig <- rbind(
    data.table(do.call(rbind, lapply(1:length(r), function(i){ 
      m <- r[[i]]$trig[, , "sup"]
      cbind(sim = i, analys = as.integer(rownames(m)), quant = "sup", m) 
    } ))), 
    
    data.table(do.call(rbind, lapply(1:length(r), function(i){ 
      m <- r[[i]]$trig[, , "inf"]
      cbind(sim = i, analys = as.integer(rownames(m)), quant = "inf", m) 
    } ))),
    
    data.table(do.call(rbind, lapply(1:length(r), function(i){ 
      m <- r[[i]]$trig[, , "fut"]
      cbind(sim = i, analys = as.integer(rownames(m)), quant = "fut", m) 
    } )))
  )
  d_trig[, `:=`(sim = as.integer(sim), analys = as.integer(analys),
                b_a_l_2 = as.logical(b_a_l_2), 
                b_a_c_2 = as.logical(b_a_c_2), 
                b_b1_l_2 = as.logical(b_b1_l_2),  
                b_b2_l_2 = as.logical(b_b2_l_2), 
                b_b1_c_2 = as.logical(b_b1_c_2), 
                b_b2_c_2 = as.logical(b_b2_c_2), 
                b_c_2 = as.logical(b_c_2))]
  d_par <- data.table(do.call(rbind, lapply(1:length(r), function(i){ 
    m <- data.table(
      par = rep(colnames(r[[i]]$par_out), 
                len = dim(r[[i]]$par_out)[1] * dim(r[[i]]$par_out)[2]), 
      smry = rep(dimnames(r[[i]]$par_out)[[3]], 
                 each = dim(r[[i]]$par_out)[2]),
      apply(r[[i]]$par_out,1,c)
    )
    m <- melt(m, id.vars = c("par", "smry"), variable.name = "analys")
    m <- dcast(m, par + analys ~ smry, value.var = "value")
    m <- cbind(sim = i, m)
    m
  } )))
  d_par[, analys := as.integer(as.character(analys))]
  d_grp <- rbindlist(lapply(1:length(r), function(i){ 
    r[[i]]$d_grp
  } ), idcol = "sim")

  l <- list(
    cfg = g_cfgsc,
    d_pr_sup = d_pr_sup, 
    d_trig = d_trig,
    d_par = d_par,
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


