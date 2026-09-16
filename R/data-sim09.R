library(data.table)
library(ggplot2)
library(patchwork)
library(fastglm)
library(parallel)
library(pbapply)
library(kableExtra)
library(here)
library(logger)
library(cmdstanr)
library(qs2)


f_log <-  here::here("logs", "log.txt")
logger::log_appender(appender_file(f_log))
# message(Sys.time(), " Log file initialised ", f_log)
logger::log_info("*** START UP - Sim09 ***")


# Command line arguments list scenario (true dose response),
# the number of simulations to run, the number of cores
# to use and the simulator to use.
args = commandArgs(trailingOnly=TRUE)

# Load cfg based on cmd line args.
if (length(args)<1) {
  log_info("Setting default run method (does nothing)")
  args[1] = "sim09_run_none"
  args[2] = "sim09/cfg-sim09-sc01-v02.yml"
} else {
  log_info("Run method ", args[1])
  log_info("Scenario config ", args[2])
}

# compile models at bottom of scipt
# if(!interactive()){
#   m_1 <- cmdstanr::cmdstan_model(cmdstanr::write_stan_file(s_mod))
# } else {
#   # based on locally stored file - for future use
#   m_1 <- cmdstanr::cmdstan_model(here::here("stan", "model-sim-09.stan"))
# }
# m_1 <- cmdstanr::cmdstan_model(here::here("stan", "model-sim-09-a.stan"))
# m_2 <- cmdstanr::cmdstan_model(here::here("stan", "model-sim-09-c.stan"))
# # hierarchical
# m_3 <- cmdstanr::cmdstan_model(here::here("stan", "model-sim-09-b.stan"))


sim09_run_trial <- function(
    l_spec,
    # initial domain state (allocation wgts)
    l_dom_state = sim09_domain_state_open(),
    # function that runs analyses and returns revised domain allocations
    # once decision threshold realised
    fn_decision = sim09_decision_fn_01,
    # data processing function - convert the row level data into binomial data
    # then put into a list format suitable for stan model
    fn_data = sim09_stan_data_01,
    # fit the model and transform the posterior into the weighted contributions
    # we need
    fn_stanfit = sim09_stan_fit_01
){
  
  log_info("sim09_run_trial: starting trial ", l_spec$ix_sim)
  # accrued data
  d_cum_dat   <- data.table()
  
  # probably not necessary, we can put everything in single result obj
  # l_state_log <- vector("list", length(l_spec$n_batch))
  
  l_res <- vector("list", length(l_spec$n_batch))
  l_dec <- NULL
  
  i <- 1
  for (i in seq_along(l_spec$n_batch)) {
    
    if(i == 1){
      # starting pt index in data
      l_spec$is <- 1
      l_spec$ie <- l_spec$is + l_spec$n_batch[i] - 1
    } else {
      l_spec$is <- nrow(d_cum_dat) + 1
      l_spec$ie <- l_spec$is + l_spec$n_batch[i] - 1
    }
    
    d_batch_dat <- sim09_batch_01(l_spec, l_dom_state = l_dom_state)
    d_batch_dat[, batch := i]
    d_cum_dat <- rbind(d_cum_dat, d_batch_dat)
    
    # everything contained in l_res now
    # state that generated batch i (starts from original/opening state setup)
    # l_state_log[[i]] <- l_dom_state         
    
    # updated for batch i+1
    l_res[[i]]   <- fn_decision(
      d_cum_dat, 
      l_dom_state, 
      i,
      l_spec, 
      l_dec,
      fn_data,
      fn_stanfit
      )  
    
    l_dom_state <- copy(l_res[[i]]$l_dom_state_new)
    # update history of decisions so we don't flip flop
    l_dec <- copy(l_res[[i]]$l_dec_new)
    
    # if all domains have been resolved in terms of decided either superiority, 
    # non inferiority or futility, then we would exit loop and return results
    # to this point under the assumption that all domains would be closed to 
    # future enrolment.
    
    # one issue is the potential for decisions made in one analysis to reverse
    # in a subsequent analysis. 
    # we don't want to flip flop turning domains on and off so we want to 
    # prevent this reversal by overriding future decisions once a domain 
    # decision has been made. for example, if we decide that d1 shows superiority
    # in the first interim then we retain that decision for all future enrolments
    # irrespective of whether the analysis suggests that we should change our
    # mind or not. 
    
    # break once done
    d_dec <- sim09_extract_dec_indicators(l_res[[i]]$l_dec_new)
    d_resolved <- d_dec[, .(resolved = any(dec)), by = domain]
    if (sum(d_resolved$resolved) == nrow(d_resolved)) {
      log_info("sim09_run_trial: All domains resolved")
      break
    }
    
  }
  
  
  log_info("sim09_run_trial: finished trial ", l_spec$ix_sim)
  
  list(
    data = d_cum_dat,
    # includes the domain state entering and after each analysis
    l_res = l_res
    )
  
  
}




# analysis and decision processing ---------
sim09_decision_fn_01 <- function(
    d_cum_dat, 
    l_dom_state, 
    batch,
    l_spec,
    l_dec, 
    fn_data = sim09_stan_data_01,
    fn_stanfit = sim09_stan_fit_01
){
  
  
  l_fit <- fn_stanfit(
    d_cum_dat, fn_data, l_spec
  )
  
  l_dec_new <- sim09_eval_all_dec(l_fit$d_rd, l_spec)
  # prevent decisions from flip-flop 
  l_dec_new <- sim09_lock_dec(l_dec, l_dec_new)          
  
  # update new state (accounts for locks)
  l_dom_state_new <- sim09_update_dom_state(l_dom_state, l_dec_new)  
  
  
  
  l_res <- list(
    # retain previous state so that we can just return a set of res objs with
    # the original state included
    l_dom_state_old = l_dom_state,
    # and the new state
    l_dom_state_new = l_dom_state_new,
    
    # posterior summary
    l_smry = list(
      d_lor_smry = l_fit$d_lor_smry,
      d_lor_std_smry = l_fit$d_lor_std_smry,
      d_rd_smry = l_fit$d_rd_smry
    ),
    # decision based on rules for each domain
    l_dec_new = l_dec_new,
    
    retn_post = l_spec$return_posterior
  )
  
  if(l_spec$return_posterior){
    l_res[["d_lor"]] <- l_fit$d_lor
    l_res[["d_lor_std"]] <- l_fit$d_lor_std
    l_res[["d_rd"]] <- l_fit$d_rd
  }
  
  # return results which may include full posterior if configured to do so
  l_res 
}



sim09_decision_fn_dummy <- function(
    d_cum_dat, 
    l_dom_state, 
    batch
){
  
  # just a dummy placeholder update on the third interim so that batch 4 and onwards
  # don't randomised d2 
  
  # in practice, this would possibly invoke the analysis from here and make the 
  # decision on the basis of the results.
  
  l_dom_state 
}

sim09_d1_alloc <- function(l_dec) {
  # safe to use isTRUE since sup is a single value
  if (isTRUE(l_dec$d1$sup)) return(c(dair = 0.0, r1 = 1/3, r2 = 2/3))   # revision superior
  if (isTRUE(l_dec$d1$fut)) return(c(dair = 1.0, r1 = 0.0, r2 = 0.0))   # revision futile
  sim09_domain_state_open()$d1
}
sim09_d2_alloc <- function(l_dec) {
  if (isTRUE(l_dec$d2$ni))  return(c(nad2 = 0.3, wk12 = 0.0, wk6 = 0.7))  # 6wk non-inferior
  if (isTRUE(l_dec$d2$fut)) return(c(nad2 = 0.3, wk12 = 0.7, wk6 = 0.0))  # 6wk futile
  sim09_domain_state_open()$d2
}
sim09_d3_alloc <- function(l_dec) {
  if (isTRUE(l_dec$d3$sup)) return(c(nad3 = 0.3, wk12 = 0.7, none = 0.0)) # wk12 superior
  if (isTRUE(l_dec$d3$fut)) return(c(nad3 = 0.3, wk12 = 0.0, none = 0.7)) # wk12 futile
  sim09_domain_state_open()$d3
}
sim09_d4_alloc <- function(l_dec) {
  if (isTRUE(l_dec$d4$sup)) return(c(nad4 = 0.3, norif = 0.0, rif = 0.7)) # rif superior
  if (isTRUE(l_dec$d4$fut)) return(c(nad4 = 0.3, norif = 0.7, rif = 0.0)) # rif futile
  sim09_domain_state_open()$d4
}

sim09_update_dom_state <- function(l_dom_state, l_dec) {
  l_dom_state$d1 <- sim09_d1_alloc(l_dec)
  l_dom_state$d2 <- sim09_d2_alloc(l_dec)
  l_dom_state$d3 <- sim09_d3_alloc(l_dec)
  l_dom_state$d4 <- sim09_d4_alloc(l_dec)
  l_dom_state
}


# Once a rule (sup/ni/fut) for a domain flips TRUE, it stays TRUE for every
# subsequent interim regardless of what a later analysis concludes. 
# NA  never overrides a locked TRUE, and never itself counts as decided.
sim09_lock_dec <- function(l_dec_prev, l_dec_new) {
  if (is.null(l_dec_prev)) return(l_dec_new)
  
  l_locked <- l_dec_new
  for (dm in names(l_dec_new)) {
    rules <- names(l_dec_new[[dm]])
    rules <- rules[!grepl("_prob$", rules)]
    for (rl in rules) {
      l_locked[[dm]][[rl]] <- isTRUE(l_dec_prev[[dm]][[rl]]) || isTRUE(l_dec_new[[dm]][[rl]])
    }
  }
  l_locked
}

# More or less generic superiority/non-inferiority + futility decision rule, 
# applied posterior samples for domain contrast theta (on whatever scale
# those draws are - log-odds or risk difference, doesn't matter to this
# function). 
# "sup" and "ni" are mechanically identical - Pr(theta > delta) >
# thresh => success - the distinction is purely in how delta/thresh are
# chosen (e.g. delta = 0 for superiority, delta = -0.05 for a non-inferiority
# margin on a risk-difference scale). 
# Futility is always the mirror check:
# Pr(theta > delta_fut) < thresh_fut => futile.
#
# rule is a list that may contain $sup and/or $ni, and/or $fut, each of the
# form list(delta = ..., thresh = ...). Domains only need to supply whichever
# of these apply to them (e.g. d2 supplies ni + fut; d1/d3/d4 supply sup + fut).
sim09_eval_rule <- function(post_draws, rule) {
  
  out <- list()
  
  if (!is.null(rule$sup)) {
    p <- mean(post_draws > rule$sup$delta, na.rm = TRUE)
    out$sup_prob <- p
    out$sup <- if (all(is.na(post_draws))) NA else p > rule$sup$thresh
  }
  
  if (!is.null(rule$ni)) {
    p <- mean(post_draws > rule$ni$delta, na.rm = TRUE)
    out$ni_prob <- p
    out$ni <- if (all(is.na(post_draws))) NA else p > rule$ni$thresh
  }
  
  if (!is.null(rule$fut)) {
    p <- mean(post_draws > rule$fut$delta, na.rm = TRUE)
    out$fut_prob <- p
    out$fut <- if (all(is.na(post_draws))) NA else p < rule$fut$thresh
  }
  
  out
}


# apply sim09_eval_rule across all four domains at once, given a data.table
# of risk-difference (or log-odds) posterior draws with columns d1, d2, d3, d4
sim09_eval_all_dec <- function(d_post, l_spec) {
  list(
    d1 = sim09_eval_rule(d_post$d1, l_spec$dec$d1),
    d2 = sim09_eval_rule(d_post$d2, l_spec$dec$d2),
    d3 = sim09_eval_rule(d_post$d3, l_spec$dec$d3),
    d4 = sim09_eval_rule(d_post$d4, l_spec$dec$d4)
  )
}






# Computing treatment contrasts ----------
sim09_std_prob <- function(
    v_b0, m_reg, m_d4, cov_grid
    ) {
  stopifnot(all(cov_grid$reg %in% colnames(m_reg)))
  stopifnot(all(cov_grid$d4  %in% colnames(m_d4)))
  n_draws <- length(v_b0)
  if (any(is.na(cov_grid$w))) return(rep(NA_real_, n_draws))
  p_acc <- numeric(n_draws)
  ii <- 1
  for (ii in seq_len(nrow(cov_grid))) {
    eta <- v_b0 + m_reg[, cov_grid$reg[ii]] + m_d4[, cov_grid$d4[ii]]
    p_acc <- p_acc + cov_grid$w[ii] * plogis(eta)
  }
  p_acc
}

# observed-proportion weights for a set of regimens; NA (not 0) if none of
# them have been observed yet, so a domain contrast with no supporting data
# comes back as NA rather than a silently as zero
sim09_reg_wgt <- function(
    regs, d_cum_dat
    ) {
  tb <- d_cum_dat[reg %in% regs, .N, keyby = reg]
  # explicitly recognise zero contribs as no informaiton
  if (sum(tb$N) == 0) return(setNames(rep(NA_real_, length(regs)), regs))
  w <- setNames(rep(0, length(regs)), regs)
  w[as.character(tb$reg)] <- tb$N / sum(tb$N)
  w[regs]
}

# Tentatively refer to as weighted conditional/model-scale log OR
sim09_compute_lor <- function(
    d_cum_dat, l_spec, f_1
    ) {
  
  # domain-level weighted contrasts, from posterior draws
  # Same weighting logic as sim09_ex_sim_2 (observed proportions across each
  # domain's contributing regimens), but applied to the full posterior of
  # b_reg/b_d4 rather than to a single point estimate
  
  # column j of b_reg <-> reg_lvls[j]
  reg_lvls <- levels(d_cum_dat$reg)   
  # column j of b_d4  <-> d4_lvls[j]
  d4_lvls  <- levels(d_cum_dat$d4)    
  
  m_reg <- as.matrix(f_1$draws(variables = "b_reg", format = "matrix"))
  colnames(m_reg) <- reg_lvls
  m_d4  <- as.matrix(f_1$draws(variables = "b_d4",  format = "matrix"))
  colnames(m_d4) <- d4_lvls
  
  w_d1 <- sim09_reg_wgt(l_spec$d1_trt_regs, d_cum_dat)
  post_d1 <- as.numeric(m_reg[, l_spec$d1_trt_regs, drop = FALSE] %*% w_d1)
  
  w_d2_12 <- sim09_reg_wgt(l_spec$d2_wk12_regs, d_cum_dat)
  w_d2_6  <- sim09_reg_wgt(l_spec$d2_wk6_regs, d_cum_dat)
  post_d2 <- as.numeric(
    m_reg[, l_spec$d2_wk6_regs,  drop = FALSE] %*% w_d2_6 -
      m_reg[, l_spec$d2_wk12_regs, drop = FALSE] %*% w_d2_12
  )
  
  w_d3_none <- sim09_reg_wgt(l_spec$d3_none_regs, d_cum_dat)
  w_d3_12   <- sim09_reg_wgt(l_spec$d3_wk12_regs, d_cum_dat)
  post_d3 <- as.numeric(
    m_reg[, l_spec$d3_wk12_regs, drop = FALSE] %*% w_d3_12 -
      m_reg[, l_spec$d3_none_regs, drop = FALSE] %*% w_d3_none
  )
  
  # d4 is unweighted - it's already a direct contrast between two b_d4 levels
  post_d4 <- as.numeric(m_d4[, "rif"] - m_d4[, "norif"])
  
  d_lor <- data.table(
    d1 = post_d1, d2 = post_d2, d3 = post_d3, d4 = post_d4
  )
  
  # d_fig <- melt(d_lor, measure.vars = names(d_lor))
  # ggplot(d_fig, aes(x = value)) +
  #   geom_density() + facet_wrap(~variable)
  
  d_lor
}



# Tentatively refer to as standardised marginal log or
sim09_compute_lor_std <- function(
    d_cum_dat, l_spec, f_1
) {
  
  v_b0 <- as.numeric(f_1$draws(variables = "b_0", format = "matrix"))
  
  reg_lvls <- levels(d_cum_dat$reg)
  d4_lvls  <- levels(d_cum_dat$d4)
  
  m_reg <- as.matrix(f_1$draws(variables = "b_reg", format = "matrix"))
  colnames(m_reg) <- reg_lvls
  
  m_d4 <- as.matrix(f_1$draws(variables = "b_d4", format = "matrix"))
  colnames(m_d4) <- d4_lvls
  
  
  # d1: revision vs DAIR
  # population is late-silo patients, standardised to their observed d4 distribution.
  # In contrast to the above compute_lor, revision represented by the observed 
  # distribution oover complete regimens and d4
  w_d1_trt <- sim09_reg_wgt(l_spec$d1_trt_regs, d_cum_dat )
  # mirrors the marginal rd approach
  d_pop_l <- d_cum_dat[ silo == "l", .N, keyby = d4 ]
  d_pop_l[, w := N / sum(N)]
  
  # tbd should probably pull some of this code out into a function that can be 
  # called from compute_rd and this function....
  # DAIR
  grid_d1_dair <- data.table(
    reg = "l_dair_nad2_nad3", d4 = as.character(d_pop_l$d4), w = d_pop_l$w
  )
  # Revision
  grid_d1_rev <- CJ(
    reg = l_spec$d1_trt_regs, d4 = as.character(d_pop_l$d4)
  )
  
  grid_d1_rev <- base::merge(
    grid_d1_rev,
    data.table(
      reg = l_spec$d1_trt_regs,
      w_trt = w_d1_trt
    ),
    by = "reg"
  )
  # now add in d4
  grid_d1_rev <- base::merge(
    grid_d1_rev,
    data.table(
      d4 = as.character(d_pop_l$d4),
      w_d4 = d_pop_l$w
    ),
    by = "d4"
  )
  # to get the final weights
  grid_d1_rev[, w := w_trt * w_d4]
  
  p_d1_dair <- sim09_std_prob(v_b0, m_reg, m_d4, grid_d1_dair)
  p_d1_rev <- sim09_std_prob(v_b0, m_reg, m_d4, grid_d1_rev)
  # marginal or
  post_lor_d1 <- qlogis(p_d1_rev) - qlogis(p_d1_dair)
  
  
  # d2: 6 weeks vs 12 weeks
  # population: observed r1 patients, standardised jointly over silo and d4.
  # The SAME target weights are used for wk6 and wk12.
  d_pop_r1 <- d_cum_dat[ d1 == "r1", .N, keyby = .(silo, d4) ]
  d_pop_r1[, w := N / sum(N)]
  d_pop_r1[, reg_wk12 := paste0(silo, "_r1_wk12_nad3") ]
  d_pop_r1[, reg_wk6 := paste0(silo, "_r1_wk6_nad3") ]
  
  grid_d2_wk12 <- data.table(
    reg = d_pop_r1$reg_wk12, d4 = as.character(d_pop_r1$d4), w = d_pop_r1$w
  )
  grid_d2_wk6 <- data.table(
    reg = d_pop_r1$reg_wk6, d4 = as.character(d_pop_r1$d4), w = d_pop_r1$w
  )
  
  p_d2_wk12 <- sim09_std_prob( v_b0, m_reg, m_d4, grid_d2_wk12 )
  p_d2_wk6 <- sim09_std_prob( v_b0, m_reg, m_d4, grid_d2_wk6 )
  post_lor_d2 <- qlogis(p_d2_wk6) - qlogis(p_d2_wk12)
  
  # d3: 12 weeks vs none
  # population: observed r2 patients, standardised jointly over silo and d4.
  
  d_pop_r2 <- d_cum_dat[ d1 == "r2", .N, keyby = .(silo, d4) ]
  d_pop_r2[, w := N / sum(N)]
  d_pop_r2[, reg_wk12 := paste0(silo, "_r2_nad2_wk12") ]
  d_pop_r2[, reg_none := paste0(silo, "_r2_nad2_none") ]
  
  grid_d3_wk12 <- data.table(
    reg = d_pop_r2$reg_wk12, d4 = as.character(d_pop_r2$d4), w = d_pop_r2$w
  )
  grid_d3_none <- data.table(
    reg = d_pop_r2$reg_none, d4 = as.character(d_pop_r2$d4), w = d_pop_r2$w
  )
  p_d3_wk12 <- sim09_std_prob( v_b0, m_reg, m_d4, grid_d3_wk12 )
  p_d3_none <- sim09_std_prob( v_b0, m_reg, m_d4, grid_d3_none )
  post_lor_d3 <- qlogis(p_d3_wk12) - qlogis(p_d3_none)
  
  
  # d4: rif vs norif
  # population: all observed patients, standardised to their observed regimen
  
  d_pop_all <- d_cum_dat[ , .N, keyby = reg ]
  d_pop_all[, w := N / sum(N)]
  
  grid_d4_rif <- data.table(
    reg = as.character(d_pop_all$reg), d4 = "rif", w = d_pop_all$w
  )
  grid_d4_norif <- data.table(
    reg = as.character(d_pop_all$reg), d4 = "norif", w = d_pop_all$w
  )
  p_d4_rif <- sim09_std_prob( v_b0, m_reg, m_d4, grid_d4_rif )
  p_d4_norif <- sim09_std_prob( v_b0, m_reg, m_d4, grid_d4_norif)
  post_lor_d4 <- qlogis(p_d4_rif) - qlogis(p_d4_norif)
  
  data.table(
    d1 = as.numeric(post_lor_d1),
    d2 = as.numeric(post_lor_d2),
    d3 = as.numeric(post_lor_d3),
    d4 = as.numeric(post_lor_d4)
  )
}

# prototype - bootstap version tbc...
sim09_reg_wgt_boot <- function(
    regs, d_cum_dat, n_draws
    ) {
  tb <- d_cum_dat[reg %in% regs, .N, keyby = reg]
  n_full <- setNames(rep(0, length(regs)), regs)
  n_full[as.character(tb$reg)] <- tb$N
  if (sum(n_full) == 0) {
    return(matrix(NA_real_, n_draws, length(regs), dimnames = list(NULL, regs)))
  }
  g <- sapply(
    n_full, function(a) {
      if (a == 0) {
        rep(0, n_draws) 
      } else { 
          rgamma(n_draws, shape = a, rate = 1) 
        }
      })
  # each row ~ Dirichlet(n_full), one row per posterior draw
  g / rowSums(g)   
}

# Producing an RD aligned with a population/regimen mix per what is observed 
# in the sample. 
# In contrast the bayesian bootstrap would try to reflect the extra uncertainty
# about whether the observed mix of regimens is a reliable estimate of the 
# population mix.
#
# d2/d3/d4 are direct. For the ACTUAL patients in the applicable population
# (r1 patients for d2, r2 patients for d3, everyone for d4), predict 
# outcome prob under each of the trt levels, holding  everything else 
# (their own silo, their own d4 for d2/d3; their own reg for
# d4) fixed at its observed value, then average the difference.
#
# d1 tricky - a DAIR patient has no observed r1/r2 (or downstream
# d2/d3) counterfactual, because they never entered that pathway.
# Therefore, treat "revision" as a probability-weighted
# mixture over the l_r1_*/l_r2_* regimens, with weights given by the
# OBSERVED proportions of actual revision patients across those cells (i.e.
# the same w_d1 already used for the log-odds jnt_d1 contrast.
# Mixture is applied uniformly to EVERY l-silo patient (not just the
# ones who actually got DAIR), so "revision" and "dair" are both counter-
# factual quantities defined the same way for the whole standardisation
# population.
sim09_comp_rd <- function(
    d_cum_dat, l_spec, f_1
    ) {
  
  v_b0  <- as.numeric(f_1$draws(variables = "b_0", format = "matrix"))
  
  reg_lvls <- levels(d_cum_dat$reg)
  d4_lvls  <- levels(d_cum_dat$d4)
  
  m_reg <- as.matrix(f_1$draws(variables = "b_reg", format = "matrix"))
  colnames(m_reg) <- reg_lvls
  m_d4  <- as.matrix(f_1$draws(variables = "b_d4",  format = "matrix"))
  colnames(m_d4) <- d4_lvls
  
  # d1: revision (mixture over r1/r2 sub-regimens) vs dair, l silo
  # same mixture weights as jnt_d1
  w_d1 <- sim09_reg_wgt(l_spec$d1_trt_regs, d_cum_dat)   
  
  # in practice i think this would need to be across the whole covariate mix, site, 
  # prognostics etc.
  
  d_pop_l <- d_cum_dat[silo == "l", .N, keyby = d4]
  d_pop_l[, w := N / sum(N)]
  
  grid_dair <- data.table(
    reg = "l_dair_nad2_nad3", d4 = as.character(d_pop_l$d4), w = d_pop_l$w)
  
  # contributions across pop
  grid_rev <- CJ(reg = l_spec$d1_trt_regs, d4 = as.character(d_pop_l$d4))
  # weight associated with each regimen
  grid_rev <- base::merge(
    grid_rev, data.table(reg = l_spec$d1_trt_regs, w_treat = w_d1), by = "reg")
  # weights acros d4
  grid_rev <- base::merge(
    grid_rev, data.table(d4 = as.character(d_pop_l$d4), w_covar = d_pop_l$w), by = "d4")
  # combined weight as product
  grid_rev[, w := w_treat * w_covar]
  
  p_dair <- sim09_std_prob(v_b0, m_reg, m_d4, cov_grid = grid_dair)
  p_rev  <- sim09_std_prob(v_b0, m_reg, m_d4, grid_rev)
  rd_d1  <- as.numeric(p_rev - p_dair)
  
  # d2: wk6 vs wk12, weigths among actual r1 patients (all silo)
  d_pop_r1 <- d_cum_dat[d1 == "r1", .N, keyby = .(silo, d4)]
  d_pop_r1[, w := N / sum(N)]
  d_pop_r1[, reg_wk12 := paste0(silo, "_r1_wk12_nad3")]
  d_pop_r1[, reg_wk6  := paste0(silo, "_r1_wk6_nad3")]
  
  grid_wk12 <- data.table(
    reg = d_pop_r1$reg_wk12, d4 = as.character(d_pop_r1$d4), w = d_pop_r1$w)
  grid_wk6  <- data.table(
    reg = d_pop_r1$reg_wk6,  d4 = as.character(d_pop_r1$d4), w = d_pop_r1$w)
  
  p_wk12 <- sim09_std_prob(v_b0, m_reg, m_d4, grid_wk12)
  p_wk6  <- sim09_std_prob(v_b0, m_reg, m_d4, grid_wk6)
  rd_d2  <- as.numeric(p_wk6 - p_wk12)
  
  # d3: wk12 vs none, as above, among actual r2 patients (any silo)
  d_pop_r2 <- d_cum_dat[d1 == "r2", .N, keyby = .(silo, d4)]
  d_pop_r2[, w := N / sum(N)]
  d_pop_r2[, reg_wk12 := paste0(silo, "_r2_nad2_wk12")]
  d_pop_r2[, reg_none := paste0(silo, "_r2_nad2_none")]
  
  grid_d3_wk12 <- data.table(
    reg = d_pop_r2$reg_wk12, d4 = as.character(d_pop_r2$d4), w = d_pop_r2$w)
  grid_d3_none <- data.table(
    reg = d_pop_r2$reg_none, d4 = as.character(d_pop_r2$d4), w = d_pop_r2$w)
  
  p_d3_wk12 <- sim09_std_prob(v_b0, m_reg, m_d4, grid_d3_wk12)
  p_d3_none <- sim09_std_prob(v_b0, m_reg, m_d4, grid_d3_none)
  rd_d3     <- as.numeric(p_d3_wk12 - p_d3_none)
  
  # d4: rif vs norif, among everyone (any reg)
  d_pop_all <- d_cum_dat[, .N, keyby = reg]
  d_pop_all[, w := N / sum(N)]
  
  grid_rif   <- data.table(
    reg = as.character(d_pop_all$reg), d4 = "rif",   w = d_pop_all$w)
  grid_norif <- data.table(
    reg = as.character(d_pop_all$reg), d4 = "norif", w = d_pop_all$w)
  
  p_rif   <- sim09_std_prob(v_b0, m_reg, m_d4, grid_rif)
  p_norif <- sim09_std_prob(v_b0, m_reg, m_d4, grid_norif)
  rd_d4   <- as.numeric(p_rif - p_norif)
  
  d_rd <- data.table(
    d1 = rd_d1, 
    d2 = rd_d2, 
    d3 = rd_d3, 
    d4 = rd_d4
  )
  d_rd
  
}

# Model fit ---------------
sim09_stan_fit_01 <- function(
    d_cum_dat,
    fn_data = sim09_stan_data_01, 
    l_spec
    ){
  
  ld <- fn_data(d_cum_dat, l_spec)
  
  foutname <- paste0(
    format(Sys.time(), format = "%Y%m%d%H%M%S"), 
    "-sim-", l_spec$ix_sim, 
    "-intrm-", max(d_cum_dat$batch))
  
  # snk <- capture.output(
  if(l_spec$mc_model == "indep1"){
    f_1 <- m_1$sample(
      ld, iter_warmup = l_spec$mc_warmup, iter_sampling = l_spec$mc_samp,
      parallel_chains = l_spec$mc_chain, chains = l_spec$mc_chain,
      refresh = 0, show_exceptions = F,
      max_treedepth = 11,
      output_dir = l_spec$mc_out_dir,
      output_basename = foutname
    )
    
    # f_1$summary(variables = c("b_reg"))
    # f_0 <- glm(y ~ reg + d4, data = d_cum_dat, family = binomial)
    # coef(f_0)
    
  } else if (l_spec$mc_model == "indep2"){
    
    f_1 <- m_2$sample(
      ld, iter_warmup = l_spec$mc_warmup, iter_sampling = l_spec$mc_samp,
      parallel_chains = l_spec$mc_chain, chains = l_spec$mc_chain,
      refresh = 0, show_exceptions = F,
      max_treedepth = 11,
      output_dir = l_spec$mc_out_dir,
      output_basename = foutname
    )
    # f_1$summary(variables = c("b_reg"))
  } else if (l_spec$mc_model == "hier"){
    f_1 <- m_3$sample(
      ld, iter_warmup = l_spec$mc_warmup, iter_sampling = l_spec$mc_samp,
      parallel_chains = l_spec$mc_chain, chains = l_spec$mc_chain,
      refresh = 0, show_exceptions = F,
      max_treedepth = 11,
      output_dir = l_spec$mc_out_dir,
      output_basename = foutname
    )
    # f_1$summary(variables = c("mu_reg", "sig_reg"))
    # f_1$summary(variables = c("b_reg"))
  } 
  
  
  # )
  
  # g-comp (standardisation) note ----------
  # In the following g-computation is used to standardise the model based 
  # parameters over a pre-specified target distribution of regimen 
  # characteristics. It is important to note that we are treating the target
  # distribution as fixed and so the posterior uncertainty is reflecting the 
  # uncertainty in the model parameters but not the uncertainty in the 
  # estimation of the target distribution as would be offered via a 
  # bayesian bootstrap. The approach adopted aligns with the common reporting
  # perspective for clinical trials.
  d_lor <- sim09_compute_lor(d_cum_dat, l_spec, f_1)
  d_lor_std <- sim09_compute_lor_std(d_cum_dat, l_spec, f_1)
  d_rd <- sim09_comp_rd(d_cum_dat, l_spec, f_1)
  
  
  par_smry <- function(dat){
    data.table(
      par = names(dat),
      mu = apply(dat, 2, mean),
      q_025 = apply(dat, 2, function(z){quantile(z, prob = 0.025)}),
      q_975 = apply(dat, 2, function(z){quantile(z, prob = 0.975)})
    )
  }
  
  d_lor_smry <- par_smry(d_lor)
  d_lor_std_smry <- par_smry(d_lor_std)
  d_rd_smry <- par_smry(d_rd)
  
  # d_fig <- melt(d_rd, measure.vars = names(d_rd))
  # ggplot(d_fig, aes(x = value)) + geom_density() + facet_wrap(~variable)
  
  list(
    f_1 = f_1,
    d_lor_smry = d_lor_smry,
    # this is the one we should use I thikn....
    d_lor_std_smry = d_lor_std_smry,  
    d_rd_smry = d_rd_smry,
    # one row per posterior draw, one column per domain contrast
    d_lor = d_lor,
    d_lor_std = d_lor_std, 
    d_rd = d_rd
  )
  
}







# Data generation ------------
# Same allocation/outcome structure as sim09_cohort_01, but d2/d3/d4 are now
# drawn directly from domain_state (a single 3-way sample() per domain,
# rather than a two-step enter/split). d1 is unchanged for now - see note
# in the accompanying discussion for how to extend it the same way.
sim09_batch_01 <- function(
    l_spec,
    l_dom_state = sim09_domain_state_open(),   
    seed = NULL
) {
  if (!is.null(seed)) set.seed(seed)
  
  stopifnot(abs(sum(l_spec$p_silo) - 1) < 1e-8)
  stopifnot(abs(sum(l_spec$p_surg_lnrd1) - 1) < 1e-8)
  stopifnot(abs(sum(l_spec$p_surg_enrd1) - 1) < 1e-8)
  stopifnot(abs(sum(l_spec$p_surg_cnrd1) - 1) < 1e-8)
  stopifnot(abs(sum(l_dom_state$d1) - 1) < 1e-8)
  stopifnot(abs(sum(l_dom_state$d2) - 1) < 1e-8)
  stopifnot(abs(sum(l_dom_state$d3) - 1) < 1e-8)
  stopifnot(abs(sum(l_dom_state$d4) - 1) < 1e-8)
  
  # all possible regime option levels
  full_reg_effect <- setNames(rep(0, length(l_spec$reg_opts)), l_spec$reg_opts)
  if (length(l_spec$reg_effect) > 0) {
    unknown <- setdiff(names(l_spec$reg_effect), l_spec$reg_opts)
    if (length(unknown) > 0) {
      stop("reg_effect has names not matching a realised silo_d1_d2_d3 combination: ",
           paste(unknown, collapse = ", "))
    }
    full_reg_effect[names(l_spec$reg_effect)] <- l_spec$reg_effect
  }
  # always fix the reference level at zero
  full_reg_effect["l_dair_nad2_nad3"] <- 0
  
  n <- length(l_spec$is:l_spec$ie)
  d <- data.table(
    id   = l_spec$is:l_spec$ie,
    t_0 = l_spec$t_0[l_spec$is:l_spec$ie],
    silo = sample(names(l_spec$p_silo), n, replace = TRUE, prob = l_spec$p_silo)
  )
  
  # d1: surgery received
  # The l silo draws from l_dom_state$d1 (the randomised comparison, which can be 
  # stopped like the other domains). The other three silos are non-randomised 
  # clinician choice and are fixed distributions.
  d[, d1 := character(.N)]
  idx_l <- d$silo == "l"
  d[idx_l, d1 := sample(names(l_dom_state$d1), .N, replace = TRUE, prob = l_dom_state$d1)]
  
  for (s in c("lnrd1", "enrd1", "cnrd1")) {
    # pick up the sampling dist for this silo
    p_vec <- switch(s, 
                    lnrd1 = l_spec$p_surg_lnrd1, 
                    enrd1 = l_spec$p_surg_enrd1, 
                    cnrd1 = l_spec$p_surg_cnrd1
    )
    # rows for this silo
    idx <- d$silo == s
    # 
    d[idx, d1 := sample(names(p_vec), .N, replace = TRUE, prob = p_vec)]
  }
  d[, d1 := factor(d1, levels = c("dair", "r1", "r2"))]
  
  
  # The d2, d3, d4 are coordinated by domain status
  d[, d2 := "nad2"]
  idx_r1 <- d$d1 == "r1"
  d[idx_r1, d2 := sample(names(l_dom_state$d2), .N, replace = TRUE, prob = l_dom_state$d2)]
  d[, d2 := factor(d2, levels = c("nad2", "wk12", "wk6"))]
  
  d[, d3 := "nad3"]
  idx_r2 <- d$d1 == "r2"
  d[idx_r2, d3 := sample(names(l_dom_state$d3), .N, replace = TRUE, prob = l_dom_state$d3)]
  d[, d3 := factor(d3, levels = c("nad3", "wk12", "none"))]
  
  d[, d4 := sample(names(l_dom_state$d4), .N, replace = TRUE, prob = l_dom_state$d4)]
  d[, d4 := factor(d4, levels = c("nad4", "norif", "rif"))]
  
  d[, silo := factor(silo, levels = c("l", "lnrd1", "enrd1", "cnrd1"))]
  
  # outcome full silo:d1:d2:d3 interaction + additive d4
  reg_key <- paste(d$silo, d$d1, d$d2, d$d3, sep = "_")
  d[, reg := factor(reg_key, levels = l_spec$reg_opts)]
  d[, d1b := copy(d1)]
  d[d1 != "dair", d1b := "rev"]
  d[, d1b := droplevels(d1b)]
  
  intercept <- qlogis(l_spec$response_p_ref)
  lp <- intercept + full_reg_effect[reg_key] + l_spec$d4_effect[as.character(d$d4)]
  
  d[, eta := lp]
  d[, p := plogis(eta)]
  d[, y := rbinom(.N, 1, p)]
  
  d[]
}

sim09_stan_data_01 <- function(d_cum_dat, l_spec){
  
  
  d_grp_dat <- d_cum_dat[, .(n = .N, y = sum(y)), keyby = .(reg, d4, d1, d2, d3, silo)]
  d_grp_dat[, ix_reg := as.integer(reg)]
  d_grp_dat[, ix_d4 := as.integer(d4)]
  
  d_grp_dat[, ix_d1 := as.integer(d1)]
  d_grp_dat[, ix_d2 := as.integer(d2)]
  d_grp_dat[, ix_d3 := as.integer(d3)]
  d_grp_dat[, ix_silo := as.integer(silo)]
  
  if(nrow(d_grp_dat[, .N, keyby = reg]) != length(l_spec$reg_opts)){
    
    reg_zero <- l_spec$reg_opts[!(l_spec$reg_opts %in% d_grp_dat$reg)]
    log_info("regs not present in data ", paste0(reg_zero, collaspse = ", "))
  }
  
  d_grp_reg <- data.table(
    reg = l_spec$reg_opts
  )
  d_grp_reg[, c("silo", "d1", "d2", "d3") := tstrsplit(reg, "_")]
  d_grp_reg[, silo := factor(silo, levels = levels(d_grp_dat$silo))]
  d_grp_reg[, d1 := factor(d1, levels = levels(d_grp_dat$d1))]
  d_grp_reg[, d2 := factor(d2, levels = levels(d_grp_dat$d2))]
  d_grp_reg[, d3 := factor(d3, levels = levels(d_grp_dat$d3))]
  
  d_grp_reg[, ix_d1 := as.integer(d1)]
  d_grp_reg[, ix_d2 := as.integer(d2)]
  d_grp_reg[, ix_d3 := as.integer(d3)]
  d_grp_reg[, ix_silo := as.integer(silo)]
  
  ld <- list(
    N = nrow(d_grp_dat),
    n = d_grp_dat$n,
    y = d_grp_dat$y,
    # reference what is observed - some regimens might not appear in our sample...
    K_reg = nrow(d_grp_reg),
    K_d4 = length(levels(d_grp_dat$d4)),
    reg = d_grp_dat$ix_reg,
    d4 = d_grp_dat$ix_d4,
    
    d1 = d_grp_dat$ix_d1,
    d2 = d_grp_dat$ix_d2,
    d3 = d_grp_dat$ix_d3,
    silo = d_grp_dat$ix_silo,
    
    reg_silo_idx = d_grp_reg$ix_silo,
    reg_d1_idx =   d_grp_reg$ix_d1, 
    reg_d2_idx =   d_grp_reg$ix_d2,
    reg_d3_idx =   d_grp_reg$ix_d3,
    
    pri_b_0 = l_spec$pri_b_0,
    pri_b_reg = l_spec$pri_b_reg,
    pri_b_d4 = l_spec$pri_b_d4,
    prior_only = l_spec$prior_only
  )
  
  ld
}



# Utils --------

sim09_report_sim_res <- function(){
  
  library(data.table)
  library(qs2)
  library(kableExtra)
  
  l <- qs2::qs_read("data/sim09/sim09-v02-20260916-091644.qs2")
  
  r = l$r
  l_spec = l$l_spec
  
  l_oc <- list(
    dec_pr = sim09_smry_pr_dec(r, l_spec),
    
    # maybe no longer necessary...
    l_dec_n = sim09_smry_dec_n(r, l_spec),
    
    # partially duplicates sim09_smry_dec_n
    l_dec_info = sim09_smry_dec_info(r, l_spec),
    
    l_effects = sim09_smry_effects(r, l_spec)
    
  )
  
  
  kableExtra::kbl(
    dcast(l_oc$dec_pr, domain + rule ~ i_anlys, value.var = "mu"),
    digits = 3, format = "simple", 
    caption = paste("Scenario ",
      l_spec$desc, "\nPr decision by domain and interim")
  )
  # sanity
  # l_oc$dec_pr[i_anlys == 5, .(pr_dec = sum(mu)), keyby = domain]
  
  kableExtra::kbl(
    l_oc$l_dec_n$d_smry,
    digits = 3, format = "simple", 
    caption = paste("Scenario ",
        l_spec$desc, "\nEnrolment at time of decision")
  )
  
  kableExtra::kbl(
    l_oc$l_dec_info$d_smry,
    digits = c(0, 0, 3, 1, 0, 0, 1, 1, 1, 3), 
    format = "simple", 
    caption = paste("Scenario ",
        l_spec$desc, "\nEnrolment and sample size informing decisions")
  )
  
  kableExtra::kbl(
    l_oc$l_dec_info$d_arms,
    digits = 1, format = "simple", 
    caption = paste("Scenario ",
        l_spec$desc, "\nSample size informing decisions by arm")
  )
  
  kableExtra::kbl(
    l_oc$l_effects$d_lor,
    digits = 3, format = "simple", 
    caption = paste("Scenario ",
        l_spec$desc, "\nWeighted conditional/model-scale log OR")
  )
  
  kableExtra::kbl(
    l_oc$l_effects$d_lor_std,
    digits = 3, format = "simple", 
    caption = paste("Scenario ",
        l_spec$desc, "\nStandardised marginal log OR")
  )
  
  kableExtra::kbl(
    l_oc$l_effects$d_rd,
    digits = 3, format = "simple", 
    caption = paste("Scenario ",
        l_spec$desc, "\nStandardised marginal RD")
  )
  
  
}

# Truth is the design-population standardised estimand.
# The fitted estimand uses realised population weights, so finite-sample
# differences between truth and analysis target are possible/likely.
sim09_true_effects_des_pop <- function(
    l_spec,
    l_dom_state = sim09_domain_state_open()
    ) {
  
  # probability under a particular regimen and d4 level
  p_reg <- function(reg, d4) {
    
    eta <- qlogis(l_spec$response_p_ref) +
      l_spec$reg_effect[reg] +
      l_spec$d4_effect[d4]
    
    plogis(eta)
  }
  
  # d4 is randomised independently of silo, so just use the domain allocation probabilities
  w_d4 <- l_dom_state$d4
  
  
  # d1: revision vs DAIR
  #
  # rev is a mixture of the six revision regimes:
  #
  #   r1_wk12_nad3
  #   r1_wk6_nad3
  #   r1_nad2_nad3
  #   r2_nad2_wk12
  #   r2_nad2_none
  #   r2_nad2_nad3
  #
  # weights are the expected proportions of revision patients in the regimes.
  # Probability of each revision regimen conditional on being in
  # the late silo and receiving revision.
  
  p_r1 <- l_dom_state$d1["r1"]
  p_r2 <- l_dom_state$d1["r2"]
  # d2 allocation conditional on r1
  p_d2 <- l_dom_state$d2
  # d3 allocation conditional on r2
  p_d3 <- l_dom_state$d3
  
  # only randomised treatment levels contribute to the target revision contrast
  # the nad2/nad3 are retained only where necessary
  w_d1 <- c(
    l_r1_wk12_nad3 = unname(p_r1 * p_d2["wk12"]),
    l_r1_wk6_nad3  = unname(p_r1 * p_d2["wk6"]),
    l_r1_nad2_nad3 = unname(p_r1 * p_d2["nad2"]),
    l_r2_nad2_wk12 = unname(p_r2 * p_d3["wk12"]),
    l_r2_nad2_none = unname(p_r2 * p_d3["none"]),
    l_r2_nad2_nad3 = unname(p_r2 * p_d3["nad3"])
  )
  # need to normalise, because the target is conditional on being a revision
  # patient rather than conditional on entering the late silo.
  w_d1 <- w_d1 / sum(w_d1)
  lor_d1 <- sum(w_d1 * l_spec$reg_effect[names(w_d1)])
  
  
  # d2: wk6 vs wk12
  # weight according to silo distribution among r1 patients
  
  silos <- names(l_spec$p_silo)
  p_r1_silo <- sapply(
    silos,
    function(s) {
      if (s == "l") {
        unname(l_spec$p_silo[s] * l_dom_state$d1["r1"])
      } else {
        unname(l_spec$p_silo[s] * l_spec[[paste0("p_surg_", s)]][["r1"]])
      }
    }
  )
  w_silo_r1 <- p_r1_silo / sum(p_r1_silo)
  
  
  # d2 log OR is the weighted contrast in the corresponding regimen effects.
  reg_wk12_d2 <- paste0(silos, "_r1_wk12_nad3")
  reg_wk6_d2  <- paste0(silos, "_r1_wk6_nad3")
  
  lor_d2 <- sum(
    w_silo_r1 * (l_spec$reg_effect[reg_wk6_d2] - l_spec$reg_effect[reg_wk12_d2]))
  
  
  # d3: wk12 vs none
  p_r2_silo <- sapply(
    silos,
    function(s) {
      if (s == "l") {
        unname(l_spec$p_silo[s] * l_dom_state$d1["r2"])
      } else {
        unname(l_spec$p_silo[s] * l_spec[[paste0("p_surg_", s)]][["r2"]])
      }
    }
  )
  w_silo_r2 <- p_r2_silo / sum(p_r2_silo)
  
  reg_wk12_d3 <- paste0(silos, "_r2_nad2_wk12")
  reg_none_d3 <- paste0(silos, "_r2_nad2_none")
  
  lor_d3 <- sum(
    w_silo_r2 *
      (l_spec$reg_effect[reg_wk12_d3] - l_spec$reg_effect[reg_none_d3])
  )
  
  # d4: rif vs norif
  lor_d4 <- l_spec$d4_effect["rif"] - l_spec$d4_effect["norif"]
  
  # Risk differences
  #
  # urgh. need to average probabilities on the probability scale,
  # rather than transform the averaged log OR.
  
  # d1
  # prob trt succss for each of d4 level
  p_dair_d1 <- sapply(
    names(w_d4),
    function(d4) {
      unname(p_reg("l_dair_nad2_nad3", d4))
    }
  )
  
  p_rev_d1 <- sapply(
    names(w_d4),
    function(d4) {
      
      sum(
        w_d1 *
          sapply(
            names(w_d1),
            function(reg) p_reg(reg, d4)
          )
      )
    }
  )
  
  rd_d1 <- sum(
    w_d4 * (p_rev_d1 - p_dair_d1)
  )
  
  
  # d2
  
  # standardise over the r1 population and d4 distribution.
  p_wk6_d2 <- numeric(length(w_silo_r1))
  p_wk12_d2 <- numeric(length(w_silo_r1))
  
  names(p_wk6_d2) <- silos
  names(p_wk12_d2) <- silos
  
  for (s in silos) {
    
    p_wk6_d2[s] <- sum(
      w_d4 *
        sapply(
          names(w_d4),
          function(d4) {
            p_reg(paste0(s, "_r1_wk6_nad3"), d4)
          }
        )
    )
    
    p_wk12_d2[s] <- sum(
      w_d4 *
        sapply(
          names(w_d4),
          function(d4) {
            p_reg(paste0(s, "_r1_wk12_nad3"), d4)
          }
        )
    )
  }
  
  rd_d2 <- sum(w_silo_r1 * (p_wk6_d2 - p_wk12_d2))
  
  
  # d3 - same deal
  
  p_wk12_d3 <- numeric(length(w_silo_r2))
  p_none_d3 <- numeric(length(w_silo_r2))
  
  names(p_wk12_d3) <- silos
  names(p_none_d3) <- silos
  
  for (s in silos) {
    
    p_wk12_d3[s] <- sum(
      w_d4 *
        sapply(
          names(w_d4),
          function(d4) {
            p_reg(paste0(s, "_r2_nad2_wk12"), d4)
          }
        )
    )
    
    p_none_d3[s] <- sum(
      w_d4 *
        sapply(
          names(w_d4),
          function(d4) {
            p_reg(paste0(s, "_r2_nad2_none"), d4)
          }
        )
    )
  }
  
  rd_d3 <- sum(w_silo_r2 * (p_wk12_d3 - p_none_d3))
  
  
  # d4
  
  # Standardise over the overall regimen distribution.
  # d4 is additive and randomised independently so this could
  # be simplified, but retaining the standardisation
  
  w_reg <- l_spec$p_silo
  
  # probability of each complete regimen in the population
  # generated under the opening domain allocation
  #
  # For d4, only the distribution of reg matters.
  
  all_regs <- l_spec$reg_opts
  
  # Expected probability of each regimen
  p_reg_pop <- setNames(numeric(length(all_regs)), all_regs)
  
  for (s in silos) {
    
    p_s <- l_spec$p_silo[s]
    
    if (s == "l") {
      
      p_d1_s <- l_dom_state$d1
      
    } else {
      
      p_d1_s <- switch(
        s,
        lnrd1 = l_spec$p_surg_lnrd1,
        enrd1 = l_spec$p_surg_enrd1,
        cnrd1 = l_spec$p_surg_cnrd1
      )
    }
    
    for (d1 in names(p_d1_s)) {
      
      if (d1 == "dair") {
        p_d2_s <- c(nad2 = 1)
        p_d3_s <- c(nad3 = 1)
      } else if (d1 == "r1") {
        p_d2_s <- l_dom_state$d2
        p_d3_s <- c(nad3 = 1)
      } else if (d1 == "r2") {
        p_d2_s <- c(nad2 = 1)
        p_d3_s <- l_dom_state$d3
      }
      
      for (d2 in names(p_d2_s)) {
        for (d3 in names(p_d3_s)) {
          
          reg <- paste(s, d1, d2, d3, sep = "_")
          
          if (reg %in% all_regs) {
            p_reg_pop[reg] <-
              p_reg_pop[reg] +
              p_s *
              p_d1_s[d1] *
              p_d2_s[d2] *
              p_d3_s[d3]
          }
        }
      }
    }
  }
  
  p_reg_pop <- p_reg_pop / sum(p_reg_pop)
  
  
  p_rif <- sum(
    p_reg_pop *
      sapply(
        all_regs,
        function(reg) p_reg(reg, "rif")
      )
  )
  
  p_norif <- sum(
    p_reg_pop *
      sapply(
        all_regs,
        function(reg) p_reg(reg, "norif")
      )
  )
  
  rd_d4 <- p_rif - p_norif
  
  
  data.table(
    domain = paste0("d", 1:4),
    lor = c(lor_d1, lor_d2, lor_d3, lor_d4),
    rd  = c(rd_d1, rd_d2, rd_d3, rd_d4)
  )
}

sim09_true_regimen_risk <- function(
    l_spec,
    l_dom_state = sim09_domain_state_open()
    ) {
  
  d4 <- names(l_spec$d4_effect)
  reg <- l_spec$reg_opts
  
  d_out <- CJ(
    reg = reg,
    d4 = d4
  )
  
  d_out[, p := plogis(
    qlogis(l_spec$response_p_ref) +
      l_spec$reg_effect[reg] +
      l_spec$d4_effect[d4]
  )]
  
  d_out
}




sim09_extract_dec_indicators <- function(l_dec){
  
  doms <- names(l_dec)
  
  d_out <- rbindlist(lapply(seq_along(l_dec), function(ii){
    
    z <- l_dec[[ii]]
    
    d_tmp <- data.table(
      domain = doms[ii],
      rule = names(z)[!(names(z) %like% "prob")]
    )
    
    d_tmp[, dec := unlist(z[rule])]
    d_tmp
  }))
  
  d_out
  
}

sim09_smry_pr_dec <- function(r, l_spec){
  
  
  d_first <- sim09_first_dec(r)
  d_first[, dec := TRUE]
  
  # incorrect coz will pick up duplicate decisions, i.e. cumulative probs
  # might be > 1 over decisions for a given domain
  # d_sims <- rbindlist(lapply(r, function(rr){
  #   
  #   rbindlist(lapply(rr$l_res, function(z){
  #     sim09_extract_dec_indicators(z$l_dec)
  #   }), idcol = "i_anlys")
  #   
  # }), idcol = "i_sim")
  # d_sims <- d_sims[order(i_sim, domain, rule, i_anlys)]
  
  d_grid <- CJ(
    i_sim = 1:l_spec$n_sim,
    i_anlys = seq_along(l_spec$n_batch),
    domain = paste0("d", 1:4)
  )
  d_grid <- base::merge(
    data.table(domain = c("d1", "d1", "d2", "d2", "d3", "d3", "d4", "d4"),
               rule = c("sup", "fut", "ni", "fut", "sup", "fut", "sup", "fut")),
    d_grid, by = "domain", all = T, allow.cartesian=TRUE)
  
  d_sims <- base::merge(
    d_grid, 
    d_first, by = c("i_sim", "i_anlys", "domain", "rule"), all.x = T)
  
  # protect against NA
  setorder(d_sims, i_sim, domain, rule, i_anlys)
  d_sims[dec == TRUE,  dec01 := 1]
  d_sims[dec == FALSE, dec01 := 0]
  d_sims[, dec01 := nafill(dec01, type = "locf"), by = .(i_sim, domain, rule)]
  # leading NAs, before any info existed
  d_sims[is.na(dec01), dec01 := 0]   
  
  # d_tmp_fut <- d_sims[domain == "d1" & rule == "fut"]
  # d_tmp_sup <- d_sims[domain == "d1" & rule == "sup"]
  # 
  # d_tmp <- merge(
  #   d_tmp_fut[, .(i_sim, i_anlys, domain, rule, dec01_fut = dec01)],
  #   d_tmp_sup[, .(i_sim, i_anlys, domain, rule, dec01_sup = dec01)],
  #   by = c("i_sim", "i_anlys", "domain"), all.x = T
  # )
  # d_tmp[dec01_fut == 1 & dec01_sup == 1]
  
  d_out <- d_sims[, .(
    mu = mean(dec01)
  ), keyby = .(i_anlys, domain, rule)]
  
  
  
  d_out
  
}


sim09_smry_dec_n <- function(r, l_spec) {
  
  d_out <- rbindlist(
    lapply(seq_along(r), function(i_sim) {
      rr <- r[[i_sim]]
      rbindlist(
        lapply(seq_along(rr$l_res), function(i_anlys) {
          z <- rr$l_res[[i_anlys]]
          if(is.null(z)) return(NULL)
          d_dec <- sim09_extract_dec_indicators(z$l_dec_new)
          d_dec <- d_dec[dec == TRUE]
          if (nrow(d_dec) == 0) {return(NULL)}
          data.table(
            domain = unique(d_dec$domain),
            i_anlys = i_anlys,
            n = sum(l_spec$n_batch[seq_len(i_anlys)])
          )
        })
      )
    }), idcol = "i_sim"
  )
  d_out[, n := as.double(n)]
  
  
  # First decision only
  setorder(d_out, i_sim, domain, i_anlys)
  
  d_first <- d_out[, .SD[1], by = .(i_sim, domain)]
  
  
  # Summary
  d_smry <- d_first[
    ,
    .(
      n_dec = .N,
      pr_dec = .N / l_spec$n_sim,
      mu_n = mean(n),
      q_025_n = quantile(n, 0.025),
      q_975_n = quantile(n, 0.975)
    ),
    by = domain
  ]
  
  list(
    d_first = d_first,
    d_smry = d_smry
  )
}


sim09_smry_effects <- function(
    r, l_spec,
    l_dom_state = sim09_domain_state_open()
) {
  
  d_true <- sim09_true_effects_des_pop(
    l_spec = l_spec,
    l_dom_state = l_dom_state
  )
  
  
  d_lor <- rbindlist(
    lapply(r, function(rr) {
      rbindlist(
        lapply(rr$l_res, function(z) {
          d <- copy(z$l_smry$d_lor_smry)
          d
        }), idcol = "i_anlys" )
    }), idcol = "i_sim")
  
  d_lor_std <- rbindlist(
    lapply(r, function(rr) {
      rbindlist(
        lapply(rr$l_res, function(z) {
          d <- copy(z$l_smry$d_lor_std_smry)
          d
        }), idcol = "i_anlys" )
    }), idcol = "i_sim")
  
  d_rd <- rbindlist(
    lapply(r, function(rr) {
      rbindlist(
        lapply(rr$l_res, function(z) {
          d <- copy(z$l_smry$d_rd_smry)
          d
        }), idcol = "i_anlys" )
    }), idcol = "i_sim")
  
  
  # Complete analysis grid and carry estimates forward
  d_grid <- CJ(
    i_sim = seq_len(l_spec$n_sim),
    i_anlys = seq_along(l_spec$n_batch),
    par = paste0("d", 1:4)
  )
  
  
  d_lor <- base::merge(
    d_grid,
    d_lor,
    by = c("i_sim", "i_anlys", "par"),
    all.x = TRUE
  )
  
  d_lor_std <- base::merge(
    d_grid,
    d_lor_std,
    by = c("i_sim", "i_anlys", "par"),
    all.x = TRUE
  )
  
  d_rd <- base::merge(
    d_grid,
    d_rd,
    by = c("i_sim", "i_anlys", "par"),
    all.x = TRUE
  )
  
  setorder(d_lor, i_sim, par, i_anlys)
  setorder(d_lor_std, i_sim, par, i_anlys)
  setorder(d_rd, i_sim, par, i_anlys)
  
  
  d_lor[, c("mu", "q_025", "q_975") :=
          lapply(.SD, nafill, type = "locf"),
        by = .(i_sim, par),
        .SDcols = c("mu", "q_025", "q_975")]
  
  d_lor_std[, c("mu", "q_025", "q_975") :=
          lapply(.SD, nafill, type = "locf"),
        by = .(i_sim, par),
        .SDcols = c("mu", "q_025", "q_975")]
  
  d_rd[, c("mu", "q_025", "q_975") :=
         lapply(.SD, nafill, type = "locf"),
       by = .(i_sim, par),
       .SDcols = c("mu", "q_025", "q_975")]
  
  
  # Add true vals
  
  d_lor <- base::merge(
    d_lor,
    d_true[, .(par = domain, truth = lor)],
    by = "par",
    all.x = TRUE
  )
  
  d_lor_std <- base::merge(
    d_lor_std,
    # note am using the same lor reference for now
    d_true[, .(par = domain, truth = lor)],
    by = "par",
    all.x = TRUE
  )
  
  d_rd <- base::merge(
    d_rd,
    d_true[, .(par = domain, truth = rd)],
    by = "par",
    all.x = TRUE
  )
  
  d_lor_out <- d_lor[, .(
    truth = data.table::first(truth),
    mean_est = mean(mu, na.rm = TRUE),
    bias = mean(mu - truth, na.rm = TRUE),
    rmse = sqrt(mean((mu - truth)^2, na.rm = TRUE)),
    coverage = mean(
      q_025 <= truth & q_975 >= truth,
      na.rm = TRUE
    )
  ), keyby = .(i_anlys, par)]
  setkey(d_lor_out, par, i_anlys)
  
  d_lor_std_out <- d_lor_std[, .(
    truth = data.table::first(truth),
    mean_est = mean(mu, na.rm = TRUE),
    bias = mean(mu - truth, na.rm = TRUE),
    rmse = sqrt(mean((mu - truth)^2, na.rm = TRUE)),
    coverage = mean(
      q_025 <= truth & q_975 >= truth,
      na.rm = TRUE
    )
  ), keyby = .(i_anlys, par)]
  setkey(d_lor_std_out, par, i_anlys)
  
  d_rd_out <- d_rd[, .(
    truth = data.table::first(truth),
    mean_est = mean(mu, na.rm = TRUE),
    bias = mean(mu - truth, na.rm = TRUE),
    rmse = sqrt(mean((mu - truth)^2, na.rm = TRUE)),
    # proportion of times truth within interval
    coverage = mean(
      q_025 <= truth & q_975 >= truth,
      na.rm = TRUE
    )
  ), keyby = .(i_anlys, par)]
  setkey(d_rd_out, par, i_anlys)
  
  list(
    d_lor = d_lor_out,
    d_lor_std = d_lor_std_out,
    d_rd = d_rd_out
  )
}


sim09_smry_dec_info <- function(r, l_spec) {
  
  d_first <- sim09_first_dec(r)
  
  if (nrow(d_first) == 0) {
    return(list(
      # for each decision
      d_inform_tot = data.table(),
      # averages by arm informing each decision
      d_arms = data.table(),
      # average totals informing each decision
      d_smry = data.table()
    ))
  }
  
  # for each first decision (by domain), obtain the number of pts in the cohort
  # that informed the decision
  # so, for example, suppose that in the first sim at the first analysis d1 
  # indicated rev was superior in the d_first decision list. the output wold 
  # then have a row for the number in the dair, r1 and r2 late silo group at 
  # the time of the first analysis
  d_out <- rbindlist(
    lapply(seq_len(nrow(d_first)), function(i) {
      i_sim   <- d_first$i_sim[i]
      domain  <- d_first$domain[i]
      i_anlys <- d_first$i_anlys[i]
      rule    <- d_first$rule[i]
      rr <- r[[i_sim]]
      d_cum <- rr$data[batch <= i_anlys]
      z <- sim09_domain_n(d_cum, domain)
      z[, `:=`(
        i_sim = i_sim,
        i_anlys = i_anlys,
        rule = rule,
        n_enrolled = nrow(d_cum)
      )]
      z
    })
  )
  
  # Number informing each contrast - ie totals
  d_inform <- d_out[, .(
    n_enrolled = as.double(first(n_enrolled)),
    n_inform = as.double(sum(n[inform]))
    ), by = .(i_sim, domain, i_anlys, rule)
  ]
  d_inform[, inform_prop := n_inform / n_enrolled]
  
  d_summary <- d_inform[
    ,
    .(
      n_dec = .N,
      pr_dec = .N / l_spec$n_sim,
      
      mu_n_enrl = mean(n_enrolled),
      q025_n_enrl = quantile(n_enrolled, 0.025),
      q975_n_enrl = quantile(n_enrolled, 0.975),
      
      mu_n_info = mean(n_inform),
      q025_n_info = quantile(n_inform, 0.025),
      q975_n_info = quantile(n_inform, 0.975),
      
      prop_info = mean(inform_prop)
    ),
    by = domain
  ]
  
  # Mean arm sizes among simulations in which the decision occurred
  d_arm_summary <- d_out[
    ,
    .(
      mu_n = mean(n),
      sd_n = sd(n),
      q025_n = quantile(n, 0.025),
      q975_n = quantile(n, 0.975)
    ),
    by = .(domain, arm, inform)
  ]
  
  list(
    # for each decision
    d_inform_tot = d_inform,
    # averages by arm informing each decision
    d_arms = d_arm_summary,
    # average totals informing each decision
    d_smry = d_summary
  )
}

sim09_first_dec <- function(r) {
  
  d_out <- rbindlist(
    lapply(seq_along(r), function(i_sim) {
      rr <- r[[i_sim]]
      rbindlist(
        lapply(seq_along(rr$l_res), function(i_anlys) {
          z <- rr$l_res[[i_anlys]]
          if (is.null(z)) return(NULL)
          d_dec <- sim09_extract_dec_indicators(z$l_dec_new)
          # isTRUE() function is not vectorized and should never be used 
          # directly inside data.table rows or columns for subsetting or 
          # creating variables
          d_dec <- d_dec[dec == TRUE]
          if (nrow(d_dec) == 0) return(NULL)
          d_dec[, `:=`(i_sim = i_sim, i_anlys = i_anlys )]
          d_dec
        })
      )
    })
  )
  
  if (nrow(d_out) == 0) {
    return(data.table())
  }
  
  # Because decisions are locked, the first TRUE for a domain is
  # the interim at which that domain was actually decided.
  setorder(d_out, i_sim, domain, i_anlys)
  
  d_out[, first_dec := seq_len(.N), by = .(i_sim, domain)]
  
  d_out[first_dec == 1, .(i_sim, domain, i_anlys, rule)]
  
}

sim09_domain_n <- function(d, domain) {
  
  if (domain == "d1") {
    z <- d[silo == "l",
           .(arm = c("dair", "r1", "r2"),
             n = c(
               sum(d1 == "dair"),
               sum(d1 == "r1"),
               sum(d1 == "r2")
             ))]
    
  } else if (domain == "d2") {
    z <- d[d1 == "r1",
           .(arm = c("wk12", "wk6", "nad2"),
             n = c(
               sum(d2 == "wk12"),
               sum(d2 == "wk6"),
               sum(d2 == "nad2")
             ))]
    
  } else if (domain == "d3") {
    z <- d[d1 == "r2",
           .(arm = c("wk12", "none", "nad3"),
             n = c(
               sum(d3 == "wk12"),
               sum(d3 == "none"),
               sum(d3 == "nad3")
             ))]
    
  } else if (domain == "d4") {
    z <- d[
      ,
      .(arm = c("rif", "norif", "nad4"),
        n = c(
          sum(d4 == "rif"),
          sum(d4 == "norif"),
          sum(d4 == "nad4")
        ))
    ]
    
  } else {
    stop("Unknown domain")
  }
  
  z[, domain := domain]
  
  # Define which arms inform the treatment contrast
  z[, inform := fcase(
    domain == "d1", arm %in% c("dair", "r1", "r2"),
    domain == "d2", arm %in% c("wk12", "wk6"),
    domain == "d3", arm %in% c("wk12", "none"),
    domain == "d4", arm %in% c("rif", "norif"),
    default = FALSE
  )]
  
  z[]
}







# Replace domains allocation with an arbitrary probability vector over any 
# subset of its levels (the rest are set to 0) 
# eg. sim09_set_domain(state, "d1", c(r1 = 1/3, r2 = 2/3)).
sim09_set_domain <- function(state, domain, vec) {
  stopifnot(abs(sum(vec) - 1) < 1e-8)
  stopifnot(all(names(vec) %in% names(state[[domain]])))
  # set up a empty vector with names aligning to the original
  full <- setNames(rep(0, length(state[[domain]])), names(state[[domain]]))
  full[names(vec)] <- vec
  state[[domain]] <- full
  state
}

# Each open domain is a simplex over its possible values, which includes the 
# non-randomised option if applicable.
# - still open, default ratio: c(nad2 = 0.3, wk12 = 0.35, wk6 = 0.35)
# - still open, re-weighted (e.g. RAR): c(nad2 = 0.3, wk12 = 0.20, wk6 = 0.50)
# - closed, reverts to standard care: c(nad2 = 1,   wk12 = 0,    wk6 = 0)
# - closed, winner becomes routine: c(nad2 = 0,   wk12 = 0,    wk6 = 1)
# Ground truth (reg_effect / d4_effect below) never changes - only the
# allocation are updated.
sim09_domain_state_open <- function() {
  list(
    # 'l' silo only - randomised dair vs revision,
    d1 = c(dair = 0.50, r1 = 1/6, r2 = 1/3),   
    # revision split r1/r2 (2/3 to r2) is just best guess
    d2 = c(nad2 = 0.3, wk12 = 0.35, wk6 = 0.35),
    d3 = c(nad3 = 0.3, wk12 = 0.35, none  = 0.35),
    d4 = c(nad4 = 0.3, norif  = 0.35, rif   = 0.35)
  )
}




sim09_enrol_time_int <- function(
    N,
    lambda = function(t, lambda_inf = 1.52, ramp_up = 90) { lambda_inf * pmin(t/ramp_up, 1) },
    lambda_max = 1.52,
    ramp_up_period  = 90
) {
  
  t <- 0
  events <- numeric(0)
  while (length(events) < N) {
    # Propose next event from homogeneous PP
    t <- t + rexp(1, rate = lambda_max)
    
    # Accept with probability lambda(t) / lambda_max
    lambda_t <- lambda(t, lambda_max, ramp_up_period)
    if (runif(1) < lambda_t / lambda_max) {
      events <- c(events, t)
    }
  }
  # event times
  events
}

sim09_reg_opts <- function() {
  combos_per_silo <- c(
    "dair_nad2_nad3", 
    "r1_wk12_nad3", 
    "r1_wk6_nad3", 
    "r1_nad2_nad3",
    "r2_nad2_wk12", 
    "r2_nad2_none", 
    "r2_nad2_nad3"
  )
  silos <- c("l", "lnrd1", "enrd1", "cnrd1")
  as.vector(outer(silos, combos_per_silo, paste, sep = "_"))
}

sim09_trt_reg_contribs <- function(){
  
  l <- list()
  
  # regimes contributing to d1 effect of interest
  l$d1_trt_regs <- c(
    "l_r1_wk12_nad3",
    "l_r1_wk6_nad3",
    "l_r1_nad2_nad3",
    "l_r2_nad2_wk12",
    "l_r2_nad2_none",
    "l_r2_nad2_nad3"
  )
  
  # regimes contributing to d2 effect of interest
  l$d2_wk12_regs <- c(
    "l_r1_wk12_nad3",
    "lnrd1_r1_wk12_nad3",
    "enrd1_r1_wk12_nad3",
    "cnrd1_r1_wk12_nad3"
  )
  l$d2_wk6_regs <- c(
    "l_r1_wk6_nad3",
    "lnrd1_r1_wk6_nad3",
    "enrd1_r1_wk6_nad3",
    "cnrd1_r1_wk6_nad3"
  )
  
  # regimes contributing to d3 effect of interest
  l$d3_wk12_regs <- c(
    "l_r2_nad2_wk12",
    "lnrd1_r2_nad2_wk12",
    "enrd1_r2_nad2_wk12",
    "cnrd1_r2_nad2_wk12"
  )
  l$d3_none_regs <- c(
    "l_r2_nad2_none",
    "lnrd1_r2_nad2_none",
    "enrd1_r2_nad2_none",
    "cnrd1_r2_nad2_none"
  )
  
  l
}

# not a fan, but basic (incomplete) way to set domain effects - take care
sim09_build_reg_effect <- function(
    reg_opts, 
    d1_trt_regs,
    d2_wk6_regs,
    d3_wk12_regs,
    d1_delta = 0, 
    d2_wk6_delta = 0, 
    d3_wk12_delta = 0
) {
  
  eff <- setNames(rep(0, length(reg_opts)), reg_opts)
  eff[d1_trt_regs]  <- eff[d1_trt_regs]  + d1_delta      # any revision vs dair
  eff[d2_wk6_regs]  <- eff[d2_wk6_regs]  + d2_wk6_delta  # wk6 vs wk12 (wk12 stays at 0)
  eff[d3_wk12_regs] <- eff[d3_wk12_regs] + d3_wk12_delta # wk12 vs none (none stays at 0)
  eff["l_dair_nad2_nad3"] <- 0
  eff
}

# helper
sim09_get_silo_contrib <- function(reg_opts, prefix = "enrd1"){
  
  reg_opts[grep(prefix, reg_opts, fixed = T)]
  
}


sim09_update_cfg <- function(l_spec){
  
  if(unname(Sys.info()[1]) == "Darwin"){
    l_spec$mc_cores <- 5
    message("On mac, resetting cores to ", l_spec$mc_cores)
  } else {
    message(paste0("Allocated ", l_spec$mc_cores, " cores"))
  }
  
  l_spec$n_batch <- unlist(l_spec$n_batch)
  l_spec$p_silo <- unlist(l_spec$p_silo)
  
  l_spec$p_surg_lnrd1 <- unlist(l_spec$p_surg_lnrd1)
  l_spec$p_surg_enrd1 <- unlist(l_spec$p_surg_enrd1)
  l_spec$p_surg_cnrd1 <- unlist(l_spec$p_surg_cnrd1)
  
  # first index will always be fixed at zero irrespective of what is put
  l_spec$reg_effect <-  unlist(l_spec$reg_effect)
  l_spec$d4_effect <- unlist(l_spec$d4_effect)
  
  l_spec$n_batch <- unlist(l_spec$n_batch)
  
  l_spec$reg_opts <- sim09_reg_opts()
  stopifnot(all(l_spec$reg_opts == names(l_spec$reg_effect)))
  
  l_tmp <- sim09_trt_reg_contribs()
  l_spec$d1_trt_regs <- l_tmp$d1_trt_regs
  l_spec$d2_wk12_regs <- l_tmp$d2_wk12_regs
  l_spec$d2_wk6_regs <- l_tmp$d2_wk6_regs
  l_spec$d3_wk12_regs <- l_tmp$d3_wk12_regs
  l_spec$d3_none_regs <- l_tmp$d3_none_regs
  
  l_spec$t_0 <- sim09_enrol_time_int(sum(l_spec$n_batch))
  
  if(l_spec$nex > 0){
    l_spec$nex <- pmin(l_spec$nex, l_spec$n_sim)
    l_spec$ex_trial_ix <- sort(sample(1:l_spec$n_sim, size = l_spec$nex, replace = F))
    l_spec$ex_trial_ix[1] <- 1
  }
  
  l_spec$pri_b_0 <- unlist(l_spec$pri_b_0)
  l_spec$pri_b_reg <- unlist(l_spec$pri_b_reg)
  l_spec$pri_b_d4 <- unlist(l_spec$pri_b_d4)
  l_spec$prior_only <- as.logical(l_spec$prior_only)
  
  # hardcoded
  l_spec$mc_out_dir <- here::here("tmp") 
  
  
  l_spec
}

sim09_default_cfg <- function(){
  l_spec <- list()
  
  if(unname(Sys.info()[1]) == "Darwin"){
    l_spec$mc_cores <- 5
    message("On mac, resetting cores to ", l_spec$mc_cores)
  } else {
    
    l_spec$mc_cores <- 60
    message(paste0("Allocated ", l_spec$mc_cores, " cores"))
  }
  
  l_spec$n_sim <- 10
  
  l_spec$seed <- 1
  l_spec$mc_model <- "indep2"
  
  l_spec$p_silo <- c(l = 0.4, lnrd1 = 0.1, enrd1 = 0.3, cnrd1 = 0.2)
  l_spec$p_surg_lnrd1 <- c(dair = 0.5, r1 = 0.2, r2 = 0.3)
  l_spec$p_surg_enrd1 <- c(dair = 0.8, r1 = 0.1, r2 = 0.1)
  l_spec$p_surg_cnrd1 <- c(dair = 0.2, r1 = 0.2, r2 = 0.6)
  # regime effects of dependent domain elements
  l_spec$reg_effect <- c(
    l_dair_nad2_nad3 = 0, 
    lnrd1_dair_nad2_nad3 = 0, 
    enrd1_dair_nad2_nad3 = 0, 
    cnrd1_dair_nad2_nad3 = 0, 
    # e.g.contrib to rand surg domain 
    l_r1_wk12_nad3 = 0, 
    lnrd1_r1_wk12_nad3 = 0,
    enrd1_r1_wk12_nad3 = 0, 
    cnrd1_r1_wk12_nad3 = 0, 
    # e.g. contrib to rand surg domain 
    l_r1_wk6_nad3 = 0, 
    lnrd1_r1_wk6_nad3 = 0, 
    enrd1_r1_wk6_nad3 = 0, 
    cnrd1_r1_wk6_nad3 = 0, 
    # e.g. contrib to rand surg domain 
    l_r1_nad2_nad3 = 0, 
    lnrd1_r1_nad2_nad3 = 0, 
    enrd1_r1_nad2_nad3 = 0, 
    cnrd1_r1_nad2_nad3 = 0, 
    # e.g. contrib to rand surg domain 
    l_r2_nad2_wk12 = 0, 
    lnrd1_r2_nad2_wk12 = 0, 
    enrd1_r2_nad2_wk12 = 0, 
    cnrd1_r2_nad2_wk12 = 0, 
    # e.g. contrib to rand surg domain 
    l_r2_nad2_none = 0, 
    lnrd1_r2_nad2_none = 0, 
    enrd1_r2_nad2_none = 0, 
    cnrd1_r2_nad2_none = 0, 
    # e.g. contrib to rand surg domain 
    l_r2_nad2_nad3 = 0, 
    lnrd1_r2_nad2_nad3 = 0, 
    enrd1_r2_nad2_nad3 = 0, 
    cnrd1_r2_nad2_nad3 = 0
    )
  l_spec$d4_effect <- c(nad4 = 0, norif = 0, rif = 0)
  l_spec$response_p_ref <- 0.6
  
  l_spec$desc <- "Default cfg"
  l_spec$nex <- 3
  l_spec$ex_trial_ix <- c(1, 2, 3)
  l_spec$n_batch <- c(500, 500, 500, 500, 500)
  
  l_spec$reg_opts <- sim09_reg_opts()
  stopifnot(all(l_spec$reg_opts == names(l_spec$reg_effect)))
  
  l_tmp <- sim09_trt_reg_contribs()
  
  l_spec$d1_trt_regs <- l_tmp$d1_trt_regs
  l_spec$d2_wk12_regs <- l_tmp$d2_wk12_regs
  l_spec$d2_wk6_regs <- l_tmp$d2_wk6_regs
  l_spec$d3_wk12_regs <- l_tmp$d3_wk12_regs
  l_spec$d3_none_regs <- l_tmp$d3_none_regs

  l_spec$t_0 <- sim09_enrol_time_int(sum(l_spec$n_batch))
  
  l_spec$pri_b_0 = c(0, 1)
  l_spec$pri_b_reg = c(0, 1)
  l_spec$pri_b_d4 = c(0, 1)
  l_spec$prior_only <- FALSE
  
  l_spec$mc_warmup <-  1000
  l_spec$mc_samp <-  1000
  l_spec$mc_chain <-  1
  
  l_spec$mc_out_dir <- here::here("tmp") 
  
  l_spec$dec <- list()
  l_spec$dec$d1 <- list()
  l_spec$dec$d2 <- list()
  l_spec$dec$d3 <- list()
  l_spec$dec$d4 <- list()
  
  # domain rules and values
  
  # theta = pr_succsess_rev - pr_succsess_dair
  # Pr(theta > delta)  > thresh => Superiority
  # Pr(theta > delta) < thresh => Futility, ie futile if Pr theta > 0.05 is < 0.3
  l_spec$dec$d1$sup <- list(delta = 0, thresh = 0.95)
  l_spec$dec$d1$fut <- list(delta = 0.05, thresh = 0.3)
  
  # theta = pr_succsess_wk6 - pr_succsess_wk12
  # Pr(theta > delta) > thresh => Non-inferior
  # Pr(theta > delta) < thresh => Futility, ie futile if Pr theta > 0.0 is < 0.1
  l_spec$dec$d2$ni <- list(delta = -0.05, thresh = 0.95)
  l_spec$dec$d2$fut <- list(delta = 0.0, thresh = 0.1)
  
  # theta = pr_succsess_wk12 - pr_succsess_wknone
  l_spec$dec$d3$sup <- list(delta = 0, thresh = 0.95)
  l_spec$dec$d3$fut <- list(delta = 0.05, thresh = 0.3)
  
  # theta = pr_succsess_rif - pr_succsess_norif
  l_spec$dec$d4$sup <- list(delta = 0, thresh = 0.99)
  l_spec$dec$d4$fut <- list(delta = 0.05, thresh = 0.3)
  
  l_spec
}

# Ex DGP figs ---------
sim09_ex_dat_1 <- function(){
  
  # CFG
  default_cfg <- T
  if(!default_cfg){
    f_cfgsc <- file.path("./etc/sim09/cfg-sim09-sc01-v01.yml")
    l_spec <- config::get(file = f_cfgsc)
    l_spec <- sim09_update_cfg(l_spec)
  } else {
    l_spec <- sim09_default_cfg()
  }
  # coordinates size of batch
  l_spec$is <- 1
  l_spec$ie <- sum(l_spec$n_batch)
  # starting state for domains
  l_dom_state = sim09_domain_state_open()
  
  d_batch <- sim09_batch_01(l_spec, l_dom_state = l_dom_state)
  
  ## figs ------------
  d_fig <- d_batch[silo == "l"]
  d_fig[, d1_b := copy(d1)]
  d_fig[d1 != "dair", d1_b := "rev"]
  p_1 <- ggplot( d_fig, aes(x = d1_b, fill = d1)) + geom_bar() + scale_x_discrete("") +
    ggtitle(
      "Surgical (late silo only)",
      subtitle = paste0("rand trt N = ", nrow(d_fig), "/", nrow(d_batch))
    )
  
  d_fig <- d_batch[d1 == "r1"]
  p_2 <- ggplot( d_fig, aes(x = d2, fill = silo)) + geom_bar() + scale_x_discrete("") +
    ggtitle(
      "Duration A (r1 rev all silo)",
      subtitle = paste0("rand trt N = ", nrow(d_fig[d2 != "nad2"]), "/", nrow(d_batch))
    )
  
  d_fig <- d_batch[d1 == "r2"]
  p_3 <- ggplot( d_fig, aes(x = d3, fill = silo)) + geom_bar() + scale_x_discrete("") +
    ggtitle(
      "Duration B (r2 rev all silo)",
      subtitle = paste0("rand trt N = ", nrow(d_fig[d3 != "nad3"]), "/", nrow(d_batch))
      )
  
  d_fig <- copy(d_batch)
  p_4 <- ggplot( d_fig, aes(x = d4, fill = silo)) + geom_bar() + scale_x_discrete("") +
    ggtitle(
      "Choice (applicable pathogen all silo)",
      subtitle = paste0("rand trt N = ", nrow(d_fig[d4 != "nad4"]), "/", nrow(d_batch))
      )
  
  p_1 + p_2 + p_3 + p_4
  
  # 
  d_fig <- copy(d_batch)
  d_fig <- d_fig[, .(count = .N), keyby = .(silo, reg)]
  d_fig[, prop := count / sum(count), keyby = silo]
  d_fig[, n := sum(count), keyby = silo]
  d_fig[, silo_lab := paste0(silo, ", N = ", n)]
  sort(unique(d_fig$silo_lab))
  d_fig[, silo_lab := factor(silo_lab, levels = c(
    "enrd1, N = 773", "l, N = 992", "lnrd1, N = 276",
    "cnrd1, N = 459"
  ))]
  p_1 <- ggplot(d_fig, aes(
    x = reg, y = prop)) + 
    geom_col() + scale_x_discrete("") +
    geom_text(aes(label = count), col = "red") +
    coord_flip() +
    facet_wrap(~ silo_lab, scales = "free_y", labeller = label_both) +
    ggtitle(
      "Proportion allocated to each regimen"
    )
  p_1
#  
  
  
  
}





# Ex prototype/minimal interim setup --------
sim09_ex_sim_1 <- function(
    # easily switch out to some different function when I build in analysis code
    fn_decision = sim09_decision_fn_dummy
  ){
  
  set.seed(1)
  
  
  # CFG
  default_cfg <- F
  if(!default_cfg){
    f_cfgsc <- file.path("./etc/sim09/cfg-sim09-sc01-v01.yml")
    l_spec <- config::get(file = f_cfgsc)
    l_spec <- sim09_update_cfg(l_spec)
  } else {
    l_spec <- sim09_default_cfg()
  }
  # starting state for domains
  l_dom_state = sim09_domain_state_open()
  
  # accrued data
  d_cum_dat   <- data.table()
  state_log <- vector("list", length(l_spec$n_batch))
  
  i <- 1
  for (i in seq_along(l_spec$n_batch)) {
    
    if(i == 1){
      # starting pt index in data
      l_spec$is <- 1
      l_spec$ie <- l_spec$is + l_spec$n_batch[i] - 1
    } else {
      l_spec$is <- nrow(d_cum_dat) + 1
      l_spec$ie <- l_spec$is + l_spec$n_batch[i] - 1
    }
    
    d_batch_dat <- sim09_batch_01(l_spec, l_dom_state = l_dom_state)
    d_batch_dat[, batch := i]
    d_cum_dat <- rbind(d_cum_dat, d_batch_dat)
    
    # state that generated batch i
    state_log[[i]] <- l_dom_state         
    # updated for batch i+1
    l_dom_state   <- fn_decision(d_cum_dat, l_dom_state, i)   
  }
  
  list(data = d_cum_dat, state_log = state_log, final_state = l_dom_state)
  
  
}


# Ex MV model vs domain level models -----
sim09_ex_sim_2 <- function(
    l_spec, l_dom_state
    ){
  
  
  d_res <- rbindlist(pbapply::pblapply(
    X=1:l_spec$n_sim, cl = l_spec$mc_cores, FUN=function(ix) {
      
      d_batch <- sim09_batch_01(l_spec, l_dom_state = l_dom_state)
      
      # multivariate (joint model) handling all domains at once.
      X <- model.matrix(~ reg + d4, data = d_batch)
      f_1 <- fastglm::fastglm(X, d_batch$y , family = binomial)
      
      # proportions with which we will weight params fur d1
      d_w_b <- d_batch[reg %in% l_spec$d1_trt_regs, .N, keyby = reg]
      d_w_b[, w := N/sum(d_w_b$N)]
      setkey(d_w_b, reg)
      wgt <- d_w_b[l_spec$d1_trt_regs, w]
      f_1_coef <- coef(f_1)[paste0("reg", l_spec$d1_trt_regs)]
      jnt_d1 <- as.numeric(wgt %*% f_1_coef)
      
      # proportions with which we will weight params fur d2
      d_w_a <- d_batch[reg %in% l_spec$d2_wk12_regs, .N, keyby = reg]
      d_w_a[, w := N/sum(d_w_a$N)]
      setkey(d_w_a, reg)
      d_w_b <- d_batch[reg %in% l_spec$d2_wk6_regs, .N, keyby = reg]
      d_w_b[, w := N/sum(d_w_b$N)]
      setkey(d_w_b, reg)
      wgt_a <- d_w_a[l_spec$d2_wk12_regs, w]
      wgt_b <- d_w_b[l_spec$d2_wk6_regs, w]
      f_1_coef_a <- coef(f_1)[paste0("reg", l_spec$d2_wk12_regs)]
      f_1_coef_b <- coef(f_1)[paste0("reg", l_spec$d2_wk6_regs)]
      jnt_d2 <- as.numeric((wgt_b %*% f_1_coef_b) - (wgt_a %*% f_1_coef_a))
      
      # proportions with which we will weight params fur d3
      d_w_a <- d_batch[reg %in% l_spec$d3_none_regs, .N, keyby = reg]
      d_w_a[, w := N/sum(d_w_a$N)]
      setkey(d_w_a, reg)
      d_w_b <- d_batch[reg %in% l_spec$d3_wk12_regs, .N, keyby = reg]
      d_w_b[, w := N/sum(d_w_b$N)]
      setkey(d_w_b, reg)
      wgt_a <- d_w_a[l_spec$d3_none_regs, w]
      wgt_b <- d_w_b[l_spec$d3_wk12_regs, w]
      f_1_coef_a <- coef(f_1)[paste0("reg", l_spec$d3_none_regs)]
      f_1_coef_b <- coef(f_1)[paste0("reg", l_spec$d3_wk12_regs)]
      jnt_d3 <- as.numeric((wgt_b %*% f_1_coef_b) - (wgt_a %*% f_1_coef_a))
      
      # get this for free
      jnt_d4 <- as.numeric(coef(f_1)["d4rif"] - coef(f_1)["d4norif"])
      
      # separate model fits which target the equivalent trt effects
      d_mod <- d_batch[silo == "l"]
      X <- model.matrix(~ d1b, data = d_mod)
      f_2 <- fastglm::fastglm(X, d_mod$y , family = binomial)
      
      uni_d1 <- coef(f_2)[2] 
      
      d_mod <- d_batch[d1 == "r1"]
      X <- model.matrix(~ d2, data = d_mod)
      f_2 <- fastglm::fastglm(X, d_mod$y , family = binomial)
      
      uni_d2 <- coef(f_2)["d2wk6"] - coef(f_2)["d2wk12"]
      
      d_mod <- d_batch[d1 == "r2"]
      X <- model.matrix(~ d3, data = d_mod)
      f_2 <- fastglm::fastglm(X, d_mod$y , family = binomial)
      
      uni_d3 <- coef(f_2)["d3wk12"] - coef(f_2)["d3none"]
      
      X <- model.matrix(~ d4, data = d_batch)
      f_2 <- fastglm::fastglm(X, d_batch$y , family = binomial)
      
      uni_d4 <- coef(f_2)["d4rif"] - coef(f_2)["d4norif"]
      
      data.table(
        jnt_d1 = jnt_d1, 
        uni_d1 = uni_d1, 
        jnt_d2 = jnt_d2, 
        uni_d2 = uni_d2, 
        jnt_d3 = jnt_d3, 
        uni_d3 = uni_d3, 
        jnt_d4 = jnt_d4, 
        uni_d4 = uni_d4
      )
      
    }
  ))
  
  d_tbl <- melt(d_res, measure.vars = names(d_res))
  d_tbl[, c("model", "domain") := tstrsplit(variable, "_", fixed = T)]
  d_smry <- dcast(
    d_tbl[, .(
      mu = mean(value), 
      sd = sd(value)
    ), keyby = .(model, domain)], domain ~ model, value.var = list("mu", "sd"))
  
  d_smry
  
}

# Ex g-comp test ------
sim09_ex_gcomp_stan_demo <- function(
){
  m_2 <- cmdstanr::cmdstan_model(here::here("stan", "model-sim-09-gcomp.stan"))
  
  d <- data.table(id = 1:1000)
  d[, trt := rbinom(.N, 1, 0.5) + 1]
  d[, cov1 := sample(1:4, size = .N, replace = T)]
  d[, cov2 := sample(1:2, size = .N, replace = T)]
  b_trt = c(0, 1)
  b_cov1 = c(0, 0.5, 0.2, 0.1)
  b_cov2 = c(0, -0.1)
  d[, eta := qlogis(0.4) + b_trt[trt] + b_cov1[cov1] + b_cov2[cov2]]
  d[, y := rbinom(.N, 1, plogis(eta))]
  
  # d[, .(y = sum(y), n = .N), by = .(trt, cov1, cov2)]
  
  ld <- list(
    N = nrow(d), y = d$y, 
    K_trt = length(unique(d$trt)),
    K_cov1 = length(unique(d$cov1 )),
    K_cov2 = length(unique(d$cov2)),
    trt = d$trt, cov1 = d$cov1, cov2 = d$cov2,
    pri_b_0 = c(0, 1),
    pri_b_trt = c(0, 1), pri_b_cov1 = c(0, 1), pri_b_cov2 = c(0, 1)
  )
  
  f_2 <- m_2$sample(
    ld, iter_warmup = 1000, iter_sampling = 1000,
    parallel_chains = 1, chains = 1,
    refresh = 0, show_exceptions = T
  )
  f_2$summary(variables = c(
    "b_0", "b_trt", "b_cov1", "b_cov2"
  ))
  f_2$summary(variables = c(
    "theta_1", "mu_1", "theta_2", "mu_2"
  ))
  
  
  
  
}


# Model fit scenarios ----
# test whether the joint model and univariate model applied to the relevant
# population can recover the same point value for the effects by weighting
# the relevatn regime parameters by their proportional representation in the 
# sample data
sim09_ex_scenarios <- function(){
  
  set.seed(1)
  default_cfg <- T
  if(!default_cfg){
    f_cfgsc <- file.path("./etc/sim09/cfg-sim09-sc01-v01.yml")
    l_spec <- config::get(file = f_cfgsc)
    l_spec <- sim09_update_cfg(l_spec)
  } else {
    l_spec <- sim09_default_cfg()
  }
  # coordinates size of batch
  l_spec$n_batch <- c(2500)
  l_spec$is <- 1
  l_spec$ie <- sum(l_spec$n_batch)
  l_spec$n_sim <- 10000
  
  # starting state for domains
  l_dom_state = sim09_domain_state_open()
  
  # Scenario - positive effect restricted to d1 ------------
  
  # contrive an effect solely attributable to revision...
  # set all revision effects, no sub domain effects
  # these are all relative to l_dair_nad2_nad3 holding d4 constant
  # l_spec$reg_opts[grep("l_", l_spec$reg_opts , fixed = T)]
  l_spec$reg_effect <- sim09_build_reg_effect(
    reg_opts = l_spec$reg_opts, 
    d1_trt_regs = l_spec$d1_trt_regs,
    d2_wk6_regs = l_spec$d2_wk6_regs,
    d3_wk12_regs = l_spec$d3_wk12_regs,
    d1_delta = 1, 
    d2_wk6_delta = 0, 
    d3_wk12_delta = 0
  )
  
  d_tbl <- sim09_ex_sim_2(l_spec, l_dom_state)
  kableExtra::kbl(
    d_tbl, digits = 3, format = "simple", 
    caption = paste0(
      "n_sim ", l_spec$n_sim, " n ", sum(l_spec$n_batch)))
  # Table: n_sim 10000 n 2500
  # 
  # domain    mu_jnt   mu_uni   sd_jnt   sd_uni
  # d1         1.025    1.003    0.148    0.145
  # d2         0.000   -0.001    0.298    0.263
  # d3         0.000    0.000    0.202    0.185
  # d4         0.001    0.001    0.103    0.099
  
  
  
  # Scenario - negative effect associ with d2 wk6 pollutes d1 ------------
  # just an effect in d2
  # l_spec$reg_opts[grep("l_", l_spec$reg_opts , fixed = T)]
  l_spec$test_fx <- c(
    
    # d2 wk6 reduces log odds trt success
    # the reference group is not diff from l_dair_nad2_nad3
    lnrd1_r1_wk12_nad3 = 0.0,
    # deleterious
    lnrd1_r1_wk6_nad3 = -1, 
    
    enrd1_r1_wk12_nad3 = 0.0, 
    enrd1_r1_wk6_nad3 = -1, 
    
    cnrd1_r1_wk12_nad3 = 0.0, 
    cnrd1_r1_wk6_nad3 = -1, 
    
    # d3 wk12 increases log odds trt success
    
    lnrd1_r2_nad2_none = 0.0,
    enrd1_r2_nad2_none = 0.0, 
    cnrd1_r2_nad2_none = 0.0,
    
    lnrd1_r2_nad2_wk12 = 0.0,
    enrd1_r2_nad2_wk12 = 0.0,
    cnrd1_r2_nad2_wk12 = 0.0,
    
    # d1 contributions
    l_r1_wk12_nad3 = 0.0,
    l_r1_wk6_nad3 = -1, 
    # assume that non-randomised trt defaults to 12wks
    l_r1_nad2_nad3 = 0.0,
    
    l_r2_nad2_wk12 = 0.0,
    l_r2_nad2_none = 0.0,
    l_r2_nad2_nad3 = 0.0
    
  )
  full_reg_effects <- setNames(rep(0, length(l_spec$reg_effect)), l_spec$reg_opts)
  full_reg_effects[names(l_spec$test_fx)] <- l_spec$test_fx
  l_spec$reg_effect <- full_reg_effects
  
  d_tbl <- sim09_ex_sim_2(l_spec, l_dom_state)
  kableExtra::kbl(
    d_tbl, digits = 3, format = "simple", 
    caption = paste0(
      "n_sim ", l_spec$n_sim, " n ", sum(l_spec$n_batch)))
  # Table: n_sim 10000 n 2500
  # 
  # domain    mu_jnt   mu_uni   sd_jnt   sd_uni
  # d1        -0.113   -0.116    0.131    0.129
  # d2        -1.036   -1.006    0.265    0.249
  # d3        -0.002   -0.002    0.178    0.175
  # d4         0.002    0.001    0.100    0.098
  
  # Scenario - positive effect associ with 12wk pollutes d1 ------------
  # just an effect in d3 
  # l_spec$reg_opts[grep("l_", l_spec$reg_opts , fixed = T)]
  l_spec$test_fx <- c(
    
    # d2 wk6 reduces log odds trt success
    # the reference group is not diff from l_dair_nad2_nad3
    lnrd1_r1_wk12_nad3 = 0.0,
    # deleterious
    lnrd1_r1_wk6_nad3 = 0.0, 
    
    enrd1_r1_wk12_nad3 = 0.0, 
    enrd1_r1_wk6_nad3 = 0.0,
    
    cnrd1_r1_wk12_nad3 = 0.0, 
    cnrd1_r1_wk6_nad3 = 0.0,
    
    # d3 wk12 increases log odds trt success
    
    lnrd1_r2_nad2_none = 0.0,
    enrd1_r2_nad2_none = 0.0, 
    cnrd1_r2_nad2_none = 0.0,
    
    lnrd1_r2_nad2_wk12 = 1.0,
    enrd1_r2_nad2_wk12 = 1.0,
    cnrd1_r2_nad2_wk12 = 1.0,
    
    # d1 contributions
    l_r1_wk12_nad3 = 0.0,
    l_r1_wk6_nad3 = 0.0,
    # assume that non-randomised trt defaults to 12wks
    l_r1_nad2_nad3 = 0.0,
    
    l_r2_nad2_wk12 = 1.0,
    l_r2_nad2_none = 0.0,
    l_r2_nad2_nad3 = 0.0
    
  )
  full_reg_effects <- setNames(rep(0, length(l_spec$reg_effect)), l_spec$reg_opts)
  full_reg_effects[names(l_spec$test_fx)] <- l_spec$test_fx
  l_spec$reg_effect <- full_reg_effects
  
  d_tbl <- sim09_ex_sim_2(l_spec, l_dom_state)
  kableExtra::kbl(
    d_tbl, digits = 3, format = "simple", 
    caption = paste0(
      "n_sim ", l_spec$n_sim, " n ", sum(l_spec$n_batch)))
  # Table: n_sim 10000 n 2500
  # 
  # domain    mu_jnt   mu_uni   sd_jnt   sd_uni
  # d1         0.240    0.201    0.136    0.131
  # d2         0.004    0.004    0.264    0.249
  # d3         1.036    1.005    0.235    0.196
  # d4        -0.002   -0.002    0.101    0.099
  
  # Scenario - positive effect associ with rif ------------
  # just an effect in d3 
  # l_spec$reg_opts[grep("l_", l_spec$reg_opts , fixed = T)]
  l_spec$test_fx <- c()
  l_spec$d4_effect["rif"]  <- 1
  full_reg_effects <- setNames(rep(0, length(l_spec$reg_effect)), l_spec$reg_opts)
  full_reg_effects[names(l_spec$test_fx)] <- l_spec$test_fx
  l_spec$reg_effect <- full_reg_effects
  
  d_tbl <- sim09_ex_sim_2(l_spec, l_dom_state)
  kableExtra::kbl(
    d_tbl, digits = 3, format = "simple", 
    caption = paste0(
      "n_sim ", l_spec$n_sim, " n ", sum(l_spec$n_batch)))
  # Table: n_sim 10000 n 2500
  # 
  # domain    mu_jnt   mu_uni   sd_jnt   sd_uni
  # d1         0.007   -0.001    0.139    0.135
  # d2         0.002    0.003    0.301    0.262
  # d3        -0.005   -0.006    0.190    0.183
  # d4         1.014    1.004    0.110    0.108
  
  
  # Scenario - opposite effects in d2/d3 ------------
  # just an effect in d3 
  # l_spec$reg_opts[grep("l_", l_spec$reg_opts , fixed = T)]
  l_spec$test_fx <- c(
    
    # d2 wk6 reduces log odds trt success
    # the reference group is not diff from l_dair_nad2_nad3
    # assume that non-randomised trt defaults to 12wks
    l_r1_wk12_nad3 = 0.0,
    l_r1_wk6_nad3 = -0.75,
    
    lnrd1_r1_wk12_nad3 = 0.0,
    # deleterious
    lnrd1_r1_wk6_nad3 = -0.75,
    
    enrd1_r1_wk12_nad3 = 0.0, 
    enrd1_r1_wk6_nad3 = -0.75,
    
    cnrd1_r1_wk12_nad3 = 0.0, 
    cnrd1_r1_wk6_nad3 = -0.75,
    
    # d3 wk12 increases log odds trt success
    
    l_r2_nad2_none = 0.0,
    l_r2_nad2_wk12 = 1.0,
    
    lnrd1_r2_nad2_none = 0.0,
    lnrd1_r2_nad2_wk12 = 1.0,
    
    enrd1_r2_nad2_none = 0.0, 
    enrd1_r2_nad2_wk12 = 1.0,
    
    cnrd1_r2_nad2_none = 0.0,
    cnrd1_r2_nad2_wk12 = 1.0
  )
  l_spec$d4_effect["rif"]  <- 0
  full_reg_effects <- setNames(rep(0, length(l_spec$reg_effect)), l_spec$reg_opts)
  full_reg_effects[names(l_spec$test_fx)] <- l_spec$test_fx
  l_spec$reg_effect <- full_reg_effects
  
  d_tbl <- sim09_ex_sim_2(l_spec, l_dom_state)
  kableExtra::kbl(
    d_tbl, digits = 3, format = "simple", 
    caption = paste0(
      "n_sim ", l_spec$n_sim, " n ", sum(l_spec$n_batch)))
  # Table: n_sim 10000 n 2500
  # 
  # domain    mu_jnt   mu_uni   sd_jnt   sd_uni
  # d1         0.154    0.111    0.135    0.129
  # d2        -0.778   -0.756    0.262    0.249
  # d3         1.034    1.005    0.230    0.195
  # d4         0.000    0.000    0.100    0.097
  
  
  
}







# more complete setup of nested design based on what is likely to happen in 
# trial and bar the inclusion of things like site, joint, prognostic factors
# etc.
# just uses a linear (not logistic) model for convenience
sim09_ex_partial_nest_01 <- function(){
  
  library(data.table)
  set.seed(1)
  
  d_cells <- rbind(
    # late silo, randomised d1 as dair vs rev but clincician selects r1/r2
    # no nad1 as some form of surgery is necessary to proceed to any other form
    # of treatment
    data.table(silo = "l", d1 = "dair", d2 = "nad2", d3 = "nad3"),
    data.table(silo = "l", d1 = "r1", d2 = "wk12", d3 = "nad3"),
    data.table(silo = "l", d1 = "r1", d2 = "wk6",  d3 = "nad3"),
    # d2 can receive non-rand trt even though r1 was revision type
    data.table(silo = "l", d1 = "r1", d2 = "nad2",  d3 = "nad3"),
    data.table(silo = "l", d1 = "r2", d2 = "nad2", d3 = "wk12"),
    data.table(silo = "l", d1 = "r2", d2 = "nad2", d3 = "none"),
    # d3 can receive non-rand trt even though r2 was revision type
    data.table(silo = "l", d1 = "r2", d2 = "nad2", d3 = "nad3"),
    
    # late silo, non-randomised d1, clincician selects dair/r1/r2
    data.table(silo = "lnrd1", d1 = "dair", d2 = "nad2", d3 = "nad3"),
    data.table(silo = "lnrd1", d1 = "r1", d2 = "wk12", d3 = "nad3"),
    data.table(silo = "lnrd1", d1 = "r1", d2 = "wk6",  d3 = "nad3"),
    data.table(silo = "lnrd1", d1 = "r1", d2 = "nad2",  d3 = "nad3"),
    data.table(silo = "lnrd1", d1 = "r2", d2 = "nad2", d3 = "wk12"),
    data.table(silo = "lnrd1", d1 = "r2", d2 = "nad2", d3 = "none"),
    data.table(silo = "lnrd1", d1 = "r2", d2 = "nad2", d3 = "nad3"),
    
    # early silo, non-rand d1
    data.table(silo = "enrd1", d1 = "dair", d2 = "nad2", d3 = "nad3"),
    data.table(silo = "enrd1", d1 = "r1", d2 = "wk12", d3 = "nad3"),
    data.table(silo = "enrd1", d1 = "r1", d2 = "wk6",  d3 = "nad3"),
    data.table(silo = "enrd1", d1 = "r1", d2 = "nad2",  d3 = "nad3"),
    data.table(silo = "enrd1", d1 = "r2", d2 = "nad2", d3 = "wk12"),
    data.table(silo = "enrd1", d1 = "r2", d2 = "nad2", d3 = "none"),
    data.table(silo = "enrd1", d1 = "r2", d2 = "nad2", d3 = "nad3"),
    
    # chronic silo, non-rand d1
    data.table(silo = "cnrd1", d1 = "dair", d2 = "nad2", d3 = "nad3"),
    data.table(silo = "cnrd1", d1 = "r1", d2 = "wk12", d3 = "nad3"),
    data.table(silo = "cnrd1", d1 = "r1", d2 = "wk6", d3 = "nad3"),
    data.table(silo = "cnrd1", d1 = "r1", d2 = "nad2", d3 = "nad3"),
    data.table(silo = "cnrd1", d1 = "r2", d2 = "nad2", d3 = "wk12"),
    data.table(silo = "cnrd1", d1 = "r2", d2 = "nad2", d3 = "none"),
    data.table(silo = "cnrd1", d1 = "r2", d2 = "nad2", d3 = "nad3")
  )
  d_cells <- d_cells[rep(1:.N, each = 3), ]
  # nad4 exists to demarcate non-randomised pts, i.e. pt inelig for d4
  d_cells[, d4 := rep(c("rif", "norif", "nad4"), length = .N)]
  
  
  dat <- d_cells[rep(seq_len(nrow(d_cells)), each = 4), ]
  dat$silo <- factor(dat$silo, levels = c("l", "lnrd1", "enrd1", "cnrd1"))
  dat$d1 <- factor(dat$d1, levels = c("dair", "r1", "r2"))
  dat$d2 <- factor(dat$d2, levels = c("nad2", "wk12", "wk6"))
  dat$d3 <- factor(dat$d3, levels = c("nad3", "wk12", "none"))
  dat$d4 <- factor(dat$d4, levels = c("nad4", "norif", "rif"))
  dat$y <- rnorm(nrow(dat))
  
  # ---- 1. Naive coding: cross Surg:d2 and Surg:d3 without respecting that
  f_1 <- lm(y ~ silo * d1 + d1:d2 + d1:d3 + d4, data = dat)
  X_1 <- model.matrix(f_1)
  ncol(X_1)
  # identifiable
  qr(X_1)$rank
  names(coef(f_1))[is.na(coef(f_1))]
  # rows indicate linear combination of cols
  alias(f_1)  
  summary(f_1)
  
  # Explicit nested coding by collapsing silo, Surg, d2, d3 into the
  # regimen combinations keeping the non-rand and rand elements separate as
  # best we can
  dat$reg <- with(dat, interaction(silo, d1, d2, d3, drop = TRUE, sep = "_"))
  levels(dat$reg)
  nlevels(dat$reg)
  
  f_2 <- lm(y ~ reg + d4, data = dat)
  X_2 <- model.matrix(f_2)
  ncol(X_2)
  qr(X_2)$rank
  summary(f_2)
  
  # Coefficients:
  #   Estimate Std. Error t value Pr(>|t|)  
  # (Intercept)              0.35426    0.28560   1.240   0.2158  
  # reglnrd1_dair_nad2_nad3  0.31421    0.39020   0.805   0.4213  
  # regenrd1_dair_nad2_nad3  0.36355    0.39020   0.932   0.3522  
  # regcnrd1_dair_nad2_nad3 -0.20377    0.39020  -0.522   0.6019  
  # regl_r1_nad2_nad3       -0.10345    0.39020  -0.265   0.7911  
  # reglnrd1_r1_nad2_nad3   -0.31784    0.39020  -0.815   0.4160  
  # regenrd1_r1_nad2_nad3   -0.02606    0.39020  -0.067   0.9468  
  # regcnrd1_r1_nad2_nad3    0.02472    0.39020   0.063   0.9495  
  # regl_r2_nad2_nad3       -0.54299    0.39020  -1.392   0.1651  
  # reglnrd1_r2_nad2_nad3   -0.02678    0.39020  -0.069   0.9453  
  # regenrd1_r2_nad2_nad3   -0.19812    0.39020  -0.508   0.6120  
  # regcnrd1_r2_nad2_nad3   -0.26283    0.39020  -0.674   0.5011  
  # regl_r1_wk12_nad3       -0.23754    0.39020  -0.609   0.5431  
  # reglnrd1_r1_wk12_nad3   -0.44692    0.39020  -1.145   0.2529  
  # regenrd1_r1_wk12_nad3   -0.39909    0.39020  -1.023   0.3072  
  # regcnrd1_r1_wk12_nad3   -0.07476    0.39020  -0.192   0.8482  
  # regl_r1_wk6_nad3        -0.37908    0.39020  -0.972   0.3321  
  # reglnrd1_r1_wk6_nad3    -0.19629    0.39020  -0.503   0.6153  
  # regenrd1_r1_wk6_nad3    -0.13091    0.39020  -0.335   0.7375  
  # regcnrd1_r1_wk6_nad3    -0.48827    0.39020  -1.251   0.2118  
  # regl_r2_nad2_wk12       -0.08503    0.39020  -0.218   0.8276  
  # reglnrd1_r2_nad2_wk12   -0.82041    0.39020  -2.103   0.0363 *
  # regenrd1_r2_nad2_wk12   -0.68172    0.39020  -1.747   0.0816 .
  # regcnrd1_r2_nad2_wk12    0.20568    0.39020   0.527   0.5985  
  # regl_r2_nad2_none        0.08783    0.39020   0.225   0.8221  
  # reglnrd1_r2_nad2_none   -0.77427    0.39020  -1.984   0.0481 *
  # regenrd1_r2_nad2_none   -0.53334    0.39020  -1.367   0.1727  
  # regcnrd1_r2_nad2_none   -0.24627    0.39020  -0.631   0.5284  
  # d4norif                 -0.12456    0.12772  -0.975   0.3302  
  # d4rif                   -0.13230    0.12772  -1.036   0.3011  
  
  # d1:
  # contrast between:
  # log odds response dair:
  # intercept (b0) is rand dair with non rand d2, non rand d3 and non rand d4
  
  # log odds response revision  (which also assumes non rand d4)
  # q %*% [ (b0 + regl_r1_nad2_nad3), (b0 + regl_r1_wk12_nad3), (b0 + regl_r1_wk6_nad3), 
  #         (b0 + regl_r2_nad2_nad3), (b0 + regl_r2_nad2_wk12), (b0 + regl_r2_nad2_none) ]'
  # q is 6x1 vector of weights based on the observed distribution of 
  # membership in each one of these regimens
  
  # could just use above as the contrast on the log odds scale or do g-comp
  # to translate into a risk diff.
  
  # the contrast on the log odds scale should match the effect of revision 
  # obtained from a d1 specific model fit to the appropriate set of data,
  # i.e. late silo pts who received randomised trt for d1
  
  # d2: 
  # Estimand: The effect of assignment to 12-week versus 6-week antibiotic 
  # duration in patients who receive a one-stage revision and are eligible for 
  # antibiotic domain on 12-month outcome (from trial entry) under their prior 
  # treatment pathway and including participant preference for downstream 
  # participation and conditionally randomised rifampicin.
  
  # contrast between:
  # wk12:
  # w [ b0 + regl_r1_wk12_nad3, b0 + reglnrd1_r1_wk12_nad3, 
  #     b0 + regenrd1_r1_wk12_nad3, b0 + regcnrd1_r1_wk12_nad3 ]
  # 
  # w is based on obs dis in relevant regs
  #
  # wk6:
  # v [ b0 + regl_r1_wk6_nad3, b0 + reglnrd1_r1_wk6_nad3, 
  #     b0 + regenrd1_r1_wk6_nad3, b0 + regcnrd1_r1_wk6_nad3 ]
  #
  # v is based on obs dis in relevant regs
  
  # the contrast on the log odds scale should match the d2 trt effect 
  # obtained from a d2 specific model fit to the appropriate subset of data,
  # e.g. r1 pts entering into rand trt for d2
  
  
}




# MAIN SIM LOOP -------------
sim09_sim_loop <- function(){
  
  log_info(paste0(match.call()[[1]]))
  
  default_cfg <- F
  if(!default_cfg){
    # load sim specification
    f_spec <- here::here("./etc", args[2])
    l_spec <- config::get(file = f_spec)
    stopifnot("Config is null" = !is.null(l_spec))
    l_spec <- sim09_update_cfg(l_spec)
  } else {
    l_spec <- sim09_default_cfg()
  }
  
  # str(l_spec)
  l_spec$return_posterior = F  ; e = NULL; ix <- 1
  log_info("Starting simulation")
  
  # temp
  l_dom_state = sim09_domain_state_open()
  
  RNGkind("L'Ecuyer-CMRG"); set.seed(l_spec$seed)
  r <- pbapply::pblapply(
    X=1:l_spec$n_sim, cl = l_spec$mc_cores, FUN=function(ix) {
      
      log_info("Simulation ", ix);
      
      l_spec$ix_sim <- ix
      
      if(ix %in% l_spec$ex_trial_ix){ l_spec$return_posterior = T  
      } else { l_spec$return_posterior = F }
      
      ll <- tryCatch({
        sim09_run_trial(
          l_spec,
          l_dom_state,
          sim09_decision_fn_01
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
  log_info("Sleep for 2 secs before processing")
  Sys.sleep(2)
  
  scen <- substr(basename(f_spec), 16, 18)
  fname <- paste0("sim09-", scen, "-", format(Sys.time(), "%Y%m%d-%H%M%S"), ".qs2")
  log_info("sim09_sim_loop: saving to file", fname)
  qs2::qs_save(
    list(
      r = r,
      l_spec = l_spec, 
      l_dom_state0 = sim09_domain_state_open()
      ),
    file = here::here("data", "sim09", fname)
  )
  
}

sim09_run_none <- function(){
  log_info("sim09_run_none: Nothing doing here bud.")
}

sim09_main <- function(){
  funcname <- paste0(args[1], "()")
  log_info("Main, invoking ", funcname)
  eval(parse(text=funcname))
}

if(!interactive()){
  sim09_main()
}




# stan models----------------
mod_a <- "
data{ 
  int N;
  array[N] int y;
  array[N] int n;
  
  int K_reg;
  int K_d4;
  
  array[N] int reg;
  array[N] int d4;
  
  vector[2] pri_b_0;
  vector[2] pri_b_reg; 
  vector[2] pri_b_d4;
  
  int prior_only;
}
transformed data{
}
parameters{
  real b_0;
  vector[K_reg-1] b_reg_raw;
  vector[K_d4-1] b_d4_raw;
}
transformed parameters{
  vector[K_reg] b_reg;
  vector[K_d4] b_d4;
  
  b_reg[1] = 0.0;
  b_reg[2:K_reg] = b_reg_raw;
  b_d4[1] = 0.0;
  b_d4[2:K_d4] = b_d4_raw;
} 
model{
  target += logistic_lpdf(b_0 | pri_b_0[1], pri_b_0[2]);
  // target += student_t_lpdf(b_reg_raw | pri_b_reg[1], pri_b_reg[2], pri_b_reg[3]);
  // target += student_t_lpdf(b_d4_raw | pri_b_d4[1], pri_b_d4[2], pri_b_d4[3]);
  target += normal_lpdf(b_reg_raw | pri_b_reg[1], pri_b_reg[2]);
  target += normal_lpdf(b_d4_raw | pri_b_d4[1], pri_b_d4[2]);
  
  if(!prior_only){
    target += binomial_logit_lpmf(y | n, b_0 + b_reg[reg] + b_d4[d4]);  
  }
}
generated quantities{
}
"


mod_b <- "
data{ 
  int N;
  array[N] int y;
  array[N] int n;
  int K_reg;
  // lookup: for reg cell k, which silo column (1..4)
  array[K_reg] int reg_silo_idx;  
   // lookup: which d1 row (1=dair, 2=r1, 3=r2)
  array[K_reg] int reg_d1_idx;   
  array[K_reg] int reg_d2_idx;
  array[K_reg] int reg_d3_idx;
  
  array[N] int silo;
  
  // each dair, r1, r2
  array[N] int  d1;
  
  // nad2, wk12, wk6
  array[N] int d2;
  // nad3, none, wk12
  array[N] int d3;
  // nad4, norif, rif
  array[N] int d4;
  
  // priors
  vector[2] pri_b_0;
  vector[2] pri_b_reg; 
  vector[2] pri_b_d4;
  
  int prior_only;
}
transformed data{
}
parameters{
  real b_0;
  vector[2] b_d1_l_raw;
  vector[3] b_d1_lnr;
  vector[3] b_d1_enr;
  vector[3] b_d1_cnr;
  
  vector[2] b_d2_raw;
  vector[2] b_d3_raw;
  vector[2] b_d4_raw;
}
transformed parameters{

  // dair, r1, r2 x l, lnr, enr, cnr
  matrix[3, 4] b_d1;
  
  b_d1[1, 1] = 0.0;
  b_d1[2:3, 1] = b_d1_l_raw;
  b_d1[, 2] = b_d1_lnr;
  b_d1[, 3] = b_d1_enr;
  b_d1[, 4] = b_d1_cnr;
  
  vector[3] b_d2;
  b_d2[1] = 0.0;
  b_d2[2:3] = b_d2_raw;
  
  vector[3] b_d3;
  b_d3[1] = 0.0;
  b_d3[2:3] = b_d3_raw;
  
  vector[3] b_d4;
  b_d4[1] = 0.0;
  b_d4[2:3] = b_d4_raw;
} 
model{
  target += logistic_lpdf(b_0 | pri_b_0[1], pri_b_0[2]);
  
  target += normal_lpdf(b_d1_l_raw | 0, 2);
  target += normal_lpdf(b_d1_lnr | 0, 2);
  target += normal_lpdf(b_d1_enr | 0, 2);
  target += normal_lpdf(b_d1_cnr | 0, 2);
  
  target += normal_lpdf(b_d2_raw | 0, 2);
  target += normal_lpdf(b_d3_raw | 0, 2);
  target += normal_lpdf(b_d4_raw | 0, 2);
  
  if(!prior_only){
    for(i in 1:N){
      target += binomial_logit_lpmf(
        y[i] | n[i], b_0 + b_d1[d1[i], silo[i]] + b_d2[d2[i]] + b_d3[d3[i]] + b_d4[d4[i]]);  
    }
  }
}
generated quantities{
  vector[K_reg] b_reg;
  for (k in 1:K_reg) {
    b_reg[k] = b_d1[reg_d1_idx[k], reg_silo_idx[k]] + b_d2[reg_d2_idx[k]] + b_d3[reg_d3_idx[k]];
  }
}
"

mod_c <- "
data{ 
  int N;
  array[N] int y;
  array[N] int n;
  int K_reg;
  int K_d4;
  array[N] int reg;
  array[N] int d4;
  
  vector[2] pri_b_0;
  vector[2] pri_b_d4;
  int prior_only;
}
transformed data{
}
parameters{
  real b_0;
  real mu_reg;
  vector[K_reg-1] z_reg;
  real<lower=0> sig_reg;
  vector[K_d4-1] b_d4_raw;
}
transformed parameters{
  vector[K_reg] b_reg;
  vector[K_d4] b_d4;
  b_d4[1] = 0.0;
  b_reg[1] = 0.0;
  b_reg[2:K_reg] = mu_reg + z_reg * sig_reg;
  b_d4[2:K_d4] = b_d4_raw;
} 
model{
  target += logistic_lpdf(b_0 | pri_b_0[1], pri_b_0[2]);
  target += normal_lpdf(mu_reg | 0, 3);
  target += normal_lpdf(z_reg | 0, 1);
  target += exponential_lpdf(sig_reg | 1);
  target += normal_lpdf(b_d4_raw | pri_b_d4[1], pri_b_d4[2]);
  if(!prior_only){
    target += binomial_logit_lpmf(y | n, b_0 + b_reg[reg] + b_d4[d4]);  
  }
}
generated quantities{
}
"


m_1 <- cmdstanr::cmdstan_model(
  cmdstanr::write_stan_file(mod_a, basename = "m_1", dir = getwd()))
m_2 <- cmdstanr::cmdstan_model(
  cmdstanr::write_stan_file(mod_b, basename = "m_2", dir = getwd()))
m_3 <- cmdstanr::cmdstan_model(
  cmdstanr::write_stan_file(mod_c, basename = "m_3", dir = getwd()))
