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
logger::log_info("*** START UP - Sim07 ***")


# Command line arguments list scenario (true dose response),
# the number of simulations to run, the number of cores
# to use and the simulator to use.
args = commandArgs(trailingOnly=TRUE)

# Load cfg based on cmd line args.
if (length(args)<1) {
  log_info("Setting default run method (does nothing)")
  args[1] = "sim07_run_none"
  args[2] = "./sim07/cfg-sim07-sc01-v01.yml"
} else {
  log_info("Run method ", args[1])
  log_info("Scenario config ", args[2])
}

# stan models----------------
mod_a <- "
// group all your data up first
data{ 
  
  int N;
  
  array[N] int y;
  array[N] int n;
  
  // Parameter indexes
  array[N] int s;
  array[N] int pref;
  // silo x intv - converted via K_d1 and silo to index b/w 1:9, see d1_ix
  array[N] int d1;
  // intv + non-rand 
  array[N] int d2;
  array[N] int d3;
  array[N] int d4;
  
  // Total indexes per covariate, e.g. number of silos, d1 interventions etc
  int K_s;
  int K_p;
  int K_d1;
  int K_d2;
  int K_d3;
  int K_d4;
  
  // G-comp related subsets
  
  // subset to specific parts of the sample in order to focus on 
  // the comparisons of interest, e.g. surgical is based on the 
  // late-acute group only and the effect of interest is obtained 
  // via g-computation
  int N_d1;
  array[N_d1] int d1_s;
  array[N_d1] int d1_p;
  // assignments to each treatment level within the domain 1 cohort
  array[N_d1] int d1_d1;
  array[N_d1] int d1_d2;
  array[N_d1] int d1_d3;
  array[N_d1] int d1_d4;
  array[N_d1] int n_d1;
  int N_d1_p1;
  int N_d1_p2;
  array[N_d1_p1] int ix_d1_p1;
  array[N_d1_p2] int ix_d1_p2;
  array[N_d1_p1] int n_d1_p1;
  array[N_d1_p2] int n_d1_p2;
  real prop_p1; // proportion with preference towards one-stage
  real prop_p2; // proportion with preference towards two-stage
  
  // d2 is evaluated for those receiving one-stage revision
  int N_d2;
  array[N_d2] int d2_s;
  array[N_d2] int d2_p;
  array[N_d2] int d2_d1;
  array[N_d2] int d2_d2;
  array[N_d2] int d2_d3;
  array[N_d2] int d2_d4;
  array[N_d2] int n_d2;
  
  // d3 is evaluated for those receiving two-stage revision
  int N_d3;
  array[N_d3] int d3_s;
  array[N_d3] int d3_p;
  array[N_d3] int d3_d1;
  array[N_d3] int d3_d2;
  array[N_d3] int d3_d3;
  array[N_d3] int d3_d4;
  array[N_d3] int n_d3;
  
  // d4 is evaluated for those receiving two-stage revision
  int N_d4;
  array[N_d4] int d4_s;
  array[N_d4] int d4_p;
  array[N_d4] int d4_d1;
  array[N_d4] int d4_d2;
  array[N_d4] int d4_d3;
  array[N_d4] int d4_d4;
  array[N_d4] int n_d4;
  
  // priors
  vector[2] pri_mu;
  vector[2] pri_bs;
  vector[2] pri_bp;
  vector[2] pri_b1;
  vector[2] pri_b2;
  vector[2] pri_b3;
  vector[2] pri_b4;
  
  int prior_only;
}
transformed data{
  array[N] int d1_ix;
  
  array[N_d2] int d2_d1_ix;
  array[N_d3] int d3_d1_ix;
  array[N_d4] int d4_d1_ix;
  
  
  for(i in 1:N){
    d1_ix[i] = d1[i] + (K_d1 * (s[i] - 1));  
  } 
  // for g-comp to pick up the correct surgical domain parameter
  for(i in 1:N_d2){
    // d2_d1 should be one-stage (2) for everything since we are conditioning on 
    // one-stage but this will convert to a silo specific index for one-stage
    d2_d1_ix[i] = d2_d1[i] + (K_d1 * (d2_s[i] - 1));  
  }
  for(i in 1:N_d3){
    // As above but for two-stage. 
    // d3_d1 should be two-stage (2) for everything since we are conditioning on 
    // two-stage but this will convert to a silo specific index for two-stage
    d3_d1_ix[i] = d3_d1[i] + (K_d1 * (d3_s[i] - 1));  
  }
  for(i in 1:N_d4){
    d4_d1_ix[i] = d4_d1[i] + (K_d1 * (d4_s[i] - 1));  
  }
}
parameters{
  real mu;
  vector[K_s-1] bs_raw;
  vector[K_p-1] bp_raw;
  vector[(K_d1 * K_s) - 1] bd1_raw;
  vector[K_d2-1] bd2_raw;
  vector[K_d3-1] bd3_raw;
  vector[K_d4-1] bd4_raw;
}
transformed parameters{
  vector[K_s] bs;
  vector[K_p] bp;
  // The surgical domain needs a parameter to account for the fact that the
  // non randomised comparisons of dair/rev(1)/rev(2) are correctly reflected
  // in the linear predictor for the duration domains.
  // For example, suppose that there is no effect of revision in the early silo
  // (for whatever reason) but in the late acute group (our randomised comparison
  // for the surgical domain) there is an effect of revision. 
  // We account for the non-randomised entry into the surgical domain via the
  // silo parameters as these are identical to an indicator for non-randomised
  // surgical intervention. 
  // Here we try to account for silo specific surgical domain effects (even if
  // simply due to the non-randomised nature of the comparison for early and 
  // chronic) and their influence on the duration domains.
  // Below I declare the silo by domain effects for surgical intervention.
  vector[K_d1 * K_s] bd1;
  vector[K_d2] bd2;
  vector[K_d3] bd3;
  vector[K_d4] bd4;
  
  vector[N] eta;
  
  bs[1] = 0.0;
  bp[1] = 0.0;
  bd1[1] = 0.0;
  bd2[1] = 0.0;
  bd3[1] = 0.0;
  bd4[1] = 0.0;
  
  bs[2:K_s] = bs_raw;
  bp[2:K_p] = bp_raw;
  bd1[2:(K_d1 * K_s)] = bd1_raw;
  bd2[2:K_d2] = bd2_raw;
  bd3[2:K_d3] = bd3_raw;
  bd4[2:K_d4] = bd4_raw;
  
  for(i in 1:N){
    // dair
    if(d1[i] == 1){
      // exclude both bd2 and bd3    
      eta[i] = mu + bs[s[i]] + bp[pref[i]] + bd1[d1_ix[i]] + bd4[d4[i]]; 
    } else if (d1[i] == 2){
      // rev(1)
      eta[i] = mu + bs[s[i]] + bp[pref[i]] + bd1[d1_ix[i]] + bd2[d2[i]] + bd4[d4[i]];  
    } else {
      // rev(2)
      eta[i] = mu + bs[s[i]] + bp[pref[i]] + bd1[d1_ix[i]] + bd3[d3[i]] + bd4[d4[i]];  
    }
  }
  

} 
model{
  target += logistic_lpdf(mu | pri_mu[1], pri_mu[2]);
  target += normal_lpdf(bs_raw | pri_bs[1], pri_bs[2]);
  target += normal_lpdf(bp_raw | pri_bp[1], pri_bp[2]);
  target += normal_lpdf(bd1_raw | pri_b1[1], pri_b1[2]);
  target += normal_lpdf(bd2_raw | pri_b2[1], pri_b2[2]);
  target += normal_lpdf(bd3_raw | pri_b3[1], pri_b3[2]);
  target += normal_lpdf(bd4_raw | pri_b4[1], pri_b4[2]);
  
  if(!prior_only){
    target += binomial_logit_lpmf(y | n, eta);  
  }

}
generated quantities{
  
  
  // Surgery domain (D1) comparisons of interest are revision relative to dair
  // restricted to the late acute group. 
  vector[N_d1] wgtsd1 = dirichlet_rng(to_vector(n_d1));
  vector[N_d1_p1] wgtsd1_p1 = dirichlet_rng(to_vector(n_d1_p1));
  vector[N_d1_p2] wgtsd1_p2 = dirichlet_rng(to_vector(n_d1_p2));
  
  // For dair and rev(2), d2 is undefined and omitted from model.
  // For dair and rev(1), d3 is undefined and omitted from model.
  
  // The comparison is dair (with whatever is recvd for d2 backbone abx) and
  // by defn no d3 (ext proph) vs rev.
  // Rev is decomposed into rev(1) and rev(2).
  // We take rev(1) as including d2 (non-randomised) and d3 is undef
  // We take rev(2) as including d3 (non-randomised) and d2 is undef

  // An alternative might be to consider bd2 and bd3 when given as usual.
  // However, this will induce effects in d1 even if d1 shows no effect but
  // either of the duration domains has an effect.
  
  vector[N_d1] mu_d1_1 = inv_logit(mu + bs[2] + bd1[4] + bp[d1_p] + bd4[d1_d4]);
  
  // Assignment to revision will either end up being one or two-stage.
  // those that had the preference for one get one and those that had pref for
  // two get two. Thus we use different weights (since these are subsets of
  // our late acute cohort).
  // I have set bp explicitly here but it would be ok to just use the data passed
  // in as all records should have been selected based on the required preference.
  vector[N_d1_p1] mu_d1_2 = inv_logit(mu + bs[2] + bp[1] + bd1[5] + bd2[1] + bd4[d1_d4[ix_d1_p1]]);
  vector[N_d1_p2] mu_d1_3 = inv_logit(mu + bs[2] + bp[2] + bd1[6] + bd3[1] + bd4[d1_d4[ix_d1_p2]]);
  real p_d1_1 = wgtsd1' * mu_d1_1   ;
  real p_d1_2 = wgtsd1_p1' * mu_d1_2   ;
  real p_d1_3 = wgtsd1_p2' * mu_d1_3   ;

  // our interpretation of revision:
  real p_d1_23 = (prop_p1 * p_d1_2) + (prop_p2 * p_d1_3) ;
  real rd_d1 = p_d1_23 - p_d1_1;
  real lor_d1 = log((p_d1_23 * (1-p_d1_1)) / ( (1 - p_d1_23) * p_d1_1));

  // AB duration domain (D2) comparisons of interest are 6 wks relative to 12 wks
  // for those that have one-stage revision
  vector[N_d2] wgtsd2 = dirichlet_rng(to_vector(n_d2));
  // In the data passed to stan, d1 needs to be set to one-stage revision
  // for all units. However, this is now a silo specific view to account for the
  // possibility of differential surgical effects (for whatever reason).
  // So, we index the bd1 using d2_d1_ix to pick up either the non-rand or rand
  // comparison from the surg domain.
  // Additionally, d3 (ext-proph) is undefined for all units here and so is 
  // omitted.
  // If we do not do these conditioning steps, then we would not be making
  // logical comparisons based on the design constraints, e.g. if you have
  // one stage revision, then it is logically impossible to receive ext proph
  // based on the design rules.
  vector[N_d2] mu_d2_2 = inv_logit(mu + bs[d2_s] + bp[d2_p] + bd1[d2_d1_ix] + bd2[2] + bd4[d2_d4] );
  vector[N_d2] mu_d2_3 = inv_logit(mu + bs[d2_s] + bp[d2_p] + bd1[d2_d1_ix] + bd2[3] + bd4[d2_d4] );
  real p_d2_2 = wgtsd2' * mu_d2_2 ;
  real p_d2_3 = wgtsd2' * mu_d2_3 ;

  real rd_d2 = p_d2_3 - p_d2_2;
  real lor_d2 = log((p_d2_3 * (1-p_d2_2)) / ( (1 - p_d2_3) * p_d2_2));

  // Ext prophylaxis domain (D3) comparisons of interest are 12 wks relative to none
  // for those that have two-stage revision
  vector[N_d3] wgtsd3 = dirichlet_rng(to_vector(n_d3));
  // In the data passed to stan, d1 needs to be set to two-stage revision
  // for all units but now this is a silo specific view so we index the bd1
  // parameter to pick up either the rand or non-rand comparison from the surg
  // domain.
  // Similarly, d2 is now underfined so omitted.
  vector[N_d3] mu_d3_2 = inv_logit(mu + bs[d3_s] + bp[d3_p] + bd1[d3_d1_ix] + bd3[2] + bd4[d3_d4]);
  vector[N_d3] mu_d3_3 = inv_logit(mu + bs[d3_s] + bp[d3_p] + bd1[d3_d1_ix] + bd3[3] + bd4[d3_d4]);
  real p_d3_2 = wgtsd3' * mu_d3_2 ;
  real p_d3_3 = wgtsd3' * mu_d3_3 ;

  real rd_d3 = p_d3_3 - p_d3_2;
  real lor_d3 = log((p_d3_3 * (1-p_d3_2)) / ( (1 - p_d3_3) * p_d3_2));

  // Ab choice domain (D4) comparisons of interest are rif to none
  vector[N_d4] wgtsd4 = dirichlet_rng(to_vector(n_d4));
  vector[N_d4] mu_d4_2;
  vector[N_d4] mu_d4_3;

  for(i in 1:N_d4){
    
    if(d4_d1[i] == 1){
      // exclude both bd2 and bd3
      mu_d4_2[i] = inv_logit(mu + bs[d4_s[i]] + bp[d4_p[i]] + bd1[d4_d1_ix[i]] + bd4[2]);
      mu_d4_3[i] = inv_logit(mu + bs[d4_s[i]] + bp[d4_p[i]] + bd1[d4_d1_ix[i]] + bd4[3]);
      
    } else if (d4_d1[i] == 2){
      // exclude the bd3 term
      mu_d4_2[i] = inv_logit(mu + bs[d4_s[i]] + bp[d4_p[i]] + bd1[d4_d1_ix[i]] + bd2[d4_d2[i]] + bd4[2]);
      mu_d4_3[i] = inv_logit(mu + bs[d4_s[i]] + bp[d4_p[i]] + bd1[d4_d1_ix[i]] + bd2[d4_d2[i]] + bd4[3]);
      
    } else {  // d4_d1 is 3
      // exclude the bd2 term
      mu_d4_2[i] = inv_logit(mu + bs[d4_s[i]] + bp[d4_p[i]] + bd1[d4_d1_ix[i]] + bd3[d4_d3[i]] + bd4[2]);
      mu_d4_3[i] = inv_logit(mu + bs[d4_s[i]] + bp[d4_p[i]] + bd1[d4_d1_ix[i]] + bd3[d4_d3[i]] + bd4[3]);
      
    }
    
  }
  real p_d4_2 = wgtsd4' * mu_d4_2   ;
  real p_d4_3 = wgtsd4' * mu_d4_3   ;

  real rd_d4 = p_d4_3 - p_d4_2;
  real lor_d4 = log((p_d4_3 * (1-p_d4_2)) / ( (1 - p_d4_3) * p_d4_2));
  
  
}
"



m_1 <- cmdstanr::cmdstan_model(cmdstanr::write_stan_file(mod_a))

l_mod <- list()
l_mod[["m_1"]] <- m_1

# Data generation script takes over from the roadmap.data pkg as of 2024-12-11.

# late silo has been put last here
g_silo = c("e", "l", "c") 
g_pr_silo <- c(0.3, 0.5, 0.2)
names(g_pr_silo) <- g_silo

# matrix of the distribution of implicated joint by silo
g_pr_jnt <- matrix(
  c(0.4, 0.6, 0.7, 0.3, 0.5, 0.5), 3, 2, byrow = T
)
colnames(g_pr_jnt) <- c("knee", "hip")
rownames(g_pr_jnt) <- g_silo

# early silo randomisation probabilities to dair/rev
g_pr_e_surg <- c(0.85, 0.15)
names(g_pr_e_surg) <- c("dair", "rev")
# early silo preference to surgical type
g_pr_e_pref <- rbind(
  # for those that receive dair, the preferences are
  dair = c(0.85, 0.1, 0.05),
  # for those that receive rev, the preferences are spread across the 
  # revision options (one-stage or two-stage)
  rev = c(0, 2/3, 1/3)
)

# late silo randomisation probabilities to dair/rev
g_pr_l_surg <- c(0.5, 0.5)
names(g_pr_l_surg) <- c("dair", "rev")
# late silo preference to surgical type
g_pr_l_pref <- rbind(
  # for those that receive dair, the preferences are
  dair = c(0.2, 0.24, 0.56),
  # for those that receive rev, the preferences are spread across the 
  # revision options (one-stage or two-stage)
  rev = c(0, 0.3, 0.7)
)

# as above but for chronic silo
g_pr_c_surg <- c(0.2, 0.8)
names(g_pr_c_surg) <- c("dair", "rev")
g_pr_c_pref <- rbind(
  dair = c(0.2, 0.2, 0.6),
  rev = c(0, 0.25, 0.75)
)

output_dir_mcmc <- here::here("tmp")


# Trial implementation ---------
run_trial <- function(
    l_spec
){
  
  log_info("Entered  run_trial for trial ", l_spec$ix_sim)
  
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
  
  if(l_spec$return_posterior){
    d_post_all <- data.table()
  }
  
  
  ## LOOP -------
  while(!stop_enrol){
    
    log_info("Trial ", l_spec$ix_sim, " analysis ", l_spec$ia)
    
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
    
    log_info("Trial ", l_spec$ix_sim, " new data generated ", l_spec$ia)
    
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
      "-sim-", l_spec$ix_sim, "-intrm-", l_spec$ia)
    
    
    
    # fit model - does it matter that I continue to fit the model after the
    # decision is made...?
    # snk <- capture.output(
    f_1 <- l_mod[[l_spec$model]]$sample(
      lsd$ld, iter_warmup = 1000, iter_sampling = 2000,
      parallel_chains = 1, chains = 1, refresh = 0, show_exceptions = F,
      max_treedepth = 11,
      output_dir = output_dir_mcmc,
      output_basename = foutname
    )
    # )
    
    log_info("Trial ", l_spec$ix_sim, " fitted models ", l_spec$ia)
    
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
    
    if(l_spec$return_posterior){
      d_post_all <- rbind(
        d_post_all,
        cbind(ia = l_spec$ia, d_post)
      )
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
    
    log_info("Trial ", l_spec$ix_sim, " extracted posterior ", l_spec$ia)
    
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
    
    log_info("Trial ", l_spec$ix_sim, " calculated decision quantities ", l_spec$ia)
    
    decision[l_spec$ia, , "sup"] <- pr_sup[l_spec$ia, ] > l_spec$thresh$sup
    decision[l_spec$ia, , "ni"] <- pr_ni[l_spec$ia, ] > l_spec$thresh$ni
    
    # futility rules
    # taken to imply negligible chance of being superior or ni
    
    # negligible chance of being superior
    decision[l_spec$ia, , "fut_sup"] <- pr_sup_fut[l_spec$ia, ] < l_spec$thresh$fut_sup
    # taken to imply negligible chance of being ni 
    decision[l_spec$ia, , "fut_ni"] <- pr_ni_fut[l_spec$ia, ] < l_spec$thresh$fut_ni
    
    log_info("Trial ", l_spec$ix_sim, " compared to thresholds ", l_spec$ia)
    
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
      log_info("Stop trial all questions addressed ", l_spec$ix_sim)
      stop_enrol <- T  
    } 
    
    log_info("Trial ", l_spec$ix_sim, " updated allocation control ", l_spec$ia)
    
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
  
  if(l_spec$return_posterior){
    l_ret$d_post_all <- copy(d_post_all)
  }
  # 
  
  
  return(l_ret)
}


# Trial data ---------
get_sim07_trial_data <- function(
    l_spec,
    dec_sup = list(surg = NA, ext_proph = NA, ab_choice = NA),
    dec_ni = list(ab_dur = NA ),
    dec_sup_fut = list(surg = NA, ext_proph = NA, ab_choice = NA),
    dec_ni_fut = list(ab_dur = NA)
    ){
  
  if(is.null(l_spec$ia)){
    ia <- 1
  } else {
    ia <- l_spec$ia
  }
  if(is.null(l_spec$is)){
    is <- 1
  } else {
    is <- l_spec$is
  }
  if(is.null(l_spec$ie)){
    ie <- 1
  } else {
    ie <- l_spec$ie
  }
  if(is.null(l_spec$t0)){
    t0 <- 1
  } else {
    t0 <- l_spec$t0
  }
  
  d <- data.table(
    ia = ia,
    id = is:ie,
    t0 = t0,
    s = sample(1:3, l_spec$N[ia], replace = T, prob = l_spec$p_s_alloc)
  )
  d[s == 1, `:=`(
    # ctl/trt allocations ignoring dependencies
    d1_alloc = rbinom(.N, 1, l_spec$l_e$p_d1_alloc),
    # 70% entery d2
    d2_entry = rbinom(.N, 1, l_spec$l_e$p_d2_entry),
    d2_alloc = rbinom(.N, 1, l_spec$l_e$p_d2_alloc),
    # 90% enter d3
    d3_entry = rbinom(.N, 1, l_spec$l_e$p_d3_entry),
    d3_alloc = rbinom(.N, 1, l_spec$l_e$p_d3_alloc),
    # 60% enter d4
    d4_entry = rbinom(.N, 1, l_spec$l_e$p_d4_entry),
    d4_alloc = rbinom(.N, 1, l_spec$l_e$p_d4_alloc),
    # preference directs type of revision (0 rev(1), 1 rev(2))
    pref = rbinom(.N, 1, l_spec$l_e$p_pref)  + 1 
  )]
  d[s == 2, `:=`(
    # ctl/trt allocations ignoring dependencies
    d1_alloc = rbinom(.N, 1, l_spec$l_l$p_d1_alloc),
    # 70% entery d2
    d2_entry = rbinom(.N, 1, l_spec$l_l$p_d2_entry),
    d2_alloc = rbinom(.N, 1, l_spec$l_l$p_d2_alloc),
    # 90% enter d3
    d3_entry = rbinom(.N, 1, l_spec$l_l$p_d3_entry),
    d3_alloc = rbinom(.N, 1, l_spec$l_l$p_d3_alloc),
    # 60% enter d4
    d4_entry = rbinom(.N, 1, l_spec$l_l$p_d4_entry),
    d4_alloc = rbinom(.N, 1, l_spec$l_l$p_d4_alloc),
    # preference directs type of revision (0 rev(1), 1 rev(2))
    pref = rbinom(.N, 1, l_spec$l_l$p_pref)  + 1 
  )]
  d[s == 3, `:=`(
    # ctl/trt allocations ignoring dependencies
    d1_alloc = rbinom(.N, 1, l_spec$l_c$p_d1_alloc),
    # 70% entery d2
    d2_entry = rbinom(.N, 1, l_spec$l_c$p_d2_entry),
    d2_alloc = rbinom(.N, 1, l_spec$l_c$p_d2_alloc),
    # 90% enter d3
    d3_entry = rbinom(.N, 1, l_spec$l_c$p_d3_entry),
    d3_alloc = rbinom(.N, 1, l_spec$l_c$p_d3_alloc),
    # 60% enter d4
    d4_entry = rbinom(.N, 1, l_spec$l_c$p_d4_entry),
    d4_alloc = rbinom(.N, 1, l_spec$l_c$p_d4_alloc),
    # preference directs type of revision (0 rev(1), 1 rev(2))
    pref = rbinom(.N, 1, l_spec$l_c$p_pref)  + 1 
  )]
  
  # surgical -----
  # assume everyone enters surgical
  if(is.na(dec_sup$surg) & is.na(dec_sup_fut$surg)){
    # dair gets dair, revision gets split
    d[d1_alloc == 0, d1 := 1]
    d[d1_alloc == 1 & pref == 1, d1 := 2]
    d[d1_alloc == 1 & pref == 2, d1 := 3]
    
  } else if (!is.na(dec_sup$surg)) {
    # Revision has been deemed superior - allocation is now either one-stage or
    # two-stage dependent on preference but only for the late silo
  
    # early and chronic as they were
    d[s != 2 & d1_alloc == 0, d1 := 1]
    d[s != 2 & d1_alloc == 1 & pref == 1, d1 := 2]
    d[s != 2 & d1_alloc == 1 & pref == 2, d1 := 3]
    
    # late now onto rev
    d[s == 2 & pref == 1, d1 := 2]
    d[s == 2 & pref == 2, d1 := 3]
    
  } else if (!is.na(dec_sup_fut$surg)) {
    # Revert to dair as the best option (superiority assessment is futile)
    
    # early and chronic as they were
    d[s != 2 & d1_alloc == 0, d1 := 1]
    d[s != 2 & d1_alloc == 1 & pref == 1, d1 := 2]
    d[s != 2 & d1_alloc == 1 & pref == 2, d1 := 3]
    
    # late now onto dair
    d[s == 2, d1 := 1]
    d[s == 2, d1 := 1]
  }
  
  # abx dur -----
  
  # undefined unless under rev(1)
  d[d1 %in% c(1, 3), d2 := NA_integer_]
  if(is.na(dec_ni$ab_dur) & is.na(dec_ni_fut$ab_dur) ){
    # default
    d[d1 == 2 & d2_entry == 0, d2 := 1]
    d[d1 == 2 & d2_entry == 1, d2 := 2 + d2_alloc]
    
  } else if (!is.na(dec_ni$ab_dur) ) {   
    # 6wks is NI to 12wks (or equivalent) so we assume that those receiving one-stage revision
    # all receive the 6wk trt
    d[d1 %in% c(2) & d2_entry == 0, d2 := 1]
    d[d1 %in% c(2) & d2_entry == 1, d2 := 3]
    
  } else if (!is.na(dec_ni_fut$ab_dur) ) {
    # 6wk ni assessment is futile or inferior, everyone now gets 12wk
    d[d1 %in% c(2) & d2_entry == 0, d2 := 1]
    d[d1 %in% c(2) & d2_entry == 1, d2 := 2]
  }
  
  # ext proph -----
  # undefined unless under rev(2)
  d[d1 %in% c(1, 2), d3 := NA_integer_]
  if(is.na(dec_sup$ext_proph) & is.na(dec_sup_fut$ext_proph)){
    # Default situation, we are allocating randomised trt
    d[d1 == 3 & d3_entry == 0, d3 := 1]
    d[d1 == 3 & d3_entry == 1, d3 := 2 + d3_alloc]
    
  } else if (!is.na(dec_sup$ext_proph)){
    # 12wks is superior to none so we assume that those receiving two-stage revision
    # all receive the 12wk trt
    d[d1 == 3 & d3_entry == 0, d3 := 1]
    d[d1 == 3 & d3_entry == 1, d3 := 3]
    
  } else if (!is.na(dec_sup_fut$ext_proph)){
    # 12wks is futile, everyone now gets none
    d[d1 == 3 & d3_entry == 0, d3 := 1]
    d[d1 == 3 & d3_entry == 1, d3 := 2]
    
  }
  
  ## Antibiotic choice ----
  
  # Choice domain is independent to others but only 60% of the population
  # enter it.
  if(is.na(dec_sup$ab_choice) & is.na(dec_sup_fut$ab_choice)){
    # Default situation
    d[d4_entry == 0, d4 := 1]
    d[d4_entry == 1, d4 := 2 + d4_alloc]
    
  } else if (!is.na(dec_sup$ab_choice)){
    # Superiority decision, all get allocated to rif
    d[d4_entry == 0, d4 := 1]
    d[d4_entry == 1, d4 := 3]
    
  } else if (!is.na(dec_sup_fut$ab_choice)){
    # rif is futile, everyone now gets none
    d[d4_entry == 0, d4 := 1]
    d[d4_entry == 1, d4 := 2]
    
  }
  
  # compute linear predictor
  
  bd1 <- c(l_spec$l_e$bd1, l_spec$l_l$bd1, l_spec$l_c$bd1)
  
  # bd1 <- c(0, rnorm(8))
  # index for d1 is a function of 
  # d1 (3 levels - dair, rev(1), rev(2)) and silo membership (1:3)
  d[, d1_ix := d1 + (3 * (s - 1))]

  # bd1 is irrelevant since d1 == 1 is the ref group, fixed at zero
  d[d1 == 1, eta := l_spec$mu + l_spec$bs[s] + l_spec$bp[pref] + bd1[d1_ix] + l_spec$bd4[d4]]
  # pref is irrelevant as d1 = 2 only occurs if pref = 0
  d[d1 == 2, eta := l_spec$mu + l_spec$bs[s] + bd1[d1_ix] + l_spec$bd2[d2] + l_spec$bd4[d4]]
  # but here pref is relevant as d1 = 3 only if pref = 1
  d[d1 == 3, eta := l_spec$mu + l_spec$bs[s] + l_spec$bp[pref] + bd1[d1_ix] + l_spec$bd3[d3] + l_spec$bd4[d4]]
  
  d[, `:=`(s = factor(s), 
           d1 = factor(d1), 
           d2 = factor(d2, levels = 1:4), 
           d3 = factor(d3, levels = 1:4), 
           d4 = factor(d4))]
  
  d[, p := plogis(eta)]
  d[, y := rbinom(.N, 1, p)]
  
  d
}






get_sim07_stan_data <- function(d_all){
  
  # convert from binary representation to binomial (successes/trials)
  d_mod <- d_all[, .(y = sum(y), n = .N, eta = round(unique(eta), 3)), 
                 keyby = .(s, pref, d1, d2, d3, d4)]
  
  
  d_mod[, `:=`(
    s = as.integer(s),
    d1 = as.integer(d1),
    d2 = as.integer(d2),
    d3 = as.integer(d3),
    d4 = as.integer(d4)
  )]
  
  d_mod[is.na(d2), d2 := 999]
  d_mod[is.na(d3), d3 := 999]
  
  K_d1 <- length(unique(d_all$d1))
  d_mod[, d1_ix := d1 + (K_d1 * (s - 1))]
  
  d_mod[, eta_obs := qlogis(y / n)]
  d_mod[, p_obs := y / n]
  
  # g-comp setup for all of the domains to provide
  # a uniform approach for generating the parameters of interest.
  
  # Surgical -----
  
  # restrict to silo 2 (late acute silo) for gcomp for d1 randomised comparisons
  # here I just use the silo assignment to imply the right subset of units
  # but in practice there may need to be an indicator variable to for this.
  d_mod_d1 <- d_mod[s == 2]
  
  
  # AB duration ----
  
  # restrict to one-stage for gcomp
  # surgery assignment being one-stage (d1 at level 2)
  # permits entry
  d_mod_d2 <- d_mod[d1 == 2]
  
  # Ext proph duration ----
  
  # restrict to two-stage for gcomp - 
  # surgery assignment being two-stage (d1 at level 3)
  # permits entry
  d_mod_d3 <- d_mod[d1 == 3]
  
  # AB choice ----
  
  # restrict to d4 randomised group - levels 2 and 3 are the randomised groups.
  # The units having d4 set to 1 were not included in ab choice
  d_mod_d4 <- d_mod[d4 %in% 2:3]
  
  ld <- list(
    # full dataset
    N = nrow(d_mod), 
    y = d_mod[, y], 
    n = d_mod[, n], 
    s = d_mod[, s], 
    pref = d_mod[, pref],
    d1 = d_mod[, d1],
    d2 = d_mod[, d2],
    d3 = d_mod[, d3],
    d4 = d_mod[, d4],
    
    # Number of levels for silos, joints, pref and each trt.
    K_s = length(unique(d_mod$s)), 
    K_p = length(unique(d_mod$pref)), 
    K_d1 = d_mod[, length(unique(d1))], 
    K_d2 = d_mod[d2 != 999, length(unique(d2))], 
    K_d3 = d_mod[d3 != 999, length(unique(d3))], 
    K_d4 = d_mod[, length(unique(d4))], 
    
    
    # cohort for surgical domain g-comp complicated due to one-stage/two-stage
    # considerations
    N_d1 = nrow(d_mod_d1),
    d1_s = d_mod_d1[, s], 
    d1_p = d_mod_d1[, pref],
    # all d1 assignments
    d1_d1 = d_mod_d1[, d1], 
    d1_d2 = d_mod_d1[, d2],
    d1_d3 = d_mod_d1[, d3],
    d1_d4 = d_mod_d1[, d4],
    # number of trials within each strata
    n_d1 = d_mod_d1[, n],
    # sample size of those where preference is for one-stage
    N_d1_p1 = d_mod_d1[pref == 1, .N],
    # sample size of those where preference is for two-stage
    N_d1_p2 = d_mod_d1[pref == 2, .N],
    # indexes for those with preference for one-stage
    ix_d1_p1 = d_mod_d1[pref == 1, which = T],
    # indexes for those with preference for two-stage
    ix_d1_p2 = d_mod_d1[pref == 2, which = T],
    # number of trials within each of these subsets
    n_d1_p1 = d_mod_d1[pref == 1, n],
    n_d1_p2 = d_mod_d1[pref == 2, n],
    
    prop_p1 = d_all[s == 2, .(wgt = .N/nrow(d_all[s == 2])), keyby = pref][pref == 1, wgt],
    prop_p2 = d_all[s == 2, .(wgt = .N/nrow(d_all[s == 2])), keyby = pref][pref == 2, wgt],
    
    N_d2 = nrow(d_mod_d2),
    d2_s = d_mod_d2[, s], 
    d2_p = d_mod_d2[, pref],
    d2_d1 = d_mod_d2[, d1], 
    d2_d2 = d_mod_d2[, d2],
    d2_d3 = d_mod_d2[, d3],
    d2_d4 = d_mod_d2[, d4],
    n_d2 = d_mod_d2[, n],
    
    N_d3 = nrow(d_mod_d3),
    d3_s = d_mod_d3[, s], 
    d3_p = d_mod_d3[, pref],
    d3_d1 = d_mod_d3[, d1], 
    d3_d2 = d_mod_d3[, d2],
    d3_d3 = d_mod_d3[, d3],
    d3_d4 = d_mod_d3[, d4],
    n_d3 = d_mod_d3[, n],
    
    N_d4 = nrow(d_mod_d4),
    d4_s = d_mod_d4[, s], 
    d4_p = d_mod_d4[, pref],
    d4_d1 = d_mod_d4[, d1], 
    d4_d2 = d_mod_d4[, d2],
    d4_d3 = d_mod_d4[, d3],
    d4_d4 = d_mod_d4[, d4],
    n_d4 = d_mod_d4[, n],
    
    prior_only = 0
  )
  
  list(
    d_mod = d_mod,
    d_mod_d1 = d_mod_d1,
    d_mod_d2 = d_mod_d2,
    d_mod_d3 = d_mod_d3,
    d_mod_d4 = d_mod_d4,
    
    ld = ld
  )
  
}


get_sim07_trial_data_subgrp <- function(
    l_spec
){
  
  
  # just demo for explanation purposes
  
  d <- get_sim07_trial_data(l_spec)
  
  # d1 -------
  # subgroups fixed to those in protocol
  # arbitrary distribution for each subgroup
  # site of infection (knee, hip) 
  # pr(site = knee| e) = 0.4
  # pr(site = knee| l) = 0.7
  # pr(site = knee| c) = 0.5
  d[s == 1, d1_sg_1 := sample(1:2, size = .N, replace = T, prob = c(0.4, 0.6))]
  d[s == 2, d1_sg_1 := sample(1:2, size = .N, replace = T, prob = c(0.7, 0.3))]
  d[s == 3, d1_sg_1 := sample(1:2, size = .N, replace = T, prob = c(0.5, 0.5))]
  
  # duration of symptom at entry (3 category)
  d[, d1_sg_2 := sample(1:3, size = .N, replace = T, prob = c(0.25, 0.5, 0.25))]

  # one of causative org is s aureus Pr(s aureus) = 0.3
  d[, d1_sg_3 := sample(1:2, size = .N, replace = T, prob = c(0.7, 0.3))]

  # crp at entry > 100 Pr(crp > 100) = 0.2
  d[, d1_sg_4 := sample(1:2, size = .N, replace = T, prob = c(0.8, 0.2))]

  
  # d2 -------
  # silo (already present in model)
  d[, d2_sg_1 := copy(s)]
  # one of causative org is s aureus Pr(s aureus) = 0.3 (same as d1)
  d[, d2_sg_2 := copy(d1_sg_3)]
  # revision all ideal
  d[, d2_sg_3 := sample(1:2, size = .N, replace = T, prob = c(0.5, 0.5))]
  
  
  # d3 ------
  # silo (already present in model)
  d[, d3_sg_1 := copy(s)]
  # one of causative org is s aureus Pr(s aureus) = 0.3 (same as d1)
  d[, d3_sg_2 := copy(d1_sg_3)]
  # Categorised duration between first-stage and reimplantation procedure
  d[, d3_sg_3 := sample(1:3, size = .N, replace = T, prob = c(0.3, 0.4, 0.3))]
  
  # d4 ------
  # type of surgery already present
  d[, d4_sg_1 := copy(d1)]
  # one of causative org is s aureus Pr(s aureus) = 0.3 (same as d1)
  d[, d4_sg_2 := copy(d1_sg_3)]
  
  
  # revise linear predictor ------
  
  # main effects for bd1 as before
  bd1_star <- c(l_spec$l_e$bd1, l_spec$l_l$bd1, l_spec$l_c$bd1)
  
  # d1 subgroups
  # silo (already dealt with as there are silo specific effects)
  # now deal with:
  # site of infection (2 level)
  # duration symptom (assume 3 level), 
  # staph (2 level), 
  # crp at baseline (2 level), 
  # duration index implant (assume 3 level), 
  # rev type (3 level main effects already present)
  
  bd1_sg_sig <- 0.3
  # 3 x 2, (dair, rev(1), rev(2)) x (knee, hip) possibilities for pr(trt success)
  bd1_sg_1 <- matrix(rnorm(6, 0, bd1_sg_sig), ncol = 2, nrow = 3)
  # 3x3, (dair, rev(1), rev(2)) x (short, med, long duration of symptoms) 
  # bd1_sg_2 <- matrix(rnorm(9, 0, bd1_sg_sig), ncol = 3, nrow = 3)
  # # stapha 
  # bd1_sg_3 <- matrix(rnorm(6, 0, bd1_sg_sig), ncol = 2, nrow = 3)
  # # crp level
  # bd1_sg_4 <- matrix(rnorm(6, 0, bd1_sg_sig), ncol = 2, nrow = 3)
  
  # d2 subgroups
  # bd2_sg_sig <- 0.2
  # bd2_sg_1 <- matrix(rnorm(9, 0, bd2_sg_sig), ncol = 3, nrow = 3)
  # bd2_sg_2 <- matrix(rnorm(6, 0, bd2_sg_sig), ncol = 2, nrow = 3)
  # bd2_sg_3 <- matrix(rnorm(9, 0, bd2_sg_sig), ncol = 3, nrow = 3)
  
  # d3 subgroups
  # bd3_sg_sig <- 0.24
  # bd3_sg_1 <- matrix(rnorm(9, 0, bd3_sg_sig), ncol = 3, nrow = 3)
  # bd3_sg_2 <- matrix(rnorm(6, 0, bd3_sg_sig), ncol = 2, nrow = 3)
  # bd3_sg_3 <- matrix(rnorm(9, 0, bd3_sg_sig), ncol = 3, nrow = 3)
  
  # d4 subgroups
  # bd4_sg_sig <- 0.1
  # # type of surgery
  # bd4_sg_1 <- matrix(rnorm(9, 0, bd4_sg_sig), ncol = 3, nrow = 3)
  # # at least one is staph
  # bd4_sg_2 <- matrix(rnorm(6, 0, bd4_sg_sig), ncol = 2, nrow = 3)
  
  d[, bd1_star := bd1_star[d1_ix] ]
  d[, bd1_sg_1 := bd1_sg_1[cbind(d1, d1_sg_1)]]
  d[, bd1 := bd1_star +  bd1_sg_1]
  # + bd1_sg_2[cbind(d1, d1_sg_2)] + bd1_sg_3[cbind(d1, d1_sg_3)] + bd1_sg_4[cbind(d1, d1_sg_4)] ]
  # d[, bd2 := l_spec$bd2[d2] + bd2_sg_1[cbind(d2, d2_sg_1)] + bd2_sg_2[cbind(d2, d2_sg_2)] + bd2_sg_3[cbind(d2, d2_sg_3)] ]
  # d[, bd3 := l_spec$bd3[d3] + bd3_sg_1[cbind(d3, d3_sg_1) ] + bd3_sg_2[cbind(d3, d3_sg_2)] + bd3_sg_3[cbind(d3, d3_sg_3)] ]
  # d[, bd4 := l_spec$bd4[d4] + bd4_sg_1[cbind(d4, d4_sg_1) ] + bd4_sg_2[cbind(d4, d4_sg_2)] ]
  
  # bd1 is irrelevant since d1 == 1 is the ref group, fixed at zero
  d[, eta := l_spec$mu + l_spec$bs[s] + l_spec$bp[pref] + bd1]
  
  # bd1 is irrelevant since d1 == 1 is the ref group, fixed at zero
  # d[d1 == 1, eta := l_spec$mu + l_spec$bs[s] + l_spec$bp[pref] + bd1 + bd4]
  # # pref is irrelevant as d1 = 2 only occurs if pref = 0
  # d[d1 == 2, eta := l_spec$mu + l_spec$bs[s] + bd1 + bd2 + bd4]
  # # but here pref is relevant as d1 = 3 only if pref = 1
  # d[d1 == 3, eta := l_spec$mu + l_spec$bs[s] + l_spec$bp[pref] + bd1 + bd3 + bd4]
  
  # bd1 is irrelevant since d1 == 1 is the ref group, fixed at zero
  # d[d1 == 1, eta := l_spec$mu + l_spec$bs[s] + l_spec$bp[pref] + bd1 + l_spec$bd4[d4]]
  # # pref is irrelevant as d1 = 2 only occurs if pref = 0
  # d[d1 == 2, eta := l_spec$mu + l_spec$bs[s] + bd1 + l_spec$bd2[d2] + l_spec$bd4[d4]]
  # # but here pref is relevant as d1 = 3 only if pref = 1
  # d[d1 == 3, eta := l_spec$mu + l_spec$bs[s] + l_spec$bp[pref] + bd1 + l_spec$bd3[d3] + l_spec$bd4[d4]]
  
  d[, `:=`(s = factor(s), 
           d1 = factor(d1), 
           d2 = factor(d2, levels = 1:4), 
           d3 = factor(d3, levels = 1:4), 
           d4 = factor(d4))]
  
  d[, p := plogis(eta)]
  d[, y := rbinom(.N, 1, p)]
  
  d_tmp <- unique(d[, .(s, pref, d1, d1_sg_1, bd1_star, bd1_sg_1, eta)])
  setkey(d_tmp, s, pref, d1, d1_sg_1)
  d_tmp[]
  d[]
  
}


get_sim07_stan_data_subgrp <- function(d_all){
  
  # convert from binary representation to binomial (successes/trials)
  d_mod <- d_all[, .(y = sum(y), n = .N, p_tru = round(unique(p), 3)), 
                 keyby = .(
                   s, pref, 
                   d1, d1_sg_1)]
  
  d_mod[, `:=`(
    s = as.integer(s),
    d1 = as.integer(d1),
    d1_sg_1 = as.integer(d1_sg_1)
    # ,
    # d1_sg_2 = as.integer(d1_sg_2),
    # d1_sg_3 = as.integer(d1_sg_3),
    # d1_sg_4 = as.integer(d1_sg_4)
  )]
  
  d_mod[, eta_obs := qlogis(y / n)]
  d_mod[, p_obs := y / n]
  
  K_d1 <- length(unique(d_all$d1))
  d_mod[, d1_ix := d1 + (K_d1 * (s - 1))]
  
  
  
  ld <- list(
    # full dataset
    N = nrow(d_mod), 
    y = d_mod[, y], 
    n = d_mod[, n], 
    s = d_mod[, s], 
    pref = d_mod[, pref],
    d1 = d_mod[, d1],
    d1_sg_1 = d_mod[, d1_sg_1],
    # d1_sg_2 = d_mod[, d1_sg_2],
    # d1_sg_3 = d_mod[, d1_sg_3],
    # d1_sg_4 = d_mod[, d1_sg_4],
    
    # Number of levels for silos, joints, pref and each trt.
    K_s = length(unique(d_mod$s)), 
    K_p = length(unique(d_mod$pref)), 
    K_d1 = d_mod[, length(unique(d1))],
    
    prior_only = 0
  )
  
  list(
    d_mod = d_mod,
    
    ld = ld
  )
  
}




get_sim07_enrol_time_int <- function(N = 2500, lambda = 1.52,
                           rho = function(t) pmin(t/360, 1)){
  
  c(0, poisson::nhpp.event.times(lambda, N - 1, rho))
}






# Util --------
# main simulation report ------
sim07_report_sim_res <- function(){
  
  library(data.table)
  library(qs2)
  library(kableExtra)
  
  get_delta <- function(l_spec, domain = "d1", rule = "sup"){
    l_spec$dec[[domain]][[rule]]$delta
  }
  get_thres <- function(l_spec, domain = "d1", rule = "sup"){
    l_spec$dec[[domain]][[rule]]$thresh
  }
  
  sim_dat_dir <-  "sim07-10"  
  
  fname <- paste0(
    sim_dat_dir, "-", format(Sys.time(), "%Y%m%d-%H%M%S"), ".md")
  f_out <- file(here::here("versions", fname), open = "w")
  
  
  f_list <- list.files(here::here("data", sim_dat_dir))
  f = f_list[1]
  for(f in f_list){
    
    message("##### ", f)
    l <- qs2::qs_read(here::here("data", sim_dat_dir, f))
    
    l_spec = l$l_spec
    # 
    # d_pr_sup = l$d_pr_sup
    # d_pr_ni = l$d_pr_ni
    # 
    # d_pr_sup_fut = l$d_pr_sup_fut
    # d_pr_ni_fut = l$d_pr_ni_fut
    # 
    # d_decision = l$d_decision
    # d_post_smry_1 = l$d_post_smry_1
    # d_post_smry_2 = l$d_post_smry_2
    # 
    # d_all = l$d_all
    # 
    # d_post_all = l$d_post_all
    
    
    writeLines(paste0("# Source data: ", f, " (", sim_dat_dir, ")"), f_out)
    
    sim07_report_report_file(l, l_spec, f_out)
    
    message("End of results for file", f)
  }
  
  message("Close file")
  close(f_out)
  
}

sim07_report_report_file <- function(
    l,
    l_spec,
    f_out
){
  
  d_pr <- sim07_decision_prob(l$d_decision)
  
  writeLines("### Decision probabilities", f_out)
  writeLines("\n", f_out)
  
  writeLines("Cumulative probability of each decision: ", f_out)
  writeLines("\n", f_out)
  tbl <- kableExtra::kbl(
    dcast(d_pr, domain + quant ~ ia, value.var = "pr_val"),
    digits = 3, format = "simple", 
    caption = paste(
      "Scenario ",
      l_spec$desc, " - Probability of decision")
  )
  writeLines(tbl, f_out)
  writeLines("\n", f_out)
}

sim07_decision_prob <- function(
    d_decision
    ){
  
  # long version of decisions
  d_dec_1 <- melt(d_decision, 
                  id.vars = c("sim", "ia", "quant"), 
                  variable.name = "domain")
  
  # Should be right, but just in case...
  if(any(is.na(d_dec_1$value))){
    message("Some of the decision values are NA in index ", i, " file ", flist[i])
    d_dec_1[is.na(value), value := FALSE]
  }
  d_dec_1[, domain := as.numeric(gsub("d", "", domain))]
    
  # Domains 1, 3 and 4 will stop for superiority or futility for superiority.
  # Domaain 2 will stop for NI or futility for NI.
  # No other stopping rules apply and so we only evaluate the operating 
  # characteristics on these, i.e. we do not care about the results for the 
  # cumualative probability of ni for domain 1, 3 and 4 because we would never
  # stop for this. 
  d_dec_1 <- rbind(
    d_dec_1[domain %in% c(1, 3, 4) & quant %in% c("sup", "fut_sup")],
    d_dec_1[domain %in% c(2) & quant %in% c("ni", "fut_ni")]
    )
    
    
  # compute the cumulative instances of a decision being made by sim, each 
  # decision type and by parameter
  d_dec_1[, value := as.logical(cumsum(value)>0), keyby = .(sim, quant, domain)]
     
  # cumulative proportion for which each decision quantity has been met by 
  # analysis and domain
  d_dec_cprob <- d_dec_1[, .(pr_val = mean(value)), keyby = .(ia, quant, domain)]
  
  d_dec_cprob
  
  
}


sim07_update_cfg <- function(l_spec){
  
  if(unname(Sys.info()[1]) == "Darwin"){
    log_info("On mac, reset cores to 5")
    l_spec$mc_cores <- 5
  } 
  
  log_info("Only model is m_1")
  l_spec$model = "m_1"
  
  l_spec$N <- l_spec$N_pt
  # These calls are just to maintain backward compatibility 
  l_spec$N_pt <- NULL
  log_info("Silo allocation is fixed")
  l_spec$p_s_alloc <- c(0.3, 0.5, 0.2)

  l_spec$mu <- l_spec$bmu
  # drop this from list only retained for backward compat
  l_spec$bmu <- NULL
  
  l_spec$bs <- unlist(l_spec$bs)
  l_spec$bp <- unlist(l_spec$bp)
  
  l_spec$l_e <- list()
  l_spec$l_l <- list()
  l_spec$l_c <- list()
  
  # dair, one, two-stage, we compare avg of one and two stage rev to dair
  l_spec$l_e$bd1 <- unlist(l_spec$bed1)
  l_spec$l_l$bd1 <- unlist(l_spec$bld1)
  l_spec$l_c$bd1 <- unlist(l_spec$bcd1)
  
  l_spec$bed1 <- NULL
  l_spec$bld1 <- NULL
  l_spec$bcd1 <- NULL
  
  l_spec$l_e$p_d1_alloc <- l_spec$e_p_d1_alloc
  l_spec$l_e$p_d2_entry <- l_spec$e_p_d2_entry
  l_spec$l_e$p_d2_alloc <- l_spec$e_p_d2_alloc
  l_spec$l_e$p_d3_entry <- l_spec$e_p_d3_entry
  l_spec$l_e$p_d3_alloc <- l_spec$e_p_d3_alloc
  l_spec$l_e$p_d4_entry <- l_spec$e_p_d4_entry
  l_spec$l_e$p_d4_alloc <- l_spec$e_p_d4_alloc
  # preference for two-stage
  l_spec$l_e$p_pref <- l_spec$e_p_pref
  
  l_spec$e_p_d1_alloc <- NULL
  l_spec$e_p_d2_entry <- NULL
  l_spec$e_p_d2_alloc <- NULL
  l_spec$e_p_d3_entry <- NULL
  l_spec$e_p_d3_alloc <- NULL
  l_spec$e_p_d4_entry <- NULL
  l_spec$e_p_d4_alloc <- NULL
  # preference for two-stage
  l_spec$e_p_pref <- NULL
  
  l_spec$l_l$p_d1_alloc <- l_spec$l_p_d1_alloc
  l_spec$l_l$p_d2_entry <- l_spec$l_p_d2_entry
  l_spec$l_l$p_d2_alloc <- l_spec$l_p_d2_alloc
  l_spec$l_l$p_d3_entry <- l_spec$l_p_d3_entry
  l_spec$l_l$p_d3_alloc <- l_spec$l_p_d3_alloc
  l_spec$l_l$p_d4_entry <- l_spec$l_p_d4_entry
  l_spec$l_l$p_d4_alloc <- l_spec$l_p_d4_alloc
  l_spec$l_l$p_pref <- l_spec$l_p_pref
  
  l_spec$l_p_d1_alloc <- NULL
  l_spec$l_p_d2_entry <- NULL
  l_spec$l_p_d2_alloc <- NULL
  l_spec$l_p_d3_entry <- NULL
  l_spec$l_p_d3_alloc <- NULL
  l_spec$l_p_d4_entry <- NULL
  l_spec$l_p_d4_alloc <- NULL
  l_spec$l_p_pref <- NULL
  
  l_spec$l_c$p_d1_alloc <- l_spec$c_p_d1_alloc
  l_spec$l_c$p_d2_entry <- l_spec$c_p_d2_entry
  l_spec$l_c$p_d2_alloc <- l_spec$c_p_d2_alloc
  l_spec$l_c$p_d3_entry <- l_spec$c_p_d3_entry
  l_spec$l_c$p_d3_alloc <- l_spec$c_p_d3_alloc
  l_spec$l_c$p_d4_entry <- l_spec$c_p_d4_entry
  l_spec$l_c$p_d4_alloc <- l_spec$c_p_d4_alloc
  # preference for two-stage
  l_spec$l_c$p_pref <- l_spec$c_p_pref
  
  l_spec$c_p_d1_alloc <- NULL
  l_spec$c_p_d2_entry <- NULL
  l_spec$c_p_d2_alloc <- NULL
  l_spec$c_p_d3_entry <- NULL
  l_spec$c_p_d3_alloc <- NULL
  l_spec$c_p_d4_entry <- NULL
  l_spec$c_p_d4_alloc <- NULL
  l_spec$c_p_pref <- NULL
  
  # always ref, 12wk, 6wk as we are assessing if 6wk ni to 12wk
  l_spec$bd2 <- unlist(l_spec$bd2)
  # always ref, 0, 12wk as we are assessing if 12wk sup to none
  l_spec$bd3 <- unlist(l_spec$bd3)
  # always ref, none, rif as we are assessing if rif is sup to none
  l_spec$bd4 <- unlist(l_spec$bd4)
  
  l_spec$prior <- list()
  # location, scale
  l_spec$prior$mu <- unlist(l_spec$pri_bmu)
  l_spec$prior$bs <- unlist(l_spec$pri_bs)
  l_spec$prior$bp <- unlist(l_spec$pri_bp)
  l_spec$prior$bd1 <- unlist(l_spec$pri_bd1)
  l_spec$prior$bd2 <- unlist(l_spec$pri_bd2)
  l_spec$prior$bd3 <- unlist(l_spec$pri_bd3)
  l_spec$prior$bd4 <- unlist(l_spec$pri_bd4)
  
  l_spec$pri_bmu <- NULL
  l_spec$pri_bs <- NULL
  l_spec$pri_bp <- NULL
  l_spec$pri_bd1 <- NULL
  l_spec$pri_bd2 <- NULL
  l_spec$pri_bd3 <- NULL
  l_spec$pri_bd4 <- NULL
  
  l_spec$delta <- list()
  l_spec$delta$sup <- l_spec$dec_delta_sup
  l_spec$delta$sup_fut <- l_spec$dec_delta_sup_fut
  l_spec$delta$ni <- l_spec$dec_delta_ni
  l_spec$delta$ni_fut <- l_spec$dec_delta_ni_fut
  
  l_spec$dec_delta_sup <- NULL
  l_spec$dec_delta_sup_fut <- NULL
  l_spec$dec_delta_ni <- NULL
  l_spec$dec_delta_ni_fut <- NULL
  
  # domain specific
  l_spec$thresh <- list()
  l_spec$thresh$sup <- unlist(l_spec$dec_thresh_sup)
  l_spec$thresh$ni <- unlist(l_spec$dec_thresh_fut_sup)
  l_spec$thresh$fut_sup <- unlist(l_spec$dec_thresh_ni)
  l_spec$thresh$fut_ni <- unlist(l_spec$dec_thresh_fut_ni)
  
  l_spec$dec_thresh_sup <- NULL
  l_spec$dec_thresh_fut_sup <- NULL
  l_spec$dec_thresh_ni <- NULL
  l_spec$dec_thresh_fut_ni <- NULL
  
  log_info("Starting simulation with following parameters:");
  log_info("N: ", paste0(l_spec$N, collapse = ", "));
  log_info("b_silo: ", paste0(l_spec$bs, collapse = ", "));
  log_info("b_pref: ", paste0(l_spec$bp, collapse = ", "));
  log_info("b_d1: ", paste0(c(l_spec$l_e$bd1, l_spec$l_l$bd1, l_spec$l_c$bd1), collapse = ", "));
  log_info("b_d2: ", paste0(l_spec$bd2, collapse = ", "));
  log_info("b_d3: ", paste0(l_spec$bd3, collapse = ", "));
  log_info("b_d4: ", paste0(l_spec$bd4, collapse = ", "));
  
  if(l_spec$nex > 0){
    l_spec$nex <- pmin(l_spec$nex, l_spec$n_sim)
    l_spec$ex_trial_ix <- sort(
      sample(1:l_spec$n_sim, size = l_spec$nex, replace = F))
    l_spec$ex_trial_ix[1] <- 1
  }
  
  l_spec
  
}


sim07_default_cfg <- function(n_sim = 5){
  
  l_spec <- list(
    desc = "RD = 0 in all domains +silo specific d1",
    n_sim = n_sim,
    mc_cores = 40,
    model = "m_1",
    N = c(500, 500, 500, 500, 500),
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
  
  if(unname(Sys.info()[1]) == "Darwin"){
    log_info("On mac, reset cores to 5")
    l_spec$mc_cores <- 5
  } 
  
  l_spec$l_e$p_d1_alloc <- 0.15
  l_spec$l_e$p_d2_entry <- 0.7
  l_spec$l_e$p_d2_alloc <- 0.5
  l_spec$l_e$p_d3_entry <- 0.9
  l_spec$l_e$p_d3_alloc <- 0.5
  l_spec$l_e$p_d4_entry <- 0.6
  l_spec$l_e$p_d4_alloc <- 0.5
  # preference for two-stage
  l_spec$l_e$p_pref <- 0.35
  
  l_spec$l_l$p_d1_alloc <- 0.5
  l_spec$l_l$p_d2_entry <- 0.7
  l_spec$l_l$p_d2_alloc <- 0.5
  l_spec$l_l$p_d3_entry <- 0.9
  l_spec$l_l$p_d3_alloc <- 0.5
  l_spec$l_l$p_d4_entry <- 0.6
  l_spec$l_l$p_d4_alloc <- 0.5
  # preference for two-stage
  l_spec$l_l$p_pref <- 0.7
  
  l_spec$l_c$p_d1_alloc <- 0.8
  l_spec$l_c$p_d2_entry <- 0.7
  l_spec$l_c$p_d2_alloc <- 0.5
  l_spec$l_c$p_d3_entry <- 0.9
  l_spec$l_c$p_d3_alloc <- 0.5
  l_spec$l_c$p_d4_entry <- 0.6
  l_spec$l_c$p_d4_alloc <- 0.5
  # preference for two-stage
  l_spec$l_c$p_pref <- 0.75
  
  # model params
  l_spec$mu <- 0.789
  l_spec$bs <- c(0.0, -0.3, -0.2)
  l_spec$bp <- c(0, -0.2)
  # dair, one, two-stage, we compare avg of one and two stage rev to dair
  l_spec$l_e$bd1 <- c(0, 0.1, 0.1)
  l_spec$l_l$bd1 <- c(-0.1, 0.62, 0.62)
  l_spec$l_c$bd1 <- c(-0.1, 0, 0.1)
  # always ref, 12wk, 6wk as we are assessing if 6wk ni to 12wk
  l_spec$bd2 <- c(0.0, 0.0, 0.0)
  # always ref, 0, 12wk as we are assessing if 12wk sup to none
  l_spec$bd3 <- c(0.0, 0.0, 0.0)
  # always ref, none, rif as we are assessing if rif is sup to none
  l_spec$bd4 <- c(0.0, 0.0, 0.0)
  
  l_spec$prior <- list()
  # location, scale
  l_spec$prior$mu <- c(0.7, 0.7)
  l_spec$prior$bs <- c(0, 1)
  l_spec$prior$bp <- c(0, 1)
  l_spec$prior$bd1 <- c(0, 1)
  l_spec$prior$bd2 <- c(0, 1)
  l_spec$prior$bd3 <- c(0, 1)
  l_spec$prior$bd4 <- c(0, 1)
  
  l_spec$delta <- list()
  l_spec$delta$sup <- 0.0
  l_spec$delta$sup_fut <- 0.05
  l_spec$delta$ni <- -0.05
  l_spec$delta$ni_fut <- 0.0
  
  # domain specific
  l_spec$thresh <- list()
  l_spec$thresh$sup <- c(0.96, 0.96, 0.96, 0.99)
  l_spec$thresh$fut_sup <- c(0.3, 0.3, 0.3, 0.3)
  l_spec$thresh$ni <- c(0.975, 0.945, 0.975, 0.975)
  l_spec$thresh$fut_ni <- c(0.25, 0.1, 0.25, 0.25)
  
  log_info("Starting simulation with following parameters:");
  log_info("N: ", paste0(l_spec$N, collapse = ", "));
  log_info("b_silo: ", paste0(l_spec$bs, collapse = ", "));
  log_info("b_pref: ", paste0(l_spec$bp, collapse = ", "));
  log_info("b_d1: ", paste0(c(l_spec$l_e$bd1, l_spec$l_l$bd1, l_spec$l_c$bd1), collapse = ", "));
  log_info("b_d2: ", paste0(l_spec$bd2, collapse = ", "));
  log_info("b_d3: ", paste0(l_spec$bd3, collapse = ", "));
  log_info("b_d4: ", paste0(l_spec$bd4, collapse = ", "));
  
  l_spec$nex <- 3
  l_spec$ex_trial_ix <- c(1, 2, 3)
  
  l_spec
}



# Example ------
sim_07_par_sim <- function(){
  
  # Rscript --vanilla ./R/sim-07.R sim_07_par_sim ./sim07/cfg-sim07-sc01-v12.yml
  
  if(unname(Sys.info()[1]) == "Darwin"){
    log_info("On mac, reset cores to 5")
    mc_cores <- 5
  } else {
    mc_cores <- 50
  }
  
  ix <- 1
  m1 <- cmdstanr::cmdstan_model("stan/model-sim-07-a.stan")
  
  output_dir_mcmc <- paste0(getwd(), "/tmp")
  
  l_spec <- list(
    N = 1,
    # silo allocation
    p_s_alloc = c(0.3, 0.5, 0.2),
    l_e = list(),
    l_l = list(),
    l_c = list()
  )
  # N by analysis
  l_spec$N <- 5000
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
  
  # str(l_spec)
  
  log_info("Starting simulation with following parameters:");
  log_info("N: ", paste0(l_spec$N, collapse = ", "));
  log_info("b_silo: ", paste0(l_spec$bs, collapse = ", "));
  log_info("b_pref: ", paste0(l_spec$bp, collapse = ", "));
  log_info("b_d1: ", paste0(c(l_spec$l_e$bd1, l_spec$l_l$bd1, l_spec$l_c$bd1), collapse = ", "));
  log_info("b_d2: ", paste0(l_spec$bd2, collapse = ", "));
  log_info("b_d3: ", paste0(l_spec$bd3, collapse = ", "));
  log_info("b_d4: ", paste0(l_spec$bd4, collapse = ", "));
  
  N_sim <- 100
  message("Start sim with ", N_sim, " iterations using the following sim params")
  message(" l_e$bd1: ", paste0(l_spec$l_e$bd1, collapse = ", "))
  message(" l_l$bd1: ", paste0(l_spec$l_l$bd1, collapse = ", "))
  message(" l_c$bd1: ", paste0(l_spec$l_c$bd1, collapse = ", "))
  
  message(" bd2: ", paste0(l_spec$bd2, collapse = ", "))
  message(" bd3: ", paste0(l_spec$bd3, collapse = ", "))
  message(" bd4: ", paste0(l_spec$bd4, collapse = ", "))
  message("")
  
  r <- pbapply::pblapply(X=1:N_sim, cl = mc_cores, FUN = function(ix){
    
    d <- get_sim07_trial_data(
      l_spec,
      dec_sup = list(surg = NA, ext_proph = NA, ab_choice = NA),
      dec_ni = list(ab_dur = NA),
      dec_sup_fut = list(surg = NA, ext_proph = NA, ab_choice = NA),
      dec_ni_fut = list(ab_dur = NA)
    )
    # d[]
    
    # create stan data format based on the relevant subsets of pt
    lsd <- get_sim07_stan_data(d)
    
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
    # only care about mean so relatively few samples
    f_1 <- m1$sample(
      lsd$ld, iter_warmup = 1000, iter_sampling = 2000,
      parallel_chains = 1, chains = 1, refresh = 0, show_exceptions = F,
      max_treedepth = 11,
      output_dir = output_dir_mcmc,
      output_basename = foutname
    )
    
    # f_1 <- m1$pathfinder(
    #   lsd$ld,
    #   num_paths=20,
    #   single_path_draws=200,
    #   history_size=50,
    #   max_lbfgs_iters=100,
    #   refresh = 0,
    #   draws = 1000,
    #   output_dir = output_dir_mcmc,
    #   output_basename = foutname
    #   )
    
    # extract posterior - marginal probability of outcome by group
    # dair vs rev
    d_post <- data.table(f_1$draws(
      variables = c(
        
        "rd_d1", 
        "rd_d2",
        "rd_d3",
        "rd_d4",
        "p_d1_1", "p_d1_2", "p_d1_3", "p_d1_23",
        "p_d2_2", "p_d2_3",
        "p_d3_2", "p_d3_3",
        "p_d4_2", "p_d4_3"
        
      ),   # risk scale
      format = "matrix"))
    d_post[, rd_d1_1 := p_d1_2 - p_d1_1]
    d_post[, rd_d1_2 := p_d1_3 - p_d1_1]
    colMeans(d_post)
    
  })
  d_fig <- data.table(do.call(rbind, r))
  
  message("Expectations for posterior means for estimated parameters:")
  d_fig <- melt(d_fig, measure.vars = names(d_fig))
  d_out <- d_fig[, .(mu = mean(value)), keyby = variable]
  for(i in 1:nrow(d_out)){
    message(sprintf("%s, %.3f", d_out[i, variable], d_out[i, mu]))
  }
  
  #ggplot(d_fig, aes(x = value, group = variable)) +
  #  geom_density() +
  #  geom_vline(data = d_fig[, .(mu = mean(value)), keyby = variable],
  #             aes(xintercept = mu), lwd = 0.2) +
  #  geom_point(data = d_fig[, .(mu = mean(value)), keyby = variable],
  #             aes(x = mu, y = 0), size = 0.001) +
  #  facet_wrap2(~variable, scales = "free_x") +
  #  geom_text_repel(
  #    data = d_fig[, .(mu = mean(value)), keyby = variable],
  #    aes(x = mu, y = 0, label = sprintf("%.3f", mu)),
  #    inherit.aes = F,
  #    force             = 0.3,
  #    nudge_y           = 0.15,
  #    direction         = "y",
  #    hjust             = 0,
  #    segment.size      = 0.2,
  #    segment.curvature = -0.1,
  #    size = 3
  #  ) +
  #  theme(
  #    axis.title = element_blank(),
  #    plot.title = element_text(size = 10),
  #    plot.subtitle = element_text(size = 8)
  #  ) 
  #
  
  # 
  # file.remove("tmp/*.csv")
}

sim_07_subgrp_dev <- function(){
  
  
  if(unname(Sys.info()[1]) == "Darwin"){
    log_info("On mac, reset cores to 5")
    mc_cores <- 5
  } else {
    mc_cores <- 50
  }
  
  # Classical setup ----
  
  set.seed(2254)
  dX <- CJ(domA = factor(1:3), s = factor(1:3), j = factor(1:2))
  X <- model.matrix(~ domA + s + j + domA:s + domA:j, data = dX)
  b <- c(
    qlogis(0.4),
    0.2,  0.1, # trt
    -0.1, -0.4, # silo
    -0.2, # joint
    
    # intxn shift rev(1) & late
    # intxn shift rev(2) & late
    # intxn shift rev(1) & chronic
    # intxn shift rev(2) & chronic
    rnorm(4, 0, 0.3), 
    
    # intxn shift rev(1) & hip
    # intxn shift rev(2) & hip
    rnorm(2, 0, 0.3) 
  )
  p_tru <- as.numeric(plogis(X %*% b))
  
  # large sample estimates approximate true values
  N <- 2e6
  d <- data.table(
    domA = factor(sample(1:3, size = N, replace = T)),
    s = factor(sample(1:3, size = N, replace = T)),
    j = factor(sample(1:2, size = N, replace = T))
  )
  X <- model.matrix(~ domA + s + j + domA:s + domA:j, data = d)
  d[, eta :=  X %*% b]
  d[, y := rbinom(.N, 1, plogis(eta))]
  
  f_1 <- glm(y ~ domA + s + j + domA:s + domA:j, data = d, family = binomial)
  
  d_out <- data.table(
    p_1 = predict(
      f_1, 
      newdata = data.table(domA = factor(1, levels = 1:3), s = d$s, j = d$j), 
      type = "response"),
    p_2 = predict(
      f_1, 
      newdata = data.table(domA = factor(2, levels = 1:3), s = d$s, j = d$j), 
      type = "response"),
    p_3 = predict(
      f_1, 
      newdata = data.table(domA = factor(3, levels = 1:3), s = d$s, j = d$j), 
      type = "response"),
    s = d$s,
    j = d$j
  )
  d_smry <- d_out[
    , .(mu_p_1 = mean(p_1), mu_p_2 = mean(p_2), mu_p_3 = mean(p_3)), keyby = s]
  d_smry[, rd_2_1 := mu_p_2 - mu_p_1]
  d_smry[, rd_3_1 := mu_p_3 - mu_p_1]
  d_smry[]
  # strata level risk in treatment arm 1, 2, 3 then rd 2 vs 1 and 3 vs 1
  # for example:
  # Key: <s>
  #   s    mu_p_1    mu_p_2    mu_p_3       rd_2_1      rd_3_1
  # <fctr>     <num>     <num>     <num>        <num>       <num>
  # 1:      1 0.3774311 0.4216500 0.4054096  0.044218902  0.02797854
  # 2:      2 0.3536127 0.3464584 0.3904160 -0.007154337  0.03680327
  # 3:      3 0.2871736 0.4192520 0.2569424  0.132078417 -0.03023120
  
  data.table(marginaleffects::avg_predictions(f_1, by = c("domA", "s")))[order(s, domA)]
  # comparisons(f_1,  by = c("domA", "s"))
  data.table(marginaleffects::comparisons(f_1, variables = "domA", by = "s"))[
    order(s), .(contrast, s, estimate)]
  # contrast      s     estimate
  # <char> <fctr>        <num>
  # 1: mean(2) - mean(1)      1  0.044218902
  # 2: mean(3) - mean(1)      1  0.027978538
  # 3: mean(2) - mean(1)      2 -0.007154337
  # 4: mean(3) - mean(1)      2  0.036803269
  # 5: mean(2) - mean(1)      3  0.132078418
  # 6: mean(3) - mean(1)      3 -0.030231201
  
  d_smry <- d_out[
    , .(mu_p_1 = mean(p_1), mu_p_2 = mean(p_2), mu_p_3 = mean(p_3)), keyby = j]
  d_smry[, rd_2_1 := mu_p_2 - mu_p_1]
  d_smry[, rd_3_1 := mu_p_3 - mu_p_1]
  d_smry[]
  # Key: <j>
  #   j    mu_p_1    mu_p_2    mu_p_3     rd_2_1      rd_3_1
  # <fctr>     <num>     <num>     <num>      <num>       <num>
  #   1:      1 0.3619614 0.4217142 0.3701458 0.05975279 0.008184397
  # 2:      2 0.3169593 0.3698630 0.3318754 0.05290368 0.014916107
  
  data.table(marginaleffects::comparisons(f_1, variables = "domA", by = "j"))[
    order(j), .(contrast, j, estimate)]
  # contrast      j    estimate
  # <char> <fctr>       <num>
  #   1: mean(2) - mean(1)      1 0.059752789
  # 2: mean(3) - mean(1)      1 0.008184397
  # 3: mean(2) - mean(1)      2 0.052903684
  # 4: mean(3) - mean(1)      2 0.014916107
  
  
  # Bayesian setup --------
  
  # fixed main effects with random effect for subgroup/interaction
  m2 <- cmdstanr::cmdstan_model("stan/model-sim-07-a-subgrp1.stan")
  
  # looking at a large sample comparison, but don't expect exact match
  N <- 1e6
  d <- data.table(
    domA = factor(sample(1:3, size = N, replace = T)),
    s = factor(sample(1:3, size = N, replace = T)),
    j = factor(sample(1:2, size = N, replace = T))
  )
  X <- model.matrix(~ domA + s + j + domA:s + domA:j, data = d)
  # use same beta from above
  d[, eta :=  X %*% b]
  d[, p_tru := plogis(eta)]
  d[, y := rbinom(.N, 1, p_tru)]
  
  # aggregate - otherwise takes for ever
  d_mod <- d[, .(y = sum(y), n = .N), keyby = .(domA, s, j)]
  
  ld <- list(
    N = nrow(d_mod),
    y = d_mod$y, n = d_mod$n,
    domA = d_mod$domA, s = d_mod$s, j = d_mod$j,
    N_s = d_mod[, .N, keyby = s][, N],
    N_j = d_mod[, .N, keyby = j][, N],
    ix_s1 = d_mod[s == 1, , which = T],
    ix_s2 = d_mod[s == 2, , which = T],
    ix_s3 = d_mod[s == 3, , which = T],
    ix_j1 = d_mod[j == 1, , which = T],
    ix_j2 = d_mod[j == 2, , which = T]
  )
  
  # m2 <- cmdstanr::cmdstan_model("stan/model-sim-07-a-subgrp1.stan")
  f_2 <- m2$sample(
    ld, iter_warmup = 1000, iter_sampling = 2000,
    parallel_chains = 1, chains = 1, refresh = 1000, show_exceptions = F,
    max_treedepth = 13
  )
  
  f_3 <- glm(cbind(y, n-y) ~ domA + s + j + domA:s + domA:j, data = d_mod, family = binomial)
  
  # results and compare to classical without shrinkage
  d_post <- data.table(
    f_2$draws(
      variables = c(
        "mu_p_dA_s1", "mu_p_dA_s2", "mu_p_dA_s3"), format = "matrix"))
  d_post <- melt(d_post, measure.vars = names(d_post))
  d_out <- data.table(marginaleffects::avg_predictions(f_3, by = c("domA", "s")))
  
  # average risk in silo 1 for treatment arm 1, 2, 3
  # average risk in silo 2 for treatment arm 1, 2, 3
  # average risk in silo 3 for treatment arm 1, 2, 3
  cbind(
    d_post[, .(mu = mean(value)), keyby = variable],
    d_out[order(s, domA), .(s, domA, estimate, conf.low, conf.high)]
  )
  
  d_post <- data.table(
    f_2$draws(
      variables = c("mu_rd_s1", "mu_rd_s2", "mu_rd_s3"), format = "matrix"))
  d_post <- melt(d_post, measure.vars = names(d_post))
  
  # silo specific risk diff trt effects (rev(1) vs dair, rev(2) vs dair)
  d_out <- data.table(marginaleffects::comparisons(f_1, variables = "domA", by = "s"))
  cbind(
    d_post[, .(mu = mean(value)), keyby = variable],
    d_out[order(s), .(contrast, s, estimate, conf.low, conf.high)]
  )
  
  
  
  # Simulation --------
  
  N <- 2500
  r <- pbapply::pblapply(X=1:2000, cl = mc_cores, FUN = function(ix){
    
    d <- data.table(
      domA = factor(sample(1:3, size = N, replace = T)),
      s = factor(sample(1:3, size = N, replace = T)),
      j = factor(sample(1:2, size = N, replace = T))
    )
    X <- model.matrix(~ domA + s + j + domA:s + domA:j, data = d)
    d[, eta :=  X %*% b]
    d[, p_tru := plogis(eta)]
    d[, y := rbinom(.N, 1, p_tru)]
    # d[, .(p_obs = mean(y), p_tru = unique(p_tru)), keyby = .(domA, s, j)]
    
    d_mod <- d[, .(y = sum(y), n = .N), keyby = .(domA, s, j)]
    
    ld <- list(
      N = nrow(d_mod),
      y = d_mod$y, n = d_mod$n,
      domA = d_mod$domA, s = d_mod$s, j = d_mod$j,
      N_s = d_mod[, .N, keyby = s][, N],
      N_j = d_mod[, .N, keyby = j][, N],
      ix_s1 = d_mod[s == 1, , which = T],
      ix_s2 = d_mod[s == 2, , which = T],
      ix_s3 = d_mod[s == 3, , which = T],
      ix_j1 = d_mod[j == 1, , which = T],
      ix_j2 = d_mod[j == 2, , which = T]
    )
    
    # m2 <- cmdstanr::cmdstan_model("stan/model-sim-07-a-subgrp1.stan")
    f_4 <- m2$sample(
      ld, iter_warmup = 1000, iter_sampling = 2000,
      parallel_chains = 1, chains = 1, refresh = 0, show_exceptions = F,
      max_treedepth = 13
    )
    
    f_5 <- glm(
      cbind(y, n-y) ~ domA + s + j + domA:s + domA:j, 
      data = d_mod, family = binomial)
    
    d_post <- data.table(
      f_4$draws(
        variables = c("mu_rd_s1", "mu_rd_s2", "mu_rd_s3",
                      "mu_rd_j1", "mu_rd_j2"), format = "matrix"))
    d_post <- melt(d_post, measure.vars = names(d_post))
    # silo specific risk diff trt effects:
    # silo 1 (rev(1) vs dair, rev(2) vs dair)
    # silo 2 (rev(1) vs dair, rev(2) vs dair)
    # silo 3 (rev(1) vs dair, rev(2) vs dair)
    # joint 1 (rev(1) vs dair, rev(2) vs dair)
    # joint 2 (rev(1) vs dair, rev(2) vs dair)
    d_post <- d_post[, .(mu = mean(value)), keyby = variable]
    
    d_out_s <- data.table(marginaleffects::comparisons(f_5, variables = "domA", by = "s"))
    d_out_j <- data.table(marginaleffects::comparisons(f_5, variables = "domA", by = "j"))
    
    cbind(
      d_post, 
      rbind(d_out_s[order(s), .(contrast, subgrp = s, estimate)], 
            d_out_j[order(j), .(contrast, subgrp = j, estimate)])
    )
    
  })
  
  d_res <- rbindlist(r, idcol = "id_sim")
  d_res[, .(mu_b = mean(mu), mu_f = mean(estimate)), keyby = variable]
  
  
  
  
  
  
  
  
  
  
  
  
  
  
}



# Run loop -------
sim07_sim_loop <- function(){
  
  log_info(paste0(match.call()[[1]]))
  
  default_cfg <- F
  if(!default_cfg){
    # load sim specification
    f_spec <- here::here("./etc", args[2])
    l_spec <- config::get(file = f_spec)
    stopifnot("Config is null" = !is.null(l_spec))
    l_spec <- sim07_update_cfg(l_spec)
  } else {
    l_spec <- sim07_default_cfg()
  }
  
  # str(l_spec)
  l_spec$return_posterior = F; e = NULL; ix <- 1
  log_info("Starting simulation")
  
  r <- parallel::mclapply(
    X=1:l_spec$n_sim, mc.cores = l_spec$mc_cores, FUN=function(ix) {
      log_info("Simulation ", ix);
      
      l_spec$ix_sim <- ix
      if(ix %in% l_spec$ex_trial_ix){
        l_spec$return_posterior = T  
      } else {
        l_spec$return_posterior = F
      }
      
      ll <- tryCatch({run_trial(l_spec)},
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
  log_info("Sleep for 1 before processing")
  Sys.sleep(1)
  
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
  
  
  d_post_all <- data.table(do.call(rbind, lapply(1:length(r), function(i){
    # if the sim contains full posterior (for example trial) then return
    if(!is.null(r[[i]]$d_post_all)){
      cbind(sim = i, r[[i]]$d_post_all)
    }
    
  } )))
  
  l <- list(
    l_spec = l_spec,
    
    # reference risk
    # d_risk_smry = d_risk_smry,
    
    d_pr_sup = d_pr_sup, 
    d_pr_ni = d_pr_ni,
    
    d_pr_sup_fut = d_pr_sup_fut, 
    d_pr_ni_fut = d_pr_ni_fut,
    
    d_decision = d_decision,
    d_post_smry_1 = d_post_smry_1,
    d_post_smry_2 = d_post_smry_2,
    
    d_all = d_all,
    
    d_post_all = d_post_all
    
    # d_n_units = d_n_units,
    # d_n_assign = d_n_assign,
    # d_grp = d_grp
  )
  
  toks <- unlist(tstrsplit(args[2], "[-.]"))
  fname <- paste0("data/sim07/sim07-", toks[4], "-", toks[5], "-", format(Sys.time(), "%Y%m%d-%H%M%S"), ".qs2")
  
  log_info("Saving results file ", fname)
  
  qs2::qs_save(l, file = fname)
}

sim07_run_none <- function(){
  log_info("sim07_run_none: Nothing doing here bud.")
}

sim07_main <- function(){
  funcname <- paste0(args[1], "()")
  log_info("Main, invoking ", funcname)
  eval(parse(text=funcname))
}



if(!interactive()){
  sim07_main()
}