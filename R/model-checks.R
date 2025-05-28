library(data.table)
library(pbapply)
library(ggplot2)
library(patchwork)




derive_parameters <- function(
  p_d1_alloc = 0.5,
  p_d2_entry = 0.7,
  p_d2_alloc = 0.5,
  p_d3_entry = 0.9,
  p_d3_alloc = 0.5,
  # 40% pref for two-stage
  p_pref = 0.4,
  # log-odds
  mu = 0.7,
  # log-or
  # different baseline risk for rev
  bp = -0.4,
  bd1 = c(0, 0, 0),
  bd2 = c(0, 0.1, 0),
  bd3 = c(0, 0, -0.2)
  
  ){
  
  # dair
  prop_tru <- c(1 - p_pref, p_pref)
  p_dair <- prop_tru[1] * plogis(mu) + prop_tru[2] * plogis(mu + bp)
  
  # rev 1
  prop_tru <- c(1 - p_d2_entry, rep(p_d2_entry/2, 2))
  p_rev_1 <- prop_tru[1] * plogis(mu + bd1[2]) + 
    # contribution of d2 12 wk
    prop_tru[2] * plogis(mu + bd1[2] + bd2[2]) + 
    # contribution of d2 6 wk
    prop_tru[3] * plogis(mu + bd1[2] + bd2[3])
  
  # rev 2
  prop_tru <- c(1-p_d3_entry, rep(p_d3_entry/2, 2))
  # contribution of non-rand, bp needs to be included here
  p_rev_2 <- prop_tru[1] * plogis(mu + bp + bd1[3]) + 
    # contribution of randomised interventions
    prop_tru[2] * plogis(mu + bp + bd1[3] + bd3[2]) + 
    prop_tru[3] * plogis(mu + bp + bd1[3] + bd3[3])
  
  # rev
  prop_tru <- c(1 - p_pref, p_pref)
  p_rev <- prop_tru[1] * p_rev_1 + prop_tru[2] * p_rev_2
  
  
  c(rd = p_rev - p_dair,
    p_dair = p_dair, p_rev = p_rev, 
    p_rev_1 = p_rev_1, p_rev_2 = p_rev_2)
  
}



test_pars <- function(){
  
  p_d1_alloc = 0.5
  p_d2_entry = 0.7
  p_d2_alloc = 0.5
  p_d3_entry = 0.9
  p_d3_alloc = 0.5
  # 40% pref for two-stage
  p_pref = 0.4
  # log-odds
  mu = 0.7
  # log-or
  # different baseline risk for rev
  bp = -0.4
  bd1 = c(0, 0, 0)
  bd2 = c(0, 0.1, 0)
  bd3 = c(0, 0, -0.2)
  
  
  obj_f1 <- function(x){
    res <- derive_parameters(bd1 = c(0, x, x))
    (res["rd"] - 0)^2
  }
  obj_f2 <- function(x){
    res <- derive_parameters(bd1 = c(0, 0, x))
    (res["rd"] - 0)^2
  }
  obj_f3 <- function(x){
    res <- derive_parameters(bd1 = c(0, x, 0))
    (res["rd"] - 0)^2
  }
  
  # for the null on the risk scale you need to have both bd1 parameters set to
  (lor <- optimize(obj_f1, c(0, 1))$minimum)
  derive_parameters(bd1 = c(0, lor, lor))
  # for the null on the risk scale you need to have bd1[2] = 0 and bd1[3] 
  optimize(obj_f2, c(0, 1))$minimum
  # for the null on the risk scale you need to have bd1[3] = 0 and bd1[2] 
  optimize(obj_f3, c(0, 1))$minimum
  
}



multi_model_approach <- function(){
  
  N_sim <- 1000
  mc_cores <- 4
  
  r <- pbapply::pblapply(X=1:N_sim, cl = mc_cores, FUN = function(ix){
    
    N <- 1e4
    d <- data.table(
      # ctl/trt allocations ignoring dependencies
      d1_alloc = rbinom(N, 1, 0.5),
      # 70% entery d2
      d2_entry = rbinom(N, 1, 0.7),
      d2_alloc = rbinom(N, 1, 0.5),
      # 90% enter d3
      d3_entry = rbinom(N, 1, 0.9),
      d3_alloc = rbinom(N, 1, 0.5),
      # sev is a proxy for some set of pt characteristics
      # sev = rnorm(N, 0.5, 0.8),
      pref = rbinom(N, 1, 0.4)
    )
    # preference directs type of revision
    # d[, pref := rbinom(N, 1, prob = plogis(sev))]
    
    # dair gets dair, revision gets split
    d[d1_alloc == 0, d1 := 1]
    d[d1_alloc == 1 & pref == 0, d1 := 2]
    d[d1_alloc == 1 & pref == 1, d1 := 3]
    
    d[d1 == 2 & d2_entry == 0, d2 := 1]
    d[d1 == 2 & d2_entry == 1, d2 := 2 + d2_alloc]
    
    d[d1 == 3 & d3_entry == 0, d3 := 1]
    d[d1 == 3 & d3_entry == 1, d3 := 2 + d3_alloc]
    
    d[, .N, keyby = .(pref, d1, d2, d3)]
    
    mu <- 0.7
    # different baseline risk for rev
    bp <- -0.4
    # actually the null case when you think about it on the risk scale
    bd1 <- c(0, 0.0188705, 0.0188705)
    bd2 <- c(0, 0.1, 0)
    bd3 <- c(0, 0, -0.2)
    
    # bd1 is irrelevant since d1 == 1 is the ref group, fixed at zero
    d[d1 == 1, eta := mu + bp*pref]
    # pref is irrelevant as d1 = 2 only occurs if pref = 0
    d[d1 == 2, eta := mu + bd1[2] + bd2[d2]]
    # but here pref is relevant as d1 = 3 only if pref = 1
    d[d1 == 3, eta := mu + bp*pref + bd1[3] + bd3[d3]]
    
    d[, `:=`(d1 = factor(d1), d2 = factor(d2), d3 = factor(d3))]
    
    d[, p := plogis(eta)]
    d[, y := rbinom(N, 1, p)]
    
    f1_1 <- glm(
      y ~ 1 + pref  , data = d, subset= d1==1, family = binomial,
      control = glm.control(epsilon = 1e-4, maxit = 100, trace = FALSE))
    # pref is zero for everyone
    f1_2 <- glm(
      y ~ 1 + d2 , data = d, subset= d1==2, family = binomial,
      control = glm.control(epsilon = 1e-4, maxit = 100, trace = FALSE))
    # effect of pref will get rolled into the intercept
    f1_3 <- glm(
      y ~ 1 + d3 , data = d, subset= d1==3, family = binomial,
      control = glm.control(epsilon = 1e-4, maxit = 100, trace = FALSE))
    
    # lapply(list(f1_1, f1_2, f1_3), function(z) coef(z))
    # dair
    d_dair <- copy(d)
    d_dair[, d1 := factor(1)]
    d_dair[, p_hat := predict(f1_1, newdata = d_dair, type = "response")]
    
    # averages over preference
    prop_tru <- c(0.6, 0.4)
    p_dair_ref <- prop_tru[1] * plogis(mu) + prop_tru[2] * plogis(mu + bp)
    d_prop <- d[d1 == 1, .(prop = .N/nrow(d[d1==1])), keyby = pref]
    p_dair_hat <- mean(d_dair$p_hat)
    
    # rev(1)
    d_rev_1 <- copy(d)
    d_rev_1[, d1 := factor(2)]
    d_rev_1[, pref := 0]
    d_rev_1[d2_entry == 1, d2 := factor(d[d2_entry == 1, d2_alloc + 2], levels = 1:3)]
    d_rev_1[d2_entry == 0, d2 := factor(1, levels = 1:3)]
    d_rev_1[, p_hat := predict(f1_2, newdata = d_rev_1, type = "response")]
    
    # average over the d2 treatments
    prop_tru <- c(0.3, 0.35, 0.35)
    # contribution of non-randomised d2
    p_rev_1_ref <- prop_tru[1] * plogis(mu + bd1[2]) + 
      # contribution of d2 12 wk
      prop_tru[2] * plogis(mu + bd1[2] + bd2[2]) + 
      # contribution of d2 6 wk
      prop_tru[3] * plogis(mu + bd1[2] + bd2[3])
    p_rev_1_hat <- mean(d_rev_1$p_hat)
    
    # # rev(2)
    d_rev_2 <- copy(d)
    d_rev_2[, d1 := factor(3)]
    # pref doesn't matter as it is rolled up into the intercept
    # d_rev_2[, pref := 1]
    d_rev_2[d3_entry == 1, d3 := factor(d[d3_entry == 1, d3_alloc + 2], levels = 1:3)]
    d_rev_2[d3_entry == 0, d3 := factor(1, levels = 1:3)]
    d_rev_2[, p_hat := predict(f1_3, newdata = d_rev_2, type = "response")]
    
    prop_tru <- c(0.1, 0.45, 0.45)
    # contribution of non-rand, bp needs to be included here
    p_rev_2_ref <- prop_tru[1] * plogis(mu + bp + bd1[3]) + 
      # contribution of randomised interventions
      prop_tru[2] * plogis(mu + bp + bd1[3] + bd3[2]) + 
      prop_tru[3] * plogis(mu + bp + bd1[3] + bd3[3])
    p_rev_2_hat <- mean(d_rev_2$p_hat)
    
    
    # need to update the 0.6 and 0.4 if the dist of pref is changed
    prop_tru <- c(0.6, 0.4)
    p_rev_ref <- prop_tru[1] * p_rev_1_ref + prop_tru[2] * p_rev_2_ref 
    # now aggregate the one and two-stage revision
    d_prop <- d[, .(prop = .N/nrow(d)), keyby = pref]
    p_rev_hat <- p_rev_1_hat * d_prop[pref == 0, prop] + p_rev_2_hat * d_prop[pref == 1, prop]
    
    rd_ref <- p_rev_ref - p_dair_ref
    rd_hat <- p_rev_hat - p_dair_hat
    
    d_uniq <- unique(d[, .(pref, d1, d2, d3, p)])
    setkey(d_uniq, pref, d1, d2, d3)
    d_uniq[d1 == 1, p_hat := predict(f1_1, newdata = d_uniq[d1 == 1], type = "response")]
    d_uniq[d1 == 2, p_hat := predict(f1_2, newdata = d_uniq[d1 == 2], type = "response")]
    d_uniq[d1 == 3, p_hat := predict(f1_3, newdata = d_uniq[d1 == 3], type = "response")]
    
    
    list(
      d_uniq = d_uniq,
      d_g = data.table(
        p_dair_ref = p_dair_ref,
        p_dair_hat = p_dair_hat,
        p_rev_1_ref = p_rev_1_ref,
        p_rev_1_hat = p_rev_1_hat,
        p_rev_2_ref = p_rev_2_ref,
        p_rev_2_hat = p_rev_2_hat,
        p_rev_ref = p_rev_ref,
        p_rev_hat = p_rev_hat,
        rd_ref = rd_ref,
        rd_hat = rd_hat
      )
    )
    
  })
  
  
  d_g <- rbindlist(lapply(r, function(z) z$d_g))
  
  d_g[, bias_rd := rd_hat - rd_ref]
  d_g[, bias_p_dair := p_dair_hat - p_dair_ref]
  d_g[, bias_p_rev := p_rev_hat - p_rev_ref]
  d_g[, bias_p_rev_1 := p_rev_1_hat - p_rev_1_ref]
  d_g[, bias_p_rev_2 := p_rev_2_hat - p_rev_2_ref]
  
  d_fig <- d_g[, .(p_dair_hat, p_rev_1_hat, p_rev_2_hat, p_rev_hat)]
  d_fig <- melt(d_fig, measure.vars = names(d_fig))
  
  p1 <- ggplot(d_fig, aes(x = value)) +
    geom_density() +
    geom_vline(data = d_fig[, .(mu = mean(value)), keyby = variable],
               aes(xintercept = mu), lwd = 0.2) + 
    facet_wrap(~variable, nrow = 1)
  #
  
  d_fig <- d_g[, .(p_dair_hat, p_rev_hat)]
  d_fig <- melt(d_fig, measure.vars = names(d_fig))
  
  p2 <- ggplot(d_fig, aes(x = value)) +
    geom_density() +
    geom_vline(data = d_fig[, .(mu = mean(value)), keyby = variable],
               aes(xintercept = mu), lwd = 0.2) +
    facet_wrap(~variable,  nrow = 1)
  
  d_fig <- d_g[, .(rd_hat)]
  p3 <- ggplot(d_fig, aes(x = rd_hat)) +
    geom_density() +
    geom_vline(data = d_fig[, .(mu = mean(rd_hat))],
               aes(xintercept = mu), lwd = 0.2) 
  
  p1 / (p2 + p3)
  
}




single_model_approach <- function(){
  
  N_sim <- 2500
  mc_cores <- 4
  
  r <- pbapply::pblapply(X=1:N_sim, cl = mc_cores, FUN = function(ix){
    
    N <- 1e6
    d <- data.table(
      # ctl/trt allocations ignoring dependencies
      d1_alloc = rbinom(N, 1, 0.5),
      # 70% entery d2
      d2_entry = rbinom(N, 1, 0.7),
      d2_alloc = rbinom(N, 1, 0.5),
      # 90% enter d3
      d3_entry = rbinom(N, 1, 0.9),
      d3_alloc = rbinom(N, 1, 0.5),
      # sev is a proxy for some set of pt characteristics
      # sev = rnorm(N, 0.5, 0.8),
      pref = rbinom(N, 1, 0.4)
    )
    # preference directs type of revision
    # d[, pref := rbinom(N, 1, prob = plogis(-sev))]
    
    # dair gets dair, revision gets split
    d[d1_alloc == 0, d1 := 1]
    d[d1_alloc == 1 & pref == 0, d1 := 2]
    d[d1_alloc == 1 & pref == 1, d1 := 3]
    
    d[d1 == 2 & d2_entry == 0, d2 := 1]
    d[d1 == 2 & d2_entry == 1, d2 := 2 + d2_alloc]
    
    d[d1 == 3 & d3_entry == 0, d3 := 1]
    d[d1 == 3 & d3_entry == 1, d3 := 2 + d3_alloc]
    
    mu <- 0.7
    # different baseline risk for rev
    bp <- -0.4
    bd1 <- c(0, 0, 0)
    bd2 <- c(0, 0.1, 0)
    bd3 <- c(0, 0, -0.2)
    
    # bp is irrelevant since it is multiplied by 0
    d[d1 == 1, eta := mu + bd1[1]]
    d[d1 == 2, eta := mu + bd1[2] + bd2[d2]]
    # bp is relevant since it is multiplied by 1
    d[d1 == 3, eta := mu + bp*pref + bd1[3] + bd3[d3]]
    
    # these are the undefined records, those for whom there is no d2 nor d3
    d[is.na(d2), d2 := 4]
    d[is.na(d3), d3 := 4]
    
    d[, `:=`(d1 = factor(d1), d2 = factor(d2), d3 = factor(d3))]
    
    d[, p := plogis(eta)]
    d[, y := rbinom(N, 1, p)]
    
    f2 <- glm(
      y ~ 1 + d1 + d2 + d3 , data = d, family = binomial,
      control = glm.control(epsilon = 1e-4, maxit = 100, trace = FALSE))
    
    summary(f2)
    # lapply(list(f1_1, f1_2, f1_3), function(z) summary(z)$coef)
    # dair
    d_dair <- copy(d)
    d_dair[, d1 := factor(1)]
    d_dair[, d2 := factor(4)]
    d_dair[, d3 := factor(4)]
    d_dair[, p_hat := predict(f2, newdata = d_dair, type = "response")]
    
    # averages over preference, the only variation for the d1 group
    d_prop <- d[d1 == 1, .(prop = .N/nrow(d[d1==1])), keyby = pref]
    p_dair_ref <- d_prop[pref == 0, prop] * plogis(mu) + d_prop[pref == 1, prop] * plogis(mu + bp)
    p_dair_hat <- mean(d_dair$p_hat)
    
    # rev(1)
    d_rev_1 <- copy(d)
    d_rev_1[, d1 := factor(2)]
    d_rev_1[, pref := 0]
    d_rev_1[d2_entry == 1, d2 := factor(d[d2_entry == 1, d2_alloc + 2], levels = 1:4)]
    d_rev_1[d2_entry == 0, d2 := factor(1)]
    d_rev_1[, d3 := factor(4)]
    d_rev_1[, p_hat := predict(f2, newdata = d_rev_1, type = "response")]
    
    d_prop <- d[d1 == 2, .(prop = .N/nrow(d[d1==2])), keyby = d2]
    p_rev_1_ref <- d_prop[d2 == 1, prop] * plogis(mu + bd1[2]) + 
      d_prop[d2 == 2, prop] * plogis(mu + bd1[2] + bd2[2]) + 
      d_prop[d2 == 3, prop] * plogis(mu + bd1[2] + bd2[3])
    p_rev_1_hat <- mean(d_rev_1$p_hat)
    
    # # rev(2)
    d_rev_2 <- copy(d)
    d_rev_2[, d1 := factor(3)]
    d_rev_2[, pref := 1]
    d_rev_2[, d2 := factor(4)]
    d_rev_2[d3_entry == 1, d3 := factor(d[d3_entry == 1, d3_alloc + 2], levels = 1:4)]
    d_rev_2[d3_entry == 0, d3 := factor(1, levels = 1:4)]
    d_rev_2[, p_hat := predict(f2, newdata = d_rev_2, type = "response")]
    
    d_prop <- d[d1 == 3, .(prop = .N/nrow(d[d1==3])), keyby = d3]
    p_rev_2_ref <- d_prop[d3 == 1, prop] * plogis(mu + bp + bd1[3]) + 
      d_prop[d3 == 2, prop] * plogis(mu + bp + bd1[3] + bd3[2]) + 
      d_prop[d3 == 3, prop] * plogis(mu + bp + bd1[3] + bd3[3])
    p_rev_2_hat <- mean(d_rev_2$p_hat)
    
    d_prop <- d[, .(prop = .N/nrow(d)), keyby = pref]
    p_rev_ref <- p_rev_1_ref * d_prop[pref == 0, prop] + p_rev_2_ref * d_prop[pref == 1, prop] 
    p_rev_hat <- p_rev_1_hat * d_prop[pref == 0, prop] + p_rev_2_hat * d_prop[pref == 1, prop]
    
    rd_ref <- p_rev_ref - p_dair_ref
    rd_hat <- p_rev_hat - p_dair_hat
    
    d_uniq <- unique(d[, .(pref, d1, d2, d3, p)])
    setkey(d_uniq, pref, d1, d2, d3)
    d_uniq[d1 == 1, p_hat := predict(f2, newdata = d_uniq[d1 == 1], type = "response")]
    d_uniq[d1 == 2, p_hat := predict(f2, newdata = d_uniq[d1 == 2], type = "response")]
    d_uniq[d1 == 3, p_hat := predict(f2, newdata = d_uniq[d1 == 3], type = "response")]
    
    
    list(
      d_uniq = d_uniq,
      d_g = data.table(
        p_dair_ref = p_dair_ref,
        p_dair_hat = p_dair_hat,
        p_rev_1_ref = p_rev_1_ref,
        p_rev_1_hat = p_rev_1_hat,
        p_rev_2_ref = p_rev_2_ref,
        p_rev_2_hat = p_rev_2_hat,
        p_rev_ref = p_rev_ref,
        p_rev_hat = p_rev_hat,
        rd_ref = rd_ref,
        rd_hat = rd_hat
      )
    )
  })
  
  d_g <- rbindlist(lapply(r, function(z) z$d_g))
  
  d_g[, bias_rd := rd_hat - rd_ref]
  d_g[, bias_p_dair := p_dair_hat - p_dair_ref]
  d_g[, bias_p_rev := p_rev_hat - p_rev_ref]
  d_g[, bias_p_rev_1 := p_rev_1_hat - p_rev_1_ref]
  d_g[, bias_p_rev_2 := p_rev_2_hat - p_rev_2_ref]
  
  d_fig <- d_g[, .(bias_rd, bias_p_dair, bias_p_rev, bias_p_rev_1, bias_p_rev_2)]
  d_fig <- melt(d_fig, measure.vars = names(d_fig))
  
  p1 <- ggplot(d_fig, aes(x = value)) +
    geom_density() +
    geom_vline(data = d_fig[, .(mu = mean(value)), keyby = variable],
               aes(xintercept = mu), lwd = 0.2) + 
    facet_wrap(~variable, nrow = 1)
  #
  
  d_fig <- d_g[, .(p_dair_hat, p_rev_hat, rd_hat)]
  d_fig <- melt(d_fig, measure.vars = names(d_fig))
  
  p2 <- ggplot(d_fig, aes(x = value)) +
    geom_density() +
    geom_vline(data = d_fig[, .(mu = mean(value)), keyby = variable],
               aes(xintercept = mu), lwd = 0.2) +
    facet_wrap(~variable, scales = "free", nrow = 1)
  
  p1 / p2
  
}

