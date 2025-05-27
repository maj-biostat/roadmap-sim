library(data.table)
library(pbapply)
library(ggplot2)

N_sim <- 2000
mc_cores <- 4

r <- pbapply::pblapply(X=1:N_sim, cl = mc_cores, FUN = function(ix){
  
  N <- 1e3
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
    sev = rnorm(N, 0.8, 0.8)
  )
  # preference directs type of revision
  d[, pref := rbinom(N, 1, prob = plogis(sev))]
  
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
  bd1 <- c(0, 0.5, 0.5)
  bd2 <- c(0, 0.1, 0)
  bd3 <- c(0, 0, -0.2)
  
  d[d1 == 1, eta := mu + bp*pref + bd1[1]]
  d[d1 == 2, eta := mu + bp*pref + bd1[2] + bd2[d2]]
  d[d1 == 3, eta := mu + bp*pref + bd1[3] + bd3[d3]]
  
  d[, `:=`(d1 = factor(d1), d2 = factor(d2), d3 = factor(d3))]
  
  # d[is.na(d2), d2 := 1]
  
  d[, p := plogis(eta)]
  d[, y := rbinom(N, 1, p)]
  
  f1_1 <- glm(
    y ~ 1 + pref  , data = d, subset= d1==1, family = binomial,
    control = glm.control(epsilon = 1e-4, maxit = 100, trace = FALSE))
  f1_2 <- glm(
    y ~ d2 , data = d, subset= d1==2, family = binomial,
    control = glm.control(epsilon = 1e-4, maxit = 100, trace = FALSE))
  f1_3 <- glm(
    y ~ d3 , data = d, subset= d1==3, family = binomial,
    control = glm.control(epsilon = 1e-4, maxit = 100, trace = FALSE))
  
  # f2 <- glm(
  #   y ~ 1 + pref + d1 + d2 + d3 , data = d, family = binomial,
  #   control = glm.control(epsilon = 1e-4, maxit = 100, trace = FALSE))
  
  # lapply(list(f1_1, f1_2, f1_3), function(z) summary(z)$coef)
  # dair
  d_dair <- copy(d)
  d_dair[, d1 := factor(1)]
  d_dair[, p_hat := predict(f1_1, newdata = d_dair, type = "response")]
  
  # averages over preference
  d_prop <- d[d1 == 1, .(prop = .N/nrow(d[d1==1])), keyby = pref]
  p_dair_ref <- d_prop[pref == 0, prop] * plogis(mu) + d_prop[pref == 1, prop] * plogis(mu + bp)
  p_dair_hat <- mean(d_dair$p_hat)
   
  # rev(1)
  d_rev_1 <- copy(d)
  d_rev_1[, d1 := factor(2)]
  d_rev_1[d2_entry == 1, d2 := factor(d[d2_entry == 1, d2_alloc + 2], levels = 1:3)]
  d_rev_1[d2_entry == 0, d2 := factor(1, levels = 1:3)]
  d_rev_1[, p_hat := predict(f1_2, newdata = d_rev_1, type = "response")]
  
  d_prop <- d[d1 == 2, .(prop = .N/nrow(d[d1==2])), keyby = d2]
  p_rev_1_ref <- d_prop[d2 == 1, prop] * plogis(mu + bd1[2]) + 
    d_prop[d2 == 2, prop] * plogis(mu + bd1[2] + bd2[2]) + 
    d_prop[d2 == 3, prop] * plogis(mu + bd1[2] + bd2[3])
  p_rev_1_hat <- mean(d_rev_1$p_hat)
   
  # # rev(2)
  d_rev_2 <- copy(d)
  d_rev_2[, d1 := factor(3)]
  d_rev_2[d3_entry == 1, d3 := factor(d[d3_entry == 1, d3_alloc + 2], levels = 1:3)]
  d_rev_2[d3_entry == 0, d3 := factor(1, levels = 1:3)]
  d_rev_2[, p_hat := predict(f1_3, newdata = d_rev_2, type = "response")]
  
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

d_fig <- d_g[, .(bias_rd, bias_p_dair, bias_p_rev, bias_p_rev_1, bias_p_rev_2)]
d_fig <- melt(d_fig, measure.vars = names(d_fig))

ggplot(d_fig, aes(x = value)) +
  geom_density() +
  geom_vline(data = d_fig[, .(mu = mean(value)), keyby = variable],
             aes(xintercept = mu), lwd = 0.3) + 
  facet_wrap(~variable)
#


