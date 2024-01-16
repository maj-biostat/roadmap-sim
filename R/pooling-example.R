library(data.table)
library(cmdstanr)
library(ggplot2)

N <- 2500
get_data_sim <- function(N = 2500){
  
  d <- data.table(id = 1:N)
  d[, silo := sample(1:3, .N, replace = T)]
  d[, x := rbinom(.N, 1, 0.5)]
  
  # silo specific treatment effects
  beta <- log(c(1.2, 0.6, 2))
  alpha <- qlogis(0.4)
  plogis(alpha + beta)
  
  # common baseline
  d[, eta := alpha]
  d[x == 1, eta := eta + beta[silo]]
  
  # sanity
  # d[, unique(eta), keyby = .(silo, x)]
  # alpha + beta
  
  d[, p := plogis(eta)]
  d[, y := rbinom(.N, 1, p)]
  
  # sanity
  # d[, .(mu = mean(y), p = unique(p)), keyby = .(silo, x)]
  
  d
}


main <- function(){
  
  
  m1 <- cmdstanr::cmdstan_model("./stan/pooling01.stan")
  m2 <- cmdstanr::cmdstan_model("./stan/pooling02.stan")
  m3 <- cmdstanr::cmdstan_model("./stan/pooling03.stan")
  
  d <- get_data_sim()
  
  ds <- d[, .(y = sum(y), n = .N), keyby = x]
  ld <- list(
    N = nrow(ds),
    y = ds$y,
    n = ds$n,
    x = ds$x
  )
  
  f1 <- m1$sample(ld,
                  iter_warmup = 1000,
                  iter_sampling = 2000,
                  parallel_chains = 2,
                  chains = 2,
                  refresh = 100, show_exceptions = F)
  
  ds <- d[, .(y = sum(y), n = .N), keyby = .(x, silo)]
  ld <- list(
    N = nrow(ds),
    y = ds$y,
    n = ds$n,
    x = ds$x,
    s = ds$silo
  )
  
  f2 <- m2$sample(ld,
                  iter_warmup = 1000,
                  iter_sampling = 2000,
                  parallel_chains = 2,
                  chains = 2,
                  refresh = 100, show_exceptions = F)
  f3 <- m3$sample(ld,
                  iter_warmup = 1000,
                  iter_sampling = 2000,
                  parallel_chains = 2,
                  chains = 2,
                  refresh = 100, show_exceptions = F)

  f1  
  f2
  f3
  
  # treatment effects log-or
  
  post1 <- data.table(f1$draws(variables = c("b"), format = "matrix"))
  post2 <- data.table(f2$draws(variables = c("b"), format = "matrix"))
  post3 <- data.table(f3$draws(variables = c("b"), format = "matrix"))
  
  d_tbl1 <- melt(post1, measure.vars = names(post1))
  d_tbl1[, value := exp(value)]
 
  d_tbl2 <- melt(post2, measure.vars = names(post2))
  d_tbl2[, value := exp(value)]
  d_tbl2[, silo := gsub("b[", "", variable, fixed = T)]
  d_tbl2[, silo := as.integer(gsub("]", "", silo, fixed = F))]
  
  d_tbl3 <- melt(post3, measure.vars = names(post2))
  d_tbl3[, value := exp(value)]
  d_tbl3[, silo := gsub("b[", "", variable, fixed = T)]
  d_tbl3[, silo := as.integer(gsub("]", "", silo, fixed = F))]
  

  
  d_fig1 <- d_tbl1[, .(mu = mean(value), 
                       q_025 = quantile(value, prob = 0.025),
                       q_975 = quantile(value, prob = 0.975)
  )]
  
  d_fig2 <- d_tbl2[, .(mu = mean(value), 
                       q_025 = quantile(value, prob = 0.025),
                       q_975 = quantile(value, prob = 0.975)
  ), keyby = .(silo)]
  
  d_fig3 <- d_tbl3[, .(mu = mean(value), 
                       q_025 = quantile(value, prob = 0.025),
                       q_975 = quantile(value, prob = 0.975)
  ), keyby = .(silo)]

  ggplot(d_fig1) + 
    geom_hline(aes(yintercept = mu)) +
    # unpooled
    geom_point(data = d_fig2, aes(x = silo, y = mu)) +
    # partial pooled
    geom_point(data = d_fig3, aes(x = silo, y = mu), col = 2)
  
  
  # pooled view of the average effect...
  
  d_fig1 <- melt(post1, measure.vars = names(post1))
  
  post3 <- data.table(f3$draws(variables = c("b_avg"), format = "matrix"))
  d_fig3 <- melt(post3, measure.vars = names(post3))
  d_fig3[, variable := "b"]
  
  ggplot(d_fig1) + 
    geom_density(aes(x = value)) +
    geom_density(data = d_fig3, aes(x = value), lty = 2)
  
  # probability
  
  post1 <- data.table(f1$draws(variables = c("p"), format = "matrix"))
  post2 <- data.table(f2$draws(variables = c("p"), format = "matrix"))
  post3 <- data.table(f3$draws(variables = c("p"), format = "matrix"))
  
  d_tbl1 <- melt(post1, measure.vars = names(post1))
  d_tbl1[variable == "p[1]", x := 0]
  d_tbl1[variable == "p[2]", x := 1]
  
  d_tbl2 <- melt(post2, measure.vars = names(post2))
  d_tbl2[variable %like% ",1]", x := 0]
  d_tbl2[variable %like% ",2]", x := 1]
  d_tbl2[, silo := gsub("p[", "", variable, fixed = T)]
  d_tbl2[, silo := as.integer(gsub(",.]", "", silo, fixed = F))]
  
  d_tbl3 <- melt(post3, measure.vars = names(post2))
  d_tbl3[variable %like% ",1]", x := 0]
  d_tbl3[variable %like% ",2]", x := 1]
  d_tbl3[, silo := gsub("p[", "", variable, fixed = T)]
  d_tbl3[, silo := as.integer(gsub(",.]", "", silo, fixed = F))]
  
  # d_tbl[, y := gsub("^p\\[.*,", "", variable)]
  
  d_fig1 <- d_tbl1[, .(mu = mean(value), 
             q_025 = quantile(value, prob = 0.025),
             q_975 = quantile(value, prob = 0.975)
             ), keyby = x]
  
  d_fig2 <- d_tbl2[, .(mu = mean(value), 
             q_025 = quantile(value, prob = 0.025),
             q_975 = quantile(value, prob = 0.975)
  ), keyby = .(x, silo)]
  
  d_fig3 <- d_tbl3[, .(mu = mean(value), 
             q_025 = quantile(value, prob = 0.025),
             q_975 = quantile(value, prob = 0.975)
  ), keyby = .(x, silo)]
  
  
  
}

