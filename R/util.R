# depends on init.R

crossjoin = function(d1, d2) {
  stopifnot(is.data.table(d1), is.data.table(d2))
  d1[, .tempCol := 1]
  d2[, .tempCol := 1]
  dd <- d1[d2, on = ".tempCol", , allow.cartesian = T]
  dd[, .tempCol := NULL]
  setcolorder(dd, names(d1)[-ncol(d1)])
  dd
}

antidiag <- function(X, offset = 0L) {
  X[col(X) + row(X) - ncol(X) - 1L == offset]
}

post_smry <- function(v){
  mu <- mean(v)
  sd <- sd(v)
  q_025 <- quantile(v, prob = 0.025)
  q_975 <- quantile(v, prob = 0.975)
  sprintf("%.3f %.3f (%.3f, %.3f)", mu, sd, q_025, q_975)
}

yn_rand <- function(N, pr_y = 0.5){
  pr_n <- 1 - pr_y
  sample(c("Y", "N"), N, replace = T, prob = c(pr_y, pr_n))
}

tbl_ex_trial <- function(d){
  
  d_b <- d[, .(y = sum(y), n = .N, p = plogis(unique(eta))), keyby = .(silo, joint, ea, a, qa, eb, b, ec, c)]
  d_b[, p_mle := y/n]
  d_b[, silo := factor(silo, levels = c("early", "late", "chronic"))]
  d_b[, joint := factor(joint, levels = c("knee", "hip"))]
  setcolorder(d_b, c("silo", "joint", "ea", "a", "qa", "eb", "b", "ec", "c", "p_mle", "p"))

  gt_tbl <- d_b[order(silo, joint)] |> 
    gt(
      groupname_col = c("silo", "joint")
    ) |> 
    cols_align(
      align = "left",
      columns = starts_with(c("silo", "joint"))
    ) |> 
    cols_align(
      align = "center",
      columns = starts_with(c("ea", "a", "qa", "eb", "b", "ec", "c"))
    ) |> 
    fmt_number(
      columns = starts_with(c("p")),
      decimals = 2
    ) |>
    tab_spanner(
      label = html("Surgical D<sub>a</sub>"),
      columns = c(ea, a, qa)
    ) |>
    tab_spanner(
      label = html("Duration D<sub>b</sub>"),
      columns = c(eb, b)
    ) |>
    tab_spanner(
      label = html("Type D<sub>c</sub>"),
      columns = c(ec, c)
    ) |>
    tab_spanner(
      label = html("Response"),
      columns = c(y, n, p_mle, p)
    ) |>
    cols_label(
      silo = html("Silo"),
      joint = html("Infection <br>site"),
      ea = html("member"),
      a = html("assigned"),
      qa = html("plan"),
      eb = html("member"),
      b = html("assigned"),
      ec = html("member"),
      c = html("assigned"),
      y = html("y"),
      n = html("n"),
      p_mle = html("MLE (p<sub>y</sub>)"),
      p = html("TRUE (p<sub>y</sub>)")
    ) |>
    summary_rows(
      groups = everything(),
      columns = c("y", "n"),
      fns = list(
        subtotal = ~ sum(., na.rm = TRUE)
      )
    ) |>
    summary_rows(
      groups = everything(),
      columns = c("p_mle",),
      fns = list(subtotal = ~ (round(sum(y) / sum(n), 2)))
    ) |>
    grand_summary_rows(
      columns = c("y", "n"),
      fns = list(
        total = ~ sum(., na.rm = TRUE)
      )
    ) |>
    grand_summary_rows(
      columns = c("p_mle",),
      fns = list(total = ~ (round(sum(y) / sum(n), 2)))
    ) |>
    tab_options(
      table.font.size = "80%"
    ) |>
    tab_footnote(
      footnote = "Transformed from the log-odds of response as used in the linear predictor to simulate data.",
      locations = cells_column_labels(columns = p)
    ) |>
    tab_footnote(
      footnote = "w06p1 = 6 weeks following one-stage procedure, w12p1 = 12 weeks following one-stage procedure etc",
      locations = cells_column_labels(columns = b)
    )
  gt_tbl
}


post_alpha <- function(post){
  
  cols <- names(post)[names(post) %like% "alpha"]
  
  alpha <- melt(post[, .SD, .SDcols = cols], measure.vars = cols)
  alpha[, idx := gsub(".*\\[", "", variable)]
  alpha[, idx := gsub("\\]", "", idx)]
  alpha[, idx := as.integer(idx)]
  
  d_fig <- alpha[, .(a_med = median(value), 
                     a_q025 = quantile(value, prob = 0.025),
                     a_q975 = quantile(value, prob = 0.975)), keyby = idx]
  
  d_fig <- merge(d_fig, sim_spec$i_a_s_u, by.x = "idx", by.y = "su")
  d_fig[, grp := paste0(silo, " (", joint, ")")]
  d_fig[, grp := factor(grp, levels = d_fig$grp)]
  
  d_fig[, a := sim_spec$a_s_u[cbind(silo, joint)]]
  d_fig
}


post_dom_a <- function(post){
  
  cols <- names(post)[names(post) %like% "b_a"]
  
  dom_a <- melt(post[, .SD, .SDcols = cols], measure.vars = cols)
  dom_a[, idx := gsub(".*\\[", "", variable)]
  dom_a[, idx := gsub("\\]", "", idx)]
  dom_a[, idx := as.integer(idx)]
  dom_a[, variable := gsub("\\[.*\\]", "", variable)]
  
  d_fig <- dom_a[, .(
    b_med = median(value), 
    b_q025 = quantile(value, prob = 0.025),
    b_q975 = quantile(value, prob = 0.975)), 
    keyby = .(variable, idx)]
  
  d_fig[variable == "b_a_l", `:=`(
    silo = "late", 
    b = sim_spec$b_a_late[idx],
    trt = names(sim_spec$b_a_late)[idx]
  )]
  
  d_fig[variable == "b_a_c", `:=`(
    silo = "chronic", 
    b = sim_spec$b_a_chronic[idx],
    trt = names(sim_spec$b_a_chronic)[idx]
  )]
  
  d_fig[, silo := factor(silo, levels = c("late", "chronic"))]
  d_fig
}


post_dom_b <- function(post){
  
  cols <- names(post)[names(post) %like% "b_b"]
  
  dom_b <- melt(post[, .SD, .SDcols = cols], measure.vars = cols)
  dom_b[, idx := gsub(".*\\[", "", variable)]
  dom_b[, idx := gsub("\\]", "", idx)]
  dom_b[, idx := as.integer(idx)]
  dom_b[, variable := gsub("\\[.*\\]", "", variable)]
  dom_b <- dom_b[idx < 3]
  
  d_fig <- dom_b[, .(
    b_med = median(value), 
    b_q025 = quantile(value, prob = 0.025),
    b_q975 = quantile(value, prob = 0.975)), 
    keyby = .(variable, idx)]
  
  d_fig[variable == "b_b1_l", `:=`(
    silo = "late", 
    qa = "one",
    b = sim_spec$b_b1_late_one[idx],
    trt = names(sim_spec$b_b1_late_one)[idx]
  )]
  
  d_fig[variable == "b_b2_l", `:=`(
    silo = "late",
    qa = "two",
    b = sim_spec$b_b2_late_two[idx],
    trt = names(sim_spec$b_b2_late_two)[idx]
  )]
  
  d_fig[variable == "b_b1_c", `:=`(
    silo = "chronic",
    # plan is same as assigned
    qa = "one",
    b = sim_spec$b_b1_chronic_one[idx],
    trt = names(sim_spec$b_b1_chronic_one)[idx]
  )]
  
  d_fig[variable == "b_b2_c", `:=`(
    silo = "chronic",
    qa = "two",
    b = sim_spec$b_b2_chronic_two[idx],
    trt = names(sim_spec$b_b2_chronic_two)[idx]
  )]
  
  d_fig[, silo := factor(silo, levels = c("late", "chronic"))]
  d_fig
}


post_dom_c <- function(post){
  
  cols <- names(post)[names(post) %like% "b_c"]
  
  dom_c <- melt(post[, .SD, .SDcols = cols], measure.vars = cols)
  dom_c[, idx := gsub(".*\\[", "", variable)]
  dom_c[, idx := gsub("\\]", "", idx)]
  dom_c[, idx := as.integer(idx)]
  dom_c[, variable := gsub("\\[.*\\]", "", variable)]
  dom_c <- dom_c[idx < 3]
  
  d_fig <- dom_c[, .(
    b_med = median(value), 
    b_q025 = quantile(value, prob = 0.025),
    b_q975 = quantile(value, prob = 0.975)), 
    keyby = .(variable, idx)]
  
  d_fig[variable == "b_c", `:=`(
    b = sim_spec$b_c[idx],
    trt = names(sim_spec$b_c)[idx]
  )]
  
  d_fig
}