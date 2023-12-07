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