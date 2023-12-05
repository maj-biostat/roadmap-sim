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
