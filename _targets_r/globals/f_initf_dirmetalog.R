initf_dirmetalog <- function() {
  deltas3 <- c(0.25, 0.25, 0.15)+runif(3, -0.02, 0.02)
  list(delta = c(deltas3, 1-sum(deltas3)))
}
