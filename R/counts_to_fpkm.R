#' Convert read counts to FPKM
#'
#' @param counts numeric.
#' @param eff_len numeric.
#'
#' @return A numeric vector of FPKM values.
#' @export

counts_to_fpkm <- function(counts, eff_len) {
  N <- sum(counts)
  exp(log(counts) + log(1e9) - log(eff_len) - log(N))
}
