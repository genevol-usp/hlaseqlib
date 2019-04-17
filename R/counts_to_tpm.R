#' Convert read counts to TPM
#'
#' @param counts numeric.
#' @param eff_len numeric.
#'
#' @return A numeric vector of TPM values.
#' @export

counts_to_tpm <- function(counts, eff_len) {
  rate <- log(counts) - log(eff_len)
  denom <- log(sum(exp(rate)))
  exp(rate - denom + log(1e6))
}
