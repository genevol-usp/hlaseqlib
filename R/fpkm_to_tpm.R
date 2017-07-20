#' Convert FPKM to TPM
#'
#' @param fpkm numeric.
#'
#' @return A numeric vector of TPM values.
#' @export

fpkm_to_tpm <- function(fpkm) 
  exp(log(fpkm) - log(sum(fpkm)) + log(1e6))
