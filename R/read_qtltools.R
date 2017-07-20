#' Read in QTLtools output from secondary chunk analysis
#'
#' @param qtltools_output character string. Path to file.
#'
#' @return data.frame
#'
#' @export

read_qtltools <- function(qtltools_output)
  readr::read_delim(qtltools_output, col_names = FALSE, delim = " ",
		    col_types = "cciiciicciiiddiiddii", progress = FALSE) %>%
  purrr::set_names(c("phen_id", "phen_chr", "phen_from", "phen_to", "strand", 
                     "n_vars", "dist", "var_id", "var_chr", "var_from", 
                     "var_to", "rank", "fwd_pval", "fwd_slope", "fwd_best",
                     "fwd_signif", "bwd_pval", "bwd_slope", "bwd_best", 
                     "bwd_signif"))