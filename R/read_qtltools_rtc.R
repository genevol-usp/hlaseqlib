#' Read in QTLtools rtc output
#'
#' @param rtc_output character string. Path to file.
#'
#' @return data.frame
#'
#' @export

read_qtltools_rtc <- function(rtc_output)
  readr::read_delim(rtc_output, col_names = FALSE, delim = " ",
		    col_types = "cccciiiiiiiiiiiiiiiddd") %>%
  purrr::set_names(c("gwas_var", "qtl_var", "gene", "gene_grp", "gwas_chr",
		     "gwas_pos", "gwas_rank", "qtl_chr", "qtl_pos", "qtl_rank",
		     "gene_chr", "gene_pos", "dist_gwas_qtl", "dist_gwas_phen",
		     "gwas_var_index", "qtl_var_index", "region_start",
		     "region_end", "n_vars", "rtc", "d_prime", "r_squared"))

