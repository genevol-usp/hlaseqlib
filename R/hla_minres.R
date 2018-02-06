#' Reduce a vector of ambiguous alleles to the maximum resolution that includes
#' all of them.
#'
#' @param allele_vec character vector.
#'
#' @return character vector. A single allele at the resolution that includes all
#' alleles in the original vector, \code{NULL} when there is no resolution that
#' includes all alleles.
#'
#' @export

hla_minres <- function(allele_vec) {

    out <- 
	lapply(1:4, function(i) hla_trimnames(allele_vec, i) %>% unique()) %>%
	purrr::keep(~length(.x) == 1) %>%
	dplyr::last()

    if (is.null(out)) return(NA)

    out
}
