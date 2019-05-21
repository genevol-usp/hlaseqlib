#' Infer HLA genotypes from RNA-seq quantifications.
#'
#' @param hla_quants data.frame.
#' @param th numeric. Threshold for expression 2nd/1st allele.
#'
#' @return A data.frame.
#'
#' @export

hla_genotype <- function(hla_quants, th = 0.05) {

    hla_quants$i <- hla_quants %>% dplyr::group_indices(locus, tpm, counts)

    hla_quants %>%
	dplyr::group_by(locus, i) %>%
	dplyr::summarise(allele = paste(allele, collapse = "-"),
		  tpm = sum(tpm),
		  counts = sum(counts)) %>%
	dplyr::group_by(locus) %>%
	dplyr::top_n(2, counts) %>%
	dplyr::filter(max(counts) == 0 | counts/max(counts) > th) %>%
	dplyr::mutate(i = ifelse(n() == 1, "1_1", "1"),
		      tpm = ifelse(i == "1_1", tpm/2L, tpm),
		      counts = ifelse(i == "1_1", counts/2L, counts)) %>%
	dplyr::ungroup() %>%
	tidyr::separate_rows(i, sep = "_") %>%
	dplyr::select(-i) %>%
	dplyr::mutate(allele = ifelse(counts == 0, NA, allele))
}
