#' Pfase HLA alleles given distances to 1000G haplotypes
#'
#' @param x data.frame. A data.frame with distances of alleles to each 1000G
#' haplotype.
#'
#' @return A data.frame with 2 rows with inferred haplotypes for each
#' individual.
#'
#' @export


phase_hla <- function(x) {

    x <- dplyr::mutate(x, locus = sub("^HLA-", "", locus))

    loci <- dplyr::syms(unique(x$locus))

    combs <- x %>%
	dplyr::group_by(locus) %>%
	dplyr::mutate(i = seq_len(n())) %>%
	dplyr::select(-diffs) %>%
	tidyr::spread(locus, allele) %>%
	dplyr::select(-subject, -i) %>%
	tidyr::expand(hap, !!! loci) %>%
	tidyr::drop_na()

    combs_1 <- combs %>%
	dplyr::filter(hap == 1) %>%
	dplyr::mutate(unique_hap_id = seq_len(n()))

    combs_2 <- combs %>%
	dplyr::filter(hap == 2) %>%
	dplyr::mutate(unique_hap_id = rep(list(combs_1$unique_hap_id), nrow(.))) %>%
	tidyr::unnest()

    min_hap <- 
        dplyr::left_join(combs_1, combs_2, by = "unique_hap_id", suffix = c(".1", ".2")) %>%
	dplyr::select(-hap.1, -hap.2, -unique_hap_id) %>%
	dplyr::mutate(geno_id = seq_len(n())) %>%
	tidyr::gather(locus, allele, -geno_id) %>%
	dplyr::group_by(geno_id) %>%
	dplyr::filter(all(unique(x$allele) %in% allele)) %>%
	tidyr::separate(locus, c("locus", "hap"), sep = "\\.", convert = TRUE) %>%
	dplyr::left_join(x, by = c("locus", "hap", "allele")) %>%
	dplyr::group_by(geno_id) %>%
	dplyr::mutate(total_diff = sum(diffs)) %>%
	dplyr::ungroup() %>%
	dplyr::filter(total_diff == min(total_diff)) %>%
	dplyr::group_by(hap, locus) %>%
	dplyr::summarise(allele = paste(unique(allele), collapse = "/")) %>%
	dplyr::ungroup() %>%
	dplyr::select(locus, hap, allele) %>%
	dplyr::arrange(locus, hap) %>%
	dplyr::mutate(locus = paste0("HLA-", locus))

    min_hap
}
