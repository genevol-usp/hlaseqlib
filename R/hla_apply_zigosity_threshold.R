#' Apply zigosity threshold to expression estimates  
#'
#' @param hla_quants data.frame.
#' @param th numeric.
#'
#' @return A data.frame
#'
#' @export

hla_apply_zigosity_threshold <- function(hla_quants, th = 0.1) {

    hla_quants <- hla_quants %>%
	dplyr::mutate(subject = as.character(subject)) %>%
	dplyr::group_by(subject, locus) %>%
	dplyr::mutate(i = as.integer(any(est_counts/max(est_counts) <= th))) %>%
	dplyr::ungroup()

    df_0 <- dplyr::filter(hla_quants, i == 0L) %>% dplyr::select(-i)

    df_1 <- hla_quants %>%
	dplyr::filter(i == 1L) %>%
	dplyr::group_by(subject, locus) %>%
	dplyr::mutate(m = as.integer(est_counts == max(est_counts)),
		      est_counts = sum(est_counts)/2L, tpm = sum(tpm)/2L) %>%
	dplyr::ungroup() %>%
	dplyr::filter(m == 1L) %>%
	dplyr::select(-i, -m)

    out <- dplyr::bind_rows(df_0, df_1, df_1) %>%
	dplyr::arrange(subject, locus, allele)

    out
}

