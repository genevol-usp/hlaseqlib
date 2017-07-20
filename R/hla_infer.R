#' Infer sequences missing on the IMGT alignments.
#'
#' @param hla_df data.frame. A data.frame such as the one returned by 
#' \code{hla_dataframe}.
#' @param .parallel logical. Whether to run in parallel. If \code{TRUE}, provide
#' a parallel backend before calling function.
#'
#' @return A data.frame with the original sequences of the complete alleles
#' plus the inferred sequences for incomplete alleles.
#'
#' @export

hla_infer <- function(hla_df, .parallel = TRUE) {

  incomplete_df <- dplyr::filter(hla_df, grepl("\\*", cds))
  
  complete_df <- dplyr::filter(hla_df, ! allele %in% incomplete_df$allele)

  if (is.data.frame(incomplete_df) && nrow(incomplete_df) == 0) {
    return(complete_df)
  }

  complete_distinct <- dplyr::distinct(complete_df, cds, .keep_all = TRUE)

  inferred_df <-
    incomplete_df %>%
    plyr::ddply(~allele, . %>%
		dplyr::mutate(closest = purrr::map_chr(cds, find_closest_allele, 
						       complete_df = complete_distinct)) %>%
		dplyr::select(allele, closest, cds) %>%
		tidyr::separate_rows(closest, sep = "/") %>%
		dplyr::left_join(complete_distinct, c("closest" = "allele")) %>%
		dplyr::mutate(inferred_seq = purrr::map2_chr(cds.x, cds.y, 
							     hla_attribute_seq)) %>%
		dplyr::select(locus, allele, cds = inferred_seq),
	      .parallel = .parallel) %>%
    tibble::as_tibble() %>%
    dplyr::distinct() %>%
    dplyr::select(locus, allele = allele_id, cds)

  dplyr::bind_rows(complete_df, inferred_df)
}

