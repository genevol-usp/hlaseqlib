#' Process HLA alignment data in *.nuc files and return a data.frame.
#'
#' @param locus character string.
#' @param imgt.database character string.
#' @param short logical. Whether to create shorter version of index.
#'
#' @return A data.frame.
#'
#' @export

hla_compile_index <- function(locus, imgt.database, short = FALSE) {

    message(paste("Processing locus", locus)) 

    # Process alignments
    hla_df <- hla_read_alignment(locus, imgt.database)

    if (all(grepl("\\*", hla_df$cds))) {
	
	message(paste0("No complete sequence for locus ", locus, "; returning NA"))
	return(NA)

    } else if (nrow(hla_df) == 1L || all(!grepl("\\*", hla_df$cds))) {
	
	# If N_alleles = 1 or all alleles are incomplete, output:
	final_df <- hla_df

    } else {

	# find closest complete allele for each incomplete allele
	distmatrix <- make_dist_matrix(hla_df) 

        closest_allele_df <- make_closest_allele_df(distmatrix)

	closest_allele_df$id <- closest_allele_df %>% 
	    dplyr::group_indices(inc_allele)

	closest_allele_df_step2 <-
	    tidyr::gather(closest_allele_df, ix, allele, 1:2) %>%
	    dplyr::distinct(id, allele) %>%
	    dplyr::left_join(hla_df, by = "allele") %>%
	    split(.$id) %>%
	    purrr::map(make_dist_matrix) %>%
	    purrr::map(make_closest_allele_df) %>%
	    dplyr::bind_rows()

	closest_within_type <- closest_allele_df_step2 %>%
	    dplyr::mutate(`1` = hla_trimnames(inc_allele, 1) == hla_trimnames(closest, 1),
			  `2` = hla_trimnames(inc_allele, 2) == hla_trimnames(closest, 2),
			  `3` = hla_trimnames(inc_allele, 3) == hla_trimnames(closest, 3),
			  `4` = hla_trimnames(inc_allele, 4) == hla_trimnames(closest, 4)) %>%
	    tidyr::gather(field, value, `1`:`4`) %>%
	    dplyr::group_by(inc_allele) %>%
	    dplyr::filter(all(value == FALSE) | value == TRUE) %>%
	    dplyr::slice(which.max(field)) %>%
	    dplyr::ungroup() %>%
	    dplyr::arrange(inc_allele) %>%
	    dplyr::select(inc_allele, closest)

	inferred_df <- closest_within_type %>%
	    dplyr::left_join(hla_df, by = c("inc_allele" = "allele")) %>%
	    dplyr::left_join(hla_df, by = c("closest" = "allele")) %>%
	    dplyr::mutate(cds = purrr::map2_chr(cds.x, cds.y, hla_attribute_seq)) %>%
	    dplyr::select(allele = inc_allele, cds)

	final_df <- hla_df %>%
	    dplyr::filter(!grepl("\\*", cds)) %>%
	    dplyr::bind_rows(inferred_df) %>%
	    dplyr::arrange(allele)
    }

    final_df %>%
	dplyr::mutate(cds = hla_format_sequence(cds)) %>%
	dplyr::mutate(allele3f = hla_trimnames(allele, 3)) %>%
	dplyr::distinct(allele3f, cds, .keep_all = TRUE) %>%
	dplyr::group_by(allele3f) %>%
	dplyr::mutate(n = dplyr::n()) %>%
	dplyr::ungroup() %>%
	dplyr::mutate(allele = ifelse(n > 1L, allele, allele3f)) %>%
	dplyr::select(allele, cds) %>%
	dplyr::arrange(allele)
}
