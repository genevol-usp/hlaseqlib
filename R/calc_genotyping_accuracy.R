#' Calculate genotyping accuracy from RNA-seq calls.
#' 
#' Calculates genotyping accuracies from RNA-seq calls using a gold standard.
#' 
#' @param data_x data.frame.
#' @param data_y data.frame. Gold standard data.
#' @param by_locus logical. Whether or not to summarize by locus.
#' 
#' @return A data.frame.
#' @export

calc_genotyping_accuracy <- function(data_x, data_y, by_locus = TRUE) {
 
    format_hla <- function(df, common) {
	dplyr::inner_join(df, common, by = c("subject", "locus")) %>%
	dplyr::arrange(subject, locus, allele) %>%
	dplyr::group_by(subject, locus) %>%
	dplyr::mutate(h = 1:n()) %>%
	dplyr::ungroup() %>%
	tidyr::separate_rows(allele, sep = "/") %>%
	tidyr::separate_rows(allele, sep = "=") %>%
	dplyr::mutate(supertype = hla_trimnames(allele, 1),
		      twofield = hla_trimnames(allele, 2),
		      threefield = hla_trimnames(allele, 3))
    }

    common <- 
	dplyr::inner_join(dplyr::select(data_x, subject, locus),
			  dplyr::select(data_y, subject, locus), 
			  by = c("subject", "locus")) %>%
	dplyr::distinct()

    data_x <- format_hla(data_x, common)
    data_y <- format_hla(data_y, common)

    x <-
	dplyr::inner_join(data_x, data_y, by = c("subject", "locus")) %>%
	dplyr::mutate(correct = allele.x == allele.y |
				allele.x == threefield.y |
				threefield.x == allele.y |
				twofield.x == allele.y | 
				supertype.x == allele.y)

    right <- 
	dplyr::filter(x, correct) %>%
	dplyr::distinct(subject, locus, h.y, .keep_all = TRUE) %>%
	dplyr::distinct(subject, locus, h.x, .keep_all = TRUE)

    wrong <-
	dplyr::anti_join(x, right, by = c("subject", "locus", "h.y")) %>%
	dplyr::anti_join(right, by = c("subject", "locus", "h.x"))

    w_same2 <-
	dplyr::filter(wrong, twofield.x == twofield.y) %>%
	dplyr::distinct(subject, locus,  allele.y, .keep_all = TRUE) %>%
	dplyr::distinct(subject, locus, h.x, .keep_all = TRUE)

    w_samet <- wrong %>%
	dplyr::filter(twofield.x != twofield.y, supertype.x == supertype.y) %>%
	dplyr::distinct(subject, locus,  allele.y, .keep_all = TRUE) %>%
	dplyr::distinct(subject, locus, h.x, .keep_all = TRUE) %>%
	dplyr::anti_join(w_same2, by = c("subject", "locus", "allele.x"))

    w_same <- dplyr::bind_rows(w_same2, w_samet)

    w_difft <-
	dplyr::anti_join(wrong, w_same, by = c("subject", "locus", "h.y")) %>%
	dplyr::anti_join(w_same, by = c("subject", "locus", "h.x"))

    df_final <- 
	dplyr::bind_rows(right, w_same2, w_samet, w_difft) %>% 
	dplyr::group_by(subject, locus, h.y) %>%
	dplyr::summarize(allele.x = paste(unique(allele.x), collapse = "/"),
			 allele.y = paste(unique(allele.y), collapse = "/"),
			 correct = unique(correct)) %>%
	dplyr::ungroup() %>%
	dplyr::select(subject, locus, allele.x, allele.y, correct) %>%
	tidyr::separate_rows(allele.x, sep = "/")

    if (by_locus) {
	df_final %>%
	    dplyr::group_by(locus) %>%
	    dplyr::summarize(accuracy = mean(correct)) %>%
	    dplyr::ungroup()
    } else {
	df_final
    }
}
