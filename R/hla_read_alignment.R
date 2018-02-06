#' Process HLA alignment data in *.nuc files and return a data.frame.
#'
#' @param nuc_file character string. Path to IMGT alignment file.
#' @param exons integer. Range of vectors to select (default = all exons).
#' @param by_exon logical. Whether to return data separated by exons.
#' @param omit_suffix character string. String of optional suffixes indicating
#' expression status. Possible values are NULL (keep all allele names with
#' suffixes) or a string of suffixes to omit (e.g., "NQL"). Default is "N".
#' @param rm_incomplete logical. Whether to remove incomplete cdss.
#' @param keep_sep logical. Whether to keep the exon separator character "|" in
#' the cdss.
#'
#' @return A data.frame.
#'
#' @export

hla_read_alignment <- function(nuc_file, exons = NULL, by_exon = FALSE, 
			       omit_suffix = "N", rm_incomplete = FALSE, 
			       keep_sep = FALSE) {
  
    locus <- sub("_nuc\\.txt", "", basename(nuc_file))
    
    nuc <-
	readLines(nuc_file) %>%
	gsub("\\s{2,}", " ", .) %>%
	trimws %>%
	.[grepl(sprintf("^%s\\d?\\*\\d{2,3}[:A-Z0-9]*\\s", locus), .)]

    hla_df <-
	tibble::tibble(allele = gsub("^(\\S+)\\s(.*)$", "\\1", nuc),
		       cds = gsub("\\s", "", gsub("^(\\S+)\\s(.*)$", "\\2", nuc))) %>%
	tibble::rownames_to_column() %>% 
	dplyr::group_by(allele) %>%
	dplyr::summarize(rowname = min(as.integer(rowname)), 
			 cds = paste(cds, collapse = "")) %>%
	dplyr::arrange(rowname) %>%
	dplyr::select(-rowname)
  
    if (!is.null(omit_suffix)) {
	suffix <- sprintf("[^%s]$", omit_suffix)
	hla_df <- dplyr::filter(hla_df, grepl(suffix, allele))
    }

    hla_df$cds <- stringr::str_split(hla_df$cds, "", simplify = TRUE) %>%
	.[, apply(., 2, function(x) !all(x == "." | x == "*"))] %>%
	apply(2, function(i) {i[i == "-"] <- i[1]; i}) %>%
	apply(1, . %>% paste(collapse = ""))

    if (rm_incomplete) {
	hla_df <- dplyr::filter(hla_df, !grepl("\\*", cds))
	if (nrow(hla_df) == 0) return(hla_df)
    }
  
    if (keep_sep && is.null(exons)) return(hla_df)

    hla_df <- hla_df %>%
	dplyr::mutate(cds = strsplit(cds, "\\|")) %>%
	tidyr::unnest() %>%
	dplyr::group_by(allele) %>%
	dplyr::mutate(exon = seq_len(n())) %>%
	dplyr::select(allele, exon, cds) 

    if (!is.null(exons) && is.numeric(exons)) {
	hla_df <- dplyr::filter(hla_df, exon %in% exons)
    }

    if (by_exon) {
	tidyr::unite(hla_df, allele, -cds)
    } else if (!is.null(exons) && is.numeric(exons) && keep_sep) {
	dplyr::summarize(hla_df, cds = paste(cds, collapse = "|")) 
    } else {
	dplyr::summarise(hla_df, cds = paste(cds, collapse = "")) 
    }
}
