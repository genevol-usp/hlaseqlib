#' Process HLA alignment data in *.gen files and return a data.frame.
#'
#' @param gen_file character string. Path to IMGT alignment file.
#' @param omit_suffix character string. String of optional suffixes indicating
#' expression status. Possible values are NULL (keep all allele names with
#' suffixes) or a string of suffixes to omit (e.g., "NQL"). Default is "N".
#'
#' @return A data.frame.
#'
#' @export

hla_read_utr <- function(gen_file, omit_suffix = "N") {
  
    locus <- sub("_gen\\.txt", "", basename(gen_file))

    gen <- readLines(gen_file) %>%
	gsub("\\s{2,}", " ", .) %>%
	trimws %>%
	.[grepl(sprintf("^%s\\d?\\*\\d{2,3}[:A-Z0-9]*\\s", locus), .)]

    hla_df <-
	tibble::tibble(allele = gsub("^(\\S+)\\s(.*)$", "\\1", gen),
		       cds = gsub("\\s", "", gsub("^(\\S+)\\s(.*)$", "\\2", gen))) %>%
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

    utr_df <- hla_df %>%
	dplyr::mutate(cds = strsplit(cds, "\\|")) %>%
	tidyr::unnest() %>%
	dplyr::group_by(allele) %>%
	dplyr::mutate(exon = seq_len(n())) %>%
	dplyr::slice(c(which.min(exon), which.max(exon))) %>%
	dplyr::ungroup() %>%
	dplyr::select(allele, exon, cds) %>%
	tidyr::spread(exon, cds) %>%
	setNames(c("allele", "utr5", "utr3")) %>%
	dplyr::mutate(utr5 = substring(utr5, nchar(utr5)-99L, nchar(utr5)),
		      utr3 = substring(utr3, 1, 200))

    utr_df
}
