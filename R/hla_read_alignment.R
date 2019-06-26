#' Process HLA alignment data in *.nuc files and return a data.frame.
#'
#' @param locus character string. The HLA locus to be processed (e.g., "A", "DRB1").
#' @param imgtdb character string. Path to the IMGTHLA directory.
#' @param exons integer. Range of exons (default is all exons).
#' @param by_exon logical. Whether to return sequences separated for each exon (default FALSE).
#'
#' @return A data.frame.
#'
#' @examples
#' \dontrun{
#' hla_read_alignment(locus = "DRB1", imgtdb = "~/IMGTHLA") 
#' }
#' \dontrun{
#' hla_read_alignment(locus = "A", imgtdb = "~/IMGTHLA", exons = c(2, 3))
#' }
#'
#' @export

hla_read_alignment <- function(locus, imgtdb, exons = NULL, by_exon = FALSE) {
 
    locus_nuc <- ifelse(grepl("DRB\\d", locus), "DRB", locus) 
    
    nuc_file <- file.path(imgtdb, "alignments", paste0(locus_nuc, "_nuc.txt"))

    alignments <- readLines(nuc_file) %>%
	gsub("\\s{2,}", " ", .) %>%
	trimws %>%
	.[grepl(sprintf("^%s\\d?\\*\\d{2,3}[:A-Z0-9]*\\s", locus), .)]

    hla_df <-
	tibble::tibble(allele = gsub("^(\\S+)\\s(.*)$", "\\1", alignments),
		       cds = gsub("\\s", "", gsub("^(\\S+)\\s(.*)$", "\\2", alignments))) %>%
	tibble::rowid_to_column() %>% 
	dplyr::group_by(allele) %>%
	dplyr::summarize(rowid = min(rowid), 
			 cds = paste(cds, collapse = "")) %>%
	dplyr::arrange(rowid) %>%
	dplyr::select(-rowid)

    hla_df$cds <- stringr::str_split(hla_df$cds, "", simplify = TRUE) %>%
	apply(2, function(i) {i[i == "-"] <- i[1]; i}) %>%
	apply(1, . %>% paste(collapse = ""))
	
    hla_df <- hla_df %>%
	dplyr::filter(sub("^([^\\*]+).+$", "\\1", allele) == locus)

    if (is.null(exons) && !by_exon) {

	hla_df <- hla_df %>% 
	    dplyr::mutate(cds = gsub("\\|", "", cds))
    
    } else if (is.null(exons) && by_exon) {

	hla_df <- hla_df %>%
	    dplyr::mutate(cds = strsplit(cds, "\\|")) %>%
	    tidyr::unnest() %>%
	    dplyr::group_by(allele) %>%
	    dplyr::mutate(exon = seq_len(dplyr::n())) %>%
	    dplyr::select(allele, exon, cds) %>% 
	    dplyr::ungroup()
    
    } else if (!is.null(exons) && is.numeric(exons) && !by_exon) {
    
	hla_df <- hla_df %>%
	    dplyr::mutate(cds = strsplit(cds, "\\|")) %>%
	    tidyr::unnest() %>%
	    dplyr::group_by(allele) %>%
	    dplyr::mutate(exon = seq_len(dplyr::n())) %>%
	    dplyr::select(allele, exon, cds) %>% 
	    dplyr::filter(exon %in% exons) %>%
	    dplyr::summarise(cds = paste(cds, collapse = "")) %>%
	    dplyr::ungroup()
    
    } else if (!is.null(exons) && is.numeric(exons) && by_exon) {
	
	hla_df <- hla_df %>%
	    dplyr::mutate(cds = strsplit(cds, "\\|")) %>%
	    tidyr::unnest() %>%
	    dplyr::group_by(allele) %>%
	    dplyr::mutate(exon = seq_len(dplyr::n())) %>%
	    dplyr::ungroup() %>%
	    dplyr::select(allele, exon, cds) %>% 
	    dplyr::filter(exon %in% exons)
    }
   
    hla_df
}
