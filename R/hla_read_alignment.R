#' Process HLA alignment data in *.nuc or *.gen files and return a data.frame.
#'
#' @param locus character string. The HLA locus to be processed (e.g., "A", "DRB1").
#' @param imgtdb character string. Path to the IMGTHLA directory.
#' @param imgtfile character string. Whether to process "nuc" or "gen" files. 
#' @param exons integer. Range of exons/introns (default is all exons/introns).
#' @param by_exon logical. Whether to return sequences separated for each exon/intron (default FALSE).
#' @param keep_sep logical. Whether to keep '|' as the separator between exons/introns.
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

hla_read_alignment <- function(locus, imgtdb, imgtfile = c("nuc", "gen"),
			       exons = NULL, by_exon = FALSE, keep_sep = FALSE) {
 
    imgtfile <- match.arg(imgtfile)

    if (imgtfile == "nuc") {

	locus_file <- ifelse(grepl("DRB\\d", locus), "DRB", locus)
	alig_file <- file.path(imgtdb, "alignments", paste0(locus_file, "_nuc.txt"))

    } else if (imgtfile == "gen") {

	locus_file <- locus
	alig_file <- file.path(imgtdb, "alignments", paste0(locus_file, "_gen.txt"))

    }

    alignments <- readLines(alig_file) %>%
	gsub("\\s{2,}", " ", .) %>%
	trimws %>%
	.[grepl(sprintf("^%s\\d?\\*\\d{2,3}[:A-Z0-9]*\\s", locus_file), .)]

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

    if (by_exon || (!is.null(exons) && is.numeric(exons))) {

	hla_df <- hla_df %>%
	    dplyr::mutate(cds = strsplit(cds, "\\|")) %>%
	    tidyr::unnest(cds) %>%
	    dplyr::group_by(allele) %>%
	    dplyr::mutate(exon = seq_len(dplyr::n())) %>%
	    dplyr::select(allele, exon, cds) %>% 
	    dplyr::ungroup()
    
	if (imgtfile == "gen") {
	    hla_df <- hla_df %>%
		dplyr::rename(idx = exon) %>%
		dplyr::mutate(feature = ifelse(idx %% 2 == 0, "exon", "intron"),
			      feature = dplyr::case_when(feature == "intron" & idx == 1 ~ "UTR",
							 feature == "intron" & idx == dplyr::last(idx) ~ "UTR",
							 TRUE ~ feature)) %>%
		dplyr::group_by(allele, feature) %>%
		dplyr::mutate(idx_grp = seq_len(dplyr::n())) %>% 
		dplyr::ungroup() %>%
		dplyr::select(allele, feature, idx, idx_grp, cds)
	} else if (imgtfile == "nuc") {
	    hla_df <- hla_df %>%
		dplyr::mutate(feature = "exon", idx = exon, idx_grp = exon) %>%
		dplyr::select(allele, feature, idx, idx_grp, cds)
	}

	if (!is.null(exons) && is.numeric(exons)) {

	    hla_df <- hla_df %>%
		dplyr::filter(feature == "exon", idx_grp %in% exons)
	}


    }
    
    if (!by_exon && !keep_sep) {

	hla_df <- hla_df %>% 
	    dplyr::mutate(cds = gsub("\\|", "", cds))
    }
   
    hla_df
}
