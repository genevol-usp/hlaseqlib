#' Process HLA alignment data in *.nuc files.
#'
#' Reads in nuc files with aligned sequences and processes them, returning a
#' matrix in which each row is the complete sequence for an allele.
#'
#' @param locus character string.
#' @param align_dir character string. Path to nuc file as a character string. 
#' Default is working directory.
#' @param omit_suffix character string. String of optional suffixes indicating
#' expression status. Possible values are NULL (default, keep all allele names
#' with suffixes) or a string of suffixes to omit (e.g., "NQL"). 
#' @param exons integer. Range of vectors to select (default = all exons).
#' @param undash logical. Should dashes be replaced with the base present in the
#' reference allele?
#' @param rm_dots_method character string. Method used to remove dots in the
#' sequence ("all" to remove only dots present in all alleles or "ref" to
#' remove those present in the first allele on the nuc file).
#' @return A matrix with the sequences from all alleles.
#' @export

hla_parsenuc <- function(locus, align_dir = getwd(), exons = NULL,
			 omit_suffix = NULL, undash = TRUE,
			 rm_dots_method = c("all", "ref")) {

  rm_dots_method <- match.arg(rm_dots_method)
  hla_pattern <- sprintf("^(%s\\d?\\*\\d{2,3}[:A-Z0-9]*[^ ])(.*)$", locus)
  
  seqs <-
    file.path(align_dir, stringr::str_c(locus, "_nuc.txt")) %>%
    readr::read_lines() %>%
    stringr::str_replace_all("[ ]{2,}", " ") %>% 
    stringr::str_trim() %>%
    stringr::str_subset(stringr::str_c("^", locus, "\\d?\\*"))

  first_allele <- stringr::str_replace(seqs[1], hla_pattern, "\\1")
 
  alleles <- stringr::str_replace_all(seqs, hla_pattern, "\\1")
  
  sequences <- 
    stringr::str_replace_all(seqs, hla_pattern, "\\2") %>%
    stringr::str_replace_all("[ ]", "")

  hla_df <-
    tibble::tibble(allele = alleles, sequence = sequences) %>%
    dplyr::group_by(allele) %>% 
    dplyr::summarize(sequence = stringr::str_c(sequence, collapse = ""))

  if (!is.null(omit_suffix)) {
    suffix <- sprintf("[^%s]$", omit_suffix)
    hla_df %<>% dplyr::filter(stringr::str_detect(allele, suffix))
  }
 
  max_len <- stringr::str_length(hla_df$sequence) %>% max()

  hla_df %<>% 
    dplyr::mutate(sequence = stringr::str_pad(sequence, max_len, "right", "."),
		  sequence = stringr::str_split(sequence, "\\|")) %>%
    tidyr::unnest()

  if (is.numeric(exons)) {  
      hla_df %<>%
	dplyr::group_by(allele) %>%
        dplyr::slice(exons)
  }
  
  hla_df %>%
    dplyr::group_by(allele) %>%
    dplyr::summarise(sequence = stringr::str_c(sequence, collapse = "")) %>%
    plyr::dlply(~allele, . %$% sequence %>% stringr::str_split("")) %>%
    purrr::flatten() %>%
    do.call(rbind, .) %>%
    hla_rmdels(method = rm_dots_method, ref_allele = first_allele) %>%
    {if (undash) hla_undash(., ref_allele = first_allele) else .} %>%
    `colnames<-`(seq_len(ncol(.)))
}
