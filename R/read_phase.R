#' Read in PHASE output file
#'
#' Reads PHASE output file and returns a data.frame with best haplotypes per
#' individual
#'
#' @param phase_out character string. Path to PHASE output file.
#' @param loci character string. Loci names in genomic order.
#'
#' @return A data.frame.
#'
#' @export

read_phase <- function(phase_out, loci) {

  x <- readLines(phase_out)

  x_best <- 
    grep("(BEGIN|END) BESTPAIRS1", x) %>% 
    {(.[[1]] + 1):(.[[2]] - 1)} %>%
    x[.] 
    
  x_list <- split(x_best, grepl("#", x_best) %>% cumsum)

  tibble::tibble(subject = purrr::map_chr(x_list, 1),
		 locus = paste(loci, collapse = " "),
		 h1 = trimws(purrr::map_chr(x_list, 2)),
		 h2 = trimws(purrr::map_chr(x_list, 3))) %>%
  dplyr::mutate(subject = gsub("0 #", "", subject)) %>%
  tidyr::separate_rows(locus, h1, h2, sep = " ") %>%
  dplyr::mutate(uncertain = as.integer(grepl("\\(", h1))) %>%
  dplyr::mutate_at(dplyr::vars(h1, h2), . %>% gsub("\\(|\\)", "", .)) %>%
  tidyr::gather(hap, allele, h1:h2) %>%
  dplyr::mutate(hap = gsub("h", "", hap),
		allele = as.integer(allele)) %>%
  dplyr::select(subject, locus, hap, allele, uncertain) %>%
  dplyr::arrange(subject, locus, hap)
}
