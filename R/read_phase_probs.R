#' Read in PHASE output file
#'
#' Reads PHASE output file and returns a data.frame with the phasing
#' probabilities.
#'
#' @param phase_out character string. Path to PHASE output file.
#'
#' @return A data.frame.
#'
#' @export

read_phase_probs <- function(phase_out) {

  x <- readr::read_lines(phase_out)

  x_inds <-
    stringr::str_detect(x, "(BEGIN|END) BESTPAIRS_SUMMARY") %>%
    which() %>%
    {(.[[1]] + 1):(.[[2]] - 1)} %>%
    x[.] %>%
    stringr::str_split(":") %>%
    purrr::map_chr(1) %>%
    stringr::str_replace_all("#", "")

  x_probs <- 
    stringr::str_detect(x, "(BEGIN|END) PHASEPROBS") %>% 
    which() %>%
    {(.[[1]] + 1):(.[[2]] - 1)} %>%
    x[.] %>%
    stringr::str_trim()

  x_df <- 
    tibble::tibble(subject = x_inds,
		   locus = stringr::str_c(gencode25_hla_coords$gene_name, 
					  collapse = " "),
		   probs = x_probs) %>%
    tidyr::separate_rows(locus, probs, sep = " ") %>%
    dplyr::arrange(subject, locus)

  x_df
}
