#' Read hlaTX output file
#'
#' @param x Character string. Path to hlatx output file.
#' @return A data.frame.

hla_readtx <- function(x) {

  tx_names <- 
    c("A", "B", "C", "DQA1", "DQB1", "DRB1") %>%
    rep(each = 2) %>%
    stringr::str_c(1:2, sep = "_") %>%
    c("sample", ., stringr::str_c(., "_est_counts"))

  hla <-
    readr::read_csv(x, col_names = FALSE) %>%
    dplyr::select(-(X14:X25)) %>%
    purrr::set_names(tx_names) %>%
    dplyr::mutate(sample = stringr::str_replace(sample, "\\.0$", "")) %>%
    dplyr::distinct() %>%
    dplyr::arrange(sample)

    hla_tidy <- 
      tibble::tibble(subject = rep(hla$sample, each = 12), 
		     allele = NA, est_counts = NA)

    counter <- 1
    for (i in seq_len(nrow(hla))) {
      hla_tidy$allele[counter:(counter+11)] <- unlist(hla[i, 2:13])
      hla_tidy$est_counts[counter:(counter+11)] <- unlist(hla[i, 14:25])
      counter %<>% `+`(12)
   }

   hla_tidy
}

