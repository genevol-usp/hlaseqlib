#' Read hlaTX output file
#'
#' @param x Character string. Path to hlatx output file.
#' @return A data.frame.

hla_readtx <- function(x) {

  tx_names <- 
    c("A", "B", "C", "DQA1", "DQB1", "DRB1") %>%
    rep(each = 2) %>%
    paste(1:2, sep = "_") %>%
    c("subject", ., paste0(., "_est_counts"))

  hla <-
    readr::read_csv(x, col_names = FALSE) %>%
    dplyr::select(-(X14:X25)) %>%
    setNames(tx_names) %>%
    dplyr::mutate(subject = sub("\\.0$", "", subject)) %>%
    dplyr::distinct() %>%
    dplyr::arrange(subject)

    hla_tidy <- tibble::tibble(subject = rep(hla$subject, each = 12), 
			       allele = NA, est_counts = NA)

    counter <- 1
    for (i in seq_len(nrow(hla))) {
      hla_tidy$allele[counter:(counter+11)] <- unlist(hla[i, 2:13])
      hla_tidy$est_counts[counter:(counter+11)] <- unlist(hla[i, 14:25])
      counter <- counter + 12
   }

   hla_tidy
}

