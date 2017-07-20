#' Convert a matrix of sequences in a list
#'
#' Convert a matrix of sequences, such that returned by hla_parsenuc, into a
#' list in which each element is the collapse sequence of each allele.
#'
#' @param hla_mat matrix. A matrix such the one returned by \code{hla_parsenuc}.
#'
#' @return a list.
#'
#' @export

hla_matrix_to_list <- function(hla_mat) {

  hla_mat %>%
  split(seq_len(nrow(.))) %>%
  `names<-`(rownames(hla_mat)) %>%
  purrr::map(~stringr::str_c(., collapse = ""))
}
