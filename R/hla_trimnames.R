#' Trim HLA allele names to the specified number of fields.
#'
#' Given a vector of HLA allele names with the format "A*01:01:01:01", returns
#' all or unique allele names at the digit resolution specified by 
#' \code{fields}.
#'
#' @param alleles_vec character string. Vector of allele names.
#' @param fields integer. Number of fields to trim names to.
#'
#' @return All or unique HLA alleles trimmed names.
#' @examples
#' hla_trimnames(c("A*01:01:01:01", "A*02:01:01:01"), fields = 2)
#' @export
hla_trimnames <- function(alleles_vec, fields = 3) {

  regex <- sprintf("^([A-Z]{1,3}\\d?\\*)((:?\\d{2,3}[NQLS]?){%d}).*$", fields)

  strsplit(alleles_vec, "/") %>%
  purrr::map(~sub(regex, "\\1\\2", .)) %>%
  purrr::map(unique) %>%
  purrr::map_chr(~paste(., collapse = "/"))
}
