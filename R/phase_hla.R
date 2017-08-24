#' Phase HLA alleles given distances to 1000G haplotypes
#'
#' @param x data.frame. A data.frame with distances of alleles to each 1000G
#' haplotype.
#'
#' @return A data.frame with 2 rows with inferred haplotypes for each
#' individual.
#'
#' @export

phase_hla <- function(x) {

  x <- dplyr::mutate(x, locus = sub("^HLA-", "", locus))

  combinations <- 
    x %>%
      dplyr::group_by(locus, hap) %>%
      dplyr::mutate(i = seq_len(n())) %>%
      tidyr::unite(allele_diff, allele, diffs) %>%
      tidyr::spread(locus, allele_diff) %>%
      tidyr::expand(hap, A, B, C, DQA1, DQB1, DRB1) %>%
      purrr::by_row(. %>% dplyr::select(A:DRB1) %>% sub("^.+_(\\d+)$", "\\1", .) %>% as.integer %>% sum,
		    .to = "S", .collate = "cols") %>%
      dplyr::mutate_at(dplyr::vars(A:DRB1), . %>% sub("_\\d+$", "", .))

  h1s <- dplyr::filter(combinations, hap == "a1") %>% dplyr::mutate(i = "A")
  h2s <- dplyr::filter(combinations, hap == "a2") %>% dplyr::mutate(i = "A")

 test.x <- dplyr::left_join(h1s, h2s, by = "i") %>% 
   dplyr::select(-i) %>%
   dplyr::mutate(r = 1:n())
  
 test.y <- test.x %>%
   dplyr::select(-S.x, -S.y) %>%
   tidyr::gather(locus, allele, -hap.x, -hap.y, -r) %>%
   dplyr::arrange(r, hap.x, hap.y, locus) %>%
   dplyr::group_by(r) %>%
   dplyr::filter(all(x$allele %in% allele)) %>%
   dplyr::ungroup()

 test.z <- dplyr::filter(test.x, r %in% test.y$r) %>%
   dplyr::select(-r) %>%
   dplyr::mutate(S = S.x + S.y) %>%
   dplyr::filter(S == min(S))

 out.x <- dplyr::select(test.z, dplyr::ends_with(".x")) %>%
   setNames(sub("\\.x", "", names(.)))
 out.y <- dplyr::select(test.z, dplyr::ends_with(".y")) %>%
   setNames(sub("\\.y", "", names(.)))

 out <- dplyr::bind_rows(out.x, out.y) %>%
   dplyr::arrange(hap)

 out %>%
    dplyr::group_by(hap) %>%
    dplyr::summarize_all(. %>% unique %>% paste(collapse = "/")) 
}
