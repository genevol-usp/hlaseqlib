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

  het_loci <-
    x %>%
      dplyr::group_by(locus) %>%
      dplyr::filter(dplyr::n_distinct(allele) == 2) %>%
      dplyr::pull(locus) %>%
      unique

  combinations <- 
    x %>%
      dplyr::group_by(locus, hap) %>%
      dplyr::mutate(i = seq_len(n())) %>%
      tidyr::unite(allele_diff, allele, diffs) %>%
      tidyr::spread(locus, allele_diff) %>%
      tidyr::fill(subject:DRB1) %>%
      tidyr::expand(hap, A, B, C, DQA1, DQB1, DRB1) %>%
      purrr::by_row(. %>% dplyr::select(A:DRB1) %>% sub("^.+_(\\d+)$", "\\1", .) %>% as.integer %>% sum,
		    .to = "S", .collate = "cols") %>%
      dplyr::mutate_at(dplyr::vars(A:DRB1), . %>% sub("_\\d+$", "", .))

  min_hap <- 
    combinations %>%
    dplyr::filter(S == min(S)) %>%
    dplyr::group_by(hap) %>%
    dplyr::mutate(n = n()) %>%
    dplyr::ungroup()

  if (any(min_hap$n == 1L)) min_hap <- dplyr::filter(min_hap, n == 1L)

  min_hap <- dplyr::select(min_hap, -n)

  if (all(c("a1", "a2") %in% min_hap$hap) && 
      all(x$allele %in% unlist(min_hap)) && 
      !any(duplicated(min_hap$hap))) {
     return(min_hap)
  }

  amb_loci <- 
    min_hap %>%
    tidyr::gather(locus, allele, -hap, -S) %>%
    dplyr::group_by(hap, locus) %>%
    dplyr::filter(dplyr::n_distinct(allele) > 1) %>%
    dplyr::pull(locus) %>%
    unique()

  min_alleles <- 
    min_hap[het_loci] %>%
    purrr::keep(~dplyr::n_distinct(.) == 1) %>% 
    unlist

  sec_hap <-
    combinations %>%
      dplyr::filter(hap != min_hap$hap,
		    ! A %in% min_alleles, ! B %in% min_alleles, 
		    ! C %in% min_alleles, ! DQA1 %in% min_alleles, 
		    ! DQB1 %in% min_alleles, ! DRB1 %in% min_alleles) %>%
      dplyr::filter(S == min(S))

  amb_in_hap2 <-
    sec_hap[amb_loci] %>%
    purrr::keep(~dplyr::n_distinct(.) == 1) %>%
    unlist

  min_hap <- min_hap %>%
    dplyr::filter(! A %in% amb_in_hap2, ! B %in% amb_in_hap2, 
		  ! C %in% amb_in_hap2, ! DQA1 %in% amb_in_hap2, 
		  ! DQB1 %in% amb_in_hap2, ! DRB1 %in% amb_in_hap2)

  dplyr::bind_rows(min_hap, sec_hap) %>%
    dplyr::group_by(hap) %>%
    dplyr::summarize_all(. %>% unique %>% paste(collapse = "/")) 
}
