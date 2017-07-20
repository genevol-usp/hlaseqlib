#' Write input for PHASE
#'
#' Takes a data.frame of HLA genotypes and writes a input file for PHASE
#'
#' @param genotypes data.frame. A data.frame with cols subject, locus, allele.
#' @param outfile character string. Path to output file.
#'
#' @return No return.
#'
#' @export

write_phase_input <- function(genotypes, outfile) {

  file.create(outfile)

  n_ind <- dplyr::n_distinct(genotypes$subject)
  n_loci <- dplyr::n_distinct(genotypes$locus)

  hla_coords <-
    gencode_chr_gene %>%
    dplyr::filter(gene_name %in% paste0("HLA-", unique(genotypes$locus))) %>%
    dplyr::mutate(gene_name = sub("HLA-", "", gene_name)) %>%
    dplyr::arrange(start)
    
  hla_pos <- c("P", hla_coords$start) %>% paste(collapse = " ")

  m_types <- rep("M", n_loci) %>% paste(collapse = "")

  genotypes_recoded <- code_alleles(genotypes)

  genotypes_recoded %>% 
    dplyr::select(locus, allele, code) %>% 
    dplyr::arrange(locus, code) %>%
    dplyr::distinct() %>%
    readr::write_tsv(stringr::str_c("codes-", outfile))

  haps_recoded <-
    genotypes_recoded %>%
    dplyr::group_by(subject, locus) %>%
    dplyr::mutate(h = seq_len(n())) %>%
    dplyr::select(subject, locus, h, code) %>%
    tidyr::spread(locus, code) %>%
    dplyr::select(subject, hla_coords$gene_name) %>%
    dplyr::ungroup() 

  write(n_ind, outfile, append = TRUE)
  write(n_loci, outfile, append = TRUE)
  write(hla_pos, outfile, append = TRUE)
  write(m_types, outfile, append = TRUE)

  haps_recoded %>%
    plyr::d_ply(~subject, 
		function(x) {
		  write(paste0("#", unique(x$subject)), outfile, append = TRUE)
		  write(unlist(x[1, -1]), outfile, ncolumns = n_loci, append = TRUE)
		  write(unlist(x[2, -1]), outfile, ncolumns = n_loci, append = TRUE)
		})
}
