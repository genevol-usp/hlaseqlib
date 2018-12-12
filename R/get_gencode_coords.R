#' Parse Gencode gtf files into a data.frame of ids and coordinates.
#'
#' @param gtf character string. Path to gtf file.
#' @param feature character string.
#'
#' @return data.frame
#' @export

get_gencode_coords <- function(gtf, feature = c("gene", "transcript", "exon")) {
  
  feature <- match.arg(feature)

  annot <-
    readr::read_tsv(gtf, comment = "##", col_names = FALSE, 
		    col_types = "c-cii-c-c", progress = FALSE) %>%
    dplyr::filter(X3 == feature) %>%
    dplyr::select(chr = X1, start = X4, end = X5, strand = X7, X9) %>%
    dplyr::mutate(i = seq_len(nrow(.)), X9 = stringr::str_split(X9, "; ")) %>%
    tidyr::unnest() %>%
    dplyr::filter(stringr::str_detect(X9, "^gene_name|^gene_id|^gene_type|^transcript_id")) %>%
    tidyr::separate(X9, c("tag", "id"), " ") %>%
    dplyr::mutate(chr = stringr::str_replace(chr, "^chr", ""),
		  id = stringr::str_replace_all(id, "\"", "")) %>%
    tidyr::spread(tag, id)

  if (feature == "gene") {
    annot %>% dplyr::select(chr, gene_id, gene_name, gene_type, start, end, strand)
  } else if (feature == "transcript") {
    annot %>% dplyr::select(chr, tx_id = transcript_id, gene_id, gene_name, 
			    start, end, strand)
  } else if (feature == "exon") {
    annot %>% dplyr::select(-i)
  } 
}
