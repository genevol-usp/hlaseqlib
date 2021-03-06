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
    tidyr::unnest(cols = c(X9)) %>%
    dplyr::filter(grepl("^gene_name|^gene_id|^gene_type|^transcript_id|^transcript_type", X9)) %>%
    tidyr::separate(X9, c("tag", "id"), " ") %>%
    dplyr::mutate(chr = sub("^chr", "", chr),
		  id = gsub("\"", "", id)) %>%
    tidyr::spread(tag, id)

  if (feature == "gene") {
    annot %>% 
	dplyr::select(chr, gene_id, gene_name, gene_type, start, end, strand)
  } else if (feature == "transcript") {
    annot %>% 
	dplyr::select(chr, gene_id, gene_name, tx_id = transcript_id, 
		      tx_type = transcript_type, start, end, strand)
  } else if (feature == "exon") {
    annot %>% 
	dplyr::select(chr, gene_id, gene_name, tx_id = transcript_id, 
		      tx_type = transcript_type, start, end, strand)

  }
}
