#' Calculate edit distances to the reference HLA allele.
#' 
#' @param hla_mat matrix. A matrix returned by \code{hla_parsenuc}.
#' @param ref_allele character string. The name of the reference allele in the 
#' same format of rownames(hla_mat).
#' @param stat character string. Return absulute number of differences ("count")
#' or proportion ("proportion") 
#' 
#' @return A vector of edit distances of each allele in respect to the reference 
#' allele.
#' 
#' @export

hla_refdist <- function(hla_mat, ref_allele, stat = c("count", "proportion")) {

  stat <- match.arg(stat)
  
  hla_mat %<>% .[apply(., 1, function(x) !any(stringr::str_detect(x, "\\*"))), ]
 
  if (stat == "count") {
    apply(hla_mat, 1, function(i) 
	      adist(hla_getseq(hla_mat[ref_allele, ]), hla_getseq(i)))
  
  } else if (stat == "proportion") {
  apply(hla_mat, 1, 
	function(i) {
	  seq_ref <- hla_getseq(hla_mat[ref_allele, ]) 
	  seq_allele <- hla_getseq(i) 
	  adist(seq_ref, seq_allele)/nchar(seq_allele)
        })
  } else {
    stop("argument stat must be either 'count' or 'proportion'")
  }
}

