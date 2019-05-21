#' Infer HLA genotypes from RNA-seq quantifications.
#'
#' @param hla_quants data.frame.
#' @param th numeric. Threshold for expression 2nd/1st allele.
#'
#' @return A data.frame
#'
#' @export

hla_genotype_dt <- function(hla_quants, th = 0.05) {

  data.table::setDT(hla_quants)
  
  hla_quants[, r := data.table::frank(-est_counts, ties.method = "dense"), by = .(subject, locus)]

  top2 <- 
    hla_quants[r == 1 | r == 2
	     ][, maxcounts := max(est_counts), by = .(subject, locus)
	     ][maxcounts == 0 | (est_counts/maxcounts > th)]
  
  out <- 
    top2[, .(allele = paste(allele, collapse = "="), tpm = mean(tpm)), 
         by = .(subject, locus, est_counts)
       ][, n := .N, by = .(subject, locus)
       ][n == 1, `:=`(est_counts = est_counts/2L, tpm = tpm/2L)
       ][, i := ifelse(n == 1, "1_1", "1")
       ][, .(i = unlist(strsplit(i, "_"))), 
         by = .(subject, locus, allele, tpm, est_counts)
       ][, .(subject, locus, allele, est_counts, tpm)
       ][order(subject, locus, allele)]

  out
}
