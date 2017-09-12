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
  
  hla_quants[, r := data.table::frank(-tpm, ties.method = "dense"), by = .(subject, locus)]

  top2 <- 
    hla_quants[r == 1 | r == 2
	     ][, maxtpm := max(tpm), by = .(subject, locus)
	     ][maxtpm == 0 | (tpm/maxtpm > th)]
  
  out <- 
    top2[, .(allele = paste(allele, collapse = "="), est_counts = mean(est_counts)), 
         by = .(subject, gene_id, locus, tpm)
       ][, n := .N, by = .(subject, gene_id, locus)
       ][n == 1, `:=`(est_counts = est_counts/2L, tpm = tpm/2L)
       ][, i := ifelse(n == 1, "1_1", "1")
       ][, .(i = unlist(strsplit(i, "_"))), 
         by = .(subject, gene_id, locus, allele, tpm, est_counts)
       ][, .(subject, gene_id, locus, allele, est_counts, tpm)
       ][order(subject, locus, allele)]

  out
}
