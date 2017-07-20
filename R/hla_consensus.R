#' Create consensus sequences for a list of HLA alleles.
#'
#' \code{hla_consensus()} takes a matrix of HLA sequences and returns a matrix 
#' of consensus sequences at the resolution specified by the argument
#' \code{field}.
#'
#' @param hla_mat matrix.
#' @param fields integer. Resolution of alleles in fields.
#' @param cores integer. Number of cores if running in parallel.
#'
#' @return A matrix of consensus sequences. Ambiguous positions contain the
#' occurring bases separated by "/".
#'
#' @export

hla_consensus <- function(hla_mat, fields, cores = 1) {

  hla_mat %>%
   `rownames<-`(hla_trimnames(rownames(.), fields)) %>%
    split.data.frame(factor(rownames(.))) %>% 
   parallel::mclapply(. %>% apply(2, . %>% unique() %>% paste0(collapse = "/")), 
		      mc.cores = cores) %>%
   do.call(rbind, .) %>%
   `colnames<-`(seq_len(ncol(.)))
}
