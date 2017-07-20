tr_to_gid <- function(tr) 
  gencode_all_tx$gene_id[match(tr, gencode_all_tx$tx_id)]

tr_to_gname <- function(tr)
  gencode_all_tx$gene_name[match(tr, gencode_all_tx$tx_id)]

gname_to_gid <- function(gname)
  gencode_all_tx$gene_id[match(gname, gencode_all_tx$gene_name)]

gid_to_gname <- function(gid)
  gencode_all_tx$gene_name[match(gid, gencode_all_tx$gene_id)]

convert_ena_ids <- function(ena_id)
  geuvadis_info$name[match(ena_id, geuvadis_info$ena_id)]

convert_kgp_to_ena <- function(kgp)
  geuvadis_info$ena_id[match(kgp, geuvadis_info$name)]

imgt_to_gname <- function(ids) {
    
  loci <- sub("^IMGT_([A-Z1-9]+)\\*.+$", "\\1", ids)
    
  ifelse(grepl("MIC|TAP|HFE", loci), loci, paste0("HLA-", loci))
}

imgt_to_gid <- function(ids) {
    
  ids %>%
    imgt_to_gname %>%
    gname_to_gid
}

targetid_to_geneid <- function(ids) {
 
  ifelse(grepl("^IMGT", ids), imgt_to_gid(ids), tr_to_gid(ids))
}

code_alleles <- function(genotypes) {

  allele_codes <-
    genotypes %>%
    dplyr::select(locus, allele) %>%
    dplyr::arrange(locus, allele) %>%
    dplyr::distinct() %>%
    dplyr::group_by(locus) %>%
    dplyr::mutate(code = seq_len(n())) 

  genotypes %>%
    dplyr::left_join(allele_codes, by = c("locus", "allele"))
}
