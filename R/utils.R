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

alleles_to_groups <- function(df, input_fields = 4L) {

  hla_groups <-
    dplyr::filter(hla_groups, locus %in% df$locus) %>%
    dplyr::mutate(allele = hla_trimnames(allele, input_fields)) %>%
    dplyr::distinct()

  df %>%  
    dplyr::group_by(subject, locus) %>%
    dplyr::mutate(e = 1:n()) %>%
    tidyr::separate_rows(allele, sep = "=") %>%
    dplyr::mutate(h = 1:n()) %>%
    tidyr::separate_rows(allele, sep = "/") %>%
    dplyr::mutate(allele = gsub("^IMGT_", "", allele) %>%
		  hla_trimnames(fields = input_fields)) %>%
    dplyr::left_join(hla_groups, by = c("locus", "allele")) %>%
    dplyr::mutate(group = ifelse(is.na(group), allele, group)) %>%
    dplyr::group_by(subject, locus, h) %>%
    dplyr::summarize(e = unique(e),
		     allele = paste(unique(group), collapse = "/")) %>%
    dplyr::group_by(subject, locus, e) %>%
    dplyr::summarize(allele = paste(unique(allele), collapse = "=")) %>%
    dplyr::ungroup() %>% 
    dplyr::select(-e) %>%
    dplyr::arrange(subject, locus, allele)
}

hla_attribute_seq <- function(incomplete_seq, complete_seq) {

  incomplete_seq %<>% strsplit("") %>% unlist
  complete_seq %<>% strsplit("") %>% unlist

  unseq <- incomplete_seq == "*"

  incomplete_seq[unseq] <- complete_seq[unseq]

  paste(incomplete_seq, collapse = "")
}

find_closest_allele <- function(incomplete, complete_df) {

  incomp <- strsplit(incomplete, "")
  
  comps <- 
    complete_df %>%
    split(.$allele) %>%
    purrr::map_chr("cds") %>%
    strsplit("")

  seqbases <- which(unlist(incomp) != "*")

  incomp_seq <- purrr::map(incomp, ~.[seqbases])
  comps_seqs <- purrr::map(comps, ~.[seqbases])
  
  purrr::map2_int(comps_seqs, incomp_seq, ~sum(.y != .x)) %>%
  which.min() %>%
  names()
}

hla_make_sequences <- function(l, infer_incomplete, n_cores) {

  paste0("/home/vitor/IMGTHLA/alignments/", l, "_nuc.txt") %>%
  hla_read_alignment(omit = "N", rm_incomplete = !infer_incomplete) %>%
  dplyr::mutate(locus = sub("^([^*]+).+$", "\\1", allele)) %>%
  dplyr::select(locus, allele, cds) %>%
  dplyr::group_by(locus) %>%
  dplyr::filter(!all(grepl("\\*", cds))) %>%
  dplyr::ungroup() %>%
  split(.$locus) %>%
  purrr::map_df(~hla_infer(., cores = n_cores)) %>%
  dplyr::select(-locus)
}

hla_format_sequence <- function(cds) {

  cds %>%
    gsub("\\.|\\|", "", .) %>%
    gsub("\\*", "N", .)
}

read_star_imgt_quants <- function(f) {

    read_tsv(f, col_types = "c--dd", progress = FALSE) %>%
        filter(grepl("^IMGT_", Name)) %>%
        mutate(locus = imgt_to_gname(Name), gene_id = gname_to_gid(locus)) %>%
        select(locus, gene_id, allele = Name, est_counts = NumReads, tpm = TPM)
}

read_star_pri_quants <- function(f) {

    read_tsv(f, col_types = "c--dd", progress = FALSE) %>%
        left_join(gencode_pri_tx, by = c("Name" = "tx_id")) %>%
        select(tx_id = Name, locus = gene_name,
               est_counts = NumReads, tpm = TPM) %>%
        filter(locus %in% hla_genes) %>%
        group_by(locus) %>%
        summarize_at(vars(tpm, est_counts), sum) %>%
        ungroup()
}

read_kallisto_imgt_quants <- function(f) {

    read_tsv(f, col_types = "c--dd", progress = FALSE) %>%
        filter(grepl("^IMGT_", target_id)) %>%
        mutate(locus = imgt_to_gname(target_id), 
	       gene_id = gname_to_gid(locus)) %>%
        select(locus, gene_id, allele = target_id, est_counts, tpm)
}

read_kallisto_pri_quants <- function(f) {

    read_tsv(f, col_types = "c--dd", progress = FALSE) %>%
        left_join(gencode_pri_tx, by = c("target_id" = "tx_id")) %>%
        select(tx_id = target_id, locus = gene_name, est_counts, tpm) %>%
        filter(locus %in% hla_genes) %>%
        group_by(locus) %>%
        summarize_at(vars(tpm, est_counts), sum) %>%
        ungroup()
}
