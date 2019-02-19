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

    readr::read_tsv(f, col_types = "cdd", progress = FALSE) %>%
        dplyr::filter(grepl("^IMGT_", Name)) %>%
        dplyr::mutate(locus = imgt_to_gname(Name), 
		      gene_id = gname_to_gid(locus)) %>%
        dplyr::select(locus, gene_id, allele = Name, 
		      est_counts = NumReads, tpm = TPM)
}

read_star_pri_quants <- function(f) {

    readr::read_tsv(f, col_types = "cdd", progress = FALSE) %>%
        dplyr::left_join(gencode_pri_tx, by = c("Name" = "tx_id")) %>%
        dplyr::select(tx_id = Name, locus = gene_name,
		      est_counts = NumReads, tpm = TPM) %>%
        dplyr::filter(locus %in% hla_genes) %>%
        dplyr::group_by(locus) %>%
        dplyr::summarize_at(vars(tpm, est_counts), sum) %>%
        dplyr::ungroup()
}

read_kallisto_imgt_quants <- function(f) {

    readr::read_tsv(f, col_types = "cdd", progress = FALSE) %>%
        dplyr::filter(grepl("^IMGT_", target_id)) %>%
        dplyr::mutate(locus = imgt_to_gname(target_id), 
		      gene_id = gname_to_gid(locus)) %>%
        dplyr::select(locus, gene_id, allele = target_id, est_counts, tpm)
}

read_kallisto_pri_quants <- function(f) {

    readr::read_tsv(f, col_types = "cdd", progress = FALSE) %>%
        dplyr::left_join(gencode_pri_tx, by = c("target_id" = "tx_id")) %>%
        dplyr::select(tx_id = target_id, locus = gene_name, est_counts, tpm) %>%
        dplyr::filter(locus %in% hla_genes) %>%
        dplyr::group_by(locus) %>%
        dplyr::summarize_at(vars(tpm, est_counts), sum) %>%
        dplyr::ungroup()
}

make_dist_matrix <- function(hla_df) {

    complete_alleles <- hla_df %>% 
	dplyr::filter(!grepl("\\*", cds)) %>%
	dplyr::pull(allele)

    incomplete_alleles <- hla_df %>% 
	dplyr::filter(grepl("\\*", cds)) %>%
	dplyr::pull(allele)

    cds_sequenced <- stringr::str_split(hla_df$cds, "", simplify = TRUE) %>%
	apply(2, function(x) !any(x == "*"))

    run <- rle(cds_sequenced)

    ends <- cumsum(run$lengths)
    starts <- ends - run$lengths + 1L

    run_df <- 
	tibble::tibble(value = run$values, 
		       start = starts, 
		       end = ends) %>%
	dplyr::filter(value == TRUE) %>%
	dplyr::slice(which.max(end - start + 1))

    hla_df_cds_common <- hla_df %>%
	dplyr::mutate(cds = substring(cds, run_df$start, run_df$end))

    hla_df_cds_common_complete <- hla_df_cds_common %>%
	dplyr::filter(allele %in% complete_alleles)

    hla_df_cds_common_incomplete <- hla_df_cds_common %>%
	dplyr::filter(allele %in% incomplete_alleles)

    cds_common_complete <- hla_df_cds_common_complete$cds
    names(cds_common_complete) <- hla_df_cds_common_complete$allele
    
    cds_common_incomplete <- hla_df_cds_common_incomplete$cds
    names(cds_common_incomplete) <- hla_df_cds_common_incomplete$allele

    stringdist::stringdistmatrix(cds_common_incomplete, cds_common_complete,
				 method = "hamming", useNames = "names", 
				 nthread = 1)
}
    
make_closest_allele_df <- function(distmatrix) {

    distmatrix %>%
	split(seq_len(nrow(distmatrix))) %>%
	setNames(rownames(distmatrix)) %>%
	lapply(function(x) which(x == min(x, na.rm = TRUE))) %>% 
	purrr::map(~tibble::tibble(closest = colnames(distmatrix)[.])) %>%
	dplyr::bind_rows(.id = "inc_allele")
}

find_closest_within_type <- function(closest_df) {

    closest_df %>% 
    dplyr::mutate(`1` = hla_trimnames(inc_allele, 1) == hla_trimnames(closest, 1),
		  `2` = hla_trimnames(inc_allele, 2) == hla_trimnames(closest, 2),
		  `3` = hla_trimnames(inc_allele, 3) == hla_trimnames(closest, 3),
		  `4` = hla_trimnames(inc_allele, 4) == hla_trimnames(closest, 4)) %>%
    tidyr::gather(field, value, `1`:`4`) %>%
    dplyr::group_by(inc_allele) %>%
    dplyr::filter(all(value == FALSE) | value == TRUE) %>%
    dplyr::slice(which.max(field)) %>%
    dplyr::ungroup() %>%
    dplyr::arrange(inc_allele) %>%
    dplyr::select(inc_allele, closest)
}
