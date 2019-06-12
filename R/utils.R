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

hla_format_sequence <- function(cds) {

    cds %>%
	gsub("\\.|\\|", "", .) %>%
	gsub("\\*", "N", .)
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

