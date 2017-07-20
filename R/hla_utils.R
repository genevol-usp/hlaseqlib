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
    dplyr::mutate(allele = gsub("^IMGT_|_s\\d+$", "", allele) %>%
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

hla_df_to_matrix <- function(hla_df) {

  hla_mat <- 
    hla_df$sequence %>%
    stringr::str_split("", simplify = TRUE) %>%
    `rownames<-`(hla_df$allele)

 hla_mat[hla_mat == ""] <- "."  

 hla_mat
}

hla_format_sequence <- function(cds) {

  cds %>%
    stringr::str_replace_all("\\.|\\|", "") %>%
    stringr::str_replace_all("\\*", "N")
}

