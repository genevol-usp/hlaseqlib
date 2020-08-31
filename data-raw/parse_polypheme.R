devtools::load_all(".")
library(tidyverse)
library(readxl)

dat <- read_excel("./polypheme_calls.xlsx") %>%
    select(subject = 3, A.1 = 4, A.2 = 5, B.1 = 6, B.2 = 7, C.1 = 8, C.2 = 9,
	   DQB1.1 = 10, DQB1.2 = 11, DRB1.1 = 12, DRB1.2 = 13) %>%
    pivot_longer(-subject, names_to = "locus", values_to = "allele") %>%
    separate(locus, c("locus", "hap"), sep = "\\.") %>%
    mutate(allele = na_if(allele, "None"),
	   allele = ifelse(grepl("\\.", allele), NA, allele)) %>%
    group_by(subject, locus) %>%
    filter(!any(is.na(allele))) %>%
    ungroup() %>%
    separate_rows(allele, sep = " ") %>%
    separate(allele, c("lineage", "allele"), sep = ":") %>%
    mutate(allele = sub("\\*$", "", allele)) %>%  
    separate_rows(allele, sep = "/") %>% 
    mutate(allele = paste0(locus, "*", lineage, ":", allele)) %>%
    group_by(subject, locus, hap) %>%
    summarise(allele = paste(allele, collapse = "/")) %>%
    ungroup() %>%
    arrange(subject, locus, hap) %>% 
    select(-hap) 

complete_genos <- anti_join(pag, dat, by = c("subject", "locus")) %>%
    bind_rows(dat) %>%
    arrange(subject, locus, allele)

write_tsv(complete_genos, "polypheme_plus_pag_calls.tsv")
