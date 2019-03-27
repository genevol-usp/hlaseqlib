devtools::load_all("/home/vitor/Libraries/hlaseqlib")
library(tidyverse)
library(readxl)

samples_phase3 <- 
    "ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/integrated_call_samples_v3.20130502.ALL.panel" %>%
    read_tsv(skip = 1, col_names = FALSE) %>%
    select(subject = X1, pop = X2)

dat <- read_excel("./polypheme_calls.xlsx") %>%
    select(subject = 3, A.1 = 4, A.2 = 5, B.1 = 6, B.2 = 7, C.1 = 8, C.2 = 9,
	   DQB1.1 = 10, DQB1.2 = 11, DRB1.1 = 12, DRB1.2 = 13) %>%
    gather(locus, allele, -1) %>%
    separate(locus, c("locus", "hap"), sep = "\\.") %>%
    inner_join(samples_phase3, by = "subject") %>%
    select(subject, pop, locus, hap, allele) %>%
    filter(pop %in% c("CEU", "FIN", "GBR", "TSI", "YRI"), !is.na(allele)) %>%
    separate(allele, c("lineage", "allele"), sep = ":") %>%
    mutate(allele = sub("\\*$", "", allele)) %>%  
    separate_rows(allele, sep = "/") %>% 
    mutate(allele = paste0(locus, "*", lineage, ":", allele)) %>%
    group_by(subject, locus, hap) %>%
    summarise(allele = paste(allele, collapse = "/")) %>%
    ungroup() %>%
    arrange(subject, locus, hap) %>% 
    complete(subject, locus, hap) %>%
    select(-hap) 

dat_na <- filter(dat, is.na(allele))

genotypes_to_complete <- 
    inner_join(pag, dat_na, by = c("subject", "locus")) %>%
    distinct() %>%
    select(subject, locus, allele = allele.x) %>%
    group_by(subject, locus) %>%
    mutate(i = ifelse(n_distinct(allele) == 1, "1_1", "1")) %>%
    ungroup() %>%
    separate_rows(i, sep = "_") %>%
    select(-i)

dat_na %>% filter(subject == "HG00108")

complete_genos <- dat %>%
    filter(!is.na(allele)) %>%
    bind_rows(genotypes_to_complete) %>%
    arrange(subject, locus, allele)

write_tsv(complete_genos, "polypheme_plus_pag_calls.tsv")

