devtools::load_all(".")
library(tidyverse)

translate_allele <- function(allele) {

    regex <- paste0("^", sub("\\*", "\\\\*", allele), "([NQLS]|:|$)")

    filter(allele_hist, grepl(regex, `3010`)) %>%
	drop_na(2) %>%
	pull(2) %>%
	unique() %>% 
	hla_minres() %>%
	paste(collapse = "/")
}

hla_readmhc <- function(file) {
  
    pag <-
	read_tsv(file, col_types = "-c---cccccccccc") %>%
	setNames(gsub("-|[ ]", "_", names(.))) %>%
	rename(subject = Subject) %>%
	gather(locus, allele, HLA_A_1:HLA_DQB1_2) %>% 
	extract(locus, c("locus", "h"), "HLA_([^_]+)+_(\\d)") %>% 
	separate_rows(allele, sep = "/") %>% 
	mutate(subject = sub("_\\d$", "", subject),
	allele = sub("XX$", "", paste(locus, allele, sep = "*")),
	allele = recode(allele, "C*0140" = "C*01:40")) %>%
	arrange(subject, locus, allele)

    pag_alleles <- tibble(allele = unique(pag$allele)) %>%
	mutate(new = map_chr(allele, translate_allele)) %>%
	select(allele, new)

    pag %>%
	left_join(pag_alleles, by = "allele") %>%
	select(subject, locus, h, allele = new) %>%
	group_by(subject, locus, h) %>%
	summarize(allele = paste(allele, collapse = "/")) %>%
	select(-h) %>%
	ungroup()
}

samples_phase3 <- 
    "ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/integrated_call_samples_v3.20130502.ALL.panel" %>%
    read_tsv(skip = 1, col_names = FALSE) %>%
    select(subject = X1, pop = X2, continent = X3, sex = X4)

allele_hist <- read_tsv("~/IMGTHLA/Allelelist_history.txt")

geuvadis_info <- 
    "http://www.ebi.ac.uk/arrayexpress/files/E-GEUV-1/E-GEUV-1.sdrf.txt" %>%
    read_tsv() %>%
    select(name = `Source Name`, ena_id = `Comment[ENA_RUN]`, 
	   assay_name = `Assay Name`, lab = Performer, 
	   lab_code = `Factor Value[laboratory]`, 
	   pop = `Characteristics[population]`, 
	   kgp_phase1 = `Comment[1000g Phase1 Genotypes]`) %>%
    distinct() %>%
    mutate(kgp_phase3 = as.integer(name %in% samples_phase3$subject)) 


gencode_chr_gene <-
    "~/gencode_data/v32/gencode.v32.annotation.gtf" %>%
    get_gencode_coords(feature = "gene")

gencode_pri_gene <-
    "~/gencode_data/v32/gencode.v32.primary_assembly.annotation.gtf.gz" %>%
    get_gencode_coords(feature = "gene")

gencode_all_gene <-
    "~/gencode_data/v32/gencode.v32.chr_patch_hapl_scaff.annotation.gtf.gz" %>%
    get_gencode_coords(feature = "gene")

gencode_chr_tx <-
    "~/gencode_data/v32/gencode.v32.annotation.gtf" %>%
    get_gencode_coords(feature = "transcript")

gencode_pri_tx <-
    "~/gencode_data/v32/gencode.v32.primary_assembly.annotation.gtf.gz" %>%
    get_gencode_coords(feature = "transcript")

gencode_all_tx <-
    "~/gencode_data/v32/gencode.v32.chr_patch_hapl_scaff.annotation.gtf.gz" %>%
    get_gencode_coords(feature = "transcript")

hla_genes <- 
    paste0("HLA-", c("A", "B", "C", "DPA1", "DPB1", "DQA1", "DQB1", "DRA", "DRB1"))

gencode_hla <- filter(gencode_chr_gene, gene_name %in% hla_genes) %>%
    select(-gene_type)

hla_groups <-
    "~/IMGTHLA/wmda/hla_nom_g.txt" %>%
    read_delim(delim = ";", col_names = FALSE, comment = "#") %>%
    separate_rows(X2, sep = "/") %>%
    transmute(allele = paste0(X1, X2), 
	      group = paste0(X1, X3),
	      locus = sub("\\*$", "", X1)) %>%
    select(locus, group, allele) %>%
    arrange(locus, group, allele)
  
pag <- hla_readmhc("./mhc.tab") 

pag_groups <- alleles_to_groups(pag) 

polypheme_pag <- read_tsv("./polypheme_plus_pag_calls.tsv")

usethis::use_data(geuvadis_info, 
		  gencode_chr_gene, gencode_chr_tx, 
		  gencode_pri_gene, gencode_pri_tx, 
		  gencode_all_gene, gencode_all_tx, 
		  gencode_hla, pag, polypheme_pag, 
		  overwrite = TRUE)

usethis::use_data(allele_hist, hla_groups, pag_groups, internal = TRUE, overwrite = TRUE)
