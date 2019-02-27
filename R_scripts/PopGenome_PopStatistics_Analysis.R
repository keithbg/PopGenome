## Analyze the outlier statistics from PopGenome analysese
## Fst < 0.05, 
## Intra-specific diversity > 0.01
## Inter-specific diversity > 0.25
## MK Test p value < 0.05

#### Libraries #################################################################
library(tidyverse)
library(PopGenome)
library(ggplot2)
################################################################################

#### FILE PATHS ################################################################
dir_input <- file.path("/Users","kbg","Documents","UC_Berkeley","CyanoMeta_NSF","Metagenomics", "Data","GenomesData", "PopGenome", "Alignments_snippy", "Alignments_snippy_format")
dir_input_annotations <- file.path("/Users","kbg","Documents","UC_Berkeley","CyanoMeta_NSF","Metagenomics", "Data","GenomesData", "PopGenome", "Annotations")
dir_output_figures <- file.path("/Users","kbg","Documents","UC_Berkeley","CyanoMeta_NSF","Metagenomics", "Data","GenomesData", "PopGenome", "Output_figures")
dir_input_table <- file.path("/Users","kbg","Documents","UC_Berkeley","CyanoMeta_NSF","Metagenomics", "Data","GenomesData", "PopGenome", "Output_tables")
################################################################################


stats <- read_tsv(file.path(dir_input_table, "PopGenome_outlier_genes.txt")) %>% 
  select(geneID, attribute, scaffold, scaff_length, outlier, pops, FST, nuc.btw.div, nuc.wtn.div, tajimaD, neutrality.index, alpha, fisher.P.value, uniref_anno, uniprot_anno, kegg_anno)

stats %>% 
  count(attribute) %>% 
  arrange(-n)


stats.anno <- stats %>% 
  count(uniref_anno) %>% 
  arrange(-n)
