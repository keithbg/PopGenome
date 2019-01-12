## Make a GFF file for each gene to be used by PopGenome
## All the GFF files will be in a single folder

#### Libraries #################################################################
library(tidyverse)
library(PopGenome)
library(ggplot2)
################################################################################

#### FILE PATHS ################################################################
dir_input <- file.path("/Users","kbg","Documents","UC_Berkeley","CyanoMeta_NSF","Metagenomics", "Data","GenomesData", "PopGenome", "Alignments_snippy", "Alignments_snippy_format")
dir_input_annotations <- file.path("/Users","kbg","Documents","UC_Berkeley","CyanoMeta_NSF","Metagenomics", "Data","GenomesData", "PopGenome", "Annotations")
dir_output_table <- file.path("/Users","kbg","Documents","UC_Berkeley","CyanoMeta_NSF","Metagenomics", "Data","GenomesData", "PopGenome", "Output_tables")
################################################################################



# Whole alignment FASTA
## Read in alignment of the genes
#mydata2 <- readData(path= file.path(dir_input, "gene_ggkbase_fasta_files"), format="FASTA", SNP.DATA = TRUE, big.data= TRUE)


## Start/end locations of each gene in GFF file
gff <- read.delim(file.path(dir_input_annotations, "Reference_ggkbase.gff"), header=F, comment.char="#", stringsAsFactors = FALSE) %>%
  as_tibble() %>%
  rename(seqid= V1, source= V2, type= V3, start= V4, end= V5, score= V6, strand= V7, phase= V8, attribute= V9) %>%
  mutate(geneID= str_c("gene_", str_pad(seq(1:nrow(.)), pad= "0", width= 4)),
         start= as.numeric(start),
         end= as.numeric(end))

gene.list <- str_replace(list.files(file.path(dir_input, "gene_ggkbase_fasta_files")), ".fa", "")

gff.subset <- gff[gff$geneID %in% gene.list, ] 

gff.subset %>% 
  rename(feature= "type") %>% 
  mutate(feature= "gene",
         attribute= str_c(geneID, attribute, sep= " ")) %>% 
  select(-geneID)
  
  

## Gene annotations
clean.up.annotation.regex <- " Tax.*$| tax.*$| evalue.*$| \\{.*$| \\(.*$|,.*$| bin.*$| n=."
anno.ref <- read.delim(file.path(dir_input_annotations, "PH2015_09S_Oscillatoriales_45_247.ql"), header=F, comment.char="#", stringsAsFactors = FALSE) %>%
  as_tibble() %>%
  select(-V5, -V6, -V7, -V8, -V11, -V13, -V14) %>%
  rename(attribute= V1, scaffold= V2, scaff_feature= V3, scaff_length= V4, start= V9, end= V10, annotation= V12) %>%
  separate(annotation, into= c("uniref_all", "uniprot_all", "kegg_all"), sep= "__ ") %>%
  mutate(uniref_anno= tolower(str_replace(.$uniref_all, clean.up.annotation.regex, "")),
         uniprot_anno= tolower(str_replace(.$uniprot_all, clean.up.annotation.regex, "")),
         kegg_anno= tolower(str_replace(.$kegg_all, clean.up.annotation.regex, ""))) %>%
  mutate_at(vars(uniref_anno:kegg_anno), funs(trimws(., which= "both"))) %>%
  mutate_at(vars(uniref_anno:kegg_anno), funs(str_replace_all(., "uncharacterized protein|protein of unknown function.*|putative uncharacterized.*|hypothetical protein", "unknown protein"))) %>%
  mutate_at(vars(uniref_anno:kegg_anno), funs(ifelse(. == "tax=cg_cyano_01", NA, .))) %>%
  mutate_at(vars(uniref_anno:kegg_anno), funs(str_replace(., " tax.*$", ""))) %>%
  select(attribute, scaffold, uniref_anno, uniprot_anno, kegg_anno)
