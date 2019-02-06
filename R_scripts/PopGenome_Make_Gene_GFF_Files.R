## Make a GFF file for each gene to be used by PopGenome
## All the GFF files will be in a single folder

make_gff_files <- function(write_table= FALSE){
#### Libraries #################################################################
require(tidyverse)
require(seqinr)
################################################################################

#### FILE PATHS ################################################################
dir_input <- file.path("/Users","kbg","Documents","UC_Berkeley","CyanoMeta_NSF","Metagenomics", "Data","GenomesData", "PopGenome", "Alignments_snippy", "Alignments_snippy_format")
dir_input_annotations <- file.path("/Users","kbg","Documents","UC_Berkeley","CyanoMeta_NSF","Metagenomics", "Data","GenomesData", "PopGenome", "Annotations")
dir_output_GFF <- file.path("/Users","kbg","Documents","UC_Berkeley","CyanoMeta_NSF","Metagenomics", "Data","GenomesData", "PopGenome", "Alignments_snippy", "Alignments_snippy_format", "GFF")
################################################################################


message("Writing GFF files for each gene to:\n", dir_output_GFF)

## Parse gene annotations
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

# Extract only attribute and annotations for future merge with GFF 
  # Annotation order: uniref, uniprot, kegg (separated by --)
 anno.for.gff <- anno.ref %>%
   mutate(anno= str_c(uniref_anno, uniprot_anno, kegg_anno, sep=" -- ")) %>% 
   select(attribute, anno) %>% 
   mutate(anno = ifelse(is.na(anno), "NA -- NA -- NA", anno),
          anno = str_replace_all(anno, "'", "")) # remove apostraphes (') otherwise PopGenome throws an error when reading GFF

## Read GFF file with start/end locations of genes
gff <- read.delim(file.path(dir_input_annotations, "Reference_ggkbase.gff"), header=F, comment.char="#", stringsAsFactors = FALSE) %>%
  as_tibble() %>%
  rename(seqid= V1, source= V2, type= V3, start= V4, end= V5, score= V6, strand= V7, phase= V8, attribute= V9) %>%
  mutate(geneID= str_c("gene_", str_pad(seq(1:nrow(.)), pad= "0", width= 4)),
         start= as.numeric(start),
         end= as.numeric(end))

# Get genes that passed filtering thresholds and are used in analyses
  # Filtering in PopGenome_Snippy_Alignment_Format.R
gene.list <- str_replace(list.files(file.path(dir_input, "gene_ggkbase_fasta_files")), "_[+|-]", "")

## Subset GFF file by genes that passed filtering thresholds and are used in analyses, 
## And merge with annotation df
gff.subset <- gff[gff$geneID %in% gene.list, ] %>% 
  left_join(., anno.for.gff, by= "attribute") %>% 
  select(seqid:phase, anno, geneID)  %>% 
  rename(attribute = anno)

if(write_table == TRUE){
  write_tsv(gff.subset, path= file.path(dir_input_annotations, "Reference_ggkbase_gff_subset.txt"), col_names = TRUE)
}

# Loop to extract GFF row for each gene and write a file
gene.files <- str_replace(list.files(file.path(dir_input, "gene_ggkbase_fasta_files")), ".fa", "")
gene.id <- str_replace(gene.files, "_[+|-]", "")

for(gene in gene.id){
  # Extract gene names
  gff.file <- gene.files[str_detect(gene.files, gene)]
  
  ## Subset GFF by geneID
  gff.gene <- gff.subset %>% 
    filter(geneID == gene) %>% 
    mutate(seqid = gene,
           end= (end - start + 1),
           start= 1) %>% 
    select(-geneID) %>% 
    #write_tsv(., path= file.path(dir_output_GFF, str_c("gff_", gff.file)), col_names = FALSE, append= FALSE) %>% 
    write_tsv(., path= file.path(dir_output_GFF, gff.file), col_names = FALSE, append= FALSE)
}


}
