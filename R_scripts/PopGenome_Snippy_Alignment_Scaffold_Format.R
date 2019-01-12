## Extract scaffolds from Snippy alignmnent

## Snippy run on contigs
## biotite location: /data5/eelriver/CYA2/PH2017/snippy/snippy_out_contigs

## This script extracts each of the scaffolds from the reference genome (PH2015_09S_Oscillatoriales_45_247)
## from the Snippy alignment, so they can be used for downstream analyses with the the PopGenome package

#### Libraries #################################################################
library(tidyverse)
library(seqinr)
################################################################################

#### FILE PATHS ################################################################
dir_input_snippy <- file.path("/Users","kbg","Documents","UC_Berkeley","CyanoMeta_NSF","Metagenomics", "Data","GenomesData", "PopGenome", "Alignments_snippy")
dir_output_snippy <- file.path("/Users","kbg","Documents","UC_Berkeley","CyanoMeta_NSF","Metagenomics", "Data","GenomesData", "PopGenome", "Alignments_snippy", "Alignments_snippy_format")
dir_output_table <- file.path("/Users","kbg","Documents","UC_Berkeley","CyanoMeta_NSF","Metagenomics", "Data","GenomesData", "PopGenome", "Output_tables")
################################################################################

# Whole alignment FASTA
alg <- read.alignment(file.path(dir_input_snippy, "core.full.aln.LO"), format= "fasta")
#alg <- read.fasta(file.path(dir_input, "core.full.aln.LO"), forceDNAtolower = FALSE)

## Scaffold lengths
scaff.df <- read.delim(file.path(dir_input_snippy, "Alignments_snippy_format", "Reference_scaffold_lengths.gff"), header=F, comment.char="#", stringsAsFactors = FALSE) %>% 
  as_tibble() %>% 
  select(V9, V4, V5) %>% 
  rename("scaffold"= V9, "start"= V4, "end"= V5) %>% 
  mutate(scaff_length= end - start + 1) %>% 
  arrange(-scaff_length, start)


## Change genome names
alg$nam <- str_replace(alg$nam, "snippy_out_LO_ctg_", "")
alg$nam <- str_replace(alg$nam, "_s25", "")
names(alg$seq) <- alg$nam

## Remove ANI species 4 from export
analysis.genomes <- alg$nam[!(alg$nam %in% c("PH2015_01D_Oscillatoriales_44_513", "PH2015_01U_Oscillatoriales_44_212"))]
analysis.genomes <- analysis.genomes[c(18, 1:17, 19:64)] # reorder, so that reference genome, PH2015_09S_Oscillatoriales_45_247, is at the top of each fasta file



## WRITE FASTA FILES FOR EACH SCAFFOLD

for(row in seq(1:nrow(scaff.df))){
  svMisc::progress(row, max.value= nrow(scaff.df), progress.bar = TRUE)
  
  ## Initialize empty list to store seauences for all genomes for each scaffold
  scaff.list <- rep(list(NULL), length(analysis.genomes))
  
  counter <- 1
  for(genome in analysis.genomes){
    ## Extract genome
    seq <- alg$seq[genome]
    genome.name <- alg$nam[alg$nam %in% genome]
  
    ## Extract scaffold sequence from genome
    scaffold_seq <- str_sub(seq, start= pull(scaff.df[row, "start"]), end= pull(scaff.df[row, "end"]))
    scaffold_name <- scaff.df[row, "scaffold"]
    
    ## Add scaffold sequence to list
    scaff.list[counter] <- scaffold_seq
    names(scaff.list)[counter] <- genome
    
    counter <- counter + 1
  }
  
  ## Write fasta including all genomes for each gene
  write.fasta(scaff.list, file.out= file.path(dir_output_snippy, "scaff_fasta_files", str_c(scaffold_name, ".fa")), names= names(scaff.list), nbchar= 80, as.string= TRUE)
  rm(seq, scaffold_seq, genome.name, counter)
}
#rm(algs)


