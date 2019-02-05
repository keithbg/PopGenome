## McDonald_Kreitman Test from PopGenome
## McDonald and Kreitman, Nature, 1991

#### Libraries #################################################################
library(tidyverse)
library(PopGenome)
library(ggplot2)
################################################################################

#### FILE PATHS ################################################################
dir_input <- file.path("/Users","kbg","Documents","UC_Berkeley","CyanoMeta_NSF","Metagenomics", "Data","GenomesData", "PopGenome", "Alignments_snippy", "Alignments_snippy_format")
dir_input_annotations <- file.path("/Users","kbg","Documents","UC_Berkeley","CyanoMeta_NSF","Metagenomics", "Data","GenomesData", "PopGenome", "Annotations")
dir_output_figures <- file.path("/Users","kbg","Documents","UC_Berkeley","CyanoMeta_NSF","Metagenomics", "Data","GenomesData", "PopGenome", "Output_figures")
dir_output_table <- file.path("/Users","kbg","Documents","UC_Berkeley", "CyanoMeta_NSF","Metagenomics", "Data","GenomesData", "PopGenome", "Output_tables")
################################################################################


#### Make GFF files for each gene that passes filtering thresholds
source("/Users/kbg/Documents/UC_Berkeley/CyanoMeta_NSF/Metagenomics/Data/GenomesData/PopGenome/R_scripts/PopGenome_Make_Gene_GFF_Files.R")
#make_gff_files(write_table= TRUE)


#### Read in FASTA and GFF files to PopGenome and get gene annotations
mydataGFF <- readData(path= file.path(dir_input, "gene_ggkbase_fasta_files"), format="FASTA", gffpath= file.path(dir_input, "GFF"))

gff <- read.delim(file.path(dir_input_annotations, "Reference_ggkbase.gff"), header=F, comment.char="#", stringsAsFactors = FALSE) %>%
  as_tibble() %>%
  rename(seqid= V1, source= V2, type= V3, start= V4, end= V5, score= V6, strand= V7, phase= V8, attribute= V9) %>%
  mutate(geneID= str_c("gene_", str_pad(seq(1:nrow(.)), pad= "0", width= 4)),
         start= as.numeric(start),
         end= as.numeric(end)) %>% 
  select(geneID, attribute) %>% 
  rename(feature= attribute) 

feature.annotation <- read_tsv(file.path(dir_input_annotations, "Reference_ggkbase_gff_subset.txt")) %>% 
  select(geneID, attribute) %>% 
  rename(anno= attribute) %>% 
  inner_join(., gff, by= "geneID")


## Check GFF data in PopGenome
# mydataGFF@Coding.region
# get.codons(mydataGFF, 735) # gene 1 in algn@region.names


#### Set populations (i.e. ANI defined species 1, 2, and 3) ####
individual.names <- sort(get.individuals(mydataGFF)[[4]]) # use region #4 because that has all samples in the FASTA

populations <- list(c(individual.names[c(-10, -12, -15, -23, -24, -25, -11, -13, -14, -16, -17, -35, -51, -36, -37, -47, -46)]),
                    c(individual.names[c(10, 12, 15, 23, 24, 25)]),
                    c(individual.names[c(11, 13, 14, 16, 17, 35, 51, 36, 37, 47, 46)]))

mydataGFF <- set.populations(mydataGFF, populations, diploid= FALSE)

#### MK Test ####
mydataGFF <- MKT(mydataGFF, new.populations= populations,
                do.fisher.test= TRUE,
                fixed.threshold.fst=FALSE)

## Extract results as list and add gene names
mkt.list <- mapply(cbind, get.MKT(mydataGFF), mydataGFF@region.names)

## Fill NULL results (where biallialec sites = 0) with a blank DF
null.elements <- which(do.call(rbind, lapply(get.MKT(mydataGFF), is.null)))
mkt.list[null.elements] <- lapply(null.elements,
                                  function(x)
                                    mkt.list[x] <- cbind(matrix(NA, nrow= 3, ncol= 9), rep(mydataGFF@region.names[x], 3))) # blank matrix of NAs, with geneID included

## Transform list into a dataframe
mkt.df <- do.call(rbind, mkt.list) %>% 
  as.data.frame() %>% 
  mutate(pops = str_replace(row.names(.), "\\.[0-9].*$", "")) %>% 
  rename("geneID" = V10) %>% 
  as_tibble() %>% 
  mutate_at(vars(P1_nonsyn:fisher.P.value), funs(as.numeric(as.character(.)))) %>% 
  mutate(geneID = str_replace(as.character(geneID), "_[+|-]$", "")) %>% 
  left_join(., feature.annotation, by= "geneID") # combine with annotation information


## MK Test significant values p < 0.05
mkt.sig <- mkt.df %>% 
  filter(fisher.P.value < 0.05) %>% 
  arrange(fisher.P.value)
write_tsv(mkt.sig, path= file.path(dir_output_table, "MK_test_sig.txt"))


#### Make Plots
#ggplot(mkt.df, aes(x= pops, y= neutrality.index)) +
#  geom_boxplot() +
#  facet_grid(.~fisher.P.value < 0.05)



