## Subset of Snippy Core alignment

## Snippy run on contigs
## biotite location: /data5/eelriver/CYA2/PH2017/snippy/snippy_out_contigs

## Many parts of the Snippy alignment had gaps in the individual genomes
## which makes certain sliding windows non-sensicle if genomes and ANI species are dominated by gaps

## This script takes gene predictions from Prokka and extracts genes that have
## <25% gaps& Ns in >50% of the genomes in each of the ANI species
## ANI species 4 is not considered because it only has 2 genomes
## Also genes were filtered to a minimum length of 225 nucleotides (75 amino acids)

##


#### Libraries #################################################################
library(tidyverse)
library(seqinr)
################################################################################

#### FILE PATHS ################################################################
dir_input_snippy <- file.path("/Users","kbg","Documents","UC_Berkeley","CyanoMeta_NSF","Metagenomics", "Data","GenomesData", "PopGenome", "Alignments_snippy")
dir_input_annotations <- file.path("/Users","kbg","Documents","UC_Berkeley","CyanoMeta_NSF","Metagenomics", "Data","GenomesData", "PopGenome", "Annotations")
dir_output_snippy <- file.path("/Users","kbg","Documents","UC_Berkeley","CyanoMeta_NSF","Metagenomics", "Data","GenomesData", "PopGenome", "Alignments_snippy", "Alignments_snippy_format")
dir_output_table <- file.path("/Users","kbg","Documents","UC_Berkeley","CyanoMeta_NSF","Metagenomics", "Data","GenomesData", "PopGenome", "Output_tables")
################################################################################


#### READ IN DATA #############################################################

## WHOLE ALIGNMENT FASTA FROM SNIPPY
alg <- read.alignment(file.path(dir_input_snippy, "core.full.aln.LO"), format= "fasta")


## GFF FILE
# gff <- read.delim(file.path(dir_input_snippy, "Reference_ggkbase.gff"), header=F, comment.char="#", stringsAsFactors = FALSE) %>%
#           slice(-seq(which(str_detect(.$V1, ">") == TRUE), nrow(.), by= 1)) %>%  # Remove fasta file sequence at end of gff file
#           as_tibble() %>%
#           rename(seqid= V1, source= V2, type= V3, start= V4, end= V5, score= V6, strand= V7, phase= V8, attribute= V9) %>%
#           filter(source == "prokka") %>%
#           mutate(geneID= str_c("gene_", str_pad(seq(1:nrow(.)), pad= "0", width= 4)))

## GGKBASE ANNOTATIONS
gff <- read.delim(file.path(dir_input_annotations, "Reference_ggkbase.gff"), header=F, comment.char="#", stringsAsFactors = FALSE) %>%
  as_tibble() %>%
  rename(seqid= V1, source= V2, type= V3, start= V4, end= V5, score= V6, strand= V7, phase= V8, attribute= V9) %>%
  mutate(geneID= str_c("gene_", str_pad(seq(1:nrow(.)), pad= "0", width= 4)),
         start= as.numeric(start),
         end= as.numeric(end))

## CALCULATE THE FRACTION OF GAPS IN EACH GENE FOR EACH GENOME ##############S##
## Takes 3-4 hours to run

gap_calcs <- function(){
## Initialize empty data frame for gene gap information
 gene.info <- data.frame(matrix(NA, nrow= length(alg$seq)*nrow(gff), ncol= 8)) %>%
                rename(genome= X1, gene= X2, length= X3, gaps= X4, gap_fraction= X5, Ns= X6, N_fraction= X7, prop_valid_sites= X8) %>%
               as_tibble()
 counter <- 1
 
 ## Run loop
 for(genome in seq(1:length(alg$seq))){
   svMisc::progress(genome, max.value= length(alg$seq), progress.bar = TRUE)

   ## Extract genome sequence and name
   seq <- alg$seq[genome]
   genome.name <- alg$nam[genome]

   for(geneID in seq(1:nrow(gff))){
     ## Extract gene from genome
     gene <- str_sub(seq, start= pull(gff[geneID, "start"]), end= pull(gff[geneID, "end"]))
     gene.number <- str_c("gene_", str_pad(geneID, pad= "0", width= 4))

     ## Calculate fraction of gaps
     gene.length <- str_length(gene)
     gap.count <- str_count(gene, "-")
     gap.fraction <- round(gap.count / gene.length, 5)

     ## Calculate fraction of Ns
     n.count <- str_count(gene, "n")
     n.fraction <- round(n.count / gene.length, 5)

     ## Proportion of valid sites
     prop.valid.sites <- 1 - round((gap.count + n.count) / gene.length, 5)

     gene.info[counter, ] <- c(alg$nam[genome], gene.number, gene.length, gap.count, gap.fraction, n.count, n.fraction, prop.valid.sites)
     counter <- counter + 1
   }
   write_tsv(gene.info, file.path(dir_output_snippy, "gene_gaps_Ns_ggkbase.tsv"), col_names= TRUE, append= FALSE)
}
}

#gene_gaps <- gap_calcs
#### FILTER GENES ####################################################################

#### Filter genes to only ones that have good coverage from all genomes in a species

gaps.df <- read_tsv(file.path(dir_output_snippy, "gene_gaps_Ns_ggkbase.tsv")) %>%
            mutate(genome= str_replace(.$genome, "snippy_out_LO_ctg_", "")) %>%
            mutate(genome= str_replace(.$genome, "_s25", "")) %>%
            #mutate(prop_valid_sites= 1-prop_valid_sites) %>%
            left_join(., read_tsv("/Users/kbg/Documents/UC_Berkeley/CyanoMeta_NSF/Metagenomics/Data/GenomesData/genome_ani_species.tsv")) %>%  # add ANI species information for each genome
            filter(ani_species != 4)

valid.sites.threshold <- 0.75 # >= 75% of sites do not have a gap or an N

#
# prop.table(table(gaps.df$prop_valid_sites > 0.5))
# ggplot(data= gaps.df) +
#   geom_histogram(aes(x= N_fraction))
# ggplot(data= gaps.df) +
#   geom_histogram(aes(x= gap_fraction))
# ggplot(data= gaps.df) +
#   geom_histogram(aes(x= prop_valid_sites, fill= prop_valid_sites > 0.75))


## Counts and proportion of genomes in each ANI species with the gene
## ANI species 1= 47 genomes, 2= 6 genomes, 3= 11 genomes, 4= 2 genomes
genes.filtered <- gaps.df %>%
                    filter(prop_valid_sites >= valid.sites.threshold) %>%
                    dplyr::count(ani_species, gene) %>%
                    mutate(n.frac= ifelse(ani_species == 1, n / 47,
                                          ifelse(ani_species == 2, n / 6,
                                                 ifelse(ani_species == 3, n / 11, n))))

# ggplot(data= genes.filtered) +
#   geom_histogram(aes(x= n), binwidth = 1, fill= "gray75", color= "black") +
#   facet_grid(.~ani_species, scales= "free_x") +
#   theme_bw()

## Identify genes that are present in >= 50% of genomes in each species
## Since Species 2 has only 6 genomes, this means all genes are present in at least 3 genomes / species
genes.filtered.complete <- genes.filtered %>%
                        select(ani_species, gene, n.frac) %>%
                        spread(key= ani_species, value= n.frac) %>%
                        filter(complete.cases(.)) %>%
                        filter(`1` >= 0.5 & `2` >= 0.5 & `3` >= 0.5) # gene present >50% of genomes in each species

# genes.filtered.complete %>%
#   gather(key= ani_species, value= n.frac, `1`:`3`) %>%
#   ggplot() +
#   geom_histogram(aes(x= n.frac), binwidth = 0.05) +
#   facet_grid(.~ani_species, scales= "free_x") +
#   theme_bw()

## Extract genes from gff table that correspond to the criteria above
## Ignore genes shorter than 75 amino acids (225 nucleotides)
gff.subset <- gff %>%
  slice(as.numeric(str_replace(genes.filtered.complete$gene, "gene_", ""))) %>%
  mutate(gene_length= end - start + 1) %>%
  filter(gene_length > 225)


#### WRITE FASTA FILES FOR EACH GENE ################################################ 
## Function to reverse and complement DNA sequence
reverse_comp_dna <- function(string){
  # split string by characters
  string_split <-  strsplit(string, split = "")
  # reverse order
  rev_string <- rev(string_split[[1]])
  # complement characters
  comp_rev_string <-  seqinr::comp(rev_string)
  # Restore the NAs back to gaps
  comp_rev_string <- ifelse(is.na(comp_rev_string), "-", comp_rev_string)
  # collapse reversed characters
  return(paste(comp_rev_string, collapse = ""))
} 

  ## Change genome names
alg$nam <- str_replace(alg$nam, "snippy_out_LO_ctg_", "")
alg$nam <- str_replace(alg$nam, "_s25", "")

## Remove ANI species 4 from export
analysis.genomes <- which(alg$nam %in% unique(gaps.df$genome))
analysis.genomes <- c(20, analysis.genomes[-18]) # reorder, so that reference genome, PH2015_09S_Oscillatoriales_45_247, is at the top of each fasta file

for(row in seq(1:nrow(gff.subset))){
  svMisc::progress(row, max.value= nrow(gff.subset), progress.bar = TRUE)

  ## Initialize empty list to store seauences for all genomes for each gene
  gene.list <- rep(list(NULL), length(analysis.genomes))

  counter <- 1
  for(genome in analysis.genomes){
    ## Extract genome
    seq <- alg$seq[genome]
    genome.name <- alg$nam[genome]

    ## Extract gene from genome
    gene <- str_sub(seq, start= pull(gff.subset[row, "start"]), end= pull(gff.subset[row, "end"]))
    gene.name <- gff.subset[row, "geneID"]
    gene.strand <- gff.subset[row, "strand"]

    # Gene on reverse strand, reverse sequence and change to complement nucleotides
    if(gff.subset[row, "strand"] == "-"){
      gene <- reverse_comp_dna(gene)
    }
    ## Calculate gap and N fraction
    gap.count <- str_count(gene, "-")
    n.count <- str_count(gene, "n")
    
    ## Proportion of valid sites
    prop.valid.sites <- 1 - round((gap.count + n.count) / str_length(gene), 5)

    ## Add sequence to list, IF valid sites are >75% of gene length
    if(prop.valid.sites >= valid.sites.threshold){
      gene.list[counter] <- gene
      names(gene.list)[counter] <- genome.name
    }
    counter <- counter + 1
  }
  
  ## Remove list elements that did not meet the gap threshold
  gene.list <- gene.list[!is.na(names(gene.list))]

  ## Write fasta including all genomes for each gene
  write.fasta(gene.list, file.out= file.path(dir_output_snippy, "gene_ggkbase_fasta_files", str_c(gene.name,"_", gene.strand, ".fa")), names= names(gene.list), nbchar= 80, as.string= TRUE)
 rm(seq, gene, genome.name, counter)
}
rm(analysis.genomes, gene.name, gene.list)


