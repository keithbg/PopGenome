## PopGenome statistics on genes extracted from Snippy alignment

## Analyses data that was formatted in PopGenome_Snippy_Alignnment_Format.R,
## which extracted gene sequences from full alignment and then filtered them based on
## gap and N thresholds

## Snippy run on contigs
## biotite location: /data5/eelriver/CYA2/PH2017/snippy/snippy_out_contigs


#### POP GENOME INFORMATION ####################################################

## Valid sites in the alignment= sites without a gap or N, these sites are excluded from calculations

## Regions are each of the nucleotide files included in the folder input into readData
## Populations are organisms within each data file defined to be in differen populations

## Run a command such as diversity.stats(mydata)
## See output list contents get.diversity(mydata)
## See contents of each list get.diversity(mydata)[[1]]


run_PopGenome<- function(gene_path, pops= c(1, 2, 3), prop_valid_sites= 0.5, export_tables= FALSE){

  require(tidyverse)
  require(PopGenome)
  
  # Whole alignment FASTA
## Read in alignment of the genes
  message("Reading in fasta files to PopGenome")
mydata <- readData(path= gene_path, format="FASTA", SNP.DATA = TRUE, big.data= TRUE)


#### SET POPULATIONS (i.e. ANI defined species 1 to 4)
#individual.names <- sort(get.individuals(mydata)[[2]])
individual.names <- suppressMessages(read_tsv("/Users/kbg/Documents/UC_Berkeley/CyanoMeta_NSF/Metagenomics/Data/GenomesData/genome_ani_species.tsv")) %>% 
  select(genome) %>% 
  #filter(genome != "PH2015_01D_Oscillatoriales_44_513" | genome != "PH2015_02D_Oscillatoriales_45_1038") %>% 
  pull()


message("Assigning populations")

if(all(pops == c(1, 2, 3))){

populations <- list(c(individual.names[c(-1, -3, -c(12:19), -c(25:27), -c(37:39), -48, -49, -53)]),
                    c(individual.names[c(12, 14, 17, 25, 26, 27)]),
                    c(individual.names[c(13, 15:16, 18:19, 37:39, 48:49, 53)]))
}

if(all(pops == c(1, 2))){
  
  populations <- list(c(individual.names[c(-1, -3, -c(12:19), -c(25:27), -c(37:39), -48, -49, -53)]),
                      c(individual.names[c(12, 14, 17, 25, 26, 27)]))
}

if(all(pops == c(1, 3))){
  
  populations <- list(c(individual.names[c(-1, -3, -c(12:19), -c(25:27), -c(37:39), -48, -49, -53)]),
                      c(NA),
                      c(individual.names[c(13, 15:16, 18:19, 37:39, 48:49, 53)]))
}

if(all(pops == c(2, 3))){
  
  populations <- list(c(NA),
                      c(individual.names[c(12, 14, 17, 25, 26, 27)]),
                      c(individual.names[c(13, 15:16, 18:19, 37:39, 48:49, 53)]))
}




mydata <- set.populations(mydata, populations, diploid= FALSE)

mydata.sum <- as.data.frame(get.sum.data(mydata)) %>%
                mutate(geneID= str_replace(row.names(.), ".fa.aln", ""),
                       prop.valid.sites= round(n.valid.sites / n.sites, 4)) %>%
                select(geneID, everything()) %>%
                as_tibble()

#### CALCULATE STATISTICS 
# show.slots(mydata) # slots are the different type of analyses that can be conducted
message("Calculating population genetic statistics")
# Calculate Neutrality statistics
mydata <- neutrality.stats(mydata, detail= TRUE)

# Calculate F_ST and Diversity statistics
mydata <- F_ST.stats(mydata, mode= "nucleotide") # this also calculates diversity statistics

# Extract FST between populations
pairwise.FST <- t(mydata@nuc.F_ST.pairwise) %>%
  as_tibble() %>%
  #select(-`pop1/pop4`, -`pop2/pop4`, -`pop3/pop4`) %>%
  mutate(geneID= str_replace(mydata@region.names, ".fa.aln", "")) %>%
  #gather(key= pops, value= FST, `pop1/pop2`:`pop2/pop3`) %>%
  gather(key= pops, value= FST, contains("pop")) %>%
  arrange(geneID)


# Extract nucleotide diversity between and within populations
nucdiv.btw <- t(mydata@nuc.diversity.between / mydata.sum$n.sites) %>%
  as_tibble() %>%
  mutate(geneID= str_replace(mydata@region.names, ".fa.aln", "")) %>%
  #gather(key= pops, value= nuc.btw.div, `pop1/pop2`:`pop2/pop3`) %>%
  gather(key= pops, value= nuc.btw.div, contains("pop")) %>%
  arrange(geneID)

nucdiv.within <- (mydata@nuc.diversity.within / mydata.sum$n.sites) %>%
  as_tibble() %>%
  mutate(geneID= str_replace(mydata@region.names, ".fa.aln", "")) %>%
  #gather(key= pops, value= nuc.wtn.div, `pop 1`:`pop 3`) %>%
  gather(key= pops, value= nuc.wtn.div, contains("pop")) %>%
  arrange(geneID)

## Extract Tajima D
tajima <- mydata@Tajima.D %>%
  as_tibble() %>%
  mutate(geneID= str_replace(mydata@region.names, ".fa.aln", "")) %>%
  gather(key= pops, value= tajimaD, contains("pop")) %>%
  arrange(geneID)


#### Combine statistics into data frames
## Remove genes with valid sites < 50% of gene length
pop.stats.btw.df <- left_join(mydata.sum, pairwise.FST) %>% # BETWEEN STATISTICS
  left_join(., nucdiv.btw) %>%
  mutate(geneID= str_replace(.$geneID, ".fa", "")) %>%
  filter(prop.valid.sites > prop_valid_sites)


pop.stats.wtn.df <- left_join(mydata.sum, nucdiv.within) %>% # WITHIN STATISTICS
  left_join(., tajima) %>%
  mutate(geneID= str_replace(.$geneID, ".fa", "")) %>%
  filter(prop.valid.sites > prop_valid_sites)

#### MK Test ####
# message("calculating MK Test")
# mydata <- MKT(mydata, new.populations= populations,
#                  do.fisher.test= TRUE,
#                  fixed.threshold.fst= FALSE)
# 
# ## Extract results as list and add gene names
# mkt.list <- mapply(cbind, get.MKT(mydata), mydata@region.names)
# 
# ## Fill NULL results (where biallialec sites = 0) with a blank DF
# null.elements <- which(do.call(rbind, lapply(get.MKT(mydata), is.null)))
# mkt.list[null.elements] <- lapply(null.elements,
#                                   function(x)
#                                     mkt.list[x] <- cbind(matrix(NA, nrow= 3, ncol= 9), rep(mydata@region.names[x], 3))) # blank matrix of NAs, with geneID included
# 
# ## Transform list into a dataframe
# mkt.df <- do.call(rbind, mkt.list) %>% 
#   as.data.frame() %>% 
#   mutate(pops = str_replace(row.names(.), "\\.[0-9].*$", "")) %>% 
#   #rename("geneID" = V10) %>% 
#   as_tibble() %>% 
#   mutate_at(vars(P1_nonsyn:fisher.P.value), funs(as.numeric(as.character(.)))) %>% 
#   #mutate(geneID = str_replace(as.character(geneID), ".fa.aln", "")) 
#   #left_join(., feature.annotation, by= "geneID") # combine with annotation information
# 
# '
## MK Test significant values p < 0.05
# mkt.sig <- mkt.df %>% 
#   filter(fisher.P.value < 0.05) %>% 
#   arrange(fisher.P.value)
# write_tsv(mkt.sig, path= file.path(export_tables, "MK_test_sig.txt"))

#### EXPORT OUTLIER TABLES 

if(export_tables != FALSE){
  message("Exporting tables")
# Fst
fst.outlier <- pop.stats.btw.df %>% 
  filter(FST < 0.5) %>% 
  mutate(outlier= "fst")
  write_tsv(fst.outlier, file.path(export_tables, "Fst_low.txt"))

# Intra-specific
intra.outlier <- pop.stats.wtn.df %>% 
  filter(nuc.wtn.div > 0.01) %>% 
  mutate(outlier= "intra")
  write_tsv(intra.outlier, file.path(export_tables, "Intra_specific_high.txt"))

# Inter-specific
inter.outlier <- pop.stats.btw.df %>% 
  filter(nuc.btw.div > 0.25) %>% 
  mutate(outlier= "inter")
  write_tsv(inter.outlier, file.path(export_tables, "Inter_specific_high.txt"))

# Tajima's D
tajima.outlier <- pop.stats.wtn.df %>% 
  filter(tajimaD > 2 | tajimaD < -2) %>% 
  mutate(outlier= "tajima")
  write_tsv(tajima.outlier, file.path(export_tables, "Tajima_2.txt"))

  
 # Read in MK Test results
 # mkt <- read_tsv(file.path(dir_output_table, "MK_test_sig.txt")) %>% 
 #   select(geneID, feature, pops, neutrality.index, alpha, fisher.P.value, anno) %>% 
 #   rename(attribute= feature) %>% 
 #   mutate(outlier= "MKT",
 #          pops= str_replace(.$pops, "\\.", "\\/")) %>% 
 #   separate(anno, into= c("uniref_anno", "uniprot_anno", "kegg_anno"), sep= " -- ")
 
 
   
 
outlier.master <- full_join(fst.outlier, intra.outlier) %>% 
                     full_join(., inter.outlier) %>% 
                     full_join(., tajima.outlier) %>% 
                     full_join(., mkt) %>% 
                     select(geneID, attribute, outlier, pops, everything()) %>% 
                     arrange(geneID)
write_tsv(outlier.master, file.path(export_tables, "PopGenome_outlier_genes.txt"))
}


#### RETURN DATA FRAMES 

output_list <- list(mydata.sum, pop.stats.wtn.df, pop.stats.btw.df, mydata)
names(output_list) <- c("mydata.sum", "pop.stats.wtn.df", "pop.stats.btw.df", "mydata")
return(output_list)
}


test <- run_PopGenome(gene_path = dir_core_genes_sp12, pops= c(1, 2))

names(test)


head(test[["mydata.sum"]])
mkt.test <- test[["mkt.df"]]




#### EXTRACT CORE GENOME ALIGNMENTS FROM ROARY OUTPUT #########################

## Extract the gene alignments of the core genes from the Roary output
## Must run Roary with the -z flag: "don't delete intermediate files"

## Gene alignments are in "pan_genome_sequences" folder

## path= folder of Roary Output
## core.prop= Proportion of genomes with gene to be considered core (0-1)

extract_core_aln <- function(roary.path, core.prop){
  require(tidyverse)
  
  ## Read in gene presence/absence table
  gene.df <- read_tsv(file.path(roary.path, "gene_presence_absence.Rtab")) %>% 
    mutate(gene_counts= rowSums(.[, -1]))
  
  ## Count genomes in table
  genome.count <- sum(str_detect(names(gene.df), "PH"))
  
  ## Filter genes in the core
  gene.core <- gene.df %>% 
    filter(gene_counts >= ceiling(genome.count * core.prop)) ## Proportion of genomes considered to be core
  
  core.genes <- gene.core %>% 
    select(Gene) %>%
    pull()
  
  ## List gene alignment sequences output by Roary
  gene.sequences <- list.files(file.path(roary.path, "pan_genome_sequences"))
  
  ## Subset the core sequences
  core.sequences <- gene.sequences[str_replace(gene.sequences, ".fa.aln", "") %in% core.genes]
  
  ## Move the core sequences to a new folder, "core_genome_sequences"
  dir.create(file.path(roary.path, "core_genome_sequences"))
  file.copy(from= file.path(roary.path, "pan_genome_sequences", core.sequences),
            to= file.path(roary.path, "core_genome_sequences"))
  
  return(gene.core)
}


getwd()
core_genes <- extract_core_aln(roary.path= "Output_PHall_bp95_c50", core.prop = 0.2)


genome_species <- read_tsv("/Users/kbg/Documents/UC_Berkeley/CyanoMeta_NSF/Metagenomics/Data/GenomesData/genome_ani_species.tsv")

core_list <- list()

for(species in 1:3){
  sp_list <- genome_species %>% 
    filter(ani_species == species)
  sp_genomes <- names(core_genes)[names(core_genes) %in% sp_list$genome]
  
  sp_genes <- core_genes %>% 
    select(Gene, sp_genomes) %>% 
    mutate(gene_counts= rowSums(.[, -1]))
  
  sp_core_genes <- sp_genes %>% 
    filter(gene_counts >= ceiling(length(sp_genomes) * 0.5))
  
  core_list[[species]] <- sp_core_genes
  
  if(species == 3){
    names(core_list) <- c("sp1", "sp2", "sp3")  
  }
}  


## Find overlapping gene names in each of the 3 species
table(core_list[["sp1"]]$Gene %in% core_list[["sp2"]]$Gene)
table(core_list[["sp1"]]$Gene %in% core_list[["sp3"]]$Gene)
table(core_list[["sp2"]]$Gene %in% core_list[["sp3"]]$Gene)


core_genes_A <- core_list[["sp1"]]$Gene[core_list[["sp1"]]$Gene %in% core_list[["sp2"]]$Gene]



dir.create(file.path(roary.path, "core_genome_sequences_sp1-2"))
file.copy(from= file.path(roary.path, "core_genome_sequences", paste0(core_genes_A, ".fa.aln")),
          to= file.path(roary.path, "core_genome_sequences_sp1-2"))

core_genes_B <- core_genes_A[core_genes_A %in% core_list[["sp3"]]$Gene]






