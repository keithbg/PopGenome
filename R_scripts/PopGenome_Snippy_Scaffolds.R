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

## Run a command such as diversity.stats(scaff.data)
## See output list contents get.diversity(scaff.data)
## See contents of each list get.diversity(scaff.data)[[1]]

#### Libraries #################################################################
library(tidyverse)
library(PopGenome)
library(ggplot2)
################################################################################

#### FILE PATHS ################################################################
dir_input_outliers <- file.path("/Users","kbg","Documents","UC_Berkeley","CyanoMeta_NSF","Metagenomics", "Data","GenomesData", "PopGenome", "Output_tables")
dir_input_scaff <- file.path("/Users","kbg","Documents","UC_Berkeley","CyanoMeta_NSF","Metagenomics", "Data","GenomesData", "PopGenome", "Alignments_snippy", "Alignments_snippy_format")
dir_input_annotations <- file.path("/Users","kbg","Documents","UC_Berkeley","CyanoMeta_NSF","Metagenomics", "Data","GenomesData", "PopGenome", "Annotations")
dir_output_figures <- file.path("/Users","kbg","Documents","UC_Berkeley","CyanoMeta_NSF","Metagenomics", "Data","GenomesData", "PopGenome", "Output_figures")
dir_output_table <- file.path("/Users","kbg","Documents","UC_Berkeley","CyanoMeta_NSF","Metagenomics", "Data","GenomesData", "PopGenome", "Output_tables")
################################################################################

## Scaffold lengths
scaff.df <- read.delim(file.path(dir_input_scaff, "Reference_scaffold_lengths.gff"), header=F, comment.char="#", stringsAsFactors = FALSE) %>% 
  as_tibble() %>% 
  select(V9, V4, V5) %>% 
  rename("scaffold"= V9, "scaf_start"= V4, "scaf_end"= V5) %>% 
  mutate(scaff_length= scaf_end - scaf_start + 1)


## Outlier information
out.df <- read_tsv(file.path(dir_input_outliers, "PopGenome_outlier_genes.txt")) %>% 
           rename(gene_start= "start", gene_end= "end") %>% 
           left_join(., scaff.df)

out.df %>% 
  dplyr::count(scaffold) 

## Remove scaffolds less than 50kb
out.df.filt <- out.df %>% 
  filter(scaff_length > 49000)

## List fasta scaffold fasta files to move to a new folder
file.list <- list.files(file.path(dir_input, "scaff_fasta_files"))[list.files(file.path(dir_input, "scaff_fasta_files")) %in% str_c(unique(out.df.filt$scaffold), ".fa")]
file.copy(from= file.path(dir_input, "scaff_fasta_files", file.list), 
          to= file.path(dir_input, "scaff_fasta_files_subset"))


# Whole alignment FASTA

## Read in alignment of the genome
scaff.data <- readData(path= file.path(dir_input_scaff, "scaff_fasta_files_subset"), format="FASTA", SNP.DATA = TRUE, big.data= FALSE)


#### SET POPULATIONS (i.e. ANI defined species 1 to 4)
individual.names <- sort(get.individuals(scaff.data)[[4]])


# populations <- list(c(individual.names[c(-1, -3, -12, -14, -17, -25, -26, -27, -13, -15, -16, -18, -19, -37, -53, -38, -39, -49, -48)]),
#                     c(individual.names[c(12, 14, 17, 25, 26, 27)]),
#                     c(individual.names[c(13, 15, 16, 18, 19, 37, 53, 38, 39, 49, 48)]),
#                     c(individual.names[c(1, 3)]))

populations <- list(c(individual.names[c(-10, -12, -15, -23, -24, -25, -11, -13, -14, -16, -17, -35, -51, -36, -37, -47, -46)]),
                    c(individual.names[c(10, 12, 15, 23, 24, 25)]),
                    c(individual.names[c(11, 13, 14, 16, 17, 35, 51, 36, 37, 47, 46)]))

scaff.data <- set.populations(scaff.data, populations, diploid= FALSE)

scaff.data.sum <- as.data.frame(get.sum.data(scaff.data)) %>%
                mutate(scaffold= row.names(.),
                       prop.valid.sites= round(n.valid.sites / n.sites, 4)) %>%
                select(scaffold, everything()) %>%
                as_tibble()


#### CALCULATE STATISTICS ######################################################
# show.slots(scaff.data) # slots are the different type of analyses that can be conducted

# Calculate Neutrality statistics
scaff.data <- neutrality.stats(scaff.data, detail= TRUE)
#get.neutrality(scaff.data)[[1]]
# scaff.data@Tajima.D
# scaff.data@n.segregating.sites

# Calculate F_ST and Diversity statistics
scaff.data <- F_ST.stats(scaff.data, mode= "nucleotide") # this also calculates diversity statistics

# scaff.data@nucleotide.F_ST
# scaff.data@nuc.F_ST.pairwise
# scaff.data@nuc.diversity.within # same as Pi
# scaff.data@nuc.diversity.between
# scaff.data@nuc.F_ST.vs.all

# Extract FST between populations
pairwise.FST <- t(scaff.data@nuc.F_ST.pairwise) %>%
  as_tibble() %>%
  #select(-`pop1/pop4`, -`pop2/pop4`, -`pop3/pop4`) %>%
  mutate(scaffold= scaff.data@region.names) %>%
  gather(key= pops, value= FST, `pop1/pop2`:`pop2/pop3`) %>%
  arrange(scaffold)

# Extract nucleotide diversity between and within populations
nucdiv.btw <- t(scaff.data@nuc.diversity.between / scaff.data.sum$n.sites) %>%
  as_tibble() %>%
  mutate(scaffold= scaff.data@region.names) %>%
  gather(key= pops, value= nuc.btw.div, `pop1/pop2`:`pop2/pop3`)%>%
  arrange(scaffold)

nucdiv.within <- (scaff.data@nuc.diversity.within / scaff.data.sum$n.sites) %>%
  as_tibble() %>%
  mutate(scaffold= scaff.data@region.names) %>%
  gather(key= pops, value= nuc.wtn.div, `pop 1`:`pop 3`)%>%
  arrange(scaffold)


## Extract Tajima D
tajima <- scaff.data@Tajima.D %>%
  as_tibble() %>%
  mutate(scaffold= scaff.data@region.names) %>%
  gather(key= pops, value= tajimaD, `pop 1`:`pop 3`) %>%
  arrange(scaffold)


#### Combine statistics into data frames
## Remove genes with valid sites < 50% of gene length
pop.stats.btw.df <- left_join(scaff.data.sum, pairwise.FST) %>% # BETWEEN STATISTICS
  left_join(., nucdiv.btw) %>%
  mutate(scaffold= str_replace(.$scaffold, ".fa", "")) #%>%
  #filter(prop.valid.sites > 0.5)
  
pop.stats.wtn.df <- left_join(scaff.data.sum, nucdiv.within) %>% # WITHIN STATISTICS
  left_join(., tajima) %>%
  mutate(scaffold= str_replace(.$scaffold, ".fa", "")) #%>%
  #filter(prop.valid.sites > 0.5) %>%
  
# rm(pairwise.FST, nucdiv.btw, nucdiv.within, tajima)



#### SLIDING WINDOW ANALYSIS ###############################################

scaff.data.slide <- sliding.window.transform(scaff.data, width= 1000, jump= 200, type= 2, whole.data= FALSE)

# Calculate F_ST and Diversity statistics
scaff.data.slide <- F_ST.stats(scaff.data.slide, mode= "nucleotide") # this also calculates diversity statistics
fst.df <- get.F_ST(scaff.data.slide, pairwise= TRUE)[[1]]
rowSums(fst.df)
show.slots(scaff.data.slide)
# Calculate Neutrality statistics
scaff.data.slide <- neutrality.stats(scaff.data.slide, detail= TRUE)
#get.neutrality(scaff.data.slide)[[1]]

scaff.data.slide@region.names
get.sum.data(scaff.data.slide[1])

fst.df.2 <- fst.df %>% 
  as_tibble() %>% 
  gather(key= pops, value= fst)

ggplot(data= fst.df.2) +
  geom_point(aes(x= pops, y= fst), size= 0.5) +
  labs(x= "Scaffold", y= expression("F"[st])) +
  #x.axis.format +
  #scale_y_continuous(limits= c(0.9, 1), expand= c(0.001, 0)) +
  #facet_wrap(~pops, nrow= 3, scales= "free_x", labeller= labeller(pops= pop.comp.facet.labels)) +
  theme_popgenome




  ##### PLOTTING PARAMETERS ######################################################
y.intercept <- geom_hline(yintercept = 0, color= "black", size= 0.25)
pop.comp.facet.labels <- as_labeller(c(`pop1/pop2` = "Species 1 & 2", `pop1/pop3` = "Species 1 & 3", `pop2/pop3` = "Species 2 & 3"))
pop.facet.labels <- as_labeller(c(`pop 1` = "Species 1", `pop 2` = "Species 2", `pop 3` = "Species 3"))

# x.axis.format.genes <- scale_x_discrete(breaks= pop.stats.btw.df$geneID[c(1, seq(300, length(pop.stats.btw.df$geneID), by= 300))],
#                                         labels= c(1, seq(300, length(pop.stats.btw.df$geneID), by= 300)/3),
#                                         expand= c(0.004, 0))
#max.break <- round(length(slide.data@region.names) * (7/7.3), 0) # Last x-axis break at 7 mega bp
#x.axis.format <-  scale_x_continuous(breaks= seq(0, max.break, length.out = 8), labels= seq(0, 7, by= 1), expand= c(0.01, 0))
x.axis.format <- scale_x_discrete(labels= NULL, expand= c(0.002, 0))



## ggplot themes
theme_popgenome <- theme(panel.grid = element_blank(),
                         plot.margin = unit(c(1, 1, 1, 1), "cm"),
                         text = element_text(size= 14),
                         plot.background = element_rect(fill = "transparent", color= "transparent"), # bg of the plot
                         panel.background = element_rect(fill= "transparent", color= "transparent"),
                         panel.border= element_rect(fill= NA, color= "black", linetype= "solid", size= 1),
                         panel.ontop = TRUE,
                         axis.text = element_text(colour="black"),
                         axis.title.x = element_text(vjust = -0.75),
                         axis.title.y = element_text(vjust = 1.5),
                         legend.background = element_rect(size=0.25, color="black", fill= "transparent"),
                         legend.key = element_blank(),
                         strip.background = element_rect(fill="transparent", color= "transparent"),
                         #axis.text.x = element_text(angle= 45, hjust= 1),
                         legend.position = "top")
##############################################################################

# Number of genes and total nucleotides inluded in analyses
length(unique(pop.stats.btw.df$geneID))
sum(pop.stats.btw.df$n.sites / 3)


#### Summary plots of number of valid sites used for analyses
ggplot(data= scaff.data.sum) +
  geom_histogram(aes(x= n.biallelic.sites), binwidth= 5, fill= "gray75", color= "black") +
  scale_y_log10(limits= c(0.95, 100), expand= c(0, 0)) +
  theme_popgenome

ggplot(data= scaff.data.sum) +
  geom_histogram(aes(x= prop.valid.sites), binwidth= 0.01, fill= "gray75", color= "black") +
  theme_popgenome

ggplot(data= pop.stats.btw.df) +
  geom_histogram(aes(x= n.valid.sites), binwidth= 100, fill= "gray75", color= "black") +
  theme_popgenome


## Fst


ggplot(data= pop.stats.btw.df) +
  geom_point(aes(x= scaffold, y= FST), size= 0.5) +
  labs(x= "Scaffold", y= expression("F"[st])) +
  x.axis.format +
  scale_y_continuous(limits= c(0.9, 1), expand= c(0.001, 0)) +
  facet_wrap(~pops, nrow= 3, scales= "free_x", labeller= labeller(pops= pop.comp.facet.labels)) +
  theme_popgenome
#ggsave(last_plot(), file= "Fst_genes.pdf", width= 8, height= 6, units= "in", path= file.path(dir_output_figures, "genes_alignment"))

ggplot(data= pop.stats.btw.df) +
  geom_boxplot(aes(x= pops, y= FST)) +
  labs(x= "", y= expression("F"[st])) +
  scale_y_continuous(limits= c(0, 1), expand= c(0.01, 0)) +
  scale_x_discrete(labels= c("Species 1 & 2", "Species 1 & 3", "Species 2 & 3")) +
  scale_color_manual(values= c("gray50", "black"), name= "Polymorphic gene") +
  theme_popgenome
#ggsave(last_plot(), file= "Fst_genes_boxplot.pdf", width= 8, height= 6, units= "in", path= file.path(dir_output_figures, "genes_alignment"))


## Nucleotide diversity BETWEEN species


ggplot(data= pop.stats.btw.df) +
  geom_point(aes(x= scaffold, y= nuc.btw.div), size= 0.5) +
  labs(x= "Scaffold", y= expression(paste("Inter-species nucleotide diversity (D"[x][y],")"))) +
  x.axis.format +
  #scale_y_continuous(limits= c(0, 1), expand= c(0.001, 0)) +
  facet_wrap(~pops, nrow= 3, scales= "free_x", labeller= labeller(pops= pop.comp.facet.labels)) +
  theme_popgenome
#ggsave(last_plot(), file= "Inter_specific_diversity_genes.pdf", width= 8, height= 6, units= "in", path= file.path(dir_output_figures, "genes_alignment"))

ggplot(data= pop.stats.btw.df) +
  geom_boxplot(aes(x= pops, y= nuc.btw.div)) +
  labs(x= "", y= expression(paste("Inter-species nucleotide diversity (D"[x][y],")"))) +
  scale_y_continuous(limits= c(0, 0.81), expand= c(0, 0.01)) +
  scale_x_discrete(labels= c("Species 1 & 2", "Species 1 & 3", "Species 2 & 3")) +
  scale_color_manual(values= c("gray50", "black"), name= "Polymorphic gene") +
  theme_popgenome
#ggsave(last_plot(), file= "Inter_specific_diversity_genes_boxplot.pdf", width= 8, height= 6, units= "in", path= file.path(dir_output_figures, "genes_alignment"))




## Nucleotide diversity WITHIN species


ggplot(data= pop.stats.wtn.df) +
  geom_point(aes(x= scaffold, y= nuc.wtn.div), size= 0.5) +
  labs(x= "Scaffold", y= expression(paste("Intra-species nucleotide diversity (D"[x][y],")"))) +
  x.axis.format +
  #scale_y_continuous(limits= c(0, 0.05), expand= c(0.03, 0)) +
  facet_wrap(~pops, nrow= 3, scales= "free_x", labeller= labeller(pops= pop.facet.labels)) +
  theme_popgenome
#ggsave(last_plot(), file= "Intra_specific_diversity_genes.pdf", width= 8, height= 6, units= "in", path= file.path(dir_output_figures, "genes_alignment"))

ggplot(data= pop.stats.wtn.df) +
  geom_boxplot(aes(x= pops, y= nuc.wtn.div)) +
  labs(x= "", y= expression("Intraspecific diversity ("~pi~")")) +
#  scale_y_continuous(limits= c(0, 0.05), expand= c(0.02, 0)) +
  scale_x_discrete(labels= c("Species 1", "Species 2", "Species 3")) +
  scale_color_manual(values= c("gray50", "black"), name= "Polymorphic gene") +
  theme_popgenome
#ggsave(last_plot(), file= "Intra_specific_diversity_genes_boxplot.pdf", width= 8, height= 6, units= "in", path= file.path(dir_output_figures, "genes_alignment"))


## Tajima's D

ggplot(data= pop.stats.wtn.df) +
  y.intercept +
  geom_hline(yintercept= c(-2, 2), linetype= "dashed", color= "black", size= 0.25) +
  geom_point(aes(x= scaffold, y= tajimaD), size= 0.5) +
  labs(x= "Scaffold", y= "Tajima's D") +
  x.axis.format +
  scale_y_continuous(limits= c(-3, 4), breaks= seq(-3, 4, by= 1), labels= c("", "-2", "", "0", "", "2", "", "4"), expand= c(0.01, 0)) +
  facet_wrap(~pops, nrow= 3, scales= "free_x", labeller= labeller(pops= pop.facet.labels)) +
  theme_popgenome
#ggsave(last_plot(), file= "TajimaD_genes.pdf", width= 8, height= 6, units= "in", path= file.path(dir_output_figures, "genes_alignment"))

ggplot(data= pop.stats.wtn.df) +
  y.intercept +
  geom_hline(yintercept= c(-2, 2), linetype= "dashed", color= "black", size= 0.25) +
  geom_boxplot(aes(x= pops, y= tajimaD)) +
  labs(x= "", y= "Tajima's D") +
  scale_y_continuous(limits= c(-3, 4), breaks= seq(-3, 4, by= 1), labels= c("", "-2", "", "0", "", "2", "", "4"), expand= c(0.01, 0)) +
  scale_x_discrete(labels= c("Species 1", "Species 2", "Species 3")) +
  scale_color_manual(values= c("gray50", "black"), name= "Polymorphic gene") +
  theme_popgenome
#ggsave(last_plot(), file= "TajimaD_genes_boxplot.pdf", width= 8, height= 6, units= "in", path= file.path(dir_output_figures, "genes_alignment"))





