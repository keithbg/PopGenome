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

#### Libraries #################################################################
library(tidyverse)
library(PopGenome)
library(ggplot2)
################################################################################

#### FILE PATHS ################################################################
#dir_core_genes <- "/Users/kbg/Documents/UC_Berkeley/CyanoMeta_NSF/Metagenomics/Data/GenomesData/Roary/Output_PH2015_sp123_Alignment_bp90_c50/core_genome_sequences"
dir_core_genes <- "/Users/kbg/Documents/UC_Berkeley/CyanoMeta_NSF/Metagenomics/Data/GenomesData/Roary/Output_PHall_bp95_c50/core_genome_sequences"
dir_core_genes_sp12 <- "/Users/kbg/Documents/UC_Berkeley/CyanoMeta_NSF/Metagenomics/Data/GenomesData/Roary/Output_PHall_bp95_c50/core_genome_sequences_sp1-2"

dir_input_annotations <- file.path("/Users","kbg","Documents","UC_Berkeley","CyanoMeta_NSF","Metagenomics", "Data","GenomesData", "PopGenome", "Annotations")
dir_output_figures <- file.path("/Users","kbg","Documents","UC_Berkeley","CyanoMeta_NSF","Metagenomics", "Data","GenomesData", "PopGenome", "Output_figures")
dir_output_table <- file.path("/Users","kbg","Documents","UC_Berkeley","CyanoMeta_NSF","Metagenomics", "Data","GenomesData", "PopGenome", "Output_tables")
################################################################################

## Run Bash script using terminal from "core_genes_sequences" folder to format fasta headers for PopGenome
"bash /Users/kbg/Documents/UC_Berkeley/CyanoMeta_NSF/Metagenomics/Data/GenomesData/Roary/Scripts/format_core_fasta_header.sh"


##### PLOTTING PARAMETERS ######################################################
y.intercept <- geom_hline(yintercept = 0, color= "black", size= 0.25)
scaffold.breaks <-   geom_vline(data= scaff.breaks.btw, aes(xintercept= end_scaf), size= 0.2, alpha= 0.5, color= "gray50")
x.axis.format.bp <- scale_x_continuous(breaks= c(1, seq(500000, 7000000, by= 500000)),
                                  labels= c("0", "", "1", "", "2", "", "3", "", "4", "", "5", "", "6", "", "7"),
                                  expand= c(0.005, 0))
pop.comp.facet.labels <- as_labeller(c(`pop1/pop2` = "Species 1 & 2", `pop1/pop3` = "Species 1 & 3", `pop2/pop3` = "Species 2 & 3"))
pop.facet.labels <- as_labeller(c(`pop 1` = "Species 1", `pop 2` = "Species 2", `pop 3` = "Species 3"))

# x.axis.format.genes <- scale_x_discrete(breaks= pop.stats.btw.df$geneID[c(1, seq(300, length(pop.stats.btw.df$geneID), by= 300))],
#                                         labels= c(1, seq(300, length(pop.stats.btw.df$geneID), by= 300)/3),
#                                         expand= c(0.004, 0))
#max.break <- round(length(slide.data@region.names) * (7/7.3), 0) # Last x-axis break at 7 mega bp
#x.axis.format <-  scale_x_continuous(breaks= seq(0, max.break, length.out = 8), labels= seq(0, 7, by= 1), expand= c(0.01, 0))
#x.axis.format <- scale_x_discrete(breaks= NULL, labels= NULL, expand= c(0.002, 0))



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
ggplot(data= mydata.sum) +
  geom_histogram(aes(x= n.biallelic.sites), binwidth= 5, fill= "gray75", color= "black") +
  scale_y_log10(limits= c(0.95, 100), expand= c(0, 0)) +
  theme_popgenome

ggplot(data= mydata.sum) +
  geom_histogram(aes(x= prop.valid.sites), binwidth= 0.01, fill= "gray75", color= "black") +
  theme_popgenome

ggplot(data= pop.stats.btw.df) +
  geom_histogram(aes(x= n.valid.sites), binwidth= 100, fill= "gray75", color= "black") +
  theme_popgenome


## Fst
ggplot(data= pop.stats.btw.df) +
  geom_point(aes(x= geneID, y= FST), size= 0.5) +
  labs(x= "Gene count", y= expression("F"[st])) +
  #x.axis.format +
  scale_x_discrete(labels= NULL) +
  scale_y_continuous(limits= c(0, 1), expand= c(0.04, 0)) +
  #scale_color_manual(values= c("gray50", "black"), name= "Polymorphic gene") +
  facet_grid(pops~.) +
  #facet_wrap(~pops, nrow= 3, labeller= labeller(pops= pop.comp.facet.labels)) +
  theme_popgenome

ggplot(data= pop.stats.btw.df) +
  #scaffold.breaks +
  geom_point(aes(x= start, y= FST), size= 0.5) +
  labs(x= "Genome location (Mbp)", y= expression("F"[st])) +
  x.axis.format.bp +
  scale_y_continuous(limits= c(0, 1), expand= c(0.001, 0)) +
  scale_color_manual(values= c("gray50", "black"), name= "Polymorphic gene") +
  facet_wrap(~pops, nrow= 3, scales= "free_x", labeller= labeller(pops= pop.comp.facet.labels)) +
  theme_popgenome
ggsave(last_plot(), file= "Fst_genes.pdf", width= 8, height= 6, units= "in", path= file.path(dir_output_figures, "genes_alignment"))

ggplot(data= pop.stats.btw.df) +
  geom_boxplot(aes(x= pops, y= FST)) +
  labs(x= "", y= expression("F"[st])) +
  scale_y_continuous(limits= c(0, 1), expand= c(0.01, 0)) +
  scale_x_discrete(labels= c("Species 1 & 2", "Species 1 & 3", "Species 2 & 3")) +
  scale_color_manual(values= c("gray50", "black"), name= "Polymorphic gene") +
  theme_popgenome
ggsave(last_plot(), file= "Fst_genes_boxplot.pdf", width= 8, height= 6, units= "in", path= file.path(dir_output_figures, "genes_alignment"))


## Nucleotide diversity BETWEEN species
# ggplot(data= pop.stats.btw.df) +
#   geom_point(aes(x= geneID, y= nuc.btw.div), size= 0.5) +
#   labs(x= "Gene count", y= expression(paste("Inter-species nucleotide diversity (D"[x][y],")"))) +
#   x.axis.format +
#   scale_y_continuous(limits= c(0, 0.81), expand= c(0.03, 0)) +
#   scale_color_manual(values= c("gray50", "black"), name= "Polymorphic gene") +
#   facet_wrap(~pops, nrow= 3, labeller= labeller(pops= pop.comp.facet.labels)) +
#   theme_popgenome

ggplot(data= pop.stats.btw.df) +
  #scaffold.breaks +
  geom_point(aes(x= start, y= nuc.btw.div), size= 0.5) +
  labs(x= "Genome location (Mbp)", y= expression(paste("Inter-species nucleotide diversity (D"[x][y],")"))) +
  x.axis.format.bp +
  scale_y_continuous(limits= c(0, 1), expand= c(0.001, 0)) +
  scale_color_manual(values= c("gray50", "black"), name= "Polymorphic gene") +
  facet_wrap(~pops, nrow= 3, scales= "free_x", labeller= labeller(pops= pop.comp.facet.labels)) +
  theme_popgenome
ggsave(last_plot(), file= "Inter_specific_diversity_genes.pdf", width= 8, height= 6, units= "in", path= file.path(dir_output_figures, "genes_alignment"))

ggplot(data= pop.stats.btw.df) +
  geom_boxplot(aes(x= pops, y= nuc.btw.div)) +
  labs(x= "", y= expression(paste("Inter-species nucleotide diversity (D"[x][y],")"))) +
  scale_y_continuous(limits= c(0, 0.9), breaks= seq(0, 0.9, by= 0.1), expand= c(0, 0.01)) +
  scale_x_discrete(labels= c("Species 1 & 2", "Species 1 & 3", "Species 2 & 3")) +
  scale_color_manual(values= c("gray50", "black"), name= "Polymorphic gene") +
  theme_popgenome
ggsave(last_plot(), file= "Inter_specific_diversity_genes_boxplot.pdf", width= 8, height= 6, units= "in", path= file.path(dir_output_figures, "genes_alignment"))




## Nucleotide diversity WITHIN species
# ggplot(data= pop.stats.wtn.df) +
#   geom_point(aes(x= geneID, y= nuc.wtn.div), size= 0.5) +
#   labs(x= "Gene count", y= expression(paste("Intra-species nucleotide diversity (",pi,")"))) +
#   x.axis.format +
#   scale_y_continuous(limits= c(0, 0.05), expand= c(0.03, 0)) +
#   scale_color_manual(values= c("gray50", "black"), name= "Polymorphic gene") +
#   facet_wrap(~pops, nrow= 3, labeller= labeller(pops= pop.facet.labels)) +
#   theme_popgenome

ggplot(data= pop.stats.wtn.df) +
  #scaffold.breaks +
  geom_point(aes(x= start, y= nuc.wtn.div), size= 0.5) +
  labs(x= "Genome location (Mbp)", y= expression(paste("Intra-species nucleotide diversity (D"[x][y],")"))) +
  x.axis.format.bp +
  scale_y_continuous(limits= c(0, 0.055), breaks= seq(0, 0.05, by= 0.01), labels= c("0.00", "", "0.02", "", "0.04", ""), expand= c(0.03, 0)) +
  scale_color_manual(values= c("gray50", "black"), name= "Polymorphic gene") +
  facet_wrap(~pops, nrow= 3, scales= "free_x", labeller= labeller(pops= pop.facet.labels)) +
  theme_popgenome
ggsave(last_plot(), file= "Intra_specific_diversity_genes.pdf", width= 8, height= 6, units= "in", path= file.path(dir_output_figures, "genes_alignment"))

ggplot(data= pop.stats.wtn.df) +
  geom_boxplot(aes(x= pops, y= nuc.wtn.div)) +
  labs(x= "", y= expression("Intraspecific diversity ("~pi~")")) +
  scale_y_continuous(limits= c(0, 0.055), expand= c(0.02, 0)) +
  scale_x_discrete(labels= c("Species 1", "Species 2", "Species 3")) +
  scale_color_manual(values= c("gray50", "black"), name= "Polymorphic gene") +
  theme_popgenome
ggsave(last_plot(), file= "Intra_specific_diversity_genes_boxplot.pdf", width= 8, height= 6, units= "in", path= file.path(dir_output_figures, "genes_alignment"))


## Tajima's D
# ggplot(data= pop.stats.wtn.df) +
#   y.intercept +
#   geom_hline(yintercept= c(-2, 2), linetype= "dashed", color= "black", size= 0.25) +
#   geom_point(aes(x= geneID, y= tajimaD), size= 0.5) +
#   labs(x= "Gene count", y= "Tajima's D") +
#   x.axis.format +
#   scale_y_continuous(limits= c(-3, 4), breaks= seq(-3, 4, by= 1), expand= c(0.01, 0)) +
#   scale_color_manual(values= c("gray50", "black"), name= "Polymorphic gene") +
#   facet_wrap(~pops, nrow= 3, labeller= labeller(pops= pop.facet.labels)) +
#   theme_popgenome


ggplot(data= pop.stats.wtn.df) +
  y.intercept +
  geom_hline(yintercept= c(-2, 2), linetype= "dashed", color= "black", size= 0.25) +
  geom_point(aes(x= start, y= tajimaD), size= 0.5) +
  labs(x= "Genome location (Mbp)", y= "Tajima's D") +
  x.axis.format.bp +
  scale_y_continuous(limits= c(-3, 4), breaks= seq(-3, 4, by= 1), labels= c("", "-2", "", "0", "", "2", "", "4"), expand= c(0.01, 0)) +
  scale_color_manual(values= c("gray50", "black"), name= "Polymorphic gene") +
  facet_wrap(~pops, nrow= 3, scales= "free_x", labeller= labeller(pops= pop.facet.labels)) +
  theme_popgenome
ggsave(last_plot(), file= "TajimaD_genes.pdf", width= 8, height= 6, units= "in", path= file.path(dir_output_figures, "genes_alignment"))

ggplot(data= pop.stats.wtn.df) +
  y.intercept +
  geom_hline(yintercept= c(-2, 2), linetype= "dashed", color= "black", size= 0.25) +
  geom_boxplot(aes(x= pops, y= tajimaD)) +
  labs(x= "", y= "Tajima's D") +
  scale_y_continuous(limits= c(-3, 4), breaks= seq(-3, 4, by= 1), labels= c("", "-2", "", "0", "", "2", "", "4"), expand= c(0.01, 0)) +
  scale_x_discrete(labels= c("Species 1", "Species 2", "Species 3")) +
  scale_color_manual(values= c("gray50", "black"), name= "Polymorphic gene") +
  theme_popgenome
ggsave(last_plot(), file= "TajimaD_genes_boxplot.pdf", width= 8, height= 6, units= "in", path= file.path(dir_output_figures, "genes_alignment"))











