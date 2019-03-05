## PopGenome statistics on genes extracted from Roary output
## Plotting the data that was formatted in PopGenome_Roary_Analysis.R,


#### LIBRARIES 
library(tidyverse)
library(ggplot2)
library(ggExtra)


#### FILE PATHS 
dir_roary <- "../Roary/Output_PHall_bp90_c90"
dir_core_genes <- "../Roary/Output_PHall_bp90_c90/core_genome_sequences"
dir_out_table <- "Output_tables"
dir_out_figs <- "Output_figures"


#### READ DATA
pg.btw.filt <- read_tsv(file.path(dir_out_table, "PopGenome_btw_filt.tsv"))
pg.wtn.filt <- read_tsv(file.path(dir_out_table, "PopGenome_wtn_filt.tsv"))

##### PLOTTING PARAMETERS 
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

#### MAKE PLOTS ##############################################################
# Number of genes and total nucleotides inluded in analyses
length(unique(pg.btw$geneID))
sum(pg.btw$n.sites / 3)


#### Summary plots of number of valid sites used for analyses
ggplot(data= mydata.sum) +
  geom_histogram(aes(x= n.biallelic.sites), binwidth= 5, fill= "gray75", color= "black") +
  scale_y_log10(limits= c(0.95, 100), expand= c(0, 0)) +
  theme_popgenome

ggplot(data= mydata.sum) +
  geom_histogram(aes(x= prop.valid.sites), binwidth= 0.01, fill= "gray75", color= "black") +
  theme_popgenome

ggplot(data= pg.btw) +
  geom_histogram(aes(x= n.valid.sites), binwidth= 100, fill= "gray75", color= "black") +
  theme_popgenome


## Fst
ggplot(data= pg.btw.filt) +
  geom_point(aes(x= geneID, y= FST), size= 0.5) +
  labs(x= "Gene count", y= expression("F"[st])) +
  #x.axis.format +
  scale_x_discrete(labels= NULL, breaks= NULL) +
  scale_y_continuous(limits= c(0, 1), expand= c(0.04, 0)) +
  facet_wrap(~pops, nrow= 3, labeller= labeller(pops= pop.comp.facet.labels)) +
  theme_popgenome
ggsave(last_plot(), file= "Fst_genes.pdf", width= 8, height= 6, units= "in", path= file.path(dir_output_figures, "genes_alignment"))

ggplot(data= pg.btw.filt) +
  geom_boxplot(aes(x= pops, y= FST)) +
  labs(x= "", y= expression("F"[st])) +
  scale_y_continuous(limits= c(0, 1), expand= c(0.01, 0)) +
  scale_x_discrete(labels= c("Species 1 & 2", "Species 1 & 3", "Species 2 & 3")) +
  theme_popgenome
#ggsave(last_plot(), file= "Fst_genes_boxplot.pdf", width= 8, height= 6, units= "in", path= file.path(dir_output_figures, "genes_alignment"))

ggplot(data= pg.btw.filt) +
  geom_histogram(aes(x= FST), binwidth = 0.01, color= "black", fill= "gray50") +
  #geom_density(aes(x= FST)) +
 # labs(y= "Count", x= expression("F"[st])) +
  #scale_x_continuous(limits= c(0, 1), expand= c(0.01, 0)) +
  #facet_grid(pops~., scales= "free_y", labeller= labeller(pops= pop.comp.facet.labels)) +
  scale_y_log10() +
  facet_grid(pops~., scales= "free_y") +
  theme_popgenome

?geom_density

## Nucleotide diversity BETWEEN species
# ggplot(data= pg.btw) +
#   geom_point(aes(x= geneID, y= nuc.btw.div), size= 0.5) +
#   labs(x= "Gene count", y= expression(paste("Inter-species nucleotide diversity (D"[x][y],")"))) +
#   x.axis.format +
#   scale_y_continuous(limits= c(0, 0.81), expand= c(0.03, 0)) +
#   scale_color_manual(values= c("gray50", "black"), name= "Polymorphic gene") +
#   facet_wrap(~pops, nrow= 3, labeller= labeller(pops= pop.comp.facet.labels)) +
#   theme_popgenome

ggplot(data= pg.btw) +
  #scaffold.breaks +
  geom_point(aes(x= geneID, y= nuc.btw.div), size= 0.5) +
  labs(x= "GeneID", y= expression(paste("Inter-species nucleotide diversity (D"[x][y],")"))) +
  #x.axis.format.bp +
  scale_x_discrete(labels= NULL) +
  scale_y_continuous(limits= c(0, 1), expand= c(0.001, 0)) +
  #scale_color_manual(values= c("gray50", "black"), name= "Polymorphic gene") +
  facet_wrap(~pops, nrow= 3, scales= "free_x", labeller= labeller(pops= pop.comp.facet.labels)) +
  theme_popgenome
#ggsave(last_plot(), file= "Inter_specific_diversity_genes.pdf", width= 8, height= 6, units= "in", path= file.path(dir_output_figures, "genes_alignment"))

ggplot(data= pg.btw) +
  geom_boxplot(aes(x= pops, y= nuc.btw.div)) +
  labs(x= "", y= expression(paste("Inter-species nucleotide diversity (D"[x][y],")"))) +
  scale_y_continuous(limits= c(0, 0.9), breaks= seq(0, 0.9, by= 0.1), expand= c(0, 0.01)) +
  scale_x_discrete(labels= c("Species 1 & 2", "Species 1 & 3", "Species 2 & 3")) +
  scale_color_manual(values= c("gray50", "black"), name= "Polymorphic gene") +
  theme_popgenome
ggsave(last_plot(), file= "Inter_specific_diversity_genes_boxplot.pdf", width= 8, height= 6, units= "in", path= file.path(dir_output_figures, "genes_alignment"))




## Nucleotide diversity WITHIN species
# ggplot(data= pg.wtn) +
#   geom_point(aes(x= geneID, y= nuc.wtn.div), size= 0.5) +
#   labs(x= "Gene count", y= expression(paste("Intra-species nucleotide diversity (",pi,")"))) +
#   x.axis.format +
#   scale_y_continuous(limits= c(0, 0.05), expand= c(0.03, 0)) +
#   scale_color_manual(values= c("gray50", "black"), name= "Polymorphic gene") +
#   facet_wrap(~pops, nrow= 3, labeller= labeller(pops= pop.facet.labels)) +
#   theme_popgenome

ggplot(data= filter(pg.wtn.filt, pops == "pop 1")) +
  #scaffold.breaks +
  geom_point(aes(x= reorder(geneID, -nuc.wtn.div), y= nuc.wtn.div), size= 0.5) +
  labs(x= "", y= expression("Intraspecific diversity ("~pi~")")) +
  #x.axis.format.bp +
  scale_x_discrete(labels= NULL) +
  #scale_y_continuous(limits= c(0, 0.055), breaks= seq(0, 0.05, by= 0.01), labels= c("0.00", "", "0.02", "", "0.04", ""), expand= c(0.03, 0)) +
  facet_wrap(~pops, nrow= 3, scales= "free_x", labeller= labeller(pops= pop.facet.labels)) +
  theme_popgenome


ggplot(data= pg.wtn.filt) +
  #scaffold.breaks +
  geom_point(aes(x= reorder(geneID, -nuc.wtn.div), y= nuc.wtn.div), size= 0.5) +
  labs(x= "", y= expression("Intraspecific diversity ("~pi~")")) +
  #x.axis.format.bp +
  scale_x_discrete(labels= NULL) +
  #scale_y_continuous(limits= c(0, 0.055), breaks= seq(0, 0.05, by= 0.01), labels= c("0.00", "", "0.02", "", "0.04", ""), expand= c(0.03, 0)) +
  facet_wrap(~pops, nrow= 3, scales= "free_x", labeller= labeller(pops= pop.facet.labels)) +
  theme_popgenome




#ggsave(last_plot(), file= "Intra_specific_diversity_genes.pdf", width= 8, height= 6, units= "in", path= file.path(dir_output_figures, "genes_alignment"))

# ggplot(data= pg.wtn) +
#   geom_boxplot(aes(x= pops, y= nuc.wtn.div)) +
#   labs(x= "", y= expression("Intraspecific diversity ("~pi~")")) +
#   scale_y_continuous(limits= c(0, 0.055), expand= c(0.02, 0)) +
#   scale_x_discrete(labels= c("Species 1", "Species 2", "Species 3")) +
#   theme_popgenome
#ggsave(last_plot(), file= "Intra_specific_diversity_genes_boxplot.pdf", width= 8, height= 6, units= "in", path= file.path(dir_output_figures, "genes_alignment"))


## Tajima's D
# ggplot(data= pg.wtn) +
#   y.intercept +
#   geom_hline(yintercept= c(-2, 2), linetype= "dashed", color= "black", size= 0.25) +
#   geom_point(aes(x= geneID, y= tajimaD), size= 0.5) +
#   labs(x= "Gene count", y= "Tajima's D") +
#   x.axis.format +
#   scale_y_continuous(limits= c(-3, 4), breaks= seq(-3, 4, by= 1), expand= c(0.01, 0)) +
#   scale_color_manual(values= c("gray50", "black"), name= "Polymorphic gene") +
#   facet_wrap(~pops, nrow= 3, labeller= labeller(pops= pop.facet.labels)) +
#   theme_popgenome


ggplot(data= pg.wtn) +
  y.intercept +
  geom_hline(yintercept= c(-2, 2), linetype= "dashed", color= "black", size= 0.25) +
  geom_point(aes(x= geneID, y= tajimaD), size= 0.5) +
  labs(x= "Genome location (Mbp)", y= "Tajima's D") +
  #x.axis.format.bp +
  scale_x_discrete(labels= NULL) +
  #scale_y_continuous(limits= c(-3, 4), breaks= seq(-3, 4, by= 1), labels= c("", "-2", "", "0", "", "2", "", "4"), expand= c(0.01, 0)) +
  #scale_color_manual(values= c("gray50", "black"), name= "Polymorphic gene") +
  facet_wrap(~pops, nrow= 3, scales= "free_x", labeller= labeller(pops= pop.facet.labels)) +
  theme_popgenome
ggsave(last_plot(), file= "TajimaD_genes.pdf", width= 8, height= 6, units= "in", path= file.path(dir_output_figures, "genes_alignment"))

ggplot(data= pg.wtn) +
  y.intercept +
  geom_hline(yintercept= c(-2, 2), linetype= "dashed", color= "black", size= 0.25) +
  geom_boxplot(aes(x= pops, y= tajimaD)) +
  labs(x= "", y= "Tajima's D") +
  scale_y_continuous(limits= c(-3, 4), breaks= seq(-3, 4, by= 1), labels= c("", "-2", "", "0", "", "2", "", "4"), expand= c(0.01, 0)) +
  scale_x_discrete(labels= c("Species 1", "Species 2", "Species 3")) +
  scale_color_manual(values= c("gray50", "black"), name= "Polymorphic gene") +
  theme_popgenome
ggsave(last_plot(), file= "TajimaD_genes_boxplot.pdf", width= 8, height= 6, units= "in", path= file.path(dir_output_figures, "genes_alignment"))











