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
pop.comp.facet.labels <- as_labeller(c(`pop1/pop2` = "Species 1 & 2", `pop1/pop3` = "Species 1 & 3", `pop2/pop3` = "Species 2 & 3"))
pop.facet.labels <- as_labeller(c(`pop 1` = "Species 1", `pop 2` = "Species 2", `pop 3` = "Species 3"))
x.axis.pops <- scale_x_discrete(labels= c("Species\n1 & 2", "Species\n1 & 3", "Species\n2 & 3"))


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
length(unique(pg.btw.filt$geneID))
sum(pg.btw.filt$n.sites / 3)


#### Summary plots of number of valid sites used for analyses
ggplot(data= pg.btw.filt) +
  geom_histogram(aes(x= n.biallelic.sites), binwidth= 5, fill= "gray75", color= "black") +
  scale_y_log10(limits= c(0.95, 100), expand= c(0, 0)) +
  theme_popgenome

ggplot(data= pg.btw.filt) +
  geom_histogram(aes(x= prop.valid.sites), binwidth= 0.01, fill= "gray75", color= "black") +
  theme_popgenome

ggplot(data= pg.btw.filt) +
  geom_histogram(aes(x= n.valid.sites), binwidth= 100, fill= "gray75", color= "black") +
  theme_popgenome


## Fst
ggplot(data= pg.btw.filt) +
  geom_point(aes(x= geneID, y= FST), size= 0.5) +
  labs(x= "Gene count", y= expression("F"[st])) +
  scale_x_discrete(labels= NULL, breaks= NULL) +
  scale_y_continuous(limits= c(0, 1), expand= c(0.04, 0)) +
  facet_wrap(~pops, nrow= 3, labeller= labeller(pops= pop.comp.facet.labels)) +
  theme_popgenome

ggplot(data= pg.btw.filt) +
  geom_boxplot(aes(x= pops, y= FST)) +
  labs(x= NULL, y= expression("F"[st])) +
  x.axis.pops +
  scale_y_continuous(limits= c(0, 1), expand= c(0.02, 0), breaks= seq(0, 1, by= 0.1), labels= c("0.0", "", "0.2", "", "0.4", "", "0.6", "", "0.8", "", "1.0")) +
  theme_popgenome
ggsave(last_plot(), file= "Fst_genes_boxplot.pdf", width= 8, height= 6, units= "in", path= filedir_out_figs, device= cairo_pdf)

# p1 <- ggplot(data= filter(pg.btw.filt, pops== "pop1/pop2")) +
#   geom_point(aes(x= geneID, y= FST), size= 0.5) +
#   labs(x= "Genes", y= expression("F"[st])) +
#   scale_y_continuous(limits= c(0, 1), expand= c(0.04, 0)) +
#   scale_x_discrete(labels= NULL, breaks= NULL) + 
#   theme_popgenome
# p1
#   ggMarginal(p1, margins= "y", type= "boxplot", size = 10)
# ?ggMarginal

  
## Nucleotide diversity BETWEEN species

ggplot(data= pg.btw.filt) +
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


ggplot(data= pg.btw.filt) +
  geom_boxplot(aes(x= pops, y= nuc.btw.div)) +
  labs(x= NULL, y= expression(paste("Inter-species nucleotide diversity (D"[x][y],")"))) +
  x.axis.pops +
  scale_y_continuous(limits= c(0, 1), expand= c(0.02, 0), breaks= seq(0, 1, by= 0.1), labels= c("0.0", "", "0.2", "", "0.4", "", "0.6", "", "0.8", "", "1.0")) +
  theme_popgenome
ggsave(last_plot(), file= "Inter_specific_diversity_genes_boxplot.pdf", width= 8, height= 6, units= "in", path= dir_out_figs, device= cairo_pdf)



## Nucleotide diversity WITHIN species

ggplot(data= pg.wtn.filt) +
  geom_boxplot(aes(x= pops, y= nuc.wtn.div)) +
  labs(x= NULL, y= expression("Intraspecific diversity ("~pi~")")) +
  x.axis.pops +
  #scale_y_continuous(limits= c(0, 1), expand= c(0.02, 0), breaks= seq(0, 1, by= 0.1), labels= c("0.0", "", "0.2", "", "0.4", "", "0.6", "", "0.8", "", "1.0")) +
  scale_y_continuous(limits= c(0, 0.15), expand= c(0.02, 0), breaks= seq(0, 0.15, by= 0.01), labels= c("0.0", rep("", 4), "0.05", rep("", 4), "0.10", rep("", 4), "0.15")) +
  theme_popgenome
ggsave(last_plot(), file= "Intra_specific_diversity_genes_boxplot.pdf", width= 8, height= 6, units= "in", path= dir_out_figs, device= cairo_pdf)


## Tajima's D
ggplot(data= pg.wtn.filt) +
  y.intercept +
  geom_hline(yintercept= c(-2, 2), linetype= "dashed", color= "black", size= 0.25) +
  geom_point(aes(x= geneID, y= tajimaD), size= 0.5) +
  labs(x= "Gene", y= "Tajima's D") +
  #x.axis.format.bp +
  scale_x_discrete(labels= NULL) +
  #scale_y_continuous(limits= c(-3, 4), breaks= seq(-3, 4, by= 1), labels= c("", "-2", "", "0", "", "2", "", "4"), expand= c(0.01, 0)) +
  #scale_color_manual(values= c("gray50", "black"), name= "Polymorphic gene") +
  facet_wrap(~pops, nrow= 3, scales= "free_x", labeller= labeller(pops= pop.facet.labels)) +
  theme_popgenome
#ggsave(last_plot(), file= "TajimaD_genes.pdf", width= 8, height= 6, units= "in", path= file.path(dir_output_figures, "genes_alignment"))

ggplot(data= pg.wtn.filt) +
  y.intercept +
  geom_hline(yintercept= c(-2, 2), linetype= "dashed", color= "black", size= 0.25) +
  geom_boxplot(aes(x= pops, y= tajimaD)) +
  labs(x= NULL, y= "Tajima's D") +
  scale_y_continuous(limits= c(-3.5, 4.5), breaks= seq(-3, 4, by= 1), labels= c("", "-2", "", "0", "", "2", "", "4"), expand= c(0.01, 0)) +
  x.axis.pops +
  theme_popgenome
ggsave(last_plot(), file= "TajimaD_genes_boxplot.pdf", width= 8, height= 6, units= "in", path= dir_out_figs, device= cairo_pdf)











