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

## Run a command such as diversity.stats(mydata2)
## See output list contents get.diversity(mydata2)
## See contents of each list get.diversity(mydata2)[[1]]

#### Libraries #################################################################
library(tidyverse)
library(PopGenome)
library(ggplot2)
################################################################################

#### FILE PATHS ################################################################
dir_input <- file.path("/Users","kbg","Documents","UC_Berkeley","CyanoMeta_NSF","Metagenomics", "Data","GenomesData", "PopGenome", "Alignments_snippy", "Alignments_snippy_format")
dir_input_annotations <- file.path("/Users","kbg","Documents","UC_Berkeley","CyanoMeta_NSF","Metagenomics", "Data","GenomesData", "PopGenome", "Annotations")
dir_output_figures <- file.path("/Users","kbg","Documents","UC_Berkeley","CyanoMeta_NSF","Metagenomics", "Data","GenomesData", "PopGenome", "Output_figures")
dir_output_table <- file.path("/Users","kbg","Documents","UC_Berkeley","CyanoMeta_NSF","Metagenomics", "Data","GenomesData", "PopGenome", "Output_tables")
################################################################################



# Whole alignment FASTA
#mydata3 <- readData(file.path(dir_input, "snippy_LO_genome_fasta"), format= "FASTA", gffpath= file.path(dir_input, "GFF"), SNP.DATA= TRUE, big.data= TRUE, FAST= TRUE)
#mydata2 <- readData(path= "snippy_genome_fasta", format="FASTA", SNP.DATA = TRUE, big.data= TRUE, FAST= TRUE)
#mydata2 <- readData(path= file.path(dir_input, "gene_fasta_files"), format="FASTA", SNP.DATA = TRUE, big.data= TRUE)

## Read in alignment of the genes
mydata2 <- readData(path= file.path(dir_input, "gene_ggkbase_fasta_files"), format="FASTA", SNP.DATA = TRUE, big.data= TRUE)



## Start/end locations of each gene in GFF file

gff <- read.delim(file.path(dir_input_annotations, "Reference_ggkbase.gff"), header=F, comment.char="#", stringsAsFactors = FALSE) %>%
  as_tibble() %>%
  rename(seqid= V1, source= V2, type= V3, start= V4, end= V5, score= V6, strand= V7, phase= V8, attribute= V9) %>%
  mutate(geneID= str_c("gene_", str_pad(seq(1:nrow(.)), pad= "0", width= 4)),
         start= as.numeric(start),
         end= as.numeric(end))

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

## Scaffold lengths
scaff.df <- read.delim(file.path(dir_input, "Reference_scaffold_lengths.gff"), header=F, comment.char="#", stringsAsFactors = FALSE) %>% 
  as_tibble() %>% 
  select(V9, V4, V5) %>% 
  rename("scaffold"= V9, "start"= V4, "end"= V5) %>% 
  mutate(scaff_length= end - start + 1)


#### SET POPULATIONS (i.e. ANI defined species 1 to 4)
individual.names <- sort(get.individuals(mydata2)[[4]])


# populations <- list(c(individual.names[c(-1, -3, -12, -14, -17, -25, -26, -27, -13, -15, -16, -18, -19, -37, -53, -38, -39, -49, -48)]),
#                     c(individual.names[c(12, 14, 17, 25, 26, 27)]),
#                     c(individual.names[c(13, 15, 16, 18, 19, 37, 53, 38, 39, 49, 48)]),
#                     c(individual.names[c(1, 3)]))

populations <- list(c(individual.names[c(-10, -12, -15, -23, -24, -25, -11, -13, -14, -16, -17, -35, -51, -36, -37, -47, -46)]),
                    c(individual.names[c(10, 12, 15, 23, 24, 25)]),
                    c(individual.names[c(11, 13, 14, 16, 17, 35, 51, 36, 37, 47, 46)]))

mydata2 <- set.populations(mydata2, populations, diploid= FALSE)

mydata2.sum <- as.data.frame(get.sum.data(mydata2)) %>%
                mutate(geneID= row.names(.),
                       prop.valid.sites= round(n.valid.sites / n.sites, 4)) %>%
                select(geneID, everything()) %>%
                as_tibble()


#### CALCULATE STATISTICS ######################################################
# show.slots(mydata2) # slots are the different type of analyses that can be conducted

# Calculate Neutrality statistics
mydata2 <- neutrality.stats(mydata2, detail= TRUE)
#get.neutrality(mydata2)[[1]]
# mydata2@Tajima.D
# mydata2@n.segregating.sites

# Calculate F_ST and Diversity statistics
mydata2 <- F_ST.stats(mydata2, mode= "nucleotide") # this also calculates diversity statistics

# mydata2@nucleotide.F_ST
# mydata2@nuc.F_ST.pairwise
# mydata2@nuc.diversity.within # same as Pi
# mydata2@nuc.diversity.between
# mydata2@nuc.F_ST.vs.all

# Extract FST between populations
pairwise.FST <- t(mydata2@nuc.F_ST.pairwise) %>%
  as_tibble() %>%
  #select(-`pop1/pop4`, -`pop2/pop4`, -`pop3/pop4`) %>%
  mutate(geneID= mydata2@region.names) %>%
  gather(key= pops, value= FST, `pop1/pop2`:`pop2/pop3`) %>%
  arrange(geneID)

# Extract nucleotide diversity between and within populations
nucdiv.btw <- t(mydata2@nuc.diversity.between / mydata2.sum$n.sites) %>%
  as_tibble() %>%
  mutate(geneID= mydata2@region.names) %>%
  gather(key= pops, value= nuc.btw.div, `pop1/pop2`:`pop2/pop3`)%>%
  arrange(geneID)

nucdiv.within <- (mydata2@nuc.diversity.within / mydata2.sum$n.sites) %>%
  as_tibble() %>%
  mutate(geneID= mydata2@region.names) %>%
  gather(key= pops, value= nuc.wtn.div, `pop 1`:`pop 3`)%>%
  arrange(geneID)

## Extract Tajima D
tajima <- mydata2@Tajima.D %>%
  as_tibble() %>%
  mutate(geneID= mydata2@region.names) %>%
  gather(key= pops, value= tajimaD, `pop 1`:`pop 3`) %>%
  arrange(geneID)


#### Combine statistics into data frames
## Remove genes with valid sites < 50% of gene length
pop.stats.btw.df <- left_join(mydata2.sum, pairwise.FST) %>% # BETWEEN STATISTICS
  left_join(., nucdiv.btw) %>%
  mutate(geneID= str_replace(.$geneID, ".fa", "")) %>%
  filter(prop.valid.sites > 0.5) %>%
  left_join(., gff) %>%
  left_join(., anno.ref) %>%
  left_join(., subset(scaff.df, select= c(scaffold, scaff_length))) %>% 
  select(geneID, attribute, scaffold, scaff_length, start, end, n.sites:nuc.btw.div, uniref_anno, uniprot_anno, kegg_anno)

pop.stats.wtn.df <- left_join(mydata2.sum, nucdiv.within) %>% # WITHIN STATISTICS
  left_join(., tajima) %>%
  mutate(geneID= str_replace(.$geneID, ".fa", "")) %>%
  filter(prop.valid.sites > 0.5) %>%
  left_join(., gff) %>%
  left_join(., anno.ref) %>%
  left_join(., subset(scaff.df, select= c(scaffold, scaff_length))) %>% 
  select(geneID, attribute, scaffold, scaff_length, start, end, n.sites:tajimaD, uniref_anno, uniprot_anno, kegg_anno)

  #rm(pairwise.FST, nucdiv.btw, nucdiv.within, tajima)
  
## Calculate the location of each scaffold for plotting
scaff.breaks.btw <- pop.stats.btw.df %>% 
  group_by(scaffold) %>% 
  summarize(start_scaf= min(start),
            end_scaf= max(end)) %>% 
  arrange(start_scaf)



#### EXPORT OUTLIER TABLES #############################################################

# Fst
fst.outlier <- pop.stats.btw.df %>% 
  filter(FST < 0.5) %>% 
  mutate(outlier= "fst")
  write_tsv(fst.outlier, file.path(dir_output_table, "Fst_low.txt"))

# Intra-specific
intra.outlier <- pop.stats.wtn.df %>% 
  filter(nuc.wtn.div > 0.007) %>% 
  mutate(outlier= "intra")
  write_tsv(intra.outlier, file.path(dir_output_table, "Intra_specific_high.txt"))

# Inter-aspecific
inter.outlier <- pop.stats.btw.df %>% 
  filter(nuc.btw.div > 0.3) %>% 
  mutate(outlier= "inter")
  write_tsv(inter.outlier, file.path(dir_output_table, "Inter_specific_high.txt"))

# Tajima's D
tajima.outlier <- pop.stats.wtn.df %>% 
  filter(tajimaD > 2 | tajimaD < -2) %>% 
  mutate(outlier= "tajima")
  write_tsv(tajima.outlier, file.path(dir_output_table, "Tajima_2.txt"))

outlier.master <- full_join(fst.outlier, intra.outlier) %>% 
                     full_join(., inter.outlier) %>% 
                     full_join(., tajima.outlier)
write_tsv(outlier.master, file.path(dir_output_table, "PopGenome_outlier_genes.txt"))
    
  
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
ggplot(data= mydata2.sum) +
  geom_histogram(aes(x= n.biallelic.sites), binwidth= 5, fill= "gray75", color= "black") +
  scale_y_log10(limits= c(0.95, 100), expand= c(0, 0)) +
  theme_popgenome

ggplot(data= mydata2.sum) +
  geom_histogram(aes(x= prop.valid.sites), binwidth= 0.01, fill= "gray75", color= "black") +
  theme_popgenome

ggplot(data= pop.stats.btw.df) +
  geom_histogram(aes(x= n.valid.sites), binwidth= 100, fill= "gray75", color= "black") +
  theme_popgenome


## Fst
# ggplot(data= pop.stats.btw.df) +
#   geom_point(aes(x= geneID, y= FST), size= 0.5) +
#   labs(x= "Gene count", y= expression("F"[st])) +
#   x.axis.format +
#   scale_y_continuous(limits= c(0, 1), expand= c(0.04, 0)) +
#   scale_color_manual(values= c("gray50", "black"), name= "Polymorphic gene") +
#   facet_wrap(~pops, nrow= 3, labeller= labeller(pops= pop.comp.facet.labels)) +
#   theme_popgenome

ggplot(data= pop.stats.btw.df) +
  scaffold.breaks +
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
  scaffold.breaks +
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
  scale_y_continuous(limits= c(0, 0.81), breaks= seq(0, 0.8, by= 0.1), expand= c(0, 0.01)) +
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
  scaffold.breaks +
  geom_point(aes(x= start, y= nuc.wtn.div), size= 0.5) +
  labs(x= "Genome location (Mbp)", y= expression(paste("Intra-species nucleotide diversity (D"[x][y],")"))) +
  x.axis.format.bp +
  scale_y_continuous(limits= c(0, 0.05), expand= c(0.03, 0)) +
  scale_color_manual(values= c("gray50", "black"), name= "Polymorphic gene") +
  facet_wrap(~pops, nrow= 3, scales= "free_x", labeller= labeller(pops= pop.facet.labels)) +
  theme_popgenome
ggsave(last_plot(), file= "Intra_specific_diversity_genes.pdf", width= 8, height= 6, units= "in", path= file.path(dir_output_figures, "genes_alignment"))

ggplot(data= pop.stats.wtn.df) +
  geom_boxplot(aes(x= pops, y= nuc.wtn.div)) +
  labs(x= "", y= expression("Intraspecific diversity ("~pi~")")) +
  scale_y_continuous(limits= c(0, 0.05), expand= c(0.02, 0)) +
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











