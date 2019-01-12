## Analysis of Snippy Core alignment

## Snippy run on contigs
## biotite location: /data5/eelriver/CYA2/PH2017/snippy/snippy_out_contigs


#### POP GENOME INFORMATION ####################################################

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
dir_input <- file.path("/Users","kbg","Documents","UC_Berkeley","CyanoMeta_NSF","Metagenomics", "Data","GenomesData", "PopGenome")
################################################################################



# Whole alignment FASTA
#mydata3 <- readData(file.path(dir_input, "snippy_LO_genome_fasta"), format= "FASTA", gffpath= file.path(dir_input, "GFF"), SNP.DATA= TRUE, big.data= TRUE, FAST= TRUE)

mydata2 <- readData(path= file.path(dir_input, "snippy_LO_genome_fasta"), format="FASTA", SNP.DATA = TRUE, big.data= TRUE, FAST= TRUE)
#mydata2 <- readData(path= "snippy_genome_fasta", format="FASTA", SNP.DATA = TRUE, big.data= TRUE, FAST= TRUE)
#mydata3 <- readData(path= "snippy_genome_fasta", format="FASTA", SNP.DATA = TRUE, big.data= TRUE) # ~8 mins to load data withou FAST= TRUE
get.sum.data(mydata2)
mydata2@region.data


#get_gff_info(mydata2, file.path(dir_input, "GFF"), chr= "Reference", position= 1)[1]


#### SET POPULATIONS
individual.names <- sort(get.individuals(mydata2)[[1]])
# populations <- list(c(individual.names[c(-1, -3, -12, -14, -17, -25, -26, -27, -13, -15, -16, -18, -19, -37, -53, -38, -39, -49, -48)]),
#                     c(individual.names[c(12, 14, 17, 25, 26, 27)]),
#                     c(individual.names[c(13, 15, 16, 18, 19, 37, 53, 38, 39, 49, 48)]),
#                     c(individual.names[c(1, 3)]))

populations <- list(c(individual.names[c(-1, -11, -13, -16, -24, -25, -26, -12, -14, -15, -17, -18, -36, -52, -37, -38, -48, -47)]),
                    c(individual.names[c(11, 13, 16, 24, 25, 26)]),
                    c(individual.names[c(12, 14, 15, 17, 18, 36, 52, 37, 38, 48, 47)]),
                    c(individual.names[c(1)]))

mydata2 <- set.populations(mydata2, populations, diploid= FALSE)
mydata2@populations
get.sum.data(mydata2)
mydata2@region.data@polyallelic.sites

#### CALCULATE STATISTICS
show.slots(mydata2) # slots are the different type of analyses that can be conducted
get.sum.data(mydata2)
# Calculate Neutrality statistics
mydata2 <- neutrality.stats(mydata2, detail= TRUE)
#get.neutrality(mydata2)[[1]]
mydata2@Tajima.D
mydata2@n.segregating.sites

# Calculate F_ST and Diversity statistics
mydata2 <- F_ST.stats(mydata2, mode= "nucleotide") # this also calculates diversity statistics


mydata2@nucleotide.F_ST
mydata2@nuc.F_ST.pairwise
mydata2@nuc.diversity.within # same as Pi
mydata2@nuc.diversity.between
mydata2@nuc.F_ST.vs.all


#### SLIDING WINDOW TRANSFORM
# type=1, scans only biallelic positions; type= 2, scans the whole genome
window.width <- 10000
window.jump <- window.width * 0.5
slide.data <- sliding.window.transform(mydata2, width= window.width, jump= window.jump, type= 2, whole.data= TRUE)

show.slots(slide.data)
head(slide.data@n.sites)

## Calculate F_ST, Diversity, and Neutrality on sliding window
slide.data <- F_ST.stats(slide.data, mode= "nucleotide")
slide.data <- neutrality.stats(slide.data, detail= TRUE)
slide.data@Tajima.D

# Extract FST between populations
pairwise.FST <- t(slide.data@nuc.F_ST.pairwise) %>%
                  as_tibble() %>%
                  select(-`pop1/pop4`, -`pop2/pop4`, -`pop3/pop4`) %>%
                  mutate(ID= seq(1, nrow(.))) %>%
                  gather(key= pops, value= FST, `pop1/pop2`:`pop2/pop3`) %>%
                  arrange(ID)

FST_0 <- pairwise.FST %>%
  filter(FST == 0)
slide.data@region.names[77]


# Extract nucleotide diversity between populations
nucdiv.btw <- t(slide.data@nuc.diversity.between / window.width) %>%
  as_tibble() %>%
  select(-`pop1/pop4`, -`pop2/pop4`, -`pop3/pop4`) %>%
  mutate(ID= seq(1, nrow(.))) %>%
  gather(key= pops, value= nuc.btw.div, `pop1/pop2`:`pop2/pop3`)

nucdiv.within <- (slide.data@nuc.diversity.within / window.width) %>%
  as_tibble() %>%
  select( -`pop 4`) %>%
  mutate(ID= seq(1, nrow(.))) %>%
  gather(key= pops, value= nuc.wtn.div, `pop 1`:`pop 3`)

## Extract Tajima D
tajima <- slide.data@Tajima.D %>%
  as_tibble() %>%
  select( -`pop 4`) %>%
  mutate(ID= seq(1, nrow(.))) %>%
  gather(key= pops, value= tajimaD, `pop 1`:`pop 3`) %>%
  arrange(ID)


##### PLOTTING PARAMETERS ######################################################
max.break <- round(length(slide.data@region.names) * (7/7.3), 0) # Last x-axis break at 7 mega bp
x.axis.format <-  scale_x_continuous(breaks= seq(0, max.break, length.out = 8), labels= seq(0, 7, by= 1), expand= c(0.01, 0))


pop.comp.facet.labels <- as_labeller(c(`pop1/pop2` = "Species 1 & 2", `pop1/pop3` = "Species 1 & 3", `pop2/pop3` = "Species 2 & 3"))
pop.facet.labels <- as_labeller(c(`pop 1` = "Species 1", `pop 2` = "Species 2", `pop 3` = "Species 3"))

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






ggplot(data= pairwise.FST) +
  geom_point(aes(x= ID, y= FST)) +
  labs(x= "Genome position (Mbp)", y= "FST") +
  x.axis.format +
  facet_wrap(~pops, nrow= 3, labeller= labeller(pops= pop.comp.facet.labels)) +
  theme_popgenome

ggplot(data= nucdiv.btw) +
  geom_point(aes(x= ID, y= nuc.btw.div)) +
  labs(x= "Genome position (Mbp)", y= "Nucleotide diversity between species") +
  x.axis.format +
  facet_wrap(~pops, nrow= 3, labeller= labeller(pops= pop.comp.facet.labels)) +
  theme_popgenome

ggplot(data= nucdiv.within) +
  geom_point(aes(x= ID, y= nuc.wtn.div)) +
  labs(x= "Genome position (Mbp)", y= "Nucleotide diversity within species") +
  x.axis.format +
  facet_wrap(~pops, nrow= 3, labeller= labeller(pops= pop.facet.labels)) +
  theme_popgenome

ggplot(data= tajima) +
  geom_hline(yintercept = 0, color= "gray", size= 0.5) +
  geom_point(aes(x= ID, y= tajimaD)) +
  labs(x= "Genome position (Mbp)", y= "Nucleotide diversity within species") +
  x.axis.format +
  ylim(-3, 3) +
  facet_wrap(~pops, nrow= 3, labeller= labeller(pops= pop.facet.labels)) +
  theme_popgenome


