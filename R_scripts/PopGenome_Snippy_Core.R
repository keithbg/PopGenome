## Analysis of Snippy Core alignment

## Snippy run on contigs
## biotite location: /data5/eelriver/CYA2/PH2017/snippy/snippy_out_contigs


#### POP GENOME INFORMATION ####################################################

## Regions are each of the nucleotide files included in the folder input into readData
## Populations are organisms within each data file defined to be in differen populations

## Run a command such as diversity.stats(mydata)
## See output list contents get.diversity(mydata)
## See contents of each list get.diversity(mydata)[[1]]

#### Libraries #################################################################
library(PopGenome)
################################################################################

#### FILE PATHS ################################################################
dir_input <- file.path("/Users","kbg","Documents","UC_Berkeley","CyanoMeta_NSF","Metagenomics", "Data","GenomesData", "PopGenome")
################################################################################


#### SNIPPY FILES

# Core SNPS VCF
mydata <- readData(path= "snippy_vcf", format="VCF", SNP.DATA = TRUE)
mydata <- readData(path= "snippy_LO_vcf", format="VCF", SNP.DATA = TRUE)
get.sum.data(mydata)
mydata@region.data

# Whole alignment FASTA
mydata2 <- readData(path= "snippy_genome_fasta", format="FASTA", SNP.DATA = TRUE, big.data= TRUE, FAST= TRUE)
get.sum.data(mydata2)
mydata2@region.data

#### SET POPULATIONS
individual.names <- sort(get.individuals(mydata2)[[1]])
populations <- list(c(individual.names[c(-1, -3, -12, -14, -17, -25, -26, -27, -13, -15, -16, -18, -19, -37, -53, -38, -39, -49, -48)]),
                    c(individual.names[c(12, 14, 17, 25, 26, 27)]),
                    c(individual.names[c(13, 15, 16, 18, 19, 37, 53, 38, 39, 49, 48)]),
                    c(individual.names[c(1, 3)]))
mydata <- set.populations(mydata, populations, diploid= FALSE)
mydata@populations
get.sum.data(mydata)


#### CALCULATE STATISTICS
show.slots(mydata) # slots are the different type of analyses that can be conducted

# Calculate Neutrality statistics
mydata <- neutrality.stats(mydata, detail= TRUE)
#get.neutrality(mydata)[[1]]
mydata@Tajima.D
mydata@n.segregating.sites

# Calculate F_ST and Diversity statistics
mydata <- F_ST.stats(mydata, mode= "nucleotide") # this also calculates diversity statistics


mydata@nucleotide.F_ST
mydata@nuc.F_ST.pairwise
mydata@nuc.diversity.within # same as Pi
mydata@nuc.diversity.between
mydata@nuc.F_ST.vs.all

#### SLIDING WINDOW TRANSFORM
#type=1, scans only biallelic positions; type= 2, scans the whole genome,
slide.data <- sliding.window.transform(mydata2, width= 5000, jump= 1000, type= 2, whole.data= TRUE)
#mydata2.slide@region.names



