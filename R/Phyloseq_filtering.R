setwd("~/Workingdirectory")

library(phyloseq)
library(ggplot2)
library(ape)
library(hilldiv)
library(decontam)
library(microbiome)
library(vegan)
library(ape)
library(fantaxtic)
library(dplyr)
library(DESeq2)
library(RColorBrewer)
library(dplyr)
library(genefilter)
library(phytools)
library(nlme)
library(mvabund)
library(ggpubr)
library(fantaxtic)


#Importing and cleaning ASV table
ASVtable <- read.csv("ASVs.csv", row.names = 1, sep=";")

ASVtable <- otu_table(ASVtable, taxa_are_rows = TRUE)
physeq <- phyloseq(ASVtable)

#Importing and filtering taxonomy table
taxtable <- read.table("Curated_Tax.csv", row.names = 1, sep=",")
taxtable$ASV <- row.names(taxtable)

taxtable<-taxtable[rownames(taxtable) %in% rownames(ASVtable), ]
taxtable <- as.matrix(taxtable)
taxtable <- tax_table(taxtable)

#Importing metadata
SamData <- read.delim("SampleData.txt", row.names = 1)
str(SamData)
SamData$Sample.type=as.factor(SamData$Sample.type)
SamData$Location=as.factor(SamData$Location)

SamData <- sample_data(SamData)

#Importing phylogenetic tree
tree <- read.tree("asvseqs.aligned.raxml.bestTree")

#Merge into phyloseq object
physeq = phyloseq(physeq, taxtable, SamData, tree)

# ##################################################################################
# ######## Decontaminate ########
# ##################################################################################
#Sequencing depth
df <- as.data.frame(sample_data(physeq)) # Put sample_data into a ggplot-friendly data.frame
df$LibrarySize <- sample_sums(physeq)
df <- df[order(df$LibrarySize),]
df$Index <- seq(nrow(df))
ggplot(data=df, aes(x=Index, y=LibrarySize, color=Sample_or_Control)) + geom_point()


sample_data(physeq)$is.neg <- sample_data(physeq)$Sample_or_Control == "Control Sample"
contamdf.prev <- isContaminant(physeq, method="prevalence", neg="is.neg", batch = "Plate")
table(contamdf.prev$contaminant)

head(which(contamdf.prev$contaminant))

contamdf.prev01 <- isContaminant(physeq, method="prevalence", neg="is.neg", threshold=0.1, batch = "Plate")
table(contamdf.prev01$contaminant)

#Make phyloseq object of presence-absence in negative controls and true samples
ps.pa <- transform_sample_counts(physeq, function(abund) 1*(abund>0))
ps.pa.neg <- prune_samples(sample_data(ps.pa)$Sample_or_Control == "Control Sample", ps.pa)
ps.pa.pos <- prune_samples(sample_data(ps.pa)$Sample_or_Control == "True Sample", ps.pa)
#Make data.frame of prevalence in positive and negative samples
df.pa <- data.frame(pa.pos=taxa_sums(ps.pa.pos), pa.neg=taxa_sums(ps.pa.neg),
                    contaminant=contamdf.prev$contaminant)
ggplot(data=df.pa, aes(x=pa.neg, y=pa.pos, color=contaminant)) + geom_point() +
xlab("Prevalence (Negative Controls)") + ylab("Prevalence (True Samples)")

physeq <- prune_taxa(!contamdf.prev$contaminant, physeq)
physeq<- subset_samples(physeq, Sample_or_Control=="True Sample")
physeq <- subset_samples(physeq, Sample.type=="Mucus" | Sample.type=="Digesta")

#Reving ASVs mapping to either Vertebrata, Chloroplast or Mitochondria
physeq <- subset_taxa(physeq, Phylum != "Vertebrata")
physeq <- subset_taxa(physeq, Order != "Chloroplast")
physeq <- subset_taxa(physeq, Family != "Mitochondria")

##################################################################################
#Separating physeq into species phyloseq objects

### Aphanius ####
Physeq_AI<- subset_samples(physeq, Species=="Aphanius iberus")

#Removing samples from location not used in the current study
Physeq_AI<- subset_samples(Physeq_AI, Location!="AI0")
Physeq_AI<- subset_samples(Physeq_AI, Location!="AI5")

#Removing sample without pairerd digesta sample
Physeq_AI <- prune_samples(!(sample_names(Physeq_AI) %in% "AI2_4_M"), Physeq_AI) 

#Removing outliers 
Physeq_AI <- prune_samples(!(sample_names(Physeq_AI) %in% "AI2_10_M"), Physeq_AI)
Physeq_AI <- prune_samples(!(sample_names(Physeq_AI) %in% "AI2_10_D"), Physeq_AI)

Physeq_AI <- prune_taxa(taxa_sums(Physeq_AI) > 0, Physeq_AI)

#Rarecurves
# cols <- c("darkred", "forestgreen", "darkblue", "pink", "orange", "blue", "yellow")
# rarecurve(t(otu_table(Physeq_AI)), step=50, cex=0.6, col = cols, xlab = "Merged reads (0-200,000)", ylab= "ASVs", label = F)
# rarecurve(t(otu_table(Physeq_AI)), step=50, cex=0.6, xlim=c(0, 20000), col = cols, xlab = "Merged reads (0-20,000)", ylab= "ASVs", label = F)
# rarecurve(t(otu_table(Physeq_AI)), step=50, cex=0.6, xlim=c(0, 5000), col = cols, xlab = "Merged reads (0-5,000)", ylab= "ASVs", label = F)

# Filtering
ASVs_AI <- as.data.frame(otu_table(Physeq_AI))
tax_AI <- as.data.frame(tax_table(Physeq_AI))
sam_AI <- as.data.frame(sample_data(Physeq_AI))
ASVs_AI = t(ASVs_AI)
## Create relative abundance table
Y_AI.rel<-ASVs_AI/rep(rowSums(ASVs_AI),times=ncol(ASVs_AI))
## Condition of being >1% abundance in at least 1 sample.
cond1=apply(Y_AI.rel,2,max)>=0.01
## Condition of being present in at least 5 samples
cond2=colSums(ASVs_AI>0)>=5
## Apply both conditions
ASVs_AI=ASVs_AI[ ,cond1&cond2]
dim(ASVs_AI)
ASVs_AI = as.data.frame(t(ASVs_AI))
tax_AI<-tax_AI[rownames(tax_AI) %in% rownames(ASVs_AI), ]

tree <- read.tree("asvseqs.aligned.raxml.bestTree")
tree <- force.ultrametric(tree,method="extend")
tree_AI <- match_data(tax_AI,tree,output="tree")
match_data(tax_AI,tree_AI)

ASVs_AI <- otu_table(ASVs_AI, taxa_are_rows = T)
tax_AI <- tax_table(as.matrix(tax_AI))
sam_AI <- sample_data(sam_AI)

Physeq_AI <- merge_phyloseq(ASVs_AI, tax_AI, sam_AI, tree_AI)

#Relative abundance
Physeq_AI <- transform_sample_counts(Physeq_AI, function(Physeq_AI) Physeq_AI/sum(Physeq_AI) )

# Creating dataframes for statistical analysis (Diversity_composition_plots.R)
write.csv(otu_table(Physeq_AI), "count_AI.csv")
write.csv(tax_table(Physeq_AI), "tax_AI.csv")
write.csv(sample_data(Physeq_AI), "sample_AI.csv")
write.tree(phy_tree(Physeq_AI), "tree_AI")

### Gambusia ####
Physeq_GH<- subset_samples(physeq, Species=="Gambusia holbrooki")

#Removing samples from location not used in the current study
Physeq_GH <- subset_samples(Physeq_GH, Location!="GH1")
Physeq_GH <- subset_samples(Physeq_GH, Location!="GH8")

#Removing sample without pairerd digesta sample
Physeq_GH <- prune_samples(!(sample_names(Physeq_GH) %in% "GH4_9_M"), Physeq_GH)

#Removing samples suspected to be switched during laboratory work
Physeq_GH <- prune_samples(!(sample_names(Physeq_GH) %in% "GH2_13_M"), Physeq_GH) #Suspected to be digesta
Physeq_GH <- prune_samples(!(sample_names(Physeq_GH) %in% "GH2_13_D"), Physeq_GH) #Suspected to be mucus
Physeq_GH <- prune_taxa(taxa_sums(Physeq_GH) > 0, Physeq_GH)

# Rarecurves
# cols <- c("darkred", "forestgreen", "darkblue", "pink", "orange", "blue", "yellow")
# rarecurve(t(otu_table(Physeq_GH)), step=50, cex=0.6, col = cols, xlab = "Merged reads (0-350,000)", ylab= "ASVs", label = F)
# rarecurve(t(otu_table(Physeq_GH)), step=50, cex=0.6, xlim=c(0, 20000), col = cols, xlab = "Merged reads (0-20,000)", ylab= "ASVs", label = F)
# rarecurve(t(otu_table(Physeq_GH)), step=50, cex=0.6, xlim=c(0, 5000), col = cols, xlab = "Merged reads (0-5,000)", ylab= "ASVs", label = F)


# Filtering
ASVs_GH <- as.data.frame(otu_table(Physeq_GH))
tax_GH <- as.data.frame(tax_table(Physeq_GH))
sam_GH <- as.data.frame(sample_data(Physeq_GH))
ASVs_GH = t(ASVs_GH)
## Create relative abundance table
Y_GH.rel<-ASVs_GH/rep(rowSums(ASVs_GH),times=ncol(ASVs_GH))
## Condition of being >1% abundance in at least 1 sample.
cond1=apply(Y_GH.rel,2,max)>=0.01
## Condition of being present in at least 5 samples
cond2=colSums(ASVs_GH>0)>=5
## Apply both conditions
ASVs_GH=ASVs_GH[ ,cond1&cond2]
dim(ASVs_GH)
ASVs_GH = as.data.frame(t(ASVs_GH))
tax_GH<-tax_GH[rownames(tax_GH) %in% rownames(ASVs_GH), ]

tree <- read.tree("asvseqs.aligned.raxml.bestTree")
tree <- force.ultrametric(tree,method="extend")
tree_GH <- match_data(tax_GH,tree,output="tree")
match_data(tax_GH,tree_GH)

ASVs_GH <- otu_table(ASVs_GH, taxa_are_rows = T)
tax_GH <- tax_table(as.matrix(tax_GH))
sam_GH <- sample_data(sam_GH)

Physeq_GH <- merge_phyloseq(ASVs_GH, tax_GH, sam_GH, tree_GH)

#Remove samples with no ASVs left after filtering 
Physeq_GH <- prune_samples(!(sample_names(Physeq_GH) %in% "GH5_5_M"), Physeq_GH)
Physeq_GH <- prune_samples(!(sample_names(Physeq_GH) %in% "GH5_5_D"), Physeq_GH)
Physeq_GH <- prune_taxa(taxa_sums(Physeq_GH) > 0, Physeq_GH)

#Relative abundance
Physeq_GH <- transform_sample_counts(Physeq_GH, function(Physeq_GH) Physeq_GH/sum(Physeq_GH) )

# Creating dataframes for statistical analysis (Diversity_composition_plots.R)
write.csv(otu_table(Physeq_GH), "count_GH.csv")
write.csv(tax_table(Physeq_GH), "tax_GH.csv")
write.csv(sample_data(Physeq_GH), "sample_GH.csv")
write.tree(phy_tree(Physeq_GH), "tree_GH")


### Valencia ####
Physeq_VH <- subset_samples(physeq, Species=="Valencia hispanica")

#Removing samples from location not used in the current study
Physeq_VH <- subset_samples(Physeq_VH, Location!="VH2")
Physeq_VH <- subset_samples(Physeq_VH, Location!="VH0")
Physeq_VH <- subset_samples(Physeq_VH, Location!="VH3")

#Removing samples identified as outliers 
Physeq_VH <- prune_samples(!(sample_names(Physeq_VH) %in% "VH6_9_M"), Physeq_VH)
Physeq_VH <- prune_samples(!(sample_names(Physeq_VH) %in% "VH6_9_D"), Physeq_VH)
Physeq_VH <- prune_samples(!(sample_names(Physeq_VH) %in% "VH6_1_M"), Physeq_VH)
Physeq_VH <- prune_samples(!(sample_names(Physeq_VH) %in% "VH6_1_D"), Physeq_VH)
Physeq_VH <- prune_taxa(taxa_sums(Physeq_VH) > 0, Physeq_VH)

# Rarecurves
# cols <- c("darkred", "forestgreen", "darkblue", "pink", "orange", "blue", "yellow")
# rarecurve(t(otu_table(Physeq_VH)), step=50, cex=0.6, col = cols, xlab = "Merged reads (0-400,000)", ylab= "ASVs", label = F)
# rarecurve(t(otu_table(Physeq_VH)), step=50, cex=0.6, xlim=c(0, 20000), col = cols, xlab = "Merged reads (0-20,000)", ylab= "ASVs", label = F)
# rarecurve(t(otu_table(Physeq_VH)), step=50, cex=0.6, xlim=c(0, 5000), col = cols, xlab = "Merged reads (0-5,000)", ylab= "ASVs", label = F)

# Filtering
ASVs_VH <- as.data.frame(otu_table(Physeq_VH))
tax_VH <- as.data.frame(tax_table(Physeq_VH))
sam_VH <- as.data.frame(sample_data(Physeq_VH))
write.csv(sam_VH, "~/Desktop/sam_VH.csv")

ASVs_VH = t(ASVs_VH)
## Create relative abundance table
Y_VH.rel<-ASVs_VH/rep(rowSums(ASVs_VH),times=ncol(ASVs_VH))
## Condition of being >1% abundance in at least 1 sample.
cond1=apply(Y_VH.rel,2,max)>=0.01
## Condition of being present in at least 5 samples
cond2=colSums(ASVs_VH>0)>=5
## Apply both conditions
ASVs_VH=ASVs_VH[ ,cond1&cond2]
dim(ASVs_VH)
ASVs_VH = as.data.frame(t(ASVs_VH))
tax_VH<-tax_VH[rownames(tax_VH) %in% rownames(ASVs_VH), ]

tree <- read.tree("asvseqs.aligned.raxml.bestTree")
tree <- force.ultrametric(tree,method="extend")
tree_VH <- match_data(tax_VH,tree,output="tree")
match_data(tax_VH,tree_VH)

ASVs_VH <- otu_table(ASVs_VH, taxa_are_rows = T)
tax_VH <- tax_table(as.matrix(tax_VH))
sam_VH <- sample_data(sam_VH)

Physeq_VH <- merge_phyloseq(ASVs_VH, tax_VH, sam_VH, tree_VH)

#Relative abundace
Physeq_VH <- transform_sample_counts(Physeq_VH, function(Physeq_VH) Physeq_VH/sum(Physeq_VH) )


# Creating dataframes for statistical analysis (Diversity_composition_plots.R)
write.csv(otu_table(Physeq_VH), "count_VH.csv")
write.csv(tax_table(Physeq_VH), "tax_VH.csv")
write.csv(sample_data(Physeq_VH), "sample_VH.csv")
write.tree(phy_tree(Physeq_VH), "tree_VH")

