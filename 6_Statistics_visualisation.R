setwd("~/Dropbox/Lasse/Bioinformatics/Amplicon/Mucus_Content/Github_code")

library(hilldiv)
library(nlme)
library(phytools)
library(ggplot2)
library(phyloseq)
library(vegan)
library(ape)


### Aphanius iberus ###
ASVs_AI <- read.csv("counts_AI.csv", row.names = 1)
tax_AI <- read.csv("taxonomy_AI.csv", row.names = 1)
meta_AI <- read.csv("sampledata_AI.csv", row.names = 1)
tree <- read.tree("tree.tree")
tree <- force.ultrametric(tree,method="extend")
tree_AI <- match_data(ASVs_AI,tree,output="tree")
match_data(ASVs_AI,tree_AI)
meta_AI_hill =meta_AI[,c("SampleID", "Location")]
row.names(meta_AI_hill) = meta_AI_hill[,c("SampleID")]
# 
# Phylooseq object
ASVs_AI.physeq <- otu_table(ASVs_AI, taxa_are_rows = TRUE)
tax_AI.physeq <- as.matrix(tax_AI)
tax_AI.physeq <- tax_table(as.matrix(tax_AI))
meta_AI.physeq <- sample_data(meta_AI)
Physeq_AI = phyloseq(ASVs_AI.physeq, tax_AI.physeq, meta_AI.physeq, tree_AI)


# Diversity profiles
 AI_div_profil_AI <- div_profile(ASVs_AI, qvalues=seq(from = 0, to = 3, by = (0.2)), hierarchy = meta_AI_hill)
 div_profile_plot(AI_div_profil_AI)

# Diversity computation
q1phy_AI <- hill_div(ASVs_AI,qvalue=1 ,tree=tree_AI)

# Linear mixed-effects models
model.q1phy_AI <- lme(q1phy_AI~1+Sample.type, random = ~1|Location/Fish_ID, meta_AI)
summary(model.q1phy_AI)
plot(model.q1phy_AI)
qqnorm(model.q1phy_AI)

# Pairwise plot
paired_AI <- read.csv("paired_AI.csv", row.names = 1, sep = ";")

pdf("Paired_AI.pdf",width=12,height=12)


  ggplot(paired_AI, aes(x = Sample.type, y = q1phy)) +
  geom_boxplot(col = "black") +
  geom_line(aes(group = Fish_ID, col = is_increasing)) +
  scale_colour_manual(values = c("green", "red"))+ylim(0,22)+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.position = "none")+
     geom_point(aes(size=10,shape = factor(Sample.type)))

dev.off()

# Beta diversity
pairdis.q1phy <- pair_dis(ASVs_AI, qvalue=1, hierarchy=meta_AI, tree = tree_AI, level = "1", metric = "U")

L1_UqN.AI <- as.data.frame(pairdis.q1phy[["L1_UqN"]])

L1_UqN.AI <- as.dist(L1_UqN.AI)

GP.ord <- ordinate(Physeq_AI, "NMDS", distance = L1_UqN.AI)

ordi_AI  <- plot_ordination(Physeq_AI, GP.ord, type="samples", color="Location", shape = "Sample.type")+
  geom_line(aes(group = Fish_ID), alpha=0.4)+
  theme_bw()+geom_point(size = 3)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank())+
  ggtitle("Aphanius iberus")+
  theme(plot.title = element_text(hjust=0.5))

# Checking assumptions of PERMANOVA
test.homogenity_AI=betadisper(L1_UqN.AI, group = meta_AI$Location)
anova(test.homogenity_AI)

# Running PERMANOVA 1 
permanova_AI<-adonis(L1_UqN.AI ~ Sample.type*Location, data=meta_AI, permutations = 999, na.rm = TRUE, strata =meta_AI$Fish_ID)

# Preparing data for PERMANOVA 2
Physeq_AI_mucus <- subset_samples(Physeq_AI, Sample.type=="Mucus")
ASVs_AI_mucus <- as.data.frame(otu_table(Physeq_AI_mucus))
tree_AI_mucus <- phy_tree(Physeq_AI_mucus)
 
Physeq_AI_digesta <- subset_samples(Physeq_AI, Sample.type=="Digesta")
ASVs_AI_digesta <- as.data.frame(otu_table(Physeq_AI_digesta))
tree_AI_digesta <- phy_tree(Physeq_AI_digesta)

# Importing metadata only for either mucus or digesta
meta_AI_mucus <- read.csv("sample_AI_mucus.csv", row.names = 1, sep = ";")
meta_AI_digesta <- read.csv("sample_AI_digesta.csv", row.names = 1, sep = ";")
 
 
pairdis.q1phy_mucus <- pair_dis(ASVs_AI_mucus, qvalue=1, hierarchy=meta_AI_mucus, tree = tree_AI_mucus, level = "1", metric = "U")
pairdis.q1phy_digesta <- pair_dis(ASVs_AI_digesta, qvalue=1, hierarchy=meta_AI_digesta, tree = tree_AI_digesta, level = "1", metric = "U")
 
L1_UqN.AI_mucus <- as.data.frame(pairdis.q1phy_mucus[["L1_UqN"]])
L1_UqN.AI_digesta <- as.data.frame(pairdis.q1phy_digesta[["L1_UqN"]])
 
L1_UqN.AI_mucus <- as.dist(L1_UqN.AI_mucus)
L1_UqN.AI_digesta <- as.dist(L1_UqN.AI_digesta)

# Checking assumptions of PERMANOVA
test.homogenity_AI_mucus=betadisper(L1_UqN.AI_mucus, group = meta_AI_mucus$Location)
test.homogenity_AI_digesta=betadisper(L1_UqN.AI_digesta, group = meta_AI_digesta$Location)
anova(test.homogenity_AI_mucus)
anova(test.homogenity_AI_digesta)

# Running PERMANOVA 2
permanova_AI_mucus<-adonis(L1_UqN.AI_mucus ~ Location, data=meta_AI_mucus, permutations = 999, na.rm = TRUE)
permanova_AI_digesta<-adonis(L1_UqN.AI_digesta ~ Location, data=meta_AI_digesta, permutations = 999, na.rm = TRUE)
 
 

### Gambusia holbrooki ###
ASVs_GH <- read.csv("counts_GH.csv", row.names = 1)
tax_GH <- read.csv("taxonomy_GH.csv", row.names = 1)
meta_GH <- read.csv("sampledata_GH.csv", row.names = 1)
tree <- read.tree("tree.tree")
tree <- force.ultrametric(tree,method="extend")
tree_GH <- match_data(ASVs_GH,tree,output="tree")
match_data(ASVs_GH,tree_GH)
meta_GH_hill =meta_GH[,c("SampleID", "Location")]
row.names(meta_GH_hill) = meta_GH_hill[,c("SampleID")]
# 
# Phylooseq object
ASVs_GH.physeq <- otu_table(ASVs_GH, taxa_are_rows = TRUE)
tax_GH.physeq <- as.matrix(tax_GH)
tax_GH.physeq <- tax_table(as.matrix(tax_GH))
meta_GH.physeq <- sample_data(meta_GH)
Physeq_GH = phyloseq(ASVs_GH.physeq, tax_GH.physeq, meta_GH.physeq, tree_GH)


# Diversity profiles
AI_div_profil_GH <- div_profile(ASVs_GH, qvalues=seq(from = 0, to = 3, by = (0.2)), hierarchy = meta_GH_hill)
div_profile_plot(AI_div_profil_GH)

# Diversity computation
q1phy_GH <- hill_div(ASVs_GH,qvalue=1 ,tree=tree_GH)

# Linear mixed-effects models
model.q1phy_GH <- lme(q1phy_GH~1+Sample.type, random = ~1|Location/Fish_ID, meta_GH)
summary(model.q1phy_GH)
plot(model.q1phy_GH)
qqnorm(model.q1phy_GH)

# Pairwise plot
paired_GH <- read.csv("paired_GH.csv", row.names = 1, sep = ";")

pdf("Paired_GH.pdf",width=12,height=12)


ggplot(paired_GH, aes(x = Sample.type, y = q1phy)) +
  geom_boxplot(col = "black") +
  geom_line(aes(group = Fish_ID, col = is_increasing)) +
  scale_colour_manual(values = c("green", "red"))+ylim(0,22)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.position = "none")+
  geom_point(aes(size=10,shape = factor(Sample.type)))

dev.off()

# Beta diversity
pairdis.q1phy <- pair_dis(ASVs_GH, qvalue=1, hierarchy=meta_GH, tree = tree_GH, level = "1", metric = "U")

L1_UqN.AI <- as.data.frame(pairdis.q1phy[["L1_UqN"]])

L1_UqN.AI <- as.dist(L1_UqN.AI)

GP.ord <- ordinate(Physeq_GH, "NMDS", distance = L1_UqN.AI)

ordi_GH  <- plot_ordination(Physeq_GH, GP.ord, type="samples", color="Location", shape = "Sample.type")+
  geom_line(aes(group = Fish_ID), alpha=0.4)+
  theme_bw()+geom_point(size = 3)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank())+
  ggtitle("Aphanius iberus")+
  theme(plot.title = element_text(hjust=0.5))

# Checking assumptions of PERMANOVA
test.homogenity_GH=betadisper(L1_UqN.AI, group = meta_GH$Location)
anova(test.homogenity_GH)

# Running PERMANOVA 1 
permanova_GH<-adonis(L1_UqN.AI ~ Sample.type*Location, data=meta_GH, permutations = 999, na.rm = TRUE, strata =meta_GH$Fish_ID)

# Preparing data for PERMANOVA 2
Physeq_GH_mucus <- subset_samples(Physeq_GH, Sample.type=="Mucus")
ASVs_GH_mucus <- as.data.frame(otu_table(Physeq_GH_mucus))
tree_GH_mucus <- phy_tree(Physeq_GH_mucus)

Physeq_GH_digesta <- subset_samples(Physeq_GH, Sample.type=="Digesta")
ASVs_GH_digesta <- as.data.frame(otu_table(Physeq_GH_digesta))
tree_GH_digesta <- phy_tree(Physeq_GH_digesta)

# Importing metadata only for either mucus or digesta
meta_GH_mucus <- read.csv("sample_GH_mucus.csv", row.names = 1, sep = ";")
meta_GH_digesta <- read.csv("sample_GH_digesta.csv", row.names = 1, sep = ";")


pairdis.q1phy_mucus <- pair_dis(ASVs_GH_mucus, qvalue=1, hierarchy=meta_GH_mucus, tree = tree_GH_mucus, level = "1", metric = "U")
pairdis.q1phy_digesta <- pair_dis(ASVs_GH_digesta, qvalue=1, hierarchy=meta_GH_digesta, tree = tree_GH_digesta, level = "1", metric = "U")

L1_UqN.AI_mucus <- as.data.frame(pairdis.q1phy_mucus[["L1_UqN"]])
L1_UqN.AI_digesta <- as.data.frame(pairdis.q1phy_digesta[["L1_UqN"]])

L1_UqN.AI_mucus <- as.dist(L1_UqN.AI_mucus)
L1_UqN.AI_digesta <- as.dist(L1_UqN.AI_digesta)

# Checking assumptions of PERMANOVA
test.homogenity_GH_mucus=betadisper(L1_UqN.AI_mucus, group = meta_GH_mucus$Location)
test.homogenity_GH_digesta=betadisper(L1_UqN.AI_digesta, group = meta_GH_digesta$Location)
anova(test.homogenity_GH_mucus)
anova(test.homogenity_GH_digesta)

# Running PERMANOVA 2
permanova_GH_mucus<-adonis(L1_UqN.AI_mucus ~ Location, data=meta_GH_mucus, permutations = 999, na.rm = TRUE)
permanova_GH_digesta<-adonis(L1_UqN.AI_digesta ~ Location, data=meta_GH_digesta, permutations = 999, na.rm = TRUE)



# Valencia hispanica # 
ASVs_VH <- read.csv("counts_VH.csv", row.names = 1)
tax_VH <- read.csv("taxonomy_VH.csv", row.names = 1)
meta_VH <- read.csv("sampledata_VH.csv", row.names = 1)
tree <- read.tree("tree.tree")
tree <- force.ultrametric(tree,method="extend")
tree_VH <- match_data(ASVs_VH,tree,output="tree")
match_data(ASVs_VH,tree_VH)
meta_VH_hill =meta_VH[,c("SampleID", "Location")]
row.names(meta_VH_hill) = meta_VH_hill[,c("SampleID")]
# 
# Phylooseq object
ASVs_VH.physeq <- otu_table(ASVs_VH, taxa_are_rows = TRUE)
tax_VH.physeq <- as.matrix(tax_VH)
tax_VH.physeq <- tax_table(as.matrix(tax_VH))
meta_VH.physeq <- sample_data(meta_VH)
Physeq_VH = phyloseq(ASVs_VH.physeq, tax_VH.physeq, meta_VH.physeq, tree_VH)


# Diversity profiles
AI_div_profil_VH <- div_profile(ASVs_VH, qvalues=seq(from = 0, to = 3, by = (0.2)), hierarchy = meta_VH_hill)
div_profile_plot(AI_div_profil_VH)

# Diversity computation
q1phy_VH <- hill_div(ASVs_VH,qvalue=1 ,tree=tree_VH)

# Linear mixed-effects models
model.q1phy_VH <- lme(q1phy_VH~1+Sample.type, random = ~1|Location/Fish_ID, meta_VH)
summary(model.q1phy_VH)
plot(model.q1phy_VH)
qqnorm(model.q1phy_VH)

# Pairwise plot
paired_VH <- read.csv("paired_VH.csv", row.names = 1, sep = ";")

pdf("Paired_VH.pdf",width=12,height=12)

ggplot(paired_VH, aes(x = Sample.type, y = q1phy)) +
  geom_boxplot(col = "black") +
  geom_line(aes(group = Fish_ID, col = is_increasing)) +
  scale_colour_manual(values = c("green", "red"))+ylim(0,22)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.position = "none")+
  geom_point(aes(size=10,shape = factor(Sample.type)))

dev.off()

# Beta diversity
pairdis.q1phy <- pair_dis(ASVs_VH, qvalue=1, hierarchy=meta_VH, tree = tree_VH, level = "1", metric = "U")

L1_UqN.AI <- as.data.frame(pairdis.q1phy[["L1_UqN"]])

L1_UqN.AI <- as.dist(L1_UqN.AI)

GP.ord <- ordinate(Physeq_VH, "NMDS", distance = L1_UqN.AI)

ordi_VH  <- plot_ordination(Physeq_VH, GP.ord, type="samples", color="Location", shape = "Sample.type")+
  geom_line(aes(group = Fish_ID), alpha=0.4)+
  theme_bw()+geom_point(size = 3)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank())+
  ggtitle("Aphanius iberus")+
  theme(plot.title = element_text(hjust=0.5))

# Checking assumptions of PERMANOVA
test.homogenity_VH=betadisper(L1_UqN.AI, group = meta_VH$Location)
anova(test.homogenity_VH)

# Running PERMANOVA 1 
permanova_VH<-adonis(L1_UqN.AI ~ Sample.type*Location, data=meta_VH, permutations = 999, na.rm = TRUE, strata =meta_VH$Fish_ID)

# Preparing data for PERMANOVA 2
Physeq_VH_mucus <- subset_samples(Physeq_VH, Sample.type=="Mucus")
ASVs_VH_mucus <- as.data.frame(otu_table(Physeq_VH_mucus))
tree_VH_mucus <- phy_tree(Physeq_VH_mucus)

Physeq_VH_digesta <- subset_samples(Physeq_VH, Sample.type=="Digesta")
ASVs_VH_digesta <- as.data.frame(otu_table(Physeq_VH_digesta))
tree_VH_digesta <- phy_tree(Physeq_VH_digesta)

# Importing metadata only for either mucus or digesta
meta_VH_mucus <- read.csv("sample_VH_mucus.csv", row.names = 1, sep = ";")
meta_VH_digesta <- read.csv("sample_VH_digesta.csv", row.names = 1, sep = ";")


pairdis.q1phy_mucus <- pair_dis(ASVs_VH_mucus, qvalue=1, hierarchy=meta_VH_mucus, tree = tree_VH_mucus, level = "1", metric = "U")
pairdis.q1phy_digesta <- pair_dis(ASVs_VH_digesta, qvalue=1, hierarchy=meta_VH_digesta, tree = tree_VH_digesta, level = "1", metric = "U")

L1_UqN.AI_mucus <- as.data.frame(pairdis.q1phy_mucus[["L1_UqN"]])
L1_UqN.AI_digesta <- as.data.frame(pairdis.q1phy_digesta[["L1_UqN"]])

L1_UqN.AI_mucus <- as.dist(L1_UqN.AI_mucus)
L1_UqN.AI_digesta <- as.dist(L1_UqN.AI_digesta)

# Checking assumptions of PERMANOVA
test.homogenity_VH_mucus=betadisper(L1_UqN.AI_mucus, group = meta_VH_mucus$Location)
test.homogenity_VH_digesta=betadisper(L1_UqN.AI_digesta, group = meta_VH_digesta$Location)
anova(test.homogenity_VH_mucus)
anova(test.homogenity_VH_digesta)

# Running PERMANOVA 2
permanova_VH_mucus<-adonis(L1_UqN.AI_mucus ~ Location, data=meta_VH_mucus, permutations = 999, na.rm = TRUE)
permanova_VH_digesta<-adonis(L1_UqN.AI_digesta ~ Location, data=meta_VH_digesta, permutations = 999, na.rm = TRUE)



# Overall diversity differences between species and sample types
### All species ###
ASVs <- read.csv("counts.csv", row.names = 1)
tax <- read.csv("taxonomy.csv", row.names = 1)
meta <- read.csv("sampledata.csv", row.names = 1)
tree <- read.tree("asvseqs.aligned.raxml.bestTree")
tree <- force.ultrametric(tree,method="extend")
tree <- match_data(ASVs,tree,output="tree")
match_data(ASVs,tree)

ASVs.physeq <- otu_table(ASVs, taxa_are_rows = TRUE)
tax.physeq <- as.matrix(tax)
tax.physeq <- tax_table(tax.physeq)
meta.physeq <- sample_data(meta)
Physeq = phyloseq(ASVs.physeq, tax.physeq, meta.physeq, tree)




Physeq_D <- subset_samples(Physeq, Sample.type=="Digesta")
Physeq_M <- subset_samples(Physeq, Sample.type=="Mucus")


# Test across species for digesta samples
ASVs_D <- as.data.frame(otu_table(Physeq_D))
tax_D <- as.data.frame(tax_table(Physeq_D))
meta_D <- as.data.frame(sample_data(Physeq_D))
tree_D <- phy_tree(Physeq_D)

meta_hill_D =meta[,c("SampleID", "Species")]
row.names(meta_hill_D) = meta_hill_D[,c("SampleID")]

div_test(ASVs_D, qvalue = 1, hierarchy = meta_hill_D, posthoc = T)

# Test across species for mucus samples
ASVs_M <- as.data.frame(otu_table(Physeq_M))
tax_M <- as.data.frame(tax_table(Physeq_M))
meta_M <- as.data.frame(sample_data(Physeq_M))
tree_M <- phy_tree(Physeq_M)

meta_hill_M =meta[,c("SampleID", "Species")]
row.names(meta_hill_M) = meta_hill_M[,c("SampleID")]

div_test(ASVs_M, qvalue = 1, hierarchy = meta_hill_M, posthoc = T)

