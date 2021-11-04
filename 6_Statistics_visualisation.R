setwd("~/Workingdirectory")

library(hilldiv)
library(nlme)
library(phytools)
library(dplyr)
library(ggpubr)
library(phyloseq)
library(fantaxtic)
library(vegan)


# Aphanius iberus: mucus and content # 
ASVs_AI <- read.csv("count_AI.csv", row.names = 1)
tax_AI <- read.csv("tax_AI.csv", row.names = 1)
meta_AI <- read.csv("sample_AI.csv", row.names = 1)
tree_AI <- read.tree("tree_AI")
meta_AI_hill =meta_AI[,c("SampleID", "Location")]
row.names(meta_AI_hill) = meta_AI_hill[,c("SampleID")]
# 
# Phylooseq object
ASVs_AI.physeq <- otu_table(ASVs_AI, taxa_are_rows = TRUE)
tax_AI.physeq <- as.matrix(tax_AI)
tax_AI.physeq <- tax_table(as.matrix(tax_AI))
meta_AI.physeq <- sample_data(meta_AI)
Physeq_AI = phyloseq(ASVs_AI.physeq, tax_AI.physeq, meta_AI.physeq, tree_AI)

#Diversity profiles
AI_div_profil_AI <- div_profile(ASVs_AI, qvalues=seq(from = 0, to = 3, by = (0.2)), hierarchy = meta_AI_hill)
div_profile_plot(AI_div_profil_AI)


#Diversity computation
q1phy_AI <- hill_div(ASVs_AI,qvalue=1 ,tree=tree_AI)

#Linear mixed-effects model
model.q1phy_AI <- lme(log(q1phy_AI)~1+Sample.type, random = ~1|Location/Fish_ID, meta_AI)
summary(model.q1phy_AI)
plot(model.q1phy_AI)
qqnorm(model.q1phy_AI)

#Compositional variation
pairdis.q1phy <- pair_dis(ASVs_AI, qvalue=1, hierarchy=meta_AI, tree = tree_AI, level = "1", metric = "U")

L1_UqN.AI <- as.data.frame(pairdis.q1phy[["L1_UqN"]])

L1_UqN.AI <- as.dist(L1_UqN.AI)

GP.ord <- ordinate(Physeq_AI, "NMDS", distance = L1_UqN.AI)

ordi_AI  <- plot_ordination(Physeq_AI, GP.ord, type="samples", color="Location", shape = "Sample.type")+
  geom_line(aes(group = Fish_ID), alpha=0.4)+
  theme_bw()+geom_point(size = 3)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank())+
  ggtitle("Aphanius iberus")+
  theme(plot.title = element_text(hjust=0.5))#+


#Checking assumptions of PERMANOVA

test.homogenity_AI=betadisper(L1_UqN.AI, group = meta_AI$Location)

anova(test.homogenity_AI)

permanova_AI<-adonis(L1_UqN.AI ~ Sample.type*Location, data=meta_AI, permutations = 999, na.rm = TRUE, strata =meta_AI$Fish_ID)


# Gambusia holbrooki # 
ASVs_GH <- read.csv("count_GH.csv", row.names = 1)
tax_GH <- read.csv("tax_GH.csv", row.names = 1)
meta_GH <- read.csv("sample_GH.csv", row.names = 1)
tree_GH <- read.tree("tree_GH")
meta_GH_hill =meta_GH[,c("SampleID", "Location")]
row.names(meta_GH_hill) = meta_GH_hill[,c("SampleID")]

# Phylooseq object
ASVs_GH.physeq <- otu_table(ASVs_GH, taxa_are_rows = TRUE)
tax_GH.physeq <- as.matrix(tax_GH)
tax_GH.physeq <- tax_table(as.matrix(tax_GH))
meta_GH.physeq <- sample_data(meta_GH)
Physeq_GH = phyloseq(ASVs_GH.physeq, tax_GH.physeq, meta_GH.physeq, tree_GH)

#Diversity profiles
div_profil_GH <- div_profile(ASVs_GH, qvalues=seq(from = 0, to = 3, by = (0.2)), hierarchy = meta_GH_hill)
div_profile_plot(div_profil_GH)

#Diversity tests
q1phy_GH <- hill_div(ASVs_GH,qvalue=1, tree=tree_GH)

# Linear mixed-effects  model
model.q1phy_GH <- lme(q1phy_GH~1+Sample.type, random = ~1|Location/Fish_ID, meta_GH)
summary(model.q1phy_GH)
plot(model.q1phy_GH)
qqnorm(model.q1phy_GH)

# Compositional variation
pairdis.q1phy_GH <- pair_dis(ASVs_GH, qvalue=1, hierarchy=meta_GH, tree = tree_GH, level = "1", metric = "U")

L1_UqN.GH <- as.data.frame(pairdis.q1phy_GH[["L1_UqN"]])

L1_UqN.GH <- as.dist(L1_UqN.GH)

GP.ord <- ordinate(Physeq_GH, "NMDS", distance = L1_UqN.GH)

ordi_GH  <- plot_ordination(Physeq_GH, GP.ord, type="samples", color="Location", shape = "Sample.type")+
  geom_line(aes(group = Fish_ID), alpha=0.4)+
  theme_bw()+geom_point(size = 3)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank())+
  ggtitle("Gambusia holbrooki")+
  theme(plot.title = element_text(hjust=0.5))

#Checking assumptions of PERMANOVA

test.homogenity_GH=betadisper(L1_UqN.GH, group = meta_GH$Sample.type)

anova(test.homogenity_GH)

permanova_GH<-adonis(L1_UqN.GH ~ Sample.type*Location, data=meta_GH, permutations = 999, na.rm = TRUE, strata =meta_GH$Fish_ID)

# Valencia hispanica # 
ASVs_VH <- read.csv("count_VH.csv", row.names = 1)
tax_VH <- read.csv("tax_VH.csv", row.names = 1)
meta_VH <- read.csv("sample_VH.csv", row.names = 1)
tree_VH <- read.tree("tree_VH")
meta_VH_hill =meta_VH[,c("SampleID", "Location")]
row.names(meta_VH_hill) = meta_VH_hill[,c("SampleID")]

# Phylooseq object
ASVs_VH.physeq <- otu_table(ASVs_VH, taxa_are_rows = TRUE)
tax_VH.physeq <- as.matrix(tax_VH)
tax_VH.physeq <- tax_table(as.matrix(tax_VH))
meta_VH.physeq <- sample_data(meta_VH)
Physeq_VH = phyloseq(ASVs_VH.physeq, tax_VH.physeq, meta_VH.physeq, tree_VH)

#Diversity profiles
div_profil_VH <- div_profile(ASVs_VH, qvalues=seq(from = 0, to = 3, by = (0.2)), hierarchy = meta_VH_hill)
div_profile_plot(div_profil_VH)
 
#Diversity computation
q1phy_VH <- hill_div(ASVs_VH,qvalue=1, tree=tree_VH)

#Linear mixed-effects model
model.q1phy_VH <- lme(q1phy_VH~1+Sample.type, random = ~1|Location/Fish_ID, meta_VH)
summary(model.q1phy_VH)
plot(model.q1phy_VH)
qqnorm(model.q1phy_VH)


#Compositional variation
pairdis.q1phy_VH <- pair_dis(ASVs_VH, qvalue=1, hierarchy=meta_VH, tree = tree_VH, level = "1", metric = "U")

L1_UqN.VH <- as.data.frame(pairdis.q1phy_VH[["L1_UqN"]])

L1_UqN.VH <- as.dist(L1_UqN.VH)

GP.ord <- ordinate(Physeq_VH, "NMDS", distance = L1_UqN.VH)

ordi_VH  <- plot_ordination(Physeq_VH, GP.ord, type="samples", color="Location", shape = "Sample.type")+
  geom_line(aes(group = Fish_ID), alpha=0.4)+
  theme_bw()+geom_point(size = 3)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank())+
  ggtitle("Valencia hispanica")+
  theme(plot.title = element_text(hjust=0.5))


#Checking assumptions of PERMANOVA

test.homogenity_GH=betadisper(L1_UqN.VH, group = meta_VH$Location)

anova(test.homogenity_GH)

permanova_VH<-adonis(L1_UqN.VH ~ Sample.type*Location, data=meta_VH, permutations = 999, na.rm = TRUE, strata =meta_VH$Fish_ID)

pdf("ordination.pdf",width=6,height=12)

ggarrange(ordi_AI, ordi_GH, ordi_VH,
          ncol = 1, nrow = 3)
dev.off()

#### Barplot for all three species
Physeq_AI_15 = get_top_taxa(Physeq_AI, 15, relative = TRUE, discard_other = FALSE,
                            other_label = "Other")

ASVs_AI_15 <- otu_table(Physeq_AI_15)
tax_AI_15 <- tax_table(Physeq_AI_15)
sam_AI_15 <- sample_data(Physeq_AI_15)

Physeq_GH_15 = get_top_taxa(Physeq_GH, 15, relative = TRUE, discard_other = FALSE,
                            other_label = "Other")

ASVs_GH_15 <- otu_table(Physeq_GH_15)
tax_GH_15 <- tax_table(Physeq_GH_15)
sam_GH_15 <- sample_data(Physeq_GH_15)

Physeq_VH_15= get_top_taxa(Physeq_VH, 15, relative = TRUE, discard_other = T,
                           other_label = "Other")


ASVs_VH_15 <- otu_table(Physeq_VH_15)
tax_VH_15 <- tax_table(Physeq_VH_15)
sam_VH_15 <- sample_data(Physeq_VH_15)

ASVs_VH_15 <- as.data.frame(otu_table(Physeq_VH_15))

colSums(otu_table(Physeq_AI_15))

Physeq_all_ASVs <- merge_phyloseq(ASVs_AI_15, ASVs_GH_15, ASVs_VH_15)
Physeq_all_tax <- merge_phyloseq(tax_AI_15, tax_GH_15, tax_VH_15)
Physeq_all_samples <- merge_phyloseq(sam_AI_15, sam_GH_15, sam_VH_15)


Physeq_all <- merge_phyloseq(Physeq_all_ASVs, Physeq_all_tax, Physeq_all_samples)

Sample= as.data.frame(sample_data(Physeq_all))

vektor = row.names(Sample)

pdf("Barplots_family.pdf",width=10,height=12)

p = plot_bar(Physeq_all, fill = "Family") + facet_wrap(~Species, scales="free_x", nrow = 3)+ theme(legend.position="left",strip.text = element_text(size=25), panel.background = element_blank())+
  labs(x="Sample", y="Relative abundance")+geom_bar(stat="identity")+labs(x="Sample", y="Relative abundance")+geom_bar(stat="identity")+theme(strip.background =element_rect(fill="white"))
p$data$Sample<- factor(p$data$Sample, levels = vektor)

dev.off()

