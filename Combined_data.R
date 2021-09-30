setwd("~/workingdirectory")

library(dada2)
library(devtools)


# Merge multiple runs (if necessary)
st1 <- readRDS("~/st1.rds")
st2 <- readRDS("~/st2.rds")
st3 <- readRDS("~/st3.rds")
st4 <- readRDS("~/st4.rds")
st5 <- readRDS("~/st5.rds")

st.all <- mergeSequenceTables(st1, st2, st3, st4, st5)

# Remove chimeras
seqtab <- removeBimeraDenovo(st.all, method="consensus", multithread=TRUE)



#Assigning taxonomy(before this step, download the database, at: https://zenodo.org/record/1172783#.XzvzQJMzbmE)
#                  There are different DADA2-formatted databases available in DADA2 website

#Assign Taxanomy
taxa <- assignTaxonomy(seqtab, "silva_nr_v132_train_set.fa.gz", tryRC=T, multithread=TRUE)
saveRDS(taxa, "taxa.RData")

taxa <- readRDS("taxa.RData")

table(nchar(getSequences(seqtab)))
#to see the breakdown of the size of these amplicons. On the top row is the size of the merged reads, and on the botton is the frequency

#Chimera identification
seqtab.nochim <- removeBimeraDenovo(seqtab, multithread=TRUE, verbose=TRUE)
saveRDS(seqtab.nochim, "seqtab.RData")

### Create a Summary
# this is one quick way to look at sequences that have been lost, to know whether they held a lot in terms of abundance
sum(seqtab.nochim)/sum(seqtab)

#Overview of counts throughout:quick way to pull out how many reads were dropped at various points of the pipeline
# set a little function
getN <- function(x) sum(getUniques(x))

row.names.remove <- c("GH3_10a_2_1_trimmed.fq.gz")

filtered_out_1 <- filtered_out[!(row.names(filtered_out) %in% row.names.remove), ]

sample.names1 <-row.names(filtered_out_1)

# making a little table
summary_tab <- data.frame(row.names=sample.names1, dada2_input=filtered_out_1[,1],
                          filtered=filtered_out_1[,2], dada_f=sapply(dada_forward, getN),
                          dada_r=sapply(dada_reverse, getN), merged=sapply(merged_amplicons, getN),
                          nonchim=rowSums(seqtab.nochim),
                          final_perc_reads_retained=round(rowSums(seqtab.nochim)/filtered_out_1[,1]*100, 1))

write.table(summary_tab, "summary_reads_table.txt", sep="\t", quote=F)


#Extracting the standard goods from DADA2
#Write out tables for further processing}
#giving to seq headers more manageable names (ASV_1, ASV_2...)
asv_seqs <- colnames(seqtab)
asv_headers <- vector(dim(seqtab)[2], mode="character")

for (i in 1:dim(seqtab)[2]) {
  asv_headers[i] <- paste(">ASV", i, sep="_")
}

#making and writing out a fasta of our final ASV seqs:
asv_fasta <- c(rbind(asv_headers, asv_seqs))
write(asv_fasta, "ASVs.fa")

#count table:
asv_tab <- t(seqtab)
row.names(asv_tab) <- sub(">", "", asv_headers)
write.table(asv_tab, "ASVs_counts.tsv", sep="\t", quote=F, col.names=NA)

#tax table:
asv_tax <- taxa
row.names(asv_tax) <- sub(">", "", asv_headers)
write.table(asv_tax, "ASVs_taxonomy.tsv", sep="\t", quote=F, col.names=NA)


