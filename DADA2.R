# Pipeline based on the online DADA2 Pipeline Tutorial (1.16) (Callahan et al. 2016)
# The first part of the pipeline deals with each sequencing run separately to infer specific error rates and
# possible overlap when merging forward and reverse reads. 


setwd("~/path/to/trimmed/data")
library(dada2)
list.files() # make sure what we think is here is actually here

## first we're setting a few variables we're going to use ##
# one with all sample names, by scanning our "samples" file we made earlier
path <- "~/path/to/trimmed/data"

# one holding the file names of all the forward reads
fnFs <- sort(list.files(path, pattern="_1_trimmed.fq.gz", full.names = TRUE))
# and one with the reverse
fnRs <- sort(list.files(path, pattern="_2_trimmed.fq.gz", full.names = TRUE))
# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(basename(fnFs), "_1_trimmed"), `[`, 1)

#Quality trimming/filtering
pdf("QC_profile_2018_03_08B.pdf", height = 15, width = 15)
plotQualityProfile(fnFs[1:6]) 
plotQualityProfile(fnRs[1:6]) 
dev.off()

# Description for filtering, using DADA2
#INPUT: forward_reads, reverse_reads
#OUTPUT: filtered_forward_reads, filtered_reverse_reads:
# maxEE=Maximum amount of estimated errors that you expect to see and read
#rm.phix: removes any reads that match the PhiX bacteriophage genome, which is typically added to Illumina sequencing runs for quality monitoring
#truncLen: parameter setting the minimum size to trim the forward and reverse reads to in order to keep the quality scores roughly above 30 overall
#minLen: is setting the minimum length reads we want to keep after trimming
#no minLen value in dada2 tutorial
#maxN: additional filtering default parameter that is removing any sequences containing any Ns
#truncq= truncates your reads at the first base pair with a quality score equal to or less than 2
#compress=TRUE: FASTQ files gzipped
#multithread=TRUE: FASTQ files process in parallel

#Filtering
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names

# truncLen is specific for each sequencing run depending on the quality of sequences

filtered_out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE, truncLen=c(223,223), compress=TRUE, multithread=TRUE, matchIDs = TRUE)
head(filtered_out)
### Store filtered reads as R dataset
saveRDS(filtered_out, "filtered_out.RData")
#r filtering based on error rates
#Generating an error model of our data by learning the specific error-signature of our dataset
err_forward_reads <- learnErrors(filtFs, multithread=TRUE)
err_reverse_reads <- learnErrors(filtRs, multithread=TRUE)

#The developers have incorporated a plotting function to visualize how well the estimated error rates match up with the observed:
#The red line is what is expected based on the quality score, the black line represents the estimate, and the black dots represent the observed.
#In generally speaking, you want the observed (black dots) to track well with the estimated (black line)
#In general, error rate decreases as Q scores increases
#QC filtered reads
pdf("plotErrors.pdf", height = 15, width = 15)
plotErrors(err_forward_reads, nominalQ=TRUE)
plotErrors(err_reverse_reads, nominalQ=TRUE)
dev.off()
#Dereplication:keep/process one, and just attach the identical sequences to it
#When DADA2 dereplicates sequences, it also generates a new quality-score profile of each unique sequence
#based on the average quality scores of each base of all of the sequences that were replicates of it.
#the sample names in these objects are initially the file names of the samples, this sets them to the sample names for the rest of the workflow

#QC filtered reads
exists <- file.exists(filtFs)
derep_forward <- derepFastq(filtFs[exists], verbose=TRUE)
derep_reverse <- derepFastq(filtRs[exists], verbose=TRUE)
names(derep_forward) <- sample.names[exists]
names(derep_reverse) <- sample.names[exists]

#Inferring of ASVs
#It does this by incorporating the consensus quality profiles and abundances of each unique sequence, and then figuring out if each sequence more likely to be of biological origin or more likely to be spurious.

#Inferring ASVs
dada_forward <- dada(derep_forward, err=err_forward_reads, multithread=TRUE)
dada_reverse <- dada(derep_reverse, err=err_reverse_reads, multithread=TRUE)

#maxMismatch is by default 0 but you can make it say 1 or 2 if you’re finding that a lot of your forward and reverse reads are not merging
#minOverlap option has been change because they were sequenced at 250 bp


#Merging forward and reverse reads. minOverlap depends on the quality of the sequncing run, but never goes below 5 bp.
merged_amplicons <- mergePairs(dada_forward, derep_forward, dada_reverse,
                               derep_reverse, minOverlap=5, verbose=TRUE)
### Store filtered amplicons as R dataset
saveRDS(merged_amplicons, "merged_amplicons.RData")

###Generate ASV table
#Generating a count table
seqtab <- makeSequenceTable(merged_amplicons)

# Extract rds file to merge with other datasets from different sequencing runs
saveRDS(seqtab, file = "~/st1.rds", ascii = FALSE, version = NULL,
        
        compress = TRUE, refhook = NULL)
