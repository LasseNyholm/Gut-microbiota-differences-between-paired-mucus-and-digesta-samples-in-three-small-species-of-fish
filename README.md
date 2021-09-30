# Gut microbiota differences between paired mucus and digesta samples in three small species of fish
16S amplicon pipeline

Here follows the pipeline used for processing and statistical analysis of data. The different sub-pipelines are chronological ordered from 1-5.


1. Demultiplexing and preprocessing
- Raw sequences are demultiplexed based on unique tags, filtered and primers are trimmed

2. DADA2
- Forward and reverse reads 

3. Decontamination, sorting and filtering (Phyloseq_filtering.sh)

4. Statistical analysis and visualization (Diversity_composition_plots.sh)

5. HMSC modelling (HMSC.sh)
