# Gut microbiota differences between paired mucus and digesta samples in three small species of fish
16S amplicon pipeline

Here follows the pipeline used for processing and statistical analysis of data. The different sub-pipelines are chronological ordered from 1-7.


1. Demultiplexing and preprocessing
- Raw sequences are demultiplexed based on unique tags, filtered and primers are trimmed

2. DADA2 - seperate rounds of sequencing
- Forward and reverse reads are filtered, ASVs are inferrred and forward and reverse reads merged. This is done for each round of sequencing to account for differences in sequencing quality between runs. 

3. DADA2 - total data
- Datasets are merged, chimeric sequences removed, taxonomy is assigned and a ASV table and taxonomy table are generated.

4. LULU

5. Decontamination, sorting and filtering (Phyloseq_filtering.sh)

6. Statistical analysis and visualization (Diversity_composition_plots.sh)

7. HMSC modelling (HMSC.sh)
