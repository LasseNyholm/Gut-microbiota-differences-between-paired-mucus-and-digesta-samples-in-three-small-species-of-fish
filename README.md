# Gut microbiota differences between paired mucus and digesta samples in three small species of fish

Here follows the pipeline used for processing and statistical analysis of data. The different sub-pipelines are chronological ordered from 1-7. Step 1-5 originally included data that are not a part of the current study, but are important for DADA2 and decontaminations. This extra data can be supplied upon request, but for convenience we supply the reader with decontaminated tables (following decontamination protocol in Phyloseq_filtering.sh) found in the "Data" directory. 


1. Demultiplexing and preprocessing
- Raw sequences are demultiplexed based on unique tags, filtered and primers are trimmed

2. DADA2 - seperate rounds of sequencing (DADA2.R)
- Forward and reverse reads are filtered, ASVs are inferrred and forward and reverse reads merged. This is done for each round of sequencing to account for differences in sequencing quality between runs. 

3. DADA2 - total data (DADA2_combined_data.R)
- Datasets are merged, chimeric sequences removed, taxonomy is assigned and a ASV table and taxonomy table are generated.

4. LULU (LULU.sh)
- Removal erroneous ASVs

5. Decontamination, sorting and filtering (Phyloseq_filtering.sh)

6. Statistical analysis and visualization (Diversity_composition_plots.sh)

7. HMSC modelling (HMSC.sh)
