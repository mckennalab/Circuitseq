library(data.table)
library(ggplot2)
library(tidyverse)
library(ggExtra)
options(scipen=1000000)

# data files for assessing plasmid length
# ----- base_directory = "/analysis/2021_07_12_PlasmidSeq_publication_96xr941/PlasmidSeq" ----- 
base_image_directory = "/analysis/2021_08_26_PlasmidSeq_paper/images/2021_08_26_plasmidseq_paper"
barcoding_summary_file = fread("/analysis/2021_07_12_PlasmidSeq_publication_96xr941/PlasmidSeq/pipelines/results/guppy_demultiplex/saved_data/barcoding_summary.txt")
read_sequencing_file = fread("/analysis/2021_07_12_PlasmidSeq_publication_96xr941/PlasmidSeq/pipelines/results/guppy/basecalling/sequencing_summary.txt")
sample_sheet = fread("/analysis/2021_08_26_PlasmidSeq_paper/fex_117_sample_sheet_with_annotations.tsv")

# ------------------------------------------------------------------------------------------------
# Plot some basic stats about the plasmids 
# ------------------------------------------------------------------------------------------------

# we used a kmer of 10
sample_sheet$unique_kmer_count = sample_sheet$kmerCount / (sample_sheet$length - 10)

# plot a plasmid length histogram 
plasmid_lengths = 
  ggplot(sample_sheet) + 
  geom_histogram(aes(length)) + 
  xlab("Known plasmid length") + 
  ylab("Count") + 
  scale_y_continuous(breaks = seq(0, 25, by = 1)) + 
  scale_x_continuous(breaks = seq(0, 15000, by = 1000)) + 
  theme_classic()
ggsave(plasmid_lengths,file=paste(base_image_directory,"2021_08_26_96x_plasmid_length_histogram.png",sep="/"),width=6,height=6)

# plot the kmer repeats per plasmid
tenmer_kmer_repeat_plot = 
  ggplot(sample_sheet) + 
  geom_histogram(aes(unique_kmer_count)) + 
  xlab("Unique 10-mer proportion") + 
  ylab("Count") + 
  scale_y_continuous(breaks = seq(0, 25, by = 1)) + 
  xlim(c(0.8,1.0)) +
  theme_classic()
ggsave(tenmer_kmer_repeat_plot,file=paste(base_image_directory,"2021_08_26_96x_plasmid_kmer_histogram.png",sep="/"),width=6,height=6)

# The GC proportion
gc_plot = 
  ggplot(sample_sheet) + 
  geom_histogram(aes(gcPercentage)) + 
  xlab("GC proportion") + 
  ylab("Count") + 
  scale_y_continuous(breaks = seq(0, 25, by = 1)) + 
  theme_classic()
ggsave(gc_plot,file=paste(base_image_directory,"2021_08_26_96x_gc_proportion.png",sep="/"),width=6,height=6)

# ------------------------------------------------------------------------------------------------
# join the reads to the barcode stats and to the sample table
# ------------------------------------------------------------------------------------------------

joined_read_file = full_join(barcoding_summary_file,read_sequencing_file,by=c("read_id"="read_id"))
joined_read_file$filtered = TRUE
joined_read_file[complete.cases(joined_read_file),]$filtered = FALSE
joined_read_file_complete = joined_read_file[complete.cases(joined_read_file),]

joined_with_samples = full_join(joined_read_file_complete,sample_sheet,by=c("barcode_arrangement"="barcodeID"))
joined_with_samples_complete = joined_with_samples[complete.cases(joined_with_samples),]

# ------------------------------------------------------------------------------------------------
# Plot stats about the reads
# ------------------------------------------------------------------------------------------------

barcode_assignment = 
  ggplot(joined_read_file_complete) + 
  geom_violin(aes(y=barcode_score,x=barcode_arrangement != "unclassified")) + 
  ylab("Barcode assignment score") +
  xlab("Read assigned to a plasmid") +
  theme_classic()
ggsave(barcode_assignment,file=paste(base_image_directory,"2021_08_26_96x_barcode_assignment_scores.png",sep="/"),width=6,height=6)

# what's the barcode assignment rate
print(dim(joined_read_file_complete))
print(summary(joined_read_file_complete$barcode_arrangement != "unclassified"))

# sequence_length_template contains the acutal read length
joined_with_samples_complete$normalized_read_length = joined_with_samples_complete$sequence_length_template / joined_with_samples_complete$length
normalized_read_length_plot = 
  ggplot(joined_with_samples_complete) + 
  geom_histogram(aes(normalized_read_length)) + 
  xlim(c(0,2.5)) +
  ylab("Count") +
  xlab("Normalized read length") +
  theme_classic()
ggsave(normalized_read_length_plot,file=paste(base_image_directory,"2021_08_26_96x_normalized_read_length.png",sep="/"),width=6,height=6)

normalized_read_length_plot_facet = 
  ggplot(joined_with_samples_complete) + 
  geom_histogram(aes(normalized_read_length)) + 
  xlim(c(0,2.5)) +
  ylab("Count") +
  xlab("Normalized read length") +
  theme_classic() + 
  facet_wrap(. ~ plasmid_size_quantile,nrow = 1, labeller = labeller(plasmid_size_quantile = supp.labs))
ggsave(normalized_read_length_plot_facet,file=paste(base_image_directory,"2021_08_26_96x_normalized_read_length_by_quantile.png",sep="/"),width=12,height=2)

filtered_read_lengths = fread("/analysis/2021_08_26_PlasmidSeq_paper/fex_117_read_lengths.tsv")
filtered_read_lengths_joined = full_join(filtered_read_lengths,sample_sheet,by=c("well"="barcodeID"))
filtered_read_lengths_joined$normalized_read_length = filtered_read_lengths_joined$readlength / filtered_read_lengths_joined$length

filtered_normalized_read_length_plot = 
  ggplot(filtered_read_lengths_joined) + 
  geom_histogram(aes(normalized_read_length)) + 
  xlim(c(0,1.5)) +
  ylab("Count") +
  xlab("Normalized read length") +
  theme_classic()
ggsave(filtered_normalized_read_length_plot,file=paste(base_image_directory,"2021_08_26_96x_normalized_filtered_read_length.png",sep="/"),width=6,height=6)


# 
final_assembly_stats = fread("/analysis/2021_07_12_PlasmidSeq_publication_96xr941/PlasmidSeq/pipelines/results/plasmid_stat/all_rotated.stats")

# classic plot :
p_global <- ggplot(final_assembly_stats, aes(x=contiguity, y=identity)) +
  geom_point() +
  theme(legend.position="none") + 
  theme_classic() + 
  xlim(c(0,1.5)) + 
  ylim(c(0,1.0))

# Set relative size of marginal plots (main plot 10x bigger than marginals)
p1_global <- ggMarginal(p_global, type="histogram", size=10)
ggsave(p1_global,file=paste(base_image_directory,"2021_08_26_96x_global_identity_contiguity.png",sep="/"),width=6,height=6)

# classic plot :
p_zoom <- ggplot(final_assembly_stats, aes(x=contiguity, y=identity)) +
  geom_point() +
  theme(legend.position="none") + 
  theme_classic() 

# Set relative size of marginal plots (main plot 10x bigger than marginals)
p1_zoom <- ggMarginal(p_zoom, type="histogram", size=10)
ggsave(p1_zoom,file=paste(base_image_directory,"2021_08_26_96x_zoom_identity_contiguity.png",sep="/"),width=6,height=6)
