library("tidyverse")
library("data.table")
library("ggplot2")
library("ggpubr")
library("ggbeeswarm")
library("OneR")

read_counts = fread("/analysis/2021_08_26_PlasmidSeq_paper/figures/sup_figure_3/2021_12_14_read_accounting.txt")
sample_sheet = fread("/analysis/2021_08_26_PlasmidSeq_paper/fex_117_sample_sheet.tsv")
repeat_content = fread("/analysis/2021_08_26_PlasmidSeq_paper/repeats/MUMmer3.23/2021_12_13_repeat_content.txt")

merged_sheet = left_join(read_counts,sample_sheet,by=c("sample"="position"))
merged_sheet_repeat = left_join(merged_sheet,repeat_content,by=c("sample"="sample"))


g = ggplot(merged_sheet_repeat,aes(length.x,guppyReads,fill=unrepeatProp)) + geom_point(color="black",size=3,shape=21,stroke =1) + 
  theme_classic() + 
  xlab("Plasmid length") + 
  ylab("Demultiplexed reads") + 
  guides(fill = guide_legend(title="Repeat content")) +
  geom_smooth(method = "lm", se=TRUE, color="black", formula = y ~ x)

ggsave(g,file="/analysis/2021_08_26_PlasmidSeq_paper/figures/sup_figure_3/2021_12_14_reads_by_plasmid_length.png",width=4,height=4)
ggsave(g,file="/analysis/2021_08_26_PlasmidSeq_paper/figures/sup_figure_3/2021_12_14_reads_by_plasmid_length.pdf",width=4,height=4)

read_sizes = fread("/analysis/2021_08_26_PlasmidSeq_paper/figures/sup_figure_3/2021_12_14_read_sizes.txt")
merged_sheet_repeat_sizes = full_join(merged_sheet_repeat,read_sizes)
merged_sheet_repeat_sizes$normalized_read_length = merged_sheet_repeat_sizes$readLength / merged_sheet_repeat_sizes$length.x
merged_sheet_repeat_sizes$size = bin(merged_sheet_repeat_sizes$length.x,nbins=3,labels=c("small","medium","large"))
merged_sheet_repeat_sizes$size_raw = bin(merged_sheet_repeat_sizes$length.x,nbins=3)

g = ggplot(merged_sheet_repeat_sizes[merged_sheet_repeat_sizes$typeof == "raw",],aes(normalized_read_length)) +
  geom_histogram() +
  theme_classic() + 
  xlab("Normalized read length") + 
  ylab("Count") +
  xlim(c(0,1.2))

ggsave(g,file="/analysis/2021_08_26_PlasmidSeq_paper/figures/sup_figure_3/2021_12_14_reads_length_filtered.png",width=4,height=4)
ggsave(g,file="/analysis/2021_08_26_PlasmidSeq_paper/figures/sup_figure_3/2021_12_14_reads_length_filtered.pdf",width=4,height=4)


g = ggplot(merged_sheet_repeat_sizes[merged_sheet_repeat_sizes$typeof == "raw",],aes(normalized_read_length)) +
  geom_histogram() +
  theme_classic() + 
  xlab("Normalized read length") + 
  ylab("Count") +
  xlim(c(0,1.2)) + 
  facet_wrap(. ~ size)

ggsave(g,file="/analysis/2021_08_26_PlasmidSeq_paper/figures/sup_figure_3/2021_12_14_reads_length_filtered_by_size.png",width=6,height=4)
ggsave(g,file="/analysis/2021_08_26_PlasmidSeq_paper/figures/sup_figure_3/2021_12_14_reads_length_filtered_by_size.pdf",width=6,height=4)


g = ggplot(merged_sheet_repeat_sizes[merged_sheet_repeat_sizes$typeof == "correct",],aes(normalized_read_length)) +
  geom_histogram() +
  theme_classic() + 
  xlab("Normalized read length") + 
  ylab("Count") +
  xlim(c(0,1.2))

ggsave(g,file="/analysis/2021_08_26_PlasmidSeq_paper/figures/sup_figure_3/2021_12_14_reads_length_correct.png",width=4,height=4)
ggsave(g,file="/analysis/2021_08_26_PlasmidSeq_paper/figures/sup_figure_3/2021_12_14_reads_length_correct.pdf",width=4,height=4)


g = ggplot(merged_sheet_repeat_sizes[merged_sheet_repeat_sizes$typeof == "correct",],aes(normalized_read_length)) +
  geom_histogram() +
  theme_classic() + 
  xlab("Normalized read length") + 
  ylab("Count") +
  xlim(c(0,1.2)) + 
  facet_wrap(. ~ size)

ggsave(g,file="/analysis/2021_08_26_PlasmidSeq_paper/figures/sup_figure_3/2021_12_14_reads_length_correct_by_size.png",width=6,height=4)
ggsave(g,file="/analysis/2021_08_26_PlasmidSeq_paper/figures/sup_figure_3/2021_12_14_reads_length_correct_by_size.pdf",width=6,height=4)

