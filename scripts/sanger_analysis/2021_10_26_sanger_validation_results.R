library("ggplot2")
library("data.table")
library("tidyverse")
library('ggbeeswarm')
setwd("/Users/aaronmck/Desktop/ab1_trial/PlasmidSeq/scripts/sanger_analysis/")

# read the summary file from the plasmid alignment
tbl = fread("summary_file_2021_10_24.txt")

# calculate some additional stats
tbl$sangerCoverage = tbl$sangerLength/tbl$plasmidLength
tbl$error_rate = tbl$sangerErrors/tbl$sangerLength
tbl$errorPerPlasmid = tbl$plasmidLength * tbl$error_rate

gg = ggplot(tbl) + geom_histogram(aes(error_rate)) + theme_classic() + xlab("Error rate") + ylab("Count")
ggsave(gg,file="plots/2021_10_25_full_error_rate.png",width=3,height=3)

# filter out the high error rate sanger results
filtered.table = tbl[tbl$error_rate < .1,]

gg = ggplot(filtered.table) + geom_histogram(aes(error_rate)) + theme_classic() + xlab("Filtered error rate") + ylab("Count")
ggsave(gg,file="plots/2021_10_25_filtered_error_rate.png",width=3,height=3)


gg = ggplot(filtered.table) + geom_histogram(aes(sangerCoverage)) + theme_classic() + xlab("Plasmid coverage") + ylab("Count")
ggsave(gg,file="plots/2021_10_25_coverage.png",width=3,height=3)

# compare to the existing assemblies for those plasmids
nanopore = fread("all_rotated.stats")
nanopore$id = sapply(nanopore$assembly,function(x) { return(as.numeric(substr(x,1,2)))})
nanopore$error_rate = 1.0 - nanopore$identity

nano_combined = nanopore[,c("id","error_rate")]
nano_combined$approach = "Nanopore"
filtered_combined = filtered.table[,c("id","error_rate")]
filtered_combined$approach = "Sanger"
combined = rbind(nano_combined,filtered_combined)

gg = ggplot(combined) + 
  geom_beeswarm(aes(y=error_rate,x=approach,fill=approach)) + 
  theme_classic() + 
  xlab("Sequencing Approach") + 
  ylab("Error rate") + 
  scale_fill_manual(values=c("#0072B2", "#56B4E9")) + 
  theme(legend.position = "none")
ggsave(gg,file="plots/2021_10_25_error_rate_by_tech.png",width=3,height=3)

combined %>% group_by(approach) %>% summarise(mean = median(error_rate), n = n())
#

nano_filtered = nano_combined[is.element(nano_combined$id,filtered_combined$id),]
combined_filtered = rbind(nano_filtered,filtered_combined)
combined_filtered %>% group_by(approach) %>% summarise(mean = median(error_rate), n = n())

gg = ggplot(combined_filtered) + 
  geom_beeswarm(aes(y=error_rate,x=approach,col=approach),size=0.3) + 
  theme_classic() + 
  xlab("Sequencing Approach") + 
  ylab("Error rate") + 
  scale_color_manual(values=c("#EE7733", "#0077BB")) + 
  theme(legend.position = "none")
ggsave(gg,file="plots/2021_10_25_filtered_error_rate_by_tech.png",width=3,height=3)

# create the contiguity and identity rate plot for the nanopore side
nanopore_melt = melt(nanopore,id.vars=c("assembly"),measure.vars = c("contiguity","identity"))

gg = ggplot(nanopore_melt,aes(y=value,x=variable,col=variable)) + 
  geom_boxplot() + 
  theme_classic() +
  xlab("") + 
  ylab("Error rate") + 
  scale_color_manual(values=c("#009E73", "#CC79A7")) + 
  theme(legend.position = "none")
ggsave(gg,file="plots/2021_10_26_full_view_contiguity_identity.pdf",width=5,height=3)

gg = ggplot(nanopore_melt,aes(y=value,x=variable,col=variable)) + 
  geom_beeswarm(size=0.5) + 
  theme_classic() +
  ylim(c(0.99,1.0)) + 
  xlab("") + 
  ylab("Error rate") + 
  scale_color_manual(values=c("#009E73", "#CC79A7")) + 
  theme(legend.position = "none")
ggsave(gg,file="plots/2021_10_26_zoom_view_contiguity_identity.pdf",width=5,height=3)
