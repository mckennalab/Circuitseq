library("tidyverse")
library("data.table")
library("ggplot2")
library("ggpubr")

run_103 = fread("/analysis/2021_08_26_PlasmidSeq_paper/run_103/10_3_all.combined",fill=TRUE)
run_103 = run_103[run_103$assembly != "56_rotated.fasta",]
run_103 = run_103[order(run_103$plasmid),]

known_contamination_samples = run_103[run_103$plasmid >= 49 & run_103$plasmid <= 56,]
known_contamination_samples$rates = c(0,5,15,30,70,85) # ,95)

cor(known_contamination_samples$contamination,known_contamination_samples$rates)
rmse(known_contamination_samples$contamination,known_contamination_samples$rates)

g = ggplot(known_contamination_samples,aes(rates,contamination)) + 
  geom_point(size=2) + 
  theme_classic() + 
  xlab("Known contamination rate") + 
  ylab("Predicted contamination rate") +
  geom_smooth(method = "lm", se=TRUE, color="black", formula = y ~ x)
  
ggsave(g,file="/analysis/2021_08_26_PlasmidSeq_paper/figures/figure_4/known_contamination_experimental_2021_12_07.pdf",width=2,height=4)
ggsave(g,file="/analysis/2021_08_26_PlasmidSeq_paper/figures/figure_4/known_contamination_experimental_2021_12_07.png",width=2,height=4)
