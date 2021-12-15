library("tidyverse")
library("data.table")
library("ggplot2")
library("ggpubr")
library(ggbeeswarm)

contam = fread("/analysis/2021_08_26_PlasmidSeq_paper/fex_117_sample_sheet.tsv.combined")
contam$complete = (1.0 - abs(1.0 - contam$contiguity)) * (1.0 - abs(1.0 - contam$identity))
contam$run = "9.4"

run_103 = fread("/analysis/2021_08_26_PlasmidSeq_paper/run_103/10_3_all.combined",fill=TRUE)
run_103 = run_103[run_103$assembly != "56_rotated.fasta",]
run_103 = run_103[order(run_103$plasmid),]
run_103$complete = (1.0 - abs(1.0 - run_103$contiguity)) * (1.0 - abs(1.0 - run_103$identity))
run_103$run = "10.3"

full_94_103_run = rbind(contam,run_103)

g = ggplot(full_94_103_run,aes(complete,contamination,fill=run)) + geom_point(color="black",size=3,shape=21,stroke =1) + 
  theme_classic() + 
  ylab("Contamination rate") + 
  xlab("Completeness") + 
  xlim(c(0,1)) + 
  ylim(c(0,100))  
  # geom_smooth(method = "lm", se=TRUE, color="black", formula = y ~ x)


ggsave(g,file="/analysis/2021_08_26_PlasmidSeq_paper/figures/figure_4/real_94_103_contamination_assembly_2021_11_07.pdf",width=4,height=4)
ggsave(g,file="/analysis/2021_08_26_PlasmidSeq_paper/figures/figure_4/real_96_103contamination_assembly_2021_11_07.png",width=4,height=4)

run_94_103_cont_ident = melt(full_94_103_run,id.vars = c("plasmid","run"), measure.vars = c("contiguity","identity"))
ggplot(run_94_103_cont_ident,aes(y=(1.0 - abs(1.0 - value)),x=variable,col=run)) + geom_quasirandom()


run_94_cont_ident = melt(contam,id.vars = c("plasmid"), measure.vars = c("contiguity","identity"))

g = ggplot(run_94_cont_ident,aes(y=(1.0 - abs(1.0 - value)),x=variable,fill=variable)) + 
  geom_quasirandom(size=3,shape=21,stroke =1,col="black") + 
  theme_classic() + 
  xlab("Value") + 
  scale_colour_manual(values = c("#0066ff","#ff5050"))

ggsave(g,file="/analysis/2021_08_26_PlasmidSeq_paper/figures/figure_1/contiguity_identity_2021_12_08.pdf",width=6,height=6)

g = ggplot(run_94_cont_ident[(1.0 - abs(1.0 - run_94_cont_ident$value)) > 0.995,],
           aes(y=(1.0 - abs(1.0 - value)),x=variable,fill=variable)) + 
  geom_quasirandom(size=3,shape=21,stroke =1,col="black") + 
  theme_classic() + 
  xlab("Value") + 
  scale_colour_manual(values = c("#0066ff","#ff5050"))

ggsave(g,file="/analysis/2021_08_26_PlasmidSeq_paper/figures/figure_1/contiguity_identity_subset_2021_12_08.pdf",width=6,height=6)


