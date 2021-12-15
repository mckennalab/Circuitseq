library("tidyverse")
library("data.table")
library("ggplot2")
tbl = fread("/analysis/2021_08_26_PlasmidSeq_paper/scripts/contamination/computational_96_well_mixture_addgene_with_contam_2021_12_06_fragment_adj.txt")

full_tbl = tbl %>% 
  group_by(plasmid,contam_set_rate) %>%
  mutate(cmedian=median(contamination, na.rm=TRUE),mad=mad(contamination, na.rm=TRUE))


full_tbl$plasmid = as.factor(full_tbl$plasmid)
sizes = c("12251" = 8890,"16337"=7397,"31815"=4290,"49792"=12598,"52961"=14873)
names = c("12251" = "medium","16337"="small","S"="XS","49792"="L","52961"="XL")
full_tbl$plasmid_id = "UNKNOWN"
full_tbl$plasmid_id[full_tbl$plasmid == "12251"] = "M"
full_tbl$plasmid_id[full_tbl$plasmid == "16337"] = "S"
full_tbl$plasmid_id[full_tbl$plasmid == "31815"] = "XS"
full_tbl$plasmid_id[full_tbl$plasmid == "49792"] = "L"
full_tbl$plasmid_id[full_tbl$plasmid == "52961"] = "XL"
full_tbl$plasmid_id = factor(full_tbl$plasmid_id, levels = c("XS","S", "M", "L","XL"))

g = ggplot(full_tbl,aes(contam_set_rate,cmedian,col=plasmid_id,group=plasmid_id)) + 
  geom_point(size=2) + 
  geom_line(size=1) + 
  theme_classic() + 
  geom_abline(slope = 1, intercept = 6) + 
  xlim(c(0,100)) + 
  ylim(c(0,100)) + 
  xlab("Simulated rate") + 
  ylab("Contamination estimate") +
  scale_colour_brewer(palette = "Dark2")

ggsave(g,file="/analysis/2021_08_26_PlasmidSeq_paper/figures/figure_4/contamination_perdicted_addgene_2021_11_07.pdf",width=5,height=5)


full_tbl = tbl %>% 
  group_by(plasmid,contam_set_rate) %>%
  mutate(cmedian=median(contamination_orig, na.rm=TRUE),mad=mad(contamination_orig, na.rm=TRUE))


full_tbl$plasmid = as.factor(full_tbl$plasmid)
sizes = c("12251" = 8890,"16337"=7397,"31815"=4290,"49792"=12598,"52961"=14873)
names = c("12251" = "medium","16337"="small","S"="XS","49792"="L","52961"="XL")
full_tbl$plasmid_id = "UNKNOWN"
full_tbl$plasmid_id[full_tbl$plasmid == "12251"] = "M"
full_tbl$plasmid_id[full_tbl$plasmid == "16337"] = "S"
full_tbl$plasmid_id[full_tbl$plasmid == "31815"] = "XS"
full_tbl$plasmid_id[full_tbl$plasmid == "49792"] = "L"
full_tbl$plasmid_id[full_tbl$plasmid == "52961"] = "XL"
full_tbl$plasmid_id = factor(full_tbl$plasmid_id, levels = c("XS","S", "M", "L","XL"))

g = ggplot(full_tbl,aes(contam_set_rate,cmedian,col=plasmid_id,group=plasmid_id)) + 
  geom_point(size=2) + 
  geom_line(size=1) + 
  theme_classic() + 
  geom_abline(slope = 1, intercept = 6) + 
  xlim(c(0,100)) + 
  ylim(c(0,100)) + 
  xlab("Simulated rate") + 
  ylab("Contamination estimate") +
  scale_colour_brewer(palette = "Dark2")

ggsave(g,file="/analysis/2021_08_26_PlasmidSeq_paper/figures/figure_4/contamination_perdicted_addgene_known_map_2021_11_07.pdf",width=5,height=5)
ggsave(g,file="/analysis/2021_08_26_PlasmidSeq_paper/figures/figure_4/contamination_perdicted_addgene_known_map_2021_11_07.png",width=5,height=5)



tbl$complete = (1.0 - abs(1.0 - tbl$contiguity)) * (1.0 - abs(1.0 - tbl$identity))
completeness = tbl %>% 
  group_by(plasmid,contam_set_rate) %>%
  summarise(completeness_median=median(complete, na.rm=TRUE),mad=mad(complete, na.rm=TRUE))
completeness$plasmid_id = "UNKNOWN"
completeness$plasmid_id[completeness$plasmid == "12251"] = "M"
completeness$plasmid_id[completeness$plasmid == "16337"] = "S"
completeness$plasmid_id[completeness$plasmid == "31815"] = "XS"
completeness$plasmid_id[completeness$plasmid == "49792"] = "L"
completeness$plasmid_id[completeness$plasmid == "52961"] = "XL"
completeness$plasmid_id = factor(completeness$plasmid_id, levels = c("XS","S", "M", "L","XL"))


g = ggplot(completeness) + 
  geom_bar(aes(x=plasmid_id,y=completeness_median,fill=plasmid_id),stat='identity') + 
  facet_wrap(. ~ completeness$contam_set_rate,ncol=12) + 
  geom_errorbar(aes(x=plasmid_id,ymin=completeness_median-mad, ymax=completeness_median+mad), width=.2,
                position=position_dodge(.9)) +
  scale_fill_brewer(palette = "Dark2") + 
  theme_classic() + 
  guides(fill=guide_legend(title="Plasmid")) +
  xlab("Plasmid") + 
  ylab("Completeness") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.border  = element_blank())

ggsave(g,file="/analysis/2021_08_26_PlasmidSeq_paper/figures/figure_4/contamination_assembly_2021_11_07.pdf",width=15,height=8)
ggsave(g,file="/analysis/2021_08_26_PlasmidSeq_paper/figures/figure_4/contamination_assembly_2021_11_07.png",width=15,height=8)
