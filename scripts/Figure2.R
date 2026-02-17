library("readxl")
library(ggplot2)
library(dplyr)
library(geosphere)
library(ggnewscale)
library(rnaturalearth)
library(stringr)
library(sf)
library(geodata)
library(raster)
library(RColorBrewer)
library(tidyr)
library(lubridate)

library(ggtree)
library(tidyverse)
library(tidytree)
library(ape)
library(treeio)
library(cowplot)


####Figure2

### Visualizing trees

#Fig2A (BA.3.x in a global subsampled tree)
timetree <- read.nexus("../data/global/global_and_ba3x/treetime_output/timetree.nexus") # ("data/timetree.nwk")
timetree$tip.label <- gsub("^['\"]|['\"]$", "", timetree$tip.label)
#timetree <- groupClade(timetree,.node=c(1359,1361))
metadata <- read_tsv("../data/global/global_and_ba3x/metadata/combined_all_metadata.tsv") 
metadata <- metadata %>% rename(label = strain) #rename sequence name column to 'label'
mrsd_val <- max(ymd(metadata$date), na.rm = TRUE)

p <- ggtree(timetree, mrsd = mrsd_val, as.Date = TRUE, size = 0.3) %<+% metadata +
  theme_tree2() +
  geom_tippoint(
    aes(fill = case_when(
      grepl("^BA\\.3\\.2", Nextclade_pango) ~ "BA.3.2",
      grepl("^BA\\.3",     Nextclade_pango) ~ "BA.3",
      TRUE                                   ~ "Other"
    )),
    shape = 21, stroke = 0.04, size = 1.8, colour = "black"
  ) +
  scale_fill_manual(
    values = c("BA.3.2" = "#59A14F", "BA.3"   = "#4E79A7", "Other"  = "grey70"
    ),
    name = NULL
  ) +
  scale_x_date(expand = expansion(mult = c(0.02, 0.10))) +
  coord_cartesian(clip = "off") +
  theme(
    legend.position   = c(0.15, 0.9),
    legend.background = element_blank(),
    plot.margin       = margin(5.5, 20, 5.5, 5.5)
  )

p

ggsave("../outputs/figure2/global_tree.png", p, width = 4, height = 6, dpi = 300, units = "in")
ggsave("../outputs/figure2/global_tree.pdf", p, width = 4, height = 6, units = "in")


#Fig2B (Root-to-tip analysis)
rtt <- read_csv("../data/global/global_and_ba3x/treetime_output/clock_results/rtt_clean.csv") %>% 
  filter(!grepl("^NODE", name))
metadata <- read_tsv("../data/global/global_and_ba3x/metadata/combined_all_metadata.tsv")
rtt_metadata <- rtt %>% rename(strain = name) %>% left_join(metadata, by = "strain")

p_rtt <- ggplot(rtt_metadata, aes(date.x, `root-to-tip distance`)) +
  geom_smooth(method="lm", se=FALSE, colour="grey60", linewidth=0.5, fullrange=TRUE) +
  geom_point(data=~filter(.x,!grepl("^BA\\.3",Nextclade_pango)),
             shape=1,size=1,stroke=0.5,colour="grey70") +
  geom_point(data=~filter(.x,grepl("^BA\\.3(\\.|$)",Nextclade_pango)&!grepl("^BA\\.3\\.2",Nextclade_pango)),
             shape=1,size=1,stroke=0.5,colour="#4E79A7") +
  geom_point(data=~filter(.x,grepl("^BA\\.3\\.2",Nextclade_pango)),
             shape=1,size=1,stroke=0.5,colour="#59A14F") +
  geom_smooth(data=~filter(.x,grepl("^BA\\.3(\\.|$)",Nextclade_pango)&!grepl("^BA\\.3\\.2",Nextclade_pango)),
              method="lm",se=FALSE,colour="#4E79A7",linewidth=0.5,fullrange=TRUE) +
  geom_smooth(data=~filter(.x,grepl("^BA\\.3\\.2",Nextclade_pango)),
              method="lm",se=FALSE,colour="#59A14F",linewidth=0.5,fullrange=TRUE) +
  theme_classic(base_size=10) +
  labs(x="Date",y="Root-to-tip divergence") +
  scale_x_continuous(limits=c(2020,2026),expand=expansion(mult=c(0.02,0.02)))

p_rtt

ggsave("../outputs/figure2/ba3_rtt.png", p_rtt, width = 4, height = 3.0, dpi = 300, units = "in")
ggsave("../outputs/figure2/ba3_rtt.pdf", p_rtt, width = 4, height = 3.0, units = "in")


# Slope/regression calculations
rtt_metadata <- rtt_metadata %>%
  mutate(lineage_group = case_when(grepl("^BA\\.3\\.2", Nextclade_pango) ~ "BA.3.2",
    grepl("^BA\\.3(\\.|$)", Nextclade_pango) ~ "BA.3", TRUE ~ "Other"))

library(broom)

rate_estimates <- rtt_metadata %>%
  filter(lineage_group %in% c("BA.3", "BA.3.2")) %>%
  group_by(lineage_group) %>%
  do(tidy(lm(`root-to-tip distance` ~ date.x, data = .))) %>%
  filter(term == "date.x") %>%
  select(lineage_group, estimate, std.error, statistic, p.value)

rate_estimates


#Fig2C (BA.3.x focused tree)
timetree_ba3 <- read.nexus("../data/global/global_and_ba3x/treetime_output/ba3_extracted_tree/ba3_extracted_tree.nexus") # ("data/timetree.nwk")
timetree_ba3$tip.label <- gsub("^['\"]|['\"]$", "", timetree_ba3$tip.label)
#timetree <- groupClade(timetree,.node=c(1359,1361))
metadata <- read_tsv("../data/global/global_and_ba3x/metadata/combined_all_metadata.tsv")
metadata <- read_tsv("../data/global/global_and_ba3x/metadata/gisaid/location_metadata.tsv") #Regional analysis
metadata <- metadata %>% rename(label = strain) #rename sequence name column to 'label'
metadata <- metadata %>% mutate(region = ifelse(region == "Africa", "Southern Africa", region))
metadata <- metadata %>% filter(!is.na(region) & region != "") # Remove N/A's
#TEMPORARY FIX (Add missing metadata manually)
missing_meta <- tibble(
  label  = c("hCoV-19/Australia/WA17182/2025",
             "hCoV-19/Mozambique/INS-PMB0734918/2025"),
  region = c("Oceania", "Southern Africa")
)
metadata <- bind_rows(metadata, missing_meta)

mrsd_val_ba3 <- max(ymd(metadata$date), na.rm = TRUE)

p_ba3 <- ggtree(timetree_ba3, mrsd = mrsd_val, as.Date = TRUE, size = 0.3) %<+% metadata +
  theme_tree2() +
  geom_tippoint(
    aes(fill = region),
    shape = 21, stroke = 0.04, size = 1.5, colour = "black"
  ) +
  scale_fill_brewer(palette = "Set3", name = "Region", direction = -1) +
  scale_x_date(expand = expansion(mult = c(0.02, 0.10))) +
  coord_cartesian(clip = "off") +
  theme(
    legend.position   = c(0.6, 0.7),
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 8),
    legend.background = element_blank(),
    plot.margin       = margin(5.5, 20, 5.5, 5.5)
  )

p_ba3

ggsave("../outputs/figure2/ba3_tree_wide.png", p_ba3, width = 4, height = 3, dpi = 300, units = "in")
ggsave("../outputs/figure2/ba3_tree.pdf", p_ba3, width = 4, height = 3, units = "in")