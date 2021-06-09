library(tidyverse)
library(reshape2)
library(RColorBrewer)

# Read in divergence data
lat_col_rm_divsum <- readr::read_delim(system("tail -n 72 RepeatMasker/Laticauda_colubrina/latCor_1.0.fasta.divsum", intern = T),
                  delim = " ", col_names = T)
lat_col_rm_divsum <- lat_col_rm_divsum[1:ncol(lat_col_rm_divsum)-1]
lat_lat_rm_divsum <- readr::read_delim(system("tail -n 72 RepeatMasker/Laticauda_laticaudata/latLat_1.0.fasta.divsum", intern = T),
                                       delim = " ", col_names = T)
lat_lat_rm_divsum <- lat_lat_rm_divsum[1:ncol(lat_lat_rm_divsum)-1]
not_scu_rm_divsum <- readr::read_delim(system("tail -n 72 RepeatMasker/Notechis_scutatus/TS10Xv2-PRI.fasta.divsum", intern = T),
                                       delim = " ", col_names = T)
not_scu_rm_divsum <- not_scu_rm_divsum[1:ncol(not_scu_rm_divsum)-1]

# Select DNA transposons
lat_col_rm_TcMar <- lat_col_rm_divsum %>%
  dplyr::select(starts_with("DNA/Tc"))
lat_col_rm_hAT <- lat_col_rm_divsum %>%
  dplyr::select(starts_with("DNA/hAT"))
lat_col_rm_PIF <- lat_col_rm_divsum %>%
  dplyr::select(starts_with("DNA/PIF"))
lat_col_rm_DNA_other <- lat_col_rm_divsum %>%
  dplyr::select(starts_with("DNA")) %>%
  dplyr::select(!starts_with("DNA/Tc")) %>%
  dplyr::select(!starts_with("DNA/hAT")) %>%
  dplyr::select(!starts_with("DNA/PIF"))

# Calculate percentages
lat_col_all_DNA <- tibble(div = 0:70,
                          `DNA/PIF-Harbinger` = rowSums(lat_col_rm_PIF),
                          `DNA/hAT` = rowSums(lat_col_rm_hAT),
                          `DNA/TcMar` = rowSums(lat_col_rm_TcMar),
                          `DNA/Other` = rowSums(lat_col_rm_DNA_other)) %>%
  reshape2::melt(id = "div") %>%
  mutate(value = value/20246879.24)

# Select DNA transposons
lat_lat_rm_TcMar <- lat_lat_rm_divsum %>%
  dplyr::select(starts_with("DNA/Tc"))
lat_lat_rm_hAT <- lat_lat_rm_divsum %>%
  dplyr::select(starts_with("DNA/hAT"))
lat_lat_rm_PIF <- lat_lat_rm_divsum %>%
  dplyr::select(starts_with("DNA/PIF"))
lat_lat_rm_DNA_other <- lat_lat_rm_divsum %>%
  dplyr::select(starts_with("DNA")) %>%
  dplyr::select(!starts_with("DNA/Tc")) %>%
  dplyr::select(!starts_with("DNA/hAT")) %>%
  dplyr::select(!starts_with("DNA/PIF"))

# Calculate percentages
lat_lat_all_DNA <- tibble(div = 0:70,
                          `DNA/PIF-Harbinger` = rowSums(lat_lat_rm_PIF),
                          `DNA/hAT` = rowSums(lat_lat_rm_hAT),
                          `DNA/TcMar` = rowSums(lat_lat_rm_TcMar),
                          `DNA/Other` = rowSums(lat_lat_rm_DNA_other)) %>%
  reshape2::melt(id = "div") %>%
  mutate(value = value/15587061.06)

# Select DNA transposons
not_scu_rm_TcMar <- not_scu_rm_divsum %>%
  dplyr::select(starts_with("DNA/Tc"))
not_scu_rm_hAT <- not_scu_rm_divsum %>%
  dplyr::select(starts_with("DNA/hAT"))
not_scu_rm_PIF <- not_scu_rm_divsum %>%
  dplyr::select(starts_with("DNA/PIF"))
not_scu_rm_DNA_other <- not_scu_rm_divsum %>%
  dplyr::select(starts_with("DNA")) %>%
  dplyr::select(!starts_with("DNA/Tc")) %>%
  dplyr::select(!starts_with("DNA/hAT")) %>%
  dplyr::select(!starts_with("DNA/PIF"))

# Calculate percentages
not_scu_all_DNA <- tibble(div = 0:70,
                          `DNA/PIF-Harbinger` = rowSums(not_scu_rm_PIF),
                          `DNA/hAT` = rowSums(not_scu_rm_hAT),
                          `DNA/TcMar` = rowSums(not_scu_rm_TcMar),
                          `DNA/Other` = rowSums(not_scu_rm_DNA_other)) %>%
  reshape2::melt(id = "div") %>%
  mutate(value = value/16655259.58)

# plot plot plot!!!
ggplot(lat_col_all_DNA, aes(x = div, y = value, fill = variable)) + geom_bar(stat = "identity") +
  scale_x_continuous(limits = c(-1, 40), expand = c(0,0), name = "Kimura substituion lebel (CpG adjusted)", breaks = (0:4*10)) +
  scale_y_continuous(limits = c(0, 3.5), name = "Perecentage of genome", expand = c(0,0)) + theme_bw() +
  ggtitle(label = "Laticauda colubrina") + theme(plot.title = element_text(face="italic", hjust = 0.5)) +
  scale_fill_brewer(name = "Transposon family", palette = "YlGnBu", direction = -1)
ggsave("Laticauda_colubrina_HT/PIF-Harbinger/divergence/latCol_plot.svg")

ggplot(lat_lat_all_DNA, aes(x = div, y = value, fill = variable)) + geom_bar(stat = "identity") +
  scale_x_continuous(limits = c(-1, 40), expand = c(0,0), name = "Kimura substituion lebel (CpG adjusted)", breaks = (0:4*10)) +
  scale_y_continuous(limits = c(0, 3.5), name = "Perecentage of genome", expand = c(0,0)) + theme_bw() +
  ggtitle(label = "Laticauda laticaudata") + theme(plot.title = element_text(face="italic", hjust = 0.5)) +
  scale_fill_brewer(name = "Transposon family", palette = "YlGnBu", direction = -1)
ggsave("Laticauda_colubrina_HT/PIF-Harbinger/divergence/latLat_plot.svg")

ggplot(not_scu_all_DNA, aes(x = div, y = value, fill = variable)) + geom_bar(stat = "identity") +
  scale_x_continuous(limits = c(-1, 40), expand = c(0,0), name = "Kimura substituion lebel (CpG adjusted)", breaks = (0:4*10)) +
  scale_y_continuous(limits = c(0, 3.5), name = "Perecentage of genome", expand = c(0,0)) + theme_bw() +
  ggtitle(label = "Notechis scutatus") + theme(plot.title = element_text(face="italic", hjust = 0.5)) +
  scale_fill_brewer(name = "Transposon family", palette = "YlGnBu", direction = -1)
ggsave("Laticauda_colubrina_HT/PIF-Harbinger/divergence/notScu_plot.svg")
