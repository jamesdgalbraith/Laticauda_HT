library(tidyverse)
library(plyranges)
library(BSgenome)

# set repeat name
repeat_name <- "all_Repbase_Harbingers"

# Read in repeat sequence and remove whotespace from names
repeat_seq <- readDNAStringSet(paste0("Laticauda_colubrina_HT/PIF-Harbinger/trees/", repeat_name, ".fasta"))
names(repeat_seq) <- sub(" .*", "", names(repeat_seq))
names(repeat_seq) <- sub("\t.*", "", names(repeat_seq))

# search for coding regions
rpst_out <- read_tsv(system(paste0("rpstblastn -query Laticauda_colubrina_HT/PIF-Harbinger/trees/", repeat_name, ".fasta -db /media/james/Samsung_T5/Databases/pfam/Pfam -outfmt \"6 qseqid sseqid qstart qend qlen sstart send slen evalue stitle\" -num_threads 16"), intern = T),
                     col_names = c("qseqid", "sseqid", "qstart", "qend", "qlen", "sstart", "send", "slen", "evalue", "stitle"))

# rearrange output
rpst_out_adjusted <- rpst_out %>%
  dplyr::filter(evalue < 0.01, send - sstart + 1 >= 0.7 * slen) %>%
  dplyr::mutate(strand = case_when(qend > qstart ~ "+", qend < qstart ~ "-"),
                seqnames = qseqid) %>%
  tidyr::separate(stitle, into = c("code", "name", "description"), sep = ", ") %>%
  dplyr::select(-description)

# Identify transposases
Tnp <- rpst_out_adjusted %>%
  filter(name %in% c("DDE_Tnp_4", "Plant_tran"), qstart < qend) %>%
  dplyr:: select(seqnames, qstart, qend, strand) %>%
  dplyr::rename(start = qstart, end = qend) %>%
  as_granges() %>%
  reduce_ranges_directed()

# identify DNA binding domains
MADF <- rpst_out_adjusted %>%
  filter(name %in% c("MADF", "MADF_DNA_bdg"), qstart > qend) %>%
  dplyr:: select(seqnames, qstart, qend, strand) %>%
  dplyr::rename(start = qend, end = qstart) %>%
  as_granges() %>%
  reduce_ranges_directed()

# Write Tnp sequences to file
tnp_seq <- Biostrings::getSeq(repeat_seq, Tnp)
names(tnp_seq) <- paste0(seqnames(Tnp),  "#", Tnp$species)
writeXStringSet(tnp_seq, paste0("Laticauda_colubrina_HT/PIF-Harbinger/trees/", repeat_name, "_tnp_nt.fasta"))
tnp_seq <- Biostrings::translate(tnp_seq, if.fuzzy.codon="solve")
writeXStringSet(tnp_seq, paste0("Laticauda_colubrina_HT/PIF-Harbinger/trees/", repeat_name, "_tnp_aa.fasta"))

# Write DNA binding sequences to file
madf_seq <- Biostrings::getSeq(repeat_seq, MADF)
names(madf_seq) <- seqnames(MADF)
writeXStringSet(madf_seq, paste0("Laticauda_colubrina_HT/PIF-Harbinger/trees/", repeat_name, "_madf_nt.fasta"))
madf_seq <- Biostrings::translate(madf_seq, if.fuzzy.codon="solve")
writeXStringSet(madf_seq, paste0("Laticauda_colubrina_HT/PIF-Harbinger/trees/", repeat_name, "_madf_aa.fasta"))
