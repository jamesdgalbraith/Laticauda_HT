# read in libraries, hiding messages
suppressMessages(library(tidyverse))
suppressMessages(library(plyranges))

## importing data
# Gene annotation (remove unknown loci)
liftoff_gff <- plyranges::read_gff3(file = "Laticauda_colubrina_HT/liftoff/latLat_1/latLat_1.0_liftoff.gff") %>%
  dplyr::filter(!str_detect(gene, "^LOC"))
# Repeat annotation
rm_gff <- plyranges::read_gff3(file = "RepeatMasker/For_Laticauda_HT/Laticauda_laticaudata/latLat_1.0_repeatmasker.gff")
# Repeat class
rm_classification <- read_tsv(system(paste0("grep \">\" RepeatMasker/custom_library.fa"), intern = T), col_names = "Target") %>%
  tidyr::separate(Target, into = c("Target", "repeat_class"), sep = "#") %>%
  dplyr::mutate(Target = sub(">", "", Target))

# Give repeats their class
rm_gff <- rm_gff %>%
  as_tibble() %>%
  dplyr::select(-width) %>%
  tidyr::separate(col = Target, into = c("Target", "qstart", "qend"), sep = " ") %>%
  inner_join(rm_classification)

# Select Harbingers
harbinger_gff <- rm_gff %>%
  filter(repeat_class == "DNA/PIF-Harbinger", grepl("Laticauda", Target)) %>%
  as_granges() %>%
  GenomicRanges::reduce(min.gapwidth = 500L)

# Select genes from liftoff
liftoff_genes <- liftoff_gff %>%
  filter(type == "gene")

# Intersect Harbingers with annotation and rearrange
overlapped_gffs <- plyranges::pair_overlaps(liftoff_gff, harbinger_gff)
overlapped_gffs_tibble <- overlapped_gffs %>%
  as_tibble() %>%
  dplyr::select(granges.x.seqnames, granges.x.start, granges.x.end, granges.x.strand, granges.y.seqnames, granges.y.start, granges.y.end, granges.y.strand,
                source, type, ID,  Name, gbkey, gene, gene_biotype, coverage, sequence_ID, extra_copy_number, copy_num_ID, Parent) %>%
  dplyr::rename(x.seqnames = granges.x.seqnames, x.start = granges.x.start, x.end = granges.x.end, x.strand = granges.x.strand,
                y.seqnames = granges.y.seqnames, y.start = granges.y.start, y.end = granges.y.end, y.strand = granges.y.strand)

# Identify genes containing Harbinger insertions
overlapped_gffs_gene <- overlapped_gffs_tibble %>%
  filter(type == "gene")

# Make list of genes containing Harbinger insertions
overlapped_gffs_genes_tbl <- tibble(gene = names(table(overlapped_gffs_gene$gene)), n = as.integer(table(overlapped_gffs_gene$gene)))
overlapped_gffs_genes_tbl %>%
  dplyr::select(gene) %>%
  write_tsv("Laticauda_colubrina_HT/PIF-Harbinger/latLat_in_gene.txt", col_names = F)

# Identify genes containing Harbinger insertions
overlapped_gffs_exon <- overlapped_gffs_tibble %>%
  dplyr::filter(type == "exon")

# Make list of exons containing Harbinger insertions
overlapped_gffs_exon_tbl <- tibble(gene = names(table(overlapped_gffs_exon$gene)), n = as.integer(table(overlapped_gffs_exon$gene)))
overlapped_gffs_exon_tbl %>%
  dplyr::select(gene) %>%
  write_tsv("Laticauda_colubrina_HT/PIF-Harbinger/latLat_in_exons.txt", col_names = F)

# Identify insertions within 5000bp of 3'UTRs
upstream_of_genes <- plyranges::pair_overlaps(liftoff_genes, decent, maxgap = 5000) %>%
  as_tibble() %>%
  mutate(granges.x.strand = as.character(granges.x.strand), granges.y.strand = as.character(granges.y.strand)) %>%
  filter((granges.x.start > granges.y.end & granges.x.strand == "+") |
           (granges.x.end < granges.y.start & granges.x.strand == "-")) %>%
  dplyr::select(gene) %>%
  base::unique()
write_tsv(upstream_of_genes, "Laticauda_colubrina_HT/PIF-Harbinger/latLat_upstream_of_gene.txt", col_names = F)
