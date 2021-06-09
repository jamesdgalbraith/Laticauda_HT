#!/bin/bash

# Convert Repeatmasker out to GFF2
tail -n +4 latLat_1.0.fasta.out | grep -v "\\*" | awk '{OFS="\t"}{if ($9=="C") {print $5,"RepeatMasker",$11,$7,$6,".","-",".","Target="$10";Divergence="$2}; if ($9=="+") {print $5,"RepeatMasker",$11,$6,$7,".","+",".","Target="$10";Divergence="$2}}' > latLat_1.0.fasta.gff

# Extract Harbinger-Snek sequences from Repeatmasker annotation
grep Harbinger latLat_1.0.fasta.gff | grep Laticauda > Harbinger-Snek_latLat_1.0.fasta.gff

# Interesect Repeat and Gene annotations and select exons
bedtools intersect -wa -a latLat_1.0_liftoff.gff -b Harbinger-Snek_latLat_1.0.fasta.gff | grep exon >  latLat_1.0_intersect.fasta.gff
