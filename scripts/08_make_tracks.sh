#!/bin/bash

filename="$1"
echo "File with exp.names: {$1}"

# Annotation tracks
echo "Getting MANE_all tracks"
cat $filename | while read p; do bedtools intersect -s -wo -a $p.bed -b /home/projects/ulitsky/igoru/data/human/MANE.GRCh38.v1.0.refseq_genomic.bed > $p.MANE_all.txt; done

echo "Getting MANE_exonic tracks"
cat $filename | while read p; do bedtools intersect -s -wo -a $p.bed -b /home/projects/ulitsky/igoru/data/human/MANE.GRCh38.v1.0.refseq_genomic.bed -split > $p.MANE_exonic.txt; done

echo "Getting MANE_exonic_antisense tracks"
cat $filename | while read p; do bedtools intersect -S -wo -a $p.bed -b /home/projects/ulitsky/igoru/data/human/MANE.GRCh38.v1.0.refseq_genomic.bed -split > $p.MANE_exonic_antisense.txt; done

echo "Getting MANE_all_antisense tracks"
cat $filename | while read p; do bedtools intersect -S -wo -a $p.bed -b /home/projects/ulitsky/igoru/data/human/MANE.GRCh38.v1.0.refseq_genomic.bed > $p.MANE_all_antisense.txt; done

echo "Getting MANE_introns_antisense tracks"
cat $filename | while read p; do bedtools intersect -S -wo -a $p.bed -b /home/projects/ulitsky/lidiag/genomic_data/GCF_000001405.40_GRCh38.p14_MANE_introns_new_format.bed > $p.MANE_introns_antisense.txt; done

echo "Getting MANE_introns_sense tracks"
cat $filename | while read p; do bedtools intersect -s -wo -a $p.bed -b /home/projects/ulitsky/lidiag/genomic_data/GCF_000001405.40_GRCh38.p14_MANE_introns_new_format.bed > $p.MANE_introns_sense.txt; done

echo "Getting gencode_noncoding_all tracks"
cat $filename | while read p; do bedtools intersect -s -wo -a $p.bed -b /home/projects/ulitsky/lidiag/genomic_data/gencode.v46.annotation.noncoding.bed > $p.gencode_noncoding_all.txt; done

# Repeats track 
echo "Getting repeats tracks"
cat $filename | while read p; do bedtools intersect -s -wo -a $p.bed -b /home/projects/ulitsky/lidiag/genomic_data/hg38.repeats.bed > $p.repeats.txt; done
