#!/bin/bash
ref=$1
reads=$2
bam=${reads//\.fastq/}.on.${ref//\.fasta}
bwa index $ref
bwa mem -t 8 $ref $reads | samtools view -bS -F 4 - | samtools sort - $bam
samtools index $bam.bam
samtools mpileup -uf $ref $bam.bam | bcftools view -cg - | vcfutils.pl vcf2fq | seqtk seq -A -l 70 - | sed -r "s/\>.*/\>$bam.consensus/g" >> $bam.consensus.fasta