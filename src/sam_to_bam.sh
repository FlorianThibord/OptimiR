#!/bin/bash
#-*- coding: utf-8 -*-
#########################################################################
### OPTIMIR pipeline: merge failed fastq, create sorted and indexed bam
### Florian THIBORD, 04/07
#########################################################################

DIR_MAPP=$1
NAME=$2
SAMTOOLS=$3

# $NAME.al.sam : output file from bowtie2
# $NAME.ok.sam : output file from filter_reads.py
# $NAME.failed.fq : reads which failed alignment with bowtie2 
# $NAME.cl.fq : output from filter_reads.py (reads discarded by filter_reads.py)

cat $DIR_MAPP/$NAME.cl.fq >> $DIR_MAPP/$NAME.failed.fq
rm $DIR_MAPP/$NAME.cl.fq

mv $DIR_MAPP/$NAME.ok.sam $DIR_MAPP/$NAME.sam
$SAMTOOLS view -bh $DIR_MAPP/$NAME.sam -o $DIR_MAPP/$NAME.bam
rm -f $DIR_MAPP/$NAME.al.sam ## aligned by bowtie2
rm -f $DIR_MAPP/$NAME.sam ## cleaned by filter_reads

# create bam
$SAMTOOLS sort $DIR_MAPP/$NAME.bam -o $DIR_MAPP/$NAME.s.bam
rm -f $DIR_MAPP/$NAME.bam
mv $DIR_MAPP/$NAME.s.bam $DIR_MAPP/$NAME.bam
$SAMTOOLS index $DIR_MAPP/$NAME.bam

# Gunzip failed fq
gzip -f $DIR_MAPP/$NAME.failed.fq
