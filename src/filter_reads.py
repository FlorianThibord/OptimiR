#!/usr/bin/env python
#-*- coding: utf-8 -*-
# Florian THIBORD (22/08/17)
########################################################
#                  OptimiR - SamCleaner                #
########################################################
# Filter reads with more than NB_MM mismatches
# BUG in Bowtie2: allow indels in seed while no mismatch is allowed in local mode

import sys

from pysam import AlignmentFile
from Bio import SeqIO

#Retrieve seqs from fastq of failed alignments
#Output: mono mapping deleted alignments
def alignments_to_seqs(reads, fastq_fn):
    fastq_dict = SeqIO.index(fastq_fn, 'fastq') # MUCH FASTER!!!!
    seqs_to_keep = []
    for read in reads:
        seq = fastq_dict[read]
        seqs_to_keep.append(seq)
    return seqs_to_keep

def delete_alignment(alignment, NB_MM):
    nb_mismatchs = alignment.get_tag('NM')
    read_len = alignment.query_length
    delete = False
    if (nb_mismatchs > NB_MM):
        delete = True
    return delete

# Delete alignment if:
#  - mm in alignment
# If deleted alignment read is mono mapping, then put it back in failed alignments .fastq
def clean_sam(sam_fn, clean_sam_fn, orig_fq_fn, NB_MM):
    sam = AlignmentFile(sam_fn, 'r')
    clean_sam = AlignmentFile(clean_sam_fn, 'wb', template=sam)
    reads_kept = set()
    reads_deleted = set()
    for alignment in sam:
        if delete_alignment(alignment, NB_MM):
            reads_deleted.add(alignment.query_name)
        else:
            clean_sam.write(alignment)
            reads_kept.add(alignment.query_name)
    # only add alignment to fastq if read involved hasn't been kept on a cross-mapping loci
    keep_as_failed_list = []
    for read in reads_deleted:
        if not(read in reads_kept):
            keep_as_failed_list.append(read)
    seqs_failed = alignments_to_seqs(keep_as_failed_list, orig_fq_fn)
    ## Print discarded alignments in fastq 
    fastq_output_fn = sam_fn[:-7] + ".cl.fq"
    SeqIO.write(seqs_failed, fastq_output_fn, 'fastq')
    ## Print in log file
    log_clean = open("{}.SamCleaner.log".format(sam_fn[:-7]), 'a')
    cpt_treated = len(reads_kept) + len(reads_deleted)
    log_clean.write('Number of alignments treated : {}\nNumber of alignment kept : {}\nNumber of alignments deleted : {}\n'.format(cpt_treated, len(reads_kept), len(reads_deleted)))
    log_clean.close()
    clean_sam.close()

def main(sam_filename, orig_fq_fn, NB_MM):
    sam_fn = sam_filename
    clean_sam_fn = "{}.ok.sam".format(sam_fn[:-7])
    clean_sam(sam_fn, clean_sam_fn, orig_fq_fn, NB_MM)
