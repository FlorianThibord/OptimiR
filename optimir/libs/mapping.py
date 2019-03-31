#!/usr/bin/env python
#-*- coding: utf-8 -*-
# Florian THIBORD (22/08/17)
########################################################
#                  OptimiR PIPELINE                    #
########################################################

# Standard libraries
import sys, os, subprocess

# Personal libraries
import optimir.libs.filter_reads as filter_reads
from optimir.libs.essentials import *

def mapping(tmpdir_mapping, fastq_in, SAMPLE_NAME, BOWTIE2, SEEDLEN, ref_library, SAMTOOLS):
    file_out_mapping = '{}/{}.al.sam'.format(tmpdir_mapping, SAMPLE_NAME)
    # Remove potential previous files
    subprocess.call("rm -f {}/{}*".format(tmpdir_mapping, SAMPLE_NAME), shell=True)
    # Bowtie2 call : very-sensitive-local, no mismatches, no reverse complement, 1 core
    command_line = ("{} --local -a -D 20 -R 3 -N 0 -L {} -i S,1,0.5 --no-1mm-upfront --norc --mp 100 --score-min G,1,8 --no-unal -p 1 --ignore-quals".format(BOWTIE2, SEEDLEN) +
                    " -x {} -q {} -S {} --un {}/{}.failed.fq 2> {}/{}.mapping_report.txt".format(ref_library, fastq_in, file_out_mapping, tmpdir_mapping, SAMPLE_NAME, tmpdir_mapping, SAMPLE_NAME))
    ret_code = subprocess.call(command_line, shell=True)
    if ret_code != 0:
        raise Except_OS_Pipe('Alignment failed: command line = {}'.format(command_line))
    # Clean sam (clean reads with indels present in seed, bug bowtie2?)
    filter_reads.main(file_out_mapping, fastq_in, 0)
    # Use samtools to convert sam to bam, and gunzip failed fq
    sam_to_bam(tmpdir_mapping, SAMPLE_NAME, SAMTOOLS)
