#!/usr/bin/env python
#-*- coding: utf-8 -*-
# Florian THIBORD (22/08/17)
########################################################
#                  OptimiR PIPELINE                    #
########################################################

# Standard libraries
import sys, os, subprocess

# Personal libraries
from . import filter_reads as CLEAN_SAM
from .essentials import *

def mapping(tmpdir_mapping, fastq_in, SAMPLE_NAME, BOWTIE2, SEEDLEN, ref_library, OptimiR_path, SAMTOOLS):
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
    CLEAN_SAM.main(file_out_mapping, fastq_in, 0)
    # Use samtools to convert sam to bam, and gunzip failed fq
    subprocess.call("{}/src/sam_to_bam.sh {} {} {}".format(OptimiR_path, tmpdir_mapping, SAMPLE_NAME, SAMTOOLS), shell=True)
    if ret_code != 0:
        raise Except_OS_Pipe("Merging step: samtools failed. Make sure it's installed and its path properly configured in PImiR_config.py")
