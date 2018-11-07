#!/usr/bin/env python
#-*- coding: utf-8 -*-
# Florian THIBORD (17/11/17)
########################################################
#                  OptimiR PIPELINE                    #
########################################################
# 1) PRE ALIGNMENT PROCESS : Trimming & Collapsing

# Standard libraries
import sys, os
import subprocess

# Personal libraries
from . import read_collapser as Collapser
from .essentials import *

## Trim reads using cutadapt
def trimming(FASTQ, SAMPLE_NAME, tmpdir_trim, CUTADAPT, ADAPT3, ADAPT5, READMIN, READMAX, BQTHRESH):
        command_line = ("{} -a {} -g {} -m {} -M {} -q {},{} -e 0.1 --max-n 0".format(CUTADAPT,
                                                                                     ADAPT3,
                                                                                     ADAPT5,
                                                                                     READMIN,
                                                                                     READMAX,
                                                                                     BQTHRESH,
                                                                                     BQTHRESH) +
                        " {} -o {}/{}.trimmed.fq > {}/{}.trimming_report.txt".format(FASTQ,
                                                                                     tmpdir_trim,
                                                                                     SAMPLE_NAME,
                                                                                     tmpdir_trim,
                                                                                     SAMPLE_NAME))
        ret_code = subprocess.call(command_line, shell=True)
        if ret_code != 0:
            raise Except_OS_Pipe('Cutadapt failed : cmd line = {}'.format(command_line))
    
# Collapse identical reads
def collapsing(SAMPLE_NAME, tmpdir_collapsed, tmpdir_trim):
    collapsed_table, report = Collapser.collapse_sample(SAMPLE_NAME, tmpdir_trim, tmpdir_collapsed)
    return collapsed_table, report
