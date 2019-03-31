#!/usr/bin/env python
# -*- coding: utf-8 -*-
########################################################
#                  OptimiR PIPELINE                    #
########################################################

# Standard libraries
import sys, os, time
import subprocess
from pysam import AlignmentFile

# Personal libraries
from optimir.libs.essentials import *
from optimir.libs.library_preparation import prepare_library
from optimir.libs.pre_process import trimming, collapsing
from optimir.libs.mapping import mapping
from optimir.libs.post_process import post_process_main

def process(args):
    ## Assign arguments to variables (all caps for variables corresponding to args)
    FASTQ = args.FASTQ
    VCF = args.VCF
    OUTDIR = args.OUTDIR
    SEEDLEN = args.SEEDLEN
    WEIGHT5 = args.WEIGHT5
    SCORE_THRESHOLD = args.SCORE_THRESH
    INCONSISTENT_THRESHOLD = args.INCONSISTENT_THRESHOLD
    RMTEMP = args.RMTEMP
    ANNOT_FILES = args.ANNOT_FILES
    WRITE_GFF = args.WRITE_GFF
    WRITE_VCF = args.WRITE_VCF
    ADAPT3 = args.ADAPT3
    ADAPT5 = args.ADAPT5    
    READMIN = args.READMIN
    READMAX = args.READMAX
    BQTHRESH = args.BQTHRESH
    TRIM_AGAIN = args.TRIM_AGAIN
    MATURES = args.MATURES
    HAIRPINS = args.HAIRPINS
    GFF3 = args.GFF3
    VERBOSE = args.VERBOSE
    CUTADAPT = args.CUTADAPT
    BOWTIE2 = args.BOWTIE2
    BOWTIE2_BUILD = args.BOWTIE2_BUILD
    SAMTOOLS = args.SAMTOOLS

    try:
        if VCF == None:
            VCF_AVAIL = False
            if WRITE_VCF:
                print("WARNING: no variant provided, thus no VCF output will be generated\n")
                WRITE_VCF = False
        else:
            check_if_file_exists(VCF)
        check_if_file_exists(MATURES)
        check_if_file_exists(HAIRPINS)
        check_if_file_exists(GFF3)
        check_if_file_exists(FASTQ)
        
    except InputError as err:
        print("ERROR during library preparation: file {} does not exists. Check input filename and try again.\n".format(err.input_name))
        sys.exit(4)
    
    if OUTDIR[-1] == "/":
        OUTDIR = OUTDIR[:-1]
    tmpdir = os.path.abspath(OUTDIR) + "/OptimiR_tmp"
    
    ########################
    ##    MAIN PROCESS    ##
    ########################

    try:
        start = time.time()
        subprocess.call("mkdir -p " + OUTDIR, shell=True)
        subprocess.call("mkdir -p " + tmpdir, shell=True)
        sample_name = os.path.basename(FASTQ)
        sample_name = sample_name.split('.')[0] ## remove extension
        fun_str_progress([sample_name], "header", VERBOSE) ## print header
        ##############################
        # LIBRARY PREPARATION
        # Create directory where OptimiR files will be stored
        out_directory = os.path.abspath(OUTDIR) + "/OptimiR_lib/"
        # Create directory for OptimiR fasta library
        fasta_dir =  out_directory + "fasta/"
        subprocess.call("mkdir -p " + fasta_dir, shell=True)
        fasta_file = fasta_dir + "optimiR_library.fa"
        # Create directory and path for bowtie2 index alignment
        index_dir =  out_directory + "bowtie2_index/"
        subprocess.call("mkdir -p " + index_dir, shell=True)
        index_path = index_dir + "optimiR_alignment_library"
        # Create directory for pickle objects
        pickle_dir = out_directory + "pickle/"
        subprocess.call("mkdir -p " + pickle_dir, shell=True)
        lib_infos_pickle_path = pickle_dir + "lib_infos.pkl"
        d_OptimiR_pickle_path = pickle_dir + "d_OptimiR.pkl"
        # in this directory, the new library is created
        prepare_library(BOWTIE2_BUILD, VCF, MATURES, HAIRPINS, GFF3, out_directory, fasta_file, index_path, lib_infos_pickle_path, d_OptimiR_pickle_path)
        VCF_AVAILABLE, GENO_AVAILABLE, sample_list, fasta_hash, vcf_hash = load_obj(lib_infos_pickle_path)
        if GENO_AVAILABLE and not(sample_name in sample_list):
            print("WARNING : Sample provided {} does not match any genotyped sample in provided vcf file".format(sample_name))
            GENO_AVAILABLE = False
        lp_time = time.time()
        fun_str_progress([VCF_AVAILABLE, GENO_AVAILABLE, round(lp_time - start, 2)], "lib_prep", VERBOSE)
        ##############################
        # PRE ALIGNMENT PROCESS : Adapter trimming / Size selection + Read collapsing
        tmpdir_trim = tmpdir + '/0_Trimming'
        if not(os.path.exists(tmpdir_trim + "/" + sample_name + ".trimmed.fq")) or TRIM_AGAIN:
            subprocess.call('mkdir -p {}'.format(tmpdir_trim), shell=True)
            trimming(FASTQ, sample_name, tmpdir_trim, CUTADAPT, ADAPT3, ADAPT5, READMIN, READMAX, BQTHRESH)
        trimming_time = time.time()
        fun_str_progress([round(trimming_time - lp_time, 2)], "trim", VERBOSE)
        tmpdir_collapsed = tmpdir + '/1_Collapsing'
        subprocess.call('mkdir -p {}'.format(tmpdir_collapsed), shell=True)
        collapse_table, collapse_report = collapsing(sample_name, tmpdir_collapsed, tmpdir_trim)
        collapsing_time = time.time()
        fun_str_progress([round(collapsing_time - trimming_time, 2)], "collaps", VERBOSE)
        ##############################
        # ALIGNMENT : local mode on index_path library
        tmpdir_mapping = tmpdir + '/2_Mapping'
        subprocess.call('mkdir -p {}'.format(tmpdir_mapping), shell=True)
        fastq_in = '{}/{}.clpsd.fq'.format(tmpdir_collapsed, sample_name)
        mapping(tmpdir_mapping, fastq_in, sample_name, BOWTIE2, SEEDLEN, index_path, SAMTOOLS)
        mapping_time = time.time()
        fun_str_progress([round(mapping_time - collapsing_time, 2)], "mapping", VERBOSE)
        ##############################
        # POST ALIGNMENT PROCESS : Read Annotation + Genotype Consistance + Alignment Scoring + Abundances generation
        tmpdir_postProcess = tmpdir + "/3_PostProcess"
        subprocess.call("mkdir -p {}".format(tmpdir_postProcess), shell=True)
        dir_results = OUTDIR + "/OptimiR_Results"
        subprocess.call("mkdir -p {}".format(dir_results), shell=True)
        sourceDB = MATURES.split('/')[-1].split('.')[0]
        post_process_main(tmpdir_mapping, tmpdir_postProcess, dir_results, collapse_table, sample_name, WEIGHT5, SCORE_THRESHOLD, INCONSISTENT_THRESHOLD, d_OptimiR_pickle_path, VCF_AVAILABLE, GENO_AVAILABLE, ANNOT_FILES, VERBOSE, WRITE_GFF, WRITE_VCF, sourceDB)
        pp_time = time.time()
        fun_str_progress([round(pp_time - mapping_time, 2)], "postproc", VERBOSE)
        end = time.time()
        if RMTEMP:
            subprocess.call("rm -r {}".format(tmpdir), shell=True)
        fun_str_progress([OUTDIR + "/OptimiR_Results/", round(end - start, 2)], "footer", VERBOSE)
    except Except_OS_Pipe as e:
        print('Error with system call: {}\nCheck that path to Bowtie2, Cutadapt, Samtools are well defined.\n'.format(e))
        exit(3)
