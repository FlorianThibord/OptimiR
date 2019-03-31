#!/usr/bin/env python
# -*- coding: utf-8 -*-
########################################################
#                  OptimiR PIPELINE                    #
########################################################

# Standard libraries
import os, sys
from argparse import ArgumentParser

# Personal libraries
from optimir.libs.process import process
from optimir.libs.library_preparation import library_preparation
from optimir.libs.summarize import summarize

### MAIN
def main():
    ## OptimiR path:
    optimiR_path = os.path.abspath(os.path.dirname(__file__))
    ## Parsing arguments
    parser = ArgumentParser(prog='optimir', description="OptimiR: A bioinformatics pipeline designed to detect and quantify miRNAs, isomiRs and polymiRs from miRSeq data, & study the impact of genetic variations on polymiRs' expression")

    subparsers = parser.add_subparsers(title='OptimiR subcommands')

    ##############
    ## Process
    ##############

    process_parser = subparsers.add_parser('process', help="Process a single fastq file with the OptimiR workflow")
    
    ## Mandatory arguments : fastq file, vcf file (optional but strongly advised), output directory
    process_parser.add_argument("-i", "--fq", dest="FASTQ", default=None, required=True, help="Full path of the sample fastq file (accepted formats and extensions: fastq, fq and fq.gz)")
    process_parser.add_argument("-o", "--dirOutput", dest="OUTDIR", default="./OptimiR_Results_Dir/", required=False, help="Full path of the directory where output files are generated [default: ./OptimiR_Results_Dir/]")
    process_parser.add_argument("-g", "--vcf", dest="VCF", default=None, required=False, help="Full path of the vcf file with genotypes")

    ## Alignment & Score Parameters : seedLen, weight5, scoreThresh
    process_parser.add_argument("--seedLen", dest="SEEDLEN", type=int, default=17, required=False, help="Choose the alignment seed length used in option '-L' by Bowtie2 [default: 17]")
    process_parser.add_argument("--w5", dest="WEIGHT5", type=int, default=4, required=False, help="Choose the weight applied on events detected on the 5' end of aligned reads [default: 4]")
    process_parser.add_argument("--scoreThresh", dest="SCORE_THRESH", type=int, default=9, required=False, help="Choose the threshold for alignment score above which alignments are discarded [default: 9]")
    process_parser.add_argument("--consistentRate", dest="INCONSISTENT_THRESHOLD", type=float, default=0.01, required=False, help="Choose the rate threshold for inconsistent reads mapped to a polymiR above which the alignment is flagged as highly suspicious [default: 0.01]")

    ## Optionnal arguments : rmTempFiles, annot, adapt 3', adapt 5', readMin, readMax, bqTresh, trimAgain
    process_parser.add_argument("--rmTempFiles", dest="RMTEMP", default=False, required=False, action='store_true', help="Add this option to remove temporary files (trimmed fastq, collapsed fastq, mapped reads, annotated bams")
    process_parser.add_argument("--annot", dest="ANNOT_FILES", default="hpics", required=False, help="Control which annotation file is produced by adding corresponding letter : 'h' for expressed_hairpins, 'p' for polymiRs_table, 'i' for consistency_table, 'c' for remaining_ambiguous, 's' for isomiRs_dist. Ex: '--annot hpics' [default] will produce all of them")
    process_parser.add_argument("--gff_out", dest="WRITE_GFF", default=False, required=False, action='store_true', help="Add this option to generate results in mirGFF3 format [default : disabled]")
    process_parser.add_argument("--vcf_out", dest="WRITE_VCF", default=False, required=False, action='store_true', help="Add this option to generate results in VCF format [default : disabled]")
    process_parser.add_argument("--adapt3", dest="ADAPT3", default="AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC -a TGGAATTCTCGGGTGCCAAGG", required=False, help="Define the 3' adaptor sequence (default is NEB & ILLUMINA: AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC -a TGGAATTCTCGGGTGCCAAGG -> hack: use -a to add adapter sequences)")
    process_parser.add_argument("--adapt5", dest="ADAPT5", default="", required=False, help="Define the 5' adaptor sequence [default: None]")
    process_parser.add_argument("--readMin", dest="READMIN", type=int, default=15, required=False, help="Define the minimum read length defined with option -m in cutadapt [default: 15]")
    process_parser.add_argument("--readMax", dest="READMAX", type=int, default=27, required=False, help="Define the maximum read length defined with option -M in cutadapt [default: 27]")
    process_parser.add_argument("--bqThresh", dest="BQTHRESH", type=int, default=28, required=False, help="Define the Base Quality threshold defined with option -q in cutadapt [default: 28]")
    process_parser.add_argument("--trimAgain", dest="TRIM_AGAIN", default=False, required=False, action='store_true', help="Add this option to trim files that have been trimmed in a previous application. By default, when temporary files are kept, trimmed files are reused. If you wish to change a paramater used in the trimming step of the workflow, this parameter is a must [default: disabled]")

    ## Optionnal arguments concerning the miRBase library used: matures.fa, hairpins.fa, miRs_coords.gff3. Usefull to switch from miRBase to miRCarta, or a new version of the miRBase
    process_parser.add_argument("--maturesFasta", dest="MATURES", default="{}/resources/fasta/hsa_matures_miRBase_v21.fa".format(optimiR_path), required=False, help="Path to the reference library containing mature miRNAs sequences [default: miRBase 21]")
    process_parser.add_argument("--hairpinsFasta", dest="HAIRPINS", default="{}/resources/fasta/hsa_hairpins_miRBase_v21.fa".format(optimiR_path), required=False, help="Path to the reference library containing pri-miRNAs sequences [default: miRBase 21]")
    process_parser.add_argument("--gff3", dest="GFF3", default="{}/resources/coordinates/hsa_miRBase_v21.gff3".format(optimiR_path), required=False, help="Path to the reference library containing miRNAs and pri-miRNAs coordinates [default: miRBase v21, GRCh38 coordinates]")

    process_parser.add_argument("--quiet", dest="VERBOSE", default=True, required=False, action='store_false', help="Add this option to remove OptimiR progression on screen [default: disabled]")

    ## Optionnal paths to cutadapt, bowtie2 and samtools (mandatory if not in $PATH)
    process_parser.add_argument("--cutadapt", dest="CUTADAPT", default="cutadapt", required=False, help="Provide path to the cutadapt binary [default: from $PATH]")
    process_parser.add_argument("--bowtie2", dest="BOWTIE2", default="bowtie2", required=False, help="Provide path to the bowtie2 binary [default: from $PATH]")
    process_parser.add_argument("--bowtie2_build", dest="BOWTIE2_BUILD", default="bowtie2-build", required=False, help="Provide path to the bowtie2 index builder binary [default: from $PATH]")
    process_parser.add_argument("--samtools", dest="SAMTOOLS", default="samtools", required=False, help="Provide path to the samtools binary [default: from $PATH]")
    process_parser.set_defaults(func=process)
    
    #######################
    ## Summarize
    #######################
    summarize_parser = subparsers.add_parser('summarize', help="Summarize results of previously processed samples")
    summarize_parser.add_argument("--dir", dest="DIR", default=None, required=True, help="Full path of the directory containing results")
    summarize_parser.set_defaults(func=summarize)

    #######################
    ## Library Preparation
    #######################
    libprep_parser = subparsers.add_parser('libprep', help="Only prepare small RNA alignment library sequences that integrates genetic variants (required before processing multiple samples in parallel)")
    libprep_parser.add_argument("-v", "--vcf", dest = "VCF", default = None, required = False, help = "Full path of the input VCF file.")
    ## Optionnal arguments concerning the miRBase library used: matures.fa, hairpins.fa, miRs_coords.gff3. Usefull to switch from miRBase to miRCarta, or a new version of the miRBase
    libprep_parser.add_argument("--maturesFasta", dest="MATURES", default="{}/resources/fasta/hsa_matures_miRBase_v21.fa".format(optimiR_path), required=False, help="Path to the reference library containing mature miRNAs sequences [default: miRBase 21]")
    libprep_parser.add_argument("--hairpinsFasta", dest="HAIRPINS", default="{}/resources/fasta/hsa_hairpins_miRBase_v21.fa".format(optimiR_path), required=False, help="Path to the reference library containing pri-miRNAs sequences [default: miRBase 21]")
    libprep_parser.add_argument("--gff3", dest="GFF3", default="{}/resources/coordinates/hsa_miRBase_v21.gff3".format(optimiR_path), required=False, help="Path to the reference library containing miRNAs and pri-miRNAs coordinates [default: miRBase v21, GRCh38 coordinates]")
    libprep_parser.add_argument("-o", "--dirOutput", dest="OUTDIR", default="./OptimiR_Results_Dir", required=False, help="Full path of the directory where output files are generated [default: ./OptimiR_Results_Dir/]")
    libprep_parser.add_argument("--bowtie2_build", dest="BOWTIE2_BUILD", default="bowtie2-build", required=False, help="Provide path to the bowtie2 index builder binary [default: from $PATH]")
    libprep_parser.set_defaults(func=library_preparation)

    args = parser.parse_args()
    args.func(args)

