#-*- coding: utf-8 -*-
# Florian THIBORD (22/08/17)
########################################################
#                  OptimiR PIPELINE                    # 
########################################################

# Standard libraries
import sys, os
from pysam import AlignmentFile

# Personal libraries
import optimir.libs.annotate as ANNOT
import optimir.libs.scoring as SCORE
import optimir.libs.consistency as CONSIST
import optimir.libs.get_counts as ABUNDANCE
from optimir.libs.essentials import *

## POST-PROCESSING (for reads mapping on miRs only)
# (Weighting, Isotyping, Genotype consistency, Abundance)
def post_process_main(tmpdir_mapping, tmpdir_postProcess, dir_results, collapse_table, sample_name, WEIGHT5, SCORE_THRESHOLD, INCONSISTENT_THRESHOLD, d_OptimiR_pickle_path, VCF_AVAILABLE, GENO_AVAILABLE, ANNOT_FILES, VERBOSE, WRITE_GFF, WRITE_VCF, sourceDB):
    # Load pickle object
    try:
        d_OptimiR = load_obj(d_OptimiR_pickle_path)
    except ImportError: ## Happens when library_preparation is called ahead
        sys.modules['essentials'] = sys.modules['src.essentials']
        d_OptimiR = load_obj(d_OptimiR_pickle_path)
    # Create dict from bam file. keys=read_name values=alignments for that read
    bam = AlignmentFile('{}/{}.bam'.format(tmpdir_mapping, sample_name), 'rb')
    bam_dict = bamFile_to_dict(bam)
    # Compute isotype, look for parental hairpin, add IT tag with isotype, PA tag with parental hairpin, and I2 tag with isotypes on different hairpins
    ANNOT.isotyper(bam_dict, d_OptimiR, WEIGHT5)
    # Compute alignment score, add SC tag
    fun_str_progress([], "alscore", VERBOSE)
    SCORE.get_alignment_score(bam_dict, WEIGHT5)
    # Discard alignments above SCORE_THRESHOLD, add XX tag with DISCARDED_HIGHSCORE
    SCORE.discard_above_score_threshold(bam_dict, SCORE_THRESHOLD)
    if VCF_AVAILABLE:
        # Discard alignment with soft clip over variant in polymIR, add XX tag with DISCARDED_SOFT_CLIPPED_VARIANT
        CONSIST.check_no_soft_clip_on_variant(bam_dict, sample_name, d_OptimiR)
        if GENO_AVAILABLE:
            # Compute genotype consistency, if discarded add XX tag with DISCARDED_GENO, add GL tag with sample genotypes and GC tag with consistency status
            fun_str_progress([], "geno", VERBOSE)
            CONSIST.compute_genotype_consistency(bam_dict, d_OptimiR, sample_name)
    # Cross mapping reads resolution : the alignment with the best score is kept, other are discarded with tag XX and DISCARDED_SCORING
    fun_str_progress([], "ambiguous", VERBOSE)
    SCORE.resolve_ambiguous_alignments(bam_dict, SCORE_THRESHOLD)
    # WEIGHT computed after possible reads discarding step
    # Compute weights, add XW tag
    SCORE.compute_weight(bam_dict)
    # Add MU tag if cross-mapping alignments, and XC tag with count of reads (retrieved from read name)
    ANNOT.add_tags(bam_dict)
    # Compute outputs: Abundances, discarded.sam
    out_fn = "{}/{}.Isotyped.bam".format(tmpdir_postProcess, sample_name)
    bam_out = AlignmentFile(out_fn, 'wb', template=bam)
    bamDict_to_file(bam_dict, bam_out, collapse_table)
    bam_out.close()
    ## Write abundances and annotation files
    fun_str_progress([], "outputs", VERBOSE)
    ABUNDANCE.compute_abundances(bam_dict, collapse_table, [sample_name], dir_results, '{}_abundances.txt'.format(sample_name), d_OptimiR, ANNOT_FILES, WRITE_GFF, sourceDB)
    ## Write polymiRs outputs (inconsistents + polymiR tables + optional VCF)
    CONSIST.write_polymiRs_outputs(bam_dict, bam, collapse_table, sample_name, dir_results, d_OptimiR, INCONSISTENT_THRESHOLD, ANNOT_FILES, WRITE_VCF)
    if "s" in ANNOT_FILES:
        write_isomiR_dist(bam_dict, sample_name, dir_results)
    bam.close()
