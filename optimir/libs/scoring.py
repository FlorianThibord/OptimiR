#-*- coding: utf-8 -*-
# Florian THIBORD (20/06/17)
########################################################
#                 OptimiR isotyper                     #
########################################################

from pysam import AlignmentFile

## Personal library
from optimir.libs.essentials import *

def compute_isotype_score(isotype, WEIGHT5):
    def get_end_nb_modif(iso_end):
        nb_trim = 0
        nb_tail = 0
        if "-" in iso_end:
            nb_trim = int(iso_end.split('-')[1].split('+')[0])
        if "+" in iso_end:
            tail = iso_end.split('+')[1]
            if tail.islower(): ## then non templated addition
                nb_tail = len(tail)
        return nb_trim + nb_tail
    iso5, iso3 = isotype.split('[')[1].split(']')[0].split(',')
    score = get_end_nb_modif(iso5) * WEIGHT5 + get_end_nb_modif(iso3)
    return score

def get_alignment_score(bam_dict, WEIGHT5):
    alignments = [a for alignments in bam_dict.values() for a in alignments]
    for alignment in alignments:
        isotype = alignment.get_tag("IT")
        score = compute_isotype_score(isotype, WEIGHT5)
        alignment.set_tag("SC", score)
      
def discard_above_score_threshold(bam_dict, SCORE_THRESHOLD):
    alignments = [a for alignments in bam_dict.values() for a in alignments]
    for alignment in alignments:
        score = int(alignment.get_tag("SC"))
        if score > SCORE_THRESHOLD:
            alignment.set_tag("XX", "DISCARDED_HIGHSCORE")


## SPECIAL CASE with polymiRs: if align on alt and reference, and one has been discarded because of inconsistent genotype, keep the other one even if score is not lower
def resolve_ambiguous_alignments(bam_dict, SCORE_THRESHOLD):
    for alignments in bam_dict.values():
        best_score = SCORE_THRESHOLD
        nb_discarded_with_bad_geno = 0
        # 1) get best score
        for alignment in alignments:
            if not(alignment.has_tag("XX")):
                score = int(alignment.get_tag("SC"))
                best_score = score if score < best_score else best_score
        # 2) discard alignment with score > best_score using tag XX and "DISCARDED_SCORING"
        for alignment in alignments:
            score = int(alignment.get_tag("SC"))
            if score > best_score:
                alignment.set_tag("XX", "DISCARDED_SCORING")
                
#Compute weight from read score
def compute_weight(bam_dict):
    for alignments in bam_dict.values():
        nb_alignments = 0
        for alignment in alignments:
            if not(alignment.has_tag("XX")): ## if alignment has not been discarded
                nb_alignments += 1
        if nb_alignments > 0:
            weight = 1. / nb_alignments
        else:
            weight = 0
        for alignment in alignments:
            if not(alignment.has_tag("XX")): ## if alignment has not been discarded
                alignment.set_tag("XW", str(weight))
            else:
                alignment.set_tag("XW", str(0))
                
