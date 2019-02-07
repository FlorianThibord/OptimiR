#-*- coding: utf-8 -*-
# Florian THIBORD (20/06/17)
########################################################
#                 OptimiR isotyper                     #
########################################################

from pysam import AlignmentFile

## Personal library
from .essentials import *
from . import scoring as SCORE

def trimming_in_5end(alignment):
    pos = alignment.reference_start
    trim5_flag = ""
    nb_5 = 0
    if pos > 0:
        x, sc_offset5 = alignment.cigartuples[0]
        # Is there soft clipping in 5' (pos 0 in CIGAR corresponds to 5' end)?
        if x == 4: ## code for soft-clip
            ## INDICATE TRIMMING FOLLOWED BY TAILING!
            trim5_flag = "-+"
            letters_tailed = alignment.seq[:sc_offset5]
            nb_trimmed = pos
            nb_5 = (nb_trimmed, letters_tailed)
        else: ## It's 5' Trimming !
            trim5_flag = "-" 
            nb_5 = pos
    return (trim5_flag, nb_5)

def trimming_in_3end(alignment, reference_seq):
    # Difficulty: we don't know if the query stops before the end of the reference, so we need to compute it (with len(ref_seq, query_seq), start_pos, soft-clip in 5' and 3')
    pos = alignment.reference_start
    trim3_flag = ""
    nb_3 = 0
    # Is there soft clipping in 3' (pos -1 in CIGAR corresponds to 3' end)?
    x3, sc_offset3 = alignment.cigartuples[-1]
    if x3 != 4:
        sc_offset3 = 0
    x5, sc_offset5 = alignment.cigartuples[0]
    if x5 != 4:
        sc_offset5 = 0
    if len(reference_seq) > len(alignment.seq) + pos - sc_offset5 - sc_offset3:
        if x3 == 4: ## code for soft-clip
            ## INDICATE TRIMMING FOLLOWED BY TAILING!
            trim3_flag = "-+"
            nb_trimmed = len(reference_seq) - (len(alignment.seq) + pos - sc_offset5 - sc_offset3)
            letters_tailed = alignment.seq[-sc_offset3:]
            nb_3 = (nb_trimmed, letters_tailed)
        else:
            trim3_flag = "-" #trim
            nb_3 = len(reference_seq) - (len(alignment.seq) + pos - sc_offset5)
    return (trim3_flag, nb_3)

def return_tail(tail_read, tail_ref): ## python3 made it weird with lambda tuples...
    return ''.join(map(lambda mytuple: mytuple[0] if mytuple[0] == mytuple[1] else mytuple[0].lower(), list(zip(tail_read, tail_ref)))) 

# Must be called after trimming
# Compute tail only (no trim + tail)
# Optional: hairpin_seq to check if tail is templated
def tailing_in_5end(alignment, reference_seq, hairpin_seq=''):
    tail5_flag = ""
    letters_5 = ''
    cigar = alignment.cigarstring
    if 'S' in cigar:
        x5, sc_offset5 = alignment.cigartuples[0]
        if (x5 == 4) and not(alignment.reference_start > 0): # else can't be tailing
            tail5_flag = "+" #tail
            letters_5 = alignment.seq[:sc_offset5]
            if hairpin_seq != '':
                hairpin_start = hairpin_seq.find(reference_seq)
                tail_in_hairpin = hairpin_seq[(hairpin_start - (len(letters_5))):hairpin_start]
                tail_in_read = alignment.seq[:len(letters_5)]
                if tail_in_read == tail_in_hairpin:
                    tail5_flag = "+[TE]"
    return (tail5_flag, letters_5)

def tailing_in_3end(alignment, reference_seq, hairpin_seq=''):
    tail3_flag = ""
    letters_3 = ''
    pos = alignment.reference_start
    cigar = alignment.cigarstring
    if 'S' in cigar:
        x3, sc_offset3 = alignment.cigartuples[-1]
        x5, sc_offset5 = alignment.cigartuples[0]
        if x5 != 4: # then it's not soft-clip
            sc_offset5 = 0
        if (x3 == 4) and (len(reference_seq) <= len(alignment.seq) + pos - sc_offset5 - sc_offset3): # else can't be tailing
            tail3_flag = "+" #tail
            letters_3 = alignment.seq[-sc_offset3:]
            if hairpin_seq != '':
                hairpin_end_tmp = hairpin_seq.find(reference_seq)
                hairpin_end = hairpin_end_tmp + len(reference_seq)
                tail_in_hairpin = hairpin_seq[hairpin_end:(hairpin_end + len(letters_3))]
                tail_in_read = alignment.seq[-len(letters_3):]
                if tail_in_read == tail_in_hairpin:
                    tail3_flag = "+[TE]"
    return (tail3_flag, letters_3)

def compute_isoform(a, reference_seq, hairpin_seq=''):
    flag5, end_5 = trimming_in_5end(a)
    flag3, end_3 = trimming_in_3end(a, reference_seq)
    if flag5 == "": #if no trimming, then search for tailing
        flag5, end_5 = tailing_in_5end(a, reference_seq, hairpin_seq)
    if flag3 == "": #if no trimming, then search for tailing
        flag3, end_3 = tailing_in_3end(a, reference_seq, hairpin_seq)
    iso5, iso3 = "", ""
    if flag5 == "-+":
        iso5 = "-{}+{}".format(end_5[0], end_5[1].lower())
    elif flag5 == "-":
        iso5 = "-{}".format(end_5)
    elif flag5 == "+":
        iso5 = "+{}".format(end_5.lower())
    elif flag5 == "+[TE]":
        iso5 = "+{}".format(end_5.upper())
    if flag3 == "-+":
        iso3 = "-{}+{}".format(end_3[0], end_3[1].lower())
    elif flag3 == "-":
        iso3 = "-{}".format(end_3)
    elif flag3 == "+":
        iso3 = "+{}".format(end_3.lower())
    elif flag3 == "+[TE]":
        iso3 = "+{}".format(end_3.upper())
    return '[{},{}]'.format(iso5, iso3)

def multi_orig_hairpin(isotype_d, WEIGHT5):
    names, isotypes = list(isotype_d.keys()), list(isotype_d.values())
    scores = [SCORE.compute_isotype_score(isotype,WEIGHT5) for isotype in isotypes]
    # if isotypes identical return only one
    identical = True
    for isotype in isotypes:
        identical = identical and (isotype == isotypes[0])
    # elif templated, resolution
    templated = {}
    for isotype, name in zip(isotypes, names):
        if any(char.isupper() for char in isotype):
            templated[name] = isotype
    # else return everything
    if identical:
        names_str = names[0]
        if len(names) > 1:
            for n in names[1:]:
                names_str += '/{}'.format(n)
        return (names_str, isotypes[0])
    elif len(templated) == 1:
        return (templated.keys()[0], templated.values()[0])
    else:
        names,isotypes,score_sorted = zip(*sorted(zip(names,isotypes,scores), key=lambda tup: tup[2]))
        names_str = names[0]
        isotypes_str = isotypes[0]
        if len(names) > 1:
            for n,i in zip(names[1:],isotypes[1:]):
                names_str += '/{}'.format(n)
                isotypes_str += '/{}'.format(i)
        return (names_str, isotypes_str)

## Update bam with isotype
def isotyper(bam_dict, d_OptimiR, WEIGHT5):
    # flatten bam_dict.values() to a list of alignments
    alignments = [a for alignments in bam_dict.values() for a in alignments]
    for alignment in alignments:
        read_name = alignment.query_name
        reference_name = alignment.reference_name
        OptimiR_obj = d_OptimiR[reference_name]
        reference_sequence = OptimiR_obj.sequence
        d_hairpins = OptimiR_obj.hairpins
        d_isotype = {}
        for hairpin_name, hairpin_sequence in d_hairpins.items():
            isotype = compute_isoform(alignment, reference_sequence, hairpin_sequence)
            d_isotype[hairpin_name] = isotype
        hairpin_names, isotypes = multi_orig_hairpin(d_isotype, WEIGHT5)
        alignment.set_tag("IT", "{}".format(isotypes.split('/')[0]))
        if len(isotypes.split('/')) > 1: ## if different isotypes on different hairpins, add them
            alignment.set_tag("I2", "{}".format(isotypes))
            alignment.set_tag("PA", "{}".format(hairpin_names.split('/')[0])) 
            alignment.set_tag("P2", "{}".format(hairpin_names)) ## add all potential hairpins parents
        else:
            alignment.set_tag("PA", "{}".format(hairpin_names)) ## add all potential hairpins parents

# Add XC tag with counts for each read, and remove it from read name
# If multiple references, add tag MU with list of references
def add_tags(bam_dict):
    for query_name, alignments in bam_dict.items():
        query_name, counts = query_name.split('_')
        for a in alignments:
            a.set_tag('XC', int(counts))
        reference_names = [a.reference_name for a in alignments if not(a.has_tag("XX"))]
        if len(reference_names) > 1:
            for a in alignments:
                a.set_tag('MU', "{}".format("/".join(reference_names)))
