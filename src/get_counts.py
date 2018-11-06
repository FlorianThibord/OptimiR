#-*- coding: utf-8 -*-
# Florian THIBORD (20/06/17)
########################################################
#                 OptimiR abundances                   #
########################################################
# Designed to also work with sample list, but only used with 1 element - list with OptimiR

from pysam import AlignmentFile
from collections import Counter

from . import essentials as ES

class annot: pass

def remove_ref_from_cross_mapp(cross_mapp, ref):
    return '+'.join([r for r in cross_mapp.split('/') if not(ref in r)])

def fill_annot_dict(annot_dict, alignment, ref):
    if alignment.has_tag('PA'):
        parent = alignment.get_tag('PA')
    else:
        parent = ''
    if alignment.has_tag('MU'):
        cross_mapp = alignment.get_tag('MU')
    else:
        cross_mapp = ''
    if alignment.has_tag('XW'):
        weight = float(alignment.get_tag('XW'))
    else:
        weight = 1
    if alignment.has_tag('XC'):
        copies = float(alignment.get_tag('XC'))
    else:
        copies = 1
    cross_mapp = remove_ref_from_cross_mapp(cross_mapp, ref)
    parent = '+'.join(parent.split('/'))
    counts = weight * copies
    if ref in annot_dict:
        if cross_mapp != '' and counts > 0:
            if cross_mapp in annot_dict[ref].cross_mapping:
                annot_dict[ref].cross_mapping[cross_mapp] += counts
            else:
                annot_dict[ref].cross_mapping[cross_mapp] = counts
        if parent != '' and counts > 0:
            if parent in annot_dict[ref].parents:
                annot_dict[ref].parents[parent] += counts
            else:
                annot_dict[ref].parents[parent] = counts
        annot_dict[ref].counts += counts
    else:
        annot_dict[ref] = annot()
        if cross_mapp != '' and counts > 0:
            annot_dict[ref].cross_mapping = {cross_mapp : counts}
        else:
            annot_dict[ref].cross_mapping = {}
        if parent != '' and counts > 0:
            annot_dict[ref].parents = {parent : counts}
        else:
            annot_dict[ref].parents = {}
        annot_dict[ref].counts = counts

def get_annot(annot_dict, ref):
    counts = annot_dict[ref].counts
    # CM is global, we only want to know where it cross-mapp and what percentage of reads cross-mapp there
    CM = annot_dict[ref].cross_mapping
    l_cross_mapp = []
    for cm_locus, nb in CM.items():
        where_cm = "+".join(cm_locus.split('/'))
        if nb > 0:
            l_cross_mapp.append('{}[{}%|{}reads]'.format(where_cm, '%.2f' % (nb*100/counts), nb))
    l_cross_mapp = '/'.join(l_cross_mapp)
    # PA same for parents hairpin
    PA = annot_dict[ref].parents
    l_parents = []
    for pa_ids, nb in PA.items():
        l_parents.append('{}[{}%|{}reads]'.format(pa_ids, '%.2f' % (nb*100/counts), nb))
    l_parents = '/'.join(l_parents)
    return l_parents + '\t' + l_cross_mapp 

def retrieve_counts_per_samples(collapse_table, sample_name_list, read, weight=1):
    counts_samples = Counter()
    sequence = collapse_table[read]
    sample_name = sample_name_list[0]
    counts_samples[sample_name] = sequence.times * weight
    return counts_samples

def compute_counts_miRs(bam_dict, collapse_table, sample_name_list):
    counts = {}
    annot_dict = {}
    alignments = [a for alignments in bam_dict.values() for a in alignments]
    for a in alignments:
        ref = a.reference_name
        if '_' in ref:
            # if polymiR, only one ref with all counts
            ref = ref.split('_')[0]
        read = a.query_sequence
        weight = float(a.get_tag('XW'))
        fill_annot_dict(annot_dict, a, ref)
        if ref not in counts:
            counts[ref] = Counter()
        counts[ref] += retrieve_counts_per_samples(collapse_table, sample_name_list, read, weight)    
    return counts, annot_dict

def compute_counts_isomiRs(bam_dict, collapse_table, sample_name_list):
    counts = {}
    annot_dict = {}
    alignments = [a for alignments in bam_dict.values() for a in alignments]
    cpt = 0
    for a in alignments:
        ref = a.reference_name
        if a.has_tag("IT"):
            cpt += 1
            iso_ref = "{}~{}".format(ref, a.get_tag('IT'))
            read = a.query_sequence
            weight = float(a.get_tag('XW'))
            fill_annot_dict(annot_dict, a, iso_ref)
            if iso_ref not in counts:
                counts[iso_ref] = Counter()
            counts[iso_ref] += retrieve_counts_per_samples(collapse_table, sample_name_list, read, weight)
    return counts, annot_dict

def compute_counts_polymiRs(bam_dict, collapse_table, sample_name_list, d_OptimiR):
    ## polymiR idents (ref and alternative)
    polymiRs = [key for key, obj in d_OptimiR.items() if len(obj.variants) > 0]
    counts = {}
    annot_dict = {}
    alignments = [a for alignments in bam_dict.values() for a in alignments]
    for a in alignments:
        ref = a.reference_name
        if ref in polymiRs:
            read = a.query_sequence
            weight = float(a.get_tag('XW'))
            fill_annot_dict(annot_dict, a, ref)
            if ref not in counts:
                counts[ref] = Counter()
            counts[ref] += retrieve_counts_per_samples(collapse_table, sample_name_list, read, weight)
    return counts, annot_dict

def compute_counts_miRs_and_polymiRs(bam_dict, collapse_table, sample_name_list):
    counts = {}
    annot_dict = {}
    alignments = [a for alignments in bam_dict.values() for a in alignments]
    for a in alignments:
        ref = a.reference_name
        read = a.query_sequence
        weight = float(a.get_tag('XW'))
        fill_annot_dict(annot_dict, a, ref)
        if ref not in counts:
            counts[ref] = Counter()
        counts[ref] += retrieve_counts_per_samples(collapse_table, sample_name_list, read, weight)
    return counts, annot_dict

def write_counts(counts, sample_name_list, out_fn, annot = {}):
    out = open(out_fn, 'w')
    # write header
    l = 'reference'
    for sample in sample_name_list:
        l += '\t{}'.format(sample)
    if annot != {}:
        l += '\tPARENTAL_HAIRPIN\tCROSS-MAPPING'
    out.write(l + '\n')
    for ref, counts_per_sample in counts.items():
        l = ref
        total_counts = 0
        for sample in sample_name_list:
            counts = counts_per_sample[sample]
            total_counts += counts
            l += '\t{}'.format(counts)
        if annot != {}:
            l += '\t' + get_annot(annot, ref)
        if total_counts > 0:
            out.write(l + '\n')
    out.close()

def write_annotation_files(sample_name_list, out_fn_cross_mapp, out_fn_hairpins, counts, annot, ANNOT_FILES):
    if annot != {}:
        if "h" in ANNOT_FILES:
            out_hp = open(out_fn_hairpins, 'w')
            # write header
            l = 'reference'
            for sample in sample_name_list:
                l += '\t{}'.format(sample)
            out_hp.write(l + '\n')
            for ref, counts_per_sample in counts.items():
                total_counts = 0
                for sample in sample_name_list:
                    count = counts_per_sample[sample]
                    total_counts += count
                l = ref
                annots = get_annot(annot, ref)
                hp, cm = annots.split('\t')
                if total_counts > 0:
                    out_hp.write(l + '\t' + hp + '\n')
            out_hp.close()
        if "c" in ANNOT_FILES:
            out_cm = open(out_fn_cross_mapp, 'w')
            # write header
            l = 'reference'
            for sample in sample_name_list:
                l += '\t{}'.format(sample)
            out_cm.write(l + '\n')
            for ref, counts_per_sample in counts.items():
                total_counts = 0
                for sample in sample_name_list:
                    count = counts_per_sample[sample]
                    total_counts += count
                l = ref
                annots = get_annot(annot, ref)
                hp, cm = annots.split('\t')
                if total_counts > 0:
                    out_cm.write(l + '\t' + cm + '\n')
            out_cm.close()

def compute_abundances(bam_dict, collapse_table, sample_name_list, out_dir, out_filename, d_OptimiR, ANNOT_FILES):
    counts, annot = compute_counts_miRs(bam_dict, collapse_table, sample_name_list)
    out_fn = out_dir + '/' + "miRs_merged." + out_filename
    write_counts(counts, sample_name_list, out_fn, annot)
    counts, annot = compute_counts_isomiRs(bam_dict, collapse_table, sample_name_list)
    out_fn = out_dir + '/' + "isomiRs_specific." + out_filename
    write_counts(counts, sample_name_list, out_fn, annot)
    counts, annot = compute_counts_miRs_and_polymiRs(bam_dict, collapse_table, sample_name_list)
    out_fn = out_dir + '/' + "miRs_and_polymiRs." + out_filename
    write_counts(counts, sample_name_list, out_fn, annot)
    out_fn_cross_mapp = out_dir + '/' + "remaining_ambiguous." + sample_name_list[0] + ".annot"
    out_fn_hairpins = out_dir + '/' + "expressed_hairpins." + sample_name_list[0] + ".annot" 
    write_annotation_files(sample_name_list, out_fn_cross_mapp, out_fn_hairpins, counts, annot, ANNOT_FILES)
    counts, annot = compute_counts_polymiRs(bam_dict, collapse_table, sample_name_list, d_OptimiR)
    out_fn = out_dir + '/' + "polymiRs_specific." + out_filename
    write_counts(counts, sample_name_list, out_fn, annot)
    
