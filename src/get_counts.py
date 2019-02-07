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

######################################
## UID managment (adapted from miRTop)
######################################
CODE2NT = {'#':'CCC','$':'GGG', '%':'TTT', '0':'TCT', '3':'GTG', '4':'GAG', '5':'GTT', '6':'GCT', '7':'TGA', '8':'GAC', '9':'TCC', '@':'AAA', 'A':'GGT', 'B':'TGT', 'C':'CGA', 'D':'CAG', 'E':'CGC', 'F':'GAT', 'G':'CGG', 'H':'CTT', 'I':'TGC', 'J':'TAG', 'K':'GGA', 'L':'TAA', 'M':'GGC', 'N':'TAC', 'O':'TTC', 'P':'TCG', 'Q':'TTA', 'R':'TTG', 'S':'CGT', 'T':'GAA', 'U':'TCA', 'V':'GCA', 'W':'GTA', 'X':'GCC', 'Y':'GTC', 'Z':'GCG', 'a':'ACC', 'b':'ATG', 'c':'AAG', 'd':'ACG', 'e':'ATC', 'f':'AAC', 'g':'ATA', 'h':'AGG', 'i':'CCT', 'j':'CTC', 'k':'AGC', 'l':'ACA', 'm':'AGA', 'n':'CAT', 'o':'AAT', 'p':'ATT', 'q':'CTG', 'r':'CTA', 's':'ACT', 't':'CAC', 'u':'TGG', 'v':'CAA', 'w':'AGT', 'x':'CCA', 'y':'CCG', 'z':'TAT'}

NT2CODE = { 'AAA':'@', 'AAC':'f', 'AAG':'c', 'AAT':'o', 'ACA':'l', 'ACC':'a', 'ACG':'d', 'ACT':'s', 'AGA':'m', 'AGC':'k', 'AGG':'h', 'AGT':'w', 'ATA':'g', 'ATC':'e', 'ATG':'b', 'ATT':'p', 'CAA':'v', 'CAC':'t', 'CAG':'D', 'CAT':'n', 'CCA':'x', 'CCC':'#', 'CCG':'y', 'CCT':'i', 'CGA':'C', 'CGC':'E', 'CGG':'G', 'CGT':'S', 'CTA':'r', 'CTC':'j', 'CTG':'q', 'CTT':'H', 'GAA':'T', 'GAC':'8', 'GAG':'4', 'GAT':'F', 'GCA':'V', 'GCC':'X', 'GCG':'Z', 'GCT':'6', 'GGA':'K', 'GGC':'M', 'GGG':'$', 'GGT':'A', 'GTA':'W', 'GTC':'Y', 'GTG':'3', 'GTT':'5', 'TAA':'L', 'TAC':'N', 'TAG':'J', 'TAT':'z', 'TCA':'U', 'TCC':'9', 'TCG':'P', 'TCT':'0', 'TGA':'7', 'TGC':'I', 'TGG':'u', 'TGT':'B', 'TTA':'Q', 'TTC':'O', 'TTG':'R', 'TTT':'%'}

def read_id(uid):
    """
    Read a unique identifier for the sequence and
    converte it to the nucleotides,
    replacing an unique character for 3 nts.
    Inspared in MINTplate: https://cm.jefferson.edu/MINTbase
    https://github.com/TJU-CMC-Org/MINTmap/tree/master/MINTplates
    Args:
        *idu(str)*: unique identifier for the sequence.

    Returns:
        *seq(str)*: nucleotides sequences.
    """
    seq = ""
    for i in uid:
        if i == "1" or i == "2":
            return seq[:-int(i)]
        else:
            seq += CODE2NT[i]
    return seq

def make_id(seq):
    """
    Create a unique identifier for the sequence from the nucleotides,
    replacing 3 nts for another unique character.
    Inspared in MINTplate: https://cm.jefferson.edu/MINTbase
    https://github.com/TJU-CMC-Org/MINTmap/tree/master/MINTplates
    Args:
        *seq(str)*: nucleotides sequences.

    Returns:
        *idName(str)*: unique identifier for the sequence.
    """
    start = 0
    uid = ""
    if not seq:
        raise ValueError("Length of sequence is Empty.")
    for i in range(0, len(seq) + 1, 3):
        if i == 0:
            continue
        trint = seq[start:i]
        uid += NT2CODE[trint]
        start = i
    if len(seq) > i:
        dummy = "A" * (3 - (len(seq) - i))
        trint = seq[i:len(seq)]
        uid += NT2CODE["%s%s" % (trint, dummy)]
        uid += str(len(dummy))
    return uid

##########################################
## GFF DEFINITIONS (todo: ugly code > clean & comment)
##########################################
def get_cigar(seq,d_snv,tail5,tail3,trim5):
    cigar = tail5.upper()
    count = 0
    for i in range(trim5 + 1, len(seq) - len(tail3) +  trim5 - len(tail5) + 1):
        if i in d_snv:
            change = d_snv[i]
            if change == seq[i - 1 - trim5 + len(tail5)]:
                if count > 0:
                    cigar += (str(count) + "M")
                    count = 0 
                cigar += change
            else: count += 1
        else:
            count += 1
    if count > 0:
        cigar += (str(count) + "M")
    cigar += tail3.upper()
    return cigar
    
def process_iso(IT, dO, seq, reference_name):
    iso5, iso3 = IT.split('[')[1].split(']')[0].split(',')
    snvs_in_ref = reference_name.split('_')[1:] if "_" in reference_name else []
    variants = []
    changes = []
    start = 1
    end = len(dO.sequence)
    trim5, trim3 = 0, 0
    tail5, tail3 = "", ""
    is_cano = True
    if "-" in iso5 and "+" in iso5:
        is_cano = False
        trim5 = int(iso5.split('-')[1].split('+')[0])
        start += trim5
        tail5 = iso5.split('+')[1]
        variants.append("iso_5p:+{}".format(trim5))
        variants.append("iso_add5p:{}".format(len(tail5)))
        changes.append("iso_add5p:{}".format(tail5.upper()))
    elif "-" in iso5:
        is_cano = False
        trim5 = int(iso5.split('-')[1])
        start += trim5
        variants.append("iso_5p:+{}".format(trim5))
    elif "+" in iso5:
        is_cano = False
        tail5 = iso5.split('+')[1]
        if tail5.isupper():
            variants.append("iso_5p:-{}".format(len(tail5)))
        else:
            variants.append("iso_add5p:{}".format(len(tail5)))
            changes.append("iso_add5p:{}".format(tail5.upper()))
    if "-" in iso3 and "+" in iso3:
        is_cano = False
        trim3 = int(iso3.split('-')[1].split('+')[0])
        end -= trim3
        tail3 = iso3.split('+')[1]
        variants.append("iso_3p:-{}".format(trim3))
        variants.append("iso_add3p:{}".format(len(tail3)))
        changes.append("iso_add3p:{}".format(tail3.upper()))
    elif "-" in iso3:
        is_cano = False
        trim3 = int(iso3.split('-')[1])
        end -= trim3
        variants.append("iso_3p:-{}".format(trim3))
    elif "+" in iso3:
        is_cano = False
        tail3 = iso3.split('+')[1]
        if tail3.isupper():
            variants.append("iso_3p:+{}".format(len(tail3)))
        else:
            variants.append("iso_add3p:{}".format(len(tail3)))
            changes.append("iso_add3p:{}".format(tail3.upper()))
    #snv (use d_O for variant pos)
    snv_list = dO.variants
    d_snv = {}
    for snv in dO.variants:
        if snv.rsID in snvs_in_ref:
            pos = snv.pos_in_rna + 1
            ref = snv.ref
            alt = snv.alt
            d_snv[pos] = alt
            var = ""
            if 2 <= pos <= 7:
                var = "iso_snv_seed"
            elif 8 == pos :
                var = "iso_snv_central_offset"
            elif 9 <= pos <= 12:
                var = "iso_snv_central"
            elif 13 <= pos <= 17:
                var = "iso_snv_central_supp"
            else:
                var = "iso_snv"
            variants.append("{}:{}".format(var, snv.rsID))
            if len(ref) != 1 or len(alt) != 1:
                changes.append("iso_snv:indel")
            else:
                changes.append("iso_snv:{}{}{}".format(pos, ref, alt))
    if len(variants) > 0:
        variants = "Variant={}".format(",".join(variants))
    else:
        variants = ""
    if len(changes) > 0:
        changes = "Change={}".format(",".join(changes))
    else:
        changes = ""
    cigar = "Cigar={}".format(get_cigar(seq,d_snv,tail5,tail3,trim5))
    return (variants, changes, start, end, cigar, is_cano)

def get_hits(alignments):
    ## retrieve only alignments where XX is HIGHSCORE OR SCORING
    hits = 0
    for a in alignments:
        if a.has_tag("XX"):
            tag = a.get_tag('XX')
            if tag != "DISCARDED_SOFT_CLIPPED_VARIANT" and tag != "DISCARDED_GENO":
                hits += 1
        else:
            hits += 1
    return str(hits)

## Todo : Add details about OptimiR Trimming + Alignment Step + VCF path/name when used
def get_gff_header(database, sample_name):
    if database == "hsa_matures_miRBase_v21":
        database = "miRBasev21 doi:10.25504/fairsharing.hmgte8"
    elif database == "hsa_matures_miRBase_v22":
        database = "miRBasev22 doi:10.25504/fairsharing.hmgte8"
    ## else : custom database
    return "## GFF3 adapted for miRNA sequencing data. VERSION 1.1\n## Source: {}\n## COLDATA: {}\n## Processed with OptimiR.\n".format(database, sample_name)

## Essentially remove polymiRs ambiguous that were discarded due to scoring instead of genotype consistency. Happens when variants are localised on extremities and are soft clipped.
## If polymiR in list (at least 2 sequences distinguished only by the presence of variant(s)), first discard the ones with inconsistent geno (if available), then keep the alignment with best score
## If same score and no geno? only keep ref for GFF
def clean_alignments(alignments):
    alignments_out = []
    ref_set = {}
    polymiRs = {} # dict in case one read aligns on several polymiRs (unlikely but with miRs you never know > actually we know, with miR-1269a/b)
    for a in alignments:
        ## If reference name is present twice, then one has rsID and one has a better score
        ref = a.reference_name.split('_')[0]
        if ref in ref_set: ## then polymiR
            if a.has_tag('GC'):
                if a.get_tag('GC') == "Genotype_INVALID":
                    continue
            if ref in polymiRs:
                polymiRs[ref].append((a, int(a.get_tag('SC'))))
            else:
                polymiRs[ref] = []
                polymiRs[ref].append((ref_set[ref], int(ref_set[ref].get_tag('SC'))))
                polymiRs[ref].append((a, int(a.get_tag('SC'))))
        else:
            ref_set[ref] = a
    for ref_name, p_l in polymiRs.items():
        p_l = sorted(p_l, key=lambda x: x[1]) ## ascending order
        best_score = p_l[0][1]
        found_ref = False
        for p in p_l:
            if not("_" in p[0].reference_name) and (p[1] == best_score):
                found_ref = True
                alignments_out.append(p[0])
        if not(found_ref): ## if reference is not found with best score, only take first
            alignments_out.append(p_l[0][0])
    ## Add the non-polymiRs
    for ref, a in ref_set.items():
        if ref not in polymiRs:
            alignments_out.append(a)
    return alignments_out
    
## Todo: rewrite and cleanup
def write_gff(bam_dict, collapse_table, sample_name, d_OptimiR, out_fn, sourceDB):
    alignments = [a for alignments in bam_dict.values() for a in alignments]
    lines = []
    with open(out_fn, "w") as out:
        out.write(get_gff_header(sourceDB, sample_name))
        for alignments in bam_dict.values():
            alignments = clean_alignments(alignments) 
            hits = "Hits={}".format(get_hits(alignments))
            ### For each alignment
            for a in alignments:
                seqid = a.reference_name.split('_')[0] ## without rsID
                variant, changes, start, end, cigar, is_cano = process_iso(a.get_tag('IT'), d_OptimiR[a.reference_name], a.seq, a.reference_name)
                if sourceDB == "hsa_matures_miRBase_v21":
                    source = "mirbase_21"
                elif sourceDB == "hsa_matures_miRBase_v22":
                    source = "mirbase_22"
                else:
                    source = sourceDB ## retrieve from provided filename
                if is_cano:
                    typ = "ref_miRNA"
                else:
                    typ = "isomiR"
                score = str(a.get_tag('SC'))
                strand = "+"
                phase = "."
                UID = "UID={}".format(make_id(a.seq))
                read = "Read={}".format(a.seq)
                parent = "Parent={}".format(a.get_tag('PA'))
                name = "Name={}".format(a.reference_name)
                alias = "Alias={}".format(",".join(list(set([i.split('_')[0] for i in d_OptimiR[seqid].ident]))))
                expression = "Expression={}".format(float(a.get_tag('XC')))
                expression_optimir = "Expression_OptimiR={}".format(float(a.get_tag('XW')) * float(a.get_tag('XC')))
                if a.has_tag('XX'):
                    tag = a.get_tag('XX')
                    if tag != "DISCARDED_SOFT_CLIPPED_VARIANT" and tag != "DISCARDED_GENO":
                        filtr = "REJECT," + tag
                    else:
                        filtr = "DELETE_GENO"
                else:
                    filtr = "PASS"
                filtr = "Filter={}".format(filtr)
                if variant == "":
                    attributes = ";".join([UID,read,parent,name,cigar,alias,expression,expression_optimir,filtr, hits])
                elif changes == "":
                    attributes = ";".join([UID,read,parent,name,variant,cigar,alias,expression,expression_optimir,filtr, hits])
                else:
                    attributes = ";".join([UID,read,parent,name,variant,changes,cigar,alias,expression,expression_optimir,filtr, hits])
                line = "\t".join([seqid, source, typ, str(start), str(end), score, strand, phase, attributes])
                if filtr != "Filter=DELETE_GENO":
                    lines.append(line)
        lines.sort()
        for line in lines:
            out.write(line + "\n")
            
def compute_abundances(bam_dict, collapse_table, sample_name_list, out_dir, out_filename, d_OptimiR, ANNOT_FILES, WRITE_GFF, sourceDB):
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
    if WRITE_GFF:
        out_fn = out_dir + "/" + sample_name_list[0] + ".gff3"
        write_gff(bam_dict, collapse_table, sample_name_list[0], d_OptimiR, out_fn, sourceDB)

