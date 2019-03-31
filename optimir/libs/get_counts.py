#-*- coding: utf-8 -*-
# Florian THIBORD (20/06/17)
########################################################
#                 OptimiR abundances                   #
########################################################
# Designed to also work with sample list, but only used with 1 element - list with OptimiR

from pysam import AlignmentFile
from collections import Counter

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
encode_hash = {'AAAAA': 'BB', 'AAAAC': 'BD', 'AAAAG': 'B0', 'AAAAT': 'BE', 'AAACA': 'BF', 'AAACC': 'B1',
               'AAACG': 'BH', 'AAACT': 'BI', 'AAAGA': 'B2', 'AAAGC': 'BJ', 'AAAGG': 'BK', 'AAAGT': 'B3',
               'AAATA': 'BL', 'AAATC': 'BM', 'AAATG': 'B4', 'AAATT': 'BN', 'AACAA': 'BO', 'AACAC': 'B5',
               'AACAG': 'BP', 'AACAT': 'BQ', 'AACCA': 'B6', 'AACCC': 'BR', 'AACCG': 'BS', 'AACCT': 'B7',
               'AACGA': 'BU', 'AACGC': 'BV', 'AACGG': 'B8', 'AACGT': 'BW', 'AACTA': 'BX', 'AACTC': 'B9',
               'AACTG': 'BY', 'AACTT': 'BZ', 'AAGAA': 'DB', 'AAGAC': 'DD', 'AAGAG': 'D0', 'AAGAT': 'DE',
               'AAGCA': 'DF', 'AAGCC': 'D1', 'AAGCG': 'DH', 'AAGCT': 'DI', 'AAGGA': 'D2', 'AAGGC': 'DJ',
               'AAGGG': 'DK', 'AAGGT': 'D3', 'AAGTA': 'DL', 'AAGTC': 'DM', 'AAGTG': 'D4', 'AAGTT': 'DN',
               'AATAA': 'DO', 'AATAC': 'D5', 'AATAG': 'DP', 'AATAT': 'DQ', 'AATCA': 'D6', 'AATCC': 'DR',
               'AATCG': 'DS', 'AATCT': 'D7', 'AATGA': 'DU', 'AATGC': 'DV', 'AATGG': 'D8', 'AATGT': 'DW',
               'AATTA': 'DX', 'AATTC': 'D9', 'AATTG': 'DY', 'AATTT': 'DZ', 'ACAAA': '0B', 'ACAAC': '0D',
               'ACAAG': '00', 'ACAAT': '0E', 'ACACA': '0F', 'ACACC': '01', 'ACACG': '0H', 'ACACT': '0I',
               'ACAGA': '02', 'ACAGC': '0J', 'ACAGG': '0K', 'ACAGT': '03', 'ACATA': '0L', 'ACATC': '0M',
               'ACATG': '04', 'ACATT': '0N', 'ACCAA': '0O', 'ACCAC': '05', 'ACCAG': '0P', 'ACCAT': '0Q',
               'ACCCA': '06', 'ACCCC': '0R', 'ACCCG': '0S', 'ACCCT': '07', 'ACCGA': '0U', 'ACCGC': '0V',
               'ACCGG': '08', 'ACCGT': '0W', 'ACCTA': '0X', 'ACCTC': '09', 'ACCTG': '0Y', 'ACCTT': '0Z',
               'ACGAA': 'EB', 'ACGAC': 'ED', 'ACGAG': 'E0', 'ACGAT': 'EE', 'ACGCA': 'EF', 'ACGCC': 'E1',
               'ACGCG': 'EH', 'ACGCT': 'EI', 'ACGGA': 'E2', 'ACGGC': 'EJ', 'ACGGG': 'EK', 'ACGGT': 'E3',
               'ACGTA': 'EL', 'ACGTC': 'EM', 'ACGTG': 'E4', 'ACGTT': 'EN', 'ACTAA': 'EO', 'ACTAC': 'E5',
               'ACTAG': 'EP', 'ACTAT': 'EQ', 'ACTCA': 'E6', 'ACTCC': 'ER', 'ACTCG': 'ES', 'ACTCT': 'E7',
               'ACTGA': 'EU', 'ACTGC': 'EV', 'ACTGG': 'E8', 'ACTGT': 'EW', 'ACTTA': 'EX', 'ACTTC': 'E9',
               'ACTTG': 'EY', 'ACTTT': 'EZ', 'AGAAA': 'FB', 'AGAAC': 'FD', 'AGAAG': 'F0', 'AGAAT': 'FE',
               'AGACA': 'FF', 'AGACC': 'F1', 'AGACG': 'FH', 'AGACT': 'FI', 'AGAGA': 'F2', 'AGAGC': 'FJ',
               'AGAGG': 'FK', 'AGAGT': 'F3', 'AGATA': 'FL', 'AGATC': 'FM', 'AGATG': 'F4', 'AGATT': 'FN',
               'AGCAA': 'FO', 'AGCAC': 'F5', 'AGCAG': 'FP', 'AGCAT': 'FQ', 'AGCCA': 'F6', 'AGCCC': 'FR',
               'AGCCG': 'FS', 'AGCCT': 'F7', 'AGCGA': 'FU', 'AGCGC': 'FV', 'AGCGG': 'F8', 'AGCGT': 'FW',
               'AGCTA': 'FX', 'AGCTC': 'F9', 'AGCTG': 'FY', 'AGCTT': 'FZ', 'AGGAA': '1B', 'AGGAC': '1D',
               'AGGAG': '10', 'AGGAT': '1E', 'AGGCA': '1F', 'AGGCC': '11', 'AGGCG': '1H', 'AGGCT': '1I',
               'AGGGA': '12', 'AGGGC': '1J', 'AGGGG': '1K', 'AGGGT': '13', 'AGGTA': '1L', 'AGGTC': '1M',
               'AGGTG': '14', 'AGGTT': '1N', 'AGTAA': '1O', 'AGTAC': '15', 'AGTAG': '1P', 'AGTAT': '1Q',
               'AGTCA': '16', 'AGTCC': '1R', 'AGTCG': '1S', 'AGTCT': '17', 'AGTGA': '1U', 'AGTGC': '1V',
               'AGTGG': '18', 'AGTGT': '1W', 'AGTTA': '1X', 'AGTTC': '19', 'AGTTG': '1Y', 'AGTTT': '1Z',
               'ATAAA': 'HB', 'ATAAC': 'HD', 'ATAAG': 'H0', 'ATAAT': 'HE', 'ATACA': 'HF', 'ATACC': 'H1',
               'ATACG': 'HH', 'ATACT': 'HI', 'ATAGA': 'H2', 'ATAGC': 'HJ', 'ATAGG': 'HK', 'ATAGT': 'H3',
               'ATATA': 'HL', 'ATATC': 'HM', 'ATATG': 'H4', 'ATATT': 'HN', 'ATCAA': 'HO', 'ATCAC': 'H5',
               'ATCAG': 'HP', 'ATCAT': 'HQ', 'ATCCA': 'H6', 'ATCCC': 'HR', 'ATCCG': 'HS', 'ATCCT': 'H7',
               'ATCGA': 'HU', 'ATCGC': 'HV', 'ATCGG': 'H8', 'ATCGT': 'HW', 'ATCTA': 'HX', 'ATCTC': 'H9',
               'ATCTG': 'HY', 'ATCTT': 'HZ', 'ATGAA': 'IB', 'ATGAC': 'ID', 'ATGAG': 'I0', 'ATGAT': 'IE',
               'ATGCA': 'IF', 'ATGCC': 'I1', 'ATGCG': 'IH', 'ATGCT': 'II', 'ATGGA': 'I2', 'ATGGC': 'IJ',
               'ATGGG': 'IK', 'ATGGT': 'I3', 'ATGTA': 'IL', 'ATGTC': 'IM', 'ATGTG': 'I4', 'ATGTT': 'IN',
               'ATTAA': 'IO', 'ATTAC': 'I5', 'ATTAG': 'IP', 'ATTAT': 'IQ', 'ATTCA': 'I6', 'ATTCC': 'IR',
               'ATTCG': 'IS', 'ATTCT': 'I7', 'ATTGA': 'IU', 'ATTGC': 'IV', 'ATTGG': 'I8', 'ATTGT': 'IW',
               'ATTTA': 'IX', 'ATTTC': 'I9', 'ATTTG': 'IY', 'ATTTT': 'IZ', 'CAAAA': '2B', 'CAAAC': '2D',
               'CAAAG': '20', 'CAAAT': '2E', 'CAACA': '2F', 'CAACC': '21', 'CAACG': '2H', 'CAACT': '2I',
               'CAAGA': '22', 'CAAGC': '2J', 'CAAGG': '2K', 'CAAGT': '23', 'CAATA': '2L', 'CAATC': '2M',
               'CAATG': '24', 'CAATT': '2N', 'CACAA': '2O', 'CACAC': '25', 'CACAG': '2P', 'CACAT': '2Q',
               'CACCA': '26', 'CACCC': '2R', 'CACCG': '2S', 'CACCT': '27', 'CACGA': '2U', 'CACGC': '2V',
               'CACGG': '28', 'CACGT': '2W', 'CACTA': '2X', 'CACTC': '29', 'CACTG': '2Y', 'CACTT': '2Z',
               'CAGAA': 'JB', 'CAGAC': 'JD', 'CAGAG': 'J0', 'CAGAT': 'JE', 'CAGCA': 'JF', 'CAGCC': 'J1',
               'CAGCG': 'JH', 'CAGCT': 'JI', 'CAGGA': 'J2', 'CAGGC': 'JJ', 'CAGGG': 'JK', 'CAGGT': 'J3',
               'CAGTA': 'JL', 'CAGTC': 'JM', 'CAGTG': 'J4', 'CAGTT': 'JN', 'CATAA': 'JO', 'CATAC': 'J5',
               'CATAG': 'JP', 'CATAT': 'JQ', 'CATCA': 'J6', 'CATCC': 'JR', 'CATCG': 'JS', 'CATCT': 'J7',
               'CATGA': 'JU', 'CATGC': 'JV', 'CATGG': 'J8', 'CATGT': 'JW', 'CATTA': 'JX', 'CATTC': 'J9',
               'CATTG': 'JY', 'CATTT': 'JZ', 'CCAAA': 'KB', 'CCAAC': 'KD', 'CCAAG': 'K0', 'CCAAT': 'KE',
               'CCACA': 'KF', 'CCACC': 'K1', 'CCACG': 'KH', 'CCACT': 'KI', 'CCAGA': 'K2', 'CCAGC': 'KJ',
               'CCAGG': 'KK', 'CCAGT': 'K3', 'CCATA': 'KL', 'CCATC': 'KM', 'CCATG': 'K4', 'CCATT': 'KN',
               'CCCAA': 'KO', 'CCCAC': 'K5', 'CCCAG': 'KP', 'CCCAT': 'KQ', 'CCCCA': 'K6', 'CCCCC': 'KR',
               'CCCCG': 'KS', 'CCCCT': 'K7', 'CCCGA': 'KU', 'CCCGC': 'KV', 'CCCGG': 'K8', 'CCCGT': 'KW',
               'CCCTA': 'KX', 'CCCTC': 'K9', 'CCCTG': 'KY', 'CCCTT': 'KZ', 'CCGAA': '3B', 'CCGAC': '3D',
               'CCGAG': '30', 'CCGAT': '3E', 'CCGCA': '3F', 'CCGCC': '31', 'CCGCG': '3H', 'CCGCT': '3I',
               'CCGGA': '32', 'CCGGC': '3J', 'CCGGG': '3K', 'CCGGT': '33', 'CCGTA': '3L', 'CCGTC': '3M',
               'CCGTG': '34', 'CCGTT': '3N', 'CCTAA': '3O', 'CCTAC': '35', 'CCTAG': '3P', 'CCTAT': '3Q',
               'CCTCA': '36', 'CCTCC': '3R', 'CCTCG': '3S', 'CCTCT': '37', 'CCTGA': '3U', 'CCTGC': '3V',
               'CCTGG': '38', 'CCTGT': '3W', 'CCTTA': '3X', 'CCTTC': '39', 'CCTTG': '3Y', 'CCTTT': '3Z',
               'CGAAA': 'LB', 'CGAAC': 'LD', 'CGAAG': 'L0', 'CGAAT': 'LE', 'CGACA': 'LF', 'CGACC': 'L1',
               'CGACG': 'LH', 'CGACT': 'LI', 'CGAGA': 'L2', 'CGAGC': 'LJ', 'CGAGG': 'LK', 'CGAGT': 'L3',
               'CGATA': 'LL', 'CGATC': 'LM', 'CGATG': 'L4', 'CGATT': 'LN', 'CGCAA': 'LO', 'CGCAC': 'L5',
               'CGCAG': 'LP', 'CGCAT': 'LQ', 'CGCCA': 'L6', 'CGCCC': 'LR', 'CGCCG': 'LS', 'CGCCT': 'L7',
               'CGCGA': 'LU', 'CGCGC': 'LV', 'CGCGG': 'L8', 'CGCGT': 'LW', 'CGCTA': 'LX', 'CGCTC': 'L9',
               'CGCTG': 'LY', 'CGCTT': 'LZ', 'CGGAA': 'MB', 'CGGAC': 'MD', 'CGGAG': 'M0', 'CGGAT': 'ME',
               'CGGCA': 'MF', 'CGGCC': 'M1', 'CGGCG': 'MH', 'CGGCT': 'MI', 'CGGGA': 'M2', 'CGGGC': 'MJ',
               'CGGGG': 'MK', 'CGGGT': 'M3', 'CGGTA': 'ML', 'CGGTC': 'MM', 'CGGTG': 'M4', 'CGGTT': 'MN',
               'CGTAA': 'MO', 'CGTAC': 'M5', 'CGTAG': 'MP', 'CGTAT': 'MQ', 'CGTCA': 'M6', 'CGTCC': 'MR',
               'CGTCG': 'MS', 'CGTCT': 'M7', 'CGTGA': 'MU', 'CGTGC': 'MV', 'CGTGG': 'M8', 'CGTGT': 'MW',
               'CGTTA': 'MX', 'CGTTC': 'M9', 'CGTTG': 'MY', 'CGTTT': 'MZ', 'CTAAA': '4B', 'CTAAC': '4D',
               'CTAAG': '40', 'CTAAT': '4E', 'CTACA': '4F', 'CTACC': '41', 'CTACG': '4H', 'CTACT': '4I',
               'CTAGA': '42', 'CTAGC': '4J', 'CTAGG': '4K', 'CTAGT': '43', 'CTATA': '4L', 'CTATC': '4M',
               'CTATG': '44', 'CTATT': '4N', 'CTCAA': '4O', 'CTCAC': '45', 'CTCAG': '4P', 'CTCAT': '4Q',
               'CTCCA': '46', 'CTCCC': '4R', 'CTCCG': '4S', 'CTCCT': '47', 'CTCGA': '4U', 'CTCGC': '4V',
               'CTCGG': '48', 'CTCGT': '4W', 'CTCTA': '4X', 'CTCTC': '49', 'CTCTG': '4Y', 'CTCTT': '4Z',
               'CTGAA': 'NB', 'CTGAC': 'ND', 'CTGAG': 'N0', 'CTGAT': 'NE', 'CTGCA': 'NF', 'CTGCC': 'N1',
               'CTGCG': 'NH', 'CTGCT': 'NI', 'CTGGA': 'N2', 'CTGGC': 'NJ', 'CTGGG': 'NK', 'CTGGT': 'N3',
               'CTGTA': 'NL', 'CTGTC': 'NM', 'CTGTG': 'N4', 'CTGTT': 'NN', 'CTTAA': 'NO', 'CTTAC': 'N5',
               'CTTAG': 'NP', 'CTTAT': 'NQ', 'CTTCA': 'N6', 'CTTCC': 'NR', 'CTTCG': 'NS', 'CTTCT': 'N7',
               'CTTGA': 'NU', 'CTTGC': 'NV', 'CTTGG': 'N8', 'CTTGT': 'NW', 'CTTTA': 'NX', 'CTTTC': 'N9',
               'CTTTG': 'NY', 'CTTTT': 'NZ', 'GAAAA': 'OB', 'GAAAC': 'OD', 'GAAAG': 'O0', 'GAAAT': 'OE',
               'GAACA': 'OF', 'GAACC': 'O1', 'GAACG': 'OH', 'GAACT': 'OI', 'GAAGA': 'O2', 'GAAGC': 'OJ',
               'GAAGG': 'OK', 'GAAGT': 'O3', 'GAATA': 'OL', 'GAATC': 'OM', 'GAATG': 'O4', 'GAATT': 'ON',
               'GACAA': 'OO', 'GACAC': 'O5', 'GACAG': 'OP', 'GACAT': 'OQ', 'GACCA': 'O6', 'GACCC': 'OR',
               'GACCG': 'OS', 'GACCT': 'O7', 'GACGA': 'OU', 'GACGC': 'OV', 'GACGG': 'O8', 'GACGT': 'OW',
               'GACTA': 'OX', 'GACTC': 'O9', 'GACTG': 'OY', 'GACTT': 'OZ', 'GAGAA': '5B', 'GAGAC': '5D',
               'GAGAG': '50', 'GAGAT': '5E', 'GAGCA': '5F', 'GAGCC': '51', 'GAGCG': '5H', 'GAGCT': '5I',
               'GAGGA': '52', 'GAGGC': '5J', 'GAGGG': '5K', 'GAGGT': '53', 'GAGTA': '5L', 'GAGTC': '5M',
               'GAGTG': '54', 'GAGTT': '5N', 'GATAA': '5O', 'GATAC': '55', 'GATAG': '5P', 'GATAT': '5Q',
               'GATCA': '56', 'GATCC': '5R', 'GATCG': '5S', 'GATCT': '57', 'GATGA': '5U', 'GATGC': '5V',
               'GATGG': '58', 'GATGT': '5W', 'GATTA': '5X', 'GATTC': '59', 'GATTG': '5Y', 'GATTT': '5Z',
               'GCAAA': 'PB', 'GCAAC': 'PD', 'GCAAG': 'P0', 'GCAAT': 'PE', 'GCACA': 'PF', 'GCACC': 'P1',
               'GCACG': 'PH', 'GCACT': 'PI', 'GCAGA': 'P2', 'GCAGC': 'PJ', 'GCAGG': 'PK', 'GCAGT': 'P3',
               'GCATA': 'PL', 'GCATC': 'PM', 'GCATG': 'P4', 'GCATT': 'PN', 'GCCAA': 'PO', 'GCCAC': 'P5',
               'GCCAG': 'PP', 'GCCAT': 'PQ', 'GCCCA': 'P6', 'GCCCC': 'PR', 'GCCCG': 'PS', 'GCCCT': 'P7',
               'GCCGA': 'PU', 'GCCGC': 'PV', 'GCCGG': 'P8', 'GCCGT': 'PW', 'GCCTA': 'PX', 'GCCTC': 'P9',
               'GCCTG': 'PY', 'GCCTT': 'PZ', 'GCGAA': 'QB', 'GCGAC': 'QD', 'GCGAG': 'Q0', 'GCGAT': 'QE',
               'GCGCA': 'QF', 'GCGCC': 'Q1', 'GCGCG': 'QH', 'GCGCT': 'QI', 'GCGGA': 'Q2', 'GCGGC': 'QJ',
               'GCGGG': 'QK', 'GCGGT': 'Q3', 'GCGTA': 'QL', 'GCGTC': 'QM', 'GCGTG': 'Q4', 'GCGTT': 'QN',
               'GCTAA': 'QO', 'GCTAC': 'Q5', 'GCTAG': 'QP', 'GCTAT': 'QQ', 'GCTCA': 'Q6', 'GCTCC': 'QR',
               'GCTCG': 'QS', 'GCTCT': 'Q7', 'GCTGA': 'QU', 'GCTGC': 'QV', 'GCTGG': 'Q8', 'GCTGT': 'QW',
               'GCTTA': 'QX', 'GCTTC': 'Q9', 'GCTTG': 'QY', 'GCTTT': 'QZ', 'GGAAA': '6B', 'GGAAC': '6D',
               'GGAAG': '60', 'GGAAT': '6E', 'GGACA': '6F', 'GGACC': '61', 'GGACG': '6H', 'GGACT': '6I',
               'GGAGA': '62', 'GGAGC': '6J', 'GGAGG': '6K', 'GGAGT': '63', 'GGATA': '6L', 'GGATC': '6M',
               'GGATG': '64', 'GGATT': '6N', 'GGCAA': '6O', 'GGCAC': '65', 'GGCAG': '6P', 'GGCAT': '6Q',
               'GGCCA': '66', 'GGCCC': '6R', 'GGCCG': '6S', 'GGCCT': '67', 'GGCGA': '6U', 'GGCGC': '6V',
               'GGCGG': '68', 'GGCGT': '6W', 'GGCTA': '6X', 'GGCTC': '69', 'GGCTG': '6Y', 'GGCTT': '6Z',
               'GGGAA': 'RB', 'GGGAC': 'RD', 'GGGAG': 'R0', 'GGGAT': 'RE', 'GGGCA': 'RF', 'GGGCC': 'R1',
               'GGGCG': 'RH', 'GGGCT': 'RI', 'GGGGA': 'R2', 'GGGGC': 'RJ', 'GGGGG': 'RK', 'GGGGT': 'R3',
               'GGGTA': 'RL', 'GGGTC': 'RM', 'GGGTG': 'R4', 'GGGTT': 'RN', 'GGTAA': 'RO', 'GGTAC': 'R5',
               'GGTAG': 'RP', 'GGTAT': 'RQ', 'GGTCA': 'R6', 'GGTCC': 'RR', 'GGTCG': 'RS', 'GGTCT': 'R7',
               'GGTGA': 'RU', 'GGTGC': 'RV', 'GGTGG': 'R8', 'GGTGT': 'RW', 'GGTTA': 'RX', 'GGTTC': 'R9',
               'GGTTG': 'RY', 'GGTTT': 'RZ', 'GTAAA': 'SB', 'GTAAC': 'SD', 'GTAAG': 'S0', 'GTAAT': 'SE',
               'GTACA': 'SF', 'GTACC': 'S1', 'GTACG': 'SH', 'GTACT': 'SI', 'GTAGA': 'S2', 'GTAGC': 'SJ',
               'GTAGG': 'SK', 'GTAGT': 'S3', 'GTATA': 'SL', 'GTATC': 'SM', 'GTATG': 'S4', 'GTATT': 'SN',
               'GTCAA': 'SO', 'GTCAC': 'S5', 'GTCAG': 'SP', 'GTCAT': 'SQ', 'GTCCA': 'S6', 'GTCCC': 'SR',
               'GTCCG': 'SS', 'GTCCT': 'S7', 'GTCGA': 'SU', 'GTCGC': 'SV', 'GTCGG': 'S8', 'GTCGT': 'SW',
               'GTCTA': 'SX', 'GTCTC': 'S9', 'GTCTG': 'SY', 'GTCTT': 'SZ', 'GTGAA': '7B', 'GTGAC': '7D',
               'GTGAG': '70', 'GTGAT': '7E', 'GTGCA': '7F', 'GTGCC': '71', 'GTGCG': '7H', 'GTGCT': '7I',
               'GTGGA': '72', 'GTGGC': '7J', 'GTGGG': '7K', 'GTGGT': '73', 'GTGTA': '7L', 'GTGTC': '7M',
               'GTGTG': '74', 'GTGTT': '7N', 'GTTAA': '7O', 'GTTAC': '75', 'GTTAG': '7P', 'GTTAT': '7Q',
               'GTTCA': '76', 'GTTCC': '7R', 'GTTCG': '7S', 'GTTCT': '77', 'GTTGA': '7U', 'GTTGC': '7V',
               'GTTGG': '78', 'GTTGT': '7W', 'GTTTA': '7X', 'GTTTC': '79', 'GTTTG': '7Y', 'GTTTT': '7Z',
               'TAAAA': 'UB', 'TAAAC': 'UD', 'TAAAG': 'U0', 'TAAAT': 'UE', 'TAACA': 'UF', 'TAACC': 'U1',
               'TAACG': 'UH', 'TAACT': 'UI', 'TAAGA': 'U2', 'TAAGC': 'UJ', 'TAAGG': 'UK', 'TAAGT': 'U3',
               'TAATA': 'UL', 'TAATC': 'UM', 'TAATG': 'U4', 'TAATT': 'UN', 'TACAA': 'UO', 'TACAC': 'U5',
               'TACAG': 'UP', 'TACAT': 'UQ', 'TACCA': 'U6', 'TACCC': 'UR', 'TACCG': 'US', 'TACCT': 'U7',
               'TACGA': 'UU', 'TACGC': 'UV', 'TACGG': 'U8', 'TACGT': 'UW', 'TACTA': 'UX', 'TACTC': 'U9',
               'TACTG': 'UY', 'TACTT': 'UZ', 'TAGAA': 'VB', 'TAGAC': 'VD', 'TAGAG': 'V0', 'TAGAT': 'VE',
               'TAGCA': 'VF', 'TAGCC': 'V1', 'TAGCG': 'VH', 'TAGCT': 'VI', 'TAGGA': 'V2', 'TAGGC': 'VJ',
               'TAGGG': 'VK', 'TAGGT': 'V3', 'TAGTA': 'VL', 'TAGTC': 'VM', 'TAGTG': 'V4', 'TAGTT': 'VN',
               'TATAA': 'VO', 'TATAC': 'V5', 'TATAG': 'VP', 'TATAT': 'VQ', 'TATCA': 'V6', 'TATCC': 'VR',
               'TATCG': 'VS', 'TATCT': 'V7', 'TATGA': 'VU', 'TATGC': 'VV', 'TATGG': 'V8', 'TATGT': 'VW',
               'TATTA': 'VX', 'TATTC': 'V9', 'TATTG': 'VY', 'TATTT': 'VZ', 'TCAAA': '8B', 'TCAAC': '8D',
               'TCAAG': '80', 'TCAAT': '8E', 'TCACA': '8F', 'TCACC': '81', 'TCACG': '8H', 'TCACT': '8I',
               'TCAGA': '82', 'TCAGC': '8J', 'TCAGG': '8K', 'TCAGT': '83', 'TCATA': '8L', 'TCATC': '8M',
               'TCATG': '84', 'TCATT': '8N', 'TCCAA': '8O', 'TCCAC': '85', 'TCCAG': '8P', 'TCCAT': '8Q',
               'TCCCA': '86', 'TCCCC': '8R', 'TCCCG': '8S', 'TCCCT': '87', 'TCCGA': '8U', 'TCCGC': '8V',
               'TCCGG': '88', 'TCCGT': '8W', 'TCCTA': '8X', 'TCCTC': '89', 'TCCTG': '8Y', 'TCCTT': '8Z',
               'TCGAA': 'WB', 'TCGAC': 'WD', 'TCGAG': 'W0', 'TCGAT': 'WE', 'TCGCA': 'WF', 'TCGCC': 'W1',
               'TCGCG': 'WH', 'TCGCT': 'WI', 'TCGGA': 'W2', 'TCGGC': 'WJ', 'TCGGG': 'WK', 'TCGGT': 'W3',
               'TCGTA': 'WL', 'TCGTC': 'WM', 'TCGTG': 'W4', 'TCGTT': 'WN', 'TCTAA': 'WO', 'TCTAC': 'W5',
               'TCTAG': 'WP', 'TCTAT': 'WQ', 'TCTCA': 'W6', 'TCTCC': 'WR', 'TCTCG': 'WS', 'TCTCT': 'W7',
               'TCTGA': 'WU', 'TCTGC': 'WV', 'TCTGG': 'W8', 'TCTGT': 'WW', 'TCTTA': 'WX', 'TCTTC': 'W9',
               'TCTTG': 'WY', 'TCTTT': 'WZ', 'TGAAA': 'XB', 'TGAAC': 'XD', 'TGAAG': 'X0', 'TGAAT': 'XE',
               'TGACA': 'XF', 'TGACC': 'X1', 'TGACG': 'XH', 'TGACT': 'XI', 'TGAGA': 'X2', 'TGAGC': 'XJ',
               'TGAGG': 'XK', 'TGAGT': 'X3', 'TGATA': 'XL', 'TGATC': 'XM', 'TGATG': 'X4', 'TGATT': 'XN',
               'TGCAA': 'XO', 'TGCAC': 'X5', 'TGCAG': 'XP', 'TGCAT': 'XQ', 'TGCCA': 'X6', 'TGCCC': 'XR',
               'TGCCG': 'XS', 'TGCCT': 'X7', 'TGCGA': 'XU', 'TGCGC': 'XV', 'TGCGG': 'X8', 'TGCGT': 'XW',
               'TGCTA': 'XX', 'TGCTC': 'X9', 'TGCTG': 'XY', 'TGCTT': 'XZ', 'TGGAA': '9B', 'TGGAC': '9D',
               'TGGAG': '90', 'TGGAT': '9E', 'TGGCA': '9F', 'TGGCC': '91', 'TGGCG': '9H', 'TGGCT': '9I',
               'TGGGA': '92', 'TGGGC': '9J', 'TGGGG': '9K', 'TGGGT': '93', 'TGGTA': '9L', 'TGGTC': '9M',
               'TGGTG': '94', 'TGGTT': '9N', 'TGTAA': '9O', 'TGTAC': '95', 'TGTAG': '9P', 'TGTAT': '9Q',
               'TGTCA': '96', 'TGTCC': '9R', 'TGTCG': '9S', 'TGTCT': '97', 'TGTGA': '9U', 'TGTGC': '9V',
               'TGTGG': '98', 'TGTGT': '9W', 'TGTTA': '9X', 'TGTTC': '99', 'TGTTG': '9Y', 'TGTTT': '9Z',
               'TTAAA': 'YB', 'TTAAC': 'YD', 'TTAAG': 'Y0', 'TTAAT': 'YE', 'TTACA': 'YF', 'TTACC': 'Y1',
               'TTACG': 'YH', 'TTACT': 'YI', 'TTAGA': 'Y2', 'TTAGC': 'YJ', 'TTAGG': 'YK', 'TTAGT': 'Y3',
               'TTATA': 'YL', 'TTATC': 'YM', 'TTATG': 'Y4', 'TTATT': 'YN', 'TTCAA': 'YO', 'TTCAC': 'Y5',
               'TTCAG': 'YP', 'TTCAT': 'YQ', 'TTCCA': 'Y6', 'TTCCC': 'YR', 'TTCCG': 'YS', 'TTCCT': 'Y7',
               'TTCGA': 'YU', 'TTCGC': 'YV', 'TTCGG': 'Y8', 'TTCGT': 'YW', 'TTCTA': 'YX', 'TTCTC': 'Y9',
               'TTCTG': 'YY', 'TTCTT': 'YZ', 'TTGAA': 'ZB', 'TTGAC': 'ZD', 'TTGAG': 'Z0', 'TTGAT': 'ZE',
               'TTGCA': 'ZF', 'TTGCC': 'Z1', 'TTGCG': 'ZH', 'TTGCT': 'ZI', 'TTGGA': 'Z2', 'TTGGC': 'ZJ',
               'TTGGG': 'ZK', 'TTGGT': 'Z3', 'TTGTA': 'ZL', 'TTGTC': 'ZM', 'TTGTG': 'Z4', 'TTGTT': 'ZN',
               'TTTAA': 'ZO', 'TTTAC': 'Z5', 'TTTAG': 'ZP', 'TTTAT': 'ZQ', 'TTTCA': 'Z6', 'TTTCC': 'ZR',
               'TTTCG': 'ZS', 'TTTCT': 'Z7', 'TTTGA': 'ZU', 'TTTGC': 'ZV', 'TTTGG': 'Z8', 'TTTGT': 'ZW',
               'TTTTA': 'ZX', 'TTTTC': 'Z9', 'TTTTG': 'ZY', 'TTTTT': 'ZZ', 'A': 'B', 'C': 'D',
               'G': '0', 'T': 'E', 'AA': 'F', 'AC': '1', 'AG': 'H', 'AT': 'I',
               'CA': '2', 'CC': 'J', 'CG': 'K', 'CT': '3', 'GA': 'L', 'GC': 'M',
               'GG': '4', 'GT': 'N', 'TA': 'O', 'TC': '5', 'TG': 'P', 'TT': 'Q',
               'AAA': '6', 'AAC': 'R', 'AAG': 'S', 'AAT': '7', 'ACA': 'U', 'ACC': 'V',
               'ACG': '8', 'ACT': 'W', 'AGA': 'X', 'AGC': '9', 'AGG': 'Y', 'AGT': 'Z',
               'ATA': 'DB', 'ATC': 'DD', 'ATG': 'D0', 'ATT': 'DE', 'CAA': 'DF', 'CAC': 'D1',
               'CAG': 'DH', 'CAT': 'DI', 'CCA': 'D2', 'CCC': 'DJ', 'CCG': 'DK', 'CCT': 'D3',
               'CGA': 'DL', 'CGC': 'DM', 'CGG': 'D4', 'CGT': 'DN', 'CTA': 'DO', 'CTC': 'D5',
               'CTG': 'DP', 'CTT': 'DQ', 'GAA': 'D6', 'GAC': 'DR', 'GAG': 'DS', 'GAT': 'D7',
               'GCA': 'DU', 'GCC': 'DV', 'GCG': 'D8', 'GCT': 'DW', 'GGA': 'DX', 'GGC': 'D9',
               'GGG': 'DY', 'GGT': 'DZ', 'GTA': '0B', 'GTC': '0D', 'GTG': '00', 'GTT': '0E',
               'TAA': '0F', 'TAC': '01', 'TAG': '0H', 'TAT': '0I', 'TCA': '02', 'TCC': '0J',
               'TCG': '0K', 'TCT': '03', 'TGA': '0L', 'TGC': '0M', 'TGG': '04', 'TGT': '0N',
               'TTA': '0O', 'TTC': '05', 'TTG': '0P', 'TTT': '0Q', 'AAAA': '06', 'AAAC': '0R',
               'AAAG': '0S', 'AAAT': '07', 'AACA': '0U', 'AACC': '0V', 'AACG': '08', 'AACT': '0W',
               'AAGA': '0X', 'AAGC': '09', 'AAGG': '0Y', 'AAGT': '0Z', 'AATA': 'EB', 'AATC': 'ED',
               'AATG': 'E0', 'AATT': 'EE', 'ACAA': 'EF', 'ACAC': 'E1', 'ACAG': 'EH', 'ACAT': 'EI',
               'ACCA': 'E2', 'ACCC': 'EJ', 'ACCG': 'EK', 'ACCT': 'E3', 'ACGA': 'EL', 'ACGC': 'EM',
               'ACGG': 'E4', 'ACGT': 'EN', 'ACTA': 'EO', 'ACTC': 'E5', 'ACTG': 'EP', 'ACTT': 'EQ',
               'AGAA': 'E6', 'AGAC': 'ER', 'AGAG': 'ES', 'AGAT': 'E7', 'AGCA': 'EU', 'AGCC': 'EV',
               'AGCG': 'E8', 'AGCT': 'EW', 'AGGA': 'EX', 'AGGC': 'E9', 'AGGG': 'EY', 'AGGT': 'EZ',
               'AGTA': 'FB', 'AGTC': 'FD', 'AGTG': 'F0', 'AGTT': 'FE', 'ATAA': 'FF', 'ATAC': 'F1',
               'ATAG': 'FH', 'ATAT': 'FI', 'ATCA': 'F2', 'ATCC': 'FJ', 'ATCG': 'FK', 'ATCT': 'F3',
               'ATGA': 'FL', 'ATGC': 'FM', 'ATGG': 'F4', 'ATGT': 'FN', 'ATTA': 'FO', 'ATTC': 'F5',
               'ATTG': 'FP', 'ATTT': 'FQ', 'CAAA': 'F6', 'CAAC': 'FR', 'CAAG': 'FS', 'CAAT': 'F7',
               'CACA': 'FU', 'CACC': 'FV', 'CACG': 'F8', 'CACT': 'FW', 'CAGA': 'FX', 'CAGC': 'F9',
               'CAGG': 'FY', 'CAGT': 'FZ', 'CATA': '1B', 'CATC': '1D', 'CATG': '10', 'CATT': '1E',
               'CCAA': '1F', 'CCAC': '11', 'CCAG': '1H', 'CCAT': '1I', 'CCCA': '12', 'CCCC': '1J',
               'CCCG': '1K', 'CCCT': '13', 'CCGA': '1L', 'CCGC': '1M', 'CCGG': '14', 'CCGT': '1N',
               'CCTA': '1O', 'CCTC': '15', 'CCTG': '1P', 'CCTT': '1Q', 'CGAA': '16', 'CGAC': '1R',
               'CGAG': '1S', 'CGAT': '17', 'CGCA': '1U', 'CGCC': '1V', 'CGCG': '18', 'CGCT': '1W',
               'CGGA': '1X', 'CGGC': '19', 'CGGG': '1Y', 'CGGT': '1Z', 'CGTA': 'HB', 'CGTC': 'HD',
               'CGTG': 'H0', 'CGTT': 'HE', 'CTAA': 'HF', 'CTAC': 'H1', 'CTAG': 'HH', 'CTAT': 'HI',
               'CTCA': 'H2', 'CTCC': 'HJ', 'CTCG': 'HK', 'CTCT': 'H3', 'CTGA': 'HL', 'CTGC': 'HM',
               'CTGG': 'H4', 'CTGT': 'HN', 'CTTA': 'HO', 'CTTC': 'H5', 'CTTG': 'HP', 'CTTT': 'HQ',
               'GAAA': 'H6', 'GAAC': 'HR', 'GAAG': 'HS', 'GAAT': 'H7', 'GACA': 'HU', 'GACC': 'HV',
               'GACG': 'H8', 'GACT': 'HW', 'GAGA': 'HX', 'GAGC': 'H9', 'GAGG': 'HY', 'GAGT': 'HZ',
               'GATA': 'IB', 'GATC': 'ID', 'GATG': 'I0', 'GATT': 'IE', 'GCAA': 'IF', 'GCAC': 'I1',
               'GCAG': 'IH', 'GCAT': 'II', 'GCCA': 'I2', 'GCCC': 'IJ', 'GCCG': 'IK', 'GCCT': 'I3',
               'GCGA': 'IL', 'GCGC': 'IM', 'GCGG': 'I4', 'GCGT': 'IN', 'GCTA': 'IO', 'GCTC': 'I5',
               'GCTG': 'IP', 'GCTT': 'IQ', 'GGAA': 'I6', 'GGAC': 'IR', 'GGAG': 'IS', 'GGAT': 'I7',
               'GGCA': 'IU', 'GGCC': 'IV', 'GGCG': 'I8', 'GGCT': 'IW', 'GGGA': 'IX', 'GGGC': 'I9',
               'GGGG': 'IY', 'GGGT': 'IZ', 'GGTA': '2B', 'GGTC': '2D', 'GGTG': '20', 'GGTT': '2E',
               'GTAA': '2F', 'GTAC': '21', 'GTAG': '2H', 'GTAT': '2I', 'GTCA': '22', 'GTCC': '2J',
               'GTCG': '2K', 'GTCT': '23', 'GTGA': '2L', 'GTGC': '2M', 'GTGG': '24', 'GTGT': '2N',
               'GTTA': '2O', 'GTTC': '25', 'GTTG': '2P', 'GTTT': '2Q', 'TAAA': '26', 'TAAC': '2R',
               'TAAG': '2S', 'TAAT': '27', 'TACA': '2U', 'TACC': '2V', 'TACG': '28', 'TACT': '2W',
               'TAGA': '2X', 'TAGC': '29', 'TAGG': '2Y', 'TAGT': '2Z', 'TATA': 'JB', 'TATC': 'JD',
               'TATG': 'J0', 'TATT': 'JE', 'TCAA': 'JF', 'TCAC': 'J1', 'TCAG': 'JH', 'TCAT': 'JI',
               'TCCA': 'J2', 'TCCC': 'JJ', 'TCCG': 'JK', 'TCCT': 'J3', 'TCGA': 'JL', 'TCGC': 'JM',
               'TCGG': 'J4', 'TCGT': 'JN', 'TCTA': 'JO', 'TCTC': 'J5', 'TCTG': 'JP', 'TCTT': 'JQ',
               'TGAA': 'J6', 'TGAC': 'JR', 'TGAG': 'JS', 'TGAT': 'J7', 'TGCA': 'JU', 'TGCC': 'JV',
               'TGCG': 'J8', 'TGCT': 'JW', 'TGGA': 'JX', 'TGGC': 'J9', 'TGGG': 'JY', 'TGGT': 'JZ',
               'TGTA': 'KB', 'TGTC': 'KD', 'TGTG': 'K0', 'TGTT': 'KE', 'TTAA': 'KF', 'TTAC': 'K1',
               'TTAG': 'KH', 'TTAT': 'KI', 'TTCA': 'K2', 'TTCC': 'KJ', 'TTCG': 'KK', 'TTCT': 'K3',
               'TTGA': 'KL', 'TTGC': 'KM', 'TTGG': 'K4', 'TTGT': 'KN', 'TTTA': 'KO', 'TTTC': 'K5',
               'TTTG': 'KP', 'TTTT': 'KQ'}

def encode_sequence(sequence, prefix):
    """
    Encodes the sequence into its corresponding license plate with given prefix (if given one)
    :param sequence: The sequence being encoded
    :param prefix: The prefix to use for the license plate
    :return: The license plate it encodes to
    """
    length = len(sequence)
    # Encode label
    if prefix is '':
        final_result = [(str(length) + '-')]
    else:
        final_result = [prefix + "-" + str(length) + "-"]

    work_sequence = sequence
    while work_sequence != '':
        try:
            final_result.append(encode_hash[work_sequence[0:5]])
        except KeyError as err:
            if not err.args:
                err.args = ('',)
            err.args = err.args + ("Error, exiting: Segment '" + work_sequence[0:5] +
                                   "' from sequence '" + sequence + "' is invalid.",)
            raise
        work_sequence = work_sequence[5:]

    return ''.join(final_result)

def convert(seq, encode, prefix):
    if prefix is None:
        prefix = ''
    else:
        if '-' in prefix or ' ' in prefix:
            sys.stderr.write("Warning: Dashes and spaces are not permitted in the license plate prefix."
                             "Program will remove all instances automatically from prefix '" + prefix + "'.\n")
            prefix = prefix.replace('-', '')
            prefix = prefix.replace(' ', '')

    if encode:
        # Encode
        if seq != '':
            cleaned = seq.upper().replace('U', 'T')
            return encode_sequence(cleaned, prefix)
            # if is_sequence(cleaned):
            #     
            # else:
            #     raise KeyError('Error, exiting: Illegal characters in line "' + seq + '"')
    # else:
    #     # Decode
    #     if seq != '':
    #         cleaned = seq.upper()
    #         return decode_sequence(cleaned)

def make_id(seq):
    """
    Create a unique identifier for the sequence from the nucleotides,
    replacing 5 nts for a unique sequence.

    It uses the code from *mirtop.mirna.keys()*.

    Inspired by MINTplate: https://cm.jefferson.edu/MINTbase
    https://github.com/TJU-CMC-Org/MINTmap/tree/master/MINTplates

    Args:
        *seq(str)*: nucleotides sequences.

    Returns:
        *idName(str)*: unique identifier for the sequence.
    """
    try:
        idu = convert(seq, True, 'iso')
    except KeyError:
        print("WARNING : miRGFF sequence issue") ## todo:log
    return idu

##########################################
## GFF DEFINITIONS (todo: ugly, clean & comment)
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
                parents_list = a.get_tag('PA')
                parent = "Parent={}".format(",".join(parents_list.split('/')))
                name = "Name={}".format(a.reference_name)
                alias = "Alias={}".format(",".join(list(set([i.split('_')[0] for i in d_OptimiR[seqid].ident]))))
                expression = "Expression={}".format(float(a.get_tag('XC')))
                expression_optimir = "expression_OptimiR={}".format(float(a.get_tag('XW')) * float(a.get_tag('XC')))
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
                line = "\t".join([seqid, source, typ, str(start), str(end), score, strand, phase, attributes]) + ";"
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

