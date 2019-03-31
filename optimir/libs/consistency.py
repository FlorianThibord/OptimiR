#-*- coding: utf-8 -*-
# Florian THIBORD (22/08/17)
########################################################
#                  OptimiR PIPELINE                    # 
########################################################

# Standard libraries
from pysam import AlignmentFile

# Personal libraries
from optimir.libs.essentials import *

## Number of inconsistent reads aligned on a polymiR above wich the alignment is flagged as suspicious
## This value is paired with INCONSISTENT_RATE_THRESHOLD provided by the user
## This thresholds allow to filter alignments of reads with sequencing errors that could replicate the variant. Alignments that go above these threshold are suspicious as they are unlikely due to sequencing errors replicating the genetic variation.
INCONSISTENT_NB_THRESHOLD = 5

## Return True if soft clip over variant (WARNING: works fine with SNPs, but not with indels...)
## po_in_rna starts at index 0
def is_variant_soft_clipped(alignment, variant, len_rna):
    pos_in_rna = variant.pos_in_rna
    isotype = alignment.get_tag("IT")
    iso5, iso3 = isotype.split('[')[1].split(']')[0].split(',')
    list_pos_soft_clipped = []
    if "-" in iso5 and "+" in iso5:
        nb_trimmed, letters_tail = iso5.split('-')[1].split('+')
        nb_trimmed = int(nb_trimmed)
        nb_tailed = len(letters_tail)
        for i in range(nb_trimmed - nb_tailed, nb_trimmed):
            list_pos_soft_clipped.append(i)
    if "-" in iso3 and "+" in iso3:
        nb_trimmed, letters_tail = iso3.split('-')[1].split('+')
        nb_trimmed = int(nb_trimmed)
        nb_tailed = len(letters_tail)
        for i in range(len_rna - nb_trimmed, len_rna - nb_trimmed + nb_tailed):
            list_pos_soft_clipped.append(i)
    return pos_in_rna in list_pos_soft_clipped

## Check if soft-clip over variant in polymiR
## If only one alignment on polymiR does not have soft clip on variant, discard the others with tag XX and "DISCARDED_SOFT_CLIPPED_VARIANT"
def check_no_soft_clip_on_variant(bam_dict, sample_name, d_OptimiR):
    polymiRs = [key for key, obj in d_OptimiR.items() if len(obj.variants) > 0]
    poly_dict = {} ## {read_name/polymiR_reference_name} -> {sc:[alignments with soft clip], others:[other alignments]}
    for read_name, alignments in bam_dict.items():
        if len(alignments) > 1: # If alone, then it did not align on another version of the polymiR, and is the most accurate mapping
            for alignment in alignments:
                reference_name = alignment.reference_name
                if reference_name in polymiRs: # it's a polymiR
                    OptimiR_obj = d_OptimiR[reference_name]
                    original_polymiR_name = OptimiR_obj.name
                    key_dict = read_name + "/" + original_polymiR_name
                    ## check if soft-clip on at least 1 variants
                    is_soft_clipped = False
                    for variant in OptimiR_obj.variants:
                        if is_variant_soft_clipped(alignment, variant, len(OptimiR_obj.sequence)):
                            is_soft_clipped = True
                    if key_dict not in poly_dict:
                        poly_dict[key_dict] = {"sc":[], "others":[]}
                    if is_soft_clipped:
                        poly_dict[key_dict]["sc"].append(reference_name)
                    else:
                        poly_dict[key_dict]["others"].append(reference_name)
    ## if only one alignment on polymiR without soft-clippin on variant, discard the soft-clipped
    for read_name, alignments in bam_dict.items():
        for alignment in alignments:
            reference_name = alignment.reference_name
            if reference_name in polymiRs: # it's a polymiR
                OptimiR_obj = d_OptimiR[reference_name]
                original_polymiR_name = OptimiR_obj.name
                key_dict = read_name + "/" + original_polymiR_name
                if key_dict in poly_dict:
                    align_dict = poly_dict[key_dict]
                    if (len(align_dict["others"]) == 1) and (reference_name in align_dict["sc"]):
                        alignment.set_tag("XX", "DISCARDED_SC")

# CASE_1 : If alignment on reference sequence, all genotypes must have at least 1 reference allele
# CASE_2 : If alignment on alternative sequence, genotypes of variant in alt sequence must have at least 1 alternative allele, while other variants (if there are other) must have at least 1 reference allele
# CASE Special : polymiR can originate from several hairpins, and sample is homozygous alternative but has reads aligned on reference sequence: do not discard reads on reference as they can originate from other hairpin that does not have the variant
def compute_genotype_consistency(bam_dict, d_OptimiR, sample_name):
    for read_name, alignments in bam_dict.items():
        for alignment in alignments:
            reference_name = alignment.reference_name
            OptimiR_obj = d_OptimiR[reference_name]
            is_reference_miRNA = reference_name == OptimiR_obj.name ## name is always the reference seq
            variants_in_sequence = set(OptimiR_obj.variants)
            has_multi_hairpin = len(OptimiR_obj.hairpins.keys()) > 1
            ## if variant list not empty, then it's a polymiR:
            if len(variants_in_sequence) > 0:
                ## retrieve OptimiR object for original sequence, which has infos for other variants in sequence. In CASE_1 it will be the same object, and the subset other_variants will be empty.
                OptimiR_obj_original = d_OptimiR[reference_name.split('_')[0]]
                variants_in_polymiR = set(OptimiR_obj_original.variants)
                other_variants = variants_in_sequence - variants_in_polymiR
                genotype_dict = {}
                variants_in_polymiR_have_alternative = True
                variants_in_polymiR_have_reference = True
                for variant in variants_in_sequence:
                    rsID = variant.rsID
                    try:
                        genotype = variant.GT_dict[sample_name]
                    except KeyError:
                        genotype = "UNAVAILABLE"
                    genotype_dict[rsID] = genotype
                    if not("1" in genotype):
                        variants_in_polymiR_have_alternative = False
                    if not("0" in genotype):
                        variants_in_polymiR_have_reference = False
                ## set(other_variants) must not contain variants with homozygous alternative
                other_variants_have_reference = True
                for variant in other_variants:
                    rsID = variant.rsID
                    try:
                        genotype = variant.GT_dict[sample_name]
                    except KeyError:
                        genotype = "UNAVAILABLE"
                    genotype_dict[rsID] = genotype
                    if not("0" in genotype):
                        other_variants_have_reference = False
                ## Save genotypes in tag GL
                geno_string = []
                for rsID, genotype in genotype_dict.items():
                    geno_string.append("{}:{}".format(rsID, genotype))
                alignment.set_tag("GL", "/".join(geno_string))
                ## Save consistency in tag GC
                ## CASE_1
                if is_reference_miRNA and variants_in_polymiR_have_reference:
                    alignment.set_tag("GC", "Genotype_OK")
                ## CASE_2
                elif not(is_reference_miRNA) and variants_in_polymiR_have_alternative and other_variants_have_reference:
                    alignment.set_tag("GC", "Genotype_OK")
                ## CASE Special
                elif is_reference_miRNA and has_multi_hairpin:
                    alignment.set_tag("GC", "Genotype_OK_MULTI_HAIRPIN")
                else:
                    alignment.set_tag("GC", "Genotype_INVALID")
                    ## Alignment to remove
                    if not(alignment.has_tag("XX")):
                        alignment.set_tag("XX", "DISCARDED_GENO")

def write_polymiRs_outputs(bam_dict, bam, collapse_table, sample_name, dir_results, d_OptimiR, INCONSISTENT_RATE_THRESHOLD, ANNOT_FILES, WRITE_VCF):
    out_inconsistents_table = dir_results + '/' + "consistency_table." + sample_name + ".annot"
    out_polymiRs_table = dir_results + '/' + "polymiRs_table." + sample_name + ".annot"
    out_inconsistent = dir_results + '/' + "inconsistents." + sample_name + ".sam"
    out_vcf = dir_results + '/' + "geno." + sample_name + ".vcf"
    bam_inconsistents = {}
    polymiRs = [key for key, obj in d_OptimiR.items() if len(obj.polymiR_list) > 0]
    poly_dict_consist = {} ## [mature name] -> (consistents, inconsistents)
    poly_dict_align = {} ## [mature_name] -> {counts on ref, counts on alt}
    for polymiR_name in polymiRs:
        poly_dict_consist[polymiR_name] = (0, 0)
        poly_dict_align[polymiR_name] = (0, 0)
    for read_name, alignments in bam_dict.items():
        reference_names = [a.reference_name for a in alignments if not(a.has_tag("XX"))]
        for alignment in alignments:
            reference_name = alignment.reference_name
            ref_name = reference_name.split('_')[0] ## reference sequence name
            if ref_name in polymiRs:
                if alignment.has_tag("XX"):
                    if alignment.get_tag("XX") == "DISCARDED_GENO":
                        if len(reference_names) == 0: ## If it doesn't map anywhere else
                            ## Save in inconsistents
                            if read_name in bam_inconsistents:
                                bam_inconsistents[read_name].append(alignment)
                            else:
                                bam_inconsistents[read_name] = [alignment]
                            ## Add in poly_dict_consist
                            counts = float(alignment.get_tag("XC"))
                            cons, incons = poly_dict_consist[ref_name]
                            poly_dict_consist[ref_name] = (cons, incons + counts)
                else:
                    counts = float(alignment.get_tag("XC"))
                    weight = float(alignment.get_tag("XW"))
                    cons, incons = poly_dict_consist[ref_name]
                    poly_dict_consist[ref_name] = (cons + (counts * weight), incons)
                    counts_ref, counts_alt = poly_dict_align[ref_name]
                    if (ref_name == reference_name): ## reference sequence
                        poly_dict_align[ref_name] = ((counts_ref + counts * weight), counts_alt)
                    else:
                        poly_dict_align[ref_name] = (counts_ref, (counts_alt + counts * weight))
    ## WRITE INCONSISTENT ALIGNMENTS
    bam_out = AlignmentFile(out_inconsistent, 'w', template=bam)
    bamDict_to_file(bam_inconsistents, bam_out, collapse_table)
    bam_out.close()
    ## WRITE INCONSISTENTS TABLE
    if "i" in ANNOT_FILES:
        with open(out_inconsistents_table, 'w') as out:
            header = "\t".join(["Sample_Name", "polymiR", "rsID_list", "Genotype_list", "Consistent_Counts", "Inconsistent_Counts", "Rate_Inconsistents", "Flag"])
            out.write(header + "\n")
            for ref_name, (cons, incons) in poly_dict_consist.items():
                total = cons + incons
                if total > 0:
                    rate_incons = float(incons) / float(total)
                    list_rsID = ",".join([v.rsID for v in d_OptimiR[ref_name].variants])
                    try:
                        list_genos= ",".join([v.GT_dict[sample_name] for v in d_OptimiR[ref_name].variants])
                    except KeyError:
                        list_genos = "UNAVAILABLE : CANNOT COMPUTE GENOTYPE CONSISTENCY"
                    l = "{}\t{}\t{}\t{}\t{}\t{}\t{}".format(sample_name, ref_name, list_rsID, list_genos, cons, incons, rate_incons)
                    if rate_incons >= INCONSISTENT_RATE_THRESHOLD and incons > INCONSISTENT_NB_THRESHOLD:
                        l += "\tSUSPICIOUS\n"
                    else:
                        l += "\t\n"
                    out.write(l)
    ## WRITE POLYMIRS TABLE
    if "p" in ANNOT_FILES:
        with open(out_polymiRs_table, 'w') as out:
            header = "\t".join(["Sample_Name", "polymiR", "rsID_list", "Genotype_list", "Counts_Reference", "Counts_Alternative", "Rate_Alternative"])
            out.write(header + "\n")
            for ref_name, (counts_ref, counts_alt) in poly_dict_align.items():
                total = counts_ref + counts_alt
                if total > 0:
                    rate_alt = float(counts_alt) / float(total)
                    list_rsID = ",".join([v.rsID for v in d_OptimiR[ref_name].variants])
                    try:
                        list_genos= ",".join([v.GT_dict[sample_name] for v in d_OptimiR[ref_name].variants])
                    except KeyError:
                        list_genos = "UNAVAILABLE"
                    l = "{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(sample_name, ref_name, list_rsID, list_genos, counts_ref, counts_alt, rate_alt)
                    out.write(l)
    if WRITE_VCF:
        ## Todo: move, clean, improve
        with open(out_vcf, 'w') as out:
            header= '##fileformat=VCFv4.1\n##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype inferred">\n##FORMAT=<ID=XC,Number=2,Type=Float,Description="Counts received on allele ref, Counts received on allele alt">\n##INFO=<ID=miR,Number=1,Type=String,Description="mature miRNA (or polymiR) from which genotype has been inferred (using expression of both alleles)">\n##INFO=<ID=Pos_miR,Number=1,Type=Integer,Description="SNP position in mature miR sequence">\n##INFO=<ID=Sens,Number=1,Type=String,Description="DNA Strand from which is transcribed the miRNA">\n##REF,ALT:Reference and alternative alleles as present in the mature miRNA. They might differ from DNA alleles, depending on miRNA sens (see INFOS)\n##Generated with OptimiR: based on the number of reads received by each polymiR allele. The decision is depends on the INCONSISTENT_RATE_THRESHOLD defined by the user (default:0.01, meaning that, to be called, an allele must gather at least 1% of the reads aligned to the polymiR)\n#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	{}\n'.format(sample_name)
            out.write(header)
            for ref_name, (counts_ref, counts_alt) in poly_dict_align.items():
                total = counts_ref + counts_alt
                if total > 0:
                    rate_alt = float(counts_alt) / float(total)
                    variant_list = d_OptimiR[ref_name].variants
                    sens = d_OptimiR[ref_name].coordinates.sens
                    for v in variant_list:
                        infos = "miR={};Pos_miR={};Sens={}".format(ref_name, int(v.pos_in_rna) + 1, sens)
                        formt = "GT:XC"
                        ## geno decided on rate alt
                        ## if rate_alt > INCONSISTENT RATE THRESHOLD and 1 - rate_alt > INCONSISTENT RATE then hetero else if rate_alt > INCONSISTENT : homo_alt else : homo_ref
                        if rate_alt >= INCONSISTENT_RATE_THRESHOLD and (1-rate_alt) >= INCONSISTENT_RATE_THRESHOLD:
                            geno = "0/1"
                            qual = rate_alt * 2 if rate_alt < 0.5 else (1-rate_alt) *2
                        elif rate_alt >= INCONSISTENT_RATE_THRESHOLD:
                            geno = "1/1"
                            qual = 1 - (1-rate_alt)*10
                        else:
                            geno = "0/0"
                            qual = 1 - rate_alt*20
                        counts = "{},{}".format(counts_ref, counts_alt)
                        geno_s = "{}:{}".format(geno, counts)
                        qual = "."
                        filtr = "."
                        line = "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(v.chromosome, v.pos, v.rsID, v.ref, v.alt, qual, filtr, infos, formt, geno_s)
                        out.write(line)
