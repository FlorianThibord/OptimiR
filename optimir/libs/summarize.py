#!/usr/bin/env python
# -*- coding: utf-8 -*-
# (08/09/17)
########################################################
#                  OPTIMIR PIPELINE                    #
########################################################
# Unify outputs from pipeline for miRs, polymiRs, and isomiRs
# Usage : python make_abundance_table.py /path/to/OptimiR/Results/dir
# Standard libraries
import sys, os

def get_sample_list(path_dir):
    file_list = os.listdir(path_dir)
    categories = set()
    samples_name = set()
    for f in file_list:
        if not('.annot' in f) and not('inconsistents.' in f) and not('.gff3' in f) and not('.vcf' in f):
            cat, name = f.split('_abundances.txt')[0].split('.')
            categories.add(cat)
            samples_name.add(name)
    return list(samples_name), list(categories)

def process_annotations(table_annot, total):
    parent_list = [p for p in table_annot['parents'] if p != ""]
    cross_mapp_list = [cm for cm in table_annot['cross_mapps'] if cm != ""]
    cross_mapp_glob = 0
    ## parents
    parent_dict = {}
    for parents in parent_list:
        parents_l = parents.split('/')
        for parent in parents_l:
            counts = float(parent.split('|')[1].split('reads')[0])
            if '+' in parent:
                name = "multiple"
            else:
                name = parent.split('[')[0]
            if name in parent_dict:
                parent_dict[name] += counts
            else:
                parent_dict[name] = counts
    parent_out = ["{}[{}%]".format(parent_name, "%.2f" % (counts*100/float(total))) for parent_name, counts in parent_dict.items()]
    ## cross_mapps
    cross_mapp_dict = {}
    for cross_mapps in cross_mapp_list:
        cross_mapps_l = cross_mapps.split('/')
        for cross_mapp in cross_mapps_l:
            counts = float(cross_mapp.split('|')[1].split('reads')[0])
            cross_mapp_glob += counts
            name = cross_mapp.split('[')[0]
            if name in cross_mapp_dict:
                cross_mapp_dict[name] += counts
            else:
                cross_mapp_dict[name] = counts
    cross_mapp_out = ["{}[{}({}%)]".format(cross_mapp_name, counts, "%.2f" % (counts*100/float(total))) for cross_mapp_name, counts in cross_mapp_dict.items()]
    return '/'.join(parent_out), '/'.join(cross_mapp_out), cross_mapp_glob

def write_abundances_table(table, sample_list, cat, path_dir):
    sample_list.sort()
    ref_list = list(table.keys())
    ref_list.sort()
    output = open('{}/../{}_abundances.txt'.format(path_dir, cat), 'w')
    header = ['REFERENCE']
    header.extend(sample_list)
    header.extend(['TOTAL', "PARENTS", "CROSS_MAPPING", "CROSS_MAPPING_COUNTS"])
    output.write('\t'.join(header) + "\n")
    for ref in ref_list:
        l = [ref]
        total = 0
        for sample in sample_list:
            try:
                count = table[ref][sample]
            except KeyError:
                count = '0'
            l.append(count)
            total += float(count)
        parent, cross_mapp, cross_mapp_glob = process_annotations(table[ref], total)
        l.extend([str(total), parent, cross_mapp, str(cross_mapp_glob)])
        output.write('\t'.join(l) + "\n")
    output.close()

def print_table(list_files, path_out, path_in):
    header_printed = False
    with open(path_out, 'w') as out:
        for f in list_files:
            with open(path_in + "/" + f, 'r') as in_f:
                header = in_f.readline()
                if not(header_printed):
                    out.write(header)
                    header_printed = True
                for l in in_f:
                    out.write(l)

def print_iso_dist(list_files, path_out, path_in):
    header = "\t".join(["Reference" ,"Canonical_PCT", "End3_Canonical_PCT", "End5_Canonical_PCT", "End3_Trim_PCT", "End5_Trim_PCT", "End3_Tail_NT_PCT", "End5_Tail_NT_PCT", "End3_Tail_TE_PCT", "End5_Tail_TE_PCT", "End3_TrimTail_PCT", "End5_TrimTail_PCT", "Total"])
    iso_dict = {}
    for f in list_files:
        with open(path_in + "/" + f, 'r') as in_f:
            h = in_f.readline()
            for l in in_f:
                elts = l.split('\t')
                reference_name, cano, i3_cano, i5_cano, i3_trim, i5_trim, i3_tail, i5_tail, i3_TE, i5_TE, i3_trimtail, i5_trimtail, total = elts
                if reference_name not in iso_dict:
                    iso_dict[reference_name] = {"cano":0, "i3_cano":0, "i5_cano":0, "i3_trim":0, "i5_trim":0, "i3_tail":0, "i5_tail":0, "i3_TE":0, "i5_TE":0, "i3_trimtail":0, "i5_trimtail":0, "total":0}
                ## Add counts
                iso_dict[reference_name]["cano"] += float(cano) * float(total) / 100.
                iso_dict[reference_name]["i3_cano"] += float(i3_cano) * float(total) / 100.
                iso_dict[reference_name]["i5_cano"] += float(i5_cano) * float(total) / 100.
                iso_dict[reference_name]["i3_trim"] += float(i3_trim) * float(total) / 100.
                iso_dict[reference_name]["i5_trim"] += float(i5_trim) * float(total) / 100.
                iso_dict[reference_name]["i3_tail"] += float(i3_tail) * float(total) / 100.
                iso_dict[reference_name]["i5_tail"] += float(i5_tail) * float(total) / 100.
                iso_dict[reference_name]["i3_TE"] += float(i3_TE) * float(total) / 100.
                iso_dict[reference_name]["i5_TE"] += float(i5_TE) * float(total) / 100.
                iso_dict[reference_name]["i3_trimtail"] += float(i3_trimtail) * float(total) / 100.
                iso_dict[reference_name]["i5_trimtail"] += float(i5_trimtail) * float(total) / 100.
                iso_dict[reference_name]["total"] += float(total)
    ref_list = list(iso_dict.keys())
    ref_list.sort()
    with open(path_out, 'w') as out:
        out.write(header + "\n")
        for reference_name in ref_list:
            total = float(iso_dict[reference_name]["total"])
            cano = float(iso_dict[reference_name]["cano"])
            i3_cano = float(iso_dict[reference_name]["i3_cano"])
            i5_cano = float(iso_dict[reference_name]["i5_cano"])
            i3_trim = float(iso_dict[reference_name]["i3_trim"])
            i5_trim = float(iso_dict[reference_name]["i5_trim"])
            i3_tail = float(iso_dict[reference_name]["i3_tail"])
            i5_tail = float(iso_dict[reference_name]["i5_tail"])
            i3_TE = float(iso_dict[reference_name]["i3_TE"])
            i5_TE = float(iso_dict[reference_name]["i5_TE"])
            i3_trimtail = float(iso_dict[reference_name]["i3_trimtail"])
            i5_trimtail = float(iso_dict[reference_name]["i5_trimtail"])
            ## This should never happen
            if total == 0:
                total = 1
                print("Essentials WARNING : NUL TOTAL in isomiRs distribution calculation")
            l = "\t".join([reference_name, str(cano * 100. / total), str(i3_cano * 100. / total), str(i5_cano * 100. / total), str(i3_trim * 100. / total), str(i5_trim * 100. / total), str(i3_tail * 100. / total), str(i5_tail * 100. / total), str(i3_TE * 100. / total), str(i5_TE * 100. / total), str(i3_trimtail * 100. / total), str(i5_trimtail * 100. / total), str(total)]) + "\n"
            out.write(l)

## GFF merge (todo: clean & comment, manage filter for variants)
def print_gff(path_dir):
    if path_dir.endswith('/'):
        path_dir = path_dir[:-1]
    out_filename = '{}/../results.gff3'.format(path_dir)
    list_files = [f for f in os.listdir(path_dir) if ".gff3" in f]
    d = {}
    header = []
    samples = []
    d_expressions = {} ## iso -> sample -> expr
    d_expressions_opt = {} ## weighted expression by optimiR
    for file_in in list_files:
        with open(path_dir + "/" + file_in, 'r') as f:
            lines = f.readlines()
            header = [lines[0], lines[1]]
            sample_name = lines[2].split('\n')[0].split(': ')[1]
            samples.append(sample_name)
            for l in lines[4:]:
                elts = l.split(';\n')[0].split('\t')
                attributes = elts[-1].split(';')
                d_attr = {}
                d_attr["line"] = "\t".join(elts[:8])
                for a in attributes:
                    attr_name, attr_val = a.split("=")
                    d_attr[attr_name] = attr_val
                entry_id = d_attr["Name"]+"_"+elts[3]+"_"+elts[4]+"_"+elts[5]+"_"+d_attr["UID"]
                if not(entry_id in d):
                    d[entry_id] = d_attr
                for attr,value in d[entry_id].items():
                    if attr == "Expression":
                        if entry_id not in d_expressions:
                            d_expressions[entry_id] = {}
                        d_expressions[entry_id][sample_name] = d_attr[attr]
                    elif attr == "expression_OptimiR":
                        if entry_id not in d_expressions_opt:
                            d_expressions_opt[entry_id] = {}
                        d_expressions_opt[entry_id][sample_name] = d_attr[attr]
                    else:
                        if value != d_attr[attr]:
                            print("Difference with {}: {} vs {}".format(entry_id, value, d_attr[attr]))
    with open(out_filename, 'w') as out:
        samples.sort()
        header = "{}{}## COLDATA: {}\n## Processed with OptimiR.\n##\n".format(header[0], header[1], ','.join(samples))
        out.write(header)
        keylist = sorted(d)
        for key in keylist:
            entry = d[key]
            attr_list = ["UID", "Read","Parent","Name","Variant","Change","Cigar","Alias","Expression","expression_OptimiR","Filter","Hits"]
            final_attr = []
            for attr in attr_list:
                if attr == "Expression":
                    expr_list = []
                    for sample_name in samples:
                        try:
                            expr_list.append(str(int(float(d_expressions[key][sample_name]))))
                        except KeyError:
                            expr_list.append("0")
                    final_attr.append("Expression={}".format(",".join(expr_list)))
                elif attr == "expression_OptimiR":
                    expr_list = []
                    for sample_name in samples:
                        try:
                            expr_list.append(d_expressions_opt[key][sample_name])
                        except KeyError:
                            expr_list.append("0.0")
                    final_attr.append("expression_OptimiR={}".format(",".join(expr_list)))
                else:
                    try:
                        value = entry[attr]
                        final_attr.append("{}={}".format(attr,value))
                    except KeyError:
                        pass
            line = entry["line"] + "\t" + ";".join(final_attr) + ";\n"
            out.write(line)
    return d_expressions

def get_vcf_header(sample_list):
    return '##fileformat=VCFv4.1\n##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype inferred">\n##FORMAT=<ID=XC,Number=2,Type=Float,Description="Counts received on allele ref, Counts received on allele alt">\n##INFO=<ID=miR,Number=1,Type=String,Description="mature miRNA (or polymiR) from which genotype has been inferred (using expression of both alleles)">\n##INFO=<ID=Pos_miR,Number=1,Type=Integer,Description="SNP position in mature miR sequence">\n##INFO=<ID=Sens,Number=1,Type=String,Description="DNA Strand from which is transcribed the miRNA">\n##REF,ALT:Reference and alternative alleles as present in the mature miRNA. They might differ from DNA alleles, depending on miRNA sens (see INFOS)\n##Generated with OptimiR: based on the number of reads received by each polymiR allele. The decision depends on the INCONSISTENT_RATE_THRESHOLD defined by the user (default:0.01, meaning that, to be called, an allele must gather at least 1% of the reads aligned to the polymiR)\n#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	{}\n'.format("\t".join(sample_list))

def print_vcf(path_dir):
    if path_dir.endswith('/'):
        path_dir = path_dir[:-1]
    out_filename = '{}/../genotypes_OptimiR.vcf'.format(path_dir)
    list_files = [f for f in os.listdir(path_dir) if ".vcf" in f]
    d_vcf = {} ##[line_elts[0:9] as str] -> [sample_name] -> geno_s
    samples = []
    for file_in in list_files:
        with open(path_dir + "/" + file_in, 'r') as f:
            lines = f.readlines()
            for line in lines:
                if line[0] == "#" and line[1] != "#":
                    sample = line.split('\n')[0].split('\t')[-1]
                    samples.append(sample)
                elif line[0] != "#":
                    elts = line.split('\n')[0].split('\t')
                    key = "\t".join(elts[0:9])
                    if key not in d_vcf:
                        d_vcf[key] = {}
                    d_vcf[key][sample] = elts[-1]
    rs_list = list(d_vcf.keys())
    rs_list.sort()
    samples.sort()
    with open(out_filename, 'w') as out:
        out.write(get_vcf_header(samples))
        for rs in rs_list:
            out_genos = []
            for sample in samples:
                try:
                    out_genos.append(d_vcf[rs][sample])
                except KeyError:
                    out_genos.append("./.:0,0")
            line = "{}\t{}\n".format(rs, "\t".join(out_genos))
            out.write(line)
            

def print_tables(path_dir):
    sample_list, categories = get_sample_list(path_dir)
    for cat in categories:
        table = {}
        for sample in sample_list:
            sample_file = open('{}/{}.{}_abundances.txt'.format(path_dir, cat, sample), 'r')
            header = sample_file.readline()
            for l in sample_file:
                elts = l.split('\n')[0].split('\t')
                ref, counts = elts[:2]
                if len(elts) == 4:#hairpin and cross_mapp availables
                    parent, cross_mapp = elts[2:]
                else:
                    parent, cross_mapp = "", ""
                if ref in table:
                    table[ref][sample] = counts
                    table[ref]["parents"].append(parent)
                    table[ref]["cross_mapps"].append(cross_mapp)
                else:
                    table[ref] = {sample : counts, "parents" : [parent], "cross_mapps" : [cross_mapp]}
            sample_file.close()
        write_abundances_table(table, sample_list, cat, path_dir)
    ## Gather polymiRs tables
    list_poly_files = [f for f in os.listdir(path_dir) if "polymiRs_table." in f]
    path_table_poly = '{}/../polymiRs_table.annot'.format(path_dir)
    if list_poly_files != []:
        print_table(list_poly_files, path_table_poly, path_dir)
    list_consist_files = [f for f in os.listdir(path_dir) if "consistency_table." in f]
    path_table_consist = '{}/../consistency_table.annot'.format(path_dir)
    if list_consist_files != []:
        print_table(list_consist_files, path_table_consist, path_dir)
    list_iso_files = [f for f in os.listdir(path_dir) if "isomiRs_dist." in f]
    path_iso_dist = '{}/../isomiRs_dist.annot'.format(path_dir)
    if list_iso_files != []:
        print_iso_dist(list_iso_files, path_iso_dist, path_dir)
    if len([f for f in os.listdir(path_dir) if ".gff3" in f]) > 0:
        print_gff(path_dir)
    if len([f for f in os.listdir(path_dir) if ".vcf" in f]) > 0:
        print_vcf(path_dir)



def summarize(args):
    print('\n#############################')
    print('#  OPTIMIR: Summary tables  #')
    print('#############################')
    path = args.DIR
    print_tables(path)
    print("OptimiR summary tables are available.")

