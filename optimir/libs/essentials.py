#!/usr/bin/env python
#-*- coding: utf-8 -*-
# (30/11/16)
########################################################
#                 OptimiR essentials                   #
########################################################
# This library contains essential classes and functions for OptimiR

from pysam import AlignmentFile
import os
import pickle
import re
from copy import deepcopy
from itertools import combinations
import hashlib
import subprocess

###########################################################
## Classes definition : Coordinates, Variant, Sequence
###########################################################

## Compatible with miRBase and miRCarta gff3 (todo: mirgenedb)
class Coordinates:
    """ Coordinates fields are defined as:
    - chromosome : str
    - start position : str
    - end position : str
    - sens : ('+' or '-')
    - ID : str (miRBase ID)
    - Alias : str (optional)
    - Name : str (optional)
    - Derives from : str (if entry_type = "miRNA")

    Coordinates infos are retrieved from gff3 file.
    entry argument must be a gff3 line
    """

    def __init__(self, entry = ""):
        self.chromosome = ""
        self.entry_type = ""
        self.start = 0
        self.end = 0
        self.sens = ""
        self.ID = ""
        self.Alias = ""
        self.Name = ""
        self.DerivesFrom = ""
        if entry != "":
            self.update_coordinates(entry)

    def update_coordinates(self, entry):
        elts = entry.split('\n')[0].split('\t')
        self.chromosome = elts[0]
        self.entry_type = elts[2]
        self.start = int(elts[3])
        self.end = int(elts[4])
        self.sens = elts[6]
        infos = elts[8].split(';')
        for info in infos:
            if "ID" in info:
                self.ID = info.split('=')[1]
            if "Alias" in info:
                self.Alias = info.split('=')[1]
            if "Name" in info:
                self.Name = info.split('=')[1]
            if "Derives_from" in info:
                self.DerivesFrom = info.split('=')[1]

    def __str__(self):
        return "ch {} : {}-{} ({}) - {} - ID={}, Alias={}, Name={}, Derives_from={}".format(self.chromosome, self.start, self.end, self.sens, self.entry_type, self.ID, self.Alias, self.Name, self.DerivesFrom)

class Variant:
    """ Variant fields are defined as:
    rsID : str
    chromosome : str
    pos : int
    ref : str
    alt : str
    pos_in_rna : int
    GT_dict : dict [sample_name] -> genotype
    infos : str
    """

    def __init__(self, line = "", rsID = '', chromosome = "", pos = 0, ref = "", alt = "", qual = 0.0, filtr=".", infos = "", genotype_format = "", GT_dict = {}, genotype_other_list = [], sample_list = []):
        self.rsID = rsID
        self.chromosome = chromosome
        self.pos = pos
        self.ref = ref
        self.alt = alt
        self.qual = qual
        self.filtr = filtr
        self.infos = infos
        self.pos_in_rna = -1 ## position starts at 0 (1st nucleotide starting 5' end)
        self.genotype_format = genotype_format
        self.GT_dict = GT_dict
        self.genotype_other_list = genotype_other_list
        if line != "":
            self.update_variant(line, sample_list)

    def update_variant(self, line, sample_list):
        elts = line.split('\n')[0].split('\t')
        self.rsID = elts[2]
        self.chromosome = elts[0].replace("chr", "")
        self.pos = int(elts[1])
        self.ref = elts[3]
        self.alt = elts[4]
        self.qual = float(elts[5])
        self.filtr = elts[6]
        self.infos = elts[7]
        if len(elts) > 8: ## genotypes availables
            self.genotype_format = elts[8].split(':')
            for elt_format, i in zip(self.genotype_format, range(0, len(self.genotype_format))):
                if elt_format == "GT":
                    GT_list = [geno.split(':')[i] for geno in elts[9:]] ## Always the first elt (i=0)
                    for GT, sample in zip(GT_list, sample_list):
                        self.GT_dict[sample] = GT
            if len(elts[9].split(':')) > 1:
                self.genotype_other_list = [":".join(geno.split(':')[1:]) for geno in elts[9:]]

    def __str__(self):
        return "{}: {}-{}, {}/{}, rsq={} infos={}, geno_nb={}, geno_other_nb={}".format(self.rsID, self.chromosome, self.pos, self.ref, self.alt, self.qual, self.infos, len(self.GT_dict), len(self.genotype_other_list))

class Sequence:
    """ Sequence fields are defined as:
    - sequence : str
    - description : str
    - coordinates : Coordinates
    - is_polymiR : bool
    - variants : Variant list (optional)
    - alternative_sequences : dict [polymiR ID] -> alternative str sequence (optional)
    - hairpins_dict : dict [mature ident / polymiR ID] -> [parental hairpin id (only one entry)] -> hairpin sequence"""
    
    def __init__(self, seq = '', description = '', coordinates = Coordinates(), is_polymiR = False):
        self.sequence = seq
        self.description = description
        self.ident = description.split(' ')[0]
        self.coordinates = coordinates
        self.is_polymiR = is_polymiR
        self.variants = []
        self.alternative_sequences = {}
        self.hairpins_dict = {}

    def add_variant(self, variant):
        self.is_polymiR = True
        self.variants.append(variant)

    def add_alternative_sequences(self):
        def make_sequence_alt(seq, variant):
            ref = variant.ref
            alt = variant.alt
            pos = variant.pos_in_rna
            if seq[pos:pos + len(ref)] == ref:
                start = seq[:pos]
                end = seq[pos + len(ref):]
                out_seq = start + alt + end
            elif seq[:pos + len(ref)] in ref: #deletion outside the sequence
                end = seq[pos + len(ref):]
                out_seq = alt + end
            elif seq[pos:] in ref: #deletion outside the sequence
                start = seq[:pos]
                out_seq = start + alt
            else:
                print("add_alternative_sequence Error : rsid : {} ref : {} alt : {} pos : {} seq : {}".format(variant.rsID, ref, alt, pos, seq))
                out_seq = ""
            return out_seq
        for variant in self.variants:
            ## 1) change other alternative sequences (combinations of indels don't work)
            alternative_sequences_copy = deepcopy(self.alternative_sequences)
            for id_alt, alt_seq in self.alternative_sequences.items():
                alternative_sequences_copy[id_alt + "_" + variant.rsID] = make_sequence_alt(alt_seq, variant)
            ## 2) change reference sequence with alternative allele
            id_alt = self.ident + "_" + variant.rsID
            alternative_sequences_copy[id_alt] = make_sequence_alt(self.sequence, variant)
            self.alternative_sequences = alternative_sequences_copy

    def variant_in_sequence(self, variant):
        seq_chrom = self.coordinates.chromosome
        seq_start = self.coordinates.start
        seq_end = self.coordinates.end
        seq_sens = self.coordinates.sens
        var_chrom = variant.chromosome
        var_pos = variant.pos
        var_ref = variant.ref
        var_alt = variant.alt
        if seq_chrom.replace('chr', '') == var_chrom.replace('chr', ''):
            if var_pos >= seq_start and var_pos <= seq_end:
                variant_copy = deepcopy(variant)
                if seq_sens == "+":
                    pos_in_rna = var_pos - seq_start
                    if not(self.sequence[pos_in_rna:(pos_in_rna + len(var_ref))] == var_ref):
                        raise VariantSequenceError(variant.rsID, self.description, self.sequence[pos_in_rna:(pos_in_rna + len(var_ref))], var_ref, pos_in_rna)
                if seq_sens == "-":
                    pos_in_rna = seq_end - var_pos - len(var_ref) + 1
                    ref = rev_compl(var_ref)
                    variant_copy.ref = ref
                    variant_copy.alt = rev_compl(var_alt)
                    if not(self.sequence[pos_in_rna:(pos_in_rna + len(var_ref))] == ref):
                        raise VariantSequenceError(variant.rsID, self.description, self.sequence[pos_in_rna:(pos_in_rna + len(var_ref))], ref, pos_in_rna)
                variant_copy.pos_in_rna = pos_in_rna
                self.add_variant(variant_copy)

    def add_hairpin_seq(self, d_hairpins, d_ident, d_coordinates):
        hairpin_id = self.coordinates.DerivesFrom
        hairpin_name = d_coordinates[hairpin_id].Name
        if hairpin_name == "":
            hairpin_name = hairpin_id ## miRCarta does not have a name
        hairpin_seq = d_hairpins[hairpin_name]
        self.hairpins_dict[self.ident] = {hairpin_name: hairpin_seq}

    def add_alternative_sequences_hairpins(self):
        ## at this point, only one sequence in hairpins_dict
        if len(self.hairpins_dict[self.ident].items()) != 1:
            print("Hairpin issue!! {}".format(hairpins_dict))
        hairpin_id, hairpin_seq = list(self.hairpins_dict[self.ident].items())[0] #todo: clean
        if hairpin_seq.count(self.sequence) != 1: ## mature seq must occur only once in hairpin
            print("Mature occurence in hairpin!! {} : {}".format(hairpin_seq, self.sequence))
        for polymiR_id, alternative_sequence in self.alternative_sequences.items():
            alternative_hairpin_id = hairpin_id + "_" + "_".join(polymiR_id.split('_')[1:])
            alternative_hairpin_sequence = hairpin_seq.replace(self.sequence, alternative_sequence)
            self.hairpins_dict[polymiR_id] = {alternative_hairpin_id: alternative_hairpin_sequence}

    def __str__(self):
        return "Sequence object: \nSequence : {}\nDescription : {}\nCoordinates : {}\nis_polymiR : {}\nVariant : {}\nAlternative sequences : {}\nhairpins dict : {}\n".format(self.sequence, self.description, self.coordinates, self.is_polymiR, str(self.variants), self.alternative_sequences, self.hairpins_dict)

## Condensed sequence object, for OptimiR use
class OptimiR:
    """ OptimiR object fields are defined as:
    - sequence : str
    - name : str, the original name of the polymiR
    - ident : str
    - coordinates : Coordinates
    - polymiR_list : str list
    - variants : Variant list (optional)
    - hairpins : dict [hairpin name] -> hairpin sequence
    """
    
    def __init__(self, seq = '', name = '', ident = '', coordinates = Coordinates(), polymiR_list = [], variants = [], hairpins = {}):
        self.sequence = seq
        self.name = name
        self.ident = ident
        self.coordinates = coordinates
        self.polymiR_list = polymiR_list
        self.variants = variants
        self.hairpins = hairpins

    
###########################################################
## Exceptions definition 
###########################################################
                
class InputError(Exception):
    """Exception raised for errors in the input.

    Attributes:
        input_name -- path given for the input file
    """

    def __init__(self, input_name):
        self.input_name = input_name

class MissingIdCoordinates(Exception):
    """Exception raised if ID or NAME missing in INFO fields of gff3 file.

    Attributes:
        line -- line with wissing ID or Name from gff3
    """

    def __init__(self, line):
        self.line = line

class VariantSequenceError(Exception):
    """Exception raised is variant does not match the fasta sequence.

    Attributes:
        rsID -- variant rsID
        seq_description -- sequence description from fasta
        nt_in_seq -- nucleotides from sequence
        var_ref -- reference allele
        pos_in_mature -- position in sequence
    """

    def __init__(self, rsID, seq_description, nt_in_seq, var_ref, pos_in_rna):
        self.rsID = rsID
        self.seq_description = seq_description
        self.nt_in_seq = nt_in_seq
        self.var_ref = var_ref
        self.pos_in_rna = pos_in_rna

class Except_OS_Pipe(Exception): pass
"""Exception raised for errors linked to subprocess calls."""


###########################################################
## Base functions 
###########################################################

def rev_compl(seq):
    d_rev = {"A":"T", "T":"A", "C":"G", "G":"C"}
    return "".join([d_rev[s] for s in seq])[::-1]

## Replace Us by Ts
def parse_fasta(fasta_file):
    """Create dict with seq description as keys and sequences as values from fasta file"""
    seq_desc = ""
    seq = ""
    d_fasta = {}
    with open(fasta_file, 'r') as fa:
        for l in fa:
            if l.startswith('>'):
                if seq != "":
                    d_fasta[seq_desc] = seq.replace("U", "T")
                seq_desc = l.split('\n')[0].replace('>', '')
                seq = ""
            else:
                seq += l.split('\n')[0]
        d_fasta[seq_desc] = seq.replace("U", "T")
    return d_fasta

def write_fasta(d_fa, out_name):
    """Write fasta file from dict with seq description and sequences"""
    with open(out_name, 'w') as out:
        for seq_desc, seq in d_fa.items():
            if "hsa" in seq_desc:
                out.write(seq_desc + "\n")
                out.write(seq + "\n")
        
def check_if_file_exists(filename):
    if filename != "":
        if not(os.path.exists(filename)):
            raise InputError(filename)
    
def bamFile_to_dict(bam):
    bam_dict = {}
    for alignment in bam:
        query_name = alignment.query_name
        if query_name in bam_dict:
            bam_dict[query_name].append(alignment)
        else:
            bam_dict[query_name] = [alignment]
    return bam_dict

def bamDict_to_file(bam_dict, bam_out, collapse_table):
    alignments = [a for alignments in bam_dict.values() for a in alignments]
    for a in alignments:
        read_name = a.query_sequence
        sequence_obj = collapse_table[read_name] ## Uncollapse here!
        for (read_name, quality) in sequence_obj.backup:
            a_copy = deepcopy(a)
            a_copy.query_name = read_name
            a_copy.qual = quality
            bam_out.write(a_copy)
            
def get_hash(f_path, mode='md5'):
    h = hashlib.new(mode)
    with open(f_path, 'rb') as file:
        block = file.read(512)
        while block:
            h.update(block)
            block = file.read(512)
    return h.hexdigest()


## PRINTINGS, AESTHETICS
## TODO : Rewrite fun_str_progress, dict in essentials [key] -> string to print
def fun_str_progress(infos_list, index, verbose=True):
    infos_list.append("")
    header = ("\n" +
              "#########################\n" +
	      "##       OPTIMIR       ##\n" +
	      "#########################\n\n" +
	      "Starting workflow for sample {}\n > Starting Library preparation...")
    footer = "\nOptimiR Run Completed.\nResults are available in {}\nTotal time: {}s\n\n"
    
    str_dict = {"header": header,
                "footer": footer,
                "lib_prep": "   Library prepared! VCF available : {}; Geno available : {}; Elapsed time: {}s\n > Starting trimming reads... (can take several minutes)",
                "trim": "   Trimming done! Elapsed time: {}s\n > Starting collapsing reads... (can take several minutes)",
		"collaps" : "   Collapsing done! Elapsed time: {}s\n > Starting mapping reads....",
                "mapping" : "   Mapping done! Elapsed time: {}s\n > Starting post-processing....",
                "postproc" : "   Post-process done! Elapsed time: {}s",
                "outputs" : "    - Generating outputs",
                "ambiguous" : "    - Resolving ambigous alignments",
                "alscore" : "    - Computing Alignment score",
                "geno" : "    - Checking genotype consistency"
    }
    if verbose:
        print(str_dict[index].format(*infos_list))



###########################################################
## PICKLE functions
###########################################################

def save_obj(obj, path):
    with open(path, 'wb') as f:
        pickle.dump(obj, f, protocol=2)
        
def load_obj(path):
    with open(path, 'rb') as f:
        return pickle.load(f)

            
###########################################################
## IsomiR distribution
###########################################################

## Write percentages of isomiR distributions for each miRNA
def write_isomiR_dist(bam_dict, sample_name, dir_results):
    out_name = dir_results + "/" + "isomiRs_dist." + sample_name + ".annot"
    alignments = [a for alignments in bam_dict.values() for a in alignments]
    iso_dict = {} # {ref_name} -> {"cano", "i3_cano", "i5_cano", "i3_trim", "i5_trim", "i3_tail", "i5_tail", "i3_TE", "i5_TE", "i3_trimtail", "i5_trimtail", "total"}
    for alignment in alignments:
        if not(alignment.has_tag("XX")):
            reference_name = alignment.reference_name
            counts = float(alignment.get_tag("XC")) * float(alignment.get_tag("XW"))            
            if reference_name not in iso_dict:
                iso_dict[reference_name] = {"cano":0, "i3_cano":0, "i5_cano":0, "i3_trim":0, "i5_trim":0, "i3_tail":0, "i5_tail":0, "i3_TE":0, "i5_TE":0, "i3_trimtail":0, "i5_trimtail":0, "total":0}
            isotype = alignment.get_tag("IT")
            iso5, iso3 = isotype.split('[')[1].split(']')[0].split(',')
            is_cano = True
            if "-" in iso5 and "+" in iso5:
                iso_dict[reference_name]["i5_trimtail"] += counts
                is_cano = False
            elif "-" in iso5:
                iso_dict[reference_name]["i5_trim"] += counts
                is_cano = False
            elif "+" in iso5:
                is_cano = False
                if iso5.split('+')[1].isupper():
                    iso_dict[reference_name]["i5_TE"] += counts
                else:
                    iso_dict[reference_name]["i5_tail"] += counts
            else:
                iso_dict[reference_name]["i5_cano"] += counts
            if "-" in iso3 and "+" in iso3:
                iso_dict[reference_name]["i3_trimtail"] += counts
                is_cano = False
            elif "-" in iso3:
                iso_dict[reference_name]["i3_trim"] += counts
                is_cano = False
            elif "+" in iso3:
                is_cano = False
                if iso3.split('+')[1].isupper():
                    iso_dict[reference_name]["i3_TE"] += counts
                else:
                    iso_dict[reference_name]["i3_tail"] += counts
            else:
                iso_dict[reference_name]["i3_cano"] += counts
            if is_cano:
                iso_dict[reference_name]["cano"] += counts
            iso_dict[reference_name]["total"] += counts
    with open(out_name, 'w') as out:
        header = "\t".join(["Reference" ,"Canonical_PCT", "End3_Canonical_PCT", "End5_Canonical_PCT", "End3_Trim_PCT", "End5_Trim_PCT", "End3_Tail_NT_PCT", "End5_Tail_NT_PCT", "End3_Tail_TE_PCT", "End5_Tail_TE_PCT", "End3_TrimTail_PCT", "End5_TrimTail_PCT", "Total"])
        out.write(header + "\n")
        RUN_cano, RUN_i3_cano, RUN_i5_cano, RUN_i3_trim, RUN_i5_trim, RUN_i3_tail, RUN_i5_tail, RUN_i3_TE, RUN_i5_TE, RUN_i3_trimtail, RUN_i5_trimtail, RUN_total = 0,0,0,0,0,0,0,0,0,0,0,0
        for reference_name, iso_d in iso_dict.items():
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
            RUN_cano = RUN_cano + iso_dict[reference_name]["cano"]
            RUN_i3_cano = RUN_i3_cano + iso_dict[reference_name]["i3_cano"]
            RUN_i5_cano = RUN_i5_cano + iso_dict[reference_name]["i5_cano"]
            RUN_i3_trim = RUN_i3_trim + iso_dict[reference_name]["i3_trim"]
            RUN_i5_trim = RUN_i5_trim + iso_dict[reference_name]["i5_trim"]
            RUN_i3_tail = RUN_i3_tail + iso_dict[reference_name]["i3_tail"]
            RUN_i5_tail = RUN_i5_tail + iso_dict[reference_name]["i5_tail"]
            RUN_i3_TE = RUN_i3_TE + iso_dict[reference_name]["i3_TE"]
            RUN_i5_TE = RUN_i5_TE + iso_dict[reference_name]["i5_TE"]
            RUN_i3_trimtail = RUN_i3_trimtail + iso_dict[reference_name]["i3_trimtail"]
            RUN_i5_trimtail = RUN_i5_trimtail + iso_dict[reference_name]["i5_trimtail"]
            RUN_total = RUN_total + iso_dict[reference_name]["total"]
        ## This should never happen
        if RUN_total == 0:
            RUN_total = 1
            print("Essentials WARNING : NUL TOTAL in isomiRs distribution calculation")
        # write total
        l = "\t".join(["TOTAL", str(RUN_cano * 100. / RUN_total), str(RUN_i3_cano * 100. / RUN_total), str(RUN_i5_cano * 100. / RUN_total), str(RUN_i3_trim * 100. / RUN_total), str(RUN_i5_trim * 100. / RUN_total), str(RUN_i3_tail * 100. / RUN_total), str(RUN_i5_tail * 100. / RUN_total), str(RUN_i3_TE * 100. / RUN_total), str(RUN_i5_TE * 100. / RUN_total), str(RUN_i3_trimtail * 100. / RUN_total), str(RUN_i5_trimtail * 100. / RUN_total), str(RUN_total)]) + "\n"
        out.write(l)

def sam_to_bam(tmpdir, sample, samtools):
    subprocess.call("cat {}/{}.cl.fq >> {}/{}.failed.fq".format(tmpdir, sample, tmpdir, sample), shell=True)
    subprocess.call("rm -f {}/{}.cl.fq".format(tmpdir, sample), shell=True)
    subprocess.call("mv {}/{}.ok.sam {}/{}.sam".format(tmpdir, sample, tmpdir, sample), shell=True)
    subprocess.call("{} view -bh {}/{}.sam -o {}/{}.bam".format(samtools, tmpdir, sample, tmpdir, sample), shell=True)
    subprocess.call("rm -f {}/{}.al.sam".format(tmpdir, sample), shell=True)## aligned by bowtie2
    subprocess.call("rm -f {}/{}.sam".format(tmpdir, sample), shell=True)## cleaned by filter_reads
    # create bam
    subprocess.call("{} sort {}/{}.bam -o {}/{}.s.bam".format(samtools, tmpdir, sample, tmpdir, sample), shell=True)
    subprocess.call("rm -f {}/{}.bam".format(tmpdir, sample), shell=True)
    subprocess.call("mv {}/{}.s.bam {}/{}.bam".format(tmpdir, sample, tmpdir, sample), shell=True)
    subprocess.call("{} index {}/{}.bam".format(samtools, tmpdir, sample), shell=True)
    # Gunzip failed fq
    subprocess.call("gzip -f {}/{}.failed.fq".format(tmpdir, sample), shell=True)

