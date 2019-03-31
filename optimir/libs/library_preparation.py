#!/usr/bin/env python
# -*- coding: utf-8 -*-
# (30/11/16)
#############################################
##         OptimiR
##  Library preparation script
#############################################

## Personal libraries
from optimir.libs.essentials import *

## Other libraries
import pickle
import subprocess
import os, sys, subprocess
import time

###########################################################
## Functions definitions 
###########################################################

## Build alignment index with bowtie2-build
def make_index(BOWTIE2_BUILD, index_path, fasta_file):
    command_line = ("{} -f {} {} -q".format(BOWTIE2_BUILD, fasta_file, index_path))
    ret_code = subprocess.call(command_line, shell=True)
    if ret_code != 0:
        raise Except_OS_Pipe('Bowtie2 index creation failed: command line = {}'.format(command_line))

## If sequence is mature, add coordinates in d_mat_coord, and append ID
def add_in_d_mat_coord(d_mat_coord, sequence, ID):
    coordinates = sequence.coordinates
    chrom = coordinates.chromosome.replace("chr", "")
    start = coordinates.start
    end = coordinates.end
    for i in range(start, end + 1):
        key = "{}:{}".format(chrom, i)
        if key not in d_mat_coord:
            d_mat_coord[key] = []
        d_mat_coord[key].append(ID)

## 1) d_coordinates dict : [Name or ID] -> Coordinates(line)
##    if not "Name" then find "ID" #case: if not miRBase then miRCarta
## 2) Make Sequence entry
## 3) Make Sequence dict : [ID] -> Sequence
## 4) Make matures coordinates dict : [chr:pos] -> miR_ID_list (to map variants)
def make_sequence_dict(matures_filename, gff3_filename):
    """Build a dict from fasta sequences file and gff3 coordinates file"""
    d_coordinates = {}
    d_ident = {}
    with open(gff3_filename, 'r') as gff3:
        for line in gff3:
            if line[0] != "#":
                coordinates = Coordinates(line)
                if coordinates.ID != "":
                    d_coordinates[coordinates.ID] = coordinates
                else:
                    raise MissingIdCoordinates(line)
                if coordinates.Name != "": ## for miRBase fasta: Name as ident description
                    ID_or_Name = coordinates.Name
                else: ## else for miRCarta : ID as ident description
                    ID_or_Name = coordinates.ID
                if ID_or_Name not in d_ident:
                    d_ident[ID_or_Name] = []
                d_ident[ID_or_Name].append(coordinates.ID)
    d_fasta_matures = parse_fasta(matures_filename)
    d_sequence = {}
    d_mat_coord = {}
    for description, sequence in d_fasta_matures.items():
        ID_or_Name = description.split(' ')[0]
        try:
            for ID in d_ident[ID_or_Name]: ## 1 seq can have several locus of origin
                coordinates = d_coordinates[ID]
                d_sequence[ID] = Sequence(sequence, description, coordinates)
                add_in_d_mat_coord(d_mat_coord, d_sequence[ID], ID)
        except KeyError:
            print("LIBRARY PREPARATION WARNING while parsing sequences: ID or Name {} does not have coordinates.".format(ID_or_Name))
    return d_sequence, d_ident, d_coordinates, d_mat_coord

def make_polymiRs(d_sequence, vcf_filename, d_mat_coord):
    """Add Variants to Sequence and generate polymiRs sequences integrating these variants"""
    sample_list = []
    GENO_AVAIL = False
    with open(vcf_filename, 'r') as vcf:
        for line in vcf:
            if line[0] == "#" and line[1] != "#": ## then header: extract sample_list
                elts = line.split('\n')[0].split('\t')
                if len(elts) > 8: ## genotypes available
                    GENO_AVAIL = True
                    sample_list = elts[9:]
            if line[0] != "#": ## variants lines
                elts = line.split('\n')[0].split('\t')
                chrom = elts[0].replace("chr", "")
                pos = int(elts[1])
                ref = elts[3]
                set_mature_ids = set()
                for i in range(pos, pos + len(ref) + 1):
                    key = "{}:{}".format(chrom, i)
                    if key in d_mat_coord:
                        ID_list = d_mat_coord[key]
                        for ID in ID_list:
                            set_mature_ids.add(ID)

                for ID in set_mature_ids:
                    sequence = d_sequence[ID]
                    try:
                        sequence.variant_in_sequence(Variant(line, sample_list=sample_list))
                    except VariantSequenceError as err:
                        print("LIBRARY PREPARATION WARNING raised while mapping {} to: {}\nnts_in_seq : {}, all_ref : {}, position in rna : {}".format(err.rsID, err.seq_description, err.nt_in_seq, err.var_ref, err.pos_in_rna))
    for sequence in d_sequence.values():
        if sequence.is_polymiR:
            sequence.add_alternative_sequences()
    return GENO_AVAIL, sample_list

def make_hairpin_seqs(d_sequence, hairpins_filename, d_ident, d_coordinates):
    d_fasta_hairpins = parse_fasta(hairpins_filename)
    d_hairpins = {}
    for description, sequence in d_fasta_hairpins.items():
        ID_or_Name = description.split(' ')[0]
        d_hairpins[ID_or_Name] = sequence
    for sequence in d_sequence.values():
        sequence.add_hairpin_seq(d_hairpins, d_ident, d_coordinates)
        if sequence.is_polymiR:
            sequence.add_alternative_sequences_hairpins()
    
def write_fasta_from_Sequences(d_sequence, out_fasta):
    """Write new fasta file integrating alternative sequences"""
    set_of_desc = set() ## do not print the same mature twice for homologous
    with open(out_fasta, 'w') as out:
        for sequence in d_sequence.values():
            if not(sequence.description in set_of_desc):
                set_of_desc.add(sequence.description)
                out.write('>' + sequence.description + '\n')
                out.write(sequence.sequence + '\n')
            if sequence.is_polymiR:
                for desc, seq in sequence.alternative_sequences.items():
                    out.write('>' + desc + '\n')
                    out.write(seq + '\n')

def write_log_from_Sequences(d_sequence, out_log):
    """Write logging informations in out_log"""
    d_pol = {} ## "rs" -> list of miRs
    pols = set()
    sequences = 0
    variants = set()
    nb_multi_snp = 0
    nb_rs_on_multi = 0
    nb_miRs = 0
    for sequence in d_sequence.values():
        if sequence.is_polymiR:
            pols.add(sequence.coordinates.ID)
            for variant in sequence.variants:
                rsID = variant.rsID
                variants.add(rsID)
                if rsID not in d_pol:
                    d_pol[rsID] = []
                d_pol[rsID].append(sequence.ident)
            if len(sequence.variants) > 1:
                nb_multi_snp += 1
            sequences += len(sequence.alternative_sequences)
    for rs, list_miRs in d_pol.items():
        if len(list_miRs) > 1:
            nb_rs_on_multi += 1
    with open(out_log, 'w') as out:
        s = "Number of new sequences : {}\nNumber of variants : {}\nNumber of polymiRs : {}\nNumber of polymiRs with multiple SNPs : {}\nNumber of SNPs on multiple polymiRs : {}".format(sequences, len(variants), len(pols), nb_multi_snp, nb_rs_on_multi)
        out.write(s)

## d_OptimiR : [seq_name] -> OptimiR_obj : mature id, is_polymiR, variant(s) (if is_polymiR), sequence, hairpin sequence(s)
def make_OptimiR_dict(d_sequence, d_ident):
    d_OptimiR = {}
    ## d_ident contains idents for both matures and hairpin. First check that ident has a key in d_sequence, to be sure it is a mature ident
    for name, ident_list in d_ident.items():
        new_hairpins_dict = {}
        is_mature = False
        polymiR_list = []
        OptimiR_variant_list = []
        ## retrieve all polymiRs, and all hairpins associated. For each polymiR, create an OptimiR obj
        ## For each sequence, store an OptimiR obj with all hairpins associated
        for ident in ident_list:
            if ident in d_sequence:
                is_mature = True
                sequence_obj = d_sequence[ident]
                coordinates = sequence_obj.coordinates
                variant_list = sequence_obj.variants
                is_polymiR = (len(variant_list) > 0)
                hairpins = sequence_obj.hairpins_dict
                if is_polymiR:
                    polymiRs_sequences = sequence_obj.alternative_sequences
                    for polymiR_id, polymiR_sequence in polymiRs_sequences.items():
                        # Check that only one hairpin per polymiR
                        d_polymiR_hairpins = {}
                        if len(hairpins[polymiR_id].items()) > 1:
                            print("More than one hairpin for polymiR {} : {}".format(polymiR_id, hairpins[polymiR_id]))
                        for hairpin_polymiR_id, hairpin_polymiR_sequence in hairpins[polymiR_id].items():
                            d_polymiR_hairpins[hairpin_polymiR_id] = hairpin_polymiR_sequence
                        # retrieve variants
                        variant_name_list = polymiR_id.split('_')[1:]
                        variants = []
                        for variant in variant_list:
                            if variant.rsID in variant_name_list:
                                variants.append(variant)
                                OptimiR_variant_list.append(variant)
                        if polymiR_id in d_OptimiR:
                            print("polymiR already exists : {}".format(polymiR_id))
                        polymiR_list.append(polymiR_id)
                        d_OptimiR[polymiR_id] = OptimiR(polymiR_sequence, name, ident, coordinates, [], list(set(variants)), d_polymiR_hairpins)
                ## Retrieve hairpins
                for hairpin_name, hairpin_seq in hairpins[name].items():
                    if hairpin_name in new_hairpins_dict:
                        print("Hairpin allready present : {} in {}".format(hairpin_name, new_hairpin_dict))
                    new_hairpins_dict[hairpin_name] = hairpin_seq
        if is_mature:
            d_OptimiR[name] = OptimiR(sequence_obj.sequence, name, ident_list, sequence_obj.coordinates, polymiR_list, list(set(OptimiR_variant_list)), new_hairpins_dict)
    return d_OptimiR
        
def check_library_modifications(matures_filename, vcf_filename, lib_infos_pickle_path, d_OptimiR_pickle_path):
    """Return True if provided informations changed between OptimiR runs. Also returns lib_infos"""
    build_index = False
    ## Retrieve hash signature from provided files
    fasta_matures_hash = get_hash(matures_filename)
    if vcf_filename != None:
        vcf_hash = get_hash(vcf_filename)
    else:
        vcf_hash = "0"
    ## pickle part: check if the library with the given files and values have already been made. If not: build it, else skip the library preparation and reuse previous files
    try:
        check_if_file_exists(lib_infos_pickle_path)
        check_if_file_exists(d_OptimiR_pickle_path)
    except InputError as err: ## Meaning this is the first run using the provided directory
        build_index = True
    if not(build_index): ## Check provided files are the same than previous run
        previous_lib_infos = load_obj(lib_infos_pickle_path)
        (VCF_AVAIL, GENO_AVAIL, sample_list, fasta_matures_previous_hash, vcf_previous_hash) = previous_lib_infos
        identical_fasta = fasta_matures_hash == fasta_matures_previous_hash
        identical_vcf = vcf_hash == vcf_previous_hash
        if not(identical_fasta) or not(identical_vcf):
            print("Rebuild library")
            build_index = True
    return build_index, vcf_hash, fasta_matures_hash
        
## Main API
## todo: index sans vcf
def prepare_library(BOWTIE2_BUILD, VCF, MATURES, HAIRPINS, GFF3, out_directory, fasta_file, index_path, lib_infos_pickle_path, d_OptimiR_pickle_path):
    """Upgrade base library with new sequences integrating alternative alleles present in vcf"""
    VCF_AVAIL = True
    GENO_AVAIL = False
    try:
        if VCF == None:
            VCF_AVAIL = False
        else:
            check_if_file_exists(VCF)
        check_if_file_exists(MATURES)
        check_if_file_exists(HAIRPINS)
        check_if_file_exists(GFF3)
    except InputError as err:
        print("ERROR during library preparation: file {} does not exists. Check input filename and try again.\n".format(err.input_name))
        sys.exit(4)
    build_index, vcf_hash, fasta_matures_hash = check_library_modifications(MATURES, VCF, lib_infos_pickle_path, d_OptimiR_pickle_path)
    if build_index:
        ## Build dicts
        d_sequence, d_ident, d_coordinates, d_mat_coord = make_sequence_dict(MATURES, GFF3)
        if VCF_AVAIL:
            GENO_AVAIL, sample_list = make_polymiRs(d_sequence, VCF, d_mat_coord)
        else:
            sample_list = []
        make_hairpin_seqs(d_sequence, HAIRPINS, d_ident, d_coordinates)
        ## Write new library in fasta file
        write_fasta_from_Sequences(d_sequence, fasta_file)
        ## Write log infos
        write_log_from_Sequences(d_sequence, out_directory + "lib_prep.log")
        ## Create d_OptimiR, a dictionary used by OptimiR
        d_OptimiR = make_OptimiR_dict(d_sequence, d_ident)
        ## save pickle objects
        lib_infos = (VCF_AVAIL, GENO_AVAIL, sample_list, fasta_matures_hash, vcf_hash)
        save_obj(d_OptimiR, d_OptimiR_pickle_path)
        save_obj(lib_infos, lib_infos_pickle_path) # matures(mirs list) and vcf (variants list, samples and genotypes)
        ## build bowtie2 index
        make_index(BOWTIE2_BUILD, index_path, fasta_file)

## Main Standalone
def library_preparation(args):
    print('\n##################################')
    print('#  OPTIMIR: Library Preparation  #')
    print('##################################')

    header = "Library preparation : \nStart process... Depending on the number of variants provided, the preparation can take a few seconds to several minutes..."
    footer = "Library built!"

    print(header)
    
    ## Assign arguments to global variables
    VCF = args.VCF
    MATURES = args.MATURES
    HAIRPINS = args.HAIRPINS
    GFF3 = args.GFF3
    OUTDIR = args.OUTDIR
    BOWTIE2_BUILD = args.BOWTIE2_BUILD

    if OUTDIR[-1] == "/":
        OUTDIR = OUTDIR[:-1]
    out_directory = os.path.abspath(OUTDIR) + "/OptimiR_lib/"
    subprocess.call("mkdir -p " + out_directory, shell=True)
    # Create directory for OptimiR fasta library
    fasta_dir =  out_directory + "fasta/"
    subprocess.call("mkdir -p " + fasta_dir, shell=True)
    fasta_file = fasta_dir + "optimiR_library.fa"
    # Create directory and path for bowtie2 index alignment
    index_dir =  out_directory + "bowtie2_index/"
    subprocess.call("mkdir -p " + index_dir, shell=True)
    index_path = index_dir + "optimiR_alignment_library"
    # Create directory for pickle objects
    pickle_dir = out_directory + "pickle/"
    subprocess.call("mkdir -p " + pickle_dir, shell=True)
    lib_infos_pickle_path = pickle_dir + "lib_infos.pkl"
    d_OptimiR_pickle_path = pickle_dir + "d_OptimiR.pkl"

    prepare_library(BOWTIE2_BUILD, VCF, MATURES, HAIRPINS, GFF3, out_directory, fasta_file, index_path, lib_infos_pickle_path, d_OptimiR_pickle_path)
    print(footer)
    
