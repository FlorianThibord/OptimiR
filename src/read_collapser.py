#-*- coding: utf-8 -*-
# Florian THIBORD (02/07/17)
########################################################
#                 READ COLLAPSER v4                    #
########################################################
# Read collapser for a single sample
# Allow uncollapsing 
# Collapse reads with identical sequences
# Compute Average BQ score between identical reads
# Add number of collapsed reads in id as "{}_{}".format(read.id, nb_reads) with nb_reads as the number of collapsed reads

from operator import add
from collections import Counter
import sys, os, time

class sequence:
    #each sequence gets a name, a total counter, an average quality and a dict with (name, qual) saved for each read
    cpt_seqs = 0
    def __init__(self, q, read_name):
        self.qual = [ord(value) for value in q]
        self.times = 1
        self.name = "readx{}".format(sequence.cpt_seqs)
        self.backup = [(read_name, q)]
        sequence.cpt_seqs += 1
    def update(self, q, read_name, counts = 1):
        now = self.qual
        self.qual = map(add, now, [ord(value) for value in q])
        self.times += counts
        self.backup.append((read_name, q))
    def get_sample_counts(self):
        try:
            return self.times
        except KeyError:
            print 'Key error in read collapser: no read {}'.format(self.name)
            return None
    def get_quality(self):
        average = map(lambda x: int(round(x/self.times)), self.qual)
        return [str(unichr(char)) for char in average]
        
def collapse(sample_name, in_dir, seq_table):
    fastq_file = "{}/{}.trimmed.fq".format(in_dir, sample_name)
    with open(fastq_file, 'r') as handle:
        for line in handle:
            if line.startswith("@"):
                name = line.split(' ')[0]
                seq = handle.next().strip()
                handle.next().strip()
                qual = handle.next().strip()
                if seq in seq_table:
                    seq_table[seq].update(qual, name)
                else:
                    seq_table[seq] = sequence(qual, name)

def write_output(out_file, seqs):
    idx = 0
    cpt = 0
    with open(out_file, 'w') as handle:
        for seq in seqs:
            idx += 1
            qual = "".join(seqs[seq].get_quality())
            name = seqs[seq].name
            counts = seqs[seq].times
            cpt += counts
            ## FILTERING BASED ON COUNTS??
            handle.write(("@{}_{}\n{}\n+\n{}\n").format(name, counts, seq, qual))
    return ' > nb of reads before : {}\n > nb of unique reads after : {}'.format(cpt, idx)

def collapse_sample(SAMPLE_NAME, input_directory, output_directory, seq_table = Counter()):
    seq_table = Counter()
    seq_table.name = SAMPLE_NAME
    collapse(SAMPLE_NAME, input_directory, seq_table)
    out_fn = '{}/{}.clpsd.fq'.format(output_directory, SAMPLE_NAME)
    report = write_output(out_fn, seq_table)
    return seq_table, report
