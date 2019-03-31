VCF=./example/genotypes.vcf
DIR_OUT=./OptimiR_Results_Dir
# MATURES=./optimir/resources/fasta/hsa_matures_miRBase_v21.fa
# HAIRPINS=./optimir/resources/fasta/hsa_hairpins_miRBase_v21.fa
# GFF3=./optimir/resources/coordinates/hsa_miRBase_v21.gff3

# First prepare library if OptimiR is launched in parallel jobs locally - or better, on a cluster. (In this example, they are launched sequentially)
optimir libprep \
       -v $VCF \
       -o $DIR_OUT # \
       # -m $MATURES \
       # -p $HAIRPINS \
       # -g $GFF3 

# Launch OptimiR on each sample (parallel computation is recommanded for many samples)
optimir process \
    --fq ./example/S1.fq.gz \
    --vcf $VCF \
    --gff_out \
    --dirOutput $DIR_OUT 

optimir process \
    --fq ./example/S2.fq.gz \
    --vcf $VCF \
    --gff_out \
    --dirOutput $DIR_OUT 

# Make summary files
optimir summarize --dir $DIR_OUT/OptimiR_Results
