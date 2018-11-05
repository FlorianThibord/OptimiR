VCF=./example/genotypes.vcf
MATURES=./resources/fasta/hsa_matures_miRBase_v21.fa
HAIRPINS=./resources/fasta//hsa_hairpins_miRBase_v21.fa
GFF3=./resources/coordinates/hsa_miRBase_v21.gff3
DIR_OUT=./Results

# First prepare library
# python ./src/library_preparation.py \
#        -v ../$VCF \
#        -m ../$MATURES \
#        -p ../$HAIRPINS \
#        -g ../$GFF3 \
#        -o ../$DIR_OUT/OptimiR_lib/

# Launch OptimiR on each sample (parallel computation is recommanded for many samples)
./OPTIMIR \
    --fq ./example/S1.fq.gz \
    --dirOutput $DIR_OUT \
    --vcf $VCF
./OPTIMIR \
    --fq ./example/S2.fq.gz \
    --dirOutput $DIR_OUT \
    --vcf $VCF

# Make summary files
./OPTIMIR_SUMMARY $DIR_OUT/OptimiR_Results
