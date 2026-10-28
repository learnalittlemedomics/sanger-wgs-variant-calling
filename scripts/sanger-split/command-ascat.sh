do_parallel[ascat]="ascat.pl \
 -o $OUTPUT_DIR/${PROTOCOL}_${NAME_MT}_vs_${NAME_WT}/ascat \
 -t $BAM_MT_TMP \
 -n $BAM_WT_TMP \
 -sg $REF_BASE/ascat/SnpGcCorrections.tsv \
 -r $REF_BASE/genome.fa \
 -q 20 \
 -g L \
 -l $REF_BASE/gender.tsv \
 -rs '$SPECIES' \
 -ra $ASSEMBLY \
 -pr $PROTOCOL \
 -pl ILLUMINA \
 -c $CPU\
 -force\
 $ASCAT_ADD_ARGS"
