do_parallel[geno]="compareBamGenotypes.pl \
  -o $OUTPUT_DIR/${PROTOCOL}_${NAME_MT}_vs_${NAME_WT}/genotyped \
  -nb $BAM_WT_TMP \
  -j $OUTPUT_DIR/${PROTOCOL}_${NAME_MT}_vs_${NAME_WT}/genotyped/result.json \
  -tb $BAM_MT_TMP \
  -s $REF_BASE/general.tsv \
  -g $REF_BASE/gender.tsv"
