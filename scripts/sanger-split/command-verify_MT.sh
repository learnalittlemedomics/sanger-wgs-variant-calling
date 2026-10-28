do_parallel[verify_MT]="verifyBamHomChk.pl -d 25 \
   -o $OUTPUT_DIR/${PROTOCOL}_$NAME_MT/contamination \
   -b $BAM_MT_TMP \
   -t $CPU \
   -a $OUTPUT_DIR/${PROTOCOL}_${NAME_MT}_vs_${NAME_WT}/ascat/${NAME_MT}.copynumber.caveman.csv \
   -j $OUTPUT_DIR/${PROTOCOL}_${NAME_MT}/contamination/result.json \
   -s $REF_BASE/verifyBamID_snps.vcf.gz"
