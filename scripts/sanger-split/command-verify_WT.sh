do_parallel[verify_WT]="verifyBamHomChk.pl -d 25 \
    -o $OUTPUT_DIR/${PROTOCOL}_${NAME_WT}/contamination \
    -b $BAM_WT_TMP \
    -t $CPU \
    -j $OUTPUT_DIR/${PROTOCOL}_${NAME_WT}/contamination/result.json \
    -s $REF_BASE/verifyBamID_snps.vcf.gz"
