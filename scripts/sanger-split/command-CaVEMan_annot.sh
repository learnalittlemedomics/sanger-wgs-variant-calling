do_parallel[CaVEMan_annot]="AnnotateVcf.pl -t -c $REF_BASE/vagrent/vagrent.cache.gz \
   -i $OUTPUT_DIR/${PROTOCOL}_${NAME_MT}_vs_${NAME_WT}/caveman/${NAME_MT}_vs_${NAME_WT}.flagged.muts.vcf.gz \
   -o $OUTPUT_DIR/${PROTOCOL}_${NAME_MT}_vs_${NAME_WT}/caveman/${NAME_MT}_vs_${NAME_WT}.annot.muts.vcf"
