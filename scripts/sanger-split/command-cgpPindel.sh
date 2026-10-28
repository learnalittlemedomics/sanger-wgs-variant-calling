# Pindel-specific parallelism; do not change shared run.params or run-stage.sh.
# The old global min(..., 8) is intentionally overridden only for this caller.
PINDEL_WORKERS=${PINDEL_WORKERS:-16}
case "$PINDEL_WORKERS" in
  ''|*[!0-9]*) echo 'PINDEL_WORKERS must be a positive integer' >&2; exit 2 ;;
esac
if [ "$PINDEL_WORKERS" -lt 1 ] || [ "$PINDEL_WORKERS" -gt 30 ]; then
  echo 'PINDEL_WORKERS must be in 1..30' >&2; exit 2
fi
PINDEL_CPU=$PINDEL_WORKERS
[ "$PINDEL_CPU" -le "$CPU" ] || PINDEL_CPU=$CPU
echo "Pindel concurrency: $PINDEL_CPU (allocated CPUs: $CPU)"
do_parallel[cgpPindel]="pindel.pl \
 -o $OUTPUT_DIR/${PROTOCOL}_${NAME_MT}_vs_${NAME_WT}/pindel \
 -r $REF_BASE/genome.fa \
 -t $BAM_MT_TMP \
 -n $BAM_WT_TMP \
 -s $REF_BASE/pindel/simpleRepeats.bed.gz \
 -u $REF_BASE/pindel/pindel_np.gff3.gz \
 -f $REF_BASE/pindel/${PROTOCOL}_Rules.lst \
 -g $REF_BASE/vagrent/codingexon_regions.indel.bed.gz \
 -st $PROTOCOL \
 -as $ASSEMBLY \
 -sp '$SPECIES' \
 -e $CONTIG_EXCLUDE \
 -b $REF_BASE/pindel/HiDepth.bed.gz \
 -c $PINDEL_CPU \
 -sf $REF_BASE/pindel/softRules.lst"
