#!/bin/bash
set -euo pipefail
SPLIT_STAGE=$1
SPLIT_CPUS=$2
shift 2
for archive in "$@"; do
  [ -f "$archive" ] && tar -xzf "$archive"
done
printf '%s\n' "options(bitmapType='cairo')" > .Rprofile
export R_PROFILE_USER="$PWD/.Rprofile"
set -- run.params
set -e

echo -e "\nStart workflow: `date`\n"

if [[ $# -eq 1 ]] ; then
  PARAM_FILE=$1
elif [ -z ${PARAM_FILE+x} ] ; then
  PARAM_FILE=$HOME/run.params
fi

echo "Loading user options from: $PARAM_FILE"
if [ ! -f $PARAM_FILE ]; then
  echo -e "\tERROR: file indicated by PARAM_FILE not found: $PARAM_FILE" 1>&2
  exit 1
fi
source $PARAM_FILE
CPU=$SPLIT_CPUS

if [ -z ${CPU+x} ]; then
  CPU=`grep -c ^processor /proc/cpuinfo`
fi

# calculate the min of user defined, host and max we should ever allow for pindel
min_cpu=`echo -e "$PINDEL_MAXCPU\n$CPU\n8" | sort -k1,1n | head -n 1`
PINDEL_CPU=$min_cpu

# create area which allows monitoring site to be started, not actively updated until after PRE-EXEC completes
#cp -r /opt/wtsi-cgp/site $OUTPUT_DIR/site

echo -e "\tBAM_MT : $BAM_MT"
echo -e "\tBAM_WT : $BAM_WT"

set -u
TMP=$OUTPUT_DIR/tmp
mkdir -p $TMP
mkdir -p $OUTPUT_DIR/timings

## get sample names from BAM headers
NAME_MT=`samtools view -H $BAM_MT | perl -ne 'chomp; if($_ =~ m/^\@RG/) {($sm) = $_ =~m/\tSM:([^\t]+)/; print "$sm\n";}' | uniq`
NAME_WT=`samtools view -H $BAM_WT | perl -ne 'chomp; if($_ =~ m/^\@RG/) {($sm) = $_ =~m/\tSM:([^\t]+)/; print "$sm\n";}' | uniq`

echo -e "\tNAME_MT : $NAME_MT"
echo -e "\tNAME_WT : $NAME_WT"

# capture index extension type (assuming same from both)
ALN_EXTN='bam'
IDX_EXTN=''
if [[ "$IDX_MT" == *.bam.bai ]]; then
  IDX_EXTN='bam.bai'
elif [[ "$IDX_MT" == *.bam.csi ]]; then
  IDX_EXTN='bam.csi'
elif [[ "$IDX_MT" == *.cram.crai ]]; then
  IDX_EXTN='cram.crai'
  ALN_EXTN='cram'
else
  echo "Alignment is not BAM or CRAM file: $" >&2
  exit 1
fi

BAM_MT_TMP=$TMP/$NAME_MT.$ALN_EXTN
IDX_MT_TMP=$TMP/$NAME_MT.$IDX_EXTN
BAM_WT_TMP=$TMP/$NAME_WT.$ALN_EXTN
IDX_WT_TMP=$TMP/$NAME_WT.$IDX_EXTN

## BAS files are not generated in this flow due to siginifcant run time overhead.

ln -fsr $BAM_MT $BAM_MT_TMP
ln -fsr $BAM_MT.bas $BAM_MT_TMP.bas
ln -fsr $IDX_MT $IDX_MT_TMP
ln -fsr $BAM_WT $BAM_WT_TMP
ln -fsr $BAM_WT.bas $BAM_WT_TMP.bas
ln -fsr $IDX_WT $IDX_WT_TMP

ASCAT_ADD_ARGS=''
# ASCAT_PURITY set
if [ ! -z ${ASCAT_PURITY+x} ]; then
  ASCAT_ADD_ARGS="$ASCAT_ADD_ARGS -pu $ASCAT_PURITY"
fi
# ASCAT_PLOIDY set
if [ ! -z ${ASCAT_PLOIDY+x} ]; then
  ASCAT_ADD_ARGS="$ASCAT_ADD_ARGS -pi $ASCAT_PLOIDY"
fi

## Make fake copynumber so we can run early steps of caveman
perl -alne 'print join(qq{\t},$F[0],0,$F[1],2);' < $REF_BASE/genome.fa.fai | tee $TMP/norm.cn.bed > $TMP/tum.cn.bed


USER_CACHE=$OUTPUT_DIR/ref_cache
if [ "$ALN_EXTN" = cram ]; then
  export REF_CACHE=$USER_CACHE/%2s/%2s/%s
  export REF_PATH=$REF_CACHE
fi
case "$SPLIT_STAGE" in
  CaVEMan|CaVEMan_scatter_prepare|CaVEMan_flag|verify_MT)
# prep ascat output for caveman and VerifyBam (Tumour):
set -x
ASCAT_CN="$OUTPUT_DIR/${PROTOCOL}_${NAME_MT}_vs_${NAME_WT}/ascat/$NAME_MT.copynumber.caveman.csv"
perl -ne '@F=(split q{,}, $_)[1,2,3,4]; $F[1]-1; print join("\t",@F)."\n";' < $ASCAT_CN > $TMP/norm.cn.bed
perl -ne '@F=(split q{,}, $_)[1,2,3,6]; $F[1]-1; print join("\t",@F)."\n";' < $ASCAT_CN > $TMP/tum.cn.bed
NORM_CONTAM=`perl -ne 'if(m/^rho\s(.+)\n/){print 1-$1;}' $OUTPUT_DIR/${PROTOCOL}_${NAME_MT}_vs_${NAME_WT}/ascat/$NAME_MT.samplestatistics.txt`
set +x


  ;;
esac
GERMLINE_BED=$OUTPUT_DIR/${PROTOCOL}_${NAME_MT}_vs_${NAME_WT}/pindel/${NAME_MT}_vs_${NAME_WT}.germline.bed
declare -A do_parallel
# Source all assignments only after their expansion variables have been defined.
NORM_CONTAM=${NORM_CONTAM:-}
# Stage-specific command input prevents unrelated caller edits invalidating this task.
COMMAND_STAGE=$SPLIT_STAGE
[ "$COMMAND_STAGE" != CaVEMan_scatter_prepare ] || COMMAND_STAGE=CaVEMan
if [ "$COMMAND_STAGE" != cache_POP ]; then
  source "split_scripts/command-$COMMAND_STAGE.sh"
fi
if [ "$SPLIT_STAGE" = cache_POP ]; then
  mkdir -p ref_cache
  if [ "$ALN_EXTN" = cram ]; then
    /usr/bin/time -v seq_cache_populate.pl -root "$USER_CACHE" "$REF_BASE/genome.fa" > "timings/${PROTOCOL}_${NAME_MT}_vs_${NAME_WT}.time.cache_POP" 2>&1
  else
    echo 'BAM input: no CRAM cache required' > "timings/${PROTOCOL}_${NAME_MT}_vs_${NAME_WT}.time.cache_POP"
  fi
  tar -czf cache_POP.stage.tar.gz timings
  exit 0
fi
case "$SPLIT_STAGE" in
  cgpPindel_annot) rm -f $OUTPUT_DIR/${PROTOCOL}_${NAME_MT}_vs_${NAME_WT}/pindel/${NAME_MT}_vs_${NAME_WT}.annot.vcf.gz* ;;
  CaVEMan_annot) rm -f $OUTPUT_DIR/${PROTOCOL}_${NAME_MT}_vs_${NAME_WT}/caveman/${NAME_MT}_vs_${NAME_WT}.annot.muts.vcf.gz* ;;
esac
if [ "$SPLIT_STAGE" = CaVEMan_scatter_prepare ]; then
  python3 split_scripts/caveman-scatter.py prepare \
    --command "${do_parallel[CaVEMan]}" \
    --cave-dir "$OUTPUT_DIR/${PROTOCOL}_${NAME_MT}_vs_${NAME_WT}/caveman" \
    --tumour "$BAM_MT" --normal "$BAM_WT" --tidx "$IDX_MT" --nidx "$IDX_WT" \
    --workers-per-job "${CAVE_WORKERS_PER_JOB:-6}" \
    --target-hours "${CAVE_TARGET_HOURS:-1.5}" \
    --mstep-minutes "${CAVE_MSTEP_MINUTES:-3.328}" \
    --estep-minutes "${CAVE_ESTEP_MINUTES:-15.83}"
  exit 0
fi
CMD="/usr/bin/time -v ${do_parallel[$SPLIT_STAGE]}"
echo "$CMD"
# Unlike the upstream parallel wrapper, this task directly propagates software exit status.
bash -c "$CMD" > "timings/${PROTOCOL}_${NAME_MT}_vs_${NAME_WT}.time.$SPLIT_STAGE" 2>&1
if [ "$SPLIT_STAGE" = cgpPindel ]; then
GERMLINE_BED=$OUTPUT_DIR/${PROTOCOL}_${NAME_MT}_vs_${NAME_WT}/pindel/${NAME_MT}_vs_${NAME_WT}.germline.bed
if [ -f $GERMLINE_BED ]; then
  # need to sort and index pindel germline
  sort -k1,1 -k2,2n -k3,3n $GERMLINE_BED | bgzip -c > $GERMLINE_BED.gz
  tabix -p bed $GERMLINE_BED.gz
  rm -f $GERMLINE_BED
fi


fi
case "$SPLIT_STAGE" in
  ascat) result_dir=${PROTOCOL}_${NAME_MT}_vs_${NAME_WT}/ascat ;;
  CaVEMan*) result_dir=${PROTOCOL}_${NAME_MT}_vs_${NAME_WT}/caveman ;;
  cgpPindel*) result_dir=${PROTOCOL}_${NAME_MT}_vs_${NAME_WT}/pindel ;;
  BRASS*) result_dir=${PROTOCOL}_${NAME_MT}_vs_${NAME_WT}/brass ;;
  geno) result_dir=${PROTOCOL}_${NAME_MT}_vs_${NAME_WT}/genotyped ;;
  verify_WT) result_dir=${PROTOCOL}_${NAME_WT}/contamination ;;
  verify_MT) result_dir=${PROTOCOL}_${NAME_MT}/contamination ;;
  *) echo "Unknown stage: $SPLIT_STAGE" >&2; exit 2 ;;
esac
# Immutable per-task snapshot; every downstream task extracts its own writable copy.
tar -czf "$SPLIT_STAGE.stage.tar.gz" "$result_dir" timings
