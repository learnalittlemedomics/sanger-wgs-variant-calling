#!/bin/bash
set -euo pipefail
PINDEL_MODE=$1
SPLIT_CPUS=$2
shift 2
PINDEL_ARGS=("$@")
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
# Fail rather than silently fall back to header/default insert sizes.
case "${PINDEL_ARGS[*]}" in
 '--index 1') test -s "$BAM_MT.bas" ;;
 '--index 2') test -s "$BAM_WT.bas" ;;
 *) echo 'Early input requires exactly --index 1 or --index 2' >&2; exit 2 ;;
esac
CPU=$SPLIT_CPUS

if [ -z ${CPU+x} ]; then
  CPU=`grep -c ^processor /proc/cpuinfo`
fi

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
if [ -f "$BAM_MT.bas" ]; then ln -fsr "$BAM_MT.bas" "$BAM_MT_TMP.bas"; fi
ln -fsr $IDX_MT $IDX_MT_TMP
ln -fsr $BAM_WT $BAM_WT_TMP
if [ -f "$BAM_WT.bas" ]; then ln -fsr "$BAM_WT.bas" "$BAM_WT_TMP.bas"; fi
ln -fsr $IDX_WT $IDX_WT_TMP

export REF_PATH="$REF_BASE/genome.fa"
export REF_CACHE="$OUTPUT_DIR/ref_cache/%2s/%2s/%s"
PINDEL_WORKERS=1
NORM_CONTAM=''
declare -A do_parallel
source split_scripts/command-cgpPindel.sh
python3 split_scripts/pindel-input16.py "$PINDEL_MODE" --command "${do_parallel[cgpPindel]}" "${PINDEL_ARGS[@]}"
