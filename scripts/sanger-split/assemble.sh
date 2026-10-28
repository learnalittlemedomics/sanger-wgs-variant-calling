#!/bin/bash
set -euo pipefail
for archive in "$@"; do tar -xzf "$archive"; done
source run.params
NAME_MT=$(samtools view -H "$BAM_MT" | perl -ne 'chomp; if($_ =~ m/^\@RG/) {($sm) = $_ =~m/	SM:([^	]+)/; print "$sm\n";}' | uniq)
NAME_WT=$(samtools view -H "$BAM_WT" | perl -ne 'chomp; if($_ =~ m/^\@RG/) {($sm) = $_ =~m/	SM:([^	]+)/; print "$sm\n";}' | uniq)
PAIR=${PROTOCOL}_${NAME_MT}_vs_${NAME_WT}
mkdir -p "$PAIR" timings
rm -rf "$PAIR"/*/logs
for role in WT MT; do
  sample=$NAME_WT
  [ "$role" = MT ] && sample=$NAME_MT
  if [ -f "timings/$PAIR.time.verify_$role" ]; then
    mv "timings/$PAIR.time.verify_$role" "timings/${PROTOCOL}_$sample.time.verify_$role"
  fi
done
results=("$PAIR")
[ ! -d "${PROTOCOL}_${NAME_MT}" ] || results+=("${PROTOCOL}_${NAME_MT}")
[ ! -d "${PROTOCOL}_${NAME_WT}" ] || results+=("${PROTOCOL}_${NAME_WT}")
tar -czf "$PAIR.result.tar.gz" "${results[@]}"
tar -czf "$PAIR.timings.tar.gz" timings
printf '%s\n' 'Split Nextflow workflow: see per-stage timings and Nextflow trace; no single-process wall time.' > "$PAIR.time"
