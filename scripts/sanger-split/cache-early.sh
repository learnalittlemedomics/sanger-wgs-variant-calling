#!/bin/bash
set -euo pipefail
source run.params
mkdir -p ref_cache timings
NAME_MT=$(samtools view -H "$BAM_MT" | perl -ne 'chomp; if(/^@RG/){($sm)=/\tSM:([^\t]+)/; print "$sm\n";}' | uniq)
NAME_WT=$(samtools view -H "$BAM_WT" | perl -ne 'chomp; if(/^@RG/){($sm)=/\tSM:([^\t]+)/; print "$sm\n";}' | uniq)
if [[ "$BAM_MT" == *.cram ]]; then
 /usr/bin/time -v seq_cache_populate.pl -root "$OUTPUT_DIR/ref_cache" "$REF_BASE/genome.fa" > "timings/${PROTOCOL}_${NAME_MT}_vs_${NAME_WT}.time.cache_POP" 2>&1
else
 echo 'BAM input: no CRAM cache required' > "timings/${PROTOCOL}_${NAME_MT}_vs_${NAME_WT}.time.cache_POP"
fi
tar -czf cache_POP.stage.tar.gz timings
