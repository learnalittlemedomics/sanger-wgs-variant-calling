# BRASS candidate-parallel remapping

This opt-in release parallelizes the per-candidate human-reference remapping
inside `filter_with_microbes_and_remapping.pl`. It preserves ordered intermediate
state, then runs the original microbial filter once on the complete state.
The step name does not imply that microbes were detected.

The three core files in `scripts/sanger-split/brass-remap-parallel/` match the
deployed release byte for byte. The root and nested workflow definitions,
shared launch scripts and default configuration are unchanged. The additional
configuration and seed utility are portable examples of deployment support;
they contain no site paths or sample identifiers.

## Execution and resource scope

- Linux/POSIX, Python 3.5 or later, and the original BRASS Perl environment,
  dependencies and reference databases are required. Only `ssearch36` is
  supported by this release.
- Consecutive candidates are divided into 256-record chunks. A 30-CPU allocation
  uses up to 15 local workers, each with two scoring threads. All workers run on
  one node. Increasing the allocation to 450 CPUs does not enable distributed
  execution or 30 threads per worker.
- `BRASS_REMAP_CHUNK_SIZE` changes chunk boundaries for validation; the production
  default is 256. A changed boundary invalidates checkpoint identity.
- The same original options, BAM, reference and database inputs must remain
  available throughout execution and recovery.

## Enable explicitly

Use an immutable release directory from this checkout. The Singularity overlay
mounts its entry point at the original helper path, and the directory at
`/opt/wtsi-cgp/brass-remap-parallel`. Invoke the existing workflow with these
additional arguments, retaining its normal inputs and configuration:

```bash
-c /absolute/checkout/conf/brass-remap-parallel.config \
--brass_remap_release_dir /absolute/checkout/scripts/sanger-split/brass-remap-parallel
```

The example targets the full process names `SangerWgs:sangerWgs:BRASS.*` and
`SangerWgs:sangerWgs:BRASSPlan`. It requests `long`, 168 hours, 30 CPUs and 108 GB
for the BRASS plan, retaining the earlier overlap-helper bind. Inspect the
effective configuration and generated submission script for the workflow
version and scheduler used at your site. Paths containing spaces, commas or
colons are rejected. Other container bind requirements must be preserved when
adapting this example.

The Nextflow controller allocation is separate: submit it through the existing
batch launcher in an activated Nextflow environment, and give it sufficient
wall time too. Where the scheduler permits shortening jobs through `time-min`,
verify both the actual time limit and minimum duration. This repository change
does not alter any running job or automatically resubmit samples.

## Ordered state and equivalence argument

Each candidate's remapping creates four arrays of scores/footprints, read
assignments and a FASTA stream. Workers return that intermediate state instead
of independently producing final filtered results. The controller joins arrays
and FASTA in original candidate order. For reads appearing in multiple chunks,
the later original candidate's assignment replaces the earlier assignment,
matching the serial hash's last-write behavior. Group rows retain their input
order and are routed to every chunk that needs the group, including repeated
group identifiers.

Only after this merge does one worker load the complete state and execute the
original microbial-filtering section against the full original input. Concatenating
independently filtered r5 files would not preserve this global behavior.

The worker adds state export/import boundaries. Removing those additions
reconstructs the original helper with SHA-256
`ff7e5d886547e2c52aded3c3f1d694139129f615a42ecdfbd9d2502142b99503`.
The public test checks that reconstruction and the exact three release hashes.
This source-level constraint supports, but does not replace, input/output tests.

## Checkpoints and failure handling

`<tmpdir>/parallel-v1/` holds the manifest, per-chunk receipts, progress and final
receipt. Only one controller can lock this directory at a time. Completed chunks
have hashed state/FASTA artifacts; failed or interrupted chunks have no valid
completion receipt and are run again. Input/group hashes, code hashes, arguments
and chunk size are part of the manifest identity. BAM, reference and database
identities use resolved path, size and modification time, **not full content
hashes**; keep those inputs immutable and retain independent checksums.

Worker failure propagates as a nonzero exit and terminates remaining children.
SIGTERM also terminates child process groups. A restart with the same paths,
inputs and options verifies completed chunks before reusing them. A changed
manifest identity or damaged receipt/output is rejected rather than silently
accepted. A complete-result cache verifies final output hashes; it does not
revalidate all intermediate chunks on that fast path.

`progress.json` reports completed records and chunks for the remapping phase.
It is not a percentage of the whole BRASS or Sanger workflow. The later global
microbial filtering and downstream workflow steps still need to finish.

## Optional reuse of completed overlap intermediates

`seed-intermediates.py` copies a separately verified BRASS directory into a new
destination via a staging directory. It requires the overlap completion marker
and an expected r4 SHA-256. Example arguments are deliberately generic:

```bash
python3 scripts/sanger-split/brass-remap-parallel/seed-intermediates.py \
  --source /analysis/previous/brass \
  --destination /analysis/new-work/brass \
  --r4-name sample.r4 \
  --expected-r4-sha256 <independently-recorded-sha256>
```

Before using this utility, stop all writers to the source/destination, establish
exclusive ownership of the destination, and independently verify that the
sample, inputs, references, versions and analysis options match. The utility
checks the named r4 and marker, not every copied artifact or provenance fact.
An existing destination is checked for those same two conditions and reported
as reused; it is not a full-directory comparison. Symlinks are preserved, so
their targets must remain available. An interrupted staging copy is retained
for inspection. This helper does not itself cancel jobs, modify Nextflow cache
metadata, or authorize reuse of incomplete upstream work.

## Validation and measured performance

Private-input results are published only as aggregate measurements in
[`benchmarks/brass-remap-parallel-2026-11-06.json`](benchmarks/brass-remap-parallel-2026-11-06.json).
No controlled sequencing inputs, sample identifiers, coordinates or read
sequences are included in the fixtures or evidence.

| Comparison | Original | Parallel | Evidence and scope |
|---|---:|---:|---|
| 48 scoring calls, 16 inputs, median of three randomized rounds | 9.108 s, 1 × 30 | 1.502 s, 15 × 2 | 6.06×; identical scores; scoring calls only |
| Remapping state, 40 real candidates | 107.242 s | 61.604 s, 15 × 2 | 1.74×; four arrays, read map and FASTA equal |
| Full filter, 40 real candidates | 315.461 s | 268.843 s, 8 × 3 | 1.17×; r5 and scores byte-identical |
| Shared-read boundary, 16 records, 63 reads crossing chunks | 9.871 s | 10.196 s, 2 × 2 | No speedup in this small case; state and FASTA equal |

Both forty-record final r5 files were empty. Equality of empty r5 alone is
insufficient: the comparison also checked 2,819 bytes of scores, nonempty
FASTA, all four arrays and the global read map. The state test's real-candidate
FASTA contained 872,453 bytes. A separate real-input interruption test sent
SIGTERM after one completed chunk, resumed, verified receipt reuse and checked
the resulting state/output against the serial reference. Corrupted checkpoints
were rejected and child errors propagated.

The entire production analysis was still running when these results were
recorded. The measurements are bounded tests, not a full-input equivalence
claim, whole-workflow speedup or completion-time guarantee. Node load, caching,
candidate distribution and the global filtering phase affect elapsed time.

## Public reproducible checks

```bash
python3 tests/scatter/test_brass_remap_parallel.py
```

The dedicated GitHub Actions workflow runs the checks on Linux. They cover
ordered merge, shared-read last ownership, repeated-group routing, actual
interrupt/resume, cached results, altered inputs, corruption, child failure,
empty input and the seed utility. A synthetic worker implements only the state
protocol: these tests validate orchestration and release scope, **not biological
alignment correctness**. POSIX orchestration tests are skipped on Windows.

For new biological validation, use an independent original helper and the same
container/inputs/options on bounded representative candidate sets. Compare
remapping arrays, global read ownership, concatenated FASTA, score output and
final r5; include cross-chunk shared reads and a case with nonempty retained
variants. Preserve original candidate order and compare full state before the
single global filter. Do not infer equality from runtime or empty r5 alone.
