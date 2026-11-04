# BRASS overlap enumeration: equivalence and validation

## Scope

`get_abs_bkpts_from_clipped_reads.pl` previously traversed every `i < j` pair,
deep-cloned both records, and immediately skipped pairs whose two chromosome
names differed. For a large candidate list this spends most of its time on
pairs that cannot affect the result.

The indexed version changes **only pair enumeration**. It builds an ordered
index keyed by the two chromosome names and visits only later records in the
same group. It keeps the original global `i` order, ascending `j` order,
deep copies, direction checks, clipping thresholds, overlap tests, read
reassignment, and output code. No parallel reduction or coordinate pruning
is introduced. Nested hashes avoid delimiter collisions in contig names.

The existing BRASSPlan script staging and container bind already select this
helper. No Nextflow graph, resource settings, or other caller is changed.

## Why equivalence does not require a multi-day whole-input reference run

For valid BRASS candidate records, let `K(i)` be the ordered pair of chromosome
names in columns 1 and 4. The original loop processes the sequence
`(0,1), (0,2), ..., (1,2), ...`.

1. If `K(i) != K(j)`, the original body creates local copies and immediately
   executes `next`. It does not mutate the input records, read collections,
   counters used by the analysis, or output.
2. The index stores every input index once in its group, in increasing order.
3. For each original outer index `i`, the new loop visits the suffix after `i`
   in that group. This is exactly the original inner sequence with only the
   guaranteed no-op pairs removed.
4. The remaining body is unchanged. By induction over the retained sequence,
   the mutable read collections and output state after each retained pair
   are the same.

The chromosome-name fields in `@lines` are never rewritten by the loop;
working record edits use deep copies. Group membership therefore cannot
become stale as the read collections change.

Preserving sequence matters: an earlier overlapping pair can reassign reads
that a later pair uses. Running groups or candidates in arbitrary order,
changing overlap thresholds, or naively parallelizing shared candidates is
not justified by this argument.

This is a semantic argument for this narrowly scoped change, **not** a claim
that all BRASS logic is correct or that an entire variant-calling workflow has
been revalidated.

## Validation layers

### Actual Perl enumeration versus a bounded exhaustive oracle

Use Python 3.7+ and Perl with its core `JSON::PP` module. The enumeration test
does not require the alignment libraries or access to sequencing data. Run:

```bash
python3 tests/scatter/test_brass_overlap_index.py
```

The test extracts the actual production index/loop headers and substitutes a
pair trace for the analysis body. It compares that trace exactly against a
separate brute-force oracle. Cases include empty/singleton input, one group,
interleaved groups, reversed chromosome pairs, delimiter-like contig names,
all binary group assignments through length six, and seeded randomized data.
It also reconstructs the original loop headers and checks the complete
script fingerprint against the validated baseline. Intentional future changes
outside the index require a new review and an updated baseline.

### Real read-level differential tests

Run the unchanged original helper and the indexed helper against the same
bounded `.r3` fixture, alignment, reference, software environment and cache:

```bash
git show 13f25b1a4aa6767297ca6ae5c5f181867d1c830b:scripts/sanger-split/get_abs_bkpts_from_clipped_reads.pl > /tmp/original-helper.pl
python3 tests/scatter/benchmark_brass_overlap.py \
  --original /tmp/original-helper.pl \
  --optimized scripts/sanger-split/get_abs_bkpts_from_clipped_reads.pl \
  --alignment tumour.cram --reference genome.fa \
  --candidates bounded-fixture.r3 --output comparison --repeats 2
```

The tool requires a new output directory, compares success codes, exact `.r4`
bytes by SHA256, and diagnostics after normalizing only script path/line
suffixes. It records both total elapsed time and the overlap-loop time from
existing phase markers. Repeated runs alternate execution order to reduce
systematic page-cache bias. Use an appropriate compute allocation; do not run
large reference comparisons on a login node.

Useful fixtures include a seeded mixture of groups, a complete chromosome-pair
group, dense candidates, intrachromosomal candidates, high-support candidates,
and controlled near-neighbour perturbations with distinct IDs to exercise
shared-read changes. Subsetting creates a test input; matching on a subset
does **not** mean its output is the corresponding slice of a whole-input run.

### Honest performance interpretation

On one inspected input, 257,020 candidates implied 33,029,511,690 unordered
pairs. Same-chromosome-pair groups contained 185,852,245 pairs, so 99.437% of
original iterations were guaranteed immediate skips. The resulting 177.7-fold
reduction in visited pairs is **not** a measured wall-clock speedup.

Reading alignments, collecting supporting reads, retained comparisons,
output generation and later BRASS stages remain. Report measured overlap-loop
speed separately from whole-helper speed, including any cases that regress.
Finite real-data tests cannot establish equivalence for every possible input.

Raw alignments, candidate coordinates, read names, patient identifiers and
private deployment paths must not be published with benchmark results.

## Measured validation (2026-11-05)

All eight paired real-input fixtures exited successfully with **byte-identical
`.r4` outputs and identical raw stderr**, including the shared-read fixture
that triggered 42 read-separation warnings in each version. The complete
largest chromosome-pair group retained every pair. The four test methods
(137 ordered-pair fixtures plus the baseline-scope guard) also passed.
The reusable benchmark harness passed four empty-input executions in the
same BRASS container, including its Python 3.5 runtime.

| Fixture | Rows | Original loop (s) | Indexed loop (s) | Loop speedup | Original total (s) | Indexed total (s) |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| Mixed groups | 2500 | 73.315 | 4.326 | 16.95x | 630.408 | 447.460 |
| Complete largest chromosome-pair group | 3472 | 1304.700 | 1301.457 | 1.00x | 1883.475 | 1932.145 |
| Dense prefix | 2000 | 60.194 | 26.597 | 2.26x | 496.808 | 405.543 |
| High support | 500 | 5.000 | 2.286 | 2.19x | 321.340 | 357.128 |
| Intrachromosomal | 1200 | 24.619 | 11.174 | 2.20x | 265.701 | 209.068 |
| Shared-read perturbations | 48 | 10.528 | 10.939 | 0.96x | 20.262 | 20.894 |
| Empty | 0 | <0.01 | <0.01 | — | 0.756 | 0.720 |
| Singleton | 1 | <0.01 | <0.01 | — | 1.070 | 1.044 |

These are single paired observations per fixture, not repeated performance
estimates. Up to four fixtures ran concurrently in one 4-CPU allocation;
variant order alternated between fixtures. Shared-storage cache state and
concurrency were not controlled, so total elapsed differences can include I/O
and system effects. Values below 0.01 seconds are not meaningful speed tests.
Cases with little or no improvement are retained in the table.

This validates the helper on bounded inputs. It does not claim that the
257,020-row original reference run completed, that the full optimized input
has completed, or that a whole-pipeline speedup has been measured. The full
input equivalence argument remains the order-preserving no-op removal above.

The [aggregate benchmark record](benchmarks/brass-overlap-index-2026-11-05.json)
contains script fingerprints and timings. Raw inputs and outputs remain private.

## Deployment boundary

Changing this repository does not hot-patch already running Perl processes or
update an independently installed server checkout. Validate in a separate
directory first. Preserve the original process and its inputs when comparing;
deploy to an active environment only with the appropriate authorization.
