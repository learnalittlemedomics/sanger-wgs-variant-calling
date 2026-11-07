# Verified local result publication

Local runs now end with `verifiedPublish`. It consumes the successful payload
channels, checks the complete `payload.files` inventory, publishes every file,
and emits a receipt. Cleanup depends on that receipt. A publication error fails
the workflow instead of permitting a successful completion marker.

## Filesystems

- Source symlinks are resolved before publication. Published results are regular
  files, never links to work-directory pathnames.
- On a filesystem that permits hard links, a temporary hard link is created in
  the destination directory and atomically renamed into place. Verification uses
  inode identity, payload size and a source MD5 check. Unlinking work files does
  not remove results. Each publication attempt reads the source once to verify
  its payload MD5, including large CRAM files; allow time for this I/O.
- On different filesystems, or when hard links are unsupported, bytes are copied
  into a unique destination-side temporary file. Source MD5, copied bytes, size,
  and a reread MD5 must match the payload before atomic replacement. Data and
  directory entries are fsynced; an unsupported durability operation fails closed.
- There is no assumption that a pathname prefix implies the same filesystem.
  An actual link operation determines whether the fast path is usable.
- Hard-linked results and work files share storage: treat them as immutable.
  To edit a result, create a new file and replace the name; do not edit in place.
- This implementation requires POSIX filesystem paths, Python 3, `flock`, hard
  links or normal file copies, and atomic rename within a destination directory.
  It is not an S3/object-storage publisher.

## Completion, recovery, and cleanup

`<publish_dir>/publication.verified.json` must have `state: VERIFIED` and the
Nextflow run must exit successfully. Individual files can appear before the
whole set is complete; consumers must use the run/receipt completion gate.
The receipt is invalidated at the beginning of each publication attempt and
records `FAILED` on an error. Missing payload groups, empty/missing files,
missing indexes, duplicate names, invalid sizes, and unsafe paths are rejected.
Existing identical files are reused; mismatching destinations are retained as
hidden `.previous` files before replacement. Temporary `.publish` files are not
completion evidence. A host crash may leave these hidden files for inspection.

The publication process disables Nextflow task caching, so `-resume` rechecks
results when source tasks remain available. It does not change caching of
analysis steps. An already running controller keeps the workflow loaded when
it started; this change applies to new or resumed controllers. Do not restart
healthy analysis solely to load this change.

Sanger verifies only enabled caller branches (ASCAT/CaVEMan/Pindel/BRASS), plus
supplement and QC. SONG/SCORE upload and remote-mode cleanup keep their original
dependencies. Local cleanup waits for successful publication. The final result
directory layout remains unchanged apart from the receipt and hidden lock.

`verifiedPublish` executes locally to the Nextflow controller without a
container. On a cluster, run the controller in a compute allocation with host
Python 3 and access to both work and output filesystems. On Dongfang, controllers
belong in the `long` partition. Publication is single-threaded I/O per workflow;
there is no new cross-sample concurrency limit.

## Tests

The tests cover resolved symlinks, same-filesystem hard links, real cross-device
copies, link-error fallback, truncated and same-size corrupt destinations,
checksum failures, missing indexes, receipt invalidation, workflow failure,
cleanup ordering, resume, and full/subset Sanger publication inventories.
Use the test scripts from a POSIX environment with Nextflow; they use synthetic
files and do not execute biological analysis.

Nextflow documents that `publishDir` is asynchronous:
[publishDir reference](https://docs.seqera.io/nextflow/reference/process/directives/publish-dir/).
That behavior motivates the explicit gate; it does not establish the underlying
cause of any previous storage short-copy incident.
