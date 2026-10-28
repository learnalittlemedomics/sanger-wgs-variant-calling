# Split Sanger workflow

The split workflow separates preparation, ASCAT, Pindel, CaVEMan, BRASS and QC into resumable tasks. CaVEMan, Pindel and BRASS include scatter/gather execution.

Enable callers with `params.sanger_steps`. CaVEMan requires ASCAT and Pindel; BRASS and tumour QC require ASCAT. Resource and scatter settings are in `conf/sanger-split.config` and `conf/sanger-scatter.config`. Values are examples, not measured guarantees for other datasets.

Set `params.sanger_split_sif` in a private configuration to a compatible Sanger container image. Set the executor, queues, container engine and reference/input/output paths in that private configuration. No site queue, account, credentials or local image path is distributed.

Run Python mock tests on Linux:

```sh
python3 tests/scatter/test_scatter.py
python3 tests/scatter/test_pindel_scatter.py
python3 tests/scatter/test_brass_scatter.py
```

`python3 tests/scatter/preflight.py` additionally requires Nextflow and exercises stub processes and a fake scheduler. It does not validate real variant calls or replace end-to-end validation with representative data and the intended container.

Before production use, validate the sanitized configuration on your environment. On managed clusters, submit the Nextflow controller according to local policy.
