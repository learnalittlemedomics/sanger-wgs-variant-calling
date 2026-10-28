# Sanger WGS Variant Calling — ChatGPT-Assisted Workflow Optimization

A personal fork of the [ICGC ARGO Sanger WGS Variant Calling Workflow](https://github.com/icgc-argo-workflows/sanger-wgs-variant-calling), with changes developed and refined with **ChatGPT assistance** to make the workflow easier to schedule, resume, configure, and maintain.

The upstream baseline is commit [817445b](https://github.com/icgc-argo-workflows/sanger-wgs-variant-calling/commit/817445baef7f6579f7240a082a8779cf2fec2b29), workflow version **2.1.0-9.7.1**. The underlying containerized tools—ASCAT, Pindel, CaVEMan, and BRASS—originate from [cancerit/dockstore-cgpwgs](https://github.com/cancerit/dockstore-cgpwgs).

ChatGPT assisted with code analysis, implementation, troubleshooting, and preparation of this public fork. The repository owner has run the modified workflow in their own environment. This is an independently maintained fork, not an official upstream release.

## Scope

The main changes concern **workflow execution and orchestration**: breaking a large caller task into smaller tasks, exposing useful parallelism, controlling resources, improving restart behavior, and fixing integration issues.

This fork does not claim to introduce new variant-calling algorithms, improve biological accuracy, or provide a measured speedup on every dataset. It continues to build on the upstream Sanger tools and result-processing modules.

## What changed and why

### 1. Split the original caller into independently scheduled tasks

**Where:** [`modules/sanger-split.nf`](modules/sanger-split.nf), [`sanger-wgs-variant-calling/main.nf`](sanger-wgs-variant-calling/main.nf), and [`scripts/sanger-split/`](scripts/sanger-split/).

The original workflow invokes the Sanger caller through one large module. This fork replaces that invocation with a split workflow containing separate preparation, reference-cache, caller, QC, and result-assembly tasks. Stage-specific helper scripts implement the individual operations.

**Practical benefits:**

- CPU, memory, and runtime requests can match individual stages instead of one allocation covering very different workloads.
- Nextflow exposes stage-level status, work directories, logs, and failures, making it easier to identify which operation needs attention.
- Successful stages can be reused with `-resume` when their cache entries, inputs, configuration, and required work files remain valid.
- Failures can be handled at a smaller task boundary, reducing the amount of work that may need to be repeated.
- Independent branches can run concurrently when their actual inputs are ready.

### 2. Add CaVEMan scatter/gather execution

**Where:** [`modules/sanger-split.nf`](modules/sanger-split.nf), [`caveman-scatter.py`](scripts/sanger-split/caveman-scatter.py), [`flag-scatter.py`](scripts/sanger-split/flag-scatter.py), and the CaVEMan command wrappers in [`scripts/sanger-split/`](scripts/sanger-split/).

CaVEMan is organized into setup, split/planning, grouped M-step execution, model merge, grouped E-step execution, gather, flagging, and optional annotation. The workflow explicitly waits for the merged model before launching dependent E-step tasks. Flagging also has planning, grouped execution, and collection stages.

Groups run bounded pools of workers rather than launching an unlimited number of processes. CPU allocations and worker counts are checked for consistency.

**Practical benefits:**

- Heavy computation can be distributed across multiple scheduler allocations.
- Grouping shards limits submission overhead compared with treating every small region as a separate scheduler job.
- A bounded worker pool gives more predictable concurrency within each allocation.
- Explicit model and gather dependencies make the execution order visible in the workflow.
- Individual failed groups can be retried without treating the complete CaVEMan computation as one indivisible task.

### 3. Separate Pindel input generation, chromosome work, and finalization

**Where:** [`modules/sanger-split.nf`](modules/sanger-split.nf), [`pindel-scatter.py`](scripts/sanger-split/pindel-scatter.py), [`pindel-input16.py`](scripts/sanger-split/pindel-input16.py), [`pindel-finalize.py`](scripts/sanger-split/pindel-finalize.py), and the Pindel wrappers and Perl helpers in [`scripts/sanger-split/`](scripts/sanger-split/).

Pindel execution is divided into input preparation, planning, chromosome-level tasks, and gather/finalization. The optimized dependency path also supports early input preparation.

**Practical benefits:**

- Chromosome-level work can run in parallel and has separate failure boundaries.
- Input generation and finalization are visible operations with their own resources and logs.
- Preparation can start earlier when the necessary inputs are available.
- Downstream CaVEMan flagging receives the required Pindel-derived dependency through explicit workflow connections.

### 4. Break BRASS into stages and parallel assembly groups

**Where:** [`modules/sanger-split.nf`](modules/sanger-split.nf), [`brass-scatter.py`](scripts/sanger-split/brass-scatter.py), [`brass-gather.py`](scripts/sanger-split/brass-gather.py), and the BRASS wrappers in [`scripts/sanger-split/`](scripts/sanger-split/).

The split workflow separates BRASS input processing, coverage, insert-size processing, grouping, copy-number normalization, planning, assembly, and final collection. Dependencies on ASCAT and other BRASS stages are expressed in the workflow graph.

**Practical benefits:**

- Preparation stages do not have to share one resource request with assembly.
- Ready branches can proceed without waiting for unrelated work.
- Assembly groups can run in parallel and be retried separately.
- The gather step centralizes collection of partial results before downstream packaging.

### 5. Configure resources by stage and handle selected resource failures

**Where:** [`conf/sanger-split.config`](conf/sanger-split.config), [`conf/sanger-scatter.config`](conf/sanger-scatter.config), and resource/retry helpers in [`modules/sanger-split.nf`](modules/sanger-split.nf).

The fork adds a per-stage resource map and scatter settings for worker counts, CPU steps, memory-per-CPU, planning estimates, and retry limits. CaVEMan retry logic recognizes selected termination codes and resource-failure messages; configured retry steps can reduce worker concurrency or adjust the requested allocation.

**Practical benefits:**

- Lightweight preparation and packaging tasks can request smaller allocations than expensive caller stages.
- Concurrency can be reduced on a retry to relieve memory pressure within an allocation.
- Resource choices are collected in configuration rather than spread across launch commands.
- Retry behavior distinguishes selected resource failures from ordinary errors that should terminate.

The included numbers are starting examples. Some implementation limits remain, including CPU and memory-per-CPU bounds in the CaVEMan helpers; this is not an unrestricted resource-autotuning system.

### 6. Start preparation earlier while retaining required dependencies

**Where:** the `sangerWgsSplit` workflow in [`modules/sanger-split.nf`](modules/sanger-split.nf), plus [`prepare-early.pl`](scripts/sanger-split/prepare-early.pl) and [`cache-early.sh`](scripts/sanger-split/cache-early.sh).

The optimized dependency path separates work that can start early from work that must wait for later reference, copy-number, or caller outputs.

**Practical benefits:**

- Preparation and cache work can overlap with other eligible tasks.
- The workflow avoids unnecessary waiting while preserving required caller dependencies.
- Dependencies can be inspected directly rather than inferred from a single large shell invocation.

Actual wall-clock savings depend on input sizes, queue availability, storage, and which stages dominate the run.

### 7. Make script changes affect only the relevant task inputs

**Where:** the `stageScripts` helper in [`modules/sanger-split.nf`](modules/sanger-split.nf).

Each stage receives the helper scripts it actually executes instead of hashing the entire helper-script directory as one shared input. Payload ID preparation also uses deep caching in the supplied split configuration.

**Practical benefits:**

- Editing one stage's script is less likely to invalidate unrelated stages.
- Iterative development can reuse more completed work when Nextflow's other cache requirements are met.
- Script dependencies are listed explicitly, making maintenance and review easier.

Cache reuse still depends on Nextflow's complete task signature and the availability of required outputs.

### 8. Select callers and adapt downstream processing

**Where:** caller-switch validation and downstream branches in [`sanger-wgs-variant-calling/main.nf`](sanger-wgs-variant-calling/main.nf), with defaults in [`conf/sanger-split.config`](conf/sanger-split.config).

The `params.sanger_steps` map controls ASCAT, Pindel, CaVEMan, BRASS, and QC. The workflow checks switch names and Boolean values, rejects an empty selection, and enforces required dependencies:

- CaVEMan requires ASCAT and Pindel.
- BRASS and tumour QC require ASCAT.

Result extraction, payload generation, supplementary outputs, upload branches, and cleanup inputs are assembled according to the enabled callers.

**Practical benefits:**

- Users can run supported subsets of the workflow.
- Invalid caller combinations fail before expensive work is submitted.
- Disabled caller branches do not require the same downstream processing as a full run.

### 9. Handle explicit scheduler submission rejections separately

**Where:** [`conf/sanger-scatter.config`](conf/sanger-scatter.config).

The configuration includes submission-rate and queue-size settings, plus retries for explicit quota or admission rejections. Submission retries are separate from task execution retries. The rejection expressions target known quota/admission messages rather than ambiguous transport timeouts.

**Practical benefits:**

- A temporary submission limit can be retried without consuming a task's resource-failure retry budget.
- Submission pressure can be tuned independently of per-task parallelism.
- Restricting retries to explicit rejection messages reduces the risk of blindly resubmitting after an uncertain transport failure.

These settings must be adapted to local scheduler policies and Nextflow behavior.

### 10. Fix publishing, cleanup staging, and payload identifiers

| Change | Location | Benefit |
| --- | --- | --- |
| Convert `publishDir.enabled` explicitly to a Boolean | [prep-sanger-qc](modules/raw.githubusercontent.com/icgc-argo/variant-calling-tools/prep-sanger-qc.0.1.3.0/tools/prep-sanger-qc/prep-sanger-qc.nf), [cleanup-workdir](wfpr_modules/github.com/icgc-argo/data-processing-utility-tools/cleanup-workdir@1.0.0/main.nf), and [payload-add-uniform-ids](wfpr_modules/github.com/icgc-argo/data-processing-utility-tools/payload-add-uniform-ids@0.1.1/main.nf) | Removes ambiguity from using a directory-string parameter as a publishing switch. |
| Stage cleanup inputs under separate numbered directories | [cleanup-workdir](wfpr_modules/github.com/icgc-argo/data-processing-utility-tools/cleanup-workdir@1.0.0/main.nf) | Avoids staging collisions when different upstream tasks produce files with identical names. |
| Use `workflow.sessionId` for the payload run identifier | [payload-gen-variant-calling](wfpr_modules/github.com/icgc-argo-workflows/data-processing-utility-tools/payload-gen-variant-calling@0.7.0/main.nf) | Uses a session UUID instead of the readable run name for this identifier. |

### 11. Separate public workflow code from private environment settings

**Where:** [`conf/`](conf/), [`sanger-wgs-variant-calling/nextflow.config`](sanger-wgs-variant-calling/nextflow.config), and the script-directory reference in [`main.nf`](sanger-wgs-variant-calling/main.nf).

The public migration removes local container-image paths, site queue overrides, account-specific quota annotations, and sample identifiers. Historical backups, installation notes, and runtime traces are excluded from the migrated changes. Shared configuration and helper-script references follow the repository layout.

**Practical benefits:**

- The repository can be shared without distributing the original deployment's private settings.
- Users can supply their own executor, queue, containers, and data paths through private configuration.
- Source changes are easier to review without historical backups and runtime artifacts mixed into the code.

### 12. Retain development checks and add resource reporting

**Where:** [`tests/scatter/`](tests/scatter/) and [`scripts/caveman_seff_report.py`](scripts/caveman_seff_report.py).

The repository retains mock integration checks for CaVEMan, Pindel, and BRASS, Nextflow stub workflows, and a fake scheduler. The reporting utility reads Nextflow trace files and collects Slurm `seff` information into job-level and summary outputs.

**Practical benefits:**

- Mock checks can exercise task coordination and error paths without real sequencing data or caller execution.
- Stub workflows provide optional checks of staging, dependencies, and retry behavior.
- Resource reports support evidence-based adjustment of CPU, memory, and grouping choices.

These checks are optional developer tools; they are not executed simply by reading or cloning the repository.

## Repository guide

| File or directory | Purpose |
| --- | --- |
| [`sanger-wgs-variant-calling/main.nf`](sanger-wgs-variant-calling/main.nf) | Main workflow entry point, caller selection, and downstream processing |
| [`sanger-wgs-variant-calling/nextflow.config`](sanger-wgs-variant-calling/nextflow.config) | Pipeline manifest and shared configuration include |
| [`modules/sanger-split.nf`](modules/sanger-split.nf) | Split processes, scatter/gather graph, resource helpers, and dependency handling |
| [`scripts/sanger-split/`](scripts/sanger-split/) | Caller wrappers, preparation helpers, scatter/gather drivers, and assembly scripts |
| [`conf/sanger-split.config`](conf/sanger-split.config) | Caller defaults and per-stage resource examples |
| [`conf/sanger-scatter.config`](conf/sanger-scatter.config) | Scatter settings, submission limits, and submission-retry configuration |
| [`tests/scatter/`](tests/scatter/) | Existing mock tests, stub workflows, fixtures, and fake scheduler |
| [`scripts/caveman_seff_report.py`](scripts/caveman_seff_report.py) | Trace and Slurm resource-report collection |
| [`docs/split-workflow.md`](docs/split-workflow.md) | Additional portable-configuration notes |

## Configuration

The pipeline declares **Nextflow >= 24.10.0**. The split module supports container execution through Singularity/Apptainer or Docker/Podman.

Before running, supply environment-specific settings in a private configuration file:

- Input alignment files and metadata, reference resources, and output locations.
- Executor and scheduler queue settings.
- Container engine configuration and compatible container images.
- Resource requests, scatter grouping, and submission limits appropriate to your environment.

For a local Sanger image, set `params.sanger_split_sif`; otherwise the split module has an upstream image fallback. Review `params.sanger_resources`, `params.sanger_scatter`, and `params.sanger_steps` before use.

Do not commit credentials, private deployment configuration, or sensitive sample metadata to a public repository. On a managed cluster, launch the Nextflow controller according to the site's scheduling policy.

## Validation and performance claims

The repository owner reports having run the modified workflow in their own environment. During migration, the existing source mock checks passed: four CaVEMan tests, one Pindel test, and one BRASS test. Python parsing and a scan of the prepared changes for known private identifiers and credential patterns also passed.

The public version includes sanitized configuration and adjusted relative paths. **Those publication-specific configuration and path changes have not been rerun end to end.**

Mock and stub checks validate orchestration behavior; they do not establish biological accuracy, numerical equivalence, or realistic memory consumption. This repository does not currently publish a controlled before/after benchmark for runtime, peak memory, compute cost, or variant-call concordance. The benefits described above explain the execution design; measured gains depend on the dataset and environment.

## Change history and acknowledgements

The migration is recorded in [PR #1](https://github.com/learnalittlemedomics/sanger-wgs-variant-calling/pull/1), with separate commits for integration fixes, the split workflow, and documentation/developer checks.

Credit for the original workflow belongs to [ICGC ARGO](https://github.com/icgc-argo-workflows/sanger-wgs-variant-calling), and credit for the underlying Sanger tools belongs to their respective developers, including [Cancer IT / Sanger](https://github.com/cancerit/dockstore-cgpwgs). ChatGPT assisted with this fork's workflow optimization and documentation; upstream attribution and applicable license terms remain in effect.
