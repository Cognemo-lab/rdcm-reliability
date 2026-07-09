# rDCM Pipeline

A MATLAB pipeline for whole-brain effective connectivity analysis using Regression Dynamic Causal Modelling (rDCM), built on the [TAPAS toolbox](https://github.com/translationalneuromodeling/tapas). The pipeline handles multi-subject resting-state fMRI data, runs rDCM estimation in parallel on a compute cluster, computes split-half reliability, performs group-level quality control, and exports per-subject effective connectivity (EC) matrices with accompanying visual QC reports.

---

## Requirements

- MATLAB R2020a or later (`exportgraphics` required for multi-page PDF output)
- [TAPAS toolbox](https://github.com/translationalneuromodeling/tapas) (rDCM module)
- SPM12
- A compute cluster with SLURM (or equivalent) — or a local workstation for small datasets

---

## Overview

The pipeline is split into two phases:

| Phase | Script | Where to run |
|---|---|---|
| **Phase 1 — Setup** | `rdcm_setup.m` | Local workstation (interactive) |
| **Phase 2 — Estimation** | `run_rdcm_pipeline.m` | Cluster or local |

**Phase 1** converts your BOLD data, lets you define ROI subsets interactively, builds a job grid, and freezes all parameters into a `.mat` file. **Phase 2** loads those frozen files and runs estimation, QC, and export non-interactively.

---

## Quick Start

### 1. Configure

Copy `set_params_TEMPLATE.m` to `set_params.m` and fill in your paths and settings. Add `set_params.m` to your `.gitignore` — it contains local paths and should not be committed.

```matlab
% Key fields to set:
params.dirs.root          = '/path/to/your/workspace';
params.filenames.roi_list = 'your_roi_labels.csv';
params.convert_input_fn   = 'rdcm_convert_input_bids';  % or your custom converter
params.params_input_fn    = 'set_params_input_bids';
params.rdcm.dt            = 0.8;   % TR in seconds
params.rdcm.est_method    = 1;     % 1 = ridge, 2 = sparse
```

### 2. Run Phase 1 (local)

```matlab
rdcm_setup('/path/to/your/workspace')
```

This produces four `.mat` files in `params.dirs.output`:

| File | Contents |
|---|---|
| `rdcm_input.mat` | BOLD data struct (`Yall`, subject IDs, ROI names) |
| `roi_subsets.mat` | ROI subset definitions from your interactive selection |
| `rdcm_grid.mat` | Job list with `run_id` hash |
| `rdcm_params.mat` | Frozen parameter snapshot |

### 3. Transfer to cluster

```bash
rsync -av /local/rdcm-output/ cluster:/scratch/youruser/rdcm-output/
```

Also transfer the pipeline code, or ensure it is on the cluster's MATLAB path.

### 4. Configure cluster paths

Copy `set_params_pipeline_TEMPLATE.m` to `set_params_pipeline.m` in the output directory on the cluster, and set the output path:

```matlab
params.dirs.output = '/scratch/youruser/rdcm-output';
params.dirs.logs   = fullfile(params.dirs.output, 'logs');
params.log_mode    = 'file';
```

### 5. Submit Phase 2

```bash
matlab -nodisplay -nosplash \
  -r "run_rdcm_pipeline('/scratch/youruser/rdcm-output'); exit"
```

On wall-time expiry, resubmit the same command unchanged. The checkpoint system automatically skips completed jobs.

---

## Input Converters

The pipeline uses pluggable input converters to accommodate different data formats. Set `params.convert_input_fn` to one of the shipped converters or your own custom function.

| Converter | Format | Param file |
|---|---|---|
| `rdcm_convert_input_bids` | BIDS-formatted NIfTI + confounds TSV | `set_params_input_bids` |
| `rdcm_convert_input_onesubmat` | One `.mat` file per subject | `set_params_input_onesubmat` |

**Custom converter signature:**

```matlab
function data = my_convert_fn(params)
% Returns data struct with fields:
%   .Yall          {1 x N_subs} cell of Y structs (Y.y, Y.dt, Y.name, Y.subj)
%   .all_roi_names {1 x R} cell of ROI name strings
%   .subject_ids   {1 x N_subs} cell of subject ID strings
%   .subs          struct with subject metadata
%   .ROIs          struct with ROI metadata
```

---

## ROI Subset Definition

`rdcm_setup` calls `rdcm_define_roisets` interactively, presenting a list of all ROIs in your dataset. You select one or more subsets (e.g., a default-mode network, a salience network) by name or network label. Each subset is independently estimated, QC'd, and exported. Subset definitions are saved in `roi_subsets.mat` and frozen with the rest of the pipeline configuration.

---

## Pipeline Architecture

```
rdcm_setup
├── rdcm_convert_input_*     ← format-specific BOLD converter
├── rdcm_diagnose_input      ← flag/trim NaN ROIs or subjects
├── rdcm_define_roisets      ← interactive ROI subset selection
└── rdcm_build_grid          ← build (roi_subset × subject) job list

run_rdcm_pipeline
└── rdcm_pipeline
    ├── rdcm_checkpoint            ← resume support
    ├── rdcm_run_single            ← per-subject rDCM via TAPAS
    ├── rdcm_fitmetrics_single     ← cosine similarity, R²
    ├── rdcm_splithalf_single      ← split-half ICC reliability
    ├── rdcm_qc_group              ← group-level QC tiers (per ROI set)
    ├── [cross-ROI-set advisory]   ← Adv E: tier discordance across ROI sets
    └── [reporting + export]
        ├── rdcm_report_qc              ← group QC summary → TXT
        ├── rdcm_report_subject         ← per-subject roster / deep-dive → TXT
        ├── rdcm_plot_qc_distributions  ← group QC distribution plots → PDF
        ├── rdcm_export_manifest        ← run manifest → TXT/JSON
        ├── rdcm_export_ec              ← EC matrices → CSV
        └── rdcm_visual_qc              ← per-subject PDF figures
```

### Results struct

`rdcm_pipeline` returns a `results` struct with the following layout:

```
results.grid        {N_roisets × N_subs}   EC output per subject per ROI set
results.fit         {N_roisets × N_subs}   fit metrics (cosine similarity, R²)
results.sh          {N_roisets × N_subs}   split-half reliability (ICC, Pearson r)
results.qc          {N_roisets × 1}        group QC summary per ROI set
results.params                              params as used (with run_id, log paths)
results.data_meta                           data.subs and data.ROIs metadata
```

---

## Quality Control

QC is applied per subject within each ROI subset using split-half reliability (ICC), not raw fit alone. The decision flow runs in stages:

**Stage 1 — Hard gates** (automatic exclusion, Tier 0):
- Gate 1a: `cosine_sim_mean < 0`
- Gate 1b: `min(ICC_diag_masked, ICC_offdiag_masked) < params.qc.icc_floor`
- Gate 1c: mean framewise displacement or scrubbed-volume proportion exceeds `params.qc.fd_threshold` / `fd_prop_threshold` (if FD data supplied)

**Stage 2 — Soft tier** (subject-level reliability, ICC_s):

| Tier | Label | Criteria |
|---|---|---|
| 3 | **Pass** | `icc_s_strong ≥ params.qc.icc_s_pass` — usable for group + individual analysis |
| 2 | **Marginal** | Below `icc_s_pass` but not hard-gated — usable for group analysis only |
| 1 | **Fail** | `icc_s_strong < params.qc.icc_floor` |
| 0 | **Excluded (hard-gate)** | Gate 1a/1b/1c triggered |

**Stage 3 — Advisory flags** (never cause exclusion, but accumulate toward `review_recommended`):
- Adv A: `logF` anomalously low (> 2 SD below group mean)
- Adv B: reliability below group median − `k × MAD` (`params.qc.mad_k`)
- Adv C: bias index (`icc_offdiag_masked − pearson_r_offdiag`) below `params.qc.bias_threshold`
- Adv D: EC matrix Mahalanobis outlier (robust PCA)
- Adv E: cross-ROI-set tier discordance, set by `rdcm_pipeline` after all ROI sets complete

Two subjects can differ in classification purely because reliability, not fit, drives inclusion — a subject with excellent cosine similarity can still be Marginal if split-half ICC is low.

QC outputs:
export/
├── ec_<label>.csv               ← EC matrix (subjects × connections)
├── ec_<label>_labels.txt        ← connection labels (ROIj->ROIi, column-major)
├── ec_<label>_meta.txt          ← provenance: run_id, N_subs, convention
├── qc_report_<run_id>.txt       ← group QC text report
├── qc_report_subject_<run_id>.txt ← per-subject roster (if enabled)
├── qc_distributions_<run_id>.pdf  ← group-level QC distribution plots
├── manifest_<run_id>.txt        ← run manifest (provenance summary)
└── qc_figures/
    └── qc_<sub_id>_<run_id>.pdf

Each visual QC page shows the full-session EC matrix (network-sorted, RdBu colormap), fit and reliability metrics, split-half EC thumbnail pair (A₁ vs A₂), and a raw self-connection panel that flags any self-connection with the wrong sign (expected: negative).

---

## Outputs

All outputs are written to `params.dirs.output/export/`:

```
export/
├── ec_<label>.csv            ← EC matrix (subjects × connections)
├── ec_<label>_labels.txt     ← connection labels (ROIj->ROIi, column-major)
├── ec_<label>_meta.txt       ← provenance: run_id, N_subs, convention
├── qc_report_<run_id>.txt    ← group QC text report
└── qc_figures/
    └── qc_<sub_id>_<run_id>.pdf
```

**EC matrix convention:** Columns are ordered column-major, source-grouped. Column `k` corresponds to label `k`: `ROIj->ROIi` where `A(i,j)` = source `j` → target `i`. This matches the vectorisation `A(:)'` applied to TAPAS output.

---

## Key Parameters

| Parameter | Default | Description |
|---|---|---|
| `params.rdcm.dt` | *(required)* | TR in seconds |
| `params.rdcm.est_method` | `1` | `1` = ridge regression, `2` = sparse rDCM |
| `params.rdcm.SNR` | `1` | SNR passed to TAPAS (no effect on empirical resting-state estimation) |
| `params.splithalf.min_timepoints` | `50` | Minimum timepoints per half for reliability |
| `params.splithalf.exclude_diagonal` | `false` | Exclude self-connections from ICC calculation |
| `params.qc.icc_floor` | `0.40` | Gate 1b hard-gate floor + Tier 1 threshold |
| `params.qc.icc_s_pass` | `0.60` | ICC_s threshold for Tier 3 (Pass) |
| `params.qc.icc_c_fail` | `0.40` | Per-connection ICC Fail threshold |
| `params.qc.icc_c_pass` | `0.60` | Per-connection ICC Pass threshold |
| `params.qc.bias_threshold` | `-0.10` | Advisory C bias index threshold |
| `params.qc.outlier_z` | `3.5` | Robust z-score flag on KL-masked ICC_s_all |
| `params.qc.mad_k` | `2` | Advisory B median-MAD multiplier |
| `params.qc.kl_mask_method` | `'top_k'` | `'top_k'` or `'top_quartile'` connection masking |
| `params.qc.kl_mask_k` | `1000` | Connections retained when `kl_mask_method = 'top_k'` |
| `params.qc.kl_mask_p` | `0.75` | Quantile retained when `kl_mask_method = 'top_quartile'` |
| `params.qc.kl_mask_min_connections` | `20` | Minimum connections retained by KL mask |
| `params.qc.fd_threshold` | `Inf` | Gate 1c mean framewise displacement (mm) |
| `params.qc.fd_prop_threshold` | `Inf` | Gate 1c proportion of scrubbed volumes |
| `params.report.skip_plots` | `false` | Skip figure generation (for headless cluster nodes) |
| `params.report.subject_roster` | `true` | Write per-subject roster via `rdcm_report_subject` |
| `params.report.deep_dive_ids` | `{}` | Subject IDs for single-subject deep-dive reports || `params.diagnose_ROIs` | `'exclude_ROI'` | `'exclude_ROI'` or `'exclude_sub'` for bad ROI handling |
| `params.include_subs` | `{}` | Subject ID whitelist; `{}` = include all subjects |
| `params.verbose` | `2` | Log verbosity level (0–3) |
| `params.log_mode` | `'both'` | `'console'`, `'file'`, or `'both'` |

---

## Cluster Resumption

The checkpoint system (`rdcm_checkpoint.m`) writes a log of completed job tokens (`rs<i>_sub<j>`) to `logs/checkpoint_<run_id>.log`. On resubmission, completed tokens are skipped automatically. To force a clean restart, delete the checkpoint file:

```bash
rm /scratch/youruser/rdcm-output/logs/checkpoint_<run_id>.log
```

The `run_id` is a hash derived from ROI subset labels and estimation method. Changing either automatically generates a new `run_id` and a fresh checkpoint, preventing stale checkpoints from silently skipping jobs in a new configuration.

---

## File Index

| File | Role |
|---|---|
| `set_params_TEMPLATE.m` | Parameter template — copy and fill in as `set_params.m` |
| `set_params_pipeline_TEMPLATE.m` | Cluster path override template |
| `rdcm_setup.m` | Phase 1 setup script (run locally, interactive) |
| `run_rdcm_pipeline.m` | Phase 2 runner (run on cluster or locally) |
| `rdcm_pipeline.m` | Main estimation loop |
| `rdcm_build_grid.m` | Job grid constructor |
| `rdcm_run_single.m` | Single-subject rDCM wrapper |
| `rdcm_fitmetrics_single.m` | Fit metric computation (cosine similarity, R²) |
| `rdcm_splithalf_single.m` | Split-half reliability (ICC, Pearson r) |
| `rdcm_qc_group.m` | Group-level QC tiering |
| `rdcm_checkpoint.m` | Checkpoint read/write for cluster resumption |
| `rdcm_export_ec.m` | EC matrix export to CSV |
| `rdcm_export_csv.m` | Generic CSV writer utility used by export functions |
| `rdcm_export_manifest.m` | Run manifest writer (provenance summary) |
| `rdcm_report_qc.m` | Group QC text report writer |
| `rdcm_report_subject.m` | Per-subject roster and deep-dive report writer |
| `rdcm_plot_qc_distributions.m` | Group-level QC distribution plots |
| `rdcm_visual_qc.m` | Per-subject visual QC PDF generator |
| `rdcm_define_roisets.m` | Interactive ROI subset selector |
| `rdcm_build_input.m` | Internal data struct constructor |
| `rdcm_diagnose_input.m` | NaN/missing ROI diagnostic and trimmer |
| `rdcm_connection_labels.m` | Connection label generator |
| `rdcm_cosine_sim.m` | Cosine similarity utility |
| `rdcm_extract_subid.m` | Subject ID extraction utility |
| `rdcm_log.m` | Logging utility |
| `rdcm_convert_input_bids.m` | BIDS input converter |
| `rdcm_convert_input_onesubmat.m` | Per-subject `.mat` input converter |
| `setup_paths.m` | MATLAB path configuration |

---

## Citation

If you use this pipeline, please cite the TAPAS rDCM toolbox papers:

- Frässle, S. et al. (2017). Regression DCM for fMRI. *NeuroImage*, 155, 406–421.
- Frässle, S. et al. (2018). A generative model of whole-brain effective connectivity. *NeuroImage*, 179, 505–529.
- Frässle, S. et al. (2020). Regression dynamic causal modeling for resting-state fMRI. *bioRxiv*. https://doi.org/10.1101/2020.08.12.247536

---

## License

This pipeline is released under the GNU General Public License v3.0, consistent with the TAPAS toolbox license terms.
