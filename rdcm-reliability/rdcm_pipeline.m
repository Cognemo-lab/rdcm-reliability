function results = rdcm_pipeline(params, data, roi_subsets, grid)
%% rdcm_pipeline
%{
Outer pipeline loop for the rDCM effective connectivity + reliability
analysis. This is the main entry point for both local and cluster runs.

The pipeline iterates over a 2-level grid:
  Level 1: ROI subsets
  Level 2: Subjects (data.Yall)

For each (ROI subset, subject) combination, the pipeline:
  1. Extracts the subject BOLD timeseries for the current ROI subset
  2. Runs full-session rDCM via rdcm_run_single [RESTING-STATE]
  3. Computes model fit metrics via rdcm_fitmetrics_single
  4. Runs split-half reliability via rdcm_splithalf_single
  5. Saves per-cell results to disk (params.dirs.output/cells/)
  6. Marks the grid index as complete in the checkpoint file

After the full grid is complete:
  Section 7 — Group QC       : rdcm_qc_group per ROI set
  Section 8 — Cross-net flags : Adv E (cross-ROI-set tier discordance)
  Section 9 — Reporting       : rdcm_report_qc, rdcm_report_subject,
                                rdcm_plot_qc_distributions
  Section 10 — Export manifest: rdcm_export_manifest
  Section 11 — Save results   : results_<run_id>.mat

NOTE — SNR: For resting-state rDCM, SNR is only used when generating
synthetic data and has no effect on empirical estimation. It is set to
a fixed value (params.rdcm.SNR = 1) and is not swept as a grid axis.
See tapas_rdcm_set_options documentation and Frässle et al. (2021).

FUTURE (task_mode): For task-based rDCM, add a driving-input extraction
step before rdcm_run_single and pass U alongside Y. SNR or other
hyperparameters may become meaningful grid axes at that point.

Cluster / wall-time resumption:
  On startup, rdcm_checkpoint checks whether a checkpoint file exists
  for this run. Completed indices are loaded and skipped. Re-submit
  naively after wall-time expiry. To force a clean restart, delete the
  checkpoint file (params.log.ck_file).
---------------------------------------------------------------------------
INPUTS
---------------------------------------------------------------------------
params := struct — pipeline configuration. Key fields:

  Logging:
    .verbose      integer 0-3 (default: 2)
    .log_mode     'console'|'file'|'both' (default: 'both')

  Directories:
    .dirs.output  char — root output directory (created if absent)
    .dirs.logs    char — log/checkpoint directory

  rDCM settings:
    .rdcm.SNR         scalar  — fixed SNR value (default: 1, not swept)
    .rdcm.est_method  integer 1=ridge / 2=sparse (default: 1)
    .rdcm.dt          scalar  — TR in seconds (must match data)

  Split-half:
    .splithalf.exclude_diagonal    logical (default: false)
    .splithalf.min_timepoints      integer (default: 50)
    .splithalf.store_half_fits     logical (default: false)
                                   When true, H1/H2 fit metrics are stored
                                   in results.fit{i_rs,i_sub}.sh_h1/.sh_h2
                                   for rdcm_report_subject Section 6.

  QC thresholds:
    .qc.icc_floor scalar (default: 0.40) — applied to min(diag_masked, offdiag_masked)
    .qc.icc_s_pass scalar (default: 0.60)
    .qc.icc_c_fail scalar (default: 0.40)
    .qc.icc_c_pass scalar (default: 0.60)
    .qc.bias_threshold scalar (default: -0.10)
    .qc.outlier_z scalar (default: 3.5) — robust z-score flag on KL-masked ICC_s_all
    .qc.kl_mask_method 'top_k'|'top_quartile' (default: 'top_k')
    .qc.kl_mask_k integer (default: 1000)
    .qc.kl_mask_p scalar in [0,1] (default: 0.75)
    .qc.kl_mask_min_connections integer (default: 20)
    .qc.kl_mask_obs_param_ratio scalar (default: 10) — sizes mask via N_TRs_half

  Reporting:
    .report.skip_plots     logical (default: false)
                           Set true on headless cluster nodes without display.
    .report.subject_roster logical (default: true)
                           Write per-subject roster via rdcm_report_subject.
    .report.deep_dive_ids  cell of char (default: {})
                           Subject IDs to run single-subject deep-dive reports.

data := struct — BOLD timeseries data loaded from rdcm_input.mat
  .Yall         {1 x N_subs} cell of subject Y structs
  .subs         struct with subject metadata
  .ROIs         struct with ROI metadata
  .subject_ids  {1 x N_subs} cell of subject ID strings

roi_subsets := struct array from rdcm_define_roisets
grid        := struct from rdcm_build_grid
---------------------------------------------------------------------------
OUTPUTS
---------------------------------------------------------------------------
results := struct with fields:
  .grid     cell(N_roisets, N_subs)  — EC matrices per cell
  .fit      cell(N_roisets, N_subs)  — fit metrics per cell
  .sh       cell(N_roisets, N_subs)  — split-half metrics per cell
  .qc       cell(N_roisets, 1)       — QC summaries from rdcm_qc_group
  .params   struct                   — params as used
  .data_meta struct                  — data.subs and data.ROIs (no timeseries)
%}

%% ===================================================================
% SECTION 1: Parameter Defaults and Validation
%% ===================================================================

if ~isfield(params, 'verbose'),  params.verbose  = 2;      end
if ~isfield(params, 'log_mode'), params.log_mode = 'both'; end

if ~isfield(params, 'dirs') || ~isfield(params.dirs, 'output')
    error('rdcm_pipeline: params.dirs.output is required.');
end
if ~isfield(params.dirs, 'logs')
    params.dirs.logs = fullfile(params.dirs.output, 'logs');
end

if ~isfield(params, 'rdcm'), params.rdcm = struct(); end
if ~isfield(params.rdcm, 'SNR'),        params.rdcm.SNR        = 1; end
if ~isfield(params.rdcm, 'est_method'), params.rdcm.est_method = 1; end
if ~isfield(params.rdcm, 'dt')
    error('rdcm_pipeline: params.rdcm.dt (TR in seconds) is required.');
end

if nargin < 3 || isempty(roi_subsets)
    error('rdcm_pipeline: roi_subsets is required (run rdcm_setup first).');
end
if nargin < 4 || isempty(grid)
    error('rdcm_pipeline: grid is required (run rdcm_setup first).');
end

if ~isfield(params, 'splithalf'), params.splithalf = struct(); end
if ~isfield(params.splithalf, 'exclude_diagonal'),  params.splithalf.exclude_diagonal  = false; end
if ~isfield(params.splithalf, 'min_timepoints'),    params.splithalf.min_timepoints    = 50;    end
if ~isfield(params.splithalf, 'store_half_fits'),   params.splithalf.store_half_fits   = false; end

if ~isfield(params, 'qc'), params.qc = struct(); end
if ~isfield(params.qc, 'icc_floor'), params.qc.icc_floor = 0.40; end
if ~isfield(params.qc, 'icc_s_pass'), params.qc.icc_s_pass = 0.60; end
if ~isfield(params.qc, 'icc_c_fail'), params.qc.icc_c_fail = 0.40; end
if ~isfield(params.qc, 'icc_c_pass'), params.qc.icc_c_pass = 0.60; end
if ~isfield(params.qc, 'bias_threshold'), params.qc.bias_threshold = -0.10; end
if ~isfield(params.qc, 'outlier_z'), params.qc.outlier_z = 3.5; end
if ~isfield(params.qc, 'kl_mask_method'), params.qc.kl_mask_method = 'top_k'; end
if ~isfield(params.qc, 'kl_mask_k'), params.qc.kl_mask_k = 1000; end
if ~isfield(params.qc, 'kl_mask_p'), params.qc.kl_mask_p = 0.75; end
if ~isfield(params.qc, 'kl_mask_min_connections'), params.qc.kl_mask_min_connections = 20; end
if ~isfield(params.qc, 'kl_mask_obs_param_ratio'), params.qc.kl_mask_obs_param_ratio = 10; end

if ~isfield(params, 'report'), params.report = struct(); end
if ~isfield(params.report, 'skip_plots'),     params.report.skip_plots     = false; end
if ~isfield(params.report, 'subject_roster'), params.report.subject_roster = true;  end
if ~isfield(params.report, 'deep_dive_ids'),  params.report.deep_dive_ids  = {};    end

%% ===================================================================
% SECTION 2: Directory Setup
%% ===================================================================

for d = {params.dirs.output, params.dirs.logs}
    if ~isfolder(d{1}), mkdir(d{1}); end
end

%% ===================================================================
% SECTION 3: Run ID and Log File Paths
%% ===================================================================

run_id = grid.run_id;

params.log.run_id   = run_id;
params.log.ck_file  = fullfile(params.dirs.logs, sprintf('checkpoint_%s.log', run_id));
params.log.log_file = fullfile(params.dirs.logs, sprintf('pipeline_%s.log',   run_id));

%% ===================================================================
% SECTION 4: Checkpoint Init
%% ===================================================================

completed = rdcm_checkpoint(params, 'init', '');

%% ===================================================================
% SECTION 5: Grid Summary
%% ===================================================================

N_roisets = grid.N_roisets;
N_subs    = grid.N_subs;

rdcm_log(params, 1, 'Grid: %d ROI sets x %d subjects = %d total jobs\n', ...
    N_roisets, N_subs, grid.N_jobs);
rdcm_log(params, 1, 'run_id: %s\n', run_id);

cell_EC  = cell(N_roisets, N_subs);
cell_fit = cell(N_roisets, N_subs);
cell_sh  = cell(N_roisets, N_subs);

saved_cells_dir = fullfile(params.dirs.output, 'cells');
if ~isfolder(saved_cells_dir), mkdir(saved_cells_dir); end

%% ===================================================================
% SECTION 6: Main Grid Loop
%% ===================================================================

rdcm_log(params, 1, '\n===== rDCM PIPELINE: MAIN LOOP START =====\n');

for i_rs = 1:N_roisets
    for i_sub = 1:N_subs

        %% --- Token and checkpoint query ----------------------------
        token = sprintf('rs%d_sub%02d', i_rs, i_sub);

        % Guard: if map was wiped by a recompile, rebuild from disk
        if isempty(completed) || completed.Count == 0
            completed = rdcm_checkpoint(params, 'init', '');
        end

        if isKey(completed, token)
            rdcm_log(params, 3, '  [skip] %s (already complete)\n', token);
            cell_file = fullfile(saved_cells_dir, [token '.mat']);
            if isfile(cell_file)
                tmp = load(cell_file, 'EC_out', 'fit_out', 'sh_out');
                cell_EC{i_rs,  i_sub} = tmp.EC_out;
                cell_fit{i_rs, i_sub} = tmp.fit_out;
                cell_sh{i_rs,  i_sub} = tmp.sh_out;
            end
            continue
        end

        %% --- Per-cell setup ----------------------------------------
        rs_i   = roi_subsets(i_rs);
        Y_full = data.Yall{i_sub};

        rdcm_log(params, 2, '\n[%s] ROI set: %s | Subject: %s\n', ...
            token, rs_i.label, Y_full.subj);

        %% --- Subset BOLD timeseries to current ROI set -------------
        [~, roi_idx] = ismember(rs_i.roi_names, Y_full.name);
        missing = rs_i.roi_names(roi_idx == 0);
        if ~isempty(missing)
            rdcm_log(params, 1, ...
                '  WARNING [%s]: %d ROI(s) not found in subject data — skipping.\n', ...
                token, numel(missing));
            rdcm_log(params, 2, '  Missing: %s\n', strjoin(missing, ', '));
            continue
        end

        Y_sub       = Y_full;
        Y_sub.y     = Y_full.y(:, roi_idx);
        Y_sub.name  = rs_i.roi_names;

        %% --- Full-session rDCM [RESTING-STATE] --------------------
        rdcm_log(params, 2, '  Running full-session rDCM...\n');
        t_start = tic;

        try
            [EC_out, fit_out] = rdcm_run_single(Y_sub, params);
        catch ME
            rdcm_log(params, 1, ...
                '  ERROR [%s] rdcm_run_single failed: %s\n', token, ME.message);
            continue
        end

        rdcm_log(params, 2, '  Full-session done (%.1fs). cosine=%.3f\n', ...
            toc(t_start), fit_out.cosine_sim_mean);

        %% --- Split-half reliability --------------------------------
        rdcm_log(params, 2, '  Running split-half reliability...\n');
        t_sh = tic;

        try
            sh_out = rdcm_splithalf_single(Y_sub, params);
        catch ME
            rdcm_log(params, 1, ...
                '  ERROR [%s] rdcm_splithalf_single failed: %s\n', token, ME.message);
            sh_out         = struct();
            sh_out.subj    = Y_sub.subj;
            sh_out.failed  = true;
            sh_out.error   = ME.message;
        end

        rdcm_log(params, 2, ' Split-half done (%.1fs). N_TRs_half=%d\n', ...
            toc(t_sh), sh_out.N_TRs_half);
        if ~isempty(sh_out.warning)
            rdcm_log(params, 2, ' Split-half warning: %s\n', sh_out.warning);
        end

        %% --- Store and save ----------------------------------------
        cell_EC{i_rs,  i_sub} = EC_out;
        cell_fit{i_rs, i_sub} = fit_out;
        cell_sh{i_rs,  i_sub} = sh_out;

        save(fullfile(saved_cells_dir, [token '.mat']), ...
            'EC_out', 'fit_out', 'sh_out', '-v7.3');

        completed = rdcm_checkpoint(params, 'mark', token, completed);

    end % i_sub
end % i_rs

rdcm_log(params, 1, '\n===== MAIN LOOP COMPLETE =====\n');

%% ===================================================================
% SECTION 7: Group QC (per ROI set)
%% ===================================================================

rdcm_log(params, 1, '\n--- Section 7: Group QC ---\n');

results.grid      = cell_EC;
results.fit       = cell_fit;
results.sh        = cell_sh;
results.qc        = cell(N_roisets, 1);
results.params    = params;
results.data_meta = struct('subs', data.subs, 'ROIs', data.ROIs);

for i_rs = 1:N_roisets
    rs_label = roi_subsets(i_rs).label;
    rdcm_log(params, 2, '  QC: %s\n', rs_label);

    qc_i = rdcm_qc_group(cell_EC(i_rs,:), cell_fit(i_rs,:), cell_sh(i_rs,:), params);

    % Attach subject IDs from data if available
    if isfield(data, 'subject_ids') && numel(data.subject_ids) == N_subs
        qc_i.subject_ids = data.subject_ids(:)';
    end

    results.qc{i_rs} = qc_i;

    rdcm_log(params, 1, ...
        '  QC [%s]: pass=%d  marginal=%d  fail=%d  hard_gate=%d\n', ...
        rs_label, qc_i.N_pass, qc_i.N_marginal, qc_i.N_fail, qc_i.N_hard_gate);
end

%% ===================================================================
% SECTION 8: Cross-network Advisory E
%  A subject whose tier is discordant across ROI sets (e.g. Pass in
%  one set, Fail in another) receives Adv E in all affected sets.
%  This is the only advisory that cannot be computed within a single
%  ROI set and must be set by the orchestrator after all QC structs
%  exist.
%% ===================================================================

rdcm_log(params, 1, '\n--- Section 8: Cross-network Advisory E ---\n');

if N_roisets > 1
    % Build [N_roisets x N_subs] tier matrix (NaN for hard-gated)
    tier_mat = NaN(N_roisets, N_subs);
    for i_rs = 1:N_roisets
        tv = results.qc{i_rs}.subject_tier(:)';
        tv_nan = double(tv);
        tv_nan(tv == 0) = NaN;
        tier_mat(i_rs, :) = tv_nan;
    end

    for i_sub = 1:N_subs
        col    = tier_mat(:, i_sub);
        finite = col(isfinite(col));
        if numel(finite) < 2, continue; end
        % Discordant if range across non-gated sets spans a tier boundary
        discordant = (max(finite) - min(finite)) >= 2;
        if discordant
            for i_rs = 1:N_roisets
                if ~results.qc{i_rs}.hard_gate_flag(i_sub)
                    results.qc{i_rs}.adv_crossnet(i_sub)       = true;
                    % Recompute review_recommended for this subject
                    n_adv = sum([ ...
                        results.qc{i_rs}.adv_logF(i_sub), ...
                        results.qc{i_rs}.adv_reliability(i_sub), ...
                        results.qc{i_rs}.adv_bias(i_sub), ...
                        results.qc{i_rs}.adv_ec_outlier(i_sub), ...
                        results.qc{i_rs}.adv_crossnet(i_sub)]);
                    results.qc{i_rs}.review_recommended(i_sub) = n_adv >= 2;
                end
            end
        end
    end

    n_adv_e = sum(cellfun(@(q) sum(q.adv_crossnet), results.qc));
    rdcm_log(params, 1, '  Adv E flags set: %d (subject-set instances)\n', n_adv_e);
else
    rdcm_log(params, 2, '  Single ROI set — Adv E not applicable.\n');
end

%% ===================================================================
% SECTION 9: Reporting
%% ===================================================================

rdcm_log(params, 1, '\n--- Section 9: Reporting ---\n');

% 9a — QC group report (text)
try
    rdcm_report_qc(results, params, roi_subsets);
catch ME
    rdcm_log(params, 1, '  WARNING: rdcm_report_qc failed: %s\n', ME.message);
end

% 9b — Subject roster
if params.report.subject_roster
    try
        rdcm_report_subject(results, params, roi_subsets);
    catch ME
        rdcm_log(params, 1, '  WARNING: rdcm_report_subject (roster) failed: %s\n', ME.message);
    end
end

% 9c — Per-subject deep-dive reports
if ~isempty(params.report.deep_dive_ids)
    for dd = 1:numel(params.report.deep_dive_ids)
        sid = params.report.deep_dive_ids{dd};
        try
            rdcm_report_subject(results, params, roi_subsets, sid);
        catch ME
            rdcm_log(params, 1, ...
                '  WARNING: rdcm_report_subject deep-dive (%s) failed: %s\n', ...
                sid, ME.message);
        end
    end
end

% 9d — QC distribution figures
if ~params.report.skip_plots
    try
        rdcm_plot_qc_distributions(results, params, roi_subsets);
    catch ME
        rdcm_log(params, 1, '  WARNING: rdcm_plot_qc_distributions failed: %s\n', ME.message);
    end
else
    rdcm_log(params, 2, '  Plots skipped (params.report.skip_plots = true).\n');
end

%% ===================================================================
% SECTION 10: Export Manifest, EC Matrices, and Visual QC
%% ===================================================================

rdcm_log(params, 1, '\n--- Section 10: Export Manifest ---\n');

try
    rdcm_export_manifest(results.qc, ...
    cellfun(@(rs) rs.label, num2cell(roi_subsets), 'UniformOutput', false), ...
    params, ...
    fullfile(params.dirs.output, 'export'));
catch ME
    rdcm_log(params, 1, ' WARNING: rdcm_export_manifest failed: %s\n', ME.message);
end

% 10a — EC matrix export
try
    rdcm_export_ec(results, params, roi_subsets);
catch ME
    rdcm_log(params, 1, ' WARNING: rdcm_export_ec failed: %s\n', ME.message);
end

% 10b — Visual QC PDFs (one per subject, one page per ROI set)
if ~params.report.skip_plots
    try
        rdcm_visual_qc(results, params, roi_subsets);
    catch ME
        rdcm_log(params, 1, ' WARNING: rdcm_visual_qc failed: %s\n', ME.message);
    end
else
    rdcm_log(params, 2, ' Visual QC skipped (params.report.skip_plots = true).\n');
end

%% ===================================================================
% SECTION 11: Save Final Results
%% ===================================================================

rdcm_log(params, 1, '\n--- Section 11: Save Results ---\n');

results_file = fullfile(params.dirs.output, sprintf('results_%s.mat', run_id));
save(results_file, 'results', '-v7.3');

rdcm_log(params, 1, 'Results saved to: %s\n', results_file);
rdcm_log(params, 1, 'Pipeline complete. Run ID: %s\n', run_id);

end % rdcm_pipeline
