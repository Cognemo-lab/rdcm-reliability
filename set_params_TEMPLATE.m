function params = set_params()
%% set_params — rDCM pipeline configuration
%{
INSTRUCTIONS: Copy set_params_TEMPLATE.m and save the copy as set_params.m
Add set_params.m to your .gitignore. Edit set_params.m with your real paths
and settings. Never edit set_params_TEMPLATE.m directly.

set_params.m is the single source of truth for a pipeline run. It is called
by rdcm_setup, which saves a frozen copy (rdcm_params.mat) into the output
directory. run_rdcm_pipeline loads that frozen copy, then optionally patches
cluster paths by looking for set_params_pipeline.m.
%}

%% ===================================================================
% SECTION 0: ROOT FOLDER & COMMON INPUTS
%% ===================================================================

% Root folder--must contain:
% - BOLD data folder
% - ROI list ordered as per BOLD data (TSV or CSV)
params.dirs.root = '/path/to/your/workspace';

% ROI list (TXT or CSV)
params.filenames.roi_list = 'your_roi_labels.csv';

% Network/community labels for each ROI, same order as roi_list.
% Used for network-sorted QC figures. Leave as '' if not available.
params.filenames.network_labels = 'your_network_labels.csv';

%% ===================================================================
% SECTION 1: INPUT CONVERSION
%% ===================================================================

% Name of the input conversion function to call during rdcm_setup.
% This function must be on the MATLAB path and must have the signature:
% data = my_convert_fn(params)
% where data has fields: Yall, all_roi_names, subject_ids, subs, ROIs.

% Shipped examples: 'rdcm_convert_input_bids'
%                    'rdcm_convert_input_onesubmat'
%                    'rdcm_convert_input_allsubsmat'
% Custom example:    'rdcm_convert_input_mylab'
params.convert_input_fn = 'rdcm_convert_input_*';
params.params_input_fn  = 'set_params_input_*';

%% ===================================================================
% SECTION 2: OUTPUT DIRECTORIES
%% ===================================================================

% Root output directory. All pipeline outputs go here.
params.dirs.output = fullfile(params.dirs.root, 'rdcm-output');

% Logs and checkpoint files (default: subdirectory of output).
params.dirs.logs = fullfile(params.dirs.output, 'logs');

%% ===================================================================
% SECTION 3: OUTPUT FILENAMES
% (relative to params.dirs.output — do not include the full path)
%% ===================================================================

params.filenames.data_mat    = 'rdcm_input.mat';
params.filenames.roi_subsets = 'roi_subsets.mat';
params.filenames.grid_mat    = 'rdcm_grid.mat';
params.filenames.params_mat  = 'rdcm_params.mat';

%% ===================================================================
% SECTION 4: rDCM ESTIMATION SETTINGS
%% ===================================================================

% TR in seconds — must match your BOLD data.
params.rdcm.dt = 0.0; % <-- SET THIS: your acquisition TR in seconds

% Estimation method:
% 1 = ridge regression rDCM (faster)
% 2 = sparse rDCM (slower, optimises sparsity hyperparameter)
params.rdcm.est_method = 1;

%% ===================================================================
% SECTION 5: SUBJECT FILTERING (applied at grid-build time)
%% ===================================================================

% Subjects to include in this pipeline run.
% Set to {} to include all subjects in data.Yall.
% Set to a cell array of subject ID strings to restrict the run:
% params.include_subs = {'sub-CMH001', 'sub-CMH002'};
% This is the soft-exclusion mechanism — it does not require re-running
% the input conversion step.
params.include_subs = {};

%% ===================================================================
% SECTION 6: ROI DIAGNOSIS STRATEGY
%% ===================================================================

% What to do when a subject has non-numeric/NaN values in one or more ROIs:
% 'exclude_ROI' Remove the problem ROI(s) from all subjects (default)
% 'exclude_sub' Remove the problem subject(s) from the dataset
params.diagnose_ROIs = 'exclude_ROI';

%% ===================================================================
% SECTION 7: SPLIT-HALF RELIABILITY SETTINGS
%% ===================================================================

params.splithalf.exclude_diagonal = false;
params.splithalf.min_timepoints   = 50;

%% ===================================================================
% SECTION 8: QC THRESHOLDS
%% ===================================================================

% --- Hard gate -------------------------------------------------------
% Subjects with full-session cosine_sim_mean < 0 are excluded entirely
% (tier 0). This threshold is fixed and not user-configurable.

% --- KL mask strategy ------------------------------------------
% The KL divergence mask determines which connections are included in ICC
% computations. Off-diagonal entries are always used as the pool;
% the mask then selects the strongest ones.

% 'top_k'        — include the top K strongest connections by abs(A).
%                  Set params.qc.kl_mask_k (default: 1000, per
%                  Frässle & Stephan 2022).
% 'top_quartile' — include the top 25% of connections by abs(A).

% For per-subject ICC (ICC_s): mask is built from the individual
% subject's full-session A matrix.
% For per-connection ICC (ICC_c): mask is built from the group-mean
% full-session A matrix (include_group subjects only). The mask is
% used only when computing a summary mean ICC_c; individual
% per-connection ICC values are computed for all connections regardless.
params.qc.kl_mask_method          = 'top_quartile'; % 'top_k' | 'top_quartile'
params.qc.kl_mask_k               = 1000;  % used only when method='top_k'
params.qc.kl_mask_p               = 0.75;
params.qc.kl_mask_min_connections = 20;
params.qc.kl_mask_obs_param_ratio = 1;

% --- Per-subject ICC thresholds (ICC_s) ------------------------------
% ICC_s is the Spearman-Brown-corrected ICC(3,1) computed on strong
% off-diagonal connections of the individual subject's A matrix.
% Three ICC_s values are computed per subject (see rdcm_qc_group docs),
% but tier assignment uses the strong off-diagonal ICC_s only.

% Tier assignment:
% ICC_s < icc_floor              -> Tier 1 (fail) / Gate 1b (hard-excluded)
% icc_floor <= ICC_s < icc_s_pass -> Tier 2 (marginal)
% ICC_s >= icc_s_pass            -> Tier 3 (pass)

% Inclusion vectors produced by rdcm_qc_group:
% include_all        — Tier 1 + 2 + 3
% include_group      — Tier 2 + 3 (group-level analysis)
% include_individual — Tier 3 only (individual-level analysis)

% Suggested defaults (Cicchetti 2001 bands):
% icc_floor  = 0.40 (below "Fair")
% icc_s_pass = 0.60 (lower boundary of "Good")
params.qc.icc_floor  = 0.40;
params.qc.icc_s_pass = 0.60;

% --- Per-connection ICC thresholds (ICC_c) ---------------------------
% ICC_c is the ICC(3,1) computed per connection across include_group
% subjects, with Spearman-Brown correction. Computed for all off-diagonal
% connections; strength mask applied only to the summary mean ICC_c.
% Tier assignment is stored and reported only — does not affect subject
% inclusion.

% ICC_c < icc_c_fail                  -> connection unreliable
% icc_c_fail <= ICC_c < icc_c_pass    -> connection marginal
% ICC_c >= icc_c_pass                 -> connection reliable
params.qc.icc_c_fail = 0.40;
params.qc.icc_c_pass = 0.60;

% --- Advisory flag thresholds (never cause exclusion) -----------------
params.qc.mad_k          = 2;     % Advisory B: median - k*MAD multiplier
params.qc.bias_threshold = -0.10; % Advisory C: bias index threshold
params.qc.outlier_z      = 3.5;   % Robust z-score flag on KL-masked ICC_s_all

% --- Motion gate (optional; requires FD data passed to rdcm_qc_group) --
params.qc.fd_threshold      = Inf; % mean framewise displacement (mm)
params.qc.fd_prop_threshold = Inf; % proportion of volumes scrubbed

%% ===================================================================
% SECTION 9: EXPORT SETTINGS
%% ===================================================================

params.export.include_marginal = true;

%% ===================================================================
% SECTION 10: REPORTING SETTINGS
%% ===================================================================

params.report.skip_plots      = false; % true on headless cluster nodes
params.report.subject_roster  = true;  % write per-subject roster
params.report.deep_dive_ids   = {};    % subject IDs for deep-dive reports

%% ===================================================================
% SECTION 11: LOGGING AND VERBOSITY
%% ===================================================================

params.verbose  = 2;
params.log_mode = 'both';

end