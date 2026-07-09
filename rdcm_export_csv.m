function rdcm_export_csv(qc_results, subject_ids, roi_set_names, out_dir, params)
%% Preamble
%{
Exports per-ROI-set subject QC reference tables as CSV files.

One CSV is written per ROI set. Each file has one row per subject and
columns ordered to follow the decision flow: identity → acquisition
context → hard gates → tier → inclusion vectors → ICC profile → fit
metrics → posterior informativeness → split-half consistency →
advisory flags.

These CSVs are the primary machine-readable output of the QC pipeline.
Downstream ML scripts should join on subject_id and filter on
include_group / include_individual — no text parsing required.
---------------------------------------------------------------------------
INPUTS
---------------------------------------------------------------------------
qc_results   := {1 x R} cell array of qc structs, one per ROI set,
                as returned by rdcm_qc_group.

subject_ids  := {N_subs x 1} cell array of subject ID strings.
                Must be indexed identically to the storage arrays passed
                to rdcm_qc_group.

roi_set_names := {1 x R} cell array of ROI set name strings,
                 e.g. {'whole_brain', 'Occ-PFr-SbC'}.
                 Used to name output files and populate the roi_set column.

out_dir      := char — output directory. Created if it does not exist.

params       := pipeline params struct. Relevant subfields:
  .n_timepoints   scalar  number of TRs in the full session
  .TR             scalar  repetition time in seconds
  (These are used to populate the acquisition context columns.
   If absent, those columns are filled with NaN.)
---------------------------------------------------------------------------
OUTPUTS
---------------------------------------------------------------------------
One CSV file per ROI set, written to out_dir:
  rdcm_qc_<roi_set_name>.csv

Column order (decision-flow):

  IDENTITY
    subject_id          char
    roi_set             char

  ACQUISITION CONTEXT
    n_timepoints        scalar
    TR_s                scalar  (seconds)
    n_rois              scalar
    n_offdiag           scalar  N*(N-1)
    obs_to_param_ratio  scalar  (N_TRs/2) / n_offdiag

  STAGE 1 — HARD GATES
    cosine_sim_mean     scalar
    gate_cosine         0/1    Gate 1a triggered
    icc_s_strong        scalar
    gate_icc_floor      0/1    Gate 1b triggered
    fd_mean             scalar (NaN if no FD data)
    fd_prop_scrubbed    scalar (NaN if no FD data)
    gate_motion         0/1    Gate 1c triggered
    hard_gate_flag      0/1    any gate triggered
    hard_gate_type      char   'none'|'cosine'|'icc_floor'|'motion'|...
    exclusion_reason    char   '' if not excluded

  STAGE 2 — TIER
    tier                0/1/2/3
    include_all         0/1
    include_group       0/1
    include_individual  0/1

  ICC PROFILE
    icc_s_offdiag       scalar
    icc_s_diag          scalar
    pearson_r_offdiag   scalar
    bias_index          scalar  icc_s_offdiag - pearson_r_offdiag

  FIT METRICS
    logF                scalar
    logF_advisory       0/1    Advisory A

  POSTERIOR INFORMATIVENESS
    R_mean_offdiag      scalar  (group-mean, reported at subject level as
    R_mean_strong       scalar   the group value — same for all subjects)
    sqrtD_mean_offdiag  scalar
    sqrtD_mean_strong   scalar

  ADVISORY FLAGS
    adv_logF            0/1    Advisory A
    adv_reliability     0/1    Advisory B
    adv_bias            0/1    Advisory C
    adv_ec_outlier      0/1    Advisory D
    adv_crossnet        0/1    Advisory E (set externally; 0 if not yet set)
    n_advisories        integer  count of A-D (E not included unless set)
    review_recommended  0/1
%}

%% --- Validate inputs -------------------------------------------------

if nargin < 5, params = struct(); end

R      = numel(qc_results);
N_subs = numel(subject_ids);

if numel(roi_set_names) ~= R
    error('rdcm_export_csv: roi_set_names must have the same length as qc_results.');
end

if ~exist(out_dir, 'dir')
    mkdir(out_dir);
end

% Optional acquisition context from params
n_tp = get_param(params, 'n_timepoints', NaN);
TR_s = get_param(params, 'TR',           NaN);

has_fd = isfield(params, 'fd_mean') && ~isempty(params.fd_mean);
if has_fd
    fd_mean_arr = params.fd_mean(:);
    fd_prop_arr = get_param(params, 'fd_prop', nan(N_subs, 1));
    fd_prop_arr = fd_prop_arr(:);
else
    fd_mean_arr = nan(N_subs, 1);
    fd_prop_arr = nan(N_subs, 1);
end

%% --- Write one CSV per ROI set --------------------------------------

for r = 1:R
    qc   = qc_results{r};
    rname = roi_set_names{r};

    % Sanitise ROI set name for use in filename
    fname_safe = regexprep(rname, '[^A-Za-z0-9_\-]', '_');
    out_path   = fullfile(out_dir, sprintf('rdcm_qc_%s.csv', fname_safe));

    N      = size(qc.icc_c, 1);           % N_ROIs
    n_off  = N * (N - 1);
    n_half = n_tp / 2;
    obs_to_param = n_half / n_off;         % NaN if n_tp is NaN

    % --- Derive per-subject gate breakdown from hard_gate_type ----------
    gate_cosine_col    = zeros(N_subs, 1);
    gate_icc_col       = zeros(N_subs, 1);
    gate_motion_col    = zeros(N_subs, 1);
    for i = 1:N_subs
        gt = qc.hard_gate_type{i};
        gate_cosine_col(i)  = contains(gt, 'cosine');
        gate_icc_col(i)     = contains(gt, 'icc_floor');
        gate_motion_col(i)  = contains(gt, 'motion');
    end

    % --- Build table row by row -----------------------------------------
    col_subject_id         = subject_ids(:);
    col_roi_set            = repmat({rname}, N_subs, 1);
    col_n_timepoints       = repmat(n_tp,    N_subs, 1);
    col_TR_s               = repmat(TR_s,    N_subs, 1);
    col_n_rois             = repmat(N,       N_subs, 1);
    col_n_offdiag          = repmat(n_off,   N_subs, 1);
    col_obs_to_param       = repmat(obs_to_param, N_subs, 1);

    col_cosine_sim_mean    = qc.cosine_sim;
    col_gate_cosine        = gate_cosine_col;
    col_icc_s_strong       = qc.icc_s_strong;
    col_gate_icc_floor     = gate_icc_col;
    col_fd_mean            = fd_mean_arr;
    col_fd_prop_scrubbed   = fd_prop_arr;
    col_gate_motion        = gate_motion_col;
    col_hard_gate_flag     = double(qc.hard_gate_flag);
    col_hard_gate_type     = qc.hard_gate_type;
    col_exclusion_reason   = qc.exclusion_reason;

    col_tier               = qc.subject_tier;
    col_include_all        = double(qc.include_all);
    col_include_group      = double(qc.include_group);
    col_include_individual = double(qc.include_individual);

    col_icc_s_offdiag      = qc.icc_s_offdiag;
    col_icc_s_diag         = qc.icc_s_diag;
    col_pearson_r_offdiag  = qc.pearson_r;
    col_bias_index         = qc.bias_index;

    col_logF               = qc.logF;

    % Posterior informativeness — group-level scalars broadcast to all rows
    col_R_mean_offdiag     = repmat(qc.R_mean_all,        N_subs, 1);
    col_R_mean_strong      = repmat(qc.R_mean_strong,     N_subs, 1);
    col_sqrtD_mean_offdiag = repmat(qc.sqrtD_mean_all,    N_subs, 1);
    col_sqrtD_mean_strong  = repmat(qc.sqrtD_mean_strong, N_subs, 1);

    col_adv_logF           = double(qc.adv_logF);
    col_adv_reliability    = double(qc.adv_reliability);
    col_adv_bias           = double(qc.adv_bias);
    col_adv_ec_outlier     = double(qc.adv_ec_outlier);
    col_adv_crossnet       = double(qc.adv_crossnet);
    col_n_advisories       = qc.n_advisories;
    col_review_recommended = double(qc.review_recommended);

    % --- Write CSV ------------------------------------------------------
    fid = fopen(out_path, 'w');
    if fid < 0
        error('rdcm_export_csv: cannot open file for writing: %s', out_path);
    end

    % Header
    header = [ ...
        'subject_id,roi_set,' ...
        'n_timepoints,TR_s,n_rois,n_offdiag,obs_to_param_ratio,' ...
        'cosine_sim_mean,gate_cosine,icc_s_strong,gate_icc_floor,' ...
        'fd_mean,fd_prop_scrubbed,gate_motion,' ...
        'hard_gate_flag,hard_gate_type,exclusion_reason,' ...
        'tier,include_all,include_group,include_individual,' ...
        'icc_s_offdiag,icc_s_diag,pearson_r_offdiag,bias_index,' ...
        'logF,adv_logF,' ...
        'R_mean_offdiag,R_mean_strong,sqrtD_mean_offdiag,sqrtD_mean_strong,' ...
        'adv_logF_flag,adv_reliability,adv_bias,adv_ec_outlier,adv_crossnet,' ...
        'n_advisories,review_recommended' ...
        ];
    fprintf(fid, '%s\n', header);

    % Rows
    for i = 1:N_subs
        fprintf(fid, '%s,%s,', ...
            csv_str(col_subject_id{i}), csv_str(col_roi_set{i}));

        fprintf(fid, '%s,%s,%d,%d,%s,', ...
            fmt_scalar(col_n_timepoints(i)), ...
            fmt_scalar(col_TR_s(i)), ...
            col_n_rois(i), ...
            col_n_offdiag(i), ...
            fmt_scalar(col_obs_to_param(i)));

        fprintf(fid, '%s,%d,%s,%d,', ...
            fmt_scalar(col_cosine_sim_mean(i)), ...
            col_gate_cosine(i), ...
            fmt_scalar(col_icc_s_strong(i)), ...
            col_gate_icc_floor(i));

        fprintf(fid, '%s,%s,%d,', ...
            fmt_scalar(col_fd_mean(i)), ...
            fmt_scalar(col_fd_prop_scrubbed(i)), ...
            col_gate_motion(i));

        fprintf(fid, '%d,%s,%s,', ...
            col_hard_gate_flag(i), ...
            csv_str(col_hard_gate_type{i}), ...
            csv_str(col_exclusion_reason{i}));

        fprintf(fid, '%d,%d,%d,%d,', ...
            col_tier(i), ...
            col_include_all(i), ...
            col_include_group(i), ...
            col_include_individual(i));

        fprintf(fid, '%s,%s,%s,%s,', ...
            fmt_scalar(col_icc_s_offdiag(i)), ...
            fmt_scalar(col_icc_s_diag(i)), ...
            fmt_scalar(col_pearson_r_offdiag(i)), ...
            fmt_scalar(col_bias_index(i)));

        fprintf(fid, '%s,%d,', ...
            fmt_scalar(col_logF(i)), ...
            col_adv_logF(i));

        fprintf(fid, '%s,%s,%s,%s,', ...
            fmt_scalar(col_R_mean_offdiag(i)), ...
            fmt_scalar(col_R_mean_strong(i)), ...
            fmt_scalar(col_sqrtD_mean_offdiag(i)), ...
            fmt_scalar(col_sqrtD_mean_strong(i)));

        fprintf(fid, '%d,%d,%d,%d,%d,%d,%d\n', ...
            col_adv_logF(i), ...
            col_adv_reliability(i), ...
            col_adv_bias(i), ...
            col_adv_ec_outlier(i), ...
            col_adv_crossnet(i), ...
            col_n_advisories(i), ...
            col_review_recommended(i));
    end

    fclose(fid);
    rdcm_log(params, 2, 'rdcm_export_csv: wrote %s\n', out_path);
end

end % main function

%% ====================================================================
%  Local helpers
%% ====================================================================

function v = get_param(s, field, default)
if isfield(s, field) && ~isempty(s.(field))
    v = s.(field);
else
    v = default;
end
end

function s = fmt_scalar(x)
% Format a scalar as a string: 6 significant figures, 'NaN' for missing.
if isnan(x)
    s = 'NaN';
else
    s = sprintf('%.6g', x);
end
end

function s = csv_str(x)
% Wrap a string in double quotes and escape internal quotes.
% Handles strings that may contain commas or newlines.
x = strrep(x, '"', '""');
s = ['"' x '"'];
end

