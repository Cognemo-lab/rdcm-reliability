function grid = rdcm_build_grid(params, data, roi_subsets)
%% rdcm_build_grid
%{
Pre-computes the full (roi_set x subject) job list for rdcm_pipeline and
saves it to params.filenames.grid_mat. The pipeline loads this at startup
to initialise its checkpoint system, enabling accurate progress reporting
and deterministic resumption after wall-time expiry.

The run_id is a short hash derived from the grid-defining parameters
(roi_subset labels, est_method, subject list). Changing any of these
automatically produces a new run_id and therefore a new checkpoint file,
preventing stale checkpoints from silently skipping new jobs.

NOTE: SNR is not a grid axis. For resting-state rDCM, SNR is only used
when generating synthetic data (type='s') and has no effect on empirical
estimation. See tapas_rdcm_set_options documentation.

Output .mat contains:
  grid.job_list   struct array — one element per (rs, sub) tuple
    .job_idx      scalar — linear index
    .rs_idx       scalar — ROI set index
    .rs_label     char   — ROI set label
    .sub_idx      scalar — subject index in data.Yall
    .sub_id       char   — subject ID string
  grid.N_jobs     scalar — total number of jobs
  grid.N_roisets  scalar
  grid.N_subs     scalar
  grid.roi_labels {1 x N_roisets} cell of ROI set label strings
  grid.run_id     char — hash string used for checkpoint/log files
  grid.created_at char — creation timestamp

INPUTS
  params      := struct from set_params / rdcm_setup
  data        := struct from rdcm_build_input (or loaded from data_mat)
  roi_subsets := struct array from rdcm_define_roisets
%}

%% --- Defaults -------------------------------------------------------

if ~isfield(params, 'verbose'),  params.verbose  = 2;      end
if ~isfield(params, 'log_mode'), params.log_mode = 'both'; end
if ~isfield(params.log, 'log_file')
    if ~isfolder(params.dirs.logs), mkdir(params.dirs.logs); end
    params.log.log_file = fullfile(params.dirs.logs, ...
        sprintf('build_grid_%s.log', datestr(now, 'yyyymmdd_HHMMSS')));
    params.log.run_id = 'build_grid';
end

rdcm_log(params, 1, '=== rdcm_build_grid ===\n');

if ~isfield(params.filenames, 'grid_mat') || isempty(params.filenames.grid_mat)
    params.filenames.grid_mat = 'rdcm_grid.mat';
end

%% --- Extract dimensions --------------------------------------------

N_roisets = numel(roi_subsets);
N_subs    = numel(data.Yall);
roi_labels = {roi_subsets.label};

rdcm_log(params, 1, 'Grid dimensions: %d ROI sets x %d subjects = %d jobs\n', ...
    N_roisets, N_subs, N_roisets * N_subs);

%% --- Validate ROI subset names against data.all_roi_names ----------

for i_rs = 1:N_roisets
    rs = roi_subsets(i_rs);
    missing = setdiff(rs.roi_names, data.all_roi_names);
    if ~isempty(missing)
        warning('rdcm_build_grid:missingROIs', ...
            ['ROI set "%s" contains %d ROI(s) not present in data.all_roi_names ' ...
             '(will be flagged at runtime):\n  %s'], ...
            rs.label, numel(missing), strjoin(missing, ', '));
    end
end

%% --- Build run_id from grid-defining parameters --------------------
% FUTURE (task_mode): include params.task_mode and any task-specific
% hyperparameter axes here so that switching modes invalidates the
% checkpoint automatically.

% Include sorted subject IDs so different conditions/datasets
% produce different run_ids even with identical ROI sets.
data_str = strjoin(sort(data.subject_ids), ',');
id_str = sprintf('%s|%d|%s', ...
    strjoin(sort(roi_labels), '+'), ...
    params.rdcm.est_method, ...
    data_str);
run_id = ['r' dec2hex(sum(double(id_str) .* (1:numel(id_str))), 6)];

rdcm_log(params, 1, 'run_id: %s\n', run_id);

%% --- Enumerate jobs ------------------------------------------------

N_jobs = N_roisets * N_subs;
job_list = struct( ...
    'job_idx',  cell(1, N_jobs), ...
    'rs_idx',   cell(1, N_jobs), ...
    'rs_label', cell(1, N_jobs), ...
    'sub_idx',  cell(1, N_jobs), ...
    'sub_id',   cell(1, N_jobs));

k = 0;
for i_rs = 1:N_roisets
    for i_sub = 1:N_subs
        k = k + 1;
        job_list(k).job_idx  = k;
        job_list(k).rs_idx   = i_rs;
        job_list(k).rs_label = roi_subsets(i_rs).label;
        job_list(k).sub_idx  = i_sub;
        job_list(k).sub_id   = data.subject_ids{i_sub};
    end
end

%% --- Assemble grid struct ------------------------------------------

grid.job_list   = job_list;
grid.N_jobs     = N_jobs;
grid.N_roisets  = N_roisets;
grid.N_subs     = N_subs;
grid.roi_labels = roi_labels;
grid.run_id     = run_id;
grid.created_at = datestr(now, 'yyyy-mm-dd HH:MM:SS');

%% --- Save ----------------------------------------------------------

out_path = fullfile(params.dirs.output, params.filenames.grid_mat);
if ~isfolder(params.dirs.output), mkdir(params.dirs.output); end
save(out_path, 'grid');

rdcm_log(params, 1, 'Grid saved to:\n  %s\n', out_path);
rdcm_log(params, 1, '%d total jobs (run_id: %s)\n', N_jobs, run_id);
rdcm_log(params, 1, '=== rdcm_build_grid complete ===\n');

end