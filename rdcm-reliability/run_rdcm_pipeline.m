function results = run_rdcm_pipeline(output_dir)
%% run_rdcm_pipeline
%{
Phase 2 pipeline runner — submit to the cluster or run locally after
rdcm_setup has completed successfully.

This script is intentionally non-interactive. It:
  1. Loads frozen params from rdcm_params.mat
  2. Applies cluster path overrides from set_params_pipeline.m (if present)
  3. Loads pre-built input data, roi_subsets, and job grid
  4. Runs rdcm_pipeline (rDCM estimation, split-half reliability, QC)
  5. Exports EC matrices and writes QC report

Usage:
  run_rdcm_pipeline('/path/to/rdcm-output')   — recommended: explicit path
  run_rdcm_pipeline()                          — fallback: uses current directory

On cluster (SLURM example):
  matlab -nodisplay -nosplash \
    -r "run_rdcm_pipeline('/scratch/peterb/rdcm-output'); exit"

On wall-time expiry: resubmit unchanged — checkpoint file ensures
completed jobs are skipped automatically.
%}

%% ===================================================================
% STEP 0: Resolve output_dir and load frozen params
%% ===================================================================

if nargin < 1 || isempty(output_dir)
    output_dir = pwd;
    dir_source = 'fallback';
else
    dir_source = 'explicit';
end

params_path = fullfile(output_dir, 'rdcm_params.mat');
if ~isfile(params_path)
    error(['run_rdcm_pipeline: rdcm_params.mat not found in:\n  %s\n' ...
           'Run rdcm_setup first and transfer the output folder here.\n' ...
           'Or call with the correct output path:\n' ...
           '  run_rdcm_pipeline(''/path/to/rdcm-output'')'], output_dir);
end

loaded = load(params_path, 'params');
params = loaded.params;
params.dirs.output = output_dir;
params.dirs.logs   = fullfile(output_dir, 'logs');
if ~isfolder(params.dirs.logs), mkdir(params.dirs.logs); end

% Bootstrap log for early messages
if ~isfield(params, 'log') || ~isfield(params.log, 'log_file')
    params.log.log_file = fullfile(params.dirs.logs, ...
        sprintf('run_pipeline_%s.log', datestr(now, 'yyyymmdd_HHMMSS')));
    params.log.run_id = 'pipeline';
end

fprintf('\n=== run_rdcm_pipeline ===\n');

% Log the directory source
if strcmp(dir_source, 'fallback')
    msg = sprintf(['WARNING: no output directory was provided to run_rdcm_pipeline.\n' ...
                   'Falling back to current directory as output path:\n' ...
                   '  %s\n' ...
                   'To suppress this warning, call: ' ...
                   'run_rdcm_pipeline(''%s'')\n'], output_dir, output_dir);
    fprintf('%s\n', msg);
    rdcm_log(params, 1, '%s', msg);
else
    rdcm_log(params, 1, 'Output directory: %s\n', output_dir);
    fprintf('Output directory: %s\n', output_dir);
end

%% Apply cluster overrides if set_params_pipeline.m is present in output_dir

set_params_pipeline_path = fullfile(output_dir, 'set_params_pipeline.m');
if isfile(set_params_pipeline_path)
    addpath(output_dir);
    params = set_params_pipeline(params);
    msg = sprintf('Cluster overrides applied from:\n  %s\n', set_params_pipeline_path);
    fprintf('%s', msg);
    rdcm_log(params, 1, '%s', msg);
else
    rdcm_log(params, 1, 'No set_params_pipeline.m found — using frozen params as-is.\n');
end

%% ===================================================================
% STEP 1: Load pre-built setup files
%% ===================================================================

required = { ...
    params.filenames.data_mat,    'Input data (run rdcm_setup first)'; ...
    params.filenames.roi_subsets, 'ROI subsets (run rdcm_setup first)'; ...
    params.filenames.grid_mat,    'Job grid (run rdcm_setup first)'};

for f = 1:size(required, 1)
    fpath = fullfile(output_dir, required{f,1});
    if ~isfile(fpath)
        error('run_rdcm_pipeline: missing file — %s\n  Expected at: %s', ...
            required{f,2}, fpath);
    end
end

loaded_data = load(fullfile(output_dir, params.filenames.data_mat),    'data');
loaded_rs   = load(fullfile(output_dir, params.filenames.roi_subsets), 'roi_subsets');
loaded_grid = load(fullfile(output_dir, params.filenames.grid_mat),    'grid');

data        = loaded_data.data;
roi_subsets = loaded_rs.roi_subsets;
grid        = loaded_grid.grid;

rdcm_log(params, 1, 'Setup files loaded.\n');
fprintf('  Subjects   : %d\n',   numel(data.Yall));
fprintf('  ROI sets   : %d\n',   numel(roi_subsets));
fprintf('  Total jobs : %d  (run_id: %s)\n\n', grid.N_jobs, grid.run_id);

%% ===================================================================
% STEP 2: Run pipeline
%% ===================================================================

results = rdcm_pipeline(params, data, roi_subsets, grid);

fprintf('\n=== run_rdcm_pipeline complete ===\n');
fprintf('Outputs written to: %s\n\n', fullfile(params.dirs.output, 'export'));

end