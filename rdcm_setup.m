function rdcm_setup(root)
%% rdcm_setup
%{
Phase 1 setup script — run this ONCE locally before cluster submission.

Steps:
  1. Convert input data (format-specific converter)
  2. Diagnose and trim ROIs/subjects (rdcm_diagnose_input)
  3. Define ROI subsets interactively (rdcm_define_roisets)
  4. Build job grid (rdcm_build_grid)

After completion, transfer the output folder to the cluster and submit
run_rdcm_pipeline.m.

Usage:
  rdcm_setup('/path/to/project')   — recommended: explicit root
  rdcm_setup()                     — fallback: uses current directory as root

Prerequisites:
  - set_params.m exists in root (copy from set_params_TEMPLATE.m)
  - The function named in params.convert_input_fn is on the MATLAB path
  - TAPAS rDCM and SPM12 are on the MATLAB path
%}

%% ===================================================================
% STEP 0: Resolve root and load parameters
%% ===================================================================

if nargin < 1 || isempty(root)
    root = pwd;
    root_source = 'fallback';
else
    root_source = 'explicit';
end

set_params_path = fullfile(root, 'set_params.m');
if ~isfile(set_params_path)
    error(['rdcm_setup: set_params.m not found in:\n  %s\n' ...
           'Copy set_params_TEMPLATE.m there, fill it in, then re-run.\n' ...
           'Or call rdcm_setup with the correct project root:\n' ...
           '  rdcm_setup(''/path/to/project'')'], root);
end

addpath(root);
params = set_params();
params.dirs.root = root;

if ~isfolder(params.dirs.output), mkdir(params.dirs.output); end
if ~isfolder(params.dirs.logs),   mkdir(params.dirs.logs);   end

% Set up a bootstrap log before rdcm_log is available
if ~isfield(params, 'verbose'),  params.verbose  = 2;      end
if ~isfield(params, 'log_mode'), params.log_mode = 'both'; end
if ~isfield(params, 'log') || ~isfield(params.log, 'log_file')
    params.log.log_file = fullfile(params.dirs.logs, ...
        sprintf('rdcm_setup_%s.log', datestr(now, 'yyyymmdd_HHMMSS')));
    params.log.run_id = 'setup';
end

fprintf('\n=== rDCM Setup ===\n');
if strcmp(root_source, 'fallback')
    msg = sprintf(['WARNING: no root directory was provided to rdcm_setup.\n' ...
                   'Falling back to current directory as project root:\n' ...
                   '  %s\n' ...
                   'To suppress this warning, call: rdcm_setup(''%s'')\n'], root, root);
    fprintf('%s\n', msg);
    rdcm_log(params, 1, '%s', msg);
else
    rdcm_log(params, 1, 'Project root: %s\n', root);
end

fprintf('Output directory: %s\n\n', params.dirs.output);
rdcm_log(params, 1, 'Output directory: %s\n', params.dirs.output);


%% ===================================================================
% STEP 1: Input conversion
%% ===================================================================

fprintf('--- Step 1: Input conversion (%s) ---\n', params.convert_input_fn);

if ~exist(params.convert_input_fn, 'file') && ...
   isempty(which(params.convert_input_fn))
    error(['rdcm_setup: convert_input_fn "%s" not found on the MATLAB path.\n' ...
           'Check params.convert_input_fn in set_params.m.'], params.convert_input_fn);
end

% STEP 1a: Convert input
data = feval(params.convert_input_fn, params);

% Promote resolved input-specific params back into params for later use
if isfield(data, 'params_used')
    if isfield(data.params_used, 'sub_id_pattern')
        params.sub_id_pattern = data.params_used.sub_id_pattern;
    end
    if isfield(data.params_used, 'cond_pattern')
        params.cond_pattern = data.params_used.cond_pattern;
    end
end

% STEP 1b: Diagnose and trim
data = rdcm_diagnose_input(data, params);

%% ===================================================================
% STEP 2: Define ROI subsets (interactive)
%% ===================================================================

fprintf('--- Step 2: Define ROI subsets ---\n');
roi_subsets = rdcm_define_roisets(params, data.all_roi_names);

%% ===================================================================
% STEP 3: Build job grid
%% ===================================================================

fprintf('--- Step 3: Build job grid ---\n');
grid = rdcm_build_grid(params, data, roi_subsets);

%% ===================================================================
% STEP 4: Save frozen params
%% ===================================================================

% Store data shape metadata for downstream reporting
params.data_meta.N_timepoints = size(data.Yall{1}.y, 1);
params.data_meta.N_subs       = numel(data.Yall);

params_path = fullfile(params.dirs.output, params.filenames.params_mat);
save(params_path, 'params');
rdcm_log(params, 1, 'Params saved to: %s\n', params_path);
fprintf('Params saved to: %s\n', params_path);

%% ===================================================================
% DONE
%% ===================================================================

fprintf('\n=== Setup complete ===\n');
fprintf('  Project root : %s\n', root);
fprintf('  Input data   : %s\n', params.filenames.data_mat);
fprintf('  ROI subsets  : %s\n', params.filenames.roi_subsets);
fprintf('  Grid         : %s\n', params.filenames.grid_mat);
fprintf('  Params       : %s\n', params.filenames.params_mat);
fprintf('  Total jobs   : %d\n', grid.N_jobs);
fprintf('  run_id       : %s\n', grid.run_id);
fprintf('\nTransfer output directory to cluster, then submit run_rdcm_pipeline.\n\n');

end