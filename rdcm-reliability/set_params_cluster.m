function params = set_params_cluster(params)
%% set_params_pipeline — cluster/remote path overrides  [NOT FOR GIT REPO]
%{
INSTRUCTIONS: Copy set_params_pipeline_TEMPLATE.m and save the copy as
set_params_pipeline.m. Add set_params_pipeline.m to your .gitignore.

This function receives the params struct saved by rdcm_setup and overrides
only the fields that differ on the cluster (typically paths and log_mode).
All rDCM estimation settings, QC thresholds, etc. are inherited unchanged.

run_rdcm_pipeline looks for this file at startup. If not found, it uses
the frozen params from rdcm_params.mat as-is (suitable for local runs).
%}

% Override output path for cluster filesystem
params.dirs.output = '/path/on/cluster/rdcm-output';
params.dirs.logs   = fullfile(params.dirs.output, 'logs');

% Cluster jobs have no console — write to file only
params.log_mode = 'file';
params.verbose  = 1;

% Optional: restrict to a subject subset for this run only
% params.include_subs = {'sub-CMH001', 'sub-CMH002'};

end