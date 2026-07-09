% run_pipeline.m
% ==========================================================================
% TOP-LEVEL rDCM RELIABILITY PIPELINE RUNNER
%
% Usage:
%   1. Copy params_template.m and fill in your paths and settings.
%   2. Edit the params_file variable below to point to your config.
%   3. Run this script. Re-run naively at any time to resume — completed
%      stages are logged and skipped automatically.
%
% Output directory layout:
%   <output.root>/
%     pipeline_log.txt                    ← checkpoint log
%     <roi_set>/<SNR_tag>/
%       <sub_id>_rdcm.mat                 ← EC parameters + raw TAPAS signals
%       <sub_id>_modelfit.mat             ← cosine similarity + free energy
%       <sub_id>_splithalf.mat            ← split-half EC reliability
%
% Checkpointing:
%   Each completed (subject × stage × ROI_set × SNR) job is logged as a
%   single line in pipeline_log.txt. On re-run, the log is loaded at
%   startup and any logged job is skipped. The log is append-only, so it
%   is safe to interrupt at any point without corrupting completed entries.
% ==========================================================================

clear; clc;
rng(0, 'twister');

%% -------------------------------------------------------------------------
%  0. Load parameters
% -------------------------------------------------------------------------

params_file = 'params_template.m';   % <-- EDIT: point to your config file
run(params_file);
fprintf('Loaded params from: %s\n', params_file);

%% -------------------------------------------------------------------------
%  1. Add paths
% -------------------------------------------------------------------------

% NEEDS TO BE CHANGED TO ONLY INCLUDE THIS REPO AND TAPAS

addpath(genpath(fullfile(params.dirs.tools, 'cognemo')));
addpath(genpath(fullfile(params.dirs.tools, 'TAPAS-rDCM', 'code')));
addpath(genpath(fullfile(params.dirs.tools, 'TAPAS-rDCM', 'misc')));
addpath(genpath(fullfile(fileparts(mfilename('fullpath')), 'pipeline')));
addpath(genpath(fullfile(fileparts(mfilename('fullpath')), 'utils')));

%% -------------------------------------------------------------------------
%  2. Load subject list
% -------------------------------------------------------------------------

cd(params.dirs.root);
warning('off', 'MATLAB:table:ModifiedAndSavedVarnames');
subs_table = readtable(params.files.sub_table);
warning('on',  'MATLAB:table:ModifiedAndSavedVarnames');

usable_mask = logical(subs_table.Usable);
subs_usable = subs_table.SubjectID(usable_mask);
N_subs      = numel(subs_usable);
fprintf('Subjects: %d usable of %d total.\n', N_subs, height(subs_table));

%% -------------------------------------------------------------------------
%  3. Prepare output root and load checkpoint log
% -------------------------------------------------------------------------
out_root = fullfile(params.dirs.root, params.dirs.output.root);
if ~isfolder(out_root)
    mkdir(out_root);
end

completed = log_load(params.log_file);
fprintf('Checkpoint log: %d completed jobs found.\n\n', completed.Count);

%% -------------------------------------------------------------------------
%  4. Report pipeline scope
% -------------------------------------------------------------------------
N_roi_sets = numel(params.roi_sets);
N_SNR      = numel(params.SNR_values);
N_stages   = params.do_build_input + ...
             (params.do_rdcm + params.do_model_fit + params.do_split_half) * N_roi_sets * N_SNR;
fprintf('Scope: %d subjects  |  %d ROI sets  |  %d SNR values\n', N_subs, N_roi_sets, N_SNR);
fprintf('Max jobs per subject: %d  |  Total: %d\n\n', N_stages, N_subs * N_stages);

%% -------------------------------------------------------------------------
%  5. Main loop
% -------------------------------------------------------------------------
for s = 1:N_subs

    sub_id = subs_usable{s};
    fprintf('[%d/%d] %s\n', s, N_subs, sub_id);

    % ---- Stage 1: Build input (once per subject) -------------------------
    if params.do_build_input
        key = make_job_key(sub_id, 'build_input');
        if ~log_check(completed, key)
            fprintf('  [1] Building input...\n');
            p1_build_input(sub_id, params);
            log_write(params.log_file, key);
            completed(key) = true;
            fprintf('  [1] Done.\n');
        else
            fprintf('  [1] build_input: skipped (logged).\n');
        end
    end

    % ---- Stages 2–4: Per ROI-set × SNR -----------------------------------
    for r = 1:N_roi_sets
        roi_set = params.roi_sets(r);

        for n = 1:N_SNR
            snr_val = params.SNR_values(n);
            snr_tag = sprintf('SNR%g', snr_val);

            fprintf('  ROI set: %-20s | %s\n', roi_set.name, snr_tag);

            % Per-condition output directory
            cond_dir = fullfile(out_root, roi_set.name, snr_tag);
            if ~isfolder(cond_dir), mkdir(cond_dir); end

            % ---- Stage 2: rDCM ------------------------------------------
            if params.do_rdcm
                key = make_job_key(sub_id, 'rdcm', roi_set.name, snr_tag);
                if ~log_check(completed, key)
                    fprintf('    [2] Running rDCM...\n');
                    p2_run_rdcm(sub_id, roi_set, snr_val, params, cond_dir);
                    log_write(params.log_file, key);
                    completed(key) = true;
                    fprintf('    [2] Done.\n');
                else
                    fprintf('    [2] rdcm: skipped (logged).\n');
                end
            end

            % ---- Stage 3: Model fit (cosine sim + free energy) ----------
            if params.do_model_fit
                key = make_job_key(sub_id, 'model_fit', roi_set.name, snr_tag);
                if ~log_check(completed, key)
                    fprintf('    [3] Computing model fit...\n');
                    p3_model_fit(sub_id, params, cond_dir);
                    log_write(params.log_file, key);
                    completed(key) = true;
                    fprintf('    [3] Done.\n');
                else
                    fprintf('    [3] model_fit: skipped (logged).\n');
                end
            end

            % ---- Stage 4: Split-half reliability ------------------------
            if params.do_split_half
                key = make_job_key(sub_id, 'split_half', roi_set.name, snr_tag);
                if ~log_check(completed, key)
                    fprintf('    [4] Computing split-half reliability...\n');
                    p4_split_half(sub_id, roi_set, snr_val, params, cond_dir);
                    log_write(params.log_file, key);
                    completed(key) = true;
                    fprintf('    [4] Done.\n');
                else
                    fprintf('    [4] split_half: skipped (logged).\n');
                end
            end

        end % SNR loop
    end % ROI set loop

    fprintf('\n');
end % subject loop

fprintf('=== Pipeline complete. Total logged jobs: %d ===\n', completed.Count);
