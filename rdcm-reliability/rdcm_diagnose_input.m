function data = rdcm_diagnose_input(data, params)
%% rdcm_diagnose_input — validate, diagnose, and trim rDCM input data
%{
Always called from rdcm_setup immediately after any rdcm_convert_input_*
function. Converter functions produce the minimal data contract; this
function handles all format-agnostic work that would otherwise have to
be duplicated across converters.

Responsibilities:
  1. Apply format-agnostic defaults (verbose, log_mode, params.log)
  2. Validate the minimal data contract from the converter
  3. Diagnose each subject's timeseries for problematic ROIs
  4. Apply the configured exclusion strategy (exclude_ROI / exclude_sub)
  5. Trim data.Yall and data.all_roi_names to retained entries
  6. Append bookkeeping structs (data.subs, data.ROIs)
  7. Save the enriched data struct to params.filenames.data_mat

MINIMAL INPUT CONTRACT (produced by rdcm_convert_input_*):
  data.Yall          {1 x N_subs} cell, each entry:
    .y               [N_t x N_ROI] double
    .dt              TR in seconds
    .name            {1 x N_ROI} ROI label cell
    .subj            subject ID string
  data.all_roi_names {1 x N_ROI} cell
  data.subject_ids   {1 x N_subs} cell

OUTPUT: enriched data struct, adding:
  data.subs.subs_all      struct with .name — all subjects before diagnosis
  data.subs.subs_problem  struct with .bool, .name — subjects with problems
  data.subs.subs_input    struct with .bool, .name — subjects kept
  data.ROIs.ROIs_all      struct with .name — all ROIs before diagnosis
  data.ROIs.ROIs_problem  struct with .bool, .name — ROIs flagged
  data.ROIs.ROIs_input    struct with .bool, .name — ROIs kept
  data.Yall               trimmed to retained subjects x retained ROIs
  data.all_roi_names      trimmed to retained ROIs
  data.subject_ids        trimmed to retained subjects

A ROI is flagged as problematic in a given subject if its column is:
  - entirely non-finite (all NaN or all Inf)
  - entirely zero
  - zero-variance
A single occurrence across any subject flags the ROI (exclude_ROI)
or that subject (exclude_sub).

params.diagnose_ROIs controls strategy:
  'exclude_ROI' — remove problem ROIs, keep all subjects (default)
  'exclude_sub' — remove subjects with any problem ROI, keep all ROIs
%}

%% ===================================================================
% FORMAT-AGNOSTIC DEFAULTS (owned here, not in converters)
%% ===================================================================

if ~isfield(params, 'diagnose_ROIs') || isempty(params.diagnose_ROIs)
    params.diagnose_ROIs = 'exclude_ROI';
end
if ~isfield(params, 'verbose'),  params.verbose  = 2;      end
if ~isfield(params, 'log_mode'), params.log_mode = 'both'; end

if ~isfield(params, 'log') || ~isfield(params.log, 'log_file')
    if ~isfolder(params.dirs.logs), mkdir(params.dirs.logs); end
    params.log.log_file = fullfile(params.dirs.logs, ...
        sprintf('diagnose_input_%s.log', datestr(now, 'yyyymmdd_HHMMSS')));
    params.log.run_id = 'diagnose_input';
end

rdcm_log(params, 1, '=== rdcm_diagnose_input ===\n');
rdcm_log(params, 1, 'diagnose_ROIs strategy: %s\n', params.diagnose_ROIs);

%% ===================================================================
% VALIDATE MINIMAL CONTRACT
%% ===================================================================

required_top = {'Yall', 'all_roi_names', 'subject_ids'};
for f = 1:numel(required_top)
    if ~isfield(data, required_top{f})
        error(['rdcm_diagnose_input: data.%s is missing.\n' ...
               'rdcm_convert_input_* must produce: Yall, all_roi_names, subject_ids.'], ...
            required_top{f});
    end
end

N_subs = numel(data.Yall);
N_ROI  = numel(data.all_roi_names);

if N_subs == 0
    error('rdcm_diagnose_input: data.Yall is empty.');
end
if N_ROI == 0
    error('rdcm_diagnose_input: data.all_roi_names is empty.');
end
if numel(data.subject_ids) ~= N_subs
    error('rdcm_diagnose_input: subject_ids has %d entries but Yall has %d cells.', ...
        numel(data.subject_ids), N_subs);
end

required_y = {'y', 'dt', 'name', 'subj'};
for i = 1:N_subs
    Y_i = data.Yall{i};
    if isempty(Y_i)
        error('rdcm_diagnose_input: data.Yall{%d} is empty (subject "%s" had no data).', ...
            i, data.subject_ids{i});
    end
    for f = 1:numel(required_y)
        if ~isfield(Y_i, required_y{f})
            error('rdcm_diagnose_input: data.Yall{%d} missing field ".%s".', ...
                i, required_y{f});
        end
    end
    if ~isnumeric(Y_i.y) || ~ismatrix(Y_i.y)
        error('rdcm_diagnose_input: data.Yall{%d}.y must be a 2-D numeric matrix.', i);
    end
    if size(Y_i.y, 2) ~= N_ROI
        error(['rdcm_diagnose_input: Yall{%d}.y has %d columns but ' ...
               'all_roi_names has %d entries.'], i, size(Y_i.y, 2), N_ROI);
    end
end

rdcm_log(params, 1, 'Found %d subjects, %d ROIs.\n', N_subs, N_ROI);

%% ===================================================================
% STEP 1: DIAGNOSE ROIs ACROSS SUBJECTS
%% ===================================================================

rdcm_log(params, 1, 'Step 1: Diagnosing ROIs across subjects...\n');

ROIs_problem_bool = false(1, N_ROI);
subs_problem_bool = false(N_subs, 1);

for i = 1:N_subs
    ts_i           = data.Yall{i}.y;
    problem_mask_i = false(1, N_ROI);

    for r = 1:N_ROI
        col = ts_i(:, r);
        if all(~isfinite(col)) || all(col == 0) || var(double(col)) == 0
            problem_mask_i(r) = true;
        end
    end

    if any(problem_mask_i)
        subs_problem_bool(i)             = true;
        ROIs_problem_bool(problem_mask_i) = true;
        rdcm_log(params, 2, '  Subject %s: %d problem ROI(s): %s\n', ...
            data.subject_ids{i}, sum(problem_mask_i), ...
            strjoin(data.all_roi_names(problem_mask_i), ', '));
    end
end

subs_all.name        = data.subject_ids;
ROIs_all.name        = data.all_roi_names;
ROIs_problem.bool    = ROIs_problem_bool;
ROIs_problem.name    = data.all_roi_names(ROIs_problem_bool);
subs_problem.bool    = subs_problem_bool;
subs_problem.name    = data.subject_ids(subs_problem_bool);

rdcm_log(params, 1, 'Diagnosis complete: %d problem ROI(s), %d problem subject(s).\n', ...
    sum(ROIs_problem_bool), sum(subs_problem_bool));

%% ===================================================================
% STEP 2: APPLY EXCLUSION STRATEGY
%% ===================================================================

switch lower(params.diagnose_ROIs)
    case 'exclude_roi'
        ROIs_input.bool = ~ROIs_problem_bool;
        ROIs_input.name = data.all_roi_names(ROIs_input.bool);
        subs_input.bool = true(N_subs, 1);
        subs_input.name = data.subject_ids;
        rdcm_log(params, 1, 'Strategy: excluding %d ROI(s), keeping all %d subjects.\n', ...
            sum(ROIs_problem_bool), N_subs);

    case 'exclude_sub'
        subs_input.bool = ~subs_problem_bool;
        subs_input.name = data.subject_ids(subs_input.bool);
        ROIs_input.bool = true(1, N_ROI);
        ROIs_input.name = data.all_roi_names;
        rdcm_log(params, 1, 'Strategy: excluding %d subject(s), keeping all %d ROIs.\n', ...
            sum(subs_problem_bool), N_ROI);

    otherwise
        error(['rdcm_diagnose_input: unknown diagnose_ROIs value "%s". ' ...
               'Use ''exclude_ROI'' or ''exclude_sub''.'], params.diagnose_ROIs);
end

if sum(subs_input.bool) == 0
    error('rdcm_diagnose_input: no subjects remain after diagnosis. Check your data.');
end
if sum(ROIs_input.bool) == 0
    error('rdcm_diagnose_input: no ROIs remain after diagnosis. Check your data.');
end

%% ===================================================================
% STEP 3: TRIM Yall AND all_roi_names
%% ===================================================================

subs_keep_idx = find(subs_input.bool);
roi_keep_idx  = find(ROIs_input.bool);
N_subs_out    = numel(subs_keep_idx);

rdcm_log(params, 1, 'Step 2: Assembling Y structs (%d subjects, %d ROIs)...\n', ...
    N_subs_out, sum(ROIs_input.bool));

Yall_out        = cell(1, N_subs_out);
subject_ids_out = cell(1, N_subs_out);

for k = 1:N_subs_out
    i     = subs_keep_idx(k);
    Y_old = data.Yall{i};

    Y_new.y    = Y_old.y(:, roi_keep_idx);
    Y_new.dt   = Y_old.dt;
    Y_new.name = ROIs_input.name;
    Y_new.subj = Y_old.subj;

    Yall_out{k}        = Y_new;
    subject_ids_out{k} = Y_old.subj;

    rdcm_log(params, 2, '  Assembled %s  [%d TRs x %d ROIs]\n', ...
        Y_old.subj, size(Y_new.y, 1), size(Y_new.y, 2));
end

%% ===================================================================
% WRITE BACK AND SAVE
%% ===================================================================

data.Yall              = Yall_out;
data.all_roi_names     = ROIs_input.name;
data.subject_ids       = subject_ids_out;
data.subs.subs_all     = subs_all;
data.subs.subs_problem = subs_problem;
data.subs.subs_input   = subs_input;
data.ROIs.ROIs_all     = ROIs_all;
data.ROIs.ROIs_problem = ROIs_problem;
data.ROIs.ROIs_input   = ROIs_input;

out_path = fullfile(params.dirs.output, params.filenames.data_mat);
if ~isfolder(params.dirs.output), mkdir(params.dirs.output); end
save(out_path, 'data', '-v7.3');

rdcm_log(params, 1, 'data struct saved to:\n  %s\n', out_path);
rdcm_log(params, 1, '  %d subjects, %d ROIs, %.1fs TR\n', ...
    N_subs_out, numel(ROIs_input.name), data.Yall{1}.dt);
rdcm_log(params, 1, '=== rdcm_diagnose_input complete ===\n\n');

end
