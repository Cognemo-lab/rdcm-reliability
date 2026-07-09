function roi_subsets = rdcm_define_roisets(params, valid_roi_names)
%% rdcm_define_roisets — interactive ROI subset definition
%{
Pure user-interaction function. Reads ROI labels from a flat CSV or TXT
file (params.filenames.roi_list), presents them interactively, and returns
a struct array of named subsets. Does not read any BOLD data files.

GUI mode (MATLAB desktop): uses listdlg / inputdlg.
Headless mode (cluster):   falls back to CLI prompts.

Required params fields:
  params.dirs.root              — study root directory
  params.filenames.roi_list     — filename of ROI label list (.csv or .txt)
  params.dirs.output            — where roi_subsets.mat is saved
  params.filenames.roi_subsets  — output filename (default: 'roi_subsets.mat')
%}

if nargin < 1, params = struct(); end
if nargin < 2, valid_roi_names = []; end
if ~isfield(params, 'filenames'), params.filenames = struct(); end
if ~isfield(params.filenames, 'roi_subsets') || isempty(params.filenames.roi_subsets)
    params.filenames.roi_subsets = 'roi_subsets.mat';
end
if ~isfield(params, 'dirs') || ~isfield(params.dirs, 'output')
    params.dirs.output = pwd;
end
if ~isfield(params.filenames, 'roi_list') || isempty(params.filenames.roi_list)
    error('rdcm_define_roisets: params.filenames.roi_list is not set in set_params.m.');
end

%% --- Detect GUI availability ----------------------------------------

use_gui = usejava('desktop') && feature('ShowFigureWindows');

if use_gui
    fprintf('GUI mode detected — using listdlg / inputdlg.\n');
else
    fprintf('Headless mode detected — using CLI prompts.\n');
end

%% --- Load ROI label list -------------------------------------------

roi_list_path = fullfile(params.dirs.root, params.filenames.roi_list);
if ~isfile(roi_list_path)
    error(['rdcm_define_roisets: ROI list file not found:\n  %s\n' ...
           'Set params.filenames.roi_list in set_params.m.'], roi_list_path);
end

[~, ~, ext] = fileparts(roi_list_path);
switch lower(ext)
    case '.csv'
        T = readtable(roi_list_path, 'ReadVariableNames', false, 'Delimiter',',');
        all_roi_names = strtrim(T{:,1})';   % first column, row-vector of strings
    case '.txt'
        fid = fopen(roi_list_path, 'r');
        raw = textscan(fid, '%s', 'Delimiter', '\n', 'Whitespace', '');
        fclose(fid);
        all_roi_names = strtrim(raw{1})';
    otherwise
        error('rdcm_define_roisets: roi_list must be .csv or .txt (got "%s").', ext);
end

all_roi_names = all_roi_names(~cellfun(@isempty, all_roi_names));  % drop blank lines
N_ROIs = numel(all_roi_names);

if N_ROIs == 0
    error('rdcm_define_roisets: no ROI names read from:\n  %s', roi_list_path);
end

fprintf('Loaded %d ROI labels from:\n  %s\n\n', N_ROIs, roi_list_path);

%% --- Intersect with valid ROIs from data diagnosis -----------------

all_roi_names_full = all_roi_names; % preserve full atlas list for label lookup

if isempty(valid_roi_names)
    excluded_roi_names = {};
else
    excluded_roi_names = setdiff(all_roi_names, valid_roi_names(:)', 'stable');
    all_roi_names      = intersect(all_roi_names, valid_roi_names(:)', 'stable');
    N_ROIs             = numel(all_roi_names);
end

if ~isempty(excluded_roi_names)
    N_excl = numel(excluded_roi_names);

    % Always print to console — gives a permanent record in session output
    fprintf('  NOTE: %d ROI(s) excluded by data diagnosis:\n', N_excl);
    for k = 1:N_excl
        fprintf('    - %s\n', excluded_roi_names{k});
    end
    fprintf('\n');

    % GUI: helpdlg before selector so user actively acknowledges
    if use_gui
        excl_msg = sprintf( ...
            ['%d ROI(s) were removed by data diagnosis and will not\n' ...
             'appear in the ROI selector below.\n\n' ...
             'Excluded ROIs:\n  %s\n\n' ...
             'These ROIs had missing, non-numeric, or all-NaN/Inf\n' ...
             'values across one or more subjects.\n\n' ...
             'Click OK to continue to the ROI selector.'], ...
            N_excl, strjoin(excluded_roi_names, '\n  '));
        uiwait(msgbox(excl_msg, ...
            sprintf('rDCM: %d ROI(s) Excluded by Diagnosis', N_excl), ...
            'help', 'modal'));
    end
end

fprintf('Selectable ROIs: %d\n\n', N_ROIs);

%% --- Load network labels (optional) --------------------------------
all_net_labels = local_load_net_labels(params, all_roi_names, all_roi_names_full);
% all_net_labels is either {1 x N_ROIs} cell of group strings, or {} (empty)

%% --- Initialise subset collection -----------------------------------

roi_subsets    = struct('label', {}, 'roi_names', {}, ...
                        'is_wholebrain', {}, ...
                        'created_by', {}, 'created_at', {});
subset_hashes  = {};   % for duplicate detection

%% ===================================================================
%  MAIN LOOP
%% ===================================================================

while true

    print_summary(roi_subsets);

    %% --- ROI selection ---------------------------------------------

    if use_gui
        if ~isempty(all_net_labels)
            roi_display = cellfun(@(n, r) sprintf('%s\t|| %s', n, r), ...
                all_net_labels(:)', all_roi_names(:)', 'UniformOutput', false);
            % Sort display entries by network label (stable — preserves within-group order)
            [~, disp_order] = sort(all_net_labels(:));
            roi_display_sorted = roi_display(disp_order);
        else
            roi_display_sorted = all_roi_names(:)';
            disp_order = 1:N_ROIs;
        end
        list_entries = [{'[ALL -- whole-brain]'}, roi_display_sorted];

        [sel_idx, ok] = listdlg( ...
            'PromptString', sprintf('Select ROIs for subset %d  (%d available, %d excluded by diagnosis)\n(Ctrl/Cmd+click to multi-select)', ...
                    numel(roi_subsets)+1, N_ROIs, numel(excluded_roi_names)), ...
            'SelectionMode', 'multiple', ...
            'ListString',    list_entries, ...
            'Name',          'rDCM: Define ROI Subset', ...
            'ListSize',      [340, 420], ...
            'OKString',      'Next: Name subset', ...
            'CancelString',  'Finish');

        if ~ok
            % User closed dialog or clicked Finish
            break
        end

        if isempty(sel_idx)
            msgbox('No ROIs selected. Please select at least one ROI or click Finish.', ...
                   'rDCM: No selection', 'warn');
            continue
        end

        % Check for [ALL] selection (index 1 in list_entries)
        if any(sel_idx == 1)
            is_wholebrain = true;
            roi_sel = all_roi_names;
        else
            is_wholebrain = false;
            % sel_idx - 1 indexes into the *sorted* display order; map back to all_roi_names
            original_idx = disp_order(sel_idx - 1);
            roi_sel = all_roi_names(original_idx);
        end

    else
        %% CLI fallback -----------------------------------------------
        fprintf('Enter ROI indices or names (comma-separated),\n');
        fprintf('"all" for whole-brain, or "done" to finish:\n');
        raw_input = strtrim(input('> ', 's'));

        if strcmpi(raw_input, 'done')
            break
        end

        if strcmpi(raw_input, 'all')
            is_wholebrain = true;
            roi_sel       = all_roi_names;
        else
            if isempty(raw_input)
                fprintf('  No input. Try again or type "done".\n\n');
                continue
            end
            is_wholebrain = false;
            tokens  = strtrim(strsplit(raw_input, ','));
            roi_sel = {};
            for t = 1:numel(tokens)
                tok = tokens{t};
                idx = str2double(tok);
                if ~isnan(idx) && idx == floor(idx) && idx >= 1 && idx <= N_ROIs
                    roi_sel{end+1} = all_roi_names{round(idx)}; %#ok<AGROW>
                else
                    match = find(strcmpi(all_roi_names, tok));
                    if ~isempty(match)
                        roi_sel{end+1} = all_roi_names{match(1)}; %#ok<AGROW>
                    else
                        fprintf('  WARNING: "%s" not recognised — skipped.\n', tok);
                    end
                end
            end
            if isempty(roi_sel)
                fprintf('  No valid ROIs. Try again.\n\n');
                continue
            end
        end
    end

    % Remove duplicates (preserve order)
    [~, ui] = unique(roi_sel, 'stable');
    roi_sel = roi_sel(ui);

    %% --- Duplicate subset detection --------------------------------

    subset_hash = make_hash(roi_sel);
    dup_idx     = find(strcmp(subset_hashes, subset_hash));

    if ~isempty(dup_idx)
        dup_label = roi_subsets(dup_idx(1)).label;
        msg = sprintf(['This ROI combination is identical to the already-defined ' ...
                       'subset "%s".\nPlease choose a different set of ROIs.'], dup_label);
        if use_gui
            msgbox(msg, 'rDCM: Duplicate subset', 'warn');
        else
            fprintf('\n  WARNING: %s\n\n', msg);
        end
        continue
    end

    %% --- Label input -----------------------------------------------

    default_label = sprintf('network_%d', numel(roi_subsets)+1);
    if is_wholebrain, default_label = 'whole_brain'; end

    if use_gui
        answer = inputdlg( ...
            sprintf('Enter a label for this subset (%d ROIs selected):', numel(roi_sel)), ...
            'rDCM: Name Subset', 1, {default_label});

        if isempty(answer)
            % User cancelled label dialog — go back to ROI selection
            fprintf('  Label cancelled. Returning to ROI selection.\n');
            continue
        end
        label_clean = strtrim(answer{1});
    else
        fprintf('\n  Selected %d ROI(s).\n', numel(roi_sel));
        if is_wholebrain, fprintf('  (whole-brain / all-to-all)\n'); end
        raw_label   = strtrim(input(sprintf('  Enter label [%s]: ', default_label), 's'));
        label_clean = raw_label;
        if isempty(label_clean), label_clean = default_label; end
    end

    if isempty(label_clean)
        label_clean = default_label;
    end

    % Warn if label already used (not a hard block — labels can repeat)
    existing_labels = {roi_subsets.label};
    if any(strcmp(existing_labels, label_clean))
        warn_msg = sprintf('Label "%s" is already used. Consider a unique label.', label_clean);
        if use_gui
            choice = questdlg([warn_msg '  Continue anyway?'], ...
                              'rDCM: Duplicate label', 'Yes', 'Rename', 'Rename');
            if ~strcmp(choice, 'Yes')
                continue   % go back to label input
            end
        else
            fprintf('  WARNING: %s\n', warn_msg);
        end
    end

    %% --- Confirm and store -----------------------------------------

    if use_gui
        confirm_msg = sprintf( ...
            'Add subset "%s"?\n  %d ROI(s)%s\n\nClick Yes to confirm, No to re-select.', ...
            label_clean, numel(roi_sel), ...
            sprintf('\n  First 5: %s', strjoin(roi_sel(1:min(5,end)), ', ')));
        choice = questdlg(confirm_msg, 'rDCM: Confirm Subset', 'Yes', 'No', 'Yes');
        if ~strcmp(choice, 'Yes')
            continue
        end
    else
        fprintf('\n  Subset "%s" (%d ROIs). Add? [y/n]: ', label_clean, numel(roi_sel));
        confirm = strtrim(input('', 's'));
        if ~strcmpi(confirm, 'y')
            fprintf('  Discarded.\n\n');
            continue
        end
    end

    n = numel(roi_subsets) + 1;
    roi_subsets(n).label        = label_clean;
    roi_subsets(n).roi_names    = roi_sel;
    roi_subsets(n).is_wholebrain = is_wholebrain;
    roi_subsets(n).created_by   = getenv('USER');
    roi_subsets(n).created_at   = datestr(now, 'yyyy-mm-dd HH:MM:SS');
    subset_hashes{end+1}        = subset_hash; %#ok<AGROW>

    fprintf('  --> Subset %d ("%s") added. %d ROI(s).\n\n', ...
            n, label_clean, numel(roi_sel));

    %% --- Ask to add another (GUI only — CLI loops naturally) -------
    if use_gui
        again = questdlg('Define another ROI subset?', ...
                         'rDCM: Add another?', ...
                         'Yes', 'No — save and finish', 'Yes');
        if strcmp(again, 'No — save and finish')
            break
        end
    end

end   % main while loop

%% ===================================================================
%  SAVE
%% ===================================================================

if isempty(roi_subsets)
    error('rdcm_define_roisets: no subsets were defined. Nothing saved.');
end

out_dir  = params.dirs.output;
if ~isempty(out_dir) && ~isfolder(out_dir)
    mkdir(out_dir);
end
out_path = fullfile(out_dir,params.filenames.roi_subsets);
save(out_path, 'roi_subsets', 'all_roi_names', 'excluded_roi_names');

fprintf('\n%d subset(s) saved to:\n  %s\n', numel(roi_subsets), out_path);
print_summary(roi_subsets);

end   % main function

%% ===================================================================
%  LOCAL HELPERS
%% ===================================================================

function print_summary(roi_subsets)
% Print the current list of defined subsets to console.
if isempty(roi_subsets)
    fprintf('  (No subsets defined yet.)\n\n');
    return
end
fprintf('--- Defined subsets so far (%d) ---\n', numel(roi_subsets));
for n = 1:numel(roi_subsets)
    wb_tag = '';
    if roi_subsets(n).is_wholebrain, wb_tag = '  [whole-brain]'; end
    fprintf('  %d. "%s"  —  %d ROIs%s\n', ...
            n, roi_subsets(n).label, numel(roi_subsets(n).roi_names), wb_tag);
end
fprintf('\n');
end

function h = make_hash(roi_names)
% Produce a deterministic string key from a sorted ROI name list.
% Used only for duplicate-subset detection within this session.
sorted = sort(roi_names(:)');
h      = strjoin(sorted, '|');
end

function net_labels = local_load_net_labels(params, roi_names, roi_names_full)
%% Load network labels for roi_names, using roi_names_full as the index key.
% roi_names      — the (possibly subset) list to return labels for
% roi_names_full — the full atlas list, positionally aligned to the
%                  network_labels file
net_labels = {};
if nargin < 3 || isempty(roi_names_full)
    roi_names_full = roi_names;  % fallback: treat as full list
end
if ~isfield(params, 'filenames') || ~isfield(params.filenames, 'network_labels')
    return
end
fname = params.filenames.network_labels;
if isempty(fname), return; end
fpath = fullfile(params.dirs.root, fname);
if ~isfile(fpath)
    fprintf(' NOTE: network_labels file not found — ROI list will show names only.\n %s\n\n', fpath);
    return
end
tbl = readtable(fpath, 'ReadVariableNames', false, 'Delimiter', ',');
raw = strtrim(table2cell(tbl(:, 1)));
if numel(raw) ~= numel(roi_names_full)
    warning('rdcm_define_roisets: network_labels length (%d) does not match full ROI list (%d) — skipping labels.', ...
        numel(raw), numel(roi_names_full));
    return
end
% Look up each roi_name in the full list and pull its label
[found, idx] = ismember(roi_names(:)', roi_names_full(:)');
if ~all(found)
    warning('rdcm_define_roisets: %d ROI(s) in selector not found in full ROI list — skipping labels.', ...
        sum(~found));
    return
end
net_labels = raw(idx)';  % row vector, positionally aligned to roi_names
end