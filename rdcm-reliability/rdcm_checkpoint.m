function [completed, ck_path] = rdcm_checkpoint(params, action, token, completed_in)
%% Preamble
%{
Machine-readable checkpoint manager for the rDCM pipeline.

Purpose
-------
Tracks which grid cells (ROI set x SNR x subject) have already completed,
so interrupted runs can resume idempotently.

Checkpoint file format
----------------------
One token per line, where each token has the form:
  rs<idx>_snr<idx>_sub<idx>
Example:
  rs2_snr3_sub07

Lines that do not match that pattern are ignored on read. This makes the
file robust to partial writes and preserves a simple append-only format.

Actions
-------
'init'  : initialize from disk or create fresh checkpoint/log files
'mark'  : append a completed token to disk and update in-memory map
'query' : legacy no-op branch; external code should use isKey(completed,...)

Important interface fix
-----------------------
The old implementation tried to update an OUTPUT variable named
'completed' in the 'mark' branch without receiving the caller's current
containers.Map as input. That broke the in-memory state and caused the
next isKey(completed, token) call in rdcm_pipeline to fail.

This revised implementation accepts the live map as the 4th input:
  completed = rdcm_checkpoint(params, 'mark', token, completed)

It is also defensive: if the 4th argument is missing or invalid, the
function rebuilds the map from disk instead of failing immediately.
---------------------------------------------------------------------------
INPUTS
---------------------------------------------------------------------------
params       struct   Pipeline params. Required fields:
                     params.log.run_id
                     params.log.log_file
                     params.log.ck_file

action       char     'init' | 'mark' | 'query'

token        char     Token string, required for 'mark'

completed_in containers.Map (optional; required for fast/clean 'mark')
---------------------------------------------------------------------------
OUTPUTS
---------------------------------------------------------------------------
completed    containers.Map  completed-token set
ck_path      char            checkpoint file path
%}

%% --- Defaults -------------------------------------------------------
if nargin < 3, token = ''; end
if nargin < 4, completed_in = []; end

%% --- Validate params ------------------------------------------------
required_fields = {'run_id', 'log_file', 'ck_file'};
for f = 1:numel(required_fields)
    if ~isfield(params, 'log') || ~isfield(params.log, required_fields{f})
        error('rdcm_checkpoint: params.log.%s is required.', required_fields{f});
    end
end

ck_path = params.log.ck_file;

%% --- Action router --------------------------------------------------
switch lower(action)

    case 'init'
        completed = local_load_or_init(params, ck_path);

    case 'mark'
        if isempty(token)
            error('rdcm_checkpoint: token must be provided for ''mark'' action.');
        end

        % Prefer the live caller map if valid; otherwise rebuild from disk.
        if isa(completed_in, 'containers.Map')
            completed = completed_in;
        else
            completed = local_read_map_from_disk(ck_path);
            rdcm_log(params, 1, ...
                'WARNING [checkpoint]: reconstructed in-memory map from disk before marking %s\n', ...
                token);
        end

        % Avoid duplicate appends if token is already present.
        if ~isKey(completed, token)
            completed(token) = true;

            fid = fopen(ck_path, 'a');
            if fid == -1
                error('rdcm_checkpoint: cannot open checkpoint file for appending: %s', ck_path);
            end
            fprintf(fid, '%s\n', token);
            fclose(fid);

            rdcm_log(params, 3, '  [checkpoint] marked: %s\n', token);
        else
            rdcm_log(params, 3, '  [checkpoint] already marked: %s\n', token);
        end

    case 'query'
        if isa(completed_in, 'containers.Map')
            completed = completed_in;
        else
            completed = local_read_map_from_disk(ck_path);
        end

    otherwise
        error('rdcm_checkpoint: unknown action ''%s''. Use init|mark|query.', action);
end

end

%% ===================================================================
function completed = local_load_or_init(params, ck_path)

if isfile(ck_path)
    completed = local_read_map_from_disk(ck_path);
    rdcm_log(params, 1, ...
        'RESUME: checkpoint loaded (%d completed indices). Run ID: %s\n', ...
        completed.Count, params.log.run_id);
    rdcm_log(params, 1, '------------------------------------------------------------\n');
    rdcm_log(params, 1, 'RESUMED at %s — %d indices already complete.\n', ...
        datestr(now, 'yyyy-mm-dd HH:MM:SS'), completed.Count);
    rdcm_log(params, 1, '------------------------------------------------------------\n');
else
    completed = containers.Map('KeyType', 'char', 'ValueType', 'logical');

    fid = fopen(ck_path, 'w');
    if fid == -1
        error('rdcm_checkpoint: cannot create checkpoint file: %s', ck_path);
    end
    fprintf(fid, '# rDCM checkpoint | run_id: %s | created: %s\n', ...
        params.log.run_id, datestr(now, 'yyyy-mm-dd HH:MM:SS'));
    fclose(fid);

    fid2 = fopen(params.log.log_file, 'w');
    if fid2 == -1
        error('rdcm_checkpoint: cannot create log file: %s', params.log.log_file);
    end
    fprintf(fid2, '# rDCM pipeline log | run_id: %s | created: %s\n', ...
        params.log.run_id, datestr(now, 'yyyy-mm-dd HH:MM:SS'));
    fclose(fid2);

    rdcm_log(params, 1, 'FRESH RUN. Run ID: %s\n', params.log.run_id);
end

end

%% ===================================================================
function completed = local_read_map_from_disk(ck_path)

completed = containers.Map('KeyType', 'char', 'ValueType', 'logical');

if ~isfile(ck_path)
    return
end

raw = fileread(ck_path);
if isempty(raw)
    return
end

lines = regexp(raw, '\r\n|\n|\r', 'split');
for i = 1:numel(lines)
    line_i = strtrim(lines{i});
    if isempty(line_i)
        continue
    end
    % Accept both 2D format (rs<i>_sub<j>) and legacy 3D format (rs<i>_snr<j>_sub<k>)
    if ~isempty(regexp(line_i, '^rs\d+_sub\d+$', 'once')) || ...
       ~isempty(regexp(line_i, '^rs\d+_snr\d+_sub\d+$', 'once'))
        completed(line_i) = true;
    end
end

end
