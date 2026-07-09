function rdcm_log(params, level, fmt, varargin)
%% Preamble
%{
Unified logging utility for the rDCM pipeline. All pipeline functions
call rdcm_log instead of fprintf directly, so that verbosity and output
routing are controlled centrally through params.

Verbosity tiers (params.verbose):
  0  Silent    — nothing written (errors still throw normally)
  1  Minimal   — pipeline start/end, step completions, warnings, summaries
  2  Normal    — per-subject progress, QC pass/fail per subject (default)
  3  Verbose   — per-ROI fit metrics, ICC per connection category, timing

Output routing (params.log_mode):
  'console'   — fprintf to command window only (local/debug runs)
  'file'      — write to log file only (cluster runs)
  'both'      — write to both (default)

The log file path is read from params.log.log_file. If params.log_mode is
'console' or params.log is not set, no file I/O is performed.

Usage:
  rdcm_log(params, 1, 'Pipeline started: %s\n', run_id)
  rdcm_log(params, 2, '  Subject %d/%d: %s\n', i, N, subj)
  rdcm_log(params, 3, '    ROI %d cosine: %.4f\n', r, c)
---------------------------------------------------------------------------
INPUTS
---------------------------------------------------------------------------
params:=    struct — pipeline params. Fields read:
              .verbose    integer 0-3 (default: 2)
              .log_mode   char 'console'|'file'|'both' (default: 'both')
              .log.log_file  char — full path to human-readable log file
                               (required when log_mode is 'file' or 'both')
level:=     integer — minimum verbosity level required to emit this message
fmt:=       char   — format string (same syntax as fprintf)
varargin:=  optional format arguments
%}

%% --- Defaults -------------------------------------------------------

if ~isfield(params, 'verbose'),  params.verbose  = 2;       end
if ~isfield(params, 'log_mode'), params.log_mode = 'both';  end

% Skip immediately if message level exceeds requested verbosity
if level > params.verbose
    return
end

%% --- Format Message -------------------------------------------------

% Always run through sprintf so escape sequences (\n, \t) are interpreted
% even when there are no format arguments.
if isempty(varargin)
    msg = sprintf(fmt);
else
    msg = sprintf(fmt, varargin{:});
end

% Prepend timestamp for file output (not console — timestamps clutter
% interactive use but are essential for post-hoc cluster log inspection)
timestamp     = datestr(now, 'yyyy-mm-dd HH:MM:SS');
msg_with_ts   = ['[' timestamp '] ' msg];

%% --- Route Output ---------------------------------------------------

write_console = strcmpi(params.log_mode, 'console') || ...
                strcmpi(params.log_mode, 'both');
write_file    = strcmpi(params.log_mode, 'file') || ...
                strcmpi(params.log_mode, 'both');

% Console output (no timestamp — matches normal MATLAB fprintf feel)
if write_console
    fprintf('%s', msg);
end

% File output (with timestamp, append mode)
if write_file
    if ~isfield(params, 'log') || ~isfield(params.log, 'log_file')
        % Warn once via console regardless of log_mode, then skip file write
        persistent warned_no_logfile
        if isempty(warned_no_logfile)
            fprintf(['[rdcm_log] WARNING: log_mode includes ''file'' but ' ...
                     'params.log.log_file is not set. File logging disabled.\n']);
            warned_no_logfile = true;
        end
        return
    end

    fid = fopen(params.log.log_file, 'a');
    if fid == -1
        fprintf(['[rdcm_log] WARNING: Cannot open log file: %s\n' ...
                 '          Check path and permissions.\n'], params.log.log_file);
        return
    end
    fprintf(fid, '%s', msg_with_ts);
    fclose(fid);
end

end