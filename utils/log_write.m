function log_write(log_file, job_key)
% LOG_WRITE  Append a completed job key and timestamp to the log file.
%
% The log file is opened in append mode so partial writes during a crash
% cannot corrupt already-completed entries. The file is closed immediately
% after each write.
%
% Usage:
%   log_write('/path/to/pipeline_log.txt', 'sub-001__build_input');

fid = fopen(log_file, 'a');
if fid == -1
    error('log_write: cannot open log file for appending: %s', log_file);
end
fprintf(fid, '%s\t%s\n', job_key, datestr(now, 'yyyy-mm-dd HH:MM:SS'));
fclose(fid);
end
