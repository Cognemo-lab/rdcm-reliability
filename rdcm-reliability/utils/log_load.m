function completed = log_load(log_file)
% LOG_LOAD  Load completed job keys from the pipeline checkpoint log.
%
% Reads a tab-delimited log file where each line is:
%   <job_key>\t<timestamp>
%
% Returns a containers.Map (char → logical) for O(1) lookup.
% If the log file does not exist, returns an empty map.
%
% Usage:
%   completed = log_load('/path/to/pipeline_log.txt');

completed = containers.Map('KeyType', 'char', 'ValueType', 'logical');

if ~isfile(log_file)
    return;
end

fid = fopen(log_file, 'r');
if fid == -1
    warning('log_load: cannot open log file: %s', log_file);
    return;
end

n_loaded = 0;
while ~feof(fid)
    line = fgetl(fid);
    if ischar(line) && ~isempty(strtrim(line))
        parts   = strsplit(line, '\t');
        job_key = strtrim(parts{1});
        if ~isempty(job_key)
            completed(job_key) = true;
            n_loaded = n_loaded + 1;
        end
    end
end
fclose(fid);
end
