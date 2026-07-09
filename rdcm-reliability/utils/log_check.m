function result = log_check(completed, job_key)
% LOG_CHECK  Return true if job_key exists in the completed log map.
%
% Usage:
%   if ~log_check(completed, key)
%       ... run job ...
%   end

result = isKey(completed, job_key);
end
