function key = make_job_key(sub_id, step, roi_name, snr_tag)
% MAKE_JOB_KEY  Build a unique, human-readable job identifier string.
%
% Key formats:
%   Stage 1 (no ROI/SNR context):
%     sub-001__build_input
%
%   Stages 2–4 (ROI set + SNR specific):
%     sub-001__roiSet-full__SNR3__rdcm
%     sub-001__roiSet-SUD_network__SNR5__model_fit
%     sub-001__roiSet-full__SNR3__split_half
%
% Usage:
%   key = make_job_key('sub-001', 'build_input')
%   key = make_job_key('sub-001', 'rdcm', 'full', 'SNR3')
%   key = make_job_key('sub-001', 'model_fit', 'SUD_network', 'SNR5')

if nargin < 3 || isempty(roi_name)
    key = sprintf('%s__%s', sub_id, step);
else
    key = sprintf('%s__roiSet-%s__%s__%s', sub_id, roi_name, snr_tag, step);
end
end
