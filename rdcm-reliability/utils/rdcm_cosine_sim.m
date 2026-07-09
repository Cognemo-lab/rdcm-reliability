function [cs_roi, cs_mean] = rdcm_cosine_sim(yd_obs_fft, yd_pred_fft)
% RDCM_COSINE_SIM  Frequency-domain cosine similarity for rDCM model fit.
%
% Computes the cosine similarity between observed and model-predicted BOLD
% temporal derivative signals in the frequency domain — the native domain
% in which rDCM performs its regression.
%
% FORMULA
%   For complex vectors a and b, the correct cosine similarity is:
%
%       cos(a, b) = Re(a' * b) / (||a|| * ||b||)
%
%   where (') is the conjugate transpose and ||.|| is the complex L2 norm:
%       ||a|| = sqrt(sum(|a_i|^2))
%
%   Taking Re() of the inner product preserves both magnitude and phase.
%   The previous implementation used real(a) and real(b) before the dot
%   product, discarding all phase information (imaginary parts), which
%   caused systematic overestimation of fit quality — worst for poorly
%   fitting subjects. This function corrects that error.
%
% INPUTS
%   yd_obs_fft   [N_freqs × N_ROIs] or [N × 1]
%                  FFT of observed BOLD temporal derivative.
%                  Typically output.signal.yd_source_fft from TAPAS rDCM.
%
%   yd_pred_fft  [N_freqs × N_ROIs] or [N × 1]
%                  FFT of model-predicted BOLD temporal derivative.
%                  Typically output.signal.yd_pred_rdcm_fft from TAPAS rDCM.
%
% OUTPUTS
%   cs_roi       [1 × N_ROIs]  Per-ROI cosine similarity. NaN for
%                              zero-variance (degenerate) ROIs.
%                              Returns a scalar if inputs are vectors.
%
%   cs_mean      Scalar. Nanmean of cs_roi across ROIs.
%
% INTERPRETATION
%   cs = 1.0  Perfect model fit for this ROI.
%   cs ~ 0    Model prediction uncorrelated with observed signal.
%   cs < 0    Anti-phase prediction (very poor fit).
%
%   Use cs_mean as a per-subject QC metric. Subjects below the group
%   distribution (e.g. > 2 MAD below median) should be flagged for review.
%   Do not apply fixed thresholds — calibrate against your dataset.
%
% EXAMPLE
%   [cs_roi, cs_mean] = rdcm_cosine_sim( ...
%       output.signal.yd_source_fft, ...
%       output.signal.yd_pred_rdcm_fft);
%
% SEE ALSO
%   p3_model_fit.m, p2_run_rdcm.m

% --- Input validation ---
if ~isequal(size(yd_obs_fft), size(yd_pred_fft))
    error('rdcm_cosine_sim: inputs must be the same size. Got [%s] vs [%s].', ...
          num2str(size(yd_obs_fft)), num2str(size(yd_pred_fft)));
end

% --- Handle vector input (whole-brain concatenated) ---
if isvector(yd_obs_fft)
    yd_obs_fft  = yd_obs_fft(:);
    yd_pred_fft = yd_pred_fft(:);
    is_vector   = true;
else
    is_vector   = false;
end

% --- Per-ROI cosine similarity ---
N_ROIs = size(yd_obs_fft, 2);
cs_roi = NaN(1, N_ROIs);

for r = 1:N_ROIs
    a   = yd_obs_fft(:, r);
    b   = yd_pred_fft(:, r);
    denom = norm(a) * norm(b);

    if denom < eps
        % Degenerate: zero-variance signal in observed or predicted
        cs_roi(r) = NaN;
    else
        % Correct complex cosine similarity — conjugate transpose
        cs_roi(r) = real(a' * b) / denom;
    end
end

cs_mean = mean(cs_roi, 'omitnan');

% Return scalar if input was a vector
if is_vector
    cs_roi = cs_roi(1);
end
end
