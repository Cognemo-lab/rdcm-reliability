function [cosine_sim_roi, cosine_sim_mean] = rdcm_cosine_sim(yd_obs, yd_pred)
%% Preamble
%{
Computes cosine similarity between observed and predicted BOLD temporal
derivatives in the frequency domain, per ROI and as a mean.

This is the corrected implementation. The previous version in
cognemo_getEC_withCosine used real(yd_obs)' * real(yd_pred), which
discards the imaginary component of the complex FFT vectors and
systematically overestimates similarity. The correct treatment for
complex vectors is to use the conjugate dot product and take the real
part of the result, since cosine similarity must be real-valued.

The formula applied per ROI r is:

    cos_r = real(yd_obs(:,r)' * yd_pred(:,r)) 
            / (norm(yd_obs(:,r)) * norm(yd_pred(:,r)))

where ' is the conjugate transpose (Hermitian), so the numerator is
the real part of the Hermitian inner product. This is equivalent to
the standard cosine similarity generalised to complex vectors, and
reduces to the familiar real formula when inputs are real-valued.

Note: norm() in MATLAB operates on the full complex vector by default,
i.e. norm(z) = sqrt(real(z'*z)) = sqrt(sum(|z_i|^2)), which is the
correct L2 norm for complex vectors.
---------------------------------------------------------------------------
INPUTS
---------------------------------------------------------------------------
yd_obs:=        [N_freqs x N_ROIs] complex matrix
                Observed BOLD temporal derivative in frequency domain.
                From output.signal.yd_source_fft, reshaped from vector.
yd_pred:=       [N_freqs x N_ROIs] complex matrix  
                Predicted BOLD temporal derivative in frequency domain.
                From output.signal.yd_pred_rdcm_fft, reshaped from vector.
---------------------------------------------------------------------------
OUTPUTS
---------------------------------------------------------------------------
cosine_sim_roi  := [1 x N_ROIs] cosine similarity per ROI, in [-1, 1]
cosine_sim_mean := scalar, mean across ROIs (NaN-safe)
fc_cosine       := scalar, FC-space cosine similarity between the
                   observed and predicted inter-regional covariance
                   structure (upper-triangle CSD-derived FC matrices).
                   Measures whether the EC solution reproduces the
                   between-ROI covariance pattern in frequency space.
                   Returns NaN if N_ROIs < 2.
%}

%% --- Input Validation -----------------------------------------------

if ~isequal(size(yd_obs), size(yd_pred))
    error('rdcm_cosine_sim: yd_obs and yd_pred must have the same size.');
end

if isempty(yd_obs)
    error('rdcm_cosine_sim: inputs are empty.');
end

[N_freqs, N_ROIs] = size(yd_obs);

if N_freqs < 2
    error('rdcm_cosine_sim: at least 2 frequency bins required (got %d).', N_freqs);
end

%% --- Per-ROI Cosine Similarity --------------------------------------

cosine_sim_roi = nan(1, N_ROIs);

for r = 1:N_ROIs
    obs_r  = yd_obs(:, r);
    pred_r = yd_pred(:, r);

    norm_obs  = norm(obs_r);
    norm_pred = norm(pred_r);

    % Skip degenerate ROIs (zero-norm signal or prediction).
    % This can occur if a ROI timeseries is flat (all-NaN after
    % preprocessing) or if rDCM produced a zero-weight prediction.
    if norm_obs < eps || norm_pred < eps
        % cosine_sim_roi(r) remains NaN — will be excluded from mean
        continue
    end

    % Hermitian (conjugate) inner product, real part taken since
    % cosine similarity is a real-valued similarity measure.
    cosine_sim_roi(r) = real(obs_r' * pred_r) / (norm_obs * norm_pred);
end

%% --- Mean (NaN-safe) ------------------------------------------------

cosine_sim_mean = mean(cosine_sim_roi, 'omitnan');


end
