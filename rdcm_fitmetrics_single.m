function fit = rdcm_fitmetrics_single(output,DCM)
%% Preamble
%{
Obtains model fit metrics for an rDCM run.
---------------------------------------------------------------------------
INPUTS
---------------------------------------------------------------------------
output:=        the parameters as calculated through rDCM
DCM:=           the DCM structure
---------------------------------------------------------------------------
OUTPUTS
---------------------------------------------------------------------------
fit:=                   model fit metrics structure
 fit.cosine_sim_roi:=   per-ROI cosine similarity
 fit.cosine_sim_mean:=  mean across ROIs
 fit.R2_roi:=           per-ROI freq-domain R^2
 fit.R2_mean:=          mean across ROIs
 fit.logF:=             negative variational free energy (avg per ROI)
 fit.N_ROIs:=           stored for bookkeeping
%}

N_ROIs = DCM.n; fit.N_ROIs = N_ROIs;
N_freqs = length(output.signal.yd_source_fft) / N_ROIs;

yd_obs  = reshape(output.signal.yd_source_fft, N_freqs, N_ROIs);
yd_pred = reshape(output.signal.yd_pred_rdcm_fft, N_freqs, N_ROIs);

% Get cosine similarity using utility
[fit.cosine_sim_roi, fit.cosine_sim_mean, fit.fc_cosine] = rdcm_cosine_sim(yd_obs, yd_pred);

% Get complementary measure, R^2:
fit.R2_roi = zeros(1,N_ROIs);
for r = 1:N_ROIs
    resid    = yd_obs(:,r) - yd_pred(:,r);
    centered = yd_obs(:,r) - mean(yd_obs(:,r));
    fit.R2_roi(r) = 1 - real(resid' * resid) / real(centered' * centered);
end
fit.R2_mean = mean(fit.R2_roi, 'omitnan');

% Get mean per-region free energy
fit.logF = output.logF / N_ROIs;

end