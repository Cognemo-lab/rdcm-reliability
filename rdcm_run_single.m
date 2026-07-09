function [EC_out,fit_out] = rdcm_run_single(Y, params, A_prior)
%% Preamble
%{
Performs rDCM on a prepped BOLD data structure Y. Wraps TAPAS functions.
Returns a clean EC struct and a fit metrics struct so that rdcm_pipeline
can store them directly in cell_EC / cell_fit without any post-processing.

The old signature was:
  [DCM, output, rdcm_options] = rdcm_run_single(Y, rdcm_options, A_prior, settings)

The new signature is:
  [EC_out, fit_out] = rdcm_run_single(Y, params, A_prior)

Migration notes
  - rdcm_options is now derived from params.rdcm (est_method, dt).
  - settings.verbose is now params.verbose.
  - A_prior is optional (default: empty → no structural prior).
  - rdcm_splithalf_single uses the same new signature; it calls this
    function with its own params struct sliced per half.
---------------------------------------------------------------------------
INPUTS
---------------------------------------------------------------------------
Y        := struct  Subject BOLD data.
             .y     [T x N_ROIs] double  BOLD timeseries
             .dt    scalar               TR in seconds
             .name  {1 x N_ROIs} cell   ROI label strings
             .subj  char                 subject identifier

params   := struct  Pipeline params (subset used here):
             .rdcm.est_method scalar  1=ridge / 2=sparse   (default: 1)
             .verbose         integer 0-3                  (default: 2)

A_prior  := [N_ROIs x N_ROIs] logical/double, or []
             Prior connectivity mask. If empty no structural prior is used.
             Self-connections (diagonal) are always forced to 1.
---------------------------------------------------------------------------
OUTPUTS
---------------------------------------------------------------------------
EC_out   := struct  Effective connectivity results.
             .A     [N_ROIs x N_ROIs] double  posterior mean EC matrix
             .Ep    struct  full TAPAS posterior means struct (Ep.A == A)
             .subj  char    propagated from Y.subj
             .N_ROIs scalar
% fit_out.cosine_sim_roi    [1 x N_ROIs]  per-ROI cosine similarity
% fit_out.cosine_sim_mean   scalar        mean across ROIs (NaN-safe)
% fit_out.cosine_sim_min    scalar        min across ROIs
% fit_out.cosine_sim_max    scalar        max across ROIs
% fit_out.N_ROIs_cosine_neg scalar        number of ROIs with cosine_sim < 0
% fit_out.logF_roi          [N_ROIs x 1]  per-ROI log model evidence
% fit_out.logF              scalar        mean log model evidence
%}

%% --- Defaults -------------------------------------------------------

if nargin < 3,  A_prior = []; end

if ~isfield(params, 'rdcm'),           params.rdcm = struct(); end
if ~isfield(params.rdcm, 'est_method'),params.rdcm.est_method = 1; end
if ~isfield(params, 'verbose'),        params.verbose = 2;     end

%% --- Build tapas_rdcm_estimate options struct -----------------------
% tapas_rdcm_estimate expects an options struct with at minimum:
%   (none)
% Other fields (filter, de-trending) use TAPAS defaults if absent.

data_type  = 'r';                  % empirical (not simulated)
est_method = params.rdcm.est_method;

%% --- Model specification -------------------------------------------

[DCM] = tapas_rdcm_model_specification(Y, [], []);

% Apply structural prior, if provided
if ~isempty(A_prior)
    DCM.a = A_prior;
    DCM.a(logical(eye(size(DCM.a)))) = 1;  % self-connections always on
end

%% --- Estimation (with timer) ----------------------------------------

if params.verbose >= 2
    fprintf('  [rdcm_run_single] Subject %s | N_ROIs=%d\n', ...
        Y.subj, size(Y.y, 2));
end

rdcm_timer = tic;

% NO need for rdcm_options when rs-fMRI mode
rdcm_options = [];

[output, ~] = tapas_rdcm_estimate(DCM, data_type, rdcm_options, est_method);
elapsed = toc(rdcm_timer);

if params.verbose >= 2
    fprintf('  [rdcm_run_single] Done (%.1fs).\n', elapsed);
end

%% --- Fit metrics ---------------------------------------------------

N_ROIs = DCM.n;
N_freqs = length(output.signal.yd_source_fft) / N_ROIs;

yd_obs  = reshape(output.signal.yd_source_fft, N_freqs, N_ROIs);
yd_pred = reshape(output.signal.yd_pred_rdcm_fft, N_freqs, N_ROIs);

% Cosine similarity
[cosine_sim_roi, cosine_sim_mean] = rdcm_cosine_sim(yd_obs, yd_pred);

% logF
logF_roi = output.logF_r;
logF_mean = output.logF/length(output.logF_r);

%% --- Posterior / prior summaries for A ------------------------------
nr = size(output.Ep.A, 1);
mu_A = output.Ep.A;
sigma2_A = nan(nr, nr);
mu0_A = nan(nr, nr);
sigma20_A = nan(nr, nr);
KL_A = nan(nr, nr);

prior_mean_A = output.priors.m0(1:nr, 1:nr);
prior_prec_A = output.priors.l0(1:nr, 1:nr);
prior_var_A = 1 ./ prior_prec_A;

for k = 1:nr
    if isempty(output.sN{k}), continue; end

    d = diag(output.sN{k});
    d = d(1:nr);

    sigma2_A(:, k) = d;
    mu0_A(:, k) = prior_mean_A(:, k);
    sigma20_A(:, k) = prior_var_A(:, k);

    ok = isfinite(mu_A(:, k)) & isfinite(mu0_A(:, k)) & ...
         isfinite(sigma2_A(:, k)) & isfinite(sigma20_A(:, k)) & ...
         sigma2_A(:, k) > 0 & sigma20_A(:, k) > 0;

    KL_A(ok, k) = 0.5 * ( ...
        log(sigma20_A(ok, k) ./ sigma2_A(ok, k)) + ...
        (sigma2_A(ok, k) + (mu_A(ok, k) - mu0_A(ok, k)).^2) ./ sigma20_A(ok, k) - 1 );
end

R_A = sigma2_A ./ prior_var_A;
D_A = (mu_A .^ 2) ./ sigma2_A;

%% --- Package EC output ---------------------------------------------

EC_out.A      = output.Ep.A;
EC_out.Ep     = output.Ep;
EC_out.subj   = Y.subj;
EC_out.N_ROIs = size(Y.y, 2);

% Cosine similarity
fit_out.cosine_sim_roi = cosine_sim_roi;
fit_out.cosine_sim_mean = cosine_sim_mean;
fit_out.cosine_sim_min     = min(cosine_sim_roi,  [], 'omitnan');
fit_out.cosine_sim_max     = max(cosine_sim_roi,  [], 'omitnan');
fit_out.N_ROIs_cosine_neg  = sum(cosine_sim_roi < 0);  % count of negative-cosine ROIs

% logF
fit_out.logF_roi = logF_roi;
fit_out.logF = logF_mean;

% posterior variance and posterior-to-prior variance ratio R_ij
fit_out.sigma2_A    = sigma2_A;    % [N x N] posterior marginal variances
fit_out.prior_var_A = prior_var_A; % [N x N] prior precisions
fit_out.R_A         = R_A;         % [N x N], bounded (0,1)
fit_out.D_A         = D_A;         % [N x N]
fit_out.KL_A        = KL_A;        % KL divergence
fit_out.mu0_A       = mu0_A;
fit_out.sigma20_A   = sigma20_A;

% Guard against zero diagonal entries in sigma2_A
fit_out.R_A(~isfinite(fit_out.D_A)) = NaN;
fit_out.D_A(~isfinite(fit_out.R_A)) = NaN;

end
