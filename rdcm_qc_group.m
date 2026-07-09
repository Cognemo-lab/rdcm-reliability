function qc = rdcm_qc_group(ec_storage, fit_storage, sh_storage, params)
%% Preamble
%{
Computes group-level QC for a single ROI set and assigns three subject-
inclusion logical vectors for downstream analyses.

All reliability metrics are computed here from the raw split-half EC
matrices in sh_storage and the full-session EC matrices in ec_storage.
No QC logic is performed inside the rDCM run loop.

Decision flow (applied per subject, in order):
  Stage 1 — Hard gates (absolute, leakage-free):
    Gate 1a: cosine_sim_mean < 0       → Tier 0 (excluded)
    Gate 1b: min(icc_s_diag_masked, icc_s_offdiag_masked) < icc_floor → Tier 1 (excluded)
  Stage 2 — Soft tier (absolute thresholds):
    icc_s_strong < icc_s_pass          → Tier 2 (marginal, group only)
    icc_s_strong >= icc_s_pass         → Tier 3 (pass, group + individual)
  Stage 3 — Advisory flags (never cause exclusion):
    Adv A: logF anomalously low (> 2 SD below group mean)
    Adv B: icc_s_strong below group median - k*MAD
    Adv C: bias index (icc_s_offdiag - pearson_r_offdiag) < -0.10
    Adv D: EC matrix Mahalanobis outlier (robust PCA, p=0.025)
    Adv E: cross-network tier discordance — populated externally by
           the pipeline orchestrator after all ROI sets are processed;
           stored in qc.adv_crossnet for each subject.
  Stage 4 — Posterior informativeness (dataset-level context only):
    R, sqrt(D), obs-to-param ratio
---------------------------------------------------------------------------
INPUTS
---------------------------------------------------------------------------
ec_storage := {1 x N_subs} cell of EC_out structs from rdcm_run_single.
  Each struct has:
    .A      [N x N]   full-session EC matrix
    .subj   char

fit_storage := {1 x N_subs} cell of fit_out structs from rdcm_run_single.
  Each struct has:
    .cosine_sim_mean  scalar
    .cosine_sim_roi   [N x 1]
    .logF             scalar
    .logF_roi         [N x 1]
    .R_A              [N x N]   posterior/prior variance ratio (optional)
    .D_A              [N x N]   mu^2/sigma^2 SNR index       (optional)

sh_storage := {1 x N_subs} cell of sh structs from rdcm_splithalf_single.
  Each struct has:
    .A1             [N x N]   EC matrix, first half
    .A2             [N x N]   EC matrix, second half
    .offdiag.icc    scalar    ICC_s off-diagonal (Spearman-Brown corrected)
    .offdiag.pearson_r scalar Pearson r, off-diagonal
    .diag.icc       scalar    ICC_s diagonal
    .subj           char
    .N_TRs_half     scalar
    .N_ROIs         scalar

params := pipeline params struct. Relevant subfields:
  .qc.strength_mask_method  'top_k' | 'top_quartile' (default: 'top_quartile')
  .qc.strength_mask_k       integer (default: 1000, used when method='top_k')
  .qc.strength_mask_p       scalar in [0,1] (default: 0.75)
  .qc.strength_mask_min_connections  integer (default: 20)
  .qc.icc_s_fail            scalar (default: 0.40)  — Gate 1b + Tier 1 floor
  .qc.icc_s_pass            scalar (default: 0.60)  — Tier 3 threshold
  .qc.icc_c_fail            scalar (default: 0.40)
  .qc.icc_c_pass            scalar (default: 0.60)
  .qc.mad_k                 scalar (default: 2)     — Advisory B multiplier
  .qc.bias_threshold        scalar (default: -0.10) — Advisory C threshold
  .qc.adv_d_pca_variance    scalar (default: 0.90)  — PCA variance for Adv D
  .qc.adv_d_alpha           scalar (default: 0.025) — Mahal. tail for Adv D
  .qc.n_perm_asymmetry      integer (default: 1000) — permutations for EC/FC check

---------------------------------------------------------------------------
OUTPUTS
---------------------------------------------------------------------------
qc := struct

  Subject-level inclusion vectors ([N_subs x 1] logical):
    .include_all          non-hard-gated (Tier >= 1)
    .include_group        group-level analysis (Tier >= 2)
    .include_individual   individual-level analysis (Tier == 3)

  Subject-level tier and flags:
    .subject_tier     [N_subs x 1]  0=hard-gate, 1=fail, 2=marginal, 3=pass
    .hard_gate_flag   [N_subs x 1]  logical — any hard gate triggered
    .hard_gate_type   {N_subs x 1}  'cosine'|'icc_floor'|'motion'|
                                    'cosine+icc_floor'|'none'
    .exclusion_reason {N_subs x 1}  human-readable string, '' if no exclusion

  Per-subject ICC (three variants, all Spearman-Brown corrected):
    .icc_s_diag     [N_subs x 1]   diagonal (self-connections)
    .icc_s_offdiag  [N_subs x 1]   all off-diagonal
    .icc_s_strong   [N_subs x 1]   strong off-diagonal (tier criterion)
    .pearson_r      [N_subs x 1]   Pearson r, off-diagonal split-half
    .bias_index     [N_subs x 1]   icc_s_offdiag - pearson_r (Advisory C)
    .kl_mask_all_s / .kl_mask_offdiag_s / .kl_mask_diag_s {N_subs x 1} [N x N] logical per-subject KL masks

  Per-connection ICC ([N x N] matrices):
    .icc_c            [N x N]   ICC(3,1) SB-corrected (NaN on diagonal)
    .icc_c_tier       [N x N]   1/2/3 per connection
    .icc_c_mean_all   scalar    mean across all off-diagonal
    .icc_c_mean_kl_masked scalar   mean across strength-masked only
    .kl_mask_c [N x N] group-level KL mask (ranked by group-mean full-run KL_A)

  Fit metrics:
    .cosine_sim   [N_subs x 1]
    .logF         [N_subs x 1]

  Advisory flags ([N_subs x 1] logical unless noted):
    .adv_logF           Adv A — logF > 2 SD below group mean
    .adv_reliability    Adv B — icc_s_strong below group median - k*MAD
    .adv_bias           Adv C — bias_index < bias_threshold
    .adv_ec_outlier     Adv D — Mahalanobis outlier in PCA-reduced EC space
    .adv_crossnet       Adv E — cross-network discordance (set externally)
    .n_advisories       [N_subs x 1]  count of Adv A-D (E added externally)
    .review_recommended [N_subs x 1]  true when n_advisories >= 2

  Posterior informativeness (group-mean, off-diagonal):
    .R_A_group        [N x N]   group-mean R per connection
    .D_A_group        [N x N]   group-mean D per connection
    .R_mean_offdiag       scalar
    .R_mean_masked    scalar
    .sqrtD_mean_offdiag   scalar    sqrt(D_mean) — posterior Z-score
    .sqrtD_mean_masked scalar

  EC vs FC informational gain (3-question checklist):
    .ec_vs_fc.self_conn_icc_gain     scalar  mean(icc_s_diag) - mean(icc_s_offdiag)
    .ec_vs_fc.self_conn_between_sub_var  scalar  var of mean A_diag across subjects
    .ec_vs_fc.self_conn_within_roi_var   scalar  mean var of A_diag across ROIs per sub
    .ec_vs_fc.self_conn_informative  logical  true if icc gain > 0 AND between > within
    .ec_vs_fc.delta_ij               [M x 1]  asymmetry indices for strength-masked pairs
    .ec_vs_fc.delta_pvalue           scalar   permutation p-value (observed vs null width)
    .ec_vs_fc.directionality_nontrivial logical
    .ec_vs_fc.uppertri_lowertri_corr [N_subs x 1]  upper/lower triangle correlation per sub
    .ec_vs_fc.symmetry_mean_r        scalar
    .ec_vs_fc.ec_preferred           logical  true if >= 2 of 3 questions are yes

  Summary counts:
    .N_subs, .N_pass, .N_marginal, .N_fail, .N_hard_gate, .prop_usable

  Parameters used:
    .params_used    struct (for manifest / reporting)
%}

%% --- Parameter defaults ----------------------------------------------

if ~isfield(params, 'qc'), params.qc = struct(); end

kl_mask_method = get_param(params.qc, 'kl_mask_method', 'top_k');
kl_mask_k = get_param(params.qc, 'kl_mask_k', 1000);
kl_mask_p = get_param(params.qc, 'kl_mask_p', 0.75);
kl_mask_min_conn = get_param(params.qc, 'kl_mask_min_connections', 20);
kl_mask_obs_param_ratio = get_param(params.qc, 'kl_mask_obs_param_ratio', 1);
icc_floor = get_param(params.qc, 'icc_floor', 0.40);
icc_s_pass = get_param(params.qc, 'icc_s_pass', 0.60);
outlier_z = get_param(params.qc, 'outlier_z', 3.5);

icc_c_fail          = get_param(params.qc, 'icc_c_fail',                  0.40);
icc_c_pass          = get_param(params.qc, 'icc_c_pass',                  0.60);
mad_k               = get_param(params.qc, 'mad_k',                       2);
bias_thr            = get_param(params.qc, 'bias_threshold',              -0.10);
adv_d_pca_var       = get_param(params.qc, 'adv_d_pca_variance',          0.90);
adv_d_alpha         = get_param(params.qc, 'adv_d_alpha',                 0.025);
n_perm_asym         = get_param(params.qc, 'n_perm_asymmetry',            1000);

N_subs = numel(fit_storage);

%% --- Determine N_ROIs from first valid entry -------------------------

N = [];
for i = 1:N_subs
    if ~isempty(ec_storage{i}) && isfield(ec_storage{i}, 'A')
        N = size(ec_storage{i}.A, 1);
        break;
    end
end
if isempty(N)
    error('rdcm_qc_group: no valid EC matrices found in ec_storage.');
end

offdiag_mask = ~logical(eye(N));
diag_idx     = logical(eye(N));

%% --- Collect fit metrics ---------------------------------------------

cosine_sim = nan(N_subs, 1);
logF       = nan(N_subs, 1);

for i = 1:N_subs
    if isempty(fit_storage{i}), continue; end
    cosine_sim(i) = fit_storage{i}.cosine_sim_mean;
    logF(i)       = fit_storage{i}.logF;
end

%% --- Stage 1: Hard gates ---------------------------------------------

% Gate 1a: spectral inversion failure
gate_cosine = cosine_sim < 0;

% Gate 1b: absolute ICC floor (computed after per-subject ICC below;
%          populated into hard_gate_flag after ICC loop)

%% --- Per-subject ICC and Pearson r -----------------------------------

icc_s_all = nan(N_subs, 1);
icc_s_offdiag = nan(N_subs, 1);
icc_s_diag = nan(N_subs, 1);
icc_s_all_masked = nan(N_subs, 1);
icc_s_offdiag_masked = nan(N_subs, 1);
icc_s_diag_masked = nan(N_subs, 1);
pearson_r = nan(N_subs, 1);
pearson_r_masked = nan(N_subs, 1);
kl_mask_all_s = cell(N_subs, 1);
kl_mask_offdiag_s = cell(N_subs, 1);
kl_mask_diag_s = cell(N_subs, 1);

for i = 1:N_subs
    if isempty(sh_storage{i}) || ~isfield(sh_storage{i}, 'A1'), continue; end
    if isempty(fit_storage{i}) || ~isfield(fit_storage{i}, 'KL_A'), continue; end
    
    A1 = sh_storage{i}.A1;
    A2 = sh_storage{i}.A2;
    
    % Ranking source: full-run KL (stable, single ranking per subject).
    % Sizing source: split-half obs/param ratio (the binding, more
    % data-constrained regime), via N_TRs_half from sh_storage.
    KL_for_rank = fit_storage{i}.KL_A;
    N_TRs_half_i = sh_storage{i}.N_TRs_half;
    
    % Unmasked ICC — all / offdiag / diag
    icc_s_all(i) = icc_sb(A1(:), A2(:));
    icc_s_offdiag(i) = icc_sb(A1(offdiag_mask), A2(offdiag_mask));
    icc_s_diag(i) = icc_sb(A1(diag_idx), A2(diag_idx));
    
    % Pearson r (off-diagonal, unmasked) — used for bias index
    r_mat = corrcoef(A1(offdiag_mask), A2(offdiag_mask));
    if ~any(isnan(r_mat(:))), pearson_r(i) = r_mat(1,2); end
    
    % Size k from the split-half obs/param regime, not the full run
    k_i = kl_mask_k_from_obsparam(N_TRs_half_i, N, kl_mask_obs_param_ratio, kl_mask_min_conn);
    
    % Diagonal and off-diagonal are ranked and sized WITHIN their own
    % partition, not pooled — a shared top-k over the full N x N matrix lets
    % the N(N-1)-large off-diagonal pool dominate and starves the N-entry
    % diagonal pool of any selected connections (this was causing
    % icc_s_diag_masked to be NaN for nearly all subjects). make_kl_mask
    % already caps k to the partition size internally, so passing k_i to
    % both calls is safe: for the diagonal it simply resolves to "select all
    % N self-connections" whenever k_i >= N, which is correct since the
    % diagonal has a far smaller parameter budget than the off-diagonal.
    mask_off = make_kl_mask(KL_for_rank, offdiag_mask, kl_mask_method, k_i, kl_mask_p, kl_mask_min_conn);
    mask_dia = make_kl_mask(KL_for_rank, diag_idx, kl_mask_method, k_i, kl_mask_p, kl_mask_min_conn);
    mask_all = mask_off | mask_dia;
    
    kl_mask_all_s{i} = mask_all;
    kl_mask_offdiag_s{i} = mask_off;
    kl_mask_diag_s{i} = mask_dia;
    
    % Masked ICC — all / offdiag / diag
    icc_s_all_masked(i) = icc_sb(A1(mask_all), A2(mask_all));
    icc_s_offdiag_masked(i) = icc_sb(A1(mask_off), A2(mask_off));
    icc_s_diag_masked(i) = icc_sb(A1(mask_dia), A2(mask_dia));
    
    % Masked Pearson r — used for bias index on masked connections
    r_mat_m = corrcoef(A1(mask_off), A2(mask_off));
    if ~any(isnan(r_mat_m(:))), pearson_r_masked(i) = r_mat_m(1,2); end
end

% Bias index computed on the masked off-diagonal set (policy: KL-masked only)
bias_index = icc_s_offdiag_masked - pearson_r_masked;

%% --- Gate 1b: absolute ICC floor (now that ICC is computed) ----------
% Gate 1b: absolute reliability floor — min(diag_masked, offdiag_masked),
% per agreed policy: self-connections must not be allowed to pass if they
% are weaker than the off-diagonal parameters.
tier_icc = min(icc_s_diag_masked, icc_s_offdiag_masked);

gate_icc_floor = tier_icc < icc_floor;
gate_icc_floor(isnan(tier_icc)) = true;
gate_icc_floor(gate_cosine) = false;

hard_gate_flag = gate_cosine | gate_icc_floor;

n_mask_counts = cellfun(@(m) sum(m(:)), kl_mask_all_s, 'UniformOutput', true);
valid_for_mask = ~hard_gate_flag & ~isnan(tier_icc);
if any(valid_for_mask)
    n_mask_median = median(n_mask_counts(valid_for_mask));
else
    n_mask_median = 0;
end
if n_mask_median < kl_mask_min_conn
    rdcm_log(params, 1, ...
    'rdcm_qc_group: KL mask median=%d connections (< min %d).\n', ...
    round(n_mask_median), kl_mask_min_conn);
end

%% --- Stage 2: Outlier flag (KL-masked values only, per policy) -------

outlier_z_score = nan(N_subs, 1);
valid_outlier = ~hard_gate_flag & ~isnan(icc_s_all_masked);
if sum(valid_outlier) >= 4
    center = median(icc_s_all_masked(valid_outlier));
    scale = 1.4826 * mad(icc_s_all_masked(valid_outlier), 1);
    if scale > 0
        outlier_z_score(valid_outlier) = ...
        (icc_s_all_masked(valid_outlier) - center) / scale;
    end
end
outlier_flag = abs(outlier_z_score) > outlier_z;
outlier_flag(isnan(outlier_flag)) = false;
outlier_flag(hard_gate_flag) = false;

%% --- Stage 3: Tier assignment ------------------------------------------
% Tier 0 = hard-gated, 1 = absolute-floor fail, 2 = outlier/marginal, 3 = pass

subject_tier = nan(N_subs, 1);
for i = 1:N_subs
    if hard_gate_flag(i)
        subject_tier(i) = 0;
    elseif gate_icc_floor(i)
        subject_tier(i) = 1;
    elseif outlier_flag(i) || tier_icc(i) < icc_s_pass
        subject_tier(i) = 2;
    else
        subject_tier(i) = 3;
    end
end

% Inclusion vectors
include_all        = subject_tier >= 1;
include_group      = subject_tier >= 2;
include_individual = subject_tier == 3;

%% --- Group-level ICC_s summary (Fisher-z averaged, include_group only) ---

ig = include_group;
icc_s_summary = struct();
icc_s_fields  = {'icc_s_all','icc_s_offdiag','icc_s_diag', ...
                 'icc_s_all_masked','icc_s_offdiag_masked','icc_s_diag_masked'};
icc_s_vectors = {icc_s_all, icc_s_offdiag, icc_s_diag, ...
                 icc_s_all_masked, icc_s_offdiag_masked, icc_s_diag_masked};

for fi = 1:numel(icc_s_fields)
    v = icc_s_vectors{fi}(ig);
    icc_s_summary.(icc_s_fields{fi}).mean = icc_fisher_mean(v);
    icc_s_summary.(icc_s_fields{fi}).min  = min(v, [], 'omitnan');
    icc_s_summary.(icc_s_fields{fi}).max  = max(v, [], 'omitnan');
end

% tier_icc and Pearson r are also correlation-bounded — Fisher-z mean
icc_s_summary.tier_icc.mean         = icc_fisher_mean(tier_icc(ig));
icc_s_summary.tier_icc.min          = min(tier_icc(ig), [], 'omitnan');
icc_s_summary.tier_icc.max          = max(tier_icc(ig), [], 'omitnan');

icc_s_summary.pearson_r_masked.mean = icc_fisher_mean(pearson_r_masked(ig));
icc_s_summary.pearson_r_masked.min  = min(pearson_r_masked(ig), [], 'omitnan');
icc_s_summary.pearson_r_masked.max  = max(pearson_r_masked(ig), [], 'omitnan');

% bias_index is a difference of two bounded quantities, not itself a
% correlation coefficient — arithmetic mean remains appropriate here
icc_s_summary.bias_index.mean = mean(bias_index(ig), 'omitnan');
icc_s_summary.bias_index.min  = min(bias_index(ig), [], 'omitnan');
icc_s_summary.bias_index.max  = max(bias_index(ig), [], 'omitnan');

%% --- Gate type strings and exclusion reasons -------------------------

hard_gate_type  = repmat({'none'}, N_subs, 1);
exclusion_reason = repmat({''}, N_subs, 1);

for i = 1:N_subs
    parts = {};
    if gate_cosine(i)
        parts{end+1} = sprintf('Gate1a: cosine_sim_mean=%.3f < 0', cosine_sim(i)); %#ok
    end
    if gate_icc_floor(i)
        parts{end+1} = sprintf('Gate1b: min(icc_diag_masked,icc_offdiag_masked)=%.3f < floor=%.2f', ...
        tier_icc(i), icc_floor); %#ok
    end
    if ~isempty(parts)
        exclusion_reason{i} = strjoin(parts, '; ');
    end

    % Gate type label
    c = gate_cosine(i); b = gate_icc_floor(i);
    if c && b, hard_gate_type{i} = 'cosine+icc_floor';
    elseif c, hard_gate_type{i} = 'cosine';
    elseif b, hard_gate_type{i} = 'icc_floor';
    end
end

%% --- Stage 3: Advisory flags -----------------------------------------
% All advisories are computed on non-hard-gated subjects only.
% Group statistics (MAD) are used only for Advisory B; this is acceptable
% because Advisory B never causes exclusion — it is a flag for human review.

valid = ~hard_gate_flag & ~isnan(icc_s_all_masked);

% Advisory A: logF anomalously low (> 2 SD below group mean)
logF_mean_grp = mean(logF(valid), 'omitnan');
logF_std_grp  = std(logF(valid),  0, 'omitnan');
adv_logF      = logF < (logF_mean_grp - 2 * logF_std_grp);
adv_logF(hard_gate_flag) = false;

% Advisory B is now redundant with the Stage 2 outlier_flag (both use
% KL-masked ICC_s_all with robust z-scoring), so it is aliased directly.
adv_reliability = outlier_flag;

% Advisory C: bias index below threshold
adv_bias = bias_index < bias_thr;
adv_bias(hard_gate_flag | isnan(bias_index)) = false;

% Advisory D: EC matrix Mahalanobis outlier in PCA-reduced space
% Uses off-diagonal entries of A_full; robust PCA via SVD on centred data;
% Mahalanobis distance in top-K PC space vs chi2 threshold.
adv_ec_outlier = false(N_subs, 1);
try
    idx_valid = find(valid);
    n_valid   = numel(idx_valid);
    n_offdiag = N * (N - 1);

    if n_valid >= 5  % minimum for meaningful PCA
        A_mat = nan(n_valid, n_offdiag);
        for ki = 1:n_valid
            Ai = ec_storage{idx_valid(ki)}.A;
            mask_i = kl_mask_offdiag_s{idx_valid(ki)};
            Ai_vals = Ai(offdiag_mask)';
            % Zero out entries not in this subject's own off-diagonal KL mask so the
            % PCA outlier check reflects only well-determined connections
            keep_i = mask_i(offdiag_mask)';
            Ai_vals(~keep_i) = NaN;
            A_mat(ki,:) = Ai_vals;
        end
        % Remove columns with no variance (constant across subjects)
        col_var = var(A_mat, 0, 1, 'omitnan');
        keep_cols = col_var > 0;
        A_sub = A_mat(:, keep_cols);

        % Centre
        A_cen = A_sub - mean(A_sub, 1);

        % SVD-based PCA
        [~, S, V] = svd(A_cen, 'econ');
        eigvals = diag(S).^2 / (n_valid - 1);
        cum_var = cumsum(eigvals) / sum(eigvals);
        K = find(cum_var >= adv_d_pca_var, 1, 'first');
        if isempty(K), K = size(V, 2); end
        K = max(K, 2);  % at least 2 PCs

        scores = A_cen * V(:, 1:K);
        mu_pc  = mean(scores, 1);
        sig_pc = cov(scores);

        % Mahalanobis distance per subject
        mah = nan(n_valid, 1);
        for ki = 1:n_valid
            diff = scores(ki,:) - mu_pc;
            mah(ki) = diff / sig_pc * diff';
        end

        % Chi-squared threshold: chi2inv(1-alpha, K)
        % Approximated without Statistics Toolbox:
        chi2_thr = chi2_approx(K, 1 - adv_d_alpha);
        outlier_flags = mah > chi2_thr;
        adv_ec_outlier(idx_valid(outlier_flags)) = true;
    end
catch ME
    rdcm_log(params, 1, ...
        'rdcm_qc_group: Advisory D (EC outlier) failed: %s\n', ME.message);
end

% Advisory E: cross-network tier discordance — placeholder, set externally
adv_crossnet = false(N_subs, 1);

% Advisory count and review flag (A-D only; E added by orchestrator)
n_advisories = double(adv_logF) + double(adv_reliability) + ...
               double(adv_bias) + double(adv_ec_outlier);
review_recommended = n_advisories >= 2;

%% --- Per-connection ICC ----------------------------------------------

icc_c           = nan(N);
icc_c_diag      = nan(N,1);
icc_c_tier_mat  = nan(N);
mean_A_group    = nan(N);
kl_mask_c = false(N);
icc_c_mean_all  = NaN;
icc_c_mean_diag = NaN;
icc_c_mean_kl_masked = NaN;

self_conn_gain_icc_c = NaN;
self_conn_gain_icc_c_pct = NaN;

idx_group = find(include_group)';
n_grp = numel(idx_group);
idx_group_no_outlier = find(include_group & ~outlier_flag)';
n_grp_no_outlier = numel(idx_group_no_outlier);

if n_grp >= 3
    A1_stack    = nan(N, N, n_grp);
    A2_stack    = nan(N, N, n_grp);
    A_full_stack = nan(N, N, n_grp);

    for ki = 1:n_grp
        ii = idx_group(ki);
        if ~isempty(sh_storage{ii}) && isfield(sh_storage{ii}, 'A1')
            A1_stack(:,:,ki)    = sh_storage{ii}.A1;
            A2_stack(:,:,ki)    = sh_storage{ii}.A2;
        end
        if ~isempty(ec_storage{ii}) && isfield(ec_storage{ii}, 'A')
            A_full_stack(:,:,ki) = ec_storage{ii}.A;
        end
    end

    for r = 1:N
        for c = 1:N
            if r == c, continue; end
            v1c = squeeze(A1_stack(r,c,:));
            v2c = squeeze(A2_stack(r,c,:));
            ok  = isfinite(v1c) & isfinite(v2c);
            if sum(ok) < 3, continue; end
            icc_c(r,c) = icc_sb(v1c(ok), v2c(ok));
        end
    end

    for r = 1:N
        v1c = squeeze(A1_stack(r,r,:));
        v2c = squeeze(A2_stack(r,r,:));
        ok = isfinite(v1c) & isfinite(v2c);
        if sum(ok) < 3, continue; end
        icc_c_diag(r) = icc_sb(v1c(ok), v2c(ok));
    end

    icc_c_tier_mat(icc_c <  icc_c_fail & offdiag_mask) = 1;
    icc_c_tier_mat(icc_c >= icc_c_fail & icc_c < icc_c_pass & offdiag_mask) = 2;
    icc_c_tier_mat(icc_c >= icc_c_pass & offdiag_mask) = 3;

    % Group-level KL mask — ranks connections by group-mean full-run KL
    % divergence (consistency with the subject-level design: "well-determined"
    % is defined by KL, not by raw connection strength). Sized via the group's
    % median N_TRs_half, the same data-constrained regime used per-subject.
    mean_A_group = mean(A_full_stack, 3, 'omitnan');
    
    KL_stack = nan(N, N, n_grp);
    N_TRs_half_grp = nan(n_grp, 1);
    for ki = 1:n_grp
    ii = idx_group(ki);
    if isfield(fit_storage{ii}, 'KL_A') && ~isempty(fit_storage{ii}.KL_A)
    KL_stack(:,:,ki) = fit_storage{ii}.KL_A;
    end
    if ~isempty(sh_storage{ii}) && isfield(sh_storage{ii}, 'N_TRs_half')
    N_TRs_half_grp(ki) = sh_storage{ii}.N_TRs_half;
    end
    end
    mean_KL_group = mean(KL_stack, 3, 'omitnan');
    
    k_grp = kl_mask_k_from_obsparam(median(N_TRs_half_grp, 'omitnan'), N, ...
        kl_mask_obs_param_ratio, kl_mask_min_conn);
    kl_mask_c = make_kl_mask(mean_KL_group, offdiag_mask, ...
        kl_mask_method, k_grp, kl_mask_p, kl_mask_min_conn);

    % Summary means
    icc_c_vals = icc_c(offdiag_mask);
    icc_c_mean_all = icc_fisher_mean(icc_c_vals);
    icc_c_mean_diag = icc_fisher_mean(icc_c_diag);
    icc_c_vals_kl_masked = icc_c(kl_mask_c);
    icc_c_mean_kl_masked = icc_fisher_mean(icc_c_vals_kl_masked);

    self_conn_gain_icc_c = icc_c_mean_diag - icc_c_mean_kl_masked;
    self_conn_gain_icc_c_pct = 100 * self_conn_gain_icc_c / icc_c_mean_kl_masked;

end

%% --- Per-connection ICC, outliers excluded (for sensitivity check) ---

icc_c_mean_all_no_outlier = NaN;
if n_grp_no_outlier >= 3
    A1_stack2 = nan(N, N, n_grp_no_outlier);
    A2_stack2 = nan(N, N, n_grp_no_outlier);
    for ki = 1:n_grp_no_outlier
        ii = idx_group_no_outlier(ki);
        if ~isempty(sh_storage{ii}) && isfield(sh_storage{ii}, 'A1')
            A1_stack2(:,:,ki) = sh_storage{ii}.A1;
            A2_stack2(:,:,ki) = sh_storage{ii}.A2;
        end
    end
    icc_c2 = nan(N);
    for r = 1:N
        for c = 1:N
            if r == c, continue; end
            v1c = squeeze(A1_stack2(r,c,:));
            v2c = squeeze(A2_stack2(r,c,:));
            ok = isfinite(v1c) & isfinite(v2c);
            if sum(ok) < 3, continue; end
            icc_c2(r,c) = icc_sb(v1c(ok), v2c(ok));
        end
    end
    vals2 = icc_c2(offdiag_mask);
    icc_c_mean_all_no_outlier = icc_fisher_mean(vals2);
end

%% --- Posterior informativeness (R and D) ----------------------------

R_stack = nan(N, N, max(n_grp,1));
D_stack = nan(N, N, max(n_grp,1));

for ki = 1:n_grp
    ii = idx_group(ki);
    if isfield(fit_storage{ii}, 'R_A') && ~isempty(fit_storage{ii}.R_A)
        R_stack(:,:,ki) = fit_storage{ii}.R_A;
    end
    if isfield(fit_storage{ii}, 'D_A') && ~isempty(fit_storage{ii}.D_A)
        D_stack(:,:,ki) = fit_storage{ii}.D_A;
    end
end

R_A_group = mean(R_stack, 3, 'omitnan');
D_A_group = mean(D_stack, 3, 'omitnan');

R_offdiag = R_A_group(offdiag_mask);
R_diag = R_A_group(diag_idx);
D_offdiag = D_A_group(offdiag_mask);
D_diag = D_A_group(diag_idx);

R_mean_offdiag     = mean(R_offdiag,                 'omitnan');
R_mean_masked      = mean(R_A_group(kl_mask_c),      'omitnan');
R_mean_diag        = mean(R_diag,                    'omitnan');
sqrtD_mean_offdiag = mean(sqrt(D_offdiag),           'omitnan');
sqrtD_mean_masked  = mean(sqrt(D_A_group(kl_mask_c)),'omitnan');
sqrtD_mean_diag    = mean(sqrt(D_diag),              'omitnan');

%% --- EC vs FC informational gain -------------

ec_vs_fc = struct();

% --- Directionality: paired t-test per connection pair, FDR-BH -------

ec_vs_fc.direction_pvals_unmasked = [];
ec_vs_fc.direction_pvals_masked = [];
ec_vs_fc.direction_pct_sig_unmasked = NaN;
ec_vs_fc.direction_pct_sig_masked = NaN;
ec_vs_fc.directionality_nontrivial = false;

if n_grp >= 5
    mean_A_sub = nan(N, N, n_grp);
    for ki = 1:n_grp
        ii = idx_group(ki);
        if ~isempty(ec_storage{ii}) && isfield(ec_storage{ii}, 'A')
            mean_A_sub(:,:,ki) = ec_storage{ii}.A;
        end
    end

    [ri, ci] = find(triu(true(N), 1));
    n_pairs_total = numel(ri);
    pvals_unmasked = nan(n_pairs_total, 1);

    for pi = 1:n_pairs_total
        a_fwd = squeeze(mean_A_sub(ri(pi), ci(pi), :));
        a_rev = squeeze(mean_A_sub(ci(pi), ri(pi), :));
        d = a_fwd - a_rev;
        d = d(isfinite(d));
        if numel(d) >= 5 && std(d) > 0
            pvals_unmasked(pi) = ttest_manual(d);
        end
    end

    valid_p = isfinite(pvals_unmasked);
    if any(valid_p)
        q_unmasked = fdr_bh(pvals_unmasked(valid_p));
        ec_vs_fc.direction_pvals_unmasked = pvals_unmasked(valid_p);
        ec_vs_fc.direction_pct_sig_unmasked = 100 * sum(q_unmasked < 0.05) / sum(valid_p);
    end

    in_mask = kl_mask_c(sub2ind([N,N], ri, ci)) | kl_mask_c(sub2ind([N,N], ci, ri));
    pvals_masked = pvals_unmasked(in_mask);
    valid_pm = isfinite(pvals_masked);
    if any(valid_pm)
        q_masked = fdr_bh(pvals_masked(valid_pm));
        ec_vs_fc.direction_pvals_masked = pvals_masked(valid_pm);
        ec_vs_fc.direction_pct_sig_masked = 100 * sum(q_masked < 0.05) / sum(valid_pm);
    end

    ec_vs_fc.directionality_nontrivial = ec_vs_fc.direction_pct_sig_masked > 0;
end

% --- Question 3: Is the EC matrix meaningfully asymmetric? -----------
% Pearson correlation between upper and lower triangles per subject.
% Low mean r across subjects → asymmetry is consistent (FC would miss it).

uplo_corr = nan(N_subs, 1);
uplo_corr_masked = nan(N_subs, 1);
for i = 1:N_subs
    if ~valid(i) || isempty(ec_storage{i}), continue; end
    A_i = ec_storage{i}.A;
    upper_idx = triu(true(N), 1);
    lower_idx = tril(true(N), -1);
    upper = A_i(upper_idx);
    lower = A_i(lower_idx);
    r_mat = corrcoef(upper, lower);
    if ~any(isnan(r_mat(:))), uplo_corr(i) = r_mat(1,2); end
    
    if ~isempty(kl_mask_offdiag_s{i})
        m = kl_mask_offdiag_s{i};
        upper_m = A_i(upper_idx & m);
        lower_full = A_i'; % transpose to align (j,i) with (i,j) mask positions
        lower_m = lower_full(upper_idx & m);
        if numel(upper_m) >= 5
            r_mat_m = corrcoef(upper_m, lower_m);
            if ~any(isnan(r_mat_m(:))), uplo_corr_masked(i) = r_mat_m(1,2); end
        end
    end
end
ec_vs_fc.uppertri_lowertri_corr = uplo_corr;
ec_vs_fc.uppertri_lowertri_corr_masked = uplo_corr_masked;
ec_vs_fc.symmetry_mean_r = mean(uplo_corr(valid), 'omitnan');
ec_vs_fc.symmetry_mean_r_masked = mean(uplo_corr_masked(valid), 'omitnan');
% Meaningful asymmetry if mean upper/lower r is substantially below 1
% and directionality_nontrivial confirmed
ec_vs_fc.asymmetry_informative = ...
    ec_vs_fc.directionality_nontrivial && ...
    (ec_vs_fc.symmetry_mean_r < 0.90);

% --- Overall recommendation ------------------------------------------
n_yes = double(ec_vs_fc.directionality_nontrivial) + ...
    double(ec_vs_fc.asymmetry_informative);
ec_vs_fc.n_yes = n_yes;
ec_vs_fc.n_questions = 2;
ec_vs_fc.ec_preferred = (n_yes >= 2);

%% --- Summary counts --------------------------------------------------

N_hard_gate = sum(subject_tier == 0);
N_fail      = sum(subject_tier == 1);
N_marginal  = sum(subject_tier == 2);
N_pass      = sum(subject_tier == 3);
prop_usable = (N_pass + N_marginal) / N_subs;

%% --- Package output --------------------------------------------------

qc.include_all        = include_all;
qc.include_group      = include_group;
qc.include_individual = include_individual;

qc.icc_s_summary = icc_s_summary;

qc.subject_tier      = subject_tier;
qc.hard_gate_flag    = hard_gate_flag;
qc.hard_gate_type    = hard_gate_type;
qc.exclusion_reason  = exclusion_reason;

qc.kl_mask_n_median = round(n_mask_median);

qc.icc_s_all = icc_s_all;
qc.icc_s_offdiag = icc_s_offdiag;
qc.icc_s_diag = icc_s_diag;
qc.icc_s_all_masked = icc_s_all_masked;
qc.icc_s_offdiag_masked = icc_s_offdiag_masked;
qc.icc_s_diag_masked = icc_s_diag_masked;
qc.tier_icc = tier_icc;
qc.pearson_r = pearson_r;
qc.pearson_r_masked = pearson_r_masked;
qc.bias_index = bias_index;
qc.kl_mask_all_s = kl_mask_all_s;
qc.kl_mask_offdiag_s = kl_mask_offdiag_s;
qc.kl_mask_diag_s = kl_mask_diag_s;
qc.outlier_flag = outlier_flag;
qc.outlier_z_score = outlier_z_score;
qc.icc_c_mean_all_no_outlier = icc_c_mean_all_no_outlier;

qc.icc_c            = icc_c;
qc.icc_c_tier       = icc_c_tier_mat;
qc.icc_c_mean_all   = icc_c_mean_all;
qc.icc_c_diag = icc_c_diag;
qc.icc_c_mean_diag = icc_c_mean_diag;
qc.icc_c_mean_kl_masked = icc_c_mean_kl_masked;
qc.self_conn_gain_icc_c = self_conn_gain_icc_c;
qc.self_conn_gain_icc_c_pct = self_conn_gain_icc_c_pct;
qc.kl_mask_c  = kl_mask_c;
qc.mean_A_group     = mean_A_group;

qc.cosine_sim = cosine_sim;
qc.logF       = logF;

qc.adv_logF          = adv_logF;
qc.adv_reliability   = adv_reliability;
qc.adv_bias          = adv_bias;
qc.adv_ec_outlier    = adv_ec_outlier;
qc.adv_crossnet      = adv_crossnet;   % placeholder; set by orchestrator
qc.n_advisories      = n_advisories;
qc.review_recommended = review_recommended;

qc.R_A_group       = R_A_group;
qc.D_A_group       = D_A_group;
qc.R_mean_offdiag  = R_mean_offdiag;
qc.R_mean_masked   = R_mean_masked;
qc.R_mean_diag     = R_mean_diag;
qc.sqrtD_mean_offdiag  = sqrtD_mean_offdiag;
qc.sqrtD_mean_masked   = sqrtD_mean_masked;
qc.sqrtD_mean_diag     = sqrtD_mean_diag;

qc.ec_vs_fc = ec_vs_fc;

qc.N_subs       = N_subs;
qc.N_pass       = N_pass;
qc.N_marginal   = N_marginal;
qc.N_fail       = N_fail;
qc.N_hard_gate  = N_hard_gate;
qc.prop_usable  = prop_usable;

qc.params_used.kl_mask_method = kl_mask_method;
qc.params_used.kl_mask_k = kl_mask_k;
qc.params_used.kl_mask_p = kl_mask_p;
qc.params_used.kl_mask_min_connections = kl_mask_min_conn;
qc.params_used.icc_floor = icc_floor;
qc.params_used.outlier_z = outlier_z;
qc.params_used.icc_s_pass                    = icc_s_pass;
qc.params_used.icc_c_fail                    = icc_c_fail;
qc.params_used.icc_c_pass                    = icc_c_pass;
qc.params_used.mad_k                         = mad_k;
qc.params_used.bias_threshold                = bias_thr;
qc.params_used.adv_d_pca_variance            = adv_d_pca_var;
qc.params_used.adv_d_alpha                   = adv_d_alpha;
qc.params_used.n_perm_asymmetry              = n_perm_asym;

end % main function

%% ====================================================================
%  Local helpers
%% ====================================================================

function v = get_param(s, field, default)
if isfield(s, field) && ~isempty(s.(field))
    v = s.(field);
else
    v = default;
end
end

% --------------------------------------------------------------------------
function icc = icc_sb(v1, v2)
% ICC(3,1) with Spearman-Brown split-half correction.
v1 = v1(:); v2 = v2(:);
ok = isfinite(v1) & isfinite(v2);
v1 = v1(ok); v2 = v2(ok);
if numel(v1) < 3
    icc = NaN; return
end
raw = icc31(v1, v2);
if isnan(raw)
    icc = NaN;
else
    icc = 2 * raw / (1 + raw);
end
end

% --------------------------------------------------------------------------
function icc = icc31(v1, v2)
% ICC(3,1) — two-way mixed effects, consistency, single measures.
% Shrout & Fleiss (1979).
n = numel(v1); k = 2;
Y = [v1(:), v2(:)];
grand_mean   = mean(Y(:));
subject_means = mean(Y, 2);
rater_means   = mean(Y, 1);
SSb = k * sum((subject_means - grand_mean).^2);
SSw = sum(sum((Y - subject_means).^2));
SSr = n * sum((rater_means - grand_mean).^2);
SSe = SSw - SSr;
MSb = SSb / (n - 1);
MSe = SSe / ((n-1) * (k-1));
denom = MSb + (k-1) * MSe;
if denom < eps
    icc = NaN;
else
    icc = (MSb - MSe) / denom;
end
end

% --------------------------------------------------------------------------
function m = icc_fisher_mean(vals)
% Fisher-z averaging for ICC/correlation-like values. Prevents bias from
% skewed distributions near the [-1,1] boundaries when taking a mean
% across many connection-wise (or subject-wise) coefficients.
vals = vals(:);
vals = vals(isfinite(vals));
if isempty(vals)
    m = NaN;
    return;
end
vals = min(max(vals, -0.999999), 0.999999); % avoid atanh(+-1) = Inf
z = atanh(vals);
m = tanh(mean(z));
end

% --------------------------------------------------------------------------
function p = ttest_manual(d)
% One-sample t-test (d vs 0), two-tailed, no Statistics Toolbox dependency.
d = d(:); d = d(isfinite(d));
n = numel(d);
if n < 2, p = NaN; return; end
m = mean(d); s = std(d);
if s == 0, p = NaN; return; end
t = m / (s / sqrt(n));
df = n - 1;
p = betainc(df / (df + t^2), df/2, 0.5); % base MATLAB, no toolbox needed
end

% --------------------------------------------------------------------------
function q = fdr_bh(p)
% Benjamini-Hochberg FDR correction.
p = p(:);
n = numel(p);
[p_sorted, idx] = sort(p);
q_sorted = p_sorted .* n ./ (1:n)';
q_sorted = min(q_sorted, 1);
for i = n-1:-1:1
    q_sorted(i) = min(q_sorted(i), q_sorted(i+1));
end
q = nan(n, 1);
q(idx) = q_sorted;
end

% --------------------------------------------------------------------------
function mask = make_kl_mask(KL_A, region_mask, method, k, p, min_conn)
% Builds a logical mask selecting the best-determined connections within
% region_mask (e.g. offdiag_mask, diag_idx, or true(N)), ranked by KL
% divergence (higher KL = posterior moved further from prior = more
% informative / better-determined).
vals = KL_A(region_mask);
n = numel(vals);
mask = false(size(KL_A));
if n == 0, return; end
switch lower(method)
    case 'top_k'
        k_use = min(k, n);
        [~, ord] = sort(vals, 'descend');
        ord = ord(isfinite(vals(ord)));
        k_use = min(k_use, numel(ord));
    case 'top_quartile'
        thr = quantile(vals(isfinite(vals)), p);
        region_idx = find(region_mask);
        keep = vals >= thr & isfinite(vals);
        mask(region_idx(keep)) = true;
        if sum(mask(:)) < min_conn
            [~, ord] = sort(vals, 'descend');
            ord = ord(isfinite(vals(ord)));
            k_use = min(min_conn, numel(ord));
            mask(:) = false;
            mask(region_idx(ord(1:k_use))) = true;
        end
        return;
    otherwise
        warning('rdcm_qc_group:unknownKLMaskMethod', ...
        'Unknown kl_mask_method "%s". Using top_k.', method);
        k_use = min(k, n);
        [~, ord] = sort(vals, 'descend');
        ord = ord(isfinite(vals(ord)));
        k_use = min(k_use, numel(ord));
    end
    region_idx = find(region_mask);
    mask(region_idx(ord(1:k_use))) = true;
    if sum(mask(:)) < min_conn && numel(ord) >= min_conn
        mask(:) = false;
        mask(region_idx(ord(1:min_conn))) = true;
    end
end

function k = kl_mask_k_from_obsparam(N_TRs_half, N_ROIs, obs_param_ratio, min_conn)
% Targets a mask size such that the retained connections satisfy the
% desired obs/param ratio under the SPLIT-HALF observation count — the
% more data-constrained regime, since that is what binds estimability
% for the split-half ICC computation.
%
% Observations available per half: N_TRs_half (per ROI, frequency-domain
% samples after rDCM's per-region regression scale similarly with N_TRs).
% Total connections in A: N_ROIs^2.
%
% k = number of connections whose combined parameter budget stays under
% the obs/param target, given N_TRs_half observations.
max_estimable = floor(N_TRs_half / obs_param_ratio);
k = max(min_conn, min(max_estimable, N_ROIs^2));
end

% --------------------------------------------------------------------------
function x = chi2_approx(k, p)
% Approximate chi2 inverse CDF (Wilson-Hilferty method).
% Avoids Statistics Toolbox dependency.
% k: degrees of freedom; p: cumulative probability (e.g. 0.975)
% Standard normal quantile via rational approximation (Abramowitz & Stegun)
z = normal_inv(p);
h = 1 - 2/(9*k);
x = k * (h + z * sqrt(2/(9*k)))^3;
x = max(x, 0);
end

function z = normal_inv(p)
% Rational approximation to the standard normal inverse CDF.
% Abramowitz & Stegun 26.2.17
if p <= 0 || p >= 1, z = Inf * sign(p - 0.5); return; end
q = min(p, 1-p);
t = sqrt(-2 * log(q));
c = [2.515517, 0.802853, 0.010328];
d = [1.432788, 0.189269, 0.001308];
z = t - (c(1) + c(2)*t + c(3)*t^2) / (1 + d(1)*t + d(2)*t^2 + d(3)*t^3);
if p < 0.5, z = -z; end
end
