function rdcm_report_qc(results, params, roi_subsets)
%% Preamble
%{
Generates a human-readable QC report summarising model fit and split-half
reliability across all ROI sets.

The report is written as a plain-text file to params.dirs.output/export/
and also printed to the console/log at verbosity level 1.

Report structure:
  HEADER         run ID, date
  HOW TO INTERPRET  decision-flow, gates 1a/1b/1c, tier table, advisories
  Per ROI set:
    1. Acquisition context
    2. Subject summary (tier counts, inclusion vectors, gate breakdown)
    3. Advisory flag summary (A-D, co-occurrence, review-recommended)
    4. Reliability summary (ICC_s six variants: {all,offdiag,diag} x {unmasked,masked}, Pearson r, bias index)
    5. Per-connection ICC summary + tier breakdown
    6. Posterior informativeness (R, sqrt(D))
    7. EC vs FC checklist (three questions + recommendation)
    8. Per-subject table (non-hard-gated)
    9. Hard-gated subject list with exclusion reasons
   10. Feature-set CSV export

The function does NOT produce figures. Visualisation is done separately
via rdcm_plot_qc_distributions using results.qc directly.
---------------------------------------------------------------------------
INPUTS
---------------------------------------------------------------------------
results     := struct from rdcm_pipeline
  .qc       {N_roisets x 1} cell of qc structs from rdcm_qc_group
  .grid     {N_roisets x N_subs} cell of EC_out structs
  .fit      {N_roisets x N_subs} cell of fit_out structs
  .params   struct
  .data_meta struct

params      := struct — pipeline params (log routing, thresholds, dirs)

roi_subsets := struct array from rdcm_define_roisets
%}

%% --- Setup -----------------------------------------------------------

if ~isfield(params, 'sub_id_pattern')
    params.sub_id_pattern = 'sub-[A-Za-z0-9]+';
end

report_dir = fullfile(params.dirs.output, 'export');
if ~isfolder(report_dir), mkdir(report_dir); end

if isfield(params, 'log') && isfield(params.log, 'run_id')
    run_id = params.log.run_id;
elseif isfield(params, 'run_id')
    run_id = params.run_id;
else
    run_id = 'unknown';
    warning('rdcm_report_qc: run_id not found; using "unknown".');
end

report_path = fullfile(report_dir, sprintf('qc_report_%s.txt', run_id));
fid = fopen(report_path, 'w');
if fid == -1
    error('rdcm_report_qc: cannot write report: %s', report_path);
end

sep_major = repmat('=', 1, 72);
sep_minor = repmat('-', 1, 72);

% Nested wlog writes to file AND console simultaneously
function wlog(fmt, varargin)
    line = sprintf(fmt, varargin{:});
    fprintf(fid, '%s', line);
    rdcm_log(params, 1, fmt, varargin{:});
end

% Read QC thresholds (used throughout)
icc_floor = get_param_safe(params, 'qc', 'icc_floor', 0.40);
icc_s_pass = get_param_safe(params, 'qc', 'icc_s_pass', 0.60);
icc_c_fail = get_param_safe(params, 'qc', 'icc_c_fail', 0.40);
icc_c_pass = get_param_safe(params, 'qc', 'icc_c_pass', 0.60);
outlier_z = get_param_safe(params, 'qc', 'outlier_z', 3.5);
bias_thr = get_param_safe(params, 'qc', 'bias_threshold', -0.10);

kmeth = get_param_safe(params, 'qc', 'kl_mask_method', 'top_k');
k_p = get_param_safe(params, 'qc', 'kl_mask_p', 0.75);
k_ratio = get_param_safe(params, 'qc', 'kl_mask_obs_param_ratio', 10);
k_min_conn = get_param_safe(params, 'qc', 'kl_mask_min_connections', 20);

%% --- Header ----------------------------------------------------------

wlog('%s\n', sep_major);
wlog('rDCM PIPELINE QC REPORT\n');
wlog('Run ID : %s\n', run_id);
wlog('Date   : %s\n', datestr(now, 'yyyy-mm-dd HH:MM:SS')); %#ok<TNOW1,DATST>
wlog('%s\n\n', sep_major);

wlog('HOW TO READ THIS REPORT\n');
wlog('%s\n', sep_minor);
wlog('\n');
wlog('This report answers two practical questions for each ROI set:\n');
wlog('1) How many subjects are usable?\n');
wlog('2) How stable are the estimated connections?\n');
wlog('Use the AT A GLANCE table first. Refer to the final section only if you need\n');
wlog('the exact meaning of terms such as reliability, better-determined connections,\n');
wlog('or the thresholds used.\n');
wlog('\n');
wlog('%s\n\n', sep_major);

%% --- AT A GLANCE -------------------------------------------------------

wlog('AT A GLANCE\n');
wlog('%s\n', sep_minor);
wlog('\n');

usable_thr = 0.70;
if isfield(params, 'qc') && isfield(params.qc, 'report') && ...
        isfield(params.qc.report, 'percent_usable_threshold') && ...
        ~isempty(params.qc.report.percent_usable_threshold)
    usable_thr = params.qc.report.percent_usable_threshold;
end

N_roisets = numel(roi_subsets);
glance_n_exclude   = nan(N_roisets, 1);
glance_n_pass      = nan(N_roisets, 1);
glance_n_marginal  = nan(N_roisets, 1);
glance_prop_usable = nan(N_roisets, 1);
glance_conn_reliab = nan(N_roisets, 1);

for i_rs = 1:N_roisets
    qc_i = results.qc{i_rs};
    glance_n_exclude(i_rs)   = qc_i.N_hard_gate + qc_i.N_fail;
    glance_n_pass(i_rs)      = qc_i.N_pass;
    glance_n_marginal(i_rs)  = qc_i.N_marginal;
    glance_prop_usable(i_rs) = qc_i.prop_usable;

    offdiag_i  = ~logical(eye(numel(roi_subsets(i_rs).roi_names)));
    n_masked_i = sum(qc_i.kl_mask_c(offdiag_i));
    if n_masked_i > 0
        n_t2_i = sum(qc_i.icc_c_tier(qc_i.kl_mask_c) == 2, 'omitnan');
        n_t3_i = sum(qc_i.icc_c_tier(qc_i.kl_mask_c) == 3, 'omitnan');
        glance_conn_reliab(i_rs) = 100 * (n_t2_i + n_t3_i) / n_masked_i;
    end
end

% Sort: EXCLUDE ascending, then PASS descending
[~, sort_idx] = sortrows([glance_n_exclude, -glance_n_pass]);

wlog(' %-20s %10s %10s %10s %14s %s\n', ...
    'ROI SET', 'PASS', 'REVIEW', 'EXCLUDE', 'C-R?', 'U-O?');

for k = 1:N_roisets
    i_rs = sort_idx(k);
    if isnan(glance_conn_reliab(i_rs))
        conn_str = 'n/a';
    else
        conn_str = sprintf('%.1f%%', glance_conn_reliab(i_rs));
    end
    if glance_prop_usable(i_rs) >= usable_thr
        ok_str = 'YES';
    else
        ok_str = 'NO';
    end
    wlog(' %-20s %10d %10d %10d %14s %s\n', ...
        roi_subsets(i_rs).label, glance_n_pass(i_rs), glance_n_marginal(i_rs), ...
        glance_n_exclude(i_rs), conn_str, ok_str);
end

wlog('\n');
wlog(' U-O -> USABLE OVERALL? = YES when the percent of usable subjects is >= %.0f%%.\n', usable_thr * 100);
wlog(' C-R -> CONNECTIONS RELIABLE? = percent of better-determined connections with\n');
wlog(' at least marginal reliability. See final section for exact definition\n');
wlog(' Table sorted by EXCLUDE (low to high), then PASS (high to low).\n');
wlog('\n');

%% --- SUBJECTS TO CHECK (top ROI set only) ------------------------

best_rs = sort_idx(1);

qc_best = results.qc{best_rs};
N_subs_best = qc_best.N_subs;
if isfield(qc_best, 'subject_ids') && numel(qc_best.subject_ids) == N_subs_best
    sub_ids_best = qc_best.subject_ids;
else
    sub_ids_best = arrayfun(@(k) sprintf('sub_%03d', k), 1:N_subs_best, ...
        'UniformOutput', false);
end

wlog('SUBJECTS TO CHECK (%s)\n', roi_subsets(best_rs).label);
wlog('%s\n', sep_minor);
wlog('\n');

wlog(' DO NOT USE\n');
wlog('   Subjects excluded from all analyses\n');
n_red = 0;
for i = 1:N_subs_best
    if qc_best.subject_tier(i) == 0
        wlog('   %-24s reason: %s\n', sub_ids_best{i}, qc_best.exclusion_reason{i});
        n_red = n_red + 1;
    elseif qc_best.subject_tier(i) == 1
        wlog('   %-24s reason: reliability below floor (tier_icc=%.3f)\n', ...
            sub_ids_best{i}, qc_best.tier_icc(i));
        n_red = n_red + 1;
    end
end
wlog('   (N=%d total)\n\n', n_red);

wlog(' REVIEW BEFORE USE\n');
wlog('   Subjects usable for group analysis, but check notes below\n');
n_yellow = 0;
for i = 1:N_subs_best
    if qc_best.subject_tier(i) == 2
        adv_parts = {};
        if qc_best.adv_logF(i), adv_parts{end+1} = plain_reason('logF low'); end %#ok
        if qc_best.adv_reliability(i), adv_parts{end+1} = plain_reason('reliability outlier'); end %#ok
        if qc_best.adv_bias(i), adv_parts{end+1} = plain_reason('bias index'); end %#ok
        if qc_best.adv_ec_outlier(i), adv_parts{end+1} = plain_reason('EC outlier'); end %#ok
        if isempty(adv_parts), adv_parts = {plain_reason('marginal tier_icc')}; end
        wlog('   %-24s reason: %s\n', sub_ids_best{i}, strjoin(adv_parts, '; '));
        n_yellow = n_yellow + 1;
    end
end
wlog('   (N=%d total)\n\n', n_yellow);

n_green = sum(qc_best.subject_tier == 3);
wlog(' CLEAR TO USE\n');
wlog('   Subjects acceptable for individual-level analysis: %d\n\n', n_green);
wlog('%s\n\n', sep_major);

%% --- APPENDIX: detailed per-ROI-set results ---------------------------

wlog('APPENDIX — DETAILED RESULTS PER ROI SET\n');
wlog('%s\n\n', sep_major);

for i_rs = 1:N_roisets

    rs       = roi_subsets(i_rs);
    rs_label = rs.label;
    rs_field = matlab.lang.makeValidName(rs_label);
    N_ROIs   = numel(rs.roi_names);
    qc_i     = results.qc{i_rs};
    N_subs_i = qc_i.N_subs;

    % Subject IDs
    if isfield(qc_i, 'subject_ids') && numel(qc_i.subject_ids) == N_subs_i
        sub_ids = qc_i.subject_ids;
    else
        sub_ids = arrayfun(@(k) sprintf('sub_%03d', k), 1:N_subs_i, ...
            'UniformOutput', false);
    end

    wlog('%s\n', sep_major);
    wlog('ROI SET: %s  (%d ROIs)\n', rs_label, N_ROIs);
    wlog('%s\n', sep_minor);

    %% 1. Acquisition context
    local_print_acq_context(@wlog, params, N_ROIs, N_subs_i, sep_minor);

    %% 2. Subject summary
    wlog(' SUBJECT SUMMARY\n');
    wlog(' %s\n', sep_minor);
    wlog('   Total subjects            : %d\n', N_subs_i);
    wlog('\n');
    wlog('   Hard-gated (Tier 0)       : %d\n', qc_i.N_hard_gate);

    % Gate-type breakdown
    n_cos = sum(cellfun(@(x) contains(x,'cosine'),    qc_i.hard_gate_type));
    n_icc = sum(cellfun(@(x) contains(x,'icc_floor'), qc_i.hard_gate_type));
    if qc_i.N_hard_gate > 0
        wlog('     Gate 1a (cosine < 0)   : %d\n', n_cos);
        wlog('     Gate 1b (ICC < floor)  : %d\n', n_icc);
    end
    wlog('\n');
    wlog(' Tier 1 — Fail : %d [tier_icc < %.2f]\n', ...
        qc_i.N_fail, icc_floor);
    wlog(' Tier 2 — Marginal : %d [outlier OR tier_icc < %.2f]\n', ...
        qc_i.N_marginal, icc_s_pass);
    wlog(' Tier 3 — Pass : %d [tier_icc >= %.2f, not outlier]\n', ...
        qc_i.N_pass, icc_s_pass);
    wlog('   Prop. usable (Tier 2+3)   : %.1f%%\n', qc_i.prop_usable * 100);
    wlog('\n');
    wlog('   include_all               : %d subjects\n', sum(qc_i.include_all));
    wlog('   include_group             : %d subjects\n', sum(qc_i.include_group));
    wlog('   include_individual        : %d subjects\n', sum(qc_i.include_individual));

    %% 3. Advisory flag summary
    wlog('\n');
    wlog(' ADVISORY FLAGS\n');
    wlog(' %s\n', sep_minor);
    valid = ~qc_i.hard_gate_flag;
    n_v   = sum(valid);
    wlog('   (non-hard-gated N = %d)\n\n', n_v);

    wlog('   Adv A  logF anomalously low           : %d  (%.1f%%)\n', ...
        sum(qc_i.adv_logF(valid)),        pct(qc_i.adv_logF(valid)));
    wlog('   Adv B  reliability below median-kMAD  : %d  (%.1f%%)\n', ...
        sum(qc_i.adv_reliability(valid)), pct(qc_i.adv_reliability(valid)));
    wlog('   Adv C  bias index below threshold      : %d  (%.1f%%)\n', ...
        sum(qc_i.adv_bias(valid)),        pct(qc_i.adv_bias(valid)));
    wlog('   Adv D  EC Mahalanobis outlier          : %d  (%.1f%%)\n', ...
        sum(qc_i.adv_ec_outlier(valid)),  pct(qc_i.adv_ec_outlier(valid)));
    wlog('   Adv E  [UPCOMING - cross-network discordance]: %d  (%.1f%%)\n', ...
        sum(qc_i.adv_crossnet(valid)),    pct(qc_i.adv_crossnet(valid)));
    wlog('   Review recommended (>= 2 advisories)  : %d  (%.1f%%)\n', ...
        sum(qc_i.review_recommended(valid)), pct(qc_i.review_recommended(valid)));

    % Advisory co-occurrence (A-D)
    adv_mat = [qc_i.adv_logF(valid), qc_i.adv_reliability(valid), ...
               qc_i.adv_bias(valid), qc_i.adv_ec_outlier(valid)];
    adv_labels = {'A','B','C','D'};
    has_cooccur = false;
    for ai = 1:4
        for aj = (ai+1):4
            n_both = sum(adv_mat(:,ai) & adv_mat(:,aj));
            if n_both > 0
                if ~has_cooccur
                    wlog('\n   Advisory co-occurrence:\n');
                    has_cooccur = true;
                end
                wlog('     Adv %s ∩ Adv %s : %d subject(s)\n', ...
                    adv_labels{ai}, adv_labels{aj}, n_both);
            end
        end
    end

    %% 4. Reliability summary
    wlog('\n');
    wlog(' SUBJECT-LEVEL STABILITY (split-half reliability)\n');
    wlog(' %s\n', sep_minor);
    ig = qc_i.include_group;
    n_ig = sum(ig);

    if n_ig > 0
        s = qc_i.icc_s_summary;
        wlog(' Per-subject ICC_s (Spearman-Brown corrected ICC(3,1), Fisher-z averaged):\n');
        wlog(' %-34s %8s %8s %8s\n',       'Variant                                  ', 'Mean', 'Min', 'Max');
        wlog(' %-34s %8.3f %8.3f %8.3f\n', 'All connections (all)                    ', s.icc_s_all.mean, s.icc_s_all.min, s.icc_s_all.max);
        wlog(' %-34s %8.3f %8.3f %8.3f\n', 'Between-ROI connections only             ', s.icc_s_offdiag.mean, s.icc_s_offdiag.min, s.icc_s_offdiag.max);
        wlog(' %-34s %8.3f %8.3f %8.3f\n', 'Self-connections only                    ', s.icc_s_diag.mean, s.icc_s_diag.min, s.icc_s_diag.max);
        wlog(' %-34s %8.3f %8.3f %8.3f\n', 'All connections (best-determined)        ', s.icc_s_all_masked.mean, s.icc_s_all_masked.min, s.icc_s_all_masked.max);
        wlog(' %-34s %8.3f %8.3f %8.3f\n', 'Between-ROI connections (best-determined)', s.icc_s_offdiag_masked.mean, s.icc_s_offdiag_masked.min, s.icc_s_offdiag_masked.max);
        wlog(' %-34s %8.3f %8.3f %8.3f\n', 'Self-connections (best-determined only)  ', s.icc_s_diag_masked.mean, s.icc_s_diag_masked.min, s.icc_s_diag_masked.max);
        wlog(' %-34s %8.3f %8.3f %8.3f\n', 'Main subject-level reliability score     ', s.tier_icc.mean, s.tier_icc.min, s.tier_icc.max);
        wlog(' %-34s %8.3f %8.3f %8.3f\n', 'Pearson r, between-ROI (best-determined) ', s.pearson_r_masked.mean, s.pearson_r_masked.min, s.pearson_r_masked.max);
        wlog(' %-34s %8.3f %8.3f %8.3f\n', 'Bias Index                               ', s.bias_index.mean, s.bias_index.min, s.bias_index.max);
        wlog('\n');
        tier_icc_mean = s.tier_icc.mean;
        [icc_grade, icc_reco] = cicchetti_grade(tier_icc_mean);
        ...
        wlog('\n');
        wlog(' Main subject-level reliability group mean : %.3f\n', tier_icc_mean);
        wlog(' Cicchetti grade : %s\n', icc_grade);
        wlog(' %s\n', icc_reco);
        wlog(' (Cicchetti benchmarks for test-retest; split-half is typically higher.)\n');
        wlog(' Outliers within include_group (|z|>%.1f) : %d / %d\n', ...
        outlier_z, sum(qc_i.outlier_flag(ig)), n_ig);
    else
        wlog(' No include_group subjects — ICC summary not computed.\n');
    end

    % Fit metrics
    wlog('\n');
    wlog('   Fit metrics (all non-hard-gated subjects, N=%d):\n', n_v);
    cos_v  = qc_i.cosine_sim(valid);
    logF_v = qc_i.logF(valid);
    wlog('   cosine_sim_mean  : mean=%.3f  min=%.3f  max=%.3f\n', ...
        mean(cos_v,'omitnan'), min(cos_v,[],'omitnan'), max(cos_v,[],'omitnan'));
    wlog('   logF             : mean=%.2f  min=%.2f  max=%.2f\n', ...
        mean(logF_v,'omitnan'), min(logF_v,[],'omitnan'), max(logF_v,[],'omitnan'));

    % KL mask note
    wlog('\n');
    wlog(' Best-determined mask: %s ', kmeth);
    wlog('(obs/param ratio target=%.1f) ', k_ratio);
    wlog('median connections selected = %d\n', qc_i.kl_mask_n_median);
    if qc_i.kl_mask_n_median < k_min_conn
        wlog(' NOTE: median mask size below configured minimum (%d) — see log warnings.\n', k_min_conn);
    end

    %% 5. Per-connection ICC
    wlog('\n');
    wlog(' CONNECTION-LEVEL STABILITY (per-connection ICC_c; diagnostic only)\n');
    wlog(' %s\n', sep_minor);

    if ~isnan(qc_i.icc_c_mean_all)
        wlog(' Mean ICC_c (between-ROI, all) : %.3f\n', qc_i.icc_c_mean_all);
        wlog(' Mean ICC_c (between-ROI, best-determined) : %.3f\n', qc_i.icc_c_mean_kl_masked);
        wlog(' Mean ICC_c (between-ROI, outliers excl.) : %.3f\n', qc_i.icc_c_mean_all_no_outlier);

        offdiag = ~logical(eye(N_ROIs));
        icc_c_v = qc_i.icc_c(offdiag);
        icc_ct  = qc_i.icc_c_tier(offdiag);
        n_off_total = sum(offdiag(:));
        n_finite    = sum(isfinite(icc_c_v));
        n_ct1 = sum(icc_ct == 1, 'omitnan');
        n_ct2 = sum(icc_ct == 2, 'omitnan');
        n_ct3 = sum(icc_ct == 3, 'omitnan');

        wlog('\n');
        wlog('   Connection breakdown (off-diagonal, N=%d):\n', n_off_total);
        wlog('   Estimable (>= 3 subjects)        : %d  (%.1f%%)\n', ...
            n_finite, 100*n_finite/n_off_total);
        if n_finite > 0
            wlog('   Tier 1 — Fail  (< %.2f)         : %d  (%.1f%% of estimable)\n', ...
                icc_c_fail, n_ct1, 100*n_ct1/n_finite);
            wlog('   Tier 2 — Marginal               : %d  (%.1f%% of estimable)\n', ...
                n_ct2, 100*n_ct2/n_finite);
            wlog('   Tier 3 — Pass  (>= %.2f)        : %d  (%.1f%% of estimable)\n', ...
                icc_c_pass, n_ct3, 100*n_ct3/n_finite);
        end

        n_masked = sum(qc_i.kl_mask_c(offdiag));
        wlog('\n');
        wlog(' KL mask connections : %d / %d\n', n_masked, n_off_total);
        if n_masked > 0
            n_usable_str = sum(qc_i.icc_c_tier(qc_i.kl_mask_c) >= 2, 'omitnan');
            wlog(' Tier 2 & 3 within mask : %d / %d (%.1f%%)\n', n_usable_str, n_masked, 100*n_usable_str/n_masked);
            n_pass_str = sum(qc_i.icc_c_tier(qc_i.kl_mask_c) == 3, 'omitnan');
            wlog(' Tier 3 within mask : %d / %d (%.1f%%)\n', ...
            n_pass_str, n_masked, 100*n_pass_str/n_masked);
        end
    else
        wlog('   Not computed (fewer than 3 include_group subjects).\n');
    end

    %% 6. Posterior informativeness
    wlog('\n');
    wlog(' HOW DATA-INFORMED THE MODEL WAS\n');
    wlog(' %s\n', sep_minor);
    wlog('   R (posterior/prior variance ratio):\n');
    wlog('     Mean (between-ROI, all)             : %.4g\n', qc_i.R_mean_offdiag);
    wlog('     Mean (between-ROI, best-determined) : %.4g\n', qc_i.R_mean_masked);
    wlog('     Mean (self-connections)             : %.4g\n', qc_i.R_mean_diag);
    wlog('   sqrt(D) (posterior SNR):\n');
    wlog('     Mean (between-ROI, all)             : %.4g\n', qc_i.sqrtD_mean_offdiag);
    wlog('     Mean (between-ROI, best-determined) : %.4g\n', qc_i.sqrtD_mean_masked);
    wlog('     Mean (self-connections)             : %.4g\n', qc_i.sqrtD_mean_diag);
    wlog('   (R ~ 1 → prior-dominated; sqrt(D) >> 1 → data-informed)\n');

    %% 6. Self-connections
    wlog('\n');
    wlog(' SELF-CONNECTION SUMMARY\n');
    wlog(' %s\n', sep_minor);
    if isfinite(qc_i.icc_c_mean_kl_masked) && abs(qc_i.icc_c_mean_kl_masked) >= 0.05
        wlog(' Difference from better-determined between-ROI connections : %+.3f (%.1f%%)\n', ...
            qc_i.self_conn_gain_icc_c, qc_i.self_conn_gain_icc_c_pct);
    else
        wlog(' Difference from better-determined between-ROI connections : %+.3f\n', ...
            qc_i.self_conn_gain_icc_c);
        wlog(' Percent difference not shown because the reference value is near zero.\n');
    end

    %% 7. Asymmetry
    wlog('\n');
    wlog(' DIRECTIONALITY / ASYMMETRY\n');
    wlog(' %s\n', sep_minor);
    wlog(' Questions: Is directionality reliably non-trivial?\n');
    wlog(' (paired t-test per connection, FDR-BH corrected)\n');
    wlog(' Pairs tested (all)                      : %d\n', numel(qc_i.ec_vs_fc.direction_pvals_unmasked));
    wlog(' Significant after FDR (all)             : %.1f%%\n', qc_i.ec_vs_fc.direction_pct_sig_unmasked);
    wlog(' Pairs tested (best-determined)          : %d\n', numel(qc_i.ec_vs_fc.direction_pvals_masked));
    wlog(' Significant after FDR (best-determined) : %.1f%%\n', qc_i.ec_vs_fc.direction_pct_sig_masked);
    wlog(' Answer                                  : %s\n', yn_str(qc_i.ec_vs_fc.directionality_nontrivial));

    %% 8. Per-subject table
    wlog('\n');
    wlog(' PER-SUBJECT TABLE (non-hard-gated)\n');
    wlog(' %s\n', sep_minor);
    wlog(' %-28s %5s %8s %8s %8s %8s %8s %6s Advisories\n', ...
        'Subject', 'Tier', 'tierICC', 'ICC_offM', 'ICC_dgM', 'cosine', 'logF', 'bias');
    wlog(' %s\n', repmat('-', 1, 100));

    tier_labels = {'FAIL','MARG','PASS'};
    for i = 1:N_subs_i
        t = qc_i.subject_tier(i);
        if t == 0, continue; end
        tlabel = tier_labels{t};

        % Advisory string: compact list of fired advisories
        adv_parts = {};
        if qc_i.adv_logF(i),        adv_parts{end+1} = 'A'; end %#ok
        if qc_i.adv_reliability(i), adv_parts{end+1} = 'B'; end %#ok
        if qc_i.adv_bias(i),        adv_parts{end+1} = 'C'; end %#ok
        if qc_i.adv_ec_outlier(i),  adv_parts{end+1} = 'D'; end %#ok
        if qc_i.adv_crossnet(i),    adv_parts{end+1} = 'E'; end %#ok
        if isempty(adv_parts)
            adv_str = '';
        else
            adv_str = ['Adv ' strjoin(adv_parts, '+')];
            if qc_i.review_recommended(i)
                adv_str = [adv_str ' REVIEW']; %#ok
            end
        end

        wlog(' %-28s %5s %8.3f %8.3f %8.3f %8.3f %8.2f %6.3f %s\n', ...
            sub_ids{i}, tlabel, ...
            qc_i.tier_icc(i), qc_i.icc_s_offdiag_masked(i), qc_i.icc_s_diag_masked(i), ...
            qc_i.cosine_sim(i), qc_i.logF(i), ...
            qc_i.bias_index(i), adv_str);
    end
    wlog(' %s\n', repmat('-', 1, 100));

    %% 9. Hard-gated subject list
    n_hg = sum(qc_i.hard_gate_flag);
    if n_hg > 0
        wlog('\n');
        wlog(' HARD-GATED SUBJECTS (Tier 0 — excluded from all analyses)\n');
        wlog(' %s\n', sep_minor);
        wlog(' %-28s %8s  Exclusion reason\n', 'Subject', 'cosine');
        for i = 1:N_subs_i
            if ~qc_i.hard_gate_flag(i), continue; end
            wlog(' %-28s %8.3f  %s\n', sub_ids{i}, ...
                qc_i.cosine_sim(i), plain_reason(qc_i.exclusion_reason{i}));
        end
    end

    %% 10. Feature-set CSV export
    wlog('\n');
    wlog(' FEATURE EXPORT\n');
    wlog(' %s\n', sep_minor);
    local_export_feature_sets(results, params, roi_subsets, i_rs, ...
        rs_label, rs_field, report_dir, sub_ids, @wlog);

    wlog('\n');

end % i_rs

%% --- Appendix: methodology -----------------------------------------

wlog('APPENDIX — HOW THIS WAS DETERMINED\n');
wlog('%s\n', sep_minor);
wlog('\n');
wlog(' DECISION FLOW\n');
wlog(' Each subject is evaluated in four stages, in this order:\n');
wlog('\n');
wlog(' Stage 1 — Hard gates (absolute, applied before ICC)\n');
wlog(' Gate 1a: cosine_sim_mean < 0 → Tier 0 (excluded)\n');
wlog(' Gate 1b: min(ICC_s_diag_masked,ICC_s_offdiag_masked) < %.2f → Tier 0\n', icc_floor);
wlog('\n');
wlog(' Stage 2 — Outlier flag + soft tier (KL-masked ICC_s_all)\n');
wlog(' tier_icc < %.2f → Tier 1 (fail; also Gate 1b)\n', icc_floor);
wlog(' Outlier: |robust z(ICC_s_all_masked)| > %.1f → Tier 2 (marginal)\n', outlier_z);
wlog(' tier_icc < %.2f (not outlier) → Tier 2 (marginal)\n', icc_s_pass);
wlog(' tier_icc >= %.2f AND not outlier → Tier 3 (pass)\n', icc_s_pass);
wlog(' (tier_icc = min(ICC_s_diag_masked, ICC_s_offdiag_masked))\n');
wlog('\n');
wlog(' Stage 3 — Advisory flags (never cause exclusion)\n');
wlog(' Adv A: logF > 2 SD below group mean\n');
wlog(' Adv B: outlier_flag — |robust z(ICC_s_all_masked)| > %.1f (alias of Stage 2)\n', outlier_z);
wlog(' Adv C: bias index (ICC_off - Pearson_r) < %.2f\n', bias_thr);
wlog(' Adv D: EC matrix Mahalanobis outlier (PCA-reduced)\n');
wlog(' Adv E: cross-network tier discordance (set by orchestrator)\n');
wlog(' review_recommended = true when >= 2 advisories fire\n');
wlog('\n');
wlog(' Stage 4 — Posterior informativeness (dataset context, not per-subject)\n');
wlog(' R (posterior/prior variance ratio) and sqrt(D) (posterior SNR)\n');
wlog(' are group-level scalars. Values near 1 for R indicate the prior\n');
wlog(' dominates; high sqrt(D) indicates data-informed connections.\n');
wlog('\n');
wlog(' INCLUSION VECTORS\n');
wlog(' include_all         — Tier 1 + 2 + 3 (all non-hard-gated)\n');
wlog(' include_group       — Tier 2 + 3 (group-level analyses)\n');
wlog(' include_individual  — Tier 3 only (individual-level analyses)\n');
wlog('\n');
wlog(' ICC VARIANTS (all Spearman-Brown corrected ICC(3,1), Fisher-z averaged)\n');
wlog(' Six values per subject: {all, offdiag, diag} x {unmasked, KL-masked}\n');
wlog(' ICC_s_all / _offdiag / _diag — unmasked\n');
wlog(' ICC_s_all_masked / _offdiag_masked / _diag_masked — KL-masked\n');
wlog(' tier_icc = min(ICC_s_diag_masked, ICC_s_offdiag_masked) — TIER CRITERION\n');
wlog(' Pearson_r_masked — Pearson r of off-diagonal split halves, KL-masked\n');
wlog(' Bias_index — ICC_offdiag_masked minus Pearson_r_masked\n');
wlog('     Negative values indicate a systematic half-session\n');
wlog('     scaling difference independent of ordinal agreement.\n');
wlog('\n');
wlog(' KL MASK\n');
if strcmpi(kmeth, 'top_k')
    wlog(' Method: top-k connections ranked by full-run KL(posterior||prior).\n');
else
    wlog(' Method: top quantile (p=%.2f) of connections ranked by full-run KL.\n', k_p);
end
wlog(' Ranking source: full-run KL_A (one stable ranking per subject).\n');
wlog(' Sizing source: split-half obs/param ratio (target=%.1f), via N_TRs_half.\n', k_ratio);
wlog(' One mask built over all connections, then sliced into diag/offdiag.\n');
wlog(' Spearman-Brown: ICC_SB = 2*ICC_raw / (1 + ICC_raw).\n');
wlog('\n');
wlog(' NOTE ON SPLIT-HALF vs TEST-RETEST ICC\n');
wlog(' These are within-session estimates. Expected to be substantially\n');
wlog(' higher than test-retest ICC (which was 0.24–0.45 in Frässle &\n');
wlog(' Stephan 2022). The two quantities are NOT directly comparable.\n');
wlog('\n');
wlog(' PER-CONNECTION ICC (ICC_c) — diagnostic only\n');
wlog(' Computed per off-diagonal connection across include_group subjects.\n');
wlog(' Spearman-Brown corrected, Fisher-z averaged for summary stats.\n');
wlog(' Ranked/masked by group-mean KL divergence (not raw |EC| strength,\n');
wlog(' per lab decision — see connection reliability note in each ROI set).\n');
wlog(' Tier thresholds:\n');
wlog('   Fail:     ICC_c < %.2f\n', icc_c_fail);
wlog('   Marginal: %.2f <= ICC_c < %.2f\n', icc_c_fail, icc_c_pass);
wlog('   Pass:     ICC_c >= %.2f\n', icc_c_pass);
wlog(' ICC_c does not affect subject inclusion.\n');
wlog('%s\n\n', sep_major);

wlog('%s\n', sep_major);
wlog('END OF REPORT\n');
wlog('%s\n', sep_major);

fclose(fid);
rdcm_log(params, 1, 'QC report written to: %s\n', report_path);

end % main function

%% ====================================================================
%  local_print_acq_context
%% ====================================================================

function local_print_acq_context(wlog, params, N_ROIs, N_subs, sep_minor)

TR  = params.rdcm.dt;
N_t = NaN;
if isfield(params, 'data_meta') && isfield(params.data_meta, 'N_timepoints')
    N_t = params.data_meta.N_timepoints;
end

wlog('\n Acquisition context:\n');
wlog('   TR               : %.2f s\n', TR);
if ~isnan(N_t)
    wlog('   N_timepoints     : %d  (%.1f min)\n', N_t, N_t * TR / 60);
else
    wlog('   N_timepoints     : unknown\n');
end
wlog('   N_ROIs           : %d\n', N_ROIs);
wlog('   N_subjects       : %d\n', N_subs);

N_offdiag = N_ROIs * (N_ROIs - 1);
wlog('   Off-diag connections : %d\n', N_offdiag);
wlog('   Self-connections     : %d\n', N_ROIs);

if ~isnan(N_t) && N_offdiag > 0
    ratio = (N_t / 2) / N_offdiag;
    wlog('   Obs-to-param ratio (per split-half) : %.2f\n', ratio);
    if ratio < 1.0
        wlog('   NOTE: < 1.0 — heavily underdetermined. Ridge prior dominates.\n');
        wlog('         Cosine similarity will be very low. ICC is the primary criterion.\n');
    elseif ratio < 2.0
        wlog('   NOTE: < 2.0 — underdetermined. Low cosine expected.\n');
    else
        wlog('   NOTE: >= 2.0 — well-determined. Moderate cosine plausible.\n');
    end
end
wlog(' %s\n', sep_minor);
end

%% ====================================================================
%  local_export_feature_sets
%% ====================================================================

function local_export_feature_sets(results, params, roi_subsets, i_rs, ...
    rs_label, rs_field, export_dir, sub_ids, wlog_fn)
%{
Writes feature-set CSVs for include_group and include_individual subjects.
  features_<set>_offdiag_<group|individual>.csv — all off-diagonal EC
  features_<set>_self_<group|individual>.csv    — diagonal (self-connections)
%}

qc_i     = results.qc{i_rs};
rs       = roi_subsets(i_rs);
roi_names = rs.roi_names;
N_ROIs   = numel(roi_names);
N_subs   = qc_i.N_subs;

sub_id_pattern = params.sub_id_pattern;

inclusion_sets = {'group', 'individual'};
for s = 1:2
    set_name = inclusion_sets{s};
    inc_vec  = qc_i.(sprintf('include_%s', set_name));
    idx      = find(inc_vec);
    N_exp    = numel(idx);

    if N_exp == 0
        wlog_fn('   [export] No %s subjects for %s — skipping.\n', ...
            set_name, rs_label);
        continue;
    end

    % Resolve subject IDs for exported rows
    exp_sub_ids = cell(N_exp, 1);
    for k = 1:N_exp
        i_sub = idx(k);
        ec_k  = results.grid{i_rs, i_sub};
        if isstruct(ec_k) && isfield(ec_k, 'subj') && ~isempty(ec_k.subj)
            tok = regexp(ec_k.subj, sub_id_pattern, 'match');
            if ~isempty(tok), exp_sub_ids{k} = tok{1};
            else,             exp_sub_ids{k} = ec_k.subj; end
        else
            exp_sub_ids{k} = sub_ids{i_sub};
        end
    end

    % --- Off-diagonal export ---
    off_labels = {};
    off_idx    = [];
    for src = 1:N_ROIs
        for tgt = 1:N_ROIs
            if src == tgt, continue; end
            off_labels{end+1} = sprintf('%s->%s', roi_names{src}, roi_names{tgt}); %#ok
            off_idx(end+1)    = (src-1)*N_ROIs + tgt; %#ok
        end
    end
    N_off   = numel(off_labels);
    off_mat = nan(N_exp, N_off);

    for k = 1:N_exp
        A = local_get_A(results.grid{i_rs, idx(k)});
        if ~isempty(A) && size(A,1) == N_ROIs
            A_vec = A(:);
            off_mat(k,:) = A_vec(off_idx)';
        end
    end

    csv_path = fullfile(export_dir, ...
        sprintf('features_%s_offdiag_%s.csv', rs_field, set_name));
    local_write_csv(csv_path, exp_sub_ids, off_labels(:), off_mat);
    wlog_fn('   [export] Off-diagonal (%s, %s): N=%d → %s\n', ...
        rs_label, set_name, N_exp, csv_path);

    % --- Self-connection export ---
    diag_labels = roi_names(:);
    diag_mat    = nan(N_exp, N_ROIs);

    for k = 1:N_exp
        A = local_get_A(results.grid{i_rs, idx(k)});
        if ~isempty(A) && size(A,1) == N_ROIs
            diag_mat(k,:) = diag(A)';
        end
    end

    csv_path = fullfile(export_dir, ...
        sprintf('features_%s_self_%s.csv', rs_field, set_name));
    local_write_csv(csv_path, exp_sub_ids, diag_labels, diag_mat);
    wlog_fn('   [export] Self-connections (%s, %s): N=%d → %s\n', ...
        rs_label, set_name, N_exp, csv_path);
end
end

%% ====================================================================
%  cicchetti_grade
%% ====================================================================

function [grade, recommendation] = cicchetti_grade(icc_val)
if isnan(icc_val)
    grade = 'Unknown';
    recommendation = 'Insufficient data to assess reliability.';
elseif icc_val >= 0.75
    grade = 'Excellent (>= 0.75)';
    recommendation = 'Supports group-level and individual-level analyses.';
elseif icc_val >= 0.60
    grade = 'Good (0.60 – 0.75)';
    recommendation = 'Supports group-level analyses. Individual-level with caution.';
elseif icc_val >= 0.40
    grade = 'Fair (0.40 – 0.60)';
    recommendation = 'Group-level analyses with caution. Avoid individual-level inference.';
else
    grade = 'Poor (< 0.40)';
    recommendation = 'Reliability is low. Null results in group analyses are expected.';
end
end

%% ====================================================================
%  Helpers
%% ====================================================================

function local_write_csv(csv_path, sub_ids, labels, mat)
fid = fopen(csv_path, 'w');
if fid == -1, error('rdcm_report_qc: cannot write %s', csv_path); end
fprintf(fid, 'subject_id');
for c = 1:numel(labels)
    fprintf(fid, ',%s', labels{c});
end
fprintf(fid, '\n');
for k = 1:size(mat, 1)
    fprintf(fid, '%s', sub_ids{k});
    fprintf(fid, ',%.6f', mat(k,:));
    fprintf(fid, '\n');
end
fclose(fid);
end

function A = local_get_A(ec_out)
A = [];
if isempty(ec_out) || ~isstruct(ec_out), return; end
if isfield(ec_out, 'A') && ~isempty(ec_out.A)
    A = ec_out.A;
elseif isfield(ec_out, 'Ep') && isstruct(ec_out.Ep) && isfield(ec_out.Ep, 'A')
    A = ec_out.Ep.A;
end
end

function v = get_param_safe(params, section, field, default)
if isfield(params, section) && isfield(params.(section), field) ...
        && ~isempty(params.(section).(field))
    v = params.(section).(field);
else
    v = default;
end
end

function p = pct(vec)
n = numel(vec);
if n == 0, p = 0; else, p = 100 * sum(vec) / n; end
end

function s = yn_str(x)
if x, s = 'YES'; else, s = 'NO'; end
end

function s = plain_reason(reason_in)
if contains(reason_in, 'marginal tier_icc')
    s = 'split-half reliability < threshold for "best"';
elseif contains(reason_in, 'logF low')
    s = 'model fit low outlier';
elseif contains(reason_in, 'EC outlier')
    s = 'EC matrix outlier';
else
    s = reason_in;
end
end