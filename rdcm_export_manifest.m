function rdcm_export_manifest(qc_results, roi_set_names, params, out_dir)
%% Preamble
%{
Writes a human-readable pipeline manifest file summarising every decision
made during a single pipeline run.

The manifest is a plain-text file structured as labelled sections so it can
be read at a glance, grepped programmatically, or diff'd between runs.  It
is NOT a substitute for the per-subject CSV (rdcm_export_csv) — it
summarises decisions at the run / group / ROI-set level.

Structure of the output file:
  ── HEADER              run timestamp, MATLAB version, git hash (if available)
  ── PIPELINE PARAMETERS all params fields, recursively formatted
  ── QC DECISION SUMMARY per ROI set: thresholds used, tier counts, gate counts
  ── EC VS FC CHECKLIST  per ROI set: three-question result + ec_preferred
  ── ADVISORY SUMMARY    per ROI set: per-advisory flag counts
  ── POSTERIOR INFO      per ROI set: R and sqrt(D) summary scalars
  ── CROSS-ROI-SET       number of subjects flagged by adv_crossnet (if set)
  ── REPRODUCIBILITY     strength mask fallback decisions, all param values
---------------------------------------------------------------------------
INPUTS
---------------------------------------------------------------------------
qc_results    := {1 x R} cell array of qc structs from rdcm_qc_group.
roi_set_names := {1 x R} cell array of ROI set name strings.
params        := pipeline params struct (full, including params.qc subfields).
out_dir       := char — output directory. Created if it does not exist.
---------------------------------------------------------------------------
OUTPUTS
---------------------------------------------------------------------------
rdcm_manifest_<timestamp>.txt written to out_dir.
Returns the output path as a char (silently if no output argument requested).
---------------------------------------------------------------------------
NOTES
---------------------------------------------------------------------------
  * rdcm_log is called at verbosity level 2 for the final "wrote" message.
  * All numeric values are formatted to 4 significant figures.
  * The manifest format is stable across pipeline versions: new sections are
    always appended; existing section headers are never renamed.  Downstream
    grep / diff scripts can key on section header lines (lines starting with
    '==').
%}

narginchk(4, 4);

if ~exist(out_dir, 'dir'), mkdir(out_dir); end

R = numel(qc_results);
if numel(roi_set_names) ~= R
    error('rdcm_export_manifest: roi_set_names length must match qc_results.');
end

%% --- Build output filename -------------------------------------------

ts       = datestr(now, 'yyyymmdd_HHMMSS'); %#ok<TNOW1,DATST>
out_path = fullfile(out_dir, sprintf('rdcm_manifest_%s.txt', ts));

fid = fopen(out_path, 'w');
if fid < 0
    error('rdcm_export_manifest: cannot open file for writing: %s', out_path);
end

try
    write_manifest(fid, qc_results, roi_set_names, params, ts, R);
catch ME
    fclose(fid);
    rethrow(ME);
end

fclose(fid);
rdcm_log(params, 2, 'rdcm_export_manifest: wrote %s\n', out_path);

if nargout > 0
    varargout{1} = out_path; %#ok<VARARG>
end

end % main function

%% ====================================================================
%  Internal writer
%% ====================================================================

function write_manifest(fid, qc_results, roi_set_names, params, ts, R)

sep_major = repmat('=', 1, 72);
sep_minor = repmat('-', 1, 72);

%% ── HEADER -----------------------------------------------------------

fprintf(fid, '%s\n', sep_major);
fprintf(fid, 'rDCM PIPELINE MANIFEST\n');
fprintf(fid, '%s\n', sep_major);
fprintf(fid, 'Generated     : %s\n', ts);
fprintf(fid, 'MATLAB version: %s\n', version());
fprintf(fid, 'Platform      : %s\n', computer());

% Git commit hash (best-effort, silent on failure)
git_hash = get_git_hash();
fprintf(fid, 'Git hash      : %s\n', git_hash);

% Pipeline script name (if in params)
if isfield(params, 'pipeline_script') && ~isempty(params.pipeline_script)
    fprintf(fid, 'Pipeline script: %s\n', params.pipeline_script);
end

fprintf(fid, 'ROI sets (%d) : %s\n', R, strjoin(roi_set_names, ' | '));
fprintf(fid, '\n');

%% ── PIPELINE PARAMETERS ---------------------------------------------

fprintf(fid, '%s\n', sep_major);
fprintf(fid, 'PIPELINE PARAMETERS\n');
fprintf(fid, '%s\n', sep_major);
print_struct(fid, params, '');
fprintf(fid, '\n');

%% ── QC DECISION SUMMARY (per ROI set) --------------------------------

fprintf(fid, '%s\n', sep_major);
fprintf(fid, 'QC DECISION SUMMARY\n');
fprintf(fid, '%s\n', sep_major);

for r = 1:R
    qc  = qc_results{r};
    rn  = roi_set_names{r};
    pu  = qc.params_used;

    fprintf(fid, '\nROI SET: %s\n', rn);
    fprintf(fid, '%s\n', sep_minor);

    % --- Thresholds used ---
    fprintf(fid, '  Thresholds\n');
    fprintf(fid, '    icc_s_fail (Gate 1b / Tier 1 floor) : %.2f\n', pu.icc_s_fail);
    fprintf(fid, '    icc_s_pass (Tier 3 threshold)       : %.2f\n', pu.icc_s_pass);
    fprintf(fid, '    icc_c_fail                          : %.2f\n', pu.icc_c_fail);
    fprintf(fid, '    icc_c_pass                          : %.2f\n', pu.icc_c_pass);
    fprintf(fid, '    mad_k (Advisory B multiplier)       : %.2f\n', pu.mad_k);
    fprintf(fid, '    bias_threshold (Advisory C)         : %.2f\n', pu.bias_threshold);
    fprintf(fid, '    adv_d_pca_variance                  : %.2f\n', pu.adv_d_pca_variance);
    fprintf(fid, '    adv_d_alpha                         : %.3f\n', pu.adv_d_alpha);
    if isfield(pu, 'fd_threshold')
        fprintf(fid, '    fd_threshold (Gate 1c)              : %.2f mm\n', pu.fd_threshold);
        fprintf(fid, '    fd_prop_threshold (Gate 1c)         : %.2f\n',    pu.fd_prop_threshold);
    end
    fprintf(fid, '    FD data supplied                    : %s\n',  tf_str(pu.has_fd_data));
    fprintf(fid, '\n');

    % --- Strength mask ---
    fprintf(fid, '  Strength mask\n');
    fprintf(fid, '    method                              : %s\n',  pu.strength_mask_method);
    if strcmpi(pu.strength_mask_method, 'top_k')
        fprintf(fid, '    k                                   : %d\n', pu.strength_mask_k);
    else
        fprintf(fid, '    p (quantile)                        : %.2f\n', pu.strength_mask_p);
    end
    fprintf(fid, '    min_connections                     : %d\n',  pu.strength_mask_min_connections);
    fprintf(fid, '    median connections selected         : %d\n',  qc.strength_mask_n_median);
    fprintf(fid, '    fallback to ICC_s_offdiag           : %s\n',  tf_str(~qc.use_strong_icc));
    fprintf(fid, '\n');

    % --- Subject counts ---
    N_subs = qc.N_subs;
    fprintf(fid, '  Subject counts (N = %d)\n', N_subs);
    fprintf(fid, '    Tier 3 — pass           : %d  (%.1f%%)\n', ...
        qc.N_pass,     100*qc.N_pass/N_subs);
    fprintf(fid, '    Tier 2 — marginal       : %d  (%.1f%%)\n', ...
        qc.N_marginal, 100*qc.N_marginal/N_subs);
    fprintf(fid, '    Tier 1 — fail           : %d  (%.1f%%)\n', ...
        qc.N_fail,     100*qc.N_fail/N_subs);
    fprintf(fid, '    Tier 0 — hard-gated     : %d  (%.1f%%)\n', ...
        qc.N_hard_gate,100*qc.N_hard_gate/N_subs);
    fprintf(fid, '    Prop usable (T2+T3)     : %.3f\n', qc.prop_usable);
    fprintf(fid, '\n');

    % --- Hard gate breakdown ---
    fprintf(fid, '  Hard gate breakdown (Gate 1a / 1b / 1c)\n');
    n_gate_cos = sum(cellfun(@(x) contains(x,'cosine'), qc.hard_gate_type));
    n_gate_icc = sum(cellfun(@(x) contains(x,'icc_floor'), qc.hard_gate_type));
    n_gate_mot = sum(cellfun(@(x) contains(x,'motion'), qc.hard_gate_type));
    fprintf(fid, '    Gate 1a (cosine < 0)    : %d subjects\n', n_gate_cos);
    fprintf(fid, '    Gate 1b (ICC < floor)   : %d subjects\n', n_gate_icc);
    fprintf(fid, '    Gate 1c (motion)        : %d subjects\n', n_gate_mot);
    fprintf(fid, '\n');

    % --- ICC summary ---
    valid = ~qc.hard_gate_flag;
    fprintf(fid, '  ICC summary (non-hard-gated subjects, N=%d)\n', sum(valid));
    fprintf(fid, '    ICC_s_strong  mean ± SD : %.4g ± %.4g\n', ...
        mean(qc.icc_s_strong(valid),'omitnan'), std(qc.icc_s_strong(valid),0,'omitnan'));
    fprintf(fid, '    ICC_s_offdiag mean ± SD : %.4g ± %.4g\n', ...
        mean(qc.icc_s_offdiag(valid),'omitnan'), std(qc.icc_s_offdiag(valid),0,'omitnan'));
    fprintf(fid, '    ICC_s_diag    mean ± SD : %.4g ± %.4g\n', ...
        mean(qc.icc_s_diag(valid),'omitnan'),   std(qc.icc_s_diag(valid),0,'omitnan'));
    fprintf(fid, '    Pearson_r     mean ± SD : %.4g ± %.4g\n', ...
        mean(qc.pearson_r(valid),'omitnan'),     std(qc.pearson_r(valid),0,'omitnan'));
    fprintf(fid, '    Bias_index    mean ± SD : %.4g ± %.4g\n', ...
        mean(qc.bias_index(valid),'omitnan'),    std(qc.bias_index(valid),0,'omitnan'));
    fprintf(fid, '\n');

    % --- Per-connection ICC ---
    fprintf(fid, '  Per-connection ICC (N_group = %d subjects)\n', sum(qc.include_group));
    fprintf(fid, '    ICC_c mean (all off-diag)   : %.4g\n', qc.icc_c_mean_all);
    fprintf(fid, '    ICC_c mean (strong-masked)  : %.4g\n', qc.icc_c_mean_strong);
    if ~all(isnan(qc.icc_c(:)))
        all_od  = qc.icc_c(~eye(size(qc.icc_c,1)));
        fin     = all_od(isfinite(all_od));
        if ~isempty(fin)
            n_t3 = sum(fin >= qc.params_used.icc_c_pass);
            n_t2 = sum(fin >= qc.params_used.icc_c_fail & fin < qc.params_used.icc_c_pass);
            n_t1 = sum(fin <  qc.params_used.icc_c_fail);
            fprintf(fid, '    Conn tier 3 (pass)          : %d  (%.1f%%)\n', n_t3, 100*n_t3/numel(fin));
            fprintf(fid, '    Conn tier 2 (marginal)      : %d  (%.1f%%)\n', n_t2, 100*n_t2/numel(fin));
            fprintf(fid, '    Conn tier 1 (fail)          : %d  (%.1f%%)\n', n_t1, 100*n_t1/numel(fin));
        end
    end
    fprintf(fid, '\n');
end

%% ── EC VS FC CHECKLIST -----------------------------------------------

fprintf(fid, '%s\n', sep_major);
fprintf(fid, 'EC VS FC INFORMATIONAL GAIN CHECKLIST\n');
fprintf(fid, '%s\n', sep_major);

for r = 1:R
    qc = qc_results{r};
    ev = qc.ec_vs_fc;
    rn = roi_set_names{r};

    fprintf(fid, '\nROI SET: %s\n', rn);
    fprintf(fid, '%s\n', sep_minor);
    fprintf(fid, '  Q1 — Self-connections informative ?\n');
    fprintf(fid, '    ICC gain (diag - offdiag)         : %.4g\n', ev.self_conn_icc_gain);
    fprintf(fid, '    Between-subject variance (diag)   : %.4g\n', ev.self_conn_between_sub_var);
    fprintf(fid, '    Within-subject variance  (diag)   : %.4g\n', ev.self_conn_within_roi_var);
    fprintf(fid, '    ANSWER                            : %s\n',   yn_str(ev.self_conn_informative));
    fprintf(fid, '\n');
    fprintf(fid, '  Q2 — Directionality reliably non-trivial ?\n');
    fprintf(fid, '    Asymmetry pairs analysed          : %d\n',   numel(ev.delta_ij));
    fprintf(fid, '    Permutation p-value               : %.4g\n', ev.delta_pvalue);
    fprintf(fid, '    ANSWER                            : %s\n',   yn_str(ev.directionality_nontrivial));
    fprintf(fid, '\n');
    fprintf(fid, '  Q3 — EC matrix meaningfully asymmetric ?\n');
    fprintf(fid, '    Mean upper/lower triangle r       : %.4g\n', ev.symmetry_mean_r);
    fprintf(fid, '    ANSWER                            : %s\n',   yn_str(ev.asymmetry_informative));
    fprintf(fid, '\n');
    fprintf(fid, '  YES count (Q1+Q2+Q3)                : %d / 3\n', ev.n_yes);
    fprintf(fid, '  EC PREFERRED over FC ?              : %s\n\n',   yn_str(ev.ec_preferred));
end

%% ── ADVISORY SUMMARY -------------------------------------------------

fprintf(fid, '%s\n', sep_major);
fprintf(fid, 'ADVISORY FLAG SUMMARY\n');
fprintf(fid, '%s\n', sep_major);

for r = 1:R
    qc  = qc_results{r};
    rn  = roi_set_names{r};
    valid = ~qc.hard_gate_flag;
    n_v   = sum(valid);

    fprintf(fid, '\nROI SET: %s  (non-gated N=%d)\n', rn, n_v);
    fprintf(fid, '%s\n', sep_minor);
    fprintf(fid, '  Adv A — logF anomalously low        : %d  (%.1f%%)\n', ...
        sum(qc.adv_logF(valid)),        pct(qc.adv_logF(valid)));
    fprintf(fid, '  Adv B — reliability below median-kMAD: %d  (%.1f%%)\n', ...
        sum(qc.adv_reliability(valid)), pct(qc.adv_reliability(valid)));
    fprintf(fid, '  Adv C — bias index below threshold  : %d  (%.1f%%)\n', ...
        sum(qc.adv_bias(valid)),        pct(qc.adv_bias(valid)));
    fprintf(fid, '  Adv D — EC outlier (Mahalanobis)    : %d  (%.1f%%)\n', ...
        sum(qc.adv_ec_outlier(valid)),  pct(qc.adv_ec_outlier(valid)));
    fprintf(fid, '  Adv E — cross-network discordance   : %d  (%.1f%%)\n', ...
        sum(qc.adv_crossnet(valid)),    pct(qc.adv_crossnet(valid)));
    fprintf(fid, '  Review recommended (>= 2 advisories): %d  (%.1f%%)\n', ...
        sum(qc.review_recommended(valid)), pct(qc.review_recommended(valid)));
    fprintf(fid, '\n');

    % Advisory co-occurrence (A-D only)
    adv_mat = [qc.adv_logF(valid), qc.adv_reliability(valid), ...
               qc.adv_bias(valid), qc.adv_ec_outlier(valid)];
    labels  = {'A','B','C','D'};
    fprintf(fid, '  Advisory co-occurrence (non-gated only):\n');
    for i = 1:4
        for j = (i+1):4
            n_both = sum(adv_mat(:,i) & adv_mat(:,j));
            if n_both > 0
                fprintf(fid, '    %s ∩ %s : %d subjects\n', labels{i}, labels{j}, n_both);
            end
        end
    end
    fprintf(fid, '\n');
end

%% ── POSTERIOR INFORMATIVENESS ----------------------------------------

fprintf(fid, '%s\n', sep_major);
fprintf(fid, 'POSTERIOR INFORMATIVENESS\n');
fprintf(fid, '%s\n', sep_major);
fprintf(fid, '  (R = posterior/prior variance ratio; sqrt(D) = posterior SNR)\n');
fprintf(fid, '  Values near 1 for R indicate prior-dominated connections.\n');
fprintf(fid, '  Values >> 1 for sqrt(D) indicate data-informed connections.\n\n');

for r = 1:R
    qc = qc_results{r};
    rn = roi_set_names{r};
    fprintf(fid, '  ROI SET: %s\n', rn);
    fprintf(fid, '    R_mean (all off-diag)    : %.4g\n', qc.R_mean_all);
    fprintf(fid, '    R_mean (strong-masked)   : %.4g\n', qc.R_mean_strong);
    fprintf(fid, '    sqrt(D)_mean (all)       : %.4g\n', qc.sqrtD_mean_all);
    fprintf(fid, '    sqrt(D)_mean (strong)    : %.4g\n', qc.sqrtD_mean_strong);
    fprintf(fid, '\n');
end

%% ── CROSS-ROI-SET ADVISORY -------------------------------------------

fprintf(fid, '%s\n', sep_major);
fprintf(fid, 'CROSS-ROI-SET ADVISORY (Adv E)\n');
fprintf(fid, '%s\n', sep_major);
fprintf(fid, '\n');
fprintf(fid, '  Adv E is set externally by the pipeline orchestrator\n');
fprintf(fid, '  after all ROI sets are processed. The counts below\n');
fprintf(fid, '  reflect the state at manifest write time.\n\n');

for r = 1:R
    qc  = qc_results{r};
    rn  = roi_set_names{r};
    n_e = sum(qc.adv_crossnet);
    fprintf(fid, '  %s : %d subject(s) flagged by Adv E\n', rn, n_e);
end
fprintf(fid, '\n');

%% ── REPRODUCIBILITY --------------------------------------------------

fprintf(fid, '%s\n', sep_major);
fprintf(fid, 'REPRODUCIBILITY\n');
fprintf(fid, '%s\n', sep_major);
fprintf(fid, '\n');
fprintf(fid, '  Permutation seed for Advisory D / EC-vs-FC asymmetry test: 42\n');
fprintf(fid, '  (fixed in rdcm_qc_group.m; reproducible across platforms)\n\n');

for r = 1:R
    qc  = qc_results{r};
    pu  = qc.params_used;
    rn  = roi_set_names{r};

    fprintf(fid, '  ROI SET: %s\n', rn);
    fprintf(fid, '    strength_mask_method          : %s\n',  pu.strength_mask_method);
    fprintf(fid, '    strength_mask_k               : %d\n',  pu.strength_mask_k);
    fprintf(fid, '    strength_mask_p               : %.2f\n',pu.strength_mask_p);
    fprintf(fid, '    strength_mask_min_connections : %d\n',  pu.strength_mask_min_connections);
    fprintf(fid, '    icc_s_fail                    : %.2f\n',pu.icc_s_fail);
    fprintf(fid, '    icc_s_pass                    : %.2f\n',pu.icc_s_pass);
    fprintf(fid, '    icc_c_fail                    : %.2f\n',pu.icc_c_fail);
    fprintf(fid, '    icc_c_pass                    : %.2f\n',pu.icc_c_pass);
    fprintf(fid, '    mad_k                         : %.2f\n',pu.mad_k);
    fprintf(fid, '    bias_threshold                : %.2f\n',pu.bias_threshold);
    fprintf(fid, '    adv_d_pca_variance            : %.2f\n',pu.adv_d_pca_variance);
    fprintf(fid, '    adv_d_alpha                   : %.3f\n',pu.adv_d_alpha);
    fprintf(fid, '    n_perm_asymmetry              : %d\n',  pu.n_perm_asymmetry);
    fprintf(fid, '    has_fd_data                   : %s\n',  tf_str(pu.has_fd_data));
    fprintf(fid, '    use_strong_icc (fallback?)    : %s\n',  tf_str(qc.use_strong_icc));
    fprintf(fid, '    strength_mask_n_median        : %d\n',  qc.strength_mask_n_median);
    fprintf(fid, '\n');
end

fprintf(fid, '%s\n', sep_major);
fprintf(fid, 'END OF MANIFEST\n');
fprintf(fid, '%s\n', sep_major);

end % write_manifest

%% ====================================================================
%  Local helpers
%% ====================================================================

function print_struct(fid, s, prefix)
% Recursively print all fields of a struct.
fields = fieldnames(s);
for i = 1:numel(fields)
    f   = fields{i};
    val = s.(f);
    key = [prefix f];
    if isstruct(val)
        print_struct(fid, val, [key '.']);
    elseif islogical(val) && isscalar(val)
        fprintf(fid, '  %-45s : %s\n', key, tf_str(val));
    elseif isnumeric(val) && isscalar(val)
        if isnan(val)
            fprintf(fid, '  %-45s : NaN\n', key);
        elseif isinf(val)
            fprintf(fid, '  %-45s : %s\n', key, num2str(val));
        else
            fprintf(fid, '  %-45s : %.4g\n', key, val);
        end
    elseif isnumeric(val) && numel(val) <= 10
        fprintf(fid, '  %-45s : [%s]\n', key, num2str(val(:)', '%.4g '));
    elseif ischar(val) || (isstring(val) && isscalar(val))
        fprintf(fid, '  %-45s : %s\n', key, char(val));
    elseif iscell(val) && numel(val) <= 8
        str_parts = cellfun(@(v) val_to_str(v), val, 'UniformOutput', false);
        fprintf(fid, '  %-45s : {%s}\n', key, strjoin(str_parts, ', '));
    else
        fprintf(fid, '  %-45s : [%s, %s]\n', key, class(val), mat2str(size(val)));
    end
end
end

function s = val_to_str(v)
if ischar(v) || (isstring(v) && isscalar(v))
    s = char(v);
elseif isnumeric(v) && isscalar(v)
    s = sprintf('%.4g', v);
elseif islogical(v) && isscalar(v)
    s = tf_str(v);
else
    s = sprintf('[%s %s]', class(v), mat2str(size(v)));
end
end

function s = tf_str(x)
if x, s = 'true'; else, s = 'false'; end
end

function s = yn_str(x)
if x, s = 'YES'; else, s = 'NO'; end
end

function p = pct(vec)
% Percentage of true values in a logical vector.
n = numel(vec);
if n == 0, p = 0; else, p = 100 * sum(vec) / n; end
end

function hash = get_git_hash()
% Return short git commit hash of the current working directory.
% Returns 'unavailable' on any failure (no git, not a repo, etc.)
hash = 'unavailable';
try
    [status, result] = system('git rev-parse --short HEAD 2>/dev/null');
    if status == 0 && ~isempty(strtrim(result))
        hash = strtrim(result);
    end
catch
end
end

