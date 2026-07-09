function rdcm_report_subject(results, params, roi_subsets, varargin)
%% Preamble
%{
Generates one plain-text subject-level QC report per ROI set, plus an
optional single-subject deep-dive report when a specific subject is
requested.

Two modes:

  MODE A — Group roster (default):
    One file per ROI set:
      subject_roster_<rs_label>_<run_id>.txt
    Columns: subject_id | tier | ICC_s_strong | ICC_s_offdiag | ICC_s_diag |
             Pearson_r | bias_index | cosine_sim | logF |
             advisories | exclusion_reason
    Sorted: hard-gated first, then Tier 1→3, then by ICC_s_strong ascending
    within each tier (most problematic subjects surface at the top).

  MODE B — Single-subject deep-dive:
    Called as rdcm_report_subject(results, params, roi_subsets, sub_id)
    One file per ROI set for the nominated subject:
      subject_report_<sub_id>_<rs_label>_<run_id>.txt
    Sections:
      1. Identification and tier assignment
      2. Decision trace (which gates and thresholds applied, and why)
      3. All QC metrics with group percentiles
      4. Advisory flag detail (which fired, which did not, thresholds)
      5. EC matrix summary (mean |A|, self-connection mean, off-diag range)
      6. Split-half half-by-half comparison (H1 vs H2 logF, cosine)
      7. Cross-ROI-set consistency (tier and ICC across all ROI sets)

---------------------------------------------------------------------------
INPUTS
---------------------------------------------------------------------------
results     := struct from rdcm_pipeline
  .qc       {N_roisets x 1} cell of qc structs
  .grid     {N_roisets x N_subs} cell of EC_out structs
  .fit      {N_roisets x N_subs} cell of fit_out structs
params      := pipeline params struct
roi_subsets := struct array from rdcm_define_roisets
varargin{1} := (optional) char — subject_id for single-subject deep-dive
%}

%% --- Setup -----------------------------------------------------------
report_dir = fullfile(params.dirs.output, 'export');
if ~isfolder(report_dir), mkdir(report_dir); end

if isfield(params,'log') && isfield(params.log,'run_id')
    run_id = params.log.run_id;
elseif isfield(params,'run_id')
    run_id = params.run_id;
else
    run_id = 'unknown';
end

% Thresholds (must match rdcm_qc_group defaults)
icc_floor  = qp(params,'icc_floor', 0.40);
icc_s_pass = qp(params,'icc_s_pass', 0.60);
bias_thr   = qp(params,'bias_threshold', -0.10);
outlier_z  = qp(params,'outlier_z', 3.5);

N_roisets   = numel(roi_subsets);

% Single-subject mode?
single_mode = nargin >= 4 && ~isempty(varargin{1});
if single_mode
    target_id = varargin{1};
else
    target_id = '';
end

%% --- Mode A: group roster per ROI set --------------------------------
if ~single_mode

    for i_rs = 1:N_roisets
        qc_i    = results.qc{i_rs};
        rs      = roi_subsets(i_rs);
        rs_lbl  = rs.label;
        rs_safe = regexprep(rs_lbl,'[^A-Za-z0-9_\-]','_');
        N_subs  = qc_i.N_subs;

        sub_ids = resolve_sub_ids(qc_i, N_subs, results.grid(i_rs, :));

        fname = fullfile(report_dir, ...
            sprintf('subject_roster_%s_%s.txt', rs_safe, run_id));
        fid   = fopen(fname, 'w');
        if fid == -1
            error('rdcm_report_subject: cannot write %s', fname);
        end

        % Header
        sep = repmat('=', 1, 120);
        fprintf(fid, '%s\n', sep);
        fprintf(fid, 'rDCM SUBJECT ROSTER — %s\n', rs_lbl);
        fprintf(fid, 'Run ID : %s\n', run_id);
        fprintf(fid, 'Date   : %s\n', datestr(now,'yyyy-mm-dd HH:MM:SS')); %#ok<TNOW1,DATST>
        fprintf(fid, 'ROI set: %s  (%d ROIs, N=%d subjects)\n', ...
            rs_lbl, numel(rs.roi_names), N_subs);
        fprintf(fid, '%s\n\n', sep);

        fprintf(fid, 'Sort order: hard-gated → Tier 1 → Tier 2 → Tier 3; ');
        fprintf(fid, 'within tier: ICC_s_strong ascending.\n\n');

        % Column header
        col_fmt = '%-28s %5s %8s %8s %8s %8s %7s %8s %8s %-20s %s\n';
        fprintf(fid, col_fmt, 'subject_id','Tier', ...
        'tierICC','ICCoffM','ICCdgM','Pearson_r','bias', ...
        'cosine','logF','Advisories','Exclusion reason');
        fprintf(fid, '(tierICC = min(ICC self-connections, ICC between-ROI), KL-masked. See appendix in group report for definitions.)\n\n');
        
        fprintf(fid, 'Tier key: EXCL=excluded from all analyses, FAIL=not usable, MARG=usable with caution (group only), PASS=usable for individual analysis.\n\n');

        % Sort: sort_key = tier (0 first), then ICC_s_strong asc within tier
        tiers = qc_i.subject_tier(:);
        icc   = qc_i.tier_icc(:);
        % remap tier 0→-1 so it sorts first
        sort_tier = tiers;
        sort_tier(tiers == 0) = -1;
        [~, ord] = sortrows([sort_tier, icc], [1, 2]);

        tier_lbl = {'FAIL','MARG','PASS'};
        for k = 1:N_subs
            i = ord(k);
            t = qc_i.subject_tier(i);
            if t == 0
                t_str = 'EXCL';
            else
                t_str = tier_lbl{t};
            end

            adv_str = build_adv_str(qc_i, i);
            ex_str  = '';
            if t == 0 && isfield(qc_i,'exclusion_reason') && ...
                    numel(qc_i.exclusion_reason) >= i
                ex_str = qc_i.exclusion_reason{i};
            end

            fprintf(fid, col_fmt, sub_ids{i}, t_str, ...
                fmt_f(qc_i.tier_icc(i)), ...
                fmt_f(qc_i.icc_s_offdiag_masked(i)), ...
                fmt_f(qc_i.icc_s_diag_masked(i)), ...
                fmt_f(qc_i.pearson_r_masked(i)), ...
                fmt_f(qc_i.bias_index(i)), ...
                fmt_f(qc_i.cosine_sim(i)), ...
                fmt_f2(qc_i.logF(i)), ...
                adv_str, ex_str);
        end

        fprintf(fid, '%s\n\n', repmat('-',1,120));

        % Summary counts
        fprintf(fid, 'Summary:\n');
        fprintf(fid, '  Hard-gated (T0) : %d\n', qc_i.N_hard_gate);
        fprintf(fid, ' Tier 1 — Fail : %d (tier_ICC < %.2f)\n', ...
            qc_i.N_fail, icc_floor);
        fprintf(fid, ' Tier 2 — Marg : %d (outlier, OR %.2f <= tier_ICC < %.2f)\n', ...
            qc_i.N_marginal, icc_floor, icc_s_pass);
        fprintf(fid, ' Tier 3 — Pass : %d (tier_ICC >= %.2f, not an outlier)\n', ...
            qc_i.N_pass, icc_s_pass);
        fprintf(fid, '  include_group   : %d\n', sum(qc_i.include_group));
        fprintf(fid, '  include_indiv   : %d\n', sum(qc_i.include_individual));

        fclose(fid);
        rdcm_log(params, 1, 'Subject roster written: %s\n', fname);
    end

%% --- Mode B: single-subject deep-dive --------------------------------
else

    for i_rs = 1:N_roisets
        qc_i   = results.qc{i_rs};
        rs     = roi_subsets(i_rs);
        rs_lbl = rs.label;
        rs_safe = regexprep(rs_lbl,'[^A-Za-z0-9_\-]','_');
        N_subs = qc_i.N_subs;

        sub_ids = resolve_sub_ids(qc_i, N_subs, results.grid(i_rs, :));

        % Find subject index
        i_sub = find(strcmp(sub_ids, target_id), 1);
        if isempty(i_sub)
            rdcm_log(params, 1, ...
                '[rdcm_report_subject] Subject "%s" not found in ROI set %s — skipping.\n', ...
                target_id, rs_lbl);
            continue
        end

        safe_id = regexprep(target_id,'[^A-Za-z0-9_\-]','_');
        fname   = fullfile(report_dir, ...
            sprintf('subject_report_%s_%s_%s.txt', safe_id, rs_safe, run_id));
        fid     = fopen(fname,'w');
        if fid == -1
            error('rdcm_report_subject: cannot write %s', fname);
        end

        sep_maj = repmat('=',1,72);
        sep_min = repmat('-',1,72);
        t       = qc_i.subject_tier(i_sub);
        tier_lbl_map = {'EXCLUDED — not usable in any analysis (Tier 0)', ...
            'FAIL — not usable in any analysis (Tier 1)', ...
            'MARGINAL — usable with caution, group only (Tier 2)', ...
            'PASS — usable for individual-level analysis (Tier 3)'};
        t_str   = tier_lbl_map{t+1};

        fprintf(fid,'%s\n',sep_maj);
        fprintf(fid,'rDCM SUBJECT DEEP-DIVE REPORT\n');
        fprintf(fid,'Subject : %s\n', target_id);
        fprintf(fid,'ROI set : %s  (%d ROIs)\n', rs_lbl, numel(rs.roi_names));
        fprintf(fid,'Run ID  : %s\n', run_id);
        fprintf(fid,'Date    : %s\n', datestr(now,'yyyy-mm-dd HH:MM:SS')); %#ok<TNOW1,DATST>
        fprintf(fid,'%s\n\n', sep_maj);

        %% Section 1 — Tier assignment
        fprintf(fid,'SECTION 1: TIER ASSIGNMENT\n');
        fprintf(fid,'%s\n', sep_min);
        fprintf(fid,'  Final tier : %s\n\n', t_str);

        %% Section 2 — Decision trace
        fprintf(fid,'SECTION 2: DECISION TRACE\n');
        fprintf(fid,'%s\n', sep_min);
        fprintf(fid,' This section shows exactly why this subject passed, was flagged, or was excluded.\n\n');

        cos_i = qc_i.cosine_sim(i_sub);
        icc_i = qc_i.tier_icc(i_sub);

        fprintf(fid,'  Gate 1a: cosine_sim_mean < 0\n');
        fprintf(fid,'    Value  : %.4f\n', cos_i);
        fprintf(fid,'    Result : %s\n\n', yn_fire(cos_i < 0));

        fprintf(fid,' Gate 1b: tier_ICC < icc_floor (%.2f)\n', icc_floor);
        fprintf(fid,' tier_ICC = min(ICC self-connections, ICC between-ROI), KL-masked.\n');
        fprintf(fid,' Value : %s\n', fmt_f(icc_i));
        fprintf(fid,' Result : %s\n\n', yn_fire(~isnan(icc_i) && icc_i < icc_floor));

        if t == 0
            ex_str = '';
            if isfield(qc_i,'exclusion_reason') && numel(qc_i.exclusion_reason) >= i_sub
                ex_str = qc_i.exclusion_reason{i_sub};
            end
            fprintf(fid,' --> Subject EXCLUDED (hard-gated). Reason: %s\n\n', ex_str);
        else
            fprintf(fid,' Hard gates: NONE fired.\n\n');
            
            fprintf(fid,' Outlier check: |robust z(ICC_all, KL-masked)| > %.1f\n', outlier_z);
            fprintf(fid,' z-score : %s\n', fmt_f(qc_i.outlier_zscore(i_sub)));
            fprintf(fid,' Result : %s\n\n', yn_fire(qc_i.outlier_flag(i_sub)));
            
            fprintf(fid,' Soft tier:\n');
            fprintf(fid,' tier_ICC = %.4f\n', icc_i);
            fprintf(fid,' icc_floor threshold = %.2f\n', icc_floor);
            fprintf(fid,' icc_s_pass threshold = %.2f\n', icc_s_pass);
            if isnan(icc_i)
                fprintf(fid,' Result: tier_ICC = NaN — treated as Tier 1 (Fail).\n\n');
            elseif qc_i.outlier_flag(i_sub) && icc_i < icc_s_pass
                fprintf(fid,' Result: outlier AND tier_ICC < pass → TIER 2 (Marginal).\n\n');
            elseif icc_i >= icc_s_pass && ~qc_i.outlier_flag(i_sub)
                fprintf(fid,' Result: tier_ICC >= pass AND not an outlier → TIER 3 (Pass).\n\n');
            else
                fprintf(fid,' Result: tier_ICC < pass → TIER 2 (Marginal).\n\n');
            end
        end

        %% Section 3 — QC metrics with group percentiles
        fprintf(fid,'SECTION 3: QC METRICS WITH GROUP PERCENTILES\n');
        fprintf(fid,'%s\n', sep_min);
        fprintf(fid,' Shows how this subject compares to the rest of the sample on each metric.\n');
        fprintf(fid,' (Percentiles computed over non-excluded subjects in this ROI set.)\n\n');
        
        valid   = ~qc_i.hard_gate_flag;
        metrics = { ...
            'tier_ICC (KL-masked)', qc_i.tier_icc, qc_i.tier_icc(i_sub); ...
            'ICC_offdiag_masked', qc_i.icc_s_offdiag_masked, qc_i.icc_s_offdiag_masked(i_sub); ...
            'ICC_diag_masked', qc_i.icc_s_diag_masked, qc_i.icc_s_diag_masked(i_sub); ...
            'ICC_offdiag_unmasked', qc_i.icc_s_offdiag, qc_i.icc_s_offdiag(i_sub); ...
            'ICC_diag_unmasked', qc_i.icc_s_diag, qc_i.icc_s_diag(i_sub); ...
            'Pearson_r_masked', qc_i.pearson_r_masked, qc_i.pearson_r_masked(i_sub); ...
            'bias_index', qc_i.bias_index, qc_i.bias_index(i_sub); ...
            'cosine_sim', qc_i.cosine_sim, qc_i.cosine_sim(i_sub); ...
            'logF', qc_i.logF, qc_i.logF(i_sub); ...
        };
        fprintf(fid,'  %-20s %10s %8s %8s %8s %8s\n', ...
            'Metric','Value','Pctl','Grp-Mean','Grp-Med','Grp-SD');
        fprintf(fid,'  %s\n', repmat('-',1,68));
        for m = 1:size(metrics,1)
            mn   = metrics{m,1};
            vec  = metrics{m,2}(valid);
            val  = metrics{m,3};
            pct  = 100 * mean(vec(isfinite(vec)) < val, 'omitnan');
            gmu  = mean(vec,'omitnan');
            gmed = median(vec,'omitnan');
            gsd  = std(vec,0,'omitnan');
            fprintf(fid,'  %-20s %10.4f %7.1f%% %8.4f %8.4f %8.4f\n', ...
                mn, val, pct, gmu, gmed, gsd);
        end
        fprintf(fid,'\n');

        %% Section 4 — Advisory flag detail
        fprintf(fid,'SECTION 4: ADVISORY FLAGS\n');
        fprintf(fid,'%s\n', sep_min);
        fprintf(fid,' Advisory flags never exclude a subject on their own — they mark cases worth a second look.\n\n');

        valid_v = ~qc_i.hard_gate_flag;
        adv_defs = {
            'A', 'adv_logF',        'logF anomalously low (> 2 SD below group mean)';
            'B', 'adv_reliability', sprintf('Outlier: |robust z(ICC_all, KL-masked)| > %.1f (same test used for marginal-tier assignment)', outlier_z);
            'C', 'adv_bias',        sprintf('bias_index < %.2f', bias_thr);
            'D', 'adv_ec_outlier',  'EC matrix Mahalanobis outlier (PCA-reduced)';
            'E', 'adv_crossnet',    'Cross-network tier discordance';
        };
        for m = 1:size(adv_defs,1)
            lbl   = adv_defs{m,1};
            fld   = adv_defs{m,2};
            desc  = adv_defs{m,3};
            fired = isfield(qc_i,fld) && qc_i.(fld)(i_sub);
            n_grp = 0;
            if isfield(qc_i,fld), n_grp = sum(qc_i.(fld)(valid_v)); end
            fprintf(fid,'  Adv %s: %s\n', lbl, desc);
            if fired
                fprintf(fid,'    Status : FIRED\n');
            else
                fprintf(fid,'    Status : not fired\n');
            end
            fprintf(fid,'    Group  : %d / %d non-gated subjects also flagged\n\n', ...
                n_grp, sum(valid_v));
        end
        rr = isfield(qc_i,'review_recommended') && qc_i.review_recommended(i_sub);
        fprintf(fid,'  review_recommended : %s\n\n', yn_str(rr));

        %% Section 5 — EC matrix summary
        fprintf(fid,'SECTION 5: EC MATRIX SUMMARY\n');
        fprintf(fid,'%s\n', sep_min);
        ec = results.grid{i_rs, i_sub};
        A  = get_A(ec);
        if ~isempty(A) && size(A,1) == numel(rs.roi_names)
            N_r = size(A,1);
            offmask = ~eye(N_r,'logical');
            A_off = A(offmask);
            A_dg  = diag(A);
            fprintf(fid,'  Off-diagonal |A|: mean=%.4f  SD=%.4f  min=%.4f  max=%.4f\n', ...
                mean(abs(A_off)), std(abs(A_off)), min(abs(A_off)), max(abs(A_off)));
            fprintf(fid,' Self-connections: mean=%.4f SD=%.4f min=%.4f max=%.4f\n', ...
                mean(A_dg), std(A_dg), min(A_dg), max(A_dg));
            grp_valid = ~qc_i.hard_gate_flag;
            fprintf(fid,' Group mean self-connection reliability (ICC_diag_masked): %.3f (context; not subject-specific)\n', ...
                mean(qc_i.icc_s_diag_masked(grp_valid),'omitnan'));
            % Top 5 strongest off-diagonal
            [~, si] = sort(abs(A_off), 'descend');
            off_idx = find(offmask);
            fprintf(fid,'  Top 5 strongest off-diagonal connections:\n');
            for kk = 1:min(5,numel(si))
                lin = off_idx(si(kk));
                [r,c] = ind2sub([N_r N_r], lin);
                fprintf(fid,'    %s -> %s  |A|=%.4f  A=%.4f\n', ...
                    rs.roi_names{c}, rs.roi_names{r}, abs(A(r,c)), A(r,c));
            end
        else
            fprintf(fid,'  EC matrix not available for this subject.\n');
        end
        fprintf(fid,'\n');

        %% Section 6 — Split-half half-by-half comparison
        fprintf(fid,'SECTION 6: SPLIT-HALF HALF-BY-HALF\n');
        fprintf(fid,'%s\n', sep_min);
        fit_i = results.fit{i_rs, i_sub};
        if isstruct(fit_i) && isfield(fit_i,'sh_h1')
            h1 = fit_i.sh_h1;  h2 = fit_i.sh_h2;
            fprintf(fid,'  %-20s %12s %12s %12s\n','Metric','H1','H2','H1-H2 diff');
            pf = {'logF','cosine_sim_mean','R2_mean'};
            for kk = 1:numel(pf)
                v1 = gf(h1,pf{kk});  v2 = gf(h2,pf{kk});
                fprintf(fid,'  %-20s %12.4f %12.4f %12.4f\n', pf{kk}, v1, v2, v1-v2);
            end
        else
            fprintf(fid,'  Half-session fit metrics not stored in results.fit.\n');
            fprintf(fid,'  (Set params.splithalf.store_half_fits = true to enable.)\n');
        end
        fprintf(fid,'\n');

        %% Section 7 — Cross-ROI-set consistency
        fprintf(fid,'SECTION 7: CROSS-ROI-SET CONSISTENCY\n');
        fprintf(fid,'%s\n', sep_min);
        fprintf(fid,' Shows whether this subject''s usability verdict changes depending on which ROI set is used.\n\n');

        fprintf(fid,' %-30s %5s %8s %8s %8s\n', ...
            'ROI set','Tier','tierICC','ICCoffM','cosine');
        fprintf(fid,'  %s\n', repmat('-',1,64));
        for j_rs = 1:N_roisets
            qc_j = results.qc{j_rs};
            ids_j = resolve_sub_ids(qc_j, qc_j.N_subs, results.grid(j_rs, :));
            i_j  = find(strcmp(ids_j, target_id), 1);
            if isempty(i_j)
                fprintf(fid,'  %-30s (not present)\n', roi_subsets(j_rs).label);
                continue
            end
            tj = qc_j.subject_tier(i_j);
            tier_lbl_map2 = {'FAIL','MARG','PASS'};
            if tj == 0
                tlj = 'GATE';
            else
                tlj = tier_lbl_map2{tj};
            end
            fprintf(fid,' %-30s %5s %8.3f %8.3f %8.3f\n', ...
                roi_subsets(j_rs).label, tlj, ...
                qc_j.tier_icc(i_j), ...
                qc_j.icc_s_offdiag_masked(i_j), ...
                qc_j.cosine_sim(i_j));
        end
        fprintf(fid,'\n');

        fprintf(fid,'%s\n',sep_maj);
        fprintf(fid,'END OF SUBJECT REPORT\n');
        fprintf(fid,'%s\n',sep_maj);

        fclose(fid);
        rdcm_log(params, 1, 'Subject report written: %s\n', fname);

    end % i_rs (single mode)
end % mode branch

end % main function

%% ====================================================================
%  Local helpers
%% ====================================================================

function ids = resolve_sub_ids(qc_i, N_subs, ec_row)
if isfield(qc_i,'subject_ids') && numel(qc_i.subject_ids) == N_subs
    ids = qc_i.subject_ids(:)';
elseif nargin >= 3 && ~isempty(ec_row)
    ids = cell(1, N_subs);
    for k = 1:N_subs
        if ~isempty(ec_row{k}) && isstruct(ec_row{k}) && isfield(ec_row{k},'subj') && ~isempty(ec_row{k}.subj)
            ids{k} = ec_row{k}.subj;
        else
            ids{k} = sprintf('sub_%03d', k);
        end
    end
else
    ids = arrayfun(@(k) sprintf('sub_%03d',k), 1:N_subs, 'UniformOutput', false);
end
end

function s = build_adv_str(qc_i, i)
parts = {};
flds  = {'adv_logF','adv_reliability','adv_bias','adv_ec_outlier','adv_crossnet'};
lbls  = {'A','B','C','D','E'};
for k = 1:5
    if isfield(qc_i, flds{k}) && qc_i.(flds{k})(i)
        parts{end+1} = lbls{k}; %#ok
    end
end
if isempty(parts)
    s = '';
else
    s = ['Adv ' strjoin(parts,'+')];
    if isfield(qc_i,'review_recommended') && qc_i.review_recommended(i)
        s = [s ' *REVIEW*'];
    end
end
end

function A = get_A(ec_out)
A = [];
if isempty(ec_out) || ~isstruct(ec_out), return; end
if isfield(ec_out,'A') && ~isempty(ec_out.A)
    A = ec_out.A;
elseif isfield(ec_out,'Ep') && isstruct(ec_out.Ep) && isfield(ec_out.Ep,'A')
    A = ec_out.Ep.A;
end
end

function v = gf(s, fld)
v = NaN;
if isstruct(s) && isfield(s, fld) && ~isempty(s.(fld)), v = s.(fld); end
end

function s = fmt_f(x)
if isnan(x), s = '     NaN';
else,         s = sprintf('%8.3f', x); end
end

function s = fmt_f2(x)
if isnan(x), s = '     NaN';
else,         s = sprintf('%8.2f', x); end
end

function s = yn_str(x)
if x, s = 'YES'; else, s = 'NO'; end
end

function s = yn_fire(x)
if x, s = 'FIRED'; else, s = 'not fired'; end
end

function v = qp(params, field, default)
if isfield(params,'qc') && isfield(params.qc, field) ...
        && ~isempty(params.qc.(field))
    v = params.qc.(field);
else
    v = default;
end
end


