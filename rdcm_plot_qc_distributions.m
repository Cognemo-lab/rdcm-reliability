function rdcm_plot_qc_distributions(results, params, roi_subsets)
%% Preamble
%{
Generates one multi-panel QC distribution figure per ROI set and saves
each as a PNG to params.dirs.output/export/qc_figures/.

Figure layout (6 panels, 2x3):
1 Histogram: tier_icc [tier criterion — icc_floor + icc_s_pass lines]
2 Histogram: ICC_all_masked [outlier/Advisory B test — median +/- outlier_z*MAD band]
3 Scatter: ICC_diag_masked vs ICC_offdiag_masked, coloured by tier [self-connection reliability vs between-ROI reliability, unity line]
4 Histogram: logF [Advisory A — 2 SD line]
5 Histogram: bias_index [Advisory C — threshold line]
6 Bar chart: advisory flag counts (A-E)

Threshold lines always reflect the absolute decision values actually used
(icc_floor, icc_s_pass, outlier_z, bias_threshold).

Cosine-based and unmasked-ICC panels were dropped: cosine_sim is expected
to be uniformly near-floor for these heavily underdetermined models (see
group report acquisition-context note), and unmasked ICC histograms are
superseded by the KL-masked variants that actually drive tier assignment.

Tier colours: Pass=green, Marginal=orange, Fail=red, Excluded=grey.
Outlier subjects (per qc.outlier_flag) are labelled by subject ID.
---------------------------------------------------------------------------
INPUTS
---------------------------------------------------------------------------
results     := struct from rdcm_pipeline
  .qc       {N_roisets x 1} cell of qc structs from rdcm_qc_group
  .grid     {N_roisets x N_subs} cell — used only for subject ID fallback
params      := pipeline params struct
roi_subsets := struct array from rdcm_define_roisets
%}

out_dir = fullfile(params.dirs.output, 'export', 'qc_figures');
if ~isfolder(out_dir), mkdir(out_dir); end

if isfield(params, 'log') && isfield(params.log, 'run_id')
    run_id = params.log.run_id;
else
    run_id = 'unknown';
end

% Read decision thresholds (must match rdcm_qc_group defaults)
icc_floor  = qp(params, 'icc_floor', 0.40);
icc_s_pass = qp(params, 'icc_s_pass', 0.60);
bias_thr   = qp(params, 'bias_threshold', -0.10);
outlier_z  = qp(params, 'outlier_z', 3.5);

% Tier palette  [Pass=green, Marginal=amber, Fail=red, Hard-gate=grey]
T_COL = [0.18 0.63 0.18;   % Tier 3 — Pass
         0.93 0.55 0.00;   % Tier 2 — Marginal
         0.80 0.12 0.12;   % Tier 1 — Fail
         0.55 0.55 0.55];  % Tier 0 — Excluded

N_roisets = numel(roi_subsets);

for i_rs = 1:N_roisets

    rs_label = roi_subsets(i_rs).label;
    N_ROIs   = numel(roi_subsets(i_rs).roi_names);
    qc_i     = results.qc{i_rs};
    N_subs   = qc_i.N_subs;

    % Subject IDs
    % Subject IDs — fall back to EC storage .subj if qc struct lacks IDs
    if isfield(qc_i, 'subject_ids') && numel(qc_i.subject_ids) == N_subs
        sub_ids = qc_i.subject_ids(:)';
    elseif isfield(results, 'grid') && size(results.grid, 1) >= i_rs
        ec_row = results.grid(i_rs, :);
        sub_ids = cell(1, N_subs);
        for k = 1:N_subs
            if ~isempty(ec_row{k}) && isstruct(ec_row{k}) && isfield(ec_row{k},'subj') && ~isempty(ec_row{k}.subj)
                sub_ids{k} = ec_row{k}.subj;
            else
                sub_ids{k} = sprintf('sub_%03d', k);
            end
        end
        else
            sub_ids = arrayfun(@(k) sprintf('sub_%03d',k), 1:N_subs, ...
                'UniformOutput', false);
    end

    % Core vectors — all [1 x N_subs]
    tier_icc_v = qc_i.tier_icc(:)';
    icc_all_m  = qc_i.icc_s_all_masked(:)';
    icc_off_m  = qc_i.icc_s_offdiag_masked(:)';
    icc_dg_m   = qc_i.icc_s_diag_masked(:)';
    logF_v     = qc_i.logF(:)';
    bias_v     = qc_i.bias_index(:)';
    tiers      = qc_i.subject_tier(:)';
    hg         = qc_i.hard_gate_flag(:)';
    outlier    = qc_i.outlier_flag(:)';

    % Per-subject tier colour matrix [N_subs x 3]
    tc = zeros(N_subs, 3);
    tc(tiers == 3, :) = repmat(T_COL(1,:), sum(tiers==3), 1);
    tc(tiers == 2, :) = repmat(T_COL(2,:), sum(tiers==2), 1);
    tc(tiers == 1, :) = repmat(T_COL(3,:), sum(tiers==1), 1);
    tc(tiers == 0, :) = repmat(T_COL(4,:), sum(tiers==0), 1);

    % Validity mask — outlier and hard-gate flags come directly from qc struct
    valid = ~hg & ~isnan(tier_icc_v);
    outlier = outlier & valid;
    
    % Outlier band on ICC_all_masked, matching Stage 2 logic exactly
    valid_all = ~hg & ~isnan(icc_all_m);
    if sum(valid_all) >= 4
        center_all = median(icc_all_m(valid_all), 'omitnan');
        scale_all  = 1.4826 * mad(icc_all_m(valid_all), 1);
    else
        center_all = NaN;
        scale_all  = NaN;
    end

    %% --- Build figure -----------------------------------------------
    fig = figure('Visible','off','Position',[0 0 1500 750]);
    sgtitle(sprintf('QC Distributions — %s  (%d ROIs, N=%d, run %s)', ...
        strrep(rs_label,'_',' '), N_ROIs, N_subs, run_id), ...
        'FontSize', 13, 'FontWeight', 'bold');

    %% Panel 1 — tier_icc histogram [TIER CRITERION]
    ax_1 = subplot(2,3,1);
    hist_panel(ax_1, tier_icc_v, 'tier\_icc = min(ICC\_diag\_masked, ICC\_offdiag\_masked)', 'Count');
    hold(ax_1,'on');
    xline(ax_1, icc_floor, 'm-', 'LineWidth', 1.8);
    xline(ax_1, icc_s_pass, 'g--', 'LineWidth', 1.5);
    yl = ylim(ax_1);
    text(ax_1, icc_floor, yl(2)*0.96, ...
    sprintf('floor=%.2f',icc_floor), ...
        'FontSize',7,'Color',[0.75 0 0.75],'HorizontalAlignment','center');
    text(ax_1, icc_s_pass, yl(2)*0.96, ...
    sprintf('pass=%.2f',icc_s_pass), ...
        'FontSize',7,'Color',[0 0.55 0],'HorizontalAlignment','center');
    n_gated = sum(hg);
    title(ax_1, sprintf('Tier reliability (N excluded = %d)', n_gated));
    hold(ax_1,'off');
    
    %% Panel 2 — ICC_all_masked histogram [OUTLIER / ADVISORY B TEST]
    ax_2 = subplot(2,3,2);
    hist_panel(ax_2, icc_all_m, 'ICC\_all\_masked [outlier test variable]', 'Count');
    hold(ax_2,'on');
    if ~isnan(center_all)
        lo_band = center_all - outlier_z * scale_all;
        hi_band = center_all + outlier_z * scale_all;
        xline(ax_2, lo_band, 'b:', 'LineWidth', 1.2);
        xline(ax_2, hi_band, 'b:', 'LineWidth', 1.2);
        yl = ylim(ax_2);
        text(ax_2, center_all, yl(2)*0.90, ...
        sprintf('median=%.2f, |z|>%.1f flagged', center_all, outlier_z), ...
            'FontSize',7,'Color',[0 0 0.8],'HorizontalAlignment','center');
    end
    title(ax_2, sprintf('Outlier check (N flagged = %d)', sum(outlier)));
    hold(ax_2,'off');
    
    %% Panel 3 — ICC_diag_masked vs ICC_offdiag_masked scatter
    ax_3 = subplot(2,3,3);
    scatter_panel(ax_3, icc_off_m, icc_dg_m, tc, outlier, sub_ids, ...
        'ICC\_offdiag\_masked', 'ICC\_diag\_masked', ...
        'Self-connection vs between-ROI reliability');
    hold(ax_3,'on');
    xl_3 = xlim(ax_3); yl_3 = ylim(ax_3);
    lo_3 = min(xl_3(1),yl_3(1)); hi_3 = max(xl_3(2),yl_3(2));
    plot(ax_3, [lo_3 hi_3], [lo_3 hi_3], 'k:', 'LineWidth', 0.8);
    hold(ax_3,'off');
    
    %% Panel 4 — logF histogram [Advisory A]
    ax_4 = subplot(2,3,4);
    hist_panel(ax_4, logF_v, 'logF [Advisory A threshold]', 'Count');
    hold(ax_4,'on');
    if any(valid)
        logF_mu = mean(logF_v(valid),'omitnan');
        logF_sd = std( logF_v(valid), 0, 'omitnan');
        adv_a_thr = logF_mu - 2 * logF_sd;
        xline(ax_4, adv_a_thr, 'r--', 'LineWidth', 1.5);
        yl = ylim(ax_4);
        text(ax_4, adv_a_thr, yl(2)*0.88, ...
        sprintf('mu-2SD=%.1f',adv_a_thr), ...
            'FontSize',7,'Color',[0.8 0 0],'HorizontalAlignment','center');
    end
    hold(ax_4,'off');
    
    %% Panel 5 — bias_index histogram [Advisory C]
    ax_5 = subplot(2,3,5);
    hist_panel(ax_5, bias_v, 'Bias index [ICC\_off\_masked - Pearson\_r\_masked]', 'Count');
    hold(ax_5,'on');
    xline(ax_5, bias_thr, 'r--', 'LineWidth', 1.5);
    xline(ax_5, 0, 'k:', 'LineWidth', 1.0);
    yl = ylim(ax_5);
    text(ax_5, bias_thr, yl(2)*0.88, ...
    sprintf('AdvC=%.2f',bias_thr), ...
        'FontSize',7,'Color',[0.8 0 0],'HorizontalAlignment','center');
    hold(ax_5,'off');
    
    %% Panel 6 — Advisory flag bar chart
    ax_6 = subplot(2,3,6);
    ax_i = ax_6;

    %% Shared legend (bottom-right corner via annotation)
    % Placed on Panel I axes to avoid layout fights
    hold(ax_i,'on');
    lg_h(1) = scatter(ax_i, NaN, NaN, 36, T_COL(1,:), 'filled');
    lg_h(2) = scatter(ax_i, NaN, NaN, 36, T_COL(2,:), 'filled');
    lg_h(3) = scatter(ax_i, NaN, NaN, 36, T_COL(3,:), 'filled');
    lg_h(4) = scatter(ax_i, NaN, NaN, 36, T_COL(4,:), 'filled');
    hold(ax_i,'off');
    legend(ax_i, lg_h, {'Pass (T3)','Marginal (T2)','Fail (T1)','Hard-gate (T0)'}, ...
        'Location','northeast','FontSize',7,'Box','off');

    %% Save
    safe_label = regexprep(rs_label, '[^A-Za-z0-9_\-]', '_');
    fname = fullfile(out_dir, sprintf('qc_dist_%s_%s.png', safe_label, run_id));
    exportgraphics(fig, fname, 'Resolution', 150);
    close(fig);
    rdcm_log(params, 2, 'rdcm_plot_qc_distributions: saved %s\n', fname);

end % i_rs
end % main function

%% ====================================================================
%  Local helpers
%% ====================================================================

function hist_panel(ax, vals, xlab, ylab)
% Clean histogram with consistent styling.
v = vals(isfinite(vals));
if isempty(v)
    title(ax, xlab, 'Interpreter','none');
    xlabel(ax, xlab, 'Interpreter','none');
    return
end
histogram(ax, v, 'NumBins', min(30, max(8, floor(numel(v)/3))), ...
    'FaceColor', [0.35 0.60 0.85], 'EdgeColor', 'white', 'FaceAlpha', 0.85);
xlabel(ax, xlab, 'Interpreter','none');
ylabel(ax, ylab);
% Subtitle with mean ± SD
hold(ax,'on');
mu = mean(v,'omitnan'); sd = std(v,0,'omitnan');
xline(ax, mu, 'k-', 'LineWidth', 1.2);
yl = ylim(ax);
text(ax, mu, yl(2)*0.70, sprintf('\\mu=%.3f\n\\sigma=%.3f', mu, sd), ...
    'FontSize', 7, 'HorizontalAlignment','center', 'Color','k');
hold(ax,'off');
grid(ax,'on'); box(ax,'off');
end

function scatter_panel(ax, xv, yv, tc, outlier, sub_ids, xlab, ylab, ttl)
% Scatter: normal points small/transparent, outliers larger/outlined.
hold(ax,'on');
non_out = ~outlier & isfinite(xv) & isfinite(yv);
out_ok  =  outlier & isfinite(xv) & isfinite(yv);
scatter(ax, xv(non_out), yv(non_out), 32, tc(non_out,:), ...
    'filled', 'MarkerFaceAlpha', 0.65);
scatter(ax, xv(out_ok), yv(out_ok), 52, tc(out_ok,:), ...
    'filled', 'MarkerEdgeColor','k','LineWidth',1.2);
% Label outliers
if any(out_ok)
    xl = xlim(ax);
    x_off = diff(xl) * 0.012;
    idx_out = find(out_ok);
    for k = 1:numel(idx_out)
        ii = idx_out(k);
        text(ax, xv(ii)+x_off, yv(ii), sub_ids{ii}, ...
            'FontSize',6,'Interpreter','none','Clipping','on');
    end
end
hold(ax,'off');
xlabel(ax, xlab, 'Interpreter','none');
ylabel(ax, ylab, 'Interpreter','none');
title(ax, ttl,   'Interpreter','none');
grid(ax,'on'); box(ax,'off');
end

function v = qp(params, field, default)
% Safe read of params.qc.<field> with fallback.
if isfield(params,'qc') && isfield(params.qc, field) ...
        && ~isempty(params.qc.(field))
    v = params.qc.(field);
else
    v = default;
end
end
