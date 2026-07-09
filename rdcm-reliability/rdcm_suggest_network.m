function mask_out = rdcm_suggest_network(results, params, roi_subsets)
%% RDCM_SUGGEST_NETWORK  Data-driven retained-connection network suggestion.
%
%  Implements the connection retention framework grounded in Fraessle &
%  Stephan (2022, Netw Neurosci) without using sparse rDCM (which performs
%  substantially worse for resting-state; r = 0.62 vs 0.96 for classical).
%
%  APPROACH (two independent criteria, both must be met)
%  -------------------------------------------------------
%  1. STATISTICAL SIGNIFICANCE  [primary — Fraessle 2022 approach]
%     One-sample t-test per connection: H0: mean_EC_ij = 0, across subjects.
%     Multiple comparison correction: FDR (default) or Bonferroni.
%     This naturally selects strong connections (Fraessle 2022, Fig 4).
%     *** With N < 10 subjects the t-test is severely underpowered.
%         The function still runs but raises a warning and reports
%         uncorrected alongside FDR results.  ***
%
%  2. SPLIT-HALF RELIABILITY  [validation — our extension of Fraessle 2022]
%     Per-connection ICC(3,1) computed from existing split-half data.
%     Requires N_pass >= 4 subjects for any estimate; >= 20 is recommended.
%     Threshold: ICC_ij > params.suggest_network.icc_threshold (default 0.5).
%
%  OUTPUT
%  ------
%  mask_out.(rs_field)
%    .mask_ttest     : NxN binary, t-test significant connections
%    .mask_icc       : NxN binary, ICC > threshold connections
%    .mask_combined  : NxN binary, BOTH criteria met (recommended)
%    .mean_EC        : NxN mean EC across pass subjects
%    .t_stat         : NxN t-statistic matrix
%    .p_fdr          : NxN FDR-corrected p-values
%    .ICC_perconn    : NxN per-connection split-half ICC
%    .n_retained     : scalar, connections kept by combined mask
%    .pct_retained   : scalar, percent of all off-diagonal connections kept
%    .N_subs_used    : number of subjects entering analysis
%
%  PDFs and .mat file are saved to <params.dirs.output>/export/suggest_network/
%
%  USAGE
%    mask_out = rdcm_suggest_network(results, params, roi_subsets);
%
%  This function is a post-hoc analysis tool — it is NOT called
%  automatically by run_rdcm_pipeline.  Call it manually after the
%  pipeline completes, ideally with N >= 20 subjects.
%
%  REFERENCE
%    Fraessle & Stephan (2022). Test-retest reliability of regression
%    dynamic causal modeling. Netw Neurosci, 6(1), 135-160.
%    doi: 10.1162/netn_a_00215

%% ── Defaults & setup ────────────────────────────────────────────────────
if ~isfield(params, 'suggest_network')
    params.suggest_network = struct();
end
sn = params.suggest_network;
if ~isfield(sn, 'icc_threshold'),   sn.icc_threshold   = 0.50; end
if ~isfield(sn, 'alpha'),           sn.alpha            = 0.05; end
if ~isfield(sn, 'correction'),      sn.correction       = 'fdr'; end  % 'fdr' | 'bonferroni' | 'none'
if ~isfield(sn, 'include_marginal'),sn.include_marginal = true;  end
if ~isfield(sn, 'exclude_diagonal'),sn.exclude_diagonal = true;  end

out_dir = fullfile(params.dirs.output, 'export', 'suggest_network');
if ~isfolder(out_dir), mkdir(out_dir); end

run_id    = local_get_run_id(params);
N_roisets = numel(roi_subsets);
SNR_grid  = params.rdcm.SNR_grid;
N_subs    = size(results.grid, 3);
mask_out  = struct();

FIG_NUM = 778;
rdcm_log(params, 1, 'rdcm_suggest_network: analysing %d ROI set(s), %d subjects\n', ...
    N_roisets, N_subs);

if N_subs < 10
    warning(['rdcm_suggest_network: N=%d subjects is too small for stable ' ...
        'per-connection statistics.\n  Results should be treated as ' ...
        'exploratory.  Recommend N >= 20 for a definitive threshold.'], N_subs);
end

%% ── Per-ROI-set analysis ─────────────────────────────────────────────────
for i_rs = 1:N_roisets

    rs       = roi_subsets(i_rs);
    rs_label = rs.label;
    rs_field = matlab.lang.makeValidName(rs_label);
    N_ROIs   = numel(rs.roi_names);
    N_conn   = N_ROIs^2 - N_ROIs;   % off-diagonal count

    if isfield(results.best, rs_field)
        best = results.best.(rs_field);
    else
        best = struct('snr_idx', 1, 'SNR', SNR_grid(1));
    end
    i_snr_best = best.snr_idx;

    qc_best   = results.qc{i_rs, i_snr_best};
    min_tier  = 2 + ~sn.include_marginal;   % 3=pass only, 2=pass+marginal
    pass_mask = qc_best.subject_tier >= min_tier;
    pass_idx  = find(pass_mask);
    N_pass    = numel(pass_idx);

    rdcm_log(params, 1, '  [%s] best SNR=%g, N_pass=%d\n', rs_label, best.SNR, N_pass);

    %% ── Collect EC and split-half matrices ──────────────────────────────
    A_all  = NaN(N_ROIs, N_ROIs, N_pass);
    A1_all = NaN(N_ROIs, N_ROIs, N_pass);
    A2_all = NaN(N_ROIs, N_ROIs, N_pass);

    for k = 1:N_pass
        i_sub = pass_idx(k);
        ec_s  = results.grid{i_rs, i_snr_best, i_sub};
        sh_s  = results.sh{i_rs,   i_snr_best, i_sub};
        A = local_get_A(ec_s);
        if ~isempty(A), A_all(:,:,k) = A; end
        if isstruct(sh_s)
            if isfield(sh_s,'A1') && ~isempty(sh_s.A1), A1_all(:,:,k) = sh_s.A1; end
            if isfield(sh_s,'A2') && ~isempty(sh_s.A2), A2_all(:,:,k) = sh_s.A2; end
        end
    end

    mean_EC = mean(A_all, 3, 'omitnan');

    %% ── Criterion 1: per-connection one-sample t-test ───────────────────
    [t_stat, p_uncorr] = local_ttest_matrix(A_all, N_ROIs, N_pass);

    % Multiple comparison correction
    p_fdr = local_fdr_matrix(p_uncorr, N_ROIs);   % always compute FDR
    p_bonf = min(p_uncorr * N_conn, 1);           % Bonferroni

    switch lower(sn.correction)
        case 'fdr',        p_thresh = p_fdr;
        case 'bonferroni', p_thresh = p_bonf;
        otherwise,         p_thresh = p_uncorr;
    end

    mask_ttest = p_thresh < sn.alpha;
    if sn.exclude_diagonal
        mask_ttest(logical(eye(N_ROIs))) = false;
    end

    %% ── Criterion 2: per-connection split-half ICC ──────────────────────
    if N_pass >= 4
        ICC_perconn = local_icc_matrix(A1_all, A2_all, N_ROIs, N_pass);
    else
        warning('rdcm_suggest_network [%s]: N_pass=%d < 4; ICC per connection skipped.', ...
            rs_label, N_pass);
        ICC_perconn = NaN(N_ROIs, N_ROIs);
    end

    mask_icc = ICC_perconn > sn.icc_threshold;
    if sn.exclude_diagonal
        mask_icc(logical(eye(N_ROIs))) = false;
    end

    %% ── Combined mask ───────────────────────────────────────────────────
    mask_combined = mask_ttest & mask_icc;
    n_ttest    = sum(mask_ttest(:));
    n_icc      = sum(mask_icc(:));
    n_combined = sum(mask_combined(:));

    rdcm_log(params, 1, ...
        '    t-test sig: %d / %d off-diag (%.1f%%)\n    ICC>%.2f: %d (%.1f%%)\n    Combined: %d (%.1f%%)\n', ...
        n_ttest, N_conn, 100*n_ttest/N_conn, ...
        sn.icc_threshold, n_icc, 100*n_icc/N_conn, ...
        n_combined, 100*n_combined/N_conn);

    %% ── Store results ───────────────────────────────────────────────────
    mask_out.(rs_field) = struct( ...
        'mask_ttest',    mask_ttest, ...
        'mask_icc',      mask_icc, ...
        'mask_combined', mask_combined, ...
        'mean_EC',       mean_EC, ...
        't_stat',        t_stat, ...
        'p_uncorr',      p_uncorr, ...
        'p_fdr',         p_fdr, ...
        'p_bonferroni',  p_bonf, ...
        'ICC_perconn',   ICC_perconn, ...
        'n_retained',    n_combined, ...
        'pct_retained',  100 * n_combined / max(N_conn, 1), ...
        'N_subs_used',   N_pass, ...
        'icc_threshold', sn.icc_threshold, ...
        'alpha',         sn.alpha, ...
        'correction',    sn.correction);

    %% ── Visualisation: four-panel figure ────────────────────────────────
    local_plot_suggest_network(FIG_NUM, mask_out.(rs_field), sn, ...
        rs_label, N_ROIs, rs.roi_names, best.SNR, N_pass, out_dir, run_id);

end  % i_rs

%% ── Save masks to .mat ──────────────────────────────────────────────────
mat_path = fullfile(out_dir, sprintf('suggest_network_%s.mat', run_id));
save(mat_path, 'mask_out', 'params');
rdcm_log(params, 1, 'Mask struct saved to: %s\n', mat_path);

if ishandle(FIG_NUM), close(FIG_NUM); end

end  % rdcm_suggest_network


%% ====================================================================
%%  VISUALISATION
%% ====================================================================

function local_plot_suggest_network(FIG_NUM, ms, sn, rs_label, N_ROIs, ...
    roi_names, best_snr, N_pass, out_dir, run_id)
%% Four-panel figure per ROI set
%  Panel A: Mean |EC| matrix (strength)
%  Panel B: Per-connection ICC matrix
%  Panel C: Combined retained network mask
%  Panel D: ICC vs |EC| scatter  (replicates Fraessle 2022 Fig 4 concept)

if ishandle(FIG_NUM), clf(FIG_NUM); end
hfig = figure(FIG_NUM);
FW = 14; FH = 10;
set(hfig, 'Visible', 'off', 'Color', [1 1 1], ...
    'Units',         'inches', 'Position',    [1 1 FW FH], ...
    'PaperUnits',    'inches', 'PaperSize',   [FW FH], ...
    'PaperPosition', [0 0 FW FH], ...
    'NumberTitle',   'off',    'Name',        'rDCM Suggest Network');

cmap_rdbu   = local_rdbu_cmap(256);
cmap_viridis = local_viridis_cmap(256);
cmap_bwr    = local_bwr_cmap(256);

% ── Reorder matrices for display ────────────────────────────────── %
[~, perm_d, blk_d, blk_s_d, blk_e_d] = local_reorder_matrix([], roi_names);
mean_EC_d  = ms.mean_EC(perm_d, perm_d);
icc_d      = ms.ICC_perconn(perm_d, perm_d);
mask_d     = double(ms.mask_combined(perm_d, perm_d));
abs_EC_d   = abs(mean_EC_d);

% ── Panel A: Mean |EC| ────────────────────────────────────────────── %
ha_A = axes('Position', [0.05, 0.47, 0.26, 0.38], 'Parent', hfig);
abs_EC = abs(ms.mean_EC);
clim_A = prctile(abs_EC_d(:), 99);
if clim_A == 0 || isnan(clim_A), clim_A = 1; end
imagesc(abs_EC_d, [0, clim_A]);
colormap(ha_A, cmap_viridis);
colorbar(ha_A, 'Location', 'southoutside', 'FontSize', 7);
axis square
title('Mean |EC| across subjects', 'FontSize', 9, 'FontWeight', 'bold');
local_add_block_labels(ha_A, blk_d, blk_s_d, blk_e_d, N_ROIs, 5);

% ── Panel B: Per-connection ICC ───────────────────────────────────── %
ha_B = axes('Position', [0.38, 0.47, 0.26, 0.38], 'Parent', hfig);
icc_disp = icc_d;
icc_disp(isnan(icc_disp)) = 0;
imagesc(icc_disp, [-1, 1]);
colormap(ha_B, cmap_bwr);
colorbar(ha_B, 'Location', 'southoutside', 'FontSize', 7);
axis square
title(sprintf('Per-connection ICC (N=%d)', N_pass), ...
    'FontSize', 9, 'FontWeight', 'bold');
local_add_block_labels(ha_B, blk_d, blk_s_d, blk_e_d, N_ROIs, 5);

% ── Panel C: Combined mask ────────────────────────────────────────── %
ha_C = axes('Position', [0.71, 0.47, 0.26, 0.38], 'Parent', hfig);
imagesc(mask_d, [0, 1]);
colormap(ha_C, [0.93 0.93 0.93; 0.15 0.50 0.15]);  % grey=dropped, green=kept
axis square
title(sprintf('Retained network  (%d connections, %.1f%%)', ...
    ms.n_retained, ms.pct_retained), 'FontSize', 9, 'FontWeight', 'bold');
local_add_block_labels(ha_C, blk_d, blk_s_d, blk_e_d, N_ROIs, 5);

% ── Panel D: ICC vs |EC| scatter (Fraessle 2022 Fig 4 concept) ──── %
ha_D = axes('Position', [0.07, 0.09, 0.42, 0.29], 'Parent', hfig);
off_diag = ~logical(eye(N_ROIs));
abs_vals = abs_EC_d(off_diag);
icc_vals = icc_d(off_diag);
valid    = ~isnan(icc_vals) & ~isnan(abs_vals);

if sum(valid) > 10
    x = abs_vals(valid);
    y = icc_vals(valid);

    % Bin means (10 bins by |EC| strength — as in Fraessle 2022)
    n_bins   = 10;
    edges    = prctile(x, linspace(0, 100, n_bins+1));
    bin_ctrs = zeros(n_bins, 1);
    bin_iccs = NaN(n_bins, 1);
    bin_sems = NaN(n_bins, 1);
    for b = 1:n_bins
        in_bin = x >= edges(b) & x < edges(b+1);
        if b == n_bins, in_bin = x >= edges(b); end
        bin_ctrs(b) = mean(x(in_bin));
        if sum(in_bin) > 1
            bin_iccs(b) = mean(y(in_bin), 'omitnan');
            bin_sems(b) = std(y(in_bin), 0) / sqrt(sum(~isnan(y(in_bin))));
        end
    end

    scatter(ha_D, x, y, 2, [0.70 0.75 0.82], 'filled', 'MarkerFaceAlpha', 0.3);
    hold(ha_D, 'on');
    errorbar(ha_D, bin_ctrs, bin_iccs, bin_sems, '-o', ...
        'Color', [0.15 0.40 0.75], 'LineWidth', 1.5, 'MarkerSize', 6, ...
        'MarkerFaceColor', [0.15 0.40 0.75], 'CapSize', 4);
    line(ha_D, ha_D.XLim, [sn.icc_threshold sn.icc_threshold], ...
        'Color', [0.75 0.44 0.00], 'LineStyle', '--', 'LineWidth', 1.0);
    line(ha_D, ha_D.XLim, [0.75 0.75], ...
        'Color', [0.12 0.52 0.12], 'LineStyle', '--', 'LineWidth', 1.0);
    hold(ha_D, 'off');
    xlabel(ha_D, '|EC weight|  (a.u.)', 'FontSize', 9);
    ylabel(ha_D, 'Split-half ICC(3,1)', 'FontSize', 9);
    title(ha_D, 'Reliability vs. Connection Strength  (Fraessle 2022)', ...
        'FontSize', 9, 'FontWeight', 'bold');
    set(ha_D, 'FontSize', 8, 'Box', 'off', 'TickDir', 'out');

    % Pearson r label
    if sum(valid) >= 3
        r_xy = corr(x(~isnan(y(~isnan(x)))), y(~isnan(y) & ~isnan(x)), ...
            'Rows', 'complete');
        text(ha_D, 0.98, 0.06, sprintf('r = %.3f', r_xy), ...
            'Units', 'normalized', 'HorizontalAlignment', 'right', ...
            'FontSize', 8, 'Color', [0.3 0.3 0.3]);
    end
else
    axis(ha_D, 'off');
    text(ha_D, 0.5, 0.5, 'Insufficient data for scatter', ...
        'Units', 'normalized', 'HorizontalAlignment', 'center', 'FontSize', 10);
end

% ── Panel E: Summary text ─────────────────────────────────────────── %
ha_E = axes('Position', [0.57, 0.09, 0.38, 0.29], 'Parent', hfig);
axis off
lines = { ...
    sprintf('ROI set:     %s  (%d ROIs)', rs_label, N_ROIs), ...
    sprintf('Best SNR:    %g', best_snr), ...
    sprintf('N subjects:  %d', N_pass), ...
    '', ...
    sprintf('Correction:  %s  (alpha=%.3f)', sn.correction, sn.alpha), ...
    sprintf('ICC thresh:  %.2f', sn.icc_threshold), ...
    '', ...
    sprintf('t-test sig:  %d  (%.1f%%)', sum(ms.mask_ttest(:)), ms.pct_retained), ...
    sprintf('ICC > %.2f:  %d', sn.icc_threshold, sum(ms.mask_icc(:))), ...
    sprintf('RETAINED:    %d  (%.1f%%)', ms.n_retained, ms.pct_retained), ...
};
if N_pass < 10
    lines{end+1} = '';
    lines{end+1} = '*** N < 10: results exploratory ***';
end

for k = 1:numel(lines)
    yp = 1.0 - (k-1)*0.092;
    col = [0 0 0];
    if contains(lines{k}, 'RETAINED'), col = [0.12 0.52 0.12]; end
    if contains(lines{k}, '***'),      col = [0.80 0.10 0.10]; end
    text(ha_E, 0, yp, lines{k}, 'Units', 'normalized', ...
        'FontName', 'Courier', 'FontSize', 8.5, 'Color', col, ...
        'VerticalAlignment', 'top', 'Interpreter', 'none');
end

% ── Figure header ─────────────────────────────────────────────────── %
annotation(hfig, 'textbox', [0.03, 0.92, 0.94, 0.06], ...
    'String', sprintf('rDCM Suggested Network — ROI set: %s  |  Run: %s', rs_label, run_id), ...
    'LineStyle', 'none', 'FontSize', 13, 'FontWeight', 'bold', ...
    'Interpreter', 'none', 'VerticalAlignment', 'middle');

% ── Export ────────────────────────────────────────────────────────── %
pdf_path = fullfile(out_dir, ...
    sprintf('suggest_network_%s_%s.pdf', matlab.lang.makeValidName(rs_label), run_id));
use_eg = ~verLessThan('matlab', '9.8');
if use_eg
    exportgraphics(hfig, pdf_path, 'ContentType', 'vector');
else
    saveas(hfig, pdf_path);
end
rdcm_log_simple(sprintf('  Network suggestion PDF: %s\n', pdf_path));

end  % local_plot_suggest_network


%% ====================================================================
%%  STATISTICS
%% ====================================================================

function [t_stat, p_uncorr] = local_ttest_matrix(A_all, N_ROIs, N_pass)
%% One-sample t-test (H0: mu = 0) for each connection across subjects.
t_stat  = zeros(N_ROIs);
p_uncorr = ones(N_ROIs);
if N_pass < 2, return; end

for ii = 1:N_ROIs
    for jj = 1:N_ROIs
        if ii == jj, continue; end
        vals = squeeze(A_all(ii, jj, :));
        vals = vals(~isnan(vals));
        if numel(vals) < 2, continue; end
        mu  = mean(vals);
        se  = std(vals) / sqrt(numel(vals));
        if se == 0, continue; end
        t   = mu / se;
        df  = numel(vals) - 1;
        p   = 2 * (1 - tcdf(abs(t), df));
        t_stat(ii,jj)   = t;
        p_uncorr(ii,jj) = p;
    end
end
end


function p_fdr = local_fdr_matrix(p_uncorr, N_ROIs)
%% Benjamini-Hochberg FDR correction on the off-diagonal p-values.
off_diag = ~logical(eye(N_ROIs));
p_vec    = p_uncorr(off_diag);

% BH procedure
[p_sorted, sort_idx] = sort(p_vec);
n  = numel(p_sorted);
bh = (1:n)' / n;   % BH critical values (alpha=1 here; scale at threshold stage)

% Cumulative minimum from the right for step-up procedure
adj = flipud(cummin(flipud(p_sorted .* n ./ (1:n)')));
adj = min(adj, 1);

p_adj = ones(n, 1);
p_adj(sort_idx) = adj;

p_fdr = ones(N_ROIs);
p_fdr(off_diag) = p_adj;
end


function ICC_mat = local_icc_matrix(A1_all, A2_all, N_ROIs, N_pass)
%% Per-connection ICC(3,1) from split-half pairs across subjects.
%  Formula for k=2 repeated measures:
%    ICC = 2*cov(y1,y2) / (var(y1) + var(y2) + (mu1-mu2)^2)
%  This is equivalent to Shrout & Fleiss ICC(3,1) for k=2.

y1 = reshape(A1_all, N_ROIs^2, N_pass);
y2 = reshape(A2_all, N_ROIs^2, N_pass);

% Replace NaN with column means for robustness
for col = 1:N_pass
    m1 = nanmean(y1(:,col)); y1(isnan(y1(:,col)),col) = m1;
    m2 = nanmean(y2(:,col)); y2(isnan(y2(:,col)),col) = m2;
end

mu1   = mean(y1, 2);
mu2   = mean(y2, 2);
cov12 = mean((y1 - mu1) .* (y2 - mu2), 2);
v1    = var(y1, 0, 2);
v2    = var(y2, 0, 2);
dmu2  = (mu1 - mu2).^2;

denom = v1 + v2 + dmu2;
icc_vec = zeros(N_ROIs^2, 1);
ok    = denom > 0;
icc_vec(ok)  = 2 * cov12(ok) ./ denom(ok);
icc_vec(~ok) = NaN;

ICC_mat = reshape(icc_vec, N_ROIs, N_ROIs);
end




%% ====================================================================
%%  MATRIX REORDERING + BLOCK LABEL HELPERS
%%
%%  Permutes both rows and columns of an EC matrix so that all ROIs of
%%  the same canonical network form a contiguous block, then draws clean
%%  block dividers and labels on both axes.
%%
%%  Resulting label order (16 blocks for 4S256):
%%    Vis-L  Vis-R  SomMot-L  SomMot-R  DAN-L  DAN-R
%%    SalVA-L  SalVA-R  Limbic-L  Limbic-R  Cont-L  Cont-R
%%    DMN-L  DMN-R  Subcort  Cereb
%%
%%  ROI name formats handled (both detected automatically):
%%    '17Networks_LH_DefaultA_PFCm_1'  (explicit prefix, any atlas copy)
%%    'LH_DefaultA_PFCm_1'             (no prefix — XCP-D 4S256 default)
%%    'Cerebellum_Left'                (cerebellar)
%%    anything else                    → 'Subcort'
%%
%%  To add another atlas: extend local_roi_to_network with a new branch.
%% ====================================================================

function [A_sorted, perm_idx, blk_lbl, blk_s, blk_e] = ...
        local_reorder_matrix(A, roi_names)
%% Permute A so that ROIs are ordered by (network, hemisphere).
%
%  Usage:
%    [A_d, perm, lbl, s, e] = local_reorder_matrix(A_win, rs.roi_names);
%    imagesc(A_d, [-c c]);
%    local_add_block_labels(ha, lbl, s, e, N_ROIs);
%
%  Pass A = [] to compute perm_idx only (avoids copying a large matrix
%  when the same permutation is reused for several matrices).

NET_BASE = {'Vis','SomMot','DAN','SalVA','Limbic','Cont','DMN'};
HEMIS    = {'L','R'};

% Build the desired label sequence
ordered = {};
for k = 1:numel(NET_BASE)
    ordered{end+1} = [NET_BASE{k} '-L'];  %#ok<AGROW>
    ordered{end+1} = [NET_BASE{k} '-R'];  %#ok<AGROW>
end
ordered{end+1} = 'Subcort';
ordered{end+1} = 'Cereb';

% Assign each ROI its combined label  ('DMN-L', 'Subcort', …)
N            = numel(roi_names);
roi_combined = cell(N, 1);
for i = 1:N
    [net, hemi]    = local_roi_to_network(roi_names{i});
    if isempty(hemi)
        roi_combined{i} = net;
    else
        roi_combined{i} = [net '-' hemi];
    end
end

% Build permutation index
perm_idx = zeros(N, 1);
blk_lbl  = {};
blk_s    = [];
blk_e    = [];
ptr      = 0;

for k = 1:numel(ordered)
    lbl = ordered{k};
    idx = find(strcmp(roi_combined, lbl));
    if isempty(idx), continue; end
    n_k = numel(idx);
    blk_s(end+1)            = ptr + 1;      %#ok<AGROW>
    perm_idx(ptr+1:ptr+n_k) = idx;
    ptr                      = ptr + n_k;
    blk_e(end+1)            = ptr;           %#ok<AGROW>
    blk_lbl{end+1}          = lbl;           %#ok<AGROW>
end

% Append any unaccounted ROIs as 'Other'
accounted = false(N, 1);
accounted(perm_idx(1:ptr)) = true;
remaining = find(~accounted);
if ~isempty(remaining)
    n_r = numel(remaining);
    blk_s(end+1)            = ptr + 1;
    perm_idx(ptr+1:ptr+n_r) = remaining;
    ptr                      = ptr + n_r;    %#ok<NASGU>
    blk_e(end+1)            = blk_s(end) + n_r - 1;
    blk_lbl{end+1}          = 'Other';
end

% Apply permutation
if isempty(A)
    A_sorted = [];
else
    A_sorted = A(perm_idx, perm_idx);
end

end  % local_reorder_matrix


function local_add_block_labels(ha, blk_lbl, blk_s, blk_e, N_ROIs, varargin)
%% Draw network block dividers and centred labels on both axes.
%  blk_lbl / blk_s / blk_e come from local_reorder_matrix.
%  Optional varargin{1} = font_size override.

N_blocks = numel(blk_lbl);
blk_ctrs = (blk_s + blk_e) / 2;

ha.XTick = [];
ha.YTick = [];
ha.XLim  = [0.5, N_ROIs + 0.5];
ha.YLim  = [0.5, N_ROIs + 0.5];

if ~isempty(varargin) && isnumeric(varargin{1})
    fs = varargin{1};
elseif N_blocks <= 10, fs = 7;
elseif N_blocks <= 16, fs = 6;
else,                  fs = 5;
end

% Rotate x-axis labels when there are many blocks
x_rot = 0;
if N_blocks > 9, x_rot = 90; end

% Divider lines
hold(ha, 'on');
for b = 1:N_blocks - 1
    bdy = blk_e(b) + 0.5;
    line(ha, [bdy bdy],         [0.5 N_ROIs+0.5], ...
         'Color', [0.55 0.55 0.55], 'LineWidth', 0.5, 'HitTest', 'off');
    line(ha, [0.5 N_ROIs+0.5], [bdy bdy], ...
         'Color', [0.55 0.55 0.55], 'LineWidth', 0.5, 'HitTest', 'off');
end
hold(ha, 'off');

% Labels
for b = 1:N_blocks
    % X-axis (source ROIs)
    text(ha, blk_ctrs(b), N_ROIs + 2, blk_lbl{b}, ...
         'HorizontalAlignment', 'center', 'VerticalAlignment', 'top', ...
         'FontSize', fs, 'Interpreter', 'none', ...
         'Rotation', x_rot, 'Clipping', 'off');
    % Y-axis (target ROIs)
    text(ha, -1, blk_ctrs(b), blk_lbl{b}, ...
         'HorizontalAlignment', 'right', 'VerticalAlignment', 'middle', ...
         'FontSize', fs, 'Interpreter', 'none', 'Clipping', 'off');
end

end  % local_add_block_labels


%% ====================================================================
%%  ATLAS PARSING
%%
%%  local_roi_to_network is the ONLY function that knows about atlas
%%  label formats.  Everything else in this file is atlas-agnostic.
%%
%%  Returns [net, hemi] where:
%%    net  — one of: 'Vis','SomMot','DAN','SalVA','Limbic','Cont','DMN',
%%                   'Subcort', 'Cereb'
%%    hemi — 'L', 'R', or '' (empty for subcortical / cerebellar)
%% ====================================================================

function [net, hemi] = local_roi_to_network(name)
%% Map one ROI label string → (canonical network, hemisphere).
%
%  Handles both XCP-D 4S256 (no prefix) and full Schaefer labels.

% ── 1. Schaefer with explicit '7Networks_' / '17Networks_' prefix ─── %
%   e.g. '17Networks_LH_DefaultA_PFCm_1'  or  '7Networks_L_Vis_1'
tok = regexp(name, '(?:7|17)Networks_([LR])[H]?_([A-Za-z]+)', ...
    'tokens', 'once');
if ~isempty(tok)
    h        = tok{1};
    net_base = regexprep(tok{2}, '[ABC]$', '');
    net      = local_schaefer_canonical(lower(net_base));
    if local_is_cortical(net)
        hemi = h;
    else
        net  = 'Subcort';
        hemi = '';
    end
    return
end

% ── 2. No prefix — XCP-D 4S256 default ───────────────────────────── %
%   'LH_DefaultA_PFCm_1' (cortical)  or  'LH_Thal_DP' (subcortical)
%
%   Both start with [LR]H?_ — the second token distinguishes them:
%   cortical network names are in local_schaefer_canonical's switch;
%   subcortical structure names (Thal, Hipp, Caud, …) are not.
tok = regexp(name, '^([LR])[H]?_([A-Za-z]+)', 'tokens', 'once');
if ~isempty(tok)
    h        = tok{1};
    net_base = regexprep(tok{2}, '[ABC]$', '');
    net      = local_schaefer_canonical(lower(net_base));
    if local_is_cortical(net)
        hemi = h;
    else
        net  = 'Subcort';
        hemi = '';
    end
    return
end

% ── 3. Cerebellar ────────────────────────────────────────────────── %
if ~isempty(regexp(name, '^[Cc]ereb', 'once'))
    net  = 'Cereb';
    hemi = '';
    return
end

% ── 4. Fallback ──────────────────────────────────────────────────── %
net  = 'Subcort';
hemi = '';

end  % local_roi_to_network


function net = local_schaefer_canonical(net_base_lower)
%% Map a lower-case, A/B/C-stripped Schaefer network token to one of
%  the seven canonical labels.  Unknown tokens pass through unchanged
%  (they will fail local_is_cortical → treated as subcortical).

switch net_base_lower
    case {'viscent','visperi','vis'}
        net = 'Vis';
    case {'sommot','sommota','sommotb'}
        net = 'SomMot';
    case {'dorsattn','dorsattna','dorsattnb'}
        net = 'DAN';
    case {'salventattn','salventattna','salventattnb', ...
          'salventatt','ventattn','ventatt'}
        net = 'SalVA';
    case {'limbic','limbica','limbicb'}
        net = 'Limbic';
    case {'cont','conta','contb','contc'}
        net = 'Cont';
    case {'default','defaulta','defaultb','defaultc', ...
          'temppar','temppara','tempparb'}
        net = 'DMN';
    otherwise
        net = net_base_lower;   % unknown → will not pass is_cortical
end

end  % local_schaefer_canonical


function tf = local_is_cortical(net)
%% True iff net is one of the seven canonical cortical network names.
tf = any(strcmp(net, {'Vis','SomMot','DAN','SalVA','Limbic','Cont','DMN'}));
end



%% ====================================================================
%%  COLORMAPS
%% ====================================================================

function cmap = local_rdbu_cmap(n)
half = ceil(n/2);
blue = [linspace(0.017,1,half)', linspace(0.443,1,half)', linspace(0.690,1,half)'];
red  = [linspace(1,0.698,n-half)', linspace(1,0.094,n-half)', linspace(1,0.168,n-half)'];
cmap = [blue; red];
end

function cmap = local_bwr_cmap(n)
half = ceil(n/2);
blue = [linspace(0,1,half)', linspace(0,1,half)', ones(half,1)];
red  = [ones(n-half,1), linspace(1,0,n-half)', linspace(1,0,n-half)'];
cmap = [blue; red];
end

function cmap = local_viridis_cmap(n)
% Approximation of viridis using control points
ctrl = [0.267 0.005 0.329; 0.283 0.141 0.458; 0.254 0.265 0.530; ...
        0.207 0.372 0.553; 0.164 0.471 0.558; 0.128 0.567 0.551; ...
        0.135 0.659 0.518; 0.267 0.749 0.441; 0.478 0.821 0.318; ...
        0.741 0.873 0.150; 0.993 0.906 0.144];
t  = linspace(0, 1, size(ctrl,1));
tt = linspace(0, 1, n);
cmap = [interp1(t, ctrl(:,1), tt)', interp1(t, ctrl(:,2), tt)', interp1(t, ctrl(:,3), tt)'];
end


%% ====================================================================
%%  MISC HELPERS
%% ====================================================================

function run_id = local_get_run_id(params)
if isfield(params,'log') && isfield(params.log,'run_id') && ~isempty(params.log.run_id)
    run_id = params.log.run_id;
elseif isfield(params,'run_id') && ~isempty(params.run_id)
    run_id = params.run_id;
else
    run_id = 'unknown';
end
end

function A = local_get_A(ec_out)
A = [];
if isempty(ec_out) || ~isstruct(ec_out), return; end
if isfield(ec_out,'A') && ~isempty(ec_out.A)
    A = ec_out.A;
elseif isfield(ec_out,'Ep') && isstruct(ec_out.Ep) && isfield(ec_out.Ep,'A')
    A = ec_out.Ep.A;
end
end

function rdcm_log_simple(msg)
% Minimal logger used when params may not be available in sub-function
fprintf('%s', msg);
end
