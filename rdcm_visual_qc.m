function rdcm_visual_qc(results, params, roi_subsets)
%% RDCM_VISUAL_QC  Per-subject visual QC PDFs for the rDCM pipeline.
%
%  One PDF per subject; one page per ROI set.  Each page shows:
% - Full-session EC matrix (RdBu diverging colormap, 99th-pct colour clip)
%   with network block labels on both axes instead of per-ROI ticks.
% - Fit & reliability metrics (single config — no SNR sweep).
% - Split-half EC thumbnail pair (A1 vs A2).
% - Self-connection panel (deviation from mean, network-ordered).
%
%  COMPATIBILITY
%    exportgraphics (R2020a+) used for multi-page PDF append.
%    Falls back to per-page PDF files on older MATLAB.

%% Setup
run_id  = local_get_run_id(params);
qc_dir  = fullfile(params.dirs.output, 'export', 'qc_figures');
if ~isfolder(qc_dir), mkdir(qc_dir); end

N_roisets = numel(roi_subsets);
N_subs    = size(results.grid, 2);

%% Load network labels and full ROI name list
net_labels_all   = local_load_network_labels(params); % net_labels_all is either an {N_ROI_total x 1} cell or {} (empty)
all_roi_names_full = local_load_roi_names(params);

rdcm_log(params, 1, 'Visual QC: generating %d subject PDF(s) -> %s\n', N_subs, qc_dir);

set(groot, 'DefaultFigureVisible', 'off');
% Create one hidden figure object for the entire run — never call figure(N)
hfig = figure('Visible', 'off', 'Color', [1 1 1], ...
    'Units', 'inches', 'Position', [1 1 13 9], ...
    'PaperUnits', 'inches', 'PaperSize', [13 9], ...
    'PaperPosition', [0 0 13 9], ...
    'NumberTitle', 'off', 'Name', 'rDCM Visual QC');

for i_sub = 1:N_subs

    sub_id = local_get_subid(results, i_sub, N_roisets);
    pdf_path   = fullfile(qc_dir, sprintf('qc_%s_%s.pdf', sub_id, run_id));
    first_page = true;

    for i_rs = 1:N_roisets

        rs       = roi_subsets(i_rs);
        rs_label = rs.label;
        rs_field = matlab.lang.makeValidName(rs_label);
        N_ROIs   = numel(rs.roi_names);

        %% Subset network labels to this ROI set, then sort
        if ~isempty(net_labels_all)
            [found, roi_idx] = ismember(rs.roi_names, all_roi_names_full);
            if ~all(found)
                warning('rdcm_visual_qc: %d ROI(s) in subset "%s" not found in network labels — skipping block labels.', ...
                    sum(~found), rs_label);
                net_labels_rs = {};
            else
                net_labels_rs = net_labels_all(roi_idx);
            end
        else
            net_labels_rs = {};
        end
        [perm_idx, blk_lbl, blk_s, blk_e] = rdcm_network_sort(rs.roi_names, net_labels_rs);

        %% Best config
        
        qc_i = results.qc{i_rs};

        tier     = local_safe_idx(qc_i.subject_tier, i_sub, 0);
        tier_str = local_tier_str(tier);
        tier_col = local_tier_color(tier);

        ec_out = local_safe_cell(results.grid, i_rs, i_sub);
        sh_out = local_safe_cell(results.sh,   i_rs, i_sub);

        %% Figure setup — reuse the single hidden figure, never raise it
        clf(hfig);
        set(hfig, 'Visible', 'off');   % clf can reset Visible on some MATLAB versions

        %% Left: Winning EC matrix
        ha_ec = axes('Position', [0.12, 0.14, 0.43, 0.71], 'Parent', hfig);
        A_win = local_get_A(ec_out);
        if ~isempty(A_win)

            A_win_d = A_win(perm_idx, perm_idx);
            
            diag_vals = diag(A_win_d);          % already in permuted/network order
            A_plot = A_win_d;
            A_plot(1:N_ROIs+1:end) = NaN;       % mask diagonal
            
            offdiag_vals = A_plot(~isnan(A_plot));
            clim_val = prctile(abs(offdiag_vals), 99);
            if isempty(clim_val) || clim_val == 0 || isnan(clim_val)
                clim_val = 1;
            end
            
            himg = imagesc(ha_ec, A_plot, [-clim_val, clim_val]);
            colormap(ha_ec, local_rdbu_cmap(256));
            
            % Draw diagonal as black squares on top
            hold(ha_ec, 'on');
            for d = 1:N_ROIs
                rectangle(ha_ec, 'Position', [d-0.5, d-0.5, 1, 1], ...
                    'FaceColor', [0 0 0], 'EdgeColor', 'none');
            end
            hold(ha_ec, 'off');

            cb = colorbar(ha_ec, 'Location', 'eastoutside');
            % cb.Position(3) = 0.018;   % narrow bar
            cb.Label.String = 'EC weight (a.u.)';
            cb.FontSize      = 8;
            axis square
            xlabel('Source ROI', 'FontSize', 9);
            ylabel('Target ROI', 'FontSize', 9);
            % Network block labels replace per-ROI ticks
            if ~isempty(blk_lbl)
                local_add_block_labels(ha_ec, blk_lbl, blk_s, blk_e, N_ROIs);
                ha_ec.XLabel.Units = 'normalized';
                ha_ec.XLabel.Position = [0.5, -0.1, 0];
                ha_ec.YLabel.Units = 'normalized';
                ha_ec.YLabel.Position = [-0.1, 0.5, 0];
            else
                % No network labels — suppress ticks but keep axis titles in default position
                ha_ec.XTick = [];
                ha_ec.YTick = [];
            end
        else
            axis off
            text(0.5, 0.5, 'EC data unavailable', 'Units', 'normalized', ...
                'HorizontalAlignment', 'center', 'FontSize', 11, ...
                'Color', [0.5 0.5 0.5]);
        end
        title('Full-session EC', ...
            'FontSize', 10, 'FontWeight', 'bold');

        %% Top-right: Metrics table
        ha_tbl = axes('Position', [0.60, 0.66, 0.38, 0.22], 'Parent', hfig);
        axis off
        title('Fit & Reliability Metrics', 'FontSize', 9, 'FontWeight', 'bold');
        
        fit_i = local_safe_cell(results.fit, i_rs, i_sub);

        cos_val = local_get_scalar(fit_i, 'cosine_sim_mean');
        tier_icc_val = local_safe_idx(qc_i.tier_icc, i_sub, NaN);
        icc_off = local_safe_idx(qc_i.icc_s_offdiag_masked, i_sub, NaN);
        icc_diag = local_safe_idx(qc_i.icc_s_diag_masked, i_sub, NaN);
        pear_off = local_safe_idx(qc_i.pearson_r_masked, i_sub, NaN);
        
        hdr_str = sprintf('%-9s %-9s %-9s %-9s %-9s', ...
            'cosine', 'tierICC', 'ICCoffM', 'ICCdgM', 'PearsonR');
        text(0, 0.90, hdr_str, 'Units', 'normalized', ...
            'FontName', 'Courier', 'FontSize', 8, 'FontWeight', 'bold', ...
            'Interpreter', 'none', 'VerticalAlignment', 'top');
        text(0, 0.80, repmat('-', 1, 52), 'Units', 'normalized', ...
            'FontName', 'Courier', 'FontSize', 8, ...
            'Interpreter', 'none', 'VerticalAlignment', 'top');
        
        row_str = sprintf('%-9.3f %-9.3f %-9.3f %-9.3f %-9.3f', ...
        cos_val, tier_icc_val, icc_off, icc_diag, pear_off);
        text(0, 0.70, row_str, 'Units', 'normalized', ...
            'FontName', 'Courier', 'FontSize', 8, ...
            'Interpreter', 'none', 'VerticalAlignment', 'top');
        text(0, 0.55, '(ICCoffM/ICCdgM/PearsonR are KL-masked, i.e.', ...
            'Units', 'normalized', 'FontName', 'Courier', 'FontSize', 6.5, ...
            'Interpreter', 'none', 'VerticalAlignment', 'top', 'Color', [0.4 0.4 0.4]);
        text(0, 0.48, ' better-determined connections only.)', ...
            'Units', 'normalized', 'FontName', 'Courier', 'FontSize', 6.5, ...
            'Interpreter', 'none', 'VerticalAlignment', 'top', 'Color', [0.4 0.4 0.4]);

        %% Middle-right: Split-half thumbnails
        ha_h1 = axes('Position', [0.60, 0.36, 0.185, 0.26], 'Parent', hfig);
        ha_h2 = axes('Position', [0.80, 0.36, 0.185, 0.26], 'Parent', hfig);
        A1 = []; A2 = [];
        if isstruct(sh_out)
            if isfield(sh_out,'A1') && ~isempty(sh_out.A1), A1 = sh_out.A1; end
            if isfield(sh_out,'A2') && ~isempty(sh_out.A2), A2 = sh_out.A2; end
        end
        if ~isempty(A1) && ~isempty(A2)
            A1_d = A1(perm_idx, perm_idx);
            A2_d = A2(perm_idx, perm_idx);
            A1_plot = A1_d; A1_plot(1:N_ROIs+1:end) = NaN;
            A2_plot = A2_d; A2_plot(1:N_ROIs+1:end) = NaN;
            combined = [A1_plot(~isnan(A1_plot)); A2_plot(~isnan(A2_plot))];
            clim_sh = prctile(abs(combined), 99);
            if clim_sh == 0 || isnan(clim_sh), clim_sh = 1; end
            cmap_sh = local_rdbu_cmap(256);
            % Half-1
            set(hfig, 'CurrentAxes', ha_h1);
            imagesc(A1_plot, [-clim_sh, clim_sh]);
            colormap(ha_h1, cmap_sh);
            hold(ha_h1, 'on');
            for d = 1:N_ROIs
                rectangle(ha_h1, 'Position', [d-0.5, d-0.5, 1, 1], ...
                'FaceColor', [0 0 0], 'EdgeColor', 'none');
            end
            hold(ha_h1, 'off');
            axis square off
            title(ha_h1, 'Half 1', 'FontSize', 8, 'FontWeight', 'bold');
            
            % Half-2
            set(hfig, 'CurrentAxes', ha_h2);
            imagesc(A2_plot, [-clim_sh, clim_sh]);
            colormap(ha_h2, cmap_sh);
            hold(ha_h2, 'on');
            for d = 1:N_ROIs
                rectangle(ha_h2, 'Position', [d-0.5, d-0.5, 1, 1], ...
                    'FaceColor', [0 0 0], 'EdgeColor', 'none');
            end
            hold(ha_h2, 'off');
            axis square off
            title(ha_h2, 'Half 2', 'FontSize', 8, 'FontWeight', 'bold');
        else
            set(hfig, 'CurrentAxes', ha_h1); axis off
            title(ha_h1, 'Half 1', 'FontSize', 8, 'FontWeight', 'bold');
            text(0.5, 0.5, 'Split-half N/A', 'Units', 'normalized', ...
            'HorizontalAlignment', 'center', 'FontSize', 8, ...
            'Color', [0.55 0.55 0.55]);
            set(hfig, 'CurrentAxes', ha_h2); axis off
            title(ha_h2, 'Half 2', 'FontSize', 8, 'FontWeight', 'bold');
        end

        %% Bottom-right: Self-connection panel (keep matrix/network order)
        ha_sc = axes('Position', [0.60, 0.14, 0.38, 0.18], 'Parent', hfig);

        % Reference is the model PRIOR for self-connections, not the empirical
        % mean — self-connections are weakly informed by data and cluster near
        % their prior, so "deviation from mean" mislabels what is really
        % "deviation from prior" and lets a handful of poorly-determined ROIs
        % stretch the axis for everyone else.
        prior_diag_vals = [];
        if isstruct(fit_i) && isfield(fit_i, 'mu0_A') && ~isempty(fit_i.mu0_A)
            mu0_d = fit_i.mu0_A(perm_idx, perm_idx);
            prior_diag_vals = diag(mu0_d);
        end
        if ~isempty(prior_diag_vals) && all(isfinite(prior_diag_vals)) && ...
            range(prior_diag_vals) < 1e-6
            prior_ref = prior_diag_vals(1); % single scalar prior, as expected
        else
            prior_ref = -0.5; % fallback: TAPAS rDCM default self-connection prior
        end

        if ~isempty(diag_vals) && numel(diag_vals) > 1
            sc_centred = diag_vals - prior_ref; % keep permuted order
            y = 1:N_ROIs;

            bh = barh(ha_sc, y, sc_centred, 0.95, 'FaceColor', 'flat', 'EdgeColor', 'none');

            cdata = zeros(N_ROIs, 3);
            cdata(sc_centred >= 0, :) = repmat([0.78 0.18 0.18], sum(sc_centred >= 0), 1);
            cdata(sc_centred < 0, :) = repmat([0.20 0.44 0.78], sum(sc_centred < 0), 1);
            bh.CData = cdata;

            xline(ha_sc, 0, 'k-', 'LineWidth', 0.8);

            % Explicit data-driven XLim with a small margin. "axis tight" was
            % unreliable on this panel (mixed barh + xline objects), and produced
            % visual clipping of the largest-magnitude self-connections. This
            % computes limits directly from the plotted values so nothing is cut off.
            sc_min = min(sc_centred);
            sc_max = max(sc_centred);
            if sc_min == sc_max
                sc_pad = max(abs(sc_min), 0.1) * 0.2;
                if sc_min > 0
                    xl_lo = 0;
                else
                    xl_lo = sc_min - sc_pad;
                end
                xl_hi = sc_max + sc_pad;
            else
                sc_range = sc_max - sc_min;
                if sc_min > 0
                    xl_lo = 0;
                else
                    xl_lo = sc_min - 0.05 * sc_range;
                end
                xl_hi = sc_max + 0.05 * sc_range;
            end
            xl_lo = min(xl_lo, 0); % always include the zero reference line
            xl_hi = max(xl_hi, 0);

            set(ha_sc, 'YDir', 'reverse', ... % top-to-bottom matches matrix
                'YLim', [0.5 N_ROIs + 0.5], ...
                'XLim', [xl_lo, xl_hi], ...
                'FontSize', 7, ...
                'TickDir', 'out', ...
                'Box', 'off');

            % Note: this is posterior mean minus prior mean in native EC units.
            % Self-connections are the most data-informed parameters in the model
            % (see sqrt(D) in the QC report), so unlike off-diagonal connections they
            % are NOT constrained to stay close to the prior — deviations larger
            % than 0.5 are expected and correct when a subject's self-connections
            % are strongly determined by their data.
            xlabel(ha_sc, 'Self-connection: posterior mean minus prior mean (EC units)', ...
                'FontSize', 7);
            title(ha_sc, sprintf('Self-connections (prior mean = %.2f; not bounded)', prior_ref), ...
                'FontSize', 8);

        if ~isempty(blk_lbl)
            local_add_sc_group_labels(ha_sc, blk_lbl, blk_s, blk_e, N_ROIs);
        else
            ha_sc.YTick = [];
        end
        else
            axis(ha_sc, 'off');
            text(0.5, 0.5, 'Self-connection data N/A', 'Units', 'normalized', ...
            'HorizontalAlignment', 'center', 'FontSize', 8, ...
            'Color', [0.55 0.55 0.55]);
        end

        % %% Bottom-right: Self-connection panel — RAW VALUES (expected sign: negative)
        % ha_sc = axes('Position', [0.60, 0.14, 0.38, 0.18], 'Parent', hfig);
        % 
        % if ~isempty(diag_vals) && numel(diag_vals) > 1
        %     y = 1:N_ROIs;
        % 
        %     % Self-connections represent self-inhibition and are expected to be
        %     % negative. Any value >= 0 is flagged as a sign violation — this is a
        %     % more direct QC signal than "deviation from prior," since it checks
        %     % the constraint that actually matters (correct sign), not distance
        %     % from a soft prior that self-connections are free to move away from.
        %     sign_violation = diag_vals >= 0;
        % 
        %     bh = barh(ha_sc, y, diag_vals, 0.95, 'FaceColor', 'flat', 'EdgeColor', 'none');
        % 
        %     cdata = zeros(N_ROIs, 3);
        %     cdata(~sign_violation, :) = repmat([0.20 0.44 0.78], sum(~sign_violation), 1); % normal (negative)
        %     cdata(sign_violation, :)  = repmat([0.85 0.10 0.10], sum(sign_violation), 1);  % flagged (>= 0)
        %     bh.CData = cdata;
        % 
        %     hold(ha_sc, 'on');
        %     xline(ha_sc, 0, 'k-', 'LineWidth', 1.0);
        %     if isfinite(prior_ref)
        %         xline(ha_sc, prior_ref, 'k--', 'LineWidth', 0.6, 'Color', [0.5 0.5 0.5]);
        %     end
        %     hold(ha_sc, 'off');
        % 
        %     % Explicit data-driven XLim so nothing is clipped
        %     sc_min = min(diag_vals);
        %     sc_max = max(diag_vals);
        %     if sc_min == sc_max
        %         sc_pad = max(abs(sc_min), 0.1) * 0.2;
        %         xl_lo = sc_min - sc_pad;
        %         xl_hi = sc_max + sc_pad;
        %     else
        %         sc_range = sc_max - sc_min;
        %         xl_lo = sc_min - 0.05 * sc_range;
        %         xl_hi = sc_max + 0.05 * sc_range;
        %     end
        %     xl_lo = min(xl_lo, 0); % always show the zero reference
        %     xl_hi = max(xl_hi, 0);
        % 
        %     set(ha_sc, 'YDir', 'reverse', ... % top-to-bottom matches matrix
        %         'YLim', [0.5 N_ROIs + 0.5], ...
        %         'XLim', [xl_lo, xl_hi], ...
        %         'FontSize', 7, ...
        %         'TickDir', 'out', ...
        %         'Box', 'off');
        % 
        %     n_flagged = sum(sign_violation);
        %     if n_flagged > 0
        %         xlabel(ha_sc, 'Self-connection value (EC units) — RED = non-negative (sign violation)', ...
        %             'FontSize', 7);
        %         title(ha_sc, sprintf('Self-connections (raw) — %d/%d flagged (>= 0)', ...
        %             n_flagged, N_ROIs), 'FontSize', 8, 'Color', [0.75 0.10 0.10]);
        %     else
        %         xlabel(ha_sc, 'Self-connection value (EC units); dashed line = prior mean', ...
        %             'FontSize', 7);
        %         title(ha_sc, sprintf('Self-connections (raw); all negative, prior = %.2f', prior_ref), ...
        %             'FontSize', 8);
        %     end
        % 
        %     if ~isempty(blk_lbl)
        %         local_add_sc_group_labels(ha_sc, blk_lbl, blk_s, blk_e, N_ROIs);
        %     else
        %         ha_sc.YTick = [];
        %     end
        % else
        %     axis(ha_sc, 'off');
        %     text(0.5, 0.5, 'Self-connection data N/A', 'Units', 'normalized', ...
        %         'HorizontalAlignment', 'center', 'FontSize', 8, ...
        %         'Color', [0.55 0.55 0.55]);
        % end

        %% Header annotations
        annotation(hfig, 'textbox', [0.04, 0.91, 0.92, 0.05], ...
            'String', sprintf('Subject: %s    ROI set: %s  (%d ROIs)    Run: %s', ...
                sub_id, rs_label, N_ROIs, run_id), ...
            'LineStyle', 'none', 'FontSize', 13, 'FontWeight', 'bold', ...
            'Interpreter', 'none', 'VerticalAlignment', 'middle');
        annotation(hfig, 'textbox', [0.04, 0.87, 0.55, 0.04], ...
            'String', sprintf('N_pass = %d  N_marginal = %d  N_fail = %d', ...
            qc_i.N_pass, qc_i.N_marginal, qc_i.N_fail), ...
            'LineStyle', 'none', 'FontSize', 10, 'Interpreter', 'none', ...
            'VerticalAlignment', 'middle');
        annotation(hfig, 'textbox', [0.60, 0.87, 0.38, 0.04], ...
            'String', sprintf('QC Status:  %s', tier_str), ...
            'LineStyle', 'none', 'FontSize', 10, 'FontWeight', 'bold', ...
            'Color', tier_col, 'Interpreter', 'none', 'VerticalAlignment', 'middle');

        local_export_page(hfig, pdf_path, first_page);
        first_page = false;

    end  % i_rs

    rdcm_log(params, 1, '  Written: %s\n', pdf_path);

end  % i_sub

if ishandle(hfig), close(hfig); end
rdcm_log(params, 1, 'Visual QC complete.\n');

end  % rdcm_visual_qc



%% ====================================================================
%  MATRIX REORDERING + BLOCK LABEL HELPERS
%
%  Permutes both rows and columns of an EC matrix so that all ROIs of
%  the same canonical network form a contiguous block, then draws clean
%  block dividers and labels on both axes.
%
%  Resulting label order (16 blocks for 4S256):
%    Vis-L  Vis-R  SomMot-L  SomMot-R  DAN-L  DAN-R
%    SalVA-L  SalVA-R  Limbic-L  Limbic-R  Cont-L  Cont-R
%    DMN-L  DMN-R  Subcort  Cereb
%
%  ROI name formats handled (both detected automatically):
%    '17Networks_LH_DefaultA_PFCm_1'  (explicit prefix, any atlas copy)
%    'LH_DefaultA_PFCm_1'             (no prefix — XCP-D 4S256 default)
%    'Cerebellum_Left'                (cerebellar)
%    anything else                    → 'Subcort'
%
%  To add another atlas: extend local_roi_to_network with a new branch.
%% ====================================================================


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
         'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
         'FontSize', fs, 'Interpreter', 'none', ...
         'Rotation', x_rot, 'Clipping', 'off');
    % Y-axis (target ROIs)
    text(ha, -1, blk_ctrs(b), blk_lbl{b}, ...
         'HorizontalAlignment', 'right', 'VerticalAlignment', 'middle', ...
         'FontSize', fs, 'Interpreter', 'none', 'Clipping', 'off');
end

end  % local_add_block_labels



function local_add_sc_group_labels(ha, blk_lbl, blk_s, blk_e, N_ROIs)
%% Draw network block dividers and centred labels on both axes.

N_blocks = numel(blk_lbl);
blk_ctrs = (blk_s + blk_e) / 2;

if N_blocks <= 10
    fs = 7;
elseif N_blocks <= 16
    fs = 6;
else
    fs = 5;
end

ha.YTick = blk_ctrs;
ha.YTickLabel = blk_lbl;
ha.TickLabelInterpreter = 'none';
ha.XGrid = 'off';
ha.YGrid = 'off';

xl = xlim(ha);
hold(ha, 'on');
for b = 1:N_blocks - 1
    yb = blk_e(b) + 0.5;
    line(ha, xl, [yb yb], ...
        'Color', [0.55 0.55 0.55], 'LineWidth', 0.5, 'HitTest', 'off');
end
hold(ha, 'off');

set(ha, 'FontSize', fs);
end % locak_add_sc_group_labels

%% ====================================================================
%  ATLAS PARSING
%
%  local_roi_to_network is the ONLY function that knows about atlas
%  label formats.  Everything else in this file is atlas-agnostic.
%
%  Returns [net, hemi] where:
%    net  — one of: 'Vis','SomMot','DAN','SalVA','Limbic','Cont','DMN',
%                   'Subcort', 'Cereb'
%    hemi — 'L', 'R', or '' (empty for subcortical / cerebellar)
%% ====================================================================



%%  PAGE EXPORT
%% ====================================================================

function local_export_page(hfig, pdf_path, is_first)
use_eg = ~verLessThan('matlab', '9.8');   % R2020a = 9.8
if is_first
    if use_eg
        exportgraphics(hfig, pdf_path, 'ContentType', 'vector');
    else
        saveas(hfig, pdf_path);
    end
else
    if use_eg
        exportgraphics(hfig, pdf_path, 'ContentType', 'vector', 'Append', true);
    else
        [fd, fn, ~] = fileparts(pdf_path);
        n = numel(dir(fullfile(fd, [fn '_page*.pdf']))) + 2;
        saveas(hfig, fullfile(fd, sprintf('%s_page%02d.pdf', fn, n)));
    end
end
end


%% ====================================================================
%%  LOCAL HELPERS (field access, colourmap, tier labels)
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

function sub_id = local_get_subid(results, i_sub, N_rs)
sub_id = sprintf('sub%02d', i_sub);
for r = 1:N_rs
    if i_sub <= size(results.grid,2) && r <= size(results.grid,1)
        ec = results.grid{r, i_sub};
        if isstruct(ec) && isfield(ec,'subj') && ~isempty(ec.subj)
            sub_id = ec.subj; return
        end
    end
end
end


function val = local_safe_cell(C, r, i)
if r<=size(C,1) && i<=size(C,2)
    val = C{r,i};
else
    val = [];
end
end

function v = local_safe_idx(arr, i, default)
if i >= 1 && i <= numel(arr)
    v = arr(i);
else
    v = default;
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

function val = local_get_scalar(s, field)
val = NaN;
if isempty(s) || ~isstruct(s) || ~isfield(s, field), return; end
v = s.(field);
if ~isempty(v), val = v(1); end
end

function val = local_get_nested(s, f1, f2)
val = NaN;
if isempty(s) || ~isstruct(s) || ~isfield(s, f1), return; end
sub = s.(f1);
if isempty(sub) || ~isstruct(sub) || ~isfield(sub, f2), return; end
v = sub.(f2);
if ~isempty(v), val = v(1); end
end

function cmap = local_rdbu_cmap(n)
half = ceil(n/2);
blue = [linspace(0.017,1,half)', linspace(0.443,1,half)', linspace(0.690,1,half)'];
red  = [linspace(1,0.698,n-half)', linspace(1,0.094,n-half)', linspace(1,0.168,n-half)'];
cmap = [blue; red];
end

function s = local_tier_str(tier)
switch tier
    case 3,    s = 'PASS';
    case 2,    s = 'REVIEW';
    case 1,    s = 'EXCLUDE (ICC < floor)';
    case 0,    s = 'EXCLUDE';
    otherwise, s = [char(string(tier)) 'UNKNOWN'];
end
end

function c = local_tier_color(tier)
switch tier
    case 3,    c = [0.12 0.52 0.12];
    case 2,    c = [0.75 0.44 0.00];
    case 1,    c = [0.80 0.10 0.10];
    case 0,    c = [0.80 0.10 0.10];
    otherwise, c = [0.40 0.40 0.40];
end
end

function net_labels = local_load_network_labels(params)
%% Load network_labels CSV if specified and present; return {} otherwise.
net_labels = {};
if ~isfield(params, 'filenames') || ~isfield(params.filenames, 'network_labels')
    return
end
fname = params.filenames.network_labels;
if isempty(fname), return; end
fpath = fullfile(params.dirs.root, fname);
if ~isfile(fpath)
    warning('rdcm_visual_qc: network_labels file not found, skipping block labels:\n  %s', fpath);
    return
end
tbl = readtable(fpath, 'ReadVariableNames', false, 'Delimiter', ',');
net_labels = table2cell(tbl(:, 1));
end

function roi_names = local_load_roi_names(params)
%% Load the full ROI name list from params.filenames.roi_list.
roi_names = {};
if ~isfield(params.filenames, 'roi_list') || isempty(params.filenames.roi_list)
    return
end
fpath = fullfile(params.dirs.root, params.filenames.roi_list);
if ~isfile(fpath), return; end
tbl = readtable(fpath, 'ReadVariableNames', false, 'Delimiter', ',');
roi_names = strtrim(table2cell(tbl(:, 1)));
end