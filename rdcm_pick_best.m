function best = rdcm_pick_best(qc_row, fit_slice, sh_slice, params, rs_i)
%% rdcm_pick_best
%{
Selects the best SNR value for a single ROI set.

Called from rdcm_pipeline after reshaping:
    fit_slice_2d = reshape(cell_fit(i_rs,:,:), N_snr, N_subs);
    sh_slice_2d  = reshape(cell_sh(i_rs,:,:),  N_snr, N_subs);
    best_rs = rdcm_pick_best(results.qc(i_rs,:), fit_slice_2d, sh_slice_2d, params, roi_subsets(i_rs));

INPUTS
    qc_row      {1 x N_snr}      QC structs from rdcm_qc_group
                                  Fields: .N_pass .N_marginal .N_subs .subject_tier
    fit_slice   {N_snr x N_subs} fit_out structs (reserved for future use)
    sh_slice    {N_snr x N_subs} sh_out structs from rdcm_splithalf_single
    params      struct            .pick.criterion = icc|cosine_sim|pearson_r
                                  .rdcm.SNR_grid
    rs_i        struct            roi_subsets(i_rs) — label used in warnings

OUTPUTS
    best        struct
        .SNR              winning SNR value
        .snr_idx          index into SNR_grid
        .criterion_value  mean criterion among pass+marginal subjects
        .N_pass           N_pass at winning SNR
        .N_marginal       N_marginal at winning SNR
        .N_subs           total subjects
        .label            ROI set label
        .criterion        criterion name used
%}

%% --- Defaults -------------------------------------------------------
if ~isfield(params, 'pick'),           params.pick = struct(); end
if ~isfield(params.pick, 'criterion'), params.pick.criterion = 'icc'; end
if ~isfield(params, 'rdcm'),           params.rdcm = struct(); end
crit = params.pick.criterion;

%% --- Safety squeeze: accept {1 x N_snr x N_subs} or {N_snr x N_subs}
sh_slice  = sh_slice(:,:);
fit_slice = fit_slice(:,:);

N_snr = numel(qc_row);
N_subs = size(sh_slice, 2);

rs_label = 'unknown';
if nargin >= 5 && isstruct(rs_i) && isfield(rs_i, 'label')
    rs_label = rs_i.label;
end

%% --- SNR grid -------------------------------------------------------
if isfield(params.rdcm, 'SNR_grid') && ~isempty(params.rdcm.SNR_grid)
    SNR_grid = params.rdcm.SNR_grid;
else
    SNR_grid = 1:N_snr;
    warning('rdcm_pick_best [%s]: SNR_grid not found; using indices.', rs_label);
end

%% --- Score each SNR -------------------------------------------------
scores     = nan(N_snr, 1);
n_pass_arr = nan(N_snr, 1);
n_marg_arr = nan(N_snr, 1);
n_subs_arr = nan(N_snr, 1);

for k = 1:N_snr
    qc_k = qc_row{k};
    if isempty(qc_k), continue; end

    n_pass_arr(k) = qc_k.N_pass;
    n_marg_arr(k) = qc_k.N_marginal;
    n_subs_arr(k) = qc_k.N_subs;

    usable    = qc_k.subject_tier >= 2;   % tier 2 (marginal) or 3 (pass)
    crit_vals = nan(N_subs, 1);

    for i = 1:N_subs
        if ~usable(i), continue; end
        sh_i = sh_slice{k, i};
        if isempty(sh_i) || ~isstruct(sh_i), continue; end

        switch crit
            case 'icc'
                if isfield(sh_i, 'offdiag') && isfield(sh_i.offdiag, 'icc')
                    crit_vals(i) = sh_i.offdiag.icc;
                elseif isfield(sh_i, 'icc')
                    crit_vals(i) = sh_i.icc;
                end
            case 'cosine_sim'
                if isfield(sh_i, 'ec_cosine_sim')
                    crit_vals(i) = sh_i.ec_cosine_sim;
                elseif isfield(sh_i, 'cosine_sim')
                    crit_vals(i) = sh_i.cosine_sim;
                end
            case 'pearson_r'
                if isfield(sh_i, 'offdiag') && isfield(sh_i.offdiag, 'pearson_r')
                    crit_vals(i) = sh_i.offdiag.pearson_r;
                elseif isfield(sh_i, 'pearson_r')
                    crit_vals(i) = sh_i.pearson_r;
                end
            otherwise
                error('rdcm_pick_best: unknown criterion "%s".', crit);
        end
    end
    scores(k) = mean(crit_vals, 'omitnan');
end

%% --- Rank: N_pass desc, criterion desc tiebreak --------------------
valid = ~isnan(scores);
if ~any(valid)
    warning('rdcm_pick_best [%s]: all NaN — defaulting to index 1.', rs_label);
    valid(1) = true; scores(1) = 0;
end
idx_v  = find(valid);
[~, o] = sortrows([n_pass_arr(idx_v), scores(idx_v)], [-1 -2]);
best_k = idx_v(o(1));

%% --- Output ---------------------------------------------------------
best.SNR             = SNR_grid(best_k);
best.snr_idx         = best_k;
best.criterion_value = scores(best_k);
best.N_pass          = n_pass_arr(best_k);
best.N_marginal      = n_marg_arr(best_k);
best.N_subs          = n_subs_arr(best_k);
best.label           = rs_label;
best.criterion       = crit;

end
