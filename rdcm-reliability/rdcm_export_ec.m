function export = rdcm_export_ec(results, params, roi_subsets)
%% Preamble
%{
Exports EC parameter matrices to CSV files, one per ROI subset.
Each CSV has subjects as rows and connections as columns, with a leading
subject-ID column and connection labels as column headers.

Vectorisation uses A(:)' (column-major, source-grouped) to match the
label order from rdcm_connection_labels. Column k of the CSV always
corresponds to label k: 'ROIj->ROIi' where A(i,j) = source j -> target i.

QC tier per subject is read from results.qc{i_rs}.
Only QC-passing subjects are included by default; marginal subjects are
optionally included (params.export.include_marginal, default: true).

Output files written to params.dirs.output/export/:
  ec_<label>.csv         EC parameter matrix (subjects x connections)
  ec_<label>_labels.txt  Connection labels, one per line
  ec_<label>_meta.txt    Provenance: run_id, SNR, N_subs, convention
---------------------------------------------------------------------------
INPUTS
---------------------------------------------------------------------------
results := struct from rdcm_pipeline (.grid, .qc, .fit, .sh, .params)
params:=        struct â€” pipeline params. Key fields:
                  .dirs.output
                  .sub_id_pattern   char regexp (default: 'sub-[A-Za-z0-9]+')
                  .export.include_marginal  logical (default: true)
roi_subsets:=   struct array from params.roi_subsets_file
---------------------------------------------------------------------------
OUTPUTS
---------------------------------------------------------------------------
export := struct array (one element per ROI subset) with fields:
            .label .csv_path .labels_path .meta_path
            .N_subs_exported .connection_labels
%}

%% --- Defaults -------------------------------------------------------

if ~isfield(params, 'sub_id_pattern')
    params.sub_id_pattern = 'sub-[A-Za-z0-9]+';
end
if ~isfield(params, 'export'),                  params.export = struct(); end
if ~isfield(params.export, 'include_marginal'), params.export.include_marginal = true; end

export_dir = fullfile(params.dirs.output, 'export');
if ~isfolder(export_dir), mkdir(export_dir); end

N_roisets = numel(roi_subsets);
export    = struct();

%% --- Loop over ROI subsets ------------------------------------------

for i_rs = 1:N_roisets

    rs       = roi_subsets(i_rs);
    rs_label = rs.label;
    rs_field = matlab.lang.makeValidName(rs_label);
    N_ROIs   = numel(rs.roi_names);

    % Connection labels: column-major, source-grouped
    conn_labels = rdcm_connection_labels(rs.roi_names);
    N_conn      = numel(conn_labels);   % N_ROIs^2

    % QC info for this cell
    qc_cell = results.qc{i_rs};

    % Which subjects to include
    N_subs = qc_cell.N_subs;
    
    % Export ALL subjects — tier column lets downstream code filter
    sub_indices = 1:N_subs;
    N_export = N_subs;
    
    EC_matrix = nan(N_export, N_conn);
    sub_ids   = cell(N_export, 1);
    qc_tiers  = zeros(N_export, 1);

    for k = 1:N_export
        i_sub  = sub_indices(k);
        ec_out = results.grid{i_rs, i_sub};

        if isempty(ec_out)
            sub_ids{k}     = sprintf('MISSING_sub%02d', i_sub);
            EC_matrix(k,:) = nan(1, N_conn);
            continue
        end

        % Vectorise A column-major (source-grouped) â€” matches label order
        A = local_get_A(ec_out);
        A_vec = A(:)';

        if numel(A_vec) ~= N_conn
            error('rdcm_export_ec: A size mismatch for ROI set "%s", sub %d. Expected %d, got %d.', ...
                  rs_label, i_sub, N_conn, numel(A_vec));
        end

        EC_matrix(k, :) = A_vec;
        sub_ids{k}  = rdcm_extract_subid(ec_out.subj, params.sub_id_pattern);
        qc_tiers(k) = qc_cell.subject_tier(i_sub);

    end

    %% --- Write CSV --------------------------------------------------

    csv_path = fullfile(export_dir, sprintf('ec_%s.csv', rs_field));
    fid = fopen(csv_path, 'w');
    if fid == -1
        error('rdcm_export_ec: cannot write CSV: %s', csv_path);
    end

    % Header row
    fprintf(fid, 'subject_id,qc_tier');
    for c = 1:N_conn
        fprintf(fid, ',%s', conn_labels{c});
    end
    fprintf(fid, '\n');

    % Data rows
    for k = 1:N_export
        fprintf(fid, '%s,%d', sub_ids{k}, qc_tiers(k));
        fprintf(fid, ',%.6f', EC_matrix(k, :));
        fprintf(fid, '\n');
    end
    fclose(fid);

    %% --- Write per-subject QC metrics CSV ---------------------------

    qc_path = fullfile(export_dir, sprintf('qc_metrics_%s.csv', rs_field));
    fid = fopen(qc_path, 'w');
    if fid == -1
        error('rdcm_export_ec: cannot write QC metrics CSV: %s', qc_path);
    end
    
    fprintf(fid, 'subject_id,qc_tier,cosine_sim,R2,ICC_off,ICC_diag,pearson_r_off,pearson_r_diag,hard_gate\n');
    
    for k = 1:N_export
        i_sub = sub_indices(k);
    
        fit_k  = results.fit{i_rs, i_sub};
        sh_k   = results.sh{i_rs,  i_sub};
    
        cos_k  = local_get_scalar(fit_k, 'cosine_sim_mean');
        r2_k   = local_get_scalar(fit_k, 'R2_mean');
        icc_o  = local_get_nested(sh_k, 'offdiag', 'icc');
        icc_d  = local_get_nested(sh_k, 'diag',    'icc');
        r_o    = local_get_nested(sh_k, 'offdiag', 'pearson_r');
        r_d    = local_get_nested(sh_k, 'diag',    'pearson_r');
        hg_k   = qc_cell.hard_gate_flag(i_sub);
    
        fprintf(fid, '%s,%d,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f,%d\n', ...
            sub_ids{k}, qc_tiers(k), cos_k, r2_k, icc_o, icc_d, r_o, r_d, hg_k);
    end
    fclose(fid);

    %% --- Write labels sidecar (.txt) --------------------------------

    labels_path = fullfile(export_dir, sprintf('ec_%s_labels.csv', rs_field));
    fid = fopen(labels_path, 'w');
    fprintf(fid, 'index,label,source_roi,target_roi\n');
    for c = 1:N_conn
        lbl = conn_labels{c};                          % format: 'ROIsrc->ROItgt'
        parts = strsplit(lbl, '->');
        fprintf(fid, '%d,%s,%s,%s\n', c, lbl, parts{1}, parts{2});
    end
    fclose(fid);

    %% --- Write metadata sidecar (JSON) ------------------------------

    meta_path = fullfile(export_dir, sprintf('ec_%s_meta.json', rs_field));
    
    run_id = local_get_run_id(params);
    meta = struct();
    meta.roi_set_label      = rs_label;
    meta.run_id             = run_id;
    meta.est_method         = params.rdcm.est_method;
    meta.N_ROIs             = N_ROIs;
    meta.N_connections      = N_conn;
    meta.N_subs_total       = N_subs;
    meta.vectorisation      = 'column-major A(:) source-grouped';
    meta.label_convention   = 'source->target [A(i,j) = source_j -> target_i]';
    meta.exported           = datestr(now, 'yyyy-mm-dd HH:MM:SS');
    meta.roi_names          = rs.roi_names(:)';   % row cell for clean JSON array
    
    fid = fopen(meta_path, 'w');
    if fid == -1
        error('rdcm_export_ec: cannot write meta JSON: %s', meta_path);
    end
    fprintf(fid, '%s\n', jsonencode(meta));
    fclose(fid);

    %% --- Store record -----------------------------------------------

    export(i_rs).label             = rs_label;
    export(i_rs).csv_path          = csv_path;
    export(i_rs).labels_path       = labels_path;
    export(i_rs).meta_path         = meta_path;
    export(i_rs).qc_metrics_path = qc_path;
    export(i_rs).N_subs_exported   = N_export;
    export(i_rs).connection_labels = conn_labels;

    rdcm_log(params, 1, 'Exported EC: %s  (%d subs, %d connections)\n', ...
        csv_path, N_export, N_conn);

end

end

function A = local_get_A(ec_out)
    if isempty(ec_out) || ~isstruct(ec_out)
        A = [];
    elseif isfield(ec_out, 'A') && ~isempty(ec_out.A)
        A = ec_out.A;
    elseif isfield(ec_out, 'Ep') && isstruct(ec_out.Ep) && isfield(ec_out.Ep, 'A')
        A = ec_out.Ep.A;
    else
        error('rdcm_export_ec: no connectivity matrix found (expected ec_out.A or ec_out.Ep.A).');
    end
end

function run_id = local_get_run_id(params)
    if isfield(params, 'log') && isfield(params.log, 'run_id')
        run_id = params.log.run_id;
    elseif isfield(params, 'run_id')
        run_id = params.run_id;
    else
        run_id = 'unknown';
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