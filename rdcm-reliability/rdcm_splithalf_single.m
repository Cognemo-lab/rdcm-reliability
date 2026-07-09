function sh = rdcm_splithalf_single(Y, params)
%% Preamble
%{
Runs a split-half reliability analysis for a single subject's BOLD
timeseries. The timeseries is split into two contiguous halves, rDCM is
estimated independently on each half, and then three reliability measures
are computed between the two resulting EC matrices:

(1) Splits the timeseries into two contiguous halves
(2) Runs rDCM independently on each half
(3) Returns the raw EC matrices and fit structs for group-level QC

All reliability metrics are computed downstream in rdcm_qc_group.
---------------------------------------------------------------------------
INPUTS
---------------------------------------------------------------------------
Y      := struct  Subject BOLD data (see rdcm_run_single for field list).

params := struct  Pipeline params struct. Relevant subfields:
           .rdcm.SNR              scalar  (default: 1)
           .rdcm.est_method       scalar  (default: 1)
           .splithalf.min_timepoints scalar (default: 50)
           .splithalf.exclude_diagonal logical (default: false)
           .verbose               integer (default: 2)
---------------------------------------------------------------------------
OUTPUTS
---------------------------------------------------------------------------
sh := struct with fields:

Raw EC matrices:
  .A1          [N_ROIs x N_ROIs]  EC matrix from first half
  .A2          [N_ROIs x N_ROIs]  EC matrix from second half

Metadata:
  .N_TRs_half  scalar
  .N_ROIs      scalar
  .exclude_diagonal  logical
  .subj        char
  .warning     char

NOTE: All reliability metrics (Pearson r, EC-space cosine similarity,
ICC(3,1) with Spearman-Brown correction) are computed in rdcm_qc_group,
which has access to the full subject pool needed for group-level ICC.


%}

%% --- Defaults -------------------------------------------------------

if ~isfield(params, 'rdcm'),           params.rdcm = struct(); end
if ~isfield(params.rdcm, 'SNR'),       params.rdcm.SNR = 1;    end
if ~isfield(params.rdcm, 'est_method'),params.rdcm.est_method = 1; end
if ~isfield(params, 'splithalf'),      params.splithalf = struct(); end
if ~isfield(params.splithalf, 'min_timepoints')
    params.splithalf.min_timepoints = 50;
end
if ~isfield(params.splithalf, 'exclude_diagonal')
    params.splithalf.exclude_diagonal = false;
end
if ~isfield(params, 'verbose'), params.verbose = 2; end

sh           = struct();
sh.subj      = Y.subj;
sh.N_ROIs    = size(Y.y, 2);
sh.exclude_diagonal = params.splithalf.exclude_diagonal;
sh.warning   = '';

%% --- Validate and Trim Timeseries -----------------------------------

N_TRs = size(Y.y, 1);

if mod(N_TRs, 2) ~= 0
    Y.y   = Y.y(1:end-1, :);
    N_TRs = N_TRs - 1;
    msg   = sprintf('Odd N_TRs: dropped last timepoint (now N=%d).', N_TRs);
    warning('rdcm_splithalf_single:oddTRs', '%s', msg);
    sh.warning = msg;
end

N_half         = N_TRs / 2;
sh.N_TRs_half  = N_half;

if N_half < params.splithalf.min_timepoints
    error(['rdcm_splithalf_single: half-length (%d TRs) is below minimum ' ...
        '(%d TRs). Increase scan length or lower ' ...
        'params.splithalf.min_timepoints.'], ...
        N_half, params.splithalf.min_timepoints);
end

%% --- Build Half-Session Y Structs -----------------------------------

Y1   = Y;  Y1.y = Y.y(1:N_half,     :);
Y2   = Y;  Y2.y = Y.y(N_half+1:end, :);

%% --- Run rDCM on Each Half ------------------------------------------
% rdcm_run_single(Y, params) now returns (EC_out, fit_out).

if params.verbose >= 2
    fprintf('  [split-half] Subject %s — half 1 of 2...\n', Y.subj);
end
[EC1,~] = rdcm_run_single(Y1, params);

if params.verbose >= 2
    fprintf('  [split-half] Subject %s — half 2 of 2...\n', Y.subj);
end
[EC2,~] = rdcm_run_single(Y2, params);

%% --- Store Raw Outputs ----------------------------------------------

sh.A1 = EC1.A;
sh.A2 = EC2.A;

%% --- Console Summary ------------------------------------------------

if params.verbose >= 2
    fprintf(' [split-half] %s done. Both halves estimated (N_half=%d TRs).\n', ...
        Y.subj, N_half);
end

end