function labels = rdcm_connection_labels(roi_names)
%% Preamble
%{
Generates connection labels for a vectorised rDCM A matrix, matching the
column order produced by A(:)' (MATLAB column-major vectorisation).

In rDCM/DCM, A(i,j) is the influence of source region j on target region i.
Column-major vectorisation groups all targets for each source together:

  A(:)' = [A(1,1), A(2,1), ..., A(N,1),   <- all targets for source 1
            A(1,2), A(2,2), ..., A(N,2),   <- all targets for source 2
            ...
            A(1,N), A(2,N), ..., A(N,N)]   <- all targets for source N

Labels are formatted as 'source->target'. Outer loop = source (column j
of A), inner loop = target (row i of A), so label k corresponds to
A(i,j) and reads 'ROIj->ROIi'.

Called by rdcm_export_ec and rdcm_report_qc to ensure CSV column headers
and reliability matrices use identical labelling.
---------------------------------------------------------------------------
INPUTS
---------------------------------------------------------------------------
roi_names:=     {1 x N_ROIs} cell array of ROI label strings
---------------------------------------------------------------------------
OUTPUTS
---------------------------------------------------------------------------
labels:=        {1 x N_ROIs^2} cell array of connection label strings
                in column-major (source-grouped) order matching A(:)'
%}

N      = numel(roi_names);
labels = cell(1, N * N);
k      = 0;
for j = 1:N          % source (column of A)
    for i = 1:N      % target (row of A)
        k         = k + 1;
        labels{k} = sprintf('%s->%s', roi_names{j}, roi_names{i});
    end
end

end