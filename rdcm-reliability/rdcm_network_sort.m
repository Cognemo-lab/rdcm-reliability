function [perm_idx, blk_lbl, blk_s, blk_e] = rdcm_network_sort(roi_names, network_labels)
%% RDCM_NETWORK_SORT  Sort ROI indices into contiguous network blocks.
%
%  [perm_idx, blk_lbl, blk_s, blk_e] = rdcm_network_sort(roi_names, network_labels)
%
%  INPUTS
%    roi_names      — {N x 1} cell of ROI name strings (from rs.roi_names)
%    network_labels — {N x 1} cell of group label strings, same order as
%                     roi_names. May be [] or {} to request no reordering.
%
%  OUTPUTS
%    perm_idx — [N x 1] permutation index (identity if no labels given)
%    blk_lbl  — {1 x B} cell of block label strings (empty if no labels)
%    blk_s    — [1 x B] block start indices in permuted order
%    blk_e    — [1 x B] block end indices in permuted order
%
%  BEHAVIOUR
%    - If network_labels is empty: returns identity permutation and empty
%      block arrays → caller suppresses axis labels entirely.
%    - If labels are already contiguous (e.g. Harvard-Oxford anatomical
%      regions): order is preserved as-is.
%    - If labels are non-contiguous (e.g. Schaefer LH/RH interleaving):
%      ROIs are stably sorted by label so all members of each group are
%      contiguous. Within each group the original ROI order is preserved.
%    - Group order in the output follows the first appearance of each
%      label in the original network_labels vector (stable sort).

N = numel(roi_names);

%% No-label fallback
if isempty(network_labels)
    perm_idx = (1:N)';
    blk_lbl  = {};
    blk_s    = [];
    blk_e    = [];
    return
end

if numel(network_labels) ~= N
    error('rdcm_network_sort: network_labels length (%d) must match roi_names length (%d).', ...
          numel(network_labels), N);
end

%% Determine group order by first appearance (preserves intended ordering)
[~, first_idx] = unique(network_labels, 'stable');
group_order    = network_labels(sort(first_idx));   % unique labels in first-appearance order

%% Build permutation: stable sort within each group
perm_idx = zeros(N, 1);
blk_lbl  = {};
blk_s    = [];
blk_e    = [];
ptr      = 0;

for k = 1:numel(group_order)
    lbl = group_order{k};
    idx = find(strcmp(network_labels, lbl));  % already in original order (stable)
    if isempty(idx), continue; end
    n_k               = numel(idx);
    blk_s(end+1)      = ptr + 1;             %#ok<AGROW>
    perm_idx(ptr+1 : ptr+n_k) = idx;
    ptr               = ptr + n_k;
    blk_e(end+1)      = ptr;                 %#ok<AGROW>
    blk_lbl{end+1}    = lbl;                 %#ok<AGROW>
end

end % rdcm_network_sort