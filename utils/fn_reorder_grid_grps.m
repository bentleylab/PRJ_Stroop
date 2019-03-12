function [grid_lab, colors] = fn_reorder_grid_grps(labels, n_row, n_col, sort_unit, direction)
%% Sort data labels by their grid locations
% INPUTS:
%   labels [cell array strs] - labels in sequential numbered order
%   n_row [int] - number of rows in grid
%   n_col [int] - number of columns in grid
%   sort_unit [str] - {'row','col'} dimension to group data along
%   direction [str] - {'up','dn'} sort in increasing or descreasing order

%% Check conditions
if numel(labels)~=n_row*n_col
    warning('WARNING!!! Grid is missing data!');
end

%% Map labels to grid and assign groupings
% Ensure they're in alphanumeric order
labels = fn_sort_labels_alphanum(labels);

% Get grid numbers for each channel
lab_n = zeros(size(labels));
for l = 1:numel(labels)
    lab_n(l) = str2num(labels{l}(regexp(labels{l},'\d')));
end

% Find sorting order
if strcmp(sort_unit,'row')
    grp_sz = n_col;
    n_grps = n_row;
elseif strcmp(sort_unit,'col')
    grp_sz = n_row;
    n_grps = n_col;
else
    error(['Unknown sort_dim ' sort_unit]);
end
if strcmp(direction,'up')
    grp_order = 1:n_grps;
elseif strcmp(direction,'dn')
    grp_order = n_grps:-1:1;
else
    error(['Unknown sorting direction: ' direction]);
end

grid_pos = reshape(1:n_row*n_col,n_row,n_col);
grid_grp = repmat(1:n_grps,grp_sz,1)';
grid_lab = cell(size(labels));
lab_ix = 1;
for grp_ix = 1:n_grps
    grp_id = grp_order(grp_ix);
    grp_nums = grid_pos(grid_grp==grp_id);
    for ch_ix = 1:numel(grp_nums)
        ch_n = grp_nums(ch_ix);
        if any(lab_n==ch_n)
            grid_lab{lab_ix} = labels{lab_n==ch_n};
            lab_ix = lab_ix+1;
        end
    end
%     grid_lab(lab_ix:lab_ix+grp_sz-1) = labels(grid_grp==grp_ix);
%     lab_ix = lab_ix+grp_sz;
    % indexing always based on grp_ix: ((grp_ix-1)*grp_sz)+1:grp_sz*grp_ix
end

%% Provide colors for each row (to plot)
colors = distinguishable_colors(n_grps);

end