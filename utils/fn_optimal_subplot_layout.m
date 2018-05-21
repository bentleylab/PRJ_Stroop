function [n_rows,n_cols] = fn_optimal_subplot_layout(n_plots)
% Returns the optimal number of rows and columns given a number of plots
% INPUTS:
%   n_plots [int] - # of plots, 1<=N<=14

switch n_plots
    case 1
        n_rows=1;n_cols=1;
    case 2
        n_rows=2;n_cols=1;  % Prefer stacking rows so time series can be compared easily
    case 3
        n_rows=3;n_cols=1;
    case 4
        n_rows=2;n_cols=2;
    case 5
        n_rows=3;n_cols=2;
    otherwise
        error('n_plots out of range! Must be 1<=n<=5 as of now...');
end

end