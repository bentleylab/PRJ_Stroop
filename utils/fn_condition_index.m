function cond_idx = fn_condition_index(cond_lab, cond_n, varargin)
% Returns binary truth array of trials matching a certain condition
% INPUTS:
%   cond_lab [str] - one of [C, N, I, MC, EQ, MI, ${trial-block}]
%       ${trial-block} can be any combination, e.g., "C-EQ"
%   cond_n [int array] - numeric condition labels, has length(n_trials)
%       C = 1-3, N = 4-6, I = 7-9
%       those are ordered: EQ, MI, MC
%       e.g., 5 is MI-N
%   varargin - 'trial_info' passes in the full struct
% OUTPUTS:
%   cond_idx [binary array] - 0/1 for identity with condition_label, length(n_trials)

% Handle variable inputs
if ~isempty(varargin)
    for v = 1:2:numel(varargin)
        if strcmp(varargin{v},'trial_info')
            trial_info = varargin{v+1};
        else
            error(['Unknown varargin ' num2str(v) ': ' varargin{v}]);
        end
    end
end

% Check trial_info was passed in for sequences
if any(strcmp({'pC','pN','pI','cI','iI'},cond_lab)) && ~exist('trial_info','var')
    error(['trial_info must be supplied to find condition: ' cond_lab]);
end

% Make sure cond_n is a column vector
if size(cond_n,2)~=1
    cond_n = cond_n';
end

switch cond_lab
    % Current Trial Congruence
    case 'C'
        cond_idx = any([cond_n==1 cond_n==2 cond_n==3],2);
    case 'N'
        cond_idx = any([cond_n==4 cond_n==5 cond_n==6],2);
    case 'I'
        cond_idx = any([cond_n==7 cond_n==8 cond_n==9],2);
    
    % Previous Trial Congruence
    case 'pC'
        C_idx = any([cond_n==1 cond_n==2 cond_n==3],2);
        block_start_ix = find(diff(trial_info.block_n)>0)+1;
        cond_idx = false(size(cond_n));
        for t = setdiff(2:numel(cond_n),block_start_ix)
            if C_idx(t-1)==1
                cond_idx(t) = 1;
            end
        end
    case 'pN'
        N_idx = any([cond_n==4 cond_n==5 cond_n==6],2);
        block_start_ix = find(diff(trial_info.block_n)>0)+1;
        cond_idx = false(size(cond_n));
        for t = setdiff(2:numel(cond_n),block_start_ix)
            if N_idx(t-1)==1
                cond_idx(t) = 1;
            end
        end
    case 'pI'
        I_idx = any([cond_n==7 cond_n==8 cond_n==9],2);
        block_start_ix = find(diff(trial_info.block_n)>0)+1;
        cond_idx = false(size(cond_n));
        for t = setdiff(2:numel(cond_n),block_start_ix)
            if I_idx(t-1)==1
                cond_idx(t) = 1;
            end
        end
        
    % Proportion Congruence Blocks
    case 'MC'
        cond_idx = any([cond_n==3 cond_n==6 cond_n==9],2);
    case 'EQ'
        cond_idx = any([cond_n==1 cond_n==4 cond_n==7],2);
    case 'MI'
        cond_idx = any([cond_n==2 cond_n==5 cond_n==8],2);

    % Trial and Block Combinations
    case 'C_EQ'
        cond_idx = [cond_n==1];
    case 'C_MI'
        cond_idx = [cond_n==2];
    case 'C_MC'
        cond_idx = [cond_n==3];
    case 'N_EQ'
        cond_idx = [cond_n==4];
    case 'N_MI'
        cond_idx = [cond_n==5];
    case 'N_MC'
        cond_idx = [cond_n==6];
    case 'I_EQ'
        cond_idx = [cond_n==7];
    case 'I_MI'
        cond_idx = [cond_n==8];
    case 'I_MC'
        cond_idx = [cond_n==9];
    
    % Congruence Sequence Effects
    case 'cI'
        C_idx = any([cond_n==1 cond_n==2 cond_n==3],2);
        I_idx = any([cond_n==7 cond_n==8 cond_n==9],2);
        block_start_ix = find(diff(trial_info.block_n)>0)+1;
        cond_idx = false(size(cond_n));
        for t = setdiff(2:numel(cond_n),block_start_ix)
            if C_idx(t-1)==1 && I_idx(t)==1
                cond_idx(t) = 1;
            end
        end
    case 'iI'
        I_idx = any([cond_n==7 cond_n==8 cond_n==9],2);
        block_start_ix = find(diff(trial_info.block_n)>0)+1;
        cond_idx = false(size(cond_n));
        for t = setdiff(2:numel(cond_n),block_start_ix)
            if I_idx(t-1)==1 && I_idx(t)==1
                cond_idx(t) = 1;
            end
        end
    otherwise
        error('Invalid condition_label');
end

end
