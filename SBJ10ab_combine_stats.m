function SBJ10ab_combine_stats(SBJ,an_id1,an_id2,stat_id1,stat_id2)
%%
% Combine output data of two stat structures
%   Currently only intended for combining S and D stat structs

%% Data Preparation
% Set up paths
[root_dir, app_dir] = fn_get_root_dir(); ft_dir = [app_dir 'fieldtrip/'];
addpath([root_dir 'PRJ_Stroop/scripts/']);
addpath([root_dir 'PRJ_Stroop/scripts/utils/']);
addpath(ft_dir);
ft_defaults

%% Load Results
SBJ_vars_cmd = ['run ' root_dir 'PRJ_Stroop/scripts/SBJ_vars/' SBJ '_vars.m'];
eval(SBJ_vars_cmd);

% Load data
an_ids = {an_id1 an_id2}; stat_ids = {stat_id1 stat_id2};
full = cell(size(stat_ids));
for st_ix = 1:2
    if strcmp(stat_ids{st_ix}(1:4),'actv')
        stat_fname = [SBJ_vars.dirs.proc SBJ '_ROI_' an_ids{st_ix} '_' stat_ids{st_ix} '.mat'];
    else
        stat_fname = [SBJ_vars.dirs.proc SBJ '_smANOVA_ROI_' stat_ids{st_ix} '_' an_ids{st_ix} '.mat'];
    end
    full{st_ix} = load(stat_fname);
end

%% Check compatibility
if ~all(strcmp(fieldnames(full{1}),fieldnames(full{2}))); error('different fields'); end
if ~strcmp(full{1}.st.evnt_lab,full{2}.st.evnt_lab); error('st.evnt_lab mismatch'); end
if full{1}.st.alpha~=full{2}.st.alpha; error('st.alpha mismatch'); end
if ~strcmp(full{1}.st.model_lab,full{2}.st.model_lab); error('st.model_lab mismatch'); end
if ~all(strcmp(full{1}.st.groups,full{2}.st.groups)); error('st.groups mismatch'); end
if full{1}.st.rt_corr~=full{2}.st.rt_corr; error('st.alpha mismatch'); end

if strcmp(full{1}.st.model_lab,'actv')
    error('actv needs implementation');
    % consider doing this with full{st_ix}.(st_name) since fileds should be
    % same names...
elseif full{1}.st.rt_corr
    error('not ready to combine rt_corr stats');
else
    if ~all(strcmp(full{1}.w2.label,full{2}.w2.label)); error('w2.label mismatch'); end
    if ~strcmp(full{1}.w2.dimord,full{2}.w2.dimord); error('w2.dimord mismatch'); end
    if ~strcmp(full{1}.w2.dimord(end-3:end),'time'); error('w2 last dim not time'); end
    for grp_ix = 1:numel(full{1}.st.groups)
        if ~all(full{1}.w2.design{grp_ix}==full{2}.w2.design{grp_ix})
            error('w2.design doesn''t match');
        end
    end
end

%% Combine stat structs
if strcmp(full{1}.st.model_lab,'actv')
    error('actv needs implementation');
else
    % Choose stat with most time points as main struct basis
    [~,main_ix] = max([numel(full{1}.w2.time),numel(full{2}.w2.time)]);
    add_ix = setdiff([1 2],main_ix);
    
    % Combine stat info
    st = full{main_ix}.st;
    st.(['st_' stat_ids{add_ix}]) = full{add_ix}.st;
    st.(['an_' stat_ids{add_ix}]) = an_ids{add_ix};
    
    % Determine order to add w2
    if numel(full{add_ix}.w2.time)>1
        error('not ready to deal with multiple add ons yet!');
    end
    if full{add_ix}.w2.time < min(full{main_ix}.w2.time)
        order_idx = [add_ix main_ix];
    elseif full{add_ix}.w2.time > max(full{main_ix}.w2.time)
        order_idx = [main_ix add_ix];
    end
    
    % Concatenate w2 structs
    w2 = full{main_ix}.w2;
    w2.time = [full{order_idx(1)}.w2.time full{order_idx(2)}.w2.time];
    w2.(['win_lim_' stat_ids{add_ix}]) = full{add_ix}.w2.win_lim;
    w2.(['win_lim_s_' stat_ids{add_ix}]) = full{add_ix}.w2.win_lim_s;
    ts_fields = {'boot','trial','pval','qval','zscore','bootmean','bootstd'};
    for f_ix = 1:numel(ts_fields)
        w2.(ts_fields{f_ix}) = cat(3,full{order_idx(1)}.w2.(ts_fields{f_ix}),...
                                     full{order_idx(2)}.w2.(ts_fields{f_ix}));
    end
    
    % Add field to adjust win_len for new data
    w2.cust_win = zeros(size(w2.time));
    w2.cust_win(w2.time==full{add_ix}.w2.time) = 1;
    w2.cust_win_len = 0;
end

%% Save stat
% Find matching substring starting from beginning
add_name = stat_ids{add_ix};
for let_ix = 1:numel(stat_ids{add_ix})
    if strcmp(stat_ids{main_ix}(1:let_ix),stat_ids{add_ix}(1:let_ix))
        add_name = add_name(2:end);
    end
end

% Save combined data
new_name = [stat_ids{main_ix} '_' add_name];
if strcmp(full{1}.st.model_lab,'actv')
    error('actv needs implementation');
else
    out_fname = [SBJ_vars.dirs.proc SBJ '_smANOVA_ROI_' new_name '_' an_id '.mat'];
    if st.rt_corr
        save(out_fname,'-v7.3','w2','stat','st');
    else
        save(out_fname,'-v7.3','w2','st');
    end
end

end