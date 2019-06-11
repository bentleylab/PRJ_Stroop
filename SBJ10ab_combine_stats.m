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
shared_fields = {'alpha','model_lab','groups','rt_corr','regress_rt','trial_cond','bad_trials'};
single_fields = {'cust_win','min_rt','win_len','win_step','n_boots','ep_lab'};
tuple_fields  = {'lim_adj','stat_lim','evnt_lab'};

% Compare field names
if ~all(strcmp(fieldnames(full{1}),fieldnames(full{2}))); error('different fields'); end
% Check shared fields are equal
for f_ix = 1:numel(shared_fields)
    if ~isequaln(full{1}.st.(shared_fields{f_ix}),full{2}.st.(shared_fields{f_ix}))
        error(['st.' shared_fields{f_ix} ' mismatch']);
    end
end
% Warn for different events
if ~strcmp(full{1}.st.evnt_lab,full{2}.st.evnt_lab);
    warning('CAREFUL: st.evnt_lab mismatch!');
end

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
    %% Determine order to add w2
    if ~all(strcmp(full{1}.st.evnt_lab,full{2}.st.evnt_lab))
        % S + R - put R last
        if any(strcmp(full{1}.st.evnt_lab,'R'))
            order_idx = [2 1];
        else
            order_idx = [1 2];
        end
    elseif numel(full{1}.w2.time)==1 || numel(full{2}.w2.time)==1
        % S + B/D - order by time
        [~, order_idx] = sort([mean(full{1}.w2.time) mean(full{2}.w2.time)]);
    else
        error('what combination is this?');
    end
    
    %% Combine stat info
    st = struct;
    st.combined_st = 1;                     % denote this is a combined st
    
    % Add new info, being careful to handle previously combined stats well
    comb_already = [isfield(full{1}.st,'combined_st') isfield(full{2}.st,'combined_st')];
    st.an_ids = {}; st.stat_ids = {};
    for st_ix = 1:2
        if comb_already(order_idx(st_ix))
            % Add in all previosuly combined an_ids and stat_ids
            st.an_ids   = [st.an_ids full{order_idx(st_ix)}.st.an_ids{:}];
            st.stat_ids = [st.stat_ids full{order_idx(st_ix)}.st.stat_ids];
        else
            st.an_ids   = [st.an_ids an_ids(order_idx(st_ix))];
            st.stat_ids = [st.stat_ids stat_ids(order_idx(st_ix))];
        end
    end
    
    % Add shared fields (already compatibility checked)
    for f_ix = 1:numel(shared_fields)
        st.(shared_fields{f_ix}) = full{1}.st.(shared_fields{f_ix});
    end
    
    % Add fields with single items
    for f_ix = 1:numel(single_fields)
        st.(single_fields{f_ix}) = [full{order_idx(1)}.st.(single_fields{f_ix}) ...
                                    full{order_idx(2)}.st.(single_fields{f_ix})];
    end
    
    % Add fields with multiple pieces of data (plus evnt_lab to avoid str cat)
%     st.bad_trial_n = {};
    for f_ix = 1:numel(tuple_fields)
        st.(tuple_fields{f_ix}) = {};
        for st_ix = 1:2
            if comb_already(order_idx(st_ix))
                st.(tuple_fields{f_ix}) = [st.(tuple_fields{f_ix}) ...
                    full{order_idx(st_ix)}.st.(tuple_fields{f_ix})];
            else
                st.(tuple_fields{f_ix}) = [st.(tuple_fields{f_ix}) ...
                    {full{order_idx(st_ix)}.st.(tuple_fields{f_ix})}];
            end
        end
    end
    
%     % Variable fields (only .bad_trial_n for now...)
%     for st_ix = 1:2
%         if isfield(full{order_idx(st_ix)}.st,'bad_trial_n')
%             st.bad_trial_n = [st.bad_trial_n {full{order_idx(st_ix)}.st.bad_trial_n}];
%         else
%             st.bad_trial_n = [st.bad_trial_n {zeros([0 1])}];
%         end
%     end
    
    % Check if anything was missed
    if ~isempty(setdiff(fieldnames(full{1}.st),fieldnames(st))) || ...
       ~isempty(setdiff(fieldnames(full{2}.st),fieldnames(st)))
       error('new st combination missed fields!');
    end
    
    %% Concatenate w2 structs
    w2 = full{order_idx(1)}.w2;
    % Keep win_lim separate to know exactly which data
    w2.win_lim = {};
    for st_ix = 1:2
        if comb_already(order_idx(st_ix))
            w2.win_lim   = [w2.win_lim full{order_idx(st_ix)}.w2.win_lim];
        else
            w2.win_lim   = [w2.win_lim {full{order_idx(st_ix)}.w2.win_lim}];
        end
    end
    
    % Concatenate time series data
    w2.time = [full{order_idx(1)}.w2.time full{order_idx(2)}.w2.time];
    w2.win_lim_s = [full{order_idx(1)}.w2.win_lim_s; full{order_idx(2)}.w2.win_lim_s];
    ts_fields = {'boot','trial','pval','qval','zscore','bootmean','bootstd'};
    for f_ix = 1:numel(ts_fields)
        w2.(ts_fields{f_ix}) = cat(3,full{order_idx(1)}.w2.(ts_fields{f_ix}),...
                                     full{order_idx(2)}.w2.(ts_fields{f_ix}));
    end
    
    % Create time series to indicate which windows are custom
    cust_win_ts = cell([2 1]);
    for st_ix = 1:2
        if isfield(full{st_ix}.w2,'cust_win')
            cust_win_ts{st_ix} = full{st_ix}.w2.cust_win;
        elseif full{st_ix}.st.cust_win==1
            cust_win_ts{st_ix} = ones(size(full{st_ix}.w2.time));
            if numel(full{st_ix}.w2.time)>1
                error('multiple custom windows???');
            end
        else
            cust_win_ts{st_ix} = zeros(size(full{st_ix}.w2.time));
        end
    end
    w2.cust_win   = [cust_win_ts{order_idx(1)} cust_win_ts{order_idx(2)}];
    
    % Add stat_id info
    w2.stat_ids = st.stat_ids; w2.an_ids = st.an_ids;
    stat_id_ts = cell([2 1]);
    for st_ix = 1:2
        if isfield(full{st_ix}.w2,'stat_id_ts')
            stat_id_ts{st_ix} = full{st_ix}.w2.stat_id_ts;
            % Re-map stat_ids to new combination
            old_st_idx = unique(full{st_ix}.w2.stat_id_ts);
            for sub_ix = 1:numel(full{st_ix}.w2.stat_ids)
                stat_id_ts{st_ix}(stat_id_ts{st_ix}==old_st_idx(sub_ix)) = find(strcmp(w2.stat_ids,full{st_ix}.w2.stat_ids{sub_ix}));
            end
        else
            stat_id_ts{st_ix} = ones(size(full{st_ix}.w2.time))*find(strcmp(stat_ids{st_ix},st.stat_ids));
        end
    end
    w2.stat_id_ts = [stat_id_ts{order_idx(1)} stat_id_ts{order_idx(2)}];
end

%% Save stat
% Find matching substring starting from beginning
add_name = stat_ids{order_idx(2)};
for let_ix = 1:numel(stat_ids{order_idx(2)})
    if strcmp(stat_ids{order_idx(1)}(1:let_ix),stat_ids{order_idx(2)}(1:let_ix))
        add_name = add_name(2:end);
    else
        break;
    end
end

% Save combined data
new_name = [stat_ids{order_idx(1)} '_' add_name];
if strcmp(full{1}.st.model_lab,'actv')
    error('actv needs implementation');
else
    out_fname = [SBJ_vars.dirs.proc SBJ '_smANOVA_ROI_' new_name '_' an_ids{order_idx(1)} '.mat'];
    if st.rt_corr
        save(out_fname,'-v7.3','w2','stat','st');
    else
        save(out_fname,'-v7.3','w2','st');
    end
end

end