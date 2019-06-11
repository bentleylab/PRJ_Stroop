function SBJ08b_HFA_plot_SR_ERPstack_cond(SBJ,conditions,an_id_s,an_id_r,stat_id_s,stat_id_r,...
                                        plt_id,save_fig,fig_vis,fig_ftype)
% Plots single trial stack for both stimulus- and response-locked HFA computed in SBJ08a_HFA_actv
%   sorts by condition, then by RT; scatter for RTs in stim-locked
% clear all; %close all;

if exist('/home/knight/hoycw/','dir');root_dir='/home/knight/hoycw/';ft_dir=[root_dir 'Apps/fieldtrip/'];
else root_dir='/Volumes/hoycw_clust/';ft_dir='/Users/colinhoy/Code/Apps/fieldtrip/';end

%% Set up paths
addpath([root_dir 'PRJ_Stroop/scripts/']);
addpath([root_dir 'PRJ_Stroop/scripts/utils/']);
addpath(ft_dir);
ft_defaults

%% Handle variable inputs
% SGE inputs
if ischar(save_fig); save_fig = str2num(save_fig); end

%% Load Results
% Load analysis parameters
SBJ_vars_cmd = ['run ' root_dir 'PRJ_Stroop/scripts/SBJ_vars/' SBJ '_vars.m'];
eval(SBJ_vars_cmd);
plt_vars_cmd = ['run ' root_dir 'PRJ_Stroop/scripts/plt_vars/' plt_id '_vars.m'];
eval(plt_vars_cmd);

% Event labels
evnt_lab = {'S', 'R'};
an_ids   = {an_id_s an_id_r}; stat_ids = {stat_id_s stat_id_r};
[cond_lab, cond_colors, cond_style] = fn_condition_label_styles(conditions);

% Load RTs
load(strcat(SBJ_vars.dirs.events,SBJ,'_trial_info_final.mat'),'trial_info');
srate = trial_info.sample_rate;

% Load data
hfa        = cell([2 1]);
evnt_check = cell([2 1]);
alpha_chk  = cell([2 1]);
grp_check  = cell([2 1]);
for sr_ix = 1:2
    % Load HFA
    hfa_fname = strcat(SBJ_vars.dirs.proc,SBJ,'_ROI_',an_ids{sr_ix},'.mat');
    tmp = load(hfa_fname,'hfa'); hfa{sr_ix} = tmp.hfa;
    tmp = load(hfa_fname,'an'); evnt_check{1} = tmp.an.evnt_lab;
    % Load Stats
    if strcmp(conditions,'actv')
        stat_fname = strcat(SBJ_vars.dirs.proc,SBJ,'_ROI_',an_ids{sr_ix},'_',stat_ids{sr_ix},'.mat');
        tmp = load(stat_fname,'actv'); actv{sr_ix} = tmp.actv;
    elseif any(strcmp(conditions,{'CNI','pCNI','PC'}))
        if ~isempty(strfind(stat_ids{sr_ix},'D1tRT')) && ~isempty(strfind(stat_ids{sr_ix},'S0tmRT'))
            an_ids{sr_ix} = strrep(an_ids{sr_ix},'S2t251','S2t151');
        end
        stat_fname = [SBJ_vars.dirs.proc SBJ '_smANOVA_ROI_' stat_ids{sr_ix} '_' an_ids{sr_ix} '.mat'];
        tmp = load(stat_fname,'w2'); w2{sr_ix} = tmp.w2;
    elseif strcmp(conditions,'RT')
        error('RT not yet implemented');
    else
        error(['unknown stat_id:' stat_id_s ' and ' stat_id_r]);
    end
    % Check events are the same
    load(stat_fname,'st');
    evnt_check{2}    = st.evnt_lab;
    alpha_chk{sr_ix} = st.alpha;
    grp_check{sr_ix} = st.groups;
    [grp_lab, ~, ~] = fn_group_label_styles(st.model_lab);
    if ~strcmp(evnt_check{1},evnt_check{2}); error('evnt_lab mismatch'); end
    clear tmp
end
% Check groups are the same
if alpha_chk{1}~=alpha_chk{2}; error('st.alpha mismatch'); end
if ~all(strcmp(grp_check{1},grp_check{2})); error('model group mismatch'); end

% Load elec
elec_fname = [SBJ_vars.dirs.recon SBJ '_elec_main_ft_pat_Dx_full.mat'];
load(elec_fname);
% Sort elecs by stat labels
cfgs = []; cfgs.channel = hfa{1}.label;
elec = fn_select_elec(cfgs,elec);
elec.roi = fn_atlas2roi_labels(elec.atlas_lab,elec.atlas_id,'ROI');

%% Quality checks
% Check channel match between analyses
if ~isempty(setdiff(hfa{1}.label,hfa{2}.label)); error('an label mismatch'); end
if strcmp(conditions,'RT');             lab_check = stat{1}.label;
elseif strcmp(conditions,'actv');       lab_check = actv{1}.label;
elseif any(strcmp(conditions,grp_lab)); lab_check = w2{1}.label; end
if any(~strcmp(hfa{1}.label,lab_check)); error('label mismatch'); end

% If singe trial format, check all time axes are the same
time_cells = [iscell(hfa{1}.time), iscell(hfa{2}.time)];
if any(time_cells)
    for sr_ix = find(time_cells)
        for t_ix = 2:numel(hfa{sr_ix}.time)
            if ~all(hfa{sr_ix}.time{1}==hfa{sr_ix}.time{t_ix})
                error('Time axes aren''t the same for all trials!');
            end
        end
    end
else
    % Confirm hfa starts before and ends after stat
    for sr_ix = 1:numel(evnt_lab)
        if strcmp(conditions,'RT')
            if ~all(hfa{sr_ix}.time==stat{sr_ix}.time)
                error('Time axes arent the same for hfa and RT stat!');
            end
        elseif strcmp(conditions,'actv')
            % Allow for some imprecision (error if off by at least 1 sample)
            if actv{sr_ix}.time(1)-hfa{sr_ix}.time(1)<=-1/srate || ...
                    actv{sr_ix}.time(end)-hfa{sr_ix}.time(end)>=-1/srate
                error('Time axes arent the same for hfa and actv!');
            end
        elseif any(strcmp(conditions,grp_lab))
            % Allow for some imprecision (error if off by at least 1 sample)
            if (w2{sr_ix}.time(1)-st.win_len/2)-hfa{sr_ix}.time(1)<=-1/srate || ...
                    (w2{sr_ix}.time(end)+st.win_len/2)-hfa{sr_ix}.time(end)>=-1/srate
                error('Time axes arent the same for hfa and ANOVA w2!');
            end
        end
    end
end

%% Prep Data
% Trim data to plotting epoch
cfg_trim = [];
cfg_trim.latency = plt_vars.plt_lim_S;
hfa{1} = ft_selectdata(cfg_trim,hfa{1});
cfg_trim.latency = plt_vars.plt_lim_R;
hfa{2} = ft_selectdata(cfg_trim,hfa{2});

% Compile cond_type, RT, trial_n
if strcmp(conditions,'actv')
    [sort_rts, sort_rt_idx] = sort(round(srate*trial_info.response_time));
    % Fake cond_mat with only sorted RTs
    cond_mat = [ones([numel(trial_info.trial_n) 1]) sort_rts sort_rt_idx];
else
    cond_mat = zeros(size(trial_info.trial_n));
    for cond_ix = 1:length(cond_lab)
        % Get binary condition index
        cond_mat(logical(fn_condition_index(cond_lab{cond_ix},...
            trial_info.condition_n,'trial_info',trial_info))) = cond_ix;
    end
    
    % Create matrix of [cond_idx, rt, trial_ix]
    cond_mat = horzcat(cond_mat,round(srate*trial_info.response_time),[1:numel(cond_mat)]');
    cond_mat = sortrows(cond_mat,[1 2]);
end

%% Prep stats
if strcmp(conditions,'RT')
    % Trim data to plotting epoch
    cfg_trim = [];
    cfg_trim.latency = plt_vars.plt_lim_S;
    stat{1} = ft_selectdata(cfg_trim,stat{1});
    cfg_trim.latency = plt_vars.plt_lim_R;
    stat{2} = ft_selectdata(cfg_trim,stat{2});
elseif strcmp(conditions,'actv')
    % Convert to mask time series
    stat = cell(size(evnt_lab));
    for sr_ix = 1:numel(evnt_lab)
        % Prepare stat mask to match hfa for specific factor
        stat{sr_ix} = rmfield(hfa{sr_ix},{'powspctrm','freq','cumtapcnt'});
        stat{sr_ix}.mask = zeros([numel(actv{sr_ix}.label) size(hfa{sr_ix}.time,2)]);
        
        % Get time offset
        [~, offset_ix] = min(abs(hfa{sr_ix}.time-actv{sr_ix}.time(1)));
        %!!! add check for cust_win for S+D analyses here
        %!!! can also switch to using .win_lim_s now!
        
        % Add mask of thresholded stats
        stat{sr_ix}.mask(:,offset_ix:offset_ix+size(actv{sr_ix}.time,2)-1) = actv{sr_ix}.mask;
    end
elseif any(strcmp(conditions,grp_lab))
    grp_ix = find(strcmp(grp_lab,conditions));
    % Convert to mask time series
    stat = cell(size(evnt_lab));
    for sr_ix = 1:numel(evnt_lab)
        % Prepare stat mask to match hfa for specific factor
        stat{sr_ix} = rmfield(hfa{sr_ix},{'powspctrm','freq','cumtapcnt'});
        stat{sr_ix}.mask = zeros([numel(w2{sr_ix}.label) size(hfa{sr_ix}.time,2)]);
        
%         % Get time offset
%         [~, offset_ix] = min(abs(hfa{sr_ix}.time-(w2{sr_ix}.time(1)-st.win_len/2)));
%         
        for ch_ix = 1:numel(stat{1}.label)
            % Add significant epochs for condition of interest to mask
            sig_chunks = fn_find_chunks(squeeze(w2{sr_ix}.qval(grp_ix,ch_ix,:))<st.alpha);
            sig_chunks(squeeze(w2{sr_ix}.qval(grp_ix,ch_ix,sig_chunks(:,1)))>st.alpha,:) = [];
            for sig_ix = 1:size(sig_chunks,1)
                sig_time_ix = zeros([1 2]); sig_time_err = zeros([1 2]);
                [sig_time_err(1),sig_time_ix(1)] = min(abs(stat{sr_ix}.time-w2{sr_ix}.win_lim_s(sig_chunks(sig_ix,1),1)));
                [sig_time_err(1),sig_time_ix(2)] = min(abs(stat{sr_ix}.time-w2{sr_ix}.win_lim_s(sig_chunks(sig_ix,2),2)));
                if any(sig_time_err>0.0001)
                    error('exact stat time not found in hfa!');
                end
                stat{sr_ix}.mask(ch_ix,sig_time_ix(1):sig_time_ix(2)) = 1;
%                 % If last window is cut off due to 10ms over, just go to end of trial
%                 if any(sig_chunks(sig_ix,:) > size(w2{sr_ix}.win_lim,1))
%                     sig_chunks(sig_ix,sig_chunks(sig_ix,:) > size(w2{sr_ix}.win_lim,1)) = size(w2{sr_ix}.win_lim,1);
%                 end
%                 stat{sr_ix}.mask(ch_ix,[w2{sr_ix}.win_lim(sig_chunks(sig_ix,1),1):...
%                                         w2{sr_ix}.win_lim(sig_chunks(sig_ix,2),2)]+offset_ix-1) = 1;
            end
        end
    end
end

%% Plot Results
fig_dir = [root_dir 'PRJ_Stroop/results/HFA/' SBJ '/ERPstack_' conditions '/' ...
            stat_id_s '-' stat_id_r '/' an_id_s '-' an_id_r '/'];
sig_ln_dir = [fig_dir 'sig_ch/'];
if ~exist(sig_ln_dir,'dir')
    [~] = mkdir(sig_ln_dir);
end

% Create a figure for each channel
for ch_ix = 1:2%numel(hfa{1}.label)
    sig_flag = 0;
    % Plot parameters
    fig_name = [SBJ '_' conditions '_SR_ERPstack_' hfa{1}.label{ch_ix}];
    figure('Name',fig_name,'units','normalized',...
        'outerposition',[0 0 0.6 0.8],'Visible',fig_vis);
    
    %% Get Single channel data and color limits
    ch_data = cell(size(hfa));
    ch_mean = cell([numel(hfa) numel(cond_lab)]);
    ch_var  = cell([numel(hfa) numel(cond_lab)]);
    trl_clims = NaN([2 2]);
    max_var = 0;
    for sr_ix = 1:2
        if isfield(hfa{sr_ix},'powspctrm')
            ch_data{sr_ix} = squeeze(hfa{sr_ix}.powspctrm(:,ch_ix,:,:));
        else
            ch_data{sr_ix} = NaN([numel(hfa{sr_ix}.trial) numel(hfa{sr_ix}.time{1})]);
            for t_ix = 1:numel(hfa{sr_ix}.trial)
                % if strcmp(evnt_lab{sr_ix},'S')
                %  [~,trl_end_ix] = min(abs(hfa{sr_ix}.time-trial_info.response_time(t_ix)));
                ch_data{sr_ix}(t_ix,:) = squeeze(hfa{sr_ix}.trial{t_ix}(ch_ix,:));
            end
        end
        % Get color limits
        trl_clims(sr_ix,1) = prctile(ch_data{sr_ix}(:),plt_vars.clim_perc(1));
        trl_clims(sr_ix,2) = prctile(ch_data{sr_ix}(:),plt_vars.clim_perc(2));
        
        % Condition specific mean, variance
        for cond_ix = 1:numel(cond_lab)
            ch_mean{sr_ix,cond_ix} = nanmean(ch_data{sr_ix}(cond_mat(cond_mat(:,1)==cond_ix,3),:),1);
            ch_var{sr_ix,cond_ix}  = squeeze(nanstd(ch_data{sr_ix}(cond_mat(cond_mat(:,1)==cond_ix,3),:),[],1)./...
                sqrt(sum(cond_mat(:,1)==cond_ix)))';
            max_var = max([max_var max(ch_var{sr_ix,cond_ix})]);
        end
        
    end
    trl_clims = [min(trl_clims(:,1)) max(trl_clims(:,2))];
    erp_ylims = [min([ch_mean{:}])-max_var max([ch_mean{:}])+max_var];
    
    %% Plot Data
    for sr_ix = 1:numel(evnt_lab)
        subplot('Position',plt_vars.subplot_pos{1,sr_ix});
        hold on;
        
        %% Plot Single Trials Per Condition
        imagesc(ch_data{sr_ix}(cond_mat(:,3),:));
        set(gca,'YDir','normal');
        if strcmp(evnt_lab{sr_ix},'S')
            x_tick_lab = plt_vars.plt_lim_S(1):plt_vars.x_step_sz:plt_vars.plt_lim_S(2);
            event_time = -plt_vars.plt_lim_S(1)*srate;
            for cond_ix = 1:numel(cond_lab)
                idx = cond_mat(:,1)==cond_ix;
                scat = scatter(cond_mat(idx,2)-plt_vars.plt_lim_S(1)*srate,find(idx),...
                    'MarkerFaceColor',[cond_colors{cond_ix}],'MarkerEdgeColor','k');
            end
        else
            x_tick_lab = plt_vars.plt_lim_R(1):plt_vars.x_step_sz:plt_vars.plt_lim_R(2);
            event_time = -plt_vars.plt_lim_R(1)*srate;
        end
        ylim([1 size(cond_mat,1)]);
        event_line = line([event_time event_time],ylim,...
            'LineWidth',plt_vars.evnt_width,'Color','k');
        
        % Stack Plotting parameters
        ax = gca;
%         ax.legend = plt_vars.legend;
        ax.Title.String  = [hfa{sr_ix}.label{ch_ix} ' (' elec.roi{ch_ix} '): ' evnt_lab{sr_ix} ' trials'];
        ax.XLim          = [0,size(ch_data{sr_ix},2)];
        ax.XTick         = 0:plt_vars.x_step_sz*srate:size(ch_data{sr_ix},2);
        ax.XTickLabel    = x_tick_lab;
        ax.XLabel.String = 'Time (s)';
        ax.YLabel.String = 'Trials';
%         legend([roi_lines{logical(roi_flag)}],lgd_lab{logical(roi_flag)},'Location',lgd_loc);
        cbar = colorbar;
        caxis(trl_clims);
        ax.Position = plt_vars.subplot_pos{1,sr_ix};    % adjust postiion after adding cbar
        
        %% Plot ERP/average trace
        subplot('Position',plt_vars.subplot_pos{2,sr_ix}); hold on;
        err_lines  = cell(size(cond_lab));
        main_lines = cell(size(cond_lab));
        for cond_ix = 1:numel(cond_lab)
            err_lines{cond_ix} = shadedErrorBar([], ch_mean{sr_ix,cond_ix}, ch_var{sr_ix,cond_ix},...
                {'Color',cond_colors{cond_ix},...
                'LineStyle',cond_style{cond_ix}},plt_vars.errbar_alpha);
            main_lines{cond_ix} = err_lines{cond_ix}.mainLine;
            
            % Plot mean RT (for S-locked)
            if strcmp(evnt_lab{sr_ix},'S')
                line([mean(cond_mat(cond_mat(:,1)==cond_ix,2))-plt_vars.plt_lim_S(1)*srate...
                    mean(cond_mat(cond_mat(:,1)==cond_ix,2))-plt_vars.plt_lim_S(1)*srate],erp_ylims,...
                    'Color',cond_colors{cond_ix},'LineStyle',plt_vars.evnt_style);
            end
        end
        
        % Find significant epochs
        sig_chunks = fn_find_chunks(squeeze(stat{sr_ix}.mask(ch_ix,:))<st.alpha);
        sig_chunks(squeeze(stat{sr_ix}.mask(ch_ix,sig_chunks(:,1)))==0,:) = [];
        
        % Plot Significance Shading
        if ~isempty(sig_chunks)
            sig_flag = 1;
            fprintf('%s %s (%s) -- %i SIGNIFICANT CLUSTERS FOUND...\n',...
                hfa{sr_ix}.label{ch_ix},conditions,evnt_lab{sr_ix},size(sig_chunks,1));
            for sig_ix = 1:size(sig_chunks,1)
                patch([sig_chunks(sig_ix,1) sig_chunks(sig_ix,1) sig_chunks(sig_ix,2) sig_chunks(sig_ix,2)],...
                    [erp_ylims(1) erp_ylims(2) erp_ylims(2) erp_ylims(1)],...
                    plt_vars.evnt_color,'FaceAlpha',plt_vars.errpatch_alpha);
            end
        end
        
        % Plot Event Line
        if strcmp(evnt_lab{sr_ix},'R')
            line([-plt_vars.plt_lim_R(1)*srate -plt_vars.plt_lim_R(1)*srate],erp_ylims,...
                'Color','k');%'LineWidth',plt_vars.evnt_width,
        end
        ax = gca;
        ax.YLim          = erp_ylims;
        ax.XLim          = [0,size(ch_data{sr_ix},2)];
        ax.XTick         = 0:plt_vars.x_step_sz*srate:size(ch_data{sr_ix},2);
        ax.XTickLabel    = x_tick_lab;
        ax.XLabel.String = 'Time (s)';
        ax.YLabel.String = 'HFA';
%         legend([main_lines{:}],cond_lab,'Location',eval(['plt_vars.legend_loc_' evnt_lab{sr_ix}]));
    end
    
    %% Save figure
    if save_fig
        fig_fname = [fig_dir fig_name '.' fig_ftype];
        fprintf('Saving %s\n',fig_fname);
        saveas(gcf,fig_fname);
        
        % Symbolic link for significant plots
        if sig_flag
            cd(sig_ln_dir);
            link_cmd = ['ln -s ../' fig_name '.' fig_ftype ' .'];
            system(link_cmd);
        end
    end
end

end
