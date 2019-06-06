function SBJ08b_HFA_plot_ERPstack_cond(SBJ,conditions,an_id,stat_id,...
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
[cond_lab, cond_colors, cond_style] = fn_condition_label_styles(conditions);

% Load RTs
load(strcat(SBJ_vars.dirs.events,SBJ,'_trial_info_final.mat'),'trial_info');
srate = trial_info.sample_rate;

% Load data
hfa_fname = strcat(SBJ_vars.dirs.proc,SBJ,'_ROI_',an_id,'.mat');
load(hfa_fname);
% Load Stats
if strcmp(conditions,'actv')
    stat_fname = strcat(SBJ_vars.dirs.proc,SBJ,'_ROI_',an_id,'_',stat_id,'.mat');
    load(stat_fname);
    if ~all(strcmp(hfa.label,actv.label)); error('label mismatch'); end
elseif any(strcmp(conditions,{'CNI','pCNI','PC'}))
    stat_fname = [SBJ_vars.dirs.proc SBJ '_smANOVA_ROI_' stat_id '_' an_id '.mat'];
    load(stat_fname);
    if ~all(strcmp(hfa.label,w2.label)); error('label mismatch'); end
elseif strcmp(conditions,'RT')
    error('RT not yet implemented');
else
    error(['unknown stat_id:' stat_id_s ' and ' stat_id_r]);
end

% Check compatibility
[grp_lab, ~, ~] = fn_group_label_styles(st.model_lab);
if ~strcmp(an.evnt_lab,st.evnt_lab); error('evnt_lab mismatch'); end

% Load elec
elec_fname = [SBJ_vars.dirs.recon SBJ '_elec_main_ft_pat_Dx_full.mat'];
load(elec_fname);
% Sort elecs by stat labels
cfgs = []; cfgs.channel = hfa.label;
elec = fn_select_elec(cfgs,elec);
elec.roi = fn_atlas2roi_labels(elec.atlas_lab,elec.atlas_id,'ROI');

% %% Quality checks
% % Check channel match between analyses
% if ~isempty(setdiff(hfa{1}.label,hfa{2}.label)); error('an label mismatch'); end
% if strcmp(conditions,'RT');             lab_check = stat{1}.label;
% elseif strcmp(conditions,'actv');       lab_check = actv{1}.label;
% elseif any(strcmp(conditions,grp_lab)); lab_check = w2{1}.label; end
% if any(~strcmp(hfa{1}.label,lab_check)); error('label mismatch'); end
% 
% % If singe trial format, check all time axes are the same
% time_cells = [iscell(hfa{1}.time), iscell(hfa{2}.time)];
% if any(time_cells)
%     for sr_ix = find(time_cells)
%         for t_ix = 2:numel(hfa.time)
%             if ~all(hfa.time{1}==hfa.time{t_ix})
%                 error('Time axes aren''t the same for all trials!');
%             end
%         end
%     end
% else
%     % Confirm hfa starts before and ends after stat
%     for sr_ix = 1:numel(evnt_lab)
%         if strcmp(conditions,'RT')
%             if ~all(hfa.time==stat.time)
%                 error('Time axes arent the same for hfa and RT stat!');
%             end
%         elseif strcmp(conditions,'actv')
%             % Allow for some imprecision (error if off by at least 1 sample)
%             if actv.time(1)-hfa.time(1)<=-1/srate || ...
%                     actv.time(end)-hfa.time(end)>=-1/srate
%                 error('Time axes arent the same for hfa and actv!');
%             end
%         elseif any(strcmp(conditions,grp_lab))
%             % Allow for some imprecision (error if off by at least 1 sample)
%             if (w2.time(1)-st.win_len/2)-hfa.time(1)<=-1/srate || ...
%                     (w2.time(end)+st.win_len/2)-hfa.time(end)>=-1/srate
%                 error('Time axes arent the same for hfa and ANOVA w2!');
%             end
%         end
%     end
% end

%% Prep Data
% Trim data to plotting epoch
cfg_trim = [];
cfg_trim.latency = eval(['plt_vars.plt_lim_' st.evnt_lab ';']);
hfa = ft_selectdata(cfg_trim,hfa);

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
    error('no RT stats');
%     % Trim data to plotting epoch
%     cfg_trim = [];
%     cfg_trim.latency = plt_vars.plt_lim_S;
%     stat{1} = ft_selectdata(cfg_trim,stat{1});
%     cfg_trim.latency = plt_vars.plt_lim_R;
%     stat{2} = ft_selectdata(cfg_trim,stat{2});
elseif strcmp(conditions,'actv')
    % Prepare stat mask to match hfa for specific factor
    stat = rmfield(hfa,{'powspctrm','freq','cumtapcnt'});
    stat.mask = zeros([numel(actv.label) size(hfa.time,2)]);
    
    % Get time offset
    [~, offset_ix] = min(abs(hfa.time-actv.time(1)));
    %!!! add check for cust_win for S+D analyses here
    %!!! can also switch to using .win_lim_s now!
    
    % Add mask of thresholded stats
    stat.mask(:,offset_ix:offset_ix+size(actv.time,2)-1) = actv.mask;
elseif any(strcmp(conditions,grp_lab))
    grp_ix = find(strcmp(grp_lab,conditions));
    % Prepare stat mask to match hfa for specific factor
    stat = rmfield(hfa,{'powspctrm','freq','cumtapcnt'});
    stat.mask = zeros([numel(w2.label) size(hfa.time,2)]);
    
    % Get time offset
    [~, offset_ix] = min(abs(hfa.time-(w2.time(1)-st.win_len/2)));
    
    for ch_ix = 1:numel(stat.label)
        % Add significant epochs for condition of interest to mask
        sig_chunks = fn_find_chunks(squeeze(w2.qval(grp_ix,ch_ix,:))<st.alpha);
        sig_chunks(squeeze(w2.qval(grp_ix,ch_ix,sig_chunks(:,1)))>st.alpha,:) = [];
        for sig_ix = 1:size(sig_chunks,1)
            % If last window is cut off due to 10ms over, just go to end of trial
            if any(sig_chunks(sig_ix,:) > size(w2.win_lim,1))
                sig_chunks(sig_ix,sig_chunks(sig_ix,:) > size(w2.win_lim,1)) = size(w2.win_lim,1);
            end
            stat.mask(ch_ix,[w2.win_lim(sig_chunks(sig_ix,1),1):w2.win_lim(sig_chunks(sig_ix,2),2)]+offset_ix-1) = 1;
        end
    end
end

%% Plot Results
fig_dir = [root_dir 'PRJ_Stroop/results/HFA/' SBJ '/ERPstack_' conditions '/' ...
            stat_id '/' an_id '/'];
sig_ln_dir = [fig_dir 'sig_ch/'];
if ~exist(sig_ln_dir,'dir')
    [~] = mkdir(sig_ln_dir);
end

% Create a figure for each channel
for ch_ix = 1:numel(hfa.label)
    sig_flag = 0;
    % Plot parameters
    fig_name = [SBJ '_' conditions '_SR_ERPstack_' hfa.label{ch_ix}];
    figure('Name',fig_name,'units','normalized',...
        'outerposition',[0 0 0.6 0.8],'Visible',fig_vis);
    
    %% Get Single channel data and color limits
    ch_mean = cell([numel(hfa) numel(cond_lab)]);
    ch_var  = cell([numel(hfa) numel(cond_lab)]);
    trl_clims = NaN([1 2]);
    max_var = 0;
    if isfield(hfa,'powspctrm')
        ch_data = squeeze(hfa.powspctrm(:,ch_ix,:,:));
    else
        ch_data = NaN([numel(hfa.trial) numel(hfa.time{1})]);
        for t_ix = 1:numel(hfa.trial)
            ch_data(t_ix,:) = squeeze(hfa.trial{t_ix}(ch_ix,:));
        end
    end
    % Get color limits
    trl_clims(1) = prctile(ch_data(:),plt_vars.clim_perc(1));
    trl_clims(2) = prctile(ch_data(:),plt_vars.clim_perc(2));
    
    % Condition specific mean, variance
    for cond_ix = 1:numel(cond_lab)
        ch_mean{cond_ix} = mean(ch_data(cond_mat(cond_mat(:,1)==cond_ix,3),:),1);
        ch_var{cond_ix}  = squeeze(std(ch_data(cond_mat(cond_mat(:,1)==cond_ix,3),:),[],1)./...
            sqrt(sum(cond_mat(:,1)==cond_ix)))';
        max_var = max([max_var max(ch_var{cond_ix})]);
    end
    
    erp_ylims = [min([ch_mean{:}])-max_var max([ch_mean{:}])+max_var];
    
    %% Plot Data
    subplot('Position',plt_vars.subplot_pos{1});
    hold on;
    
    %% Plot Single Trials Per Condition
    imagesc(ch_data(cond_mat(:,3),:));
    set(gca,'YDir','normal');
    if strcmp(st.evnt_lab,'S')
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
    ax.Title.String  = [hfa.label{ch_ix} ' (' elec.roi{ch_ix} '): ' st.evnt_lab ' trials'];
    ax.XLim          = [0,size(ch_data,2)];
    ax.XTick         = 0:plt_vars.x_step_sz*srate:size(ch_data,2);
    ax.XTickLabel    = x_tick_lab;
    ax.XLabel.String = 'Time (s)';
    ax.YLabel.String = 'Trials';
    %         legend([roi_lines{logical(roi_flag)}],lgd_lab{logical(roi_flag)},'Location',lgd_loc);
    cbar = colorbar;
    caxis(trl_clims);
    ax.Position = plt_vars.subplot_pos{1};    % adjust postiion after adding cbar
    
    %% Plot ERP/average trace
    subplot('Position',plt_vars.subplot_pos{2}); hold on;
    err_lines  = cell(size(cond_lab));
    main_lines = cell(size(cond_lab));
    for cond_ix = 1:numel(cond_lab)
        err_lines{cond_ix} = shadedErrorBar([], ch_mean{cond_ix}, ch_var{cond_ix},...
            {'Color',cond_colors{cond_ix},...
            'LineStyle',cond_style{cond_ix}},plt_vars.errbar_alpha);
        main_lines{cond_ix} = err_lines{cond_ix}.mainLine;
        
        % Plot mean RT (for S-locked)
        if strcmp(st.evnt_lab,'S')
            line([mean(cond_mat(cond_mat(:,1)==cond_ix,2))-plt_vars.plt_lim_S(1)*srate...
                mean(cond_mat(cond_mat(:,1)==cond_ix,2))-plt_vars.plt_lim_S(1)*srate],erp_ylims,...
                'Color',cond_colors{cond_ix},'LineStyle',plt_vars.evnt_style);
        end
    end
    
    % Find significant epochs
    sig_chunks = fn_find_chunks(squeeze(stat.mask(ch_ix,:))<st.alpha);
    sig_chunks(squeeze(stat.mask(ch_ix,sig_chunks(:,1)))==0,:) = [];
    
    % Plot Significance Shading
    if ~isempty(sig_chunks)
        sig_flag = 1;
        fprintf('%s %s (%s) -- %i SIGNIFICANT CLUSTERS FOUND...\n',...
            hfa.label{ch_ix},conditions,st.evnt_lab,size(sig_chunks,1));
        for sig_ix = 1:size(sig_chunks,1)
            patch([sig_chunks(sig_ix,1) sig_chunks(sig_ix,1) sig_chunks(sig_ix,2) sig_chunks(sig_ix,2)],...
                [erp_ylims(1) erp_ylims(2) erp_ylims(2) erp_ylims(1)],...
                plt_vars.evnt_color,'FaceAlpha',plt_vars.errpatch_alpha);
        end
    end
    
    % Plot Event Line
    if strcmp(st.evnt_lab,'R')
        line([-plt_vars.plt_lim_R(1)*srate -plt_vars.plt_lim_R(1)*srate],erp_ylims,...
            'Color','k');%'LineWidth',plt_vars.evnt_width,
    end
    ax = gca;
    ax.YLim          = erp_ylims;
    ax.XLim          = [0,size(ch_data,2)];
    ax.XTick         = 0:plt_vars.x_step_sz*srate:size(ch_data,2);
    ax.XTickLabel    = x_tick_lab;
    ax.XLabel.String = 'Time (s)';
    ax.YLabel.String = 'HFA';
    legend([main_lines{:}],cond_lab,'Location',eval(['plt_vars.legend_loc_' st.evnt_lab]));

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
