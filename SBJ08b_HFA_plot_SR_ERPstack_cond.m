function SBJ08b_HFA_plot_SR_ERPstack_cond(SBJ,conditions,an_id_s,an_id_r,...
                                        plt_id,save_fig,fig_vis,fig_ftype,varargin)
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

if ~isempty(varargin)
    for v = 1:2:numel(varargin)
        if strcmp(varargin{v},'actv_win')
            actv_win = varargin{v+1};
            if isnumeric(actv_win); actv_win = num2str(actv_win); end
        else
            error(['Unknown varargin ' num2str(v) ': ' varargin{v}]);
        end
    end
end

%% Load Results
SBJ_vars_cmd = ['run ' root_dir 'PRJ_Stroop/scripts/SBJ_vars/' SBJ '_vars.m'];
eval(SBJ_vars_cmd);
plt_vars_cmd = ['run ' root_dir 'PRJ_Stroop/scripts/plt_vars/' plt_id '_vars.m'];
eval(plt_vars_cmd);

% Load RTs
load(strcat(SBJ_vars.dirs.events,SBJ,'_trial_info_final.mat'),'trial_info');

% Load data
hfa_filename1 = strcat(SBJ_vars.dirs.proc,SBJ,'_ROI_',an_id_s,'.mat');
hfa_filename2 = strcat(SBJ_vars.dirs.proc,SBJ,'_ROI_',an_id_r,'.mat');
tmp = load(hfa_filename1,'hfa'); hfa{1} = tmp.hfa;
tmp = load(hfa_filename2,'hfa'); hfa{2} = tmp.hfa;
clear tmp

% Load elec
elec_fname = [SBJ_vars.dirs.recon SBJ '_elec_main_ft_pat_Dx_full.mat'];
load(elec_fname);
% Sort elecs by stat labels
cfgs = []; cfgs.channel = hfa{1}.label;
elec = fn_select_elec(cfgs,elec);
elec.roi = fn_atlas2roi_labels(elec.atlas_lab,elec.atlas_id,'ROI');

% Load stats
if strcmp(conditions,'actv')
    stat_fname1 = strcat(SBJ_vars.dirs.proc,SBJ,'_ROI_',an_id_s,'_actv_mn',actv_win,'.mat');
    stat_fname2 = strcat(SBJ_vars.dirs.proc,SBJ,'_ROI_',an_id_r,'_actv_mn',actv_win,'.mat');
    tmp = load(stat_fname1,'actv_ch'); actv_ch{1} = tmp.actv_ch;
    tmp = load(stat_fname2,'actv_ch'); actv_ch{2} = tmp.actv_ch;
    tmp = load(stat_fname1,'actv_ch_epochs'); actv_ch_epochs{1} = tmp.actv_ch_epochs;
    tmp = load(stat_fname2,'actv_ch_epochs'); actv_ch_epochs{2} = tmp.actv_ch_epochs;
    clear tmp
end

%% Quality checks
% Check channel match between analyses
if ~isempty(setdiff(hfa{1}.label,hfa{2}.label))
    error('ERROR: channels do not match between the two analyses!');
end

% If singe trial format, check all time axes are the same
time_cells = [iscell(hfa{1}.time), iscell(hfa{2}.time)];
if any(time_cells)
    sample_rate = zeros([1 2]);
    for sr_ix = find(time_cells)
        for t_ix = 2:numel(hfa{sr_ix}.time)
            if ~all(hfa{sr_ix}.time{1}==hfa{sr_ix}.time{t_ix})
                error('Time axes arent the same for all trials!');
            end
        end
        sample_rate(sr_ix) = (numel(hfa{sr_ix}.time{1})-1)/(hfa{sr_ix}.time{1}(end)-hfa{sr_ix}.time{1}(1));
    end
    if diff(sample_rate) > 1
        error('sample rates dont match between analyses!');
    else
        sample_rate = sample_rate(1);
    end
else
    sample_rate = (numel(hfa{1}.time)-1)/(hfa{1}.time(end)-hfa{1}.time(1));
end

%% Prep Data
evnt_lab = {'S', 'R'};

% Trim data to plotting epoch
cfg_trim = [];
cfg_trim.latency = plt_vars.plt_lim_S;
hfa{1} = ft_selectdata(cfg_trim,hfa{1});
cfg_trim.latency = plt_vars.plt_lim_R;
hfa{2} = ft_selectdata(cfg_trim,hfa{2});

% Compile cond_type, RT, trial_n
[cond_lab, cond_colors, cond_style] = fn_condition_label_styles(conditions);
if strcmp(conditions,'actv')
    [sort_rts, sort_rt_idx] = sort(round(sample_rate*trial_info.response_time));
    % Fake cond_mat with only sorted RTs
    cond_mat = [ones([numel(trial_info.trial_n) 1]) sort_rts sort_rt_idx];
else
    cond_idx = false([length(cond_lab) length(trial_info.trial_n)]);
    for cond_ix = 1:length(cond_lab)
        % Get binary condition index
        cond_idx(cond_ix,:) = logical(fn_condition_index(cond_lab{cond_ix},...
            trial_info.condition_n));
    end
    
    % Create matrix of [cond_idx, rt, trial_ix]
    [cond_mat,~] = find(cond_idx);
    cond_mat = horzcat(cond_mat,round(sample_rate*trial_info.response_time),[1:numel(cond_mat)]');
    cond_mat = sortrows(cond_mat,[1 2]);
end

%% Plot Results
fig_dir = [root_dir 'PRJ_Stroop/results/HFA/' SBJ '/ERPstack_' conditions '/' an_id_s '-' an_id_r '/'];
if ~exist(fig_dir,'dir')
    [~] = mkdir(fig_dir);
end

% Create a figure for each channel
for ch_ix = 1:numel(hfa{1}.label)
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
                ch_data{sr_ix}(t_ix,:) = squeeze(hfa{sr_ix}.trial{t_ix}(ch_ix,:));
            end
        end
        % Get color limits
        trl_clims(sr_ix,1) = prctile(ch_data{sr_ix}(:),plt_vars.clim_perc(1));
        trl_clims(sr_ix,2) = prctile(ch_data{sr_ix}(:),plt_vars.clim_perc(2));
        
        % Condition specific mean, variance
        for cond_ix = 1:numel(cond_lab)
            ch_mean{sr_ix,cond_ix} = mean(ch_data{sr_ix}(cond_mat(:,1)==cond_ix,:),1);
            ch_var{sr_ix,cond_ix}  = squeeze(std(ch_data{sr_ix}(cond_mat(:,1)==cond_ix,:),[],1)./...
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
        ax     = gca;
        
        %% Plot Single Trials Per Condition
        imagesc(ch_data{sr_ix}(cond_mat(:,3),:));
        set(gca,'YDir','normal');
        if strcmp(evnt_lab{sr_ix},'S')
            x_tick_lab = plt_vars.plt_lim_S(1):plt_vars.x_step_sz:plt_vars.plt_lim_S(2);
            event_time = -plt_vars.plt_lim_S(1)*sample_rate;
            for cond_ix = 1:numel(cond_lab)
                idx = cond_mat(:,1)==cond_ix;
                scat = scatter(cond_mat(idx,2)-plt_vars.plt_lim_S(1)*sample_rate,find(idx),...
                    'MarkerFaceColor',[cond_colors{cond_ix}],'MarkerEdgeColor','k');
            end
        else
            x_tick_lab = plt_vars.plt_lim_R(1):plt_vars.x_step_sz:plt_vars.plt_lim_R(2);
            event_time = -plt_vars.plt_lim_R(1)*sample_rate;
        end
        ylim([1 size(cond_mat,1)]);
        event_line = line([event_time event_time],ylim,...
            'LineWidth',plt_vars.evnt_width,'Color','k');
        
        % Stack Plotting parameters
        ax = gca;
%         ax.legend = plt_vars.legend;
        ax.Title.String  = [hfa{sr_ix}.label{ch_ix} ' (' elec.roi{ch_ix} '): ' evnt_lab{sr_ix} ' trials'];
        ax.XLim          = [0,size(ch_data{sr_ix},2)];
        ax.XTick         = 0:plt_vars.x_step_sz*sample_rate:size(ch_data{sr_ix},2);
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
                line([mean(cond_mat(cond_mat(:,1)==cond_ix,2))-plt_vars.plt_lim_S(1)*sample_rate...
                    mean(cond_mat(cond_mat(:,1)==cond_ix,2))-plt_vars.plt_lim_S(1)*sample_rate],erp_ylims,...
                    'Color',cond_colors{cond_ix},'LineStyle',plt_vars.evnt_style);
            end
            
            % Plot significance shading
            if strcmp(conditions,'actv')
                % Find significant time periods
                if any(strcmp(hfa{sr_ix}.label{ch_ix},actv_ch{sr_ix}))
                    % Find significant epoch indices
                    actv_ch_ix = strcmp(hfa{sr_ix}.label{ch_ix},actv_ch{sr_ix});
                    sig_chunks = NaN(size(actv_ch_epochs{sr_ix}{actv_ch_ix}));
                    for win_ix = 1:size(actv_ch_epochs{sr_ix}{actv_ch_ix},1)
                        sig_chunks(win_ix,1) = find(hfa{sr_ix}.time==actv_ch_epochs{sr_ix}{actv_ch_ix}(win_ix,1));
                        sig_chunks(win_ix,2) = find(hfa{sr_ix}.time==actv_ch_epochs{sr_ix}{actv_ch_ix}(win_ix,2));
                    end
                    ylims = ylim;
                    for sig_ix = 1:size(sig_chunks,1)
                        patch([sig_chunks(sig_ix,1) sig_chunks(sig_ix,1) sig_chunks(sig_ix,2) sig_chunks(sig_ix,2)],...
                              [ylims(1) ylims(2) ylims(2) ylims(1)],...
                              plt_vars.evnt_color,'FaceAlpha',plt_vars.errpatch_alpha);
                    end
                end
            end
        end
        if strcmp(evnt_lab{sr_ix},'R')
            line([-plt_vars.plt_lim_R(1)*sample_rate -plt_vars.plt_lim_R(1)*sample_rate],erp_ylims,...
                'Color','k');%'LineWidth',plt_vars.evnt_width,
        end
        ax = gca;
        ax.YLim          = erp_ylims;
        ax.XLim          = [0,size(ch_data{sr_ix},2)];
        ax.XTick         = 0:plt_vars.x_step_sz*sample_rate:size(ch_data{sr_ix},2);
        ax.XTickLabel    = x_tick_lab;
        ax.XLabel.String = 'Time (s)';
        ax.YLabel.String = 'HFA';
        legend([main_lines{:}],cond_lab,'Location',eval(['plt_vars.legend_loc_' evnt_lab{sr_ix}]));
    end
    
    %% Save figure
    if save_fig
        fig_fname = [fig_dir fig_name '.' fig_ftype];
        fprintf('Saving %s\n',fig_fname);
        saveas(gcf,fig_fname);
        %eval(['export_fig ' fig_filename]);
    end
end

end
