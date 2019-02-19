% function [ output_args ] = untitled3(SBJ,elec_lab,conditions,an_id,stat_id,plt_id,save_fig,fig_vis)
%Plot HFA time series averaged by RT splits
%   Detailed explanation goes here

SBJ = 'IR35';
conditions = 'CI';
an_id = 'HGm_S_zbtS_trl2to151_sm0_wn100_stat15';
elec_lab = {'LAC2-3'};
pipeline_id = 'main_ft';
save_fig = 0;
fig_vis = 'on';
plt_id = 'ts_errbr_evnt';

RT_split = 4;    % 2 = median, 3 = thirds, 4 = quartiles

fig_filetype = 'png';
if ischar(save_fig); save_fig = str2num(save_fig); end

%% Data Preparation
% Set up paths
if exist('/home/knight/hoycw/','dir');root_dir='/home/knight/hoycw/';ft_dir=[root_dir 'Apps/fieldtrip/'];
else root_dir='/Volumes/hoycw_clust/';ft_dir='/Users/colinhoy/Code/Apps/fieldtrip/';end
addpath([root_dir 'PRJ_Stroop/scripts/']);
addpath([root_dir 'PRJ_Stroop/scripts/utils/']);
addpath(ft_dir);
ft_defaults

%% Load Results
SBJ_vars_cmd = ['run ' root_dir 'PRJ_Stroop/scripts/SBJ_vars/' SBJ '_vars.m'];
eval(SBJ_vars_cmd);
plt_vars_cmd = ['run ' root_dir 'PRJ_Stroop/scripts/plt_vars/' plt_id '_vars.m'];
eval(plt_vars_cmd);

[cond_lab, cond_colors, cond_style] = fn_condition_label_styles(conditions);
event_lab = {'stim', 'resp'};
% Load RTs
load(strcat(SBJ_vars.dirs.events,SBJ,'_trial_info_final.mat'),'trial_info');

elec_fname = [SBJ_vars.dirs.recon SBJ '_elec_' pipeline_id '_pat_Dx_tis.mat'];
load(elec_fname);
% roi = fn_atlas2roi_labels(elec.atlas_label(strcmp(elec.label,elec_lab)),'Dx','ROI');
% roi = roi{1};

stats_fname = strcat(SBJ_vars.dirs.proc,SBJ,'_',conditions,'_ROI_',an_id,'.mat');
load(stats_fname);

sample_rate = (numel(hfa{1}.time)-1)/(hfa{1}.time(end)-hfa{1}.time(1));

%% Prep Data
% Trim data to plotting epoch
cfg_trim = [];
cfg_trim.channel = elec_lab;
stat = ft_selectdata(cfg_trim,stat);
% cfg_trim.latency = plt_vars.plt_lim_S;
hfa{1} = ft_selectdata(cfg_trim,hfa{1});
hfa{2} = ft_selectdata(cfg_trim,hfa{2});
% cfg_trim.latency = plt_vars.plt_lim_R;
% hfa{2,1} = ft_selectdata(cfg_trim,hfa{2,1});

% Compute mean RT per condition
RTs = round(1000*trial_info.response_time); % converts sec to ms
cond_idx = zeros(size(RTs));
rt_cond = cell(size(cond_lab));
for cond_ix = 1:numel(cond_lab)
    cond_idx(fn_condition_index([cond_lab{cond_ix}], trial_info.condition_n)==1) = cond_ix;
    rt_cond{cond_ix} = RTs(fn_condition_index([cond_lab{cond_ix}], trial_info.condition_n)==1);
%     RT_means{cond_ix} = mean(RTs(fn_condition_index([cond_lab{cond_ix}], trial_info.condition_n)==1));
%     % Add in the baseline offset to plot correctly
%     RT_means{cond_ix} = RT_means{cond_ix}-plt_vars.plt_lim_S(1)*1000;
end

%% Split trials by RT
% split_cut = zeros([1 RT_split+1]);
split_idx = cell(size(cond_lab));
split_cnt = zeros([numel(cond_lab) RT_split+1]);
split_avg = zeros([numel(cond_lab) RT_split]);

split_cut = quantile(RTs,0:1/RT_split:1);
split_cut(end) = split_cut(end)+1;  % make sure longest RT gets included in last bin (histc will put an equal value in singleton bin)
[full_cnt, full_idx] = histc(RTs,split_cut);
if full_cnt(end)~=0
    error('histc messed up and put something below the lowest RT!');
end
for cond_ix = 1:numel(cond_lab)
    [split_cnt(cond_ix,:), split_idx{cond_ix}] = histc(rt_cond{cond_ix},split_cut);
    if split_cnt(cond_ix,end)~=0
        error('histc messed up and put something below the lowest RT!');
    end
    for split_ix = 1:RT_split
        split_avg(cond_ix,split_ix) = mean(rt_cond{cond_ix}(split_idx{cond_ix}==split_ix));
    end
end
split_cnt(:,end) = [];

%% Plot Results
fig_dir = [root_dir 'PRJ_Stroop/results/HFA/' SBJ '/' conditions '/RT_split' num2str(RT_split)];
if ~exist(fig_dir,'dir')
    mkdir(fig_dir);
end

% Create a figure for each channel
ch_ix = 1; % for loop here
% Plot parameters
fig_name = strcat(SBJ,'_',conditions,'_',elec_lab{ch_ix});%roi,'_',
figure('Name',fig_name,'units','normalized',...
    'outerposition',[0 0 0.5 1],'Visible',fig_vis); %twice as wide for the double plot
plot_info.fig        = gcf;
plot_info.x_step     = plt_vars.x_step_sz*sample_rate;
plot_info.legend_loc = plt_vars.legend_loc;
plot_info.sig_alpha  = plt_vars.sig_alpha;
plot_info.sig_color  = plt_vars.sig_color;
% Condition plotting params
cond_info.style      = cond_style;
cond_info.color      = cond_colors;
cond_info.alpha      = repmat(plt_vars.errbar_alpha,[1 numel(cond_lab)]);

for split_ix = 1:RT_split
    subplot(RT_split,1,split_ix);
    plot_info.ax     = gca;
    plot_info.title  = ['spilt ix = ' num2str(split_ix) ', rt ratio = ' num2str(split_cnt(1,split_ix)/split_cnt(2,split_ix))];%[event_lab{set_ix} '-Locked'];
    plot_info.legend = 1;
    plot_info.x_lab  = hfa{1}.time(1):plt_vars.x_step_sz:hfa{1}.time(end);
    % Stimulus plotting params
    event_info.name  = {'stim', cond_lab{:}};
    event_info.color = {[0 0 0], cond_colors{:}};
    event_info.width = repmat(plt_vars.evnt_width,[1 numel(event_info.name)]);
    event_info.style = repmat({plt_vars.evnt_style},[1 numel(event_info.name)]);
    event_info.style{1} = '-';
    event_info.time  = [find(hfa{1}.time == 0), split_avg(:,split_ix)'];
    
    % Compute means and variance
    means = NaN([numel(cond_lab) size(hfa{1}.powspctrm,4)]);
    var   = NaN([numel(cond_lab) size(hfa{1}.powspctrm,4)]);
    for cond_ix = 1:numel(cond_lab)
        means(cond_ix,:) = squeeze(mean(hfa{cond_ix}.powspctrm(split_idx{cond_ix}==split_ix,1,:,:),1));
        var(cond_ix,:) = squeeze(std(hfa{cond_ix}.powspctrm(split_idx{cond_ix}==split_ix,1,:,:),[],1)./...
            sqrt(size(hfa{cond_ix}.powspctrm,1)))';
        % add trial counts to condition labels
        cond_info.name{cond_ix} = [cond_lab{cond_ix} ' (' num2str(split_cnt(cond_ix,split_ix)) ')'];
    end
    
    fn_plot_ts_error_bar(plot_info,means,var,event_info,cond_info);
    
    % Plot RT histogram underneath
%     histogram(rt_cond{cond_ix}(split_idx{cond_ix}==split_ix),5);
end

%% Save figure
% if save_fig
%     fig_filename = [fig_dir fig_name '.' fig_filetype];
%     fprintf('Saving %s\n',fig_filename,'svg');
%     saveas(gcf,fig_filename);
%     %eval(['export_fig ' fig_filename]);
% end

