function SBJ08c_HFA_summary_plot_actv(SBJ,an_id,actv_win,plt_id,save_fig,fig_vis)
% Load HFA analysis results for active ONLY (no ROI info)
%   epochs, plot a summary of those time period per electrode
% clear all; %close all;
fig_filetype = 'png';
if ischar(save_fig); save_fig = str2num(save_fig); end
if isnumeric(actv_win); actv_win = num2str(actv_win); end

%% Data Preparation
% Set up paths
addpath('/home/knight/hoycw/PRJ_Stroop/scripts/');
addpath('/home/knight/hoycw/PRJ_Stroop/scripts/utils/');
addpath('/home/knight/hoycw/Apps/fieldtrip/');
ft_defaults

%% Load Results
% Load variables
SBJ_vars_cmd = ['run /home/knight/hoycw/PRJ_Stroop/scripts/SBJ_vars/' SBJ '_vars.m'];
eval(SBJ_vars_cmd);
an_vars_cmd = ['run /home/knight/hoycw/PRJ_Stroop/scripts/an_vars/' an_id '_vars.m'];
eval(an_vars_cmd);
plt_vars_cmd = ['run /home/knight/hoycw/PRJ_Stroop/scripts/plt_vars/' plt_id '_vars.m'];
eval(plt_vars_cmd);

% Load data
actv_filename = strcat(SBJ_vars.dirs.proc,SBJ,'_actv_ROI_',an_id,'_mn',actv_win,'.mat');
load(actv_filename);

%% Process parameters
% Get event timing
if strcmp(an_id(1:5),'HGm_S')
    event = 'stim';
    % Compute mean RT
    load(strcat(SBJ_vars.dirs.events,SBJ,'_trial_info_final.mat'),'trial_info');
    RT_mean = mean(trial_info.response_time); % converts sec to ms
    % Add in the baseline offset to plot correctly
    event_time = RT_mean-plt_vars.plt_lim(1);
elseif strcmp(an_id(1:5),'HGm_R')
    event = 'resp';
    event_time = 0;
end

%% Aggregate results per channel
actv_epochs = {};
actv_sign   = {};
actv_code   = zeros([1, numel(hfa.label)]);
for ch_ix = 1:numel(hfa.label)
    % Check if active, get epochs
    if ~isempty(strmatch(hfa.label{ch_ix},actv_ch,'exact'))
        % Find significant epoch indices
        actv_ix = strmatch(hfa.label{ch_ix},actv_ch,'exact');
        actv_epochs{ch_ix} = actv_ch_epochs{actv_ix};
        actv_ep_sign = NaN([1 size(actv_ch_epochs{actv_ix},1)]);
        sig_chunk_ix = NaN([1 2]);
        for ep_ix = 1:size(actv_ch_epochs{actv_ix},1)
            sig_chunk_ix = [find(hfa.time==actv_ch_epochs{actv_ix}(ep_ix,1))...
                            find(hfa.time==actv_ch_epochs{actv_ix}(ep_ix,2))];
            % Find sign of (de)activation
            if 0<=squeeze(mean(mean(hfa.powspctrm(:,ch_ix,1,sig_chunk_ix(1):sig_chunk_ix(2)),1),4))
                actv_ep_sign(ep_ix) = 1;
            else
                actv_ep_sign(ep_ix) = -1;
            end
        end
        actv_sign{ch_ix} = actv_ep_sign;
        % Code significance and sign: 0=none, -1=deactive, 1=active, 2=both
        if any(actv_ep_sign==1) && any(actv_ep_sign==-1)
            actv_code(ch_ix) = 2;
        elseif any(actv_ep_sign==1)
            actv_code(ch_ix) = 1;
        elseif any(actv_ep_sign==-1)
            actv_code(ch_ix) = -1;
        end
    else
        actv_epochs{ch_ix} = [];
        actv_sign{ch_ix}   = 0;
    end
end

%% Plot Results
% Creat and format the plot
fig_name = [SBJ '_HFA_summary_actv_' event];
figure('Name',fig_name,'units','normalized',...
        'outerposition',[0 0 1 1],'Visible',fig_vis);
hold on;

% Plot labels
ax = gca;
ax.XLabel.String   = 'Time (s)';
ax.XLabel.FontSize = 14;
ax.XLim   = plt_vars.plt_lim;
ax.XTick  = plt_vars.plt_lim(1):plt_vars.x_step_sz:plt_vars.plt_lim(2);

ax.YLabel.String   = 'Channels';
ax.YLabel.FontSize = 14;
ax.YLim  = [0 numel(hfa.label)+1];
ax.YTick = [1:numel(hfa.label)];
ylabels  = hfa.label;
for lab_ix = 2:2:numel(ylabels) % To space out the labels and make them readable
    ylabels{lab_ix} = ['     ' ylabels{lab_ix}];
end
ax.YTickLabel         = hfa.label;
ax.YTickLabelRotation = 45;
% plot_info.legend_loc = plt_vars.legend_loc;

% Plot Summary info in title
perc_actv      = numel(find(actv_code~=0))/numel(hfa.label);
perc_actv_pos  = numel(find(actv_code==1))/numel(hfa.label);
perc_actv_neg  = numel(find(actv_code==-1))/numel(hfa.label);
perc_actv_both = numel(find(actv_code==2))/numel(hfa.label);

ax.Title.String = sprintf('Active: %.2f (+=%.2f;-=%.2f;+-=%.2f)',...
    perc_actv,perc_actv_pos,perc_actv_neg,perc_actv_both);
ax.Title.FontSize = 16;

% Plot Event
if strcmp(an_id(1:5),'HGm_S')
    event_time = RT_mean;
elseif strcmp(an_id(1:5),'HGm_R')
    event_time = 0;
end
line([event_time event_time], [0 numel(hfa.label)+1],...
    'LineWidth',plt_vars.evnt_width, 'Color',plt_vars.evnt_color,...
    'LineStyle',plt_vars.evnt_style);

% Plot Summary Epochs
for ch_ix = 1:numel(hfa.label)
    % Plot active epochs
    for ep_ix = 1:size(actv_epochs{ch_ix},1)
        % [top_L top_R bot_R bot_L]
        x = [actv_epochs{ch_ix}(ep_ix,1) actv_epochs{ch_ix}(ep_ix,2) ...
             actv_epochs{ch_ix}(ep_ix,2) actv_epochs{ch_ix}(ep_ix,1)];
        y = [ch_ix+plt_vars.actv_height/2 ch_ix+plt_vars.actv_height/2 ...
             ch_ix-plt_vars.actv_height/2 ch_ix-plt_vars.actv_height/2];
        if actv_sign{ch_ix}(ep_ix)==1
            patch(x,y,plt_vars.actv_color{1},'FaceAlpha',plt_vars.actv_alpha);
        else
            patch(x,y,plt_vars.actv_color{2},'FaceAlpha',plt_vars.actv_alpha);
        end
    end
end

%% Save figure
if save_fig
    fig_dir = ['/home/knight/hoycw/PRJ_Stroop/results/HFA/' SBJ '/summary_actv/'...
        an_id '_mn' actv_win '/'];
    if ~exist(fig_dir,'dir')
        mkdir(fig_dir);
    end
    
    fig_filename = [fig_dir fig_name '.' fig_filetype];
    fprintf('Saving %s\n',fig_filename);
    saveas(gcf,fig_filename);
    %eval(['export_fig ' fig_filename]);
end

end
