function SBJ08c_HFA_summary_plot_actv_cond(SBJ,conditions,an_id,actv_win,plt_id,save_fig,fig_vis)
% Load HFA analysis results for active and condition-differentiating
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
load(strcat(SBJ_vars.dirs.proc,SBJ,'_',conditions,'_ROI_',an_id,'.mat'));
actv_filename = strcat(SBJ_vars.dirs.proc,SBJ,'_actv_ROI_',an_id,'_mn',actv_win,'.mat');
load(actv_filename,'actv_ch','actv_ch_epochs');
tmp = load(actv_filename,'hfa'); hfa_actv = tmp.hfa;

%% Load ROI and GM/WM info
einfo_filename = [SBJ_vars.dirs.proc SBJ '_elec_ROIs.tsv'];
% roi_filename = strcat(SBJ_vars.dirs.proc,SBJ,'_elec_ROIs.csv');
% einfo = textread(roi_filename,'%s','delimiter',',');
% einfo = reshape(einfo,[3, numel(stat.label)]); % elec_name,ROI,tissue(GG,GW,GW,WW,BS)
% % gm_elec_idx = logical(zeros([1 size(einfo,2)]));
% % for e_ix = 1:size(einfo,2)
% %     if ~isempty(strfind(einfo{3,e_ix},'G'))
% %         gm_elec_idx(e_ix) = true;
% %     end
% % end
% % einfo_gm = einfo(:,gm_elec_idx);
% roi_list = unique(einfo(2,:));

%% Process parameters
%!!! is this the best way to do this??? Maybe not...
sample_rate = (numel(hfa{1}.time)-1)/(hfa{1}.time(end)-hfa{1}.time(1));

% Confirm channels and time axis are the same
%   All the rounding and such is because some stupid rounding errors...
same_start = round(hfa_actv.time(1)*uint8(sample_rate))==round(hfa{1}.time(1)*uint8(sample_rate));
same_end   = round(hfa_actv.time(end)*uint8(sample_rate))==round(hfa{1}.time(end)*uint8(sample_rate));
same_numel = size(hfa_actv.time,2)==size(hfa{1}.time,2);
if ~same_start || ~same_end || ~same_numel
    error('time axes are not the same across hfa analyses!');
end
if ~isempty(setdiff(hfa{1}.label,hfa_actv.label));
    error('Different electrodes across hfa analyses!');
end

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
actv_code   = zeros([1, numel(stat.label)]);
cond_epochs = {};
cond_count  = 0;
both_count  = 0;
for ch_ix = 1:numel(stat.label)
    both_ix = 0;
    % Check if active, get epochs
    if ~isempty(strmatch(stat.label{ch_ix},actv_ch,'exact'))
        % Find significant epoch indices
        actv_ix = strmatch(stat.label{ch_ix},actv_ch,'exact');
        actv_epochs{ch_ix} = actv_ch_epochs{actv_ix};
        actv_ep_sign = NaN([1 size(actv_ch_epochs{actv_ix},1)]);
        sig_chunk_ix = NaN([1 2]);
        for ep_ix = 1:size(actv_ch_epochs{actv_ix},1)
            sig_chunk_ix = [find(hfa_actv.time==actv_ch_epochs{actv_ix}(ep_ix,1))...
                            find(hfa_actv.time==actv_ch_epochs{actv_ix}(ep_ix,2))];
            % Find sign of (de)activation
            if 0<=squeeze(mean(mean(hfa_actv.powspctrm(:,ch_ix,1,sig_chunk_ix(1):sig_chunk_ix(2)),1),4))
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
        both_ix = both_ix + 1;
    else
        actv_epochs{ch_ix} = [];
        actv_sign{ch_ix}   = 0;
    end
    
    % Check for condition differences, get epochs
    if sum(squeeze(stat.mask(ch_ix,1,:)))>0
        mask_chunks = fn_find_chunks(squeeze(stat.mask(ch_ix,1,:)));
        mask_chunks(squeeze(stat.mask(ch_ix,1,mask_chunks(:,1)))==0,:) = [];
        % Convert to time
        for ep_ix = 1:size(mask_chunks,1)
            mask_chunks(ep_ix,1) = stat.time(mask_chunks(ep_ix,1));
            mask_chunks(ep_ix,2) = stat.time(mask_chunks(ep_ix,2));
        end
        cond_epochs{ch_ix} = mask_chunks;
        cond_count = cond_count + 1;
        both_ix = both_ix + 1;
    else
        cond_epochs{ch_ix} = [];
    end
    
    % Check for overlap
    if both_ix==2
        both_count = both_count+1;
    end
end

%% Plot Results
% Creat and format the plot
fig_name = [SBJ '_HFA_summary_actv_cond_' event];
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
ax.YLim  = [0 numel(stat.label)+1];
ax.YTick = [1:numel(stat.label)];
ylabels  = stat.label;
for lab_ix = 2:2:numel(ylabels) % To space out the labels and make them readable
    ylabels{lab_ix} = ['     ' ylabels{lab_ix}];
end
ax.YTickLabel         = stat.label;
ax.YTickLabelRotation = 45;
% plot_info.legend_loc = plt_vars.legend_loc;

% Plot Summary info in title
perc_actv      = numel(find(actv_code~=0))/numel(stat.label);
perc_actv_pos  = numel(find(actv_code==1))/numel(stat.label);
perc_actv_neg  = numel(find(actv_code==-1))/numel(stat.label);
perc_actv_both = numel(find(actv_code==2))/numel(stat.label);
perc_cond      = cond_count/numel(stat.label);
perc_both      = both_count/numel(stat.label);

ax.Title.String = sprintf('Condition: %.2f; Active: %.2f (+=%.2f;-=%.2f;+-=%.2f); Both: %.2f',...
    perc_cond,perc_actv,perc_actv_pos,perc_actv_neg,perc_actv_both,perc_both);
ax.Title.FontSize = 16;

% Plot Event
if strcmp(an_id(1:5),'HGm_S')
    event_time = RT_mean;
elseif strcmp(an_id(1:5),'HGm_R')
    event_time = 0;
end
line([event_time event_time], [0 numel(stat.label)+1],...
    'LineWidth',plt_vars.evnt_width, 'Color',plt_vars.evnt_color,...
    'LineStyle',plt_vars.evnt_style);

% Plot Summary Epochs
for ch_ix = 1:numel(stat.label)
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
    
    % Plot epochs with significant condition differences
    for ep_ix = 1:size(cond_epochs{ch_ix},1)
        plot(cond_epochs{ch_ix}(ep_ix,1):(1/sample_rate):cond_epochs{ch_ix}(ep_ix,2),...
            ones([1 numel(cond_epochs{ch_ix}(ep_ix,1):(1/sample_rate):cond_epochs{ch_ix}(ep_ix,2))])*ch_ix,...
            'LineWidth',plt_vars.cond_width,'Color',plt_vars.cond_color);
    end
end

%% Save figure
if save_fig
    fig_dir = ['/home/knight/hoycw/PRJ_Stroop/results/HFA/' SBJ '/summary_' conditions '_actv/'...
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
