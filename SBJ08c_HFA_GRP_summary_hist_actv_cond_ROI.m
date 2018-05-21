function SBJ08c_HFA_GRP_summary_hist_actv_cond_ROI(SBJs,conditions,pipeline_id,an_id,actv_win,plt_id,save_fig,fig_vis)
% Load HFA analysis results for active and condition-differentiating
%   epochs, plot a summary of those time period per electrode
% clear all; %close all;
fig_filetype = 'png';
label_spacer = 0;
groi_label_spacer = '      ';
if ischar(save_fig); save_fig = str2num(save_fig); end
if isnumeric(actv_win); actv_win = num2str(actv_win); end

%% Data Preparation
% Set up paths
addpath('/home/knight/hoycw/PRJ_Stroop/scripts/');
addpath('/home/knight/hoycw/PRJ_Stroop/scripts/utils/');
addpath('/home/knight/hoycw/Apps/fieldtrip/');
ft_defaults

%% Prep variables
an_vars_cmd = ['run /home/knight/hoycw/PRJ_Stroop/scripts/an_vars/' an_id '_vars.m'];
eval(an_vars_cmd);
plt_vars_cmd = ['run /home/knight/hoycw/PRJ_Stroop/scripts/plt_vars/' plt_id '_vars.m'];
eval(plt_vars_cmd);

% Load all ROI info
load('~/PRJ_Stroop/data/full_roi_lists.mat');
roi_list = all_rois;
roi_colors = cell(size(roi_list));
for roi_ix = 1:numel(roi_list)
    roi_colors{roi_ix} = fn_roi2color(roi_list{roi_ix});
end
groi_list = all_grois;
groi_colors = cell(size(groi_list));
for groi_ix = 1:numel(groi_list)
    groi_colors{groi_ix} = fn_roi2color(groi_list{groi_ix});
end

% Get event timing
if strcmp(an_id(1:5),'HGm_S')
    event = 'stim';
elseif strcmp(an_id(1:5),'HGm_R')
    event = 'resp';
    event_time = 0;
end

% Set up onset counts
actv_onsets = cell(size(roi_list));
dact_onsets = cell(size(roi_list));
actv_code   = cell(size(roi_list));
cond_onsets = cell(size(roi_list));
cond_count  = zeros([1 numel(roi_list)]);
both_count  = zeros([1 numel(roi_list)]);
actv_g_onsets = cell(size(groi_list));
dact_g_onsets = cell(size(groi_list));
actv_g_code   = cell(size(groi_list));
cond_g_onsets = cell(size(groi_list));
cond_g_count  = zeros([1 numel(groi_list)]);
both_g_count  = zeros([1 numel(groi_list)]);

mean_RTs = NaN(size(SBJs));

%% Load Results
for sbj_ix = 1:numel(SBJs)
    SBJ = SBJs{sbj_ix};
    % Load variables
    SBJ_vars_cmd = ['run /home/knight/hoycw/PRJ_Stroop/scripts/SBJ_vars/' SBJ '_vars.m'];
    eval(SBJ_vars_cmd);
    
    % Compute mean RT
    if strcmp(event,'stim')
        load(strcat(SBJ_vars.dirs.events,SBJ,'_trial_info_final.mat'),'trial_info');
        mean_RTs(sbj_ix) = mean(trial_info.response_time); % converts sec to ms
    end
    % Load data
    load(strcat(SBJ_vars.dirs.proc,SBJ,'_',conditions,'_ROI_',an_id,'.mat'));
    actv_filename = strcat(SBJ_vars.dirs.proc,SBJ,'_actv_ROI_',an_id,'_mn',actv_win,'.mat');
    load(actv_filename,'actv_ch','actv_ch_epochs');
    tmp = load(actv_filename,'hfa'); hfa_actv = tmp.hfa;
    
    %% Load ROI and GM/WM info
    einfo_filename = [SBJ_vars.dirs.preproc SBJ '_einfo_' pipeline_id '.mat'];
    load(einfo_filename);
    % Electrode Info Table:
    %   label- name of electrode
    %   ROI- specific region
    %   gROI- general region (LPFC, MPFC, OFC, FWM=frontal white matter)
    %   ROI2- specific region of second electrode
    %   tissue- primary tissue type
    %   GM weight- percentage of electrode pair in GM
    %   Out- 0/1 flag for whether this is partially out of the brain
    
    % Sort by gROI, then ROI
    einfo = sortrows(einfo,[3,2]);
    if ~isempty(setdiff(einfo(:,1),hfa{1}.label))
        error('ERROR: Electrodes do not match between hfa and einfo!');
    end
    
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
    
    %% Aggregate results per ROI
    for ch_ix = 1:numel(stat.label)
        both_ix = 0;
        einfo_ix = strmatch(stat.label(ch_ix),einfo(:,1),'exact');
        roi_ix = strmatch(einfo(einfo_ix,2),roi_list,'exact');
        groi_ix = strmatch(einfo(einfo_ix,3),groi_list,'exact');
        % Check if active, get epochs
        if ~isempty(strmatch(stat.label{ch_ix},actv_ch,'exact'))
            % Find significant epoch indices
            actv_ix = strmatch(stat.label{ch_ix},actv_ch,'exact');
            actv_epochs = actv_ch_epochs{actv_ix};
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
            % Only take onsets for activations (not deactivations)
            if ~isempty(find(actv_ep_sign==1,1))
                actv_onsets{roi_ix} = [actv_onsets{roi_ix} actv_epochs(find(actv_ep_sign==1,1),1)];
                actv_g_onsets{groi_ix} = [actv_g_onsets{groi_ix} actv_epochs(find(actv_ep_sign==1,1),1)];
            end
            if ~isempty(find(actv_ep_sign==-1,1))
                dact_onsets{roi_ix} = [dact_onsets{roi_ix} actv_epochs(find(actv_ep_sign==-1,1),1)];
                dact_g_onsets{groi_ix} = [dact_g_onsets{groi_ix} actv_epochs(find(actv_ep_sign==-1,1),1)];
            end
            
            % Code significance and sign: 0=none, -1=deactive, 1=active, 2=both
            if any(actv_ep_sign==1) && any(actv_ep_sign==-1)
                actv_code{roi_ix}    = [actv_code{roi_ix} 2];
                actv_g_code{groi_ix} = [actv_g_code{groi_ix} 2];
            elseif any(actv_ep_sign==1)
                actv_code{roi_ix}    = [actv_code{roi_ix} 1];
                actv_g_code{groi_ix} = [actv_g_code{groi_ix} 1];
            elseif any(actv_ep_sign==-1)
                actv_code{roi_ix}    = [actv_code{roi_ix} -1];
                actv_g_code{groi_ix} = [actv_g_code{groi_ix} -1];
            end
            both_ix = both_ix + 1;
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
            cond_onsets{roi_ix} = [cond_onsets{roi_ix} min(mask_chunks(:,1))];
            cond_g_onsets{groi_ix} = [cond_g_onsets{groi_ix} min(mask_chunks(:,1))];
            cond_count(roi_ix) = cond_count(roi_ix) + 1;
            cond_g_count(groi_ix) = cond_g_count(groi_ix) + 1;
            both_ix = both_ix + 1;
        end
        
        % Check for overlap
        if both_ix==2
            both_count(roi_ix) = both_count(roi_ix)+1;
            both_g_count(groi_ix) = both_g_count(groi_ix)+1;
        end
    end
    
    clear SBJ SBJ_vars hfa hfa_actv stat einfo actv_ch actv_ch_epochs tmp trial_info
end
%% Compile into histograms
actv_onsets_hist = zeros([numel(roi_list) numel(plt_vars.x_data)-1]);
dact_onsets_hist = zeros([numel(roi_list) numel(plt_vars.x_data)-1]);
cond_onsets_hist = zeros([numel(roi_list) numel(plt_vars.x_data)-1]);
for roi_ix = 1:numel(roi_list)
    actv_onsets_hist(roi_ix,:) = histcounts(actv_onsets{roi_ix},plt_vars.x_data);
    dact_onsets_hist(roi_ix,:) = histcounts(dact_onsets{roi_ix},plt_vars.x_data);
    cond_onsets_hist(roi_ix,:) = histcounts(cond_onsets{roi_ix},plt_vars.x_data);
end
actv_g_onsets_hist = zeros([numel(groi_list) numel(plt_vars.x_data)-1]);
dact_g_onsets_hist = zeros([numel(groi_list) numel(plt_vars.x_data)-1]);
cond_g_onsets_hist = zeros([numel(groi_list) numel(plt_vars.x_data)-1]);
for groi_ix = 1:numel(groi_list)
    actv_g_onsets_hist(groi_ix,:) = histcounts(actv_g_onsets{groi_ix},plt_vars.x_data);
    dact_g_onsets_hist(groi_ix,:) = histcounts(dact_g_onsets{groi_ix},plt_vars.x_data);
    cond_g_onsets_hist(groi_ix,:) = histcounts(cond_g_onsets{groi_ix},plt_vars.x_data);
end

% Add in the baseline offset to plot correctly
if strcmp(event,'stim')
    event_time = mean(mean_RTs)-plt_vars.plt_lim(1);
end

%% Plot Results
% Create and format the plot
fig_name = ['GRP_HFA_ROI_hist_onset_actv_cond_ROI_' event];
figure('Name',fig_name,'units','normalized',...
        'outerposition',[0 0 1 1],'Visible',fig_vis);

% Plot Activations by ROI
bars = {};
ax   = {};
titles = {'Activation','Deactivation','Condition Differences'};
ymaxs = [max(sum(actv_onsets_hist,1)) max(sum(dact_onsets_hist,1)) max(sum(cond_onsets_hist,1))];
for plot_ix = [1 3 2]   %plot middle last because it has a long legend
    subplot(3,1,plot_ix);
    if plot_ix==1
        bars{plot_ix} = bar(plt_vars.x_data(1:end-1)+(diff(plt_vars.x_data(1:2))/2),actv_onsets_hist',1,'stack');
    elseif plot_ix==2
        bars{plot_ix} = bar(plt_vars.x_data(1:end-1)+(diff(plt_vars.x_data(1:2))/2),dact_onsets_hist',1,'stack');
    else
        bars{plot_ix} = bar(plt_vars.x_data(1:end-1)+(diff(plt_vars.x_data(1:2))/2),cond_onsets_hist',1,'stack');
    end
    for roi_ix = 1:numel(roi_list)
        set(bars{plot_ix}(roi_ix), 'FaceColor', roi_colors{roi_ix});
        set(bars{plot_ix}(roi_ix), 'EdgeColor', 'k');
    end
    if plot_ix==2
        legend(roi_list,'Location',plt_vars.legend_loc);
    end
    ax{plot_ix} = gca;
    
    % Plot Event
    line([event_time event_time], ax{1}.YLim,...
        'LineWidth',plt_vars.evnt_width, 'Color',plt_vars.evnt_color,...
        'LineStyle',plt_vars.evnt_style);
    
    % Plot labels
    ax{plot_ix}.XLabel.String   = 'Time (s)';
    ax{plot_ix}.XLabel.FontSize = 14;
    ax{plot_ix}.XLim    = plt_vars.plt_lim;
    ax{plot_ix}.XTick   = plt_vars.plt_lim(1):plt_vars.x_step_sz:plt_vars.plt_lim(2);
    ax{plot_ix}.XColor  = 'k';
    
    ax{plot_ix}.YLabel.String   = '# Onsets';
    ax{plot_ix}.YLabel.FontSize = 14;
    ax{plot_ix}.YLim            = [0 ymaxs(plot_ix)];
%     ax{plot_ix}.YTick           = ax{plot_ix}.YLim(1):1:ax{plot_ix}.YLim(2);        % Ticks anywhere non-zero in plot_idx
    % ax.YTickLabel      = roi_list;
    % ax.YTickLabelRotation = 45;
    ax{plot_ix}.YColor  = 'k';
    
    ax{plot_ix}.Title.String = titles{plot_ix};
%     ax.Title.String = sprintf('Condition: %.2f; Active: %.2f (+=%.2f;-=%.2f;+-=%.2f); Both: %.2f',...
%         perc_cond,perc_actv,perc_actv_pos,perc_actv_neg,perc_actv_both,perc_both);
    ax{plot_ix}.Title.FontSize = 16;
end

%% Save figure
if save_fig
    fig_dir = ['/home/knight/hoycw/PRJ_Stroop/results/HFA/GRP_summary_' conditions '_actv_ROI/'...
        an_id '_mn' actv_win '/'];
    if ~exist(fig_dir,'dir')
        mkdir(fig_dir);
    end
    
    fig_filename = [fig_dir fig_name '.' fig_filetype];
    fprintf('Saving %s\n',fig_filename);
    saveas(gcf,fig_filename);
    %eval(['export_fig ' fig_filename]);
end

%% Plot for General ROIs
% Create and format the plot
fig_name = ['GRP_HFA_ROI_hist_onset_actv_cond_gROI_' event];
figure('Name',fig_name,'units','normalized',...
        'outerposition',[0 0 1 1],'Visible',fig_vis);

% Plot Activations by ROI
bars = {};
ax   = {};
titles = {'Activation','Deactivation','Condition Differences'};
ymaxs = [max(sum(actv_g_onsets_hist,1)) max(sum(dact_g_onsets_hist,1)) max(sum(cond_g_onsets_hist,1))];
for plot_ix = 1:3
    subplot(3,1,plot_ix);
    if plot_ix==1
        bars{plot_ix} = bar(plt_vars.x_data(1:end-1)+(diff(plt_vars.x_data(1:2))/2),actv_g_onsets_hist',1,'stack');
    elseif plot_ix==2
        bars{plot_ix} = bar(plt_vars.x_data(1:end-1)+(diff(plt_vars.x_data(1:2))/2),dact_g_onsets_hist',1,'stack');
    else
        bars{plot_ix} = bar(plt_vars.x_data(1:end-1)+(diff(plt_vars.x_data(1:2))/2),cond_g_onsets_hist',1,'stack');
    end
    for groi_ix = 1:numel(groi_list)
        set(bars{plot_ix}(groi_ix), 'FaceColor', groi_colors{groi_ix});
        set(bars{plot_ix}(groi_ix), 'EdgeColor', 'k');
    end
    legend(groi_list,'Location',plt_vars.legend_loc);
    ax{plot_ix} = gca;
    
    % Plot Event
    line([event_time event_time], ax{1}.YLim,...
        'LineWidth',plt_vars.evnt_width, 'Color',plt_vars.evnt_color,...
        'LineStyle',plt_vars.evnt_style);
    
    % Plot labels
    ax{plot_ix}.XLabel.String   = 'Time (s)';
    ax{plot_ix}.XLabel.FontSize = 14;
    ax{plot_ix}.XLim    = plt_vars.plt_lim;
    ax{plot_ix}.XTick   = plt_vars.plt_lim(1):plt_vars.x_step_sz:plt_vars.plt_lim(2);
    ax{plot_ix}.XColor  = 'k';
    
    ax{plot_ix}.YLabel.String   = '# Onsets';
    ax{plot_ix}.YLabel.FontSize = 14;
    ax{plot_ix}.YLim            = [0 ymaxs(plot_ix)];
%     ax{plot_ix}.YTick           = ax{plot_ix}.YLim(1):1:ax{plot_ix}.YLim(2);        % Ticks anywhere non-zero in plot_idx
    % ax.YTickLabel      = roi_list;
    % ax.YTickLabelRotation = 45;
    ax{plot_ix}.YColor  = 'k';
    
    ax{plot_ix}.Title.String = titles{plot_ix};
%     ax.Title.String = sprintf('Condition: %.2f; Active: %.2f (+=%.2f;-=%.2f;+-=%.2f); Both: %.2f',...
%         perc_cond,perc_actv,perc_actv_pos,perc_actv_neg,perc_actv_both,perc_both);
    ax{plot_ix}.Title.FontSize = 16;
end

%% Save figure
if save_fig
    fig_dir = ['/home/knight/hoycw/PRJ_Stroop/results/HFA/GRP_summary_' conditions '_actv_ROI/'...
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
