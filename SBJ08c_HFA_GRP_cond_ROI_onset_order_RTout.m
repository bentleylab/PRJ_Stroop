function SBJ08c_HFA_GRP_cond_ROI_onset_order_RTout(SBJs,conditions,pipeline_id,an_id,plt_id,save_fig,fig_vis)
% Load HFA analysis results for active and condition-differentiating
%   epochs, plot a summary of those time period per electrode
% clear all; %close all;
fig_filetype = 'png';
label_spacer = 0;
groi_label_spacer = '      ';
if ischar(save_fig); save_fig = str2num(save_fig); end
% if isnumeric(actv_win); actv_win = num2str(actv_win); end

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

% Get event timing
[cond_lab, cond_colors, cond_style] = fn_condition_label_styles(conditions);
mean_RTs = NaN([numel(cond_lab)+1 numel(SBJs)]);
if strcmp(an_id(1:5),'HGm_S')
    event = 'stim';
elseif strcmp(an_id(1:5),'HGm_R')
    event = 'resp';
    event_time = 0;
end

% Load all ROI info
load('~/PRJ_Stroop/data/full_roi_lists.mat');
roi_list = all_rois;
roi_list(strcmp(roi_list,'')) = [];
roi_list(strcmp(roi_list,'FWM')) = [];
roi_list(strcmp(roi_list,'OUT')) = [];
roi_colors = cell(size(roi_list));
for roi_ix = 1:numel(roi_list)
    roi_colors{roi_ix} = fn_roi2color(roi_list{roi_ix});
end
% groi_list = all_grois;
groi_list = {'LPFC','MPFC','INS','OFC'};
groi_colors = cell(size(groi_list));
for groi_ix = 1:numel(groi_list)
    groi_colors{groi_ix} = fn_roi2color(groi_list{groi_ix});
end

% Set up onset counts
cond_onsets   = cell([numel(SBJs) numel(roi_list)]);
cond_g_onsets = cell([numel(SBJs) numel(groi_list)]);
% cond_g_rois   = cell([numel(SBJs) 1]);

%% Load Results
for sbj_ix = 1:numel(SBJs)
    SBJ = SBJs{sbj_ix};
    % Load variables
    SBJ_vars_cmd = ['run /home/knight/hoycw/PRJ_Stroop/scripts/SBJ_vars/' SBJ '_vars.m'];
    eval(SBJ_vars_cmd);
    
    % Compute mean RT
    if strcmp(event,'stim')
        load(strcat(SBJ_vars.dirs.events,SBJ,'_trial_info_final.mat'),'trial_info');
        % Compute mean RT per condition
        for cond_ix = 1:numel(cond_lab)
            mean_RTs(cond_ix,sbj_ix) = mean(trial_info.response_time(fn_condition_index([cond_lab{cond_ix}],...
                                                                    trial_info.condition_n)==1))-plt_vars.plt_lim(1);
        end
        mean_RTs(numel(cond_lab)+1,sbj_ix) = mean(trial_info.response_time)-plt_vars.plt_lim(1);
    end
    % Load data
    load(strcat(SBJ_vars.dirs.proc,SBJ,'_',conditions,'_ROI_',an_id,'.mat'));
%     actv_filename = strcat(SBJ_vars.dirs.proc,SBJ,'_actv_ROI_',an_id,'_mn',actv_win,'.mat');
%     load(actv_filename,'actv_ch','actv_ch_epochs');
%     tmp = load(actv_filename,'hfa'); hfa_actv = tmp.hfa;
    
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
%     same_start = round(hfa_actv.time(1)*uint8(sample_rate))==round(hfa{1}.time(1)*uint8(sample_rate));
%     same_end   = round(hfa_actv.time(end)*uint8(sample_rate))==round(hfa{1}.time(end)*uint8(sample_rate));
%     same_numel = size(hfa_actv.time,2)==size(hfa{1}.time,2);
%     if ~same_start || ~same_end || ~same_numel
%         error('time axes are not the same across hfa analyses!');
%     end
%     if ~isempty(setdiff(hfa{1}.label,hfa_actv.label));
%         error('Different electrodes across hfa analyses!');
%     end
    
    %% Aggregate results per ROI
    for ch_ix = 1:numel(stat.label)
        groi_ix = [];
        roi_ix = [];
        einfo_ix = strmatch(stat.label(ch_ix),einfo(:,1),'exact');
        % Check for condition differences, get epochs
        if sum(squeeze(stat.mask(ch_ix,1,:)))>0
            mask_chunks = fn_find_chunks(squeeze(stat.mask(ch_ix,1,:)));
            mask_chunks(squeeze(stat.mask(ch_ix,1,mask_chunks(:,1)))==0,:) = [];
            % Convert to time
            for ep_ix = 1:size(mask_chunks,1)
                mask_chunks(ep_ix,1) = stat.time(mask_chunks(ep_ix,1));
                mask_chunks(ep_ix,2) = stat.time(mask_chunks(ep_ix,2));
            end
            % Exclude differences after the mean RT for this SBJ
            if min(mask_chunks(:,1))<mean_RTs(numel(cond_lab)+1,sbj_ix)
                roi_ix = strmatch(einfo(einfo_ix,2),roi_list,'exact');
                if ~isempty(roi_ix)
                    cond_onsets{sbj_ix,roi_ix} = [cond_onsets{sbj_ix,roi_ix} min(mask_chunks(:,1))];
                end
                groi_ix = strmatch(einfo(einfo_ix,3),groi_list,'exact');
                if ~isempty(groi_ix)
                    cond_g_onsets{sbj_ix,groi_ix} = [cond_g_onsets{sbj_ix,groi_ix} min(mask_chunks(:,1))];
                end
            end
        end
    end
    
%     clear SBJ SBJ_vars hfa hfa_actv stat einfo actv_ch actv_ch_epochs tmp
end
%% Compute medians per gROI
median_onsets   = NaN([numel(SBJs) numel(roi_list)]);
median_g_onsets = NaN([numel(SBJs) numel(groi_list)]);
for sbj_ix = 1:numel(SBJs)
    for roi_ix = 1:numel(roi_list)
        median_onsets(sbj_ix,roi_ix) = median(cond_onsets{sbj_ix,roi_ix});
    end
    for groi_ix = 1:numel(groi_list)
        median_g_onsets(sbj_ix,groi_ix) = median(cond_g_onsets{sbj_ix,groi_ix});
%         fprintf('%s , %s: %f (N=%i)\n',SBJs{sbj_ix},groi_list{groi_ix},...
%             median_onsets(sbj_ix,groi_ix),numel(cond_g_onsets{sbj_ix,groi_ix}));
%         disp(cond_g_onsets{sbj_ix,groi_ix});
%         fprintf('\n');
    end
end

%% Plot GROI Results
% Create and format the plot
fig_name = ['GRP_HFA_cond_onsets_gROI_' event '_meanRTout'];
figure('Name',fig_name,'units','normalized',...
        'outerposition',[0 0 0.9 1],'Visible',fig_vis);

% Plot Condition Difference Onsets per SBJ
ax = subplot(2,1,1);
hold on;
median_line_height = 0.3;
SBJ_y_spacer       = 0.35;
SBJ_ys             = SBJ_y_spacer*[1:numel(SBJs)]-median_line_height/3;
rt_marker_width    = 3;
rt_marker_size     = 250;
onset_marker_size  = 150;
marker_offset      = 0.08;
groi_y_offset      = linspace(-marker_offset,marker_offset,numel(groi_list));
for sbj_ix = 1:numel(SBJs)
    s = [];
    for groi_ix = 1:numel(groi_list)
        s(groi_ix) = scatter(cond_g_onsets{sbj_ix,groi_ix},...
            repmat(SBJ_ys(sbj_ix)+groi_y_offset(groi_ix),[1 numel(cond_g_onsets{sbj_ix,groi_ix})]),...
            onset_marker_size,'o','filled');
        set(s(groi_ix),'MarkerFaceColor',groi_colors{groi_ix},'MarkerEdgeColor','k');
    end
    r = [];
    for cond_ix = 1:numel(cond_lab)
        r(cond_ix) = scatter(mean_RTs(cond_ix,sbj_ix),SBJ_ys(sbj_ix),rt_marker_size,'+');
        set(r(cond_ix),'MarkerEdgeColor',cond_colors{cond_ix},'LineWidth',rt_marker_width);
    end
    for groi_ix = 1:numel(groi_list)
        line([median_g_onsets(sbj_ix,groi_ix) median_g_onsets(sbj_ix,groi_ix)],...
            [SBJ_ys(sbj_ix)-median_line_height/2 SBJ_ys(sbj_ix)+median_line_height/2],...
            'Color',groi_colors{groi_ix},'LineStyle','-','LineWidth',3);
    end
    if sbj_ix==4    % IR35 has all 4 main gROIs
        legend([s r],groi_list{:},['Mean RT: ' cond_lab{1}],['Mean RT: ' cond_lab{2}],'Location','northeast');
    end
    line([plt_vars.plt_lim(1) plt_vars.plt_lim(2)],[SBJ_ys(sbj_ix) SBJ_ys(sbj_ix)],...
        'Color',[0.2 0.2 0.2],'LineStyle',':');
end

% ax = gca;
% Plot labels
% ax.XLabel.String   = 'Time (s)';
% ax.XLabel.FontSize = 14;
ax.XLim    = plt_vars.plt_lim;
ax.XTick   = plt_vars.plt_lim(1):plt_vars.x_step_sz:plt_vars.plt_lim(2);
ax.XColor  = 'k';

ax.YLabel.String   = 'Subject';
ax.YLabel.FontSize = 14;
ax.YLim            = [0 max(SBJ_ys)+median_line_height/2];
ax.YTick           = SBJ_ys;        % Ticks anywhere non-zero in plot_idx
ax.YTickLabel      = SBJs;
% ax.YTickLabelRotation = 45;
ax.YColor  = 'k';

ax.Title.String = 'Condition Difference Onsets by gROI in Individuals';
ax.Title.FontSize = 16;

% Plot histogram of the median onsets by ROI across subjects
ax2 = subplot(2,1,2);
hold on;
% ax2 = gca;
b = [];
l = [];
% hist_alpha = 0.6;
hist_data = zeros([numel(groi_list) numel(plt_vars.x_data)-1]);
for groi_ix = 1:numel(groi_list)
    hist_data(groi_ix,:)  = histcounts(median_g_onsets(:,groi_ix),plt_vars.x_data);
end
b = bar(plt_vars.x_data(1:end-1)+(diff(plt_vars.x_data(1:2))/2),hist_data',1,'stacked');
for groi_ix = numel(groi_list):-1:1 %backwards so OFC doesn't overwrite LPFC mean onset
    set(b(groi_ix),'FaceColor',groi_colors{groi_ix},'EdgeColor','k');
%     l(groi_ix) = line([nanmean(median_g_onsets(:,groi_ix)) nanmean(median_g_onsets(:,groi_ix))],...
%         ax2.YLim,'Color',groi_colors{groi_ix},'LineStyle','--','LineWidth',3);
end
legend(b,groi_list{:},'Location','northeast');
% Plot labels
ax2.XLabel.String   = 'Time (s)';
ax2.XLabel.FontSize = 14;
ax2.XLim    = plt_vars.plt_lim;
ax2.XTick   = plt_vars.plt_lim(1):plt_vars.x_step_sz:plt_vars.plt_lim(2);
ax2.XColor  = 'k';

ax2.YLabel.String   = '# Median Onsets';
ax2.YLabel.FontSize = 14;
% ax2.YLim            = [0 numel(SBJs)+1];
ax2.YTick           = 1:1:ax2.YLim(2);        % Ticks anywhere non-zero in plot_idx
% ax2.YTickLabel      = SBJs;
% ax2.YTickLabelRotation = 45;
ax2.YColor  = 'k';

ax2.Title.String = 'Group-Level Median Condition Difference Onsets by gROI';
ax2.Title.FontSize = 16;

% Reposition axes
set(ax, 'Position',[0.05 0.33 0.9 0.63]);    %[left bottom width height]
set(ax2,'Position',[0.05 0.05 0.9 0.21]);    %[left bottom width height]

%% Save figure
if save_fig
    fig_dir = ['/home/knight/hoycw/PRJ_Stroop/results/HFA/GRP_' conditions '_onsets_ROI/'...
        an_id '/'];
    if ~exist(fig_dir,'dir')
        mkdir(fig_dir);
    end
    
    fig_filename = [fig_dir fig_name '.' fig_filetype];
    fprintf('Saving %s\n',fig_filename);
    saveas(gcf,fig_filename);
    %eval(['export_fig ' fig_filename]);
end

%% Plot ROI Results
% Create and format the plot
fig_name = ['GRP_HFA_cond_onsets_ROI_' event '_meanRTout'];
figure('Name',fig_name,'units','normalized',...
        'outerposition',[0 0 0.9 1],'Visible',fig_vis);

% Plot Condition Difference Onsets per SBJ
ax = subplot(2,1,1);
hold on;
roi_y_offset = linspace(-marker_offset,marker_offset,numel(roi_list));
for sbj_ix = 1:numel(SBJs)
    s = [];
    for roi_ix = 1:numel(roi_list)
        s(roi_ix) = scatter(cond_onsets{sbj_ix,roi_ix},...
            repmat(SBJ_ys(sbj_ix)+roi_y_offset(roi_ix),[1 numel(cond_onsets{sbj_ix,roi_ix})]),...
            onset_marker_size,'o','filled');
        set(s(roi_ix),'MarkerFaceColor',roi_colors{roi_ix},'MarkerEdgeColor','k');
    end
    r = [];
    for cond_ix = 1:numel(cond_lab)
        r(cond_ix) = scatter(mean_RTs(cond_ix,sbj_ix),SBJ_ys(sbj_ix),rt_marker_size,'+');
        set(r(cond_ix),'MarkerEdgeColor',cond_colors{cond_ix},'LineWidth',rt_marker_width);
    end
    for roi_ix = 1:numel(roi_list)
        line([median_onsets(sbj_ix,roi_ix) median_onsets(sbj_ix,roi_ix)],...
            [SBJ_ys(sbj_ix)-median_line_height/2 SBJ_ys(sbj_ix)+median_line_height/2],...
            'Color',roi_colors{roi_ix},'LineStyle','-','LineWidth',3);
    end
    line([plt_vars.plt_lim(1) plt_vars.plt_lim(2)],[SBJ_ys(sbj_ix) SBJ_ys(sbj_ix)],...
        'Color',[0.2 0.2 0.2],'LineStyle',':');
end
scat_obj = findobj(gca,'Type','Scatter');
roi_scat_idx = logical([1 numel(scat_obj)]);
roi_scat_lgd = zeros([1 numel(scat_obj)]);
roi_rt_lab = {roi_list{:} ['Mean RT: ' cond_lab{1}],['Mean RT: ' cond_lab{2}]};
for scat_ix = 1:numel(scat_obj)
    if ~isempty(scat_obj(scat_ix).DisplayName)
        roi_scat_idx(scat_ix) = 1;
        roi_scat_lgd(scat_ix) = strmatch(scat_obj(scat_ix).DisplayName,roi_rt_lab);
    end
end
legend(scat_obj(roi_scat_idx),roi_rt_lab{roi_scat_lgd(roi_scat_idx)},'Location','northeast');

% ax = gca;
% Plot labels
% ax.XLabel.String   = 'Time (s)';
% ax.XLabel.FontSize = 14;
ax.XLim    = plt_vars.plt_lim;
ax.XTick   = plt_vars.plt_lim(1):plt_vars.x_step_sz:plt_vars.plt_lim(2);
ax.XColor  = 'k';

ax.YLabel.String   = 'Subject';
ax.YLabel.FontSize = 14;
ax.YLim            = [0 max(SBJ_ys)+median_line_height/2];
ax.YTick           = SBJ_ys;        % Ticks anywhere non-zero in plot_idx
ax.YTickLabel      = SBJs;
% ax.YTickLabelRotation = 45;
ax.YColor  = 'k';

ax.Title.String = 'Condition Difference Onsets by ROI in Individuals';
ax.Title.FontSize = 16;

% Plot histogram of the median onsets by ROI across subjects
ax2 = subplot(2,1,2);
hold on;
% ax2 = gca;
b = [];
l = [];
% hist_alpha = 0.6;
hist_data = zeros([numel(roi_list) numel(plt_vars.x_data)-1]);
for roi_ix = 1:numel(roi_list)
    hist_data(roi_ix,:)  = histcounts(median_onsets(:,roi_ix),plt_vars.x_data);
end
b = bar(plt_vars.x_data(1:end-1)+(diff(plt_vars.x_data(1:2))/2),hist_data',1,'stacked');
for roi_ix = numel(roi_list):-1:1 %backwards so OFC doesn't overwrite LPFC mean onset
    set(b(roi_ix),'FaceColor',roi_colors{roi_ix},'EdgeColor','k');
%     l(roi_ix) = line([nanmean(median_onsets(:,roi_ix)) nanmean(median_onsets(:,roi_ix))],...
%         ax2.YLim,'Color',roi_colors{roi_ix},'LineStyle','--','LineWidth',3);
end
legend(b,roi_list{:},'Location','northeast');
% Plot labels
ax2.XLabel.String   = 'Time (s)';
ax2.XLabel.FontSize = 14;
ax2.XLim    = plt_vars.plt_lim;
ax2.XTick   = plt_vars.plt_lim(1):plt_vars.x_step_sz:plt_vars.plt_lim(2);
ax2.XColor  = 'k';

ax2.YLabel.String   = '# Median Onsets';
ax2.YLabel.FontSize = 14;
% ax2.YLim            = [0 numel(SBJs)+1];
ax2.YTick           = 1:1:ax2.YLim(2);        % Ticks anywhere non-zero in plot_idx
% ax2.YTickLabel      = SBJs;
% ax2.YTickLabelRotation = 45;
ax2.YColor  = 'k';

ax2.Title.String = 'Group-Level Median Condition Difference Onsets by ROI';
ax2.Title.FontSize = 16;

% Reposition axes
set(ax, 'Position',[0.05 0.33 0.9 0.63]);    %[left bottom width height]
set(ax2,'Position',[0.05 0.05 0.9 0.21]);    %[left bottom width height]

%% Save figure
if save_fig
    fig_filename = [fig_dir fig_name '.' fig_filetype];
    fprintf('Saving %s\n',fig_filename);
    saveas(gcf,fig_filename);
    %eval(['export_fig ' fig_filename]);
end

end
