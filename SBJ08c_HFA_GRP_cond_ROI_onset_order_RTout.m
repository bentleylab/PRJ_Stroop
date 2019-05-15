function SBJ10c_HFA_GRP_ANOVA_ROI_onsets_RTout(SBJs,stat_id,pipeline_id,an_id,roi_id,plt_id,save_fig,fig_vis)
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
stat_vars_cmd = ['run /home/knight/hoycw/PRJ_Stroop/scripts/stat_vars/' stat_id '_vars.m'];
eval(stat_vars_cmd);

% Get condition info
[grp_lab, grp_colors, grp_style] = fn_group_label_styles(model_lab);
cond_lab = {'RT', grp_lab{:}};

% Get event timing
mean_RTs = zeros(size(SBJs));
if strcmp(an_id(1:5),'HGm_S')
    event_lab = 'stim';
elseif strcmp(an_id(1:5),'HGm_R')
    event_lab = 'resp';
end

% Load all ROI info
load('~/PRJ_Stroop/data/full_roi_lists.mat');
[roi_list, roi_colors, einfo_roi_col] = fn_roi_label_styles(roi_id);

% Set up onset counts
all_onsets  = cell([numel(SBJs) numel(roi_list) numel(grp_lab)+1]);

%% Load Results
for sbj_ix = 1:numel(SBJs)
    SBJ = SBJs{sbj_ix};
    % Load variables
    SBJ_vars_cmd = ['run /home/knight/hoycw/PRJ_Stroop/scripts/SBJ_vars/' SBJ '_vars.m'];
    eval(SBJ_vars_cmd);
    
    % Compute mean RT
    if strcmp(event_lab,'stim')
        load(strcat(SBJ_vars.dirs.events,SBJ,'_trial_info_final.mat'),'trial_info');
        % Compute mean RT
        mean_RTs(sbj_ix) = mean(trial_info.response_time);
    end
    % Load data
    load(strcat(SBJ_vars.dirs.proc,SBJ,'_ANOVA_ROI_',model_lab,'_',an_id,'.mat'));
    
    %% Process parameters
    %!!! is this the best way to do this??? Maybe not...
    sample_rate = (numel(hfa{1}.time)-1)/(hfa{1}.time(end)-hfa{1}.time(1));
    
    % FDR correct pvalues for ANOVA
    qvals = NaN(size(w2.pval));
    for ch_ix = 1:numel(stat.label)
        [~, ~, ~, qvals(:,ch_ix,:)] = fdr_bh(squeeze(w2.pval(:,ch_ix,:)));%,0.05,'pdep','yes');
    end
    
    % Get Sliding Window Parameters
    win_lim    = fn_sliding_window_lim(stat{1}.time,win_len*sample_rate,win_step*sample_rate);
    win_center = round(mean(win_lim,2));
    
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
    if ~isempty(setdiff(hfa{1}.label,einfo(:,1)))
        error('ERROR: Electrodes do not match between hfa and einfo!');
    end
    
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
        roi_ix = [];
        einfo_ix = strmatch(stat.label(ch_ix),einfo(:,1),'exact');
        % Check for condition differences, get epochs
        if strmatch(einfo(einfo_ix,einfo_roi_col),roi_list,'exact')
            roi_ix = strmatch(einfo(einfo_ix,einfo_roi_col),roi_list,'exact');
            
            % Get RT correlation onset
            if sum(squeeze(stat.mask(ch_ix,1,:)))>0
                mask_chunks = fn_find_chunks(squeeze(stat.mask(ch_ix,1,:)));
                mask_chunks(squeeze(stat.mask(ch_ix,1,mask_chunks(:,1)))==0,:) = [];
                % Convert the first onset of significance to time
                onset_time = stat.time(mask_chunks(1,1));
                % Exclude differences after the mean RT for this SBJ
                if strcmp(event_lab,'resp')
                    all_onsets{sbj_ix,roi_ix,1} = [all_onsets{sbj_ix,roi_ix,1} onset_time];
                elseif strcmp(event_lab,'stim') && (min(mask_chunks(:,1))<mean_RTs(sbj_ix))
                    all_onsets{sbj_ix,roi_ix,1} = [all_onsets{sbj_ix,roi_ix,1} onset_time];
                end
            end
            
            % Get ANOVA group onsets
            for grp_ix = 1:numel(grp_lab)
                if any(squeeze(qvals(grp_ix,ch_ix,:))<0.05)
                    sig_onsets = win_lim(squeeze(qvals(grp_ix,ch_ix,:))<0.05,1);
                    all_onsets{sbj_ix,roi_ix,grp_ix+1} = [all_onsets{sbj_ix,roi_ix,grp_ix+1} sig_onsets(1)];
                end
            end
        end
    end
    
    clear SBJ SBJ_vars hfa hfa_actv stat einfo w2
end
%% Compute medians per gROI
median_onsets = NaN([numel(SBJs) numel(roi_list) numel(grp_lab)+1]);
for sbj_ix = 1:numel(SBJs)
    for roi_ix = 1:numel(roi_list)
        median_onsets(sbj_ix,roi_ix,1) = median(all_onsets{sbj_ix,roi_ix,1});
        for grp_ix = 1:numel(grp_lab)
            median_onsets(sbj_ix,roi_ix,grp_ix+1) = median(all_onsets{sbj_ix,roi_ix,grp_ix+1});
        end
%         fprintf('%s , %s: %f (N=%i)\n',SBJs{sbj_ix},groi_list{groi_ix},...
%             median_onsets(sbj_ix,groi_ix),numel(cond_g_onsets{sbj_ix,groi_ix}));
%         disp(cond_g_onsets{sbj_ix,groi_ix});
%         fprintf('\n');
    end
end

%% Plot GROI Results
for cond_ix = 1:numel(cond_lab)
    % Create and format the plot
    fig_name = ['GRP_HFA_onsets_' stat_id '_' roi_id '_' event_lab '_meanRTout'];
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
    groi_y_offset      = linspace(-marker_offset,marker_offset,numel(roi_list));
    for sbj_ix = 1:numel(SBJs)
        s = [];
        for roi_ix = 1:numel(roi_list)
            s(roi_ix) = scatter(all_onsets{sbj_ix,roi_ix},...
                repmat(SBJ_ys(sbj_ix)+groi_y_offset(roi_ix),[1 numel(all_onsets{sbj_ix,roi_ix})]),...
                onset_marker_size,'o','filled');
            set(s(roi_ix),'MarkerFaceColor',roi_colors{roi_ix},'MarkerEdgeColor','k');
        end
        r = [];
        for cond_ix = 1:numel(grp_lab)
            r(cond_ix) = scatter(mean_RTs(cond_ix,sbj_ix)-plt_vars.plt_lim(1),SBJ_ys(sbj_ix),rt_marker_size,'+');
            set(r(cond_ix),'MarkerEdgeColor',cond_colors{cond_ix},'LineWidth',rt_marker_width);
        end
        for roi_ix = 1:numel(roi_list)
            line([median_onsets(sbj_ix,roi_ix) median_onsets(sbj_ix,roi_ix)],...
                [SBJ_ys(sbj_ix)-median_line_height/2 SBJ_ys(sbj_ix)+median_line_height/2],...
                'Color',roi_colors{roi_ix},'LineStyle','-','LineWidth',3);
        end
        if sbj_ix==4    % IR35 has all 4 main gROIs
            legend([s r],roi_list{:},['Mean RT: ' grp_lab{1}],['Mean RT: ' grp_lab{2}],'Location','northeast');
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
    hist_data = zeros([numel(roi_list) numel(plt_vars.x_data)-1]);
    for roi_ix = 1:numel(roi_list)
        hist_data(roi_ix,:)  = histcounts(median_onsets(:,roi_ix),plt_vars.x_data);
    end
    b = bar(plt_vars.x_data(1:end-1)+(diff(plt_vars.x_data(1:2))/2),hist_data',1,'stacked');
    for roi_ix = numel(roi_list):-1:1 %backwards so OFC doesn't overwrite LPFC mean onset
        set(b(roi_ix),'FaceColor',roi_colors{roi_ix},'EdgeColor','k');
        %     l(groi_ix) = line([nanmean(median_g_onsets(:,groi_ix)) nanmean(median_g_onsets(:,groi_ix))],...
        %         ax2.YLim,'Color',groi_colors{groi_ix},'LineStyle','--','LineWidth',3);
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
end

end
