function SBJ08c_HFA_GRP_summary_bar_perc_actv_cond_ROI(SBJs,conditions,pipeline_id,an_id,actv_win,plt_id,save_fig,fig_vis)
% Load HFA analysis results for active and condition-differentiating epochs
%   Plot bar chart with % active, % deactivated, % condition sensitive
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
% roi_list = all_rois;
% roi_colors = cell(size(roi_list));
% for roi_ix = 1:numel(roi_list)
%     roi_colors{roi_ix} = fn_roi2color(roi_list{roi_ix});
% end
groi_list = {'LPFC','MPFC','INS','OFC'};
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
% actv_onsets = cell(size(roi_list));
% dact_onsets = cell(size(roi_list));
% actv_code   = cell(size(roi_list));
% cond_onsets = cell(size(roi_list));
% cond_count  = zeros([1 numel(roi_list)]);
% both_count  = zeros([1 numel(roi_list)]);
% actv_g_onsets = cell(size(groi_list));
% dact_g_onsets = cell(size(groi_list));
actv_g_count  = zeros([numel(SBJs) numel(groi_list)]);
dact_g_count  = zeros([numel(SBJs) numel(groi_list)]);
% cond_g_onsets = cell(size(groi_list));
cond_g_count  = zeros([numel(SBJs) numel(groi_list)]);
both_g_count  = zeros([numel(SBJs) numel(groi_list)]);
elec_g_count  = zeros([numel(SBJs) numel(groi_list)]);

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
        %         roi_ix = strmatch(einfo(einfo_ix,2),roi_list,'exact');
        if ~isempty(strmatch(einfo(einfo_ix,3),groi_list,'exact'))
            groi_ix = strmatch(einfo(einfo_ix,3),groi_list,'exact');
            elec_g_count(sbj_ix,groi_ix) = elec_g_count(sbj_ix,groi_ix)+1;
            % Check if active, get epochs
            if ~isempty(strmatch(stat.label{ch_ix},actv_ch,'exact'))
                % Find significant epoch indices
                actv_ix = strmatch(stat.label{ch_ix},actv_ch,'exact');
                actv_epochs = actv_ch_epochs{actv_ix};
                % Toss late epcohs if S-locked
                if strcmp('stim',event)
                    actv_epochs(actv_epochs(:,1)<mean_RTs(sbj_ix),:) = [];
                end
                actv_ep_sign = NaN([1 size(actv_epochs,1)]);
                sig_chunk_ix = NaN([1 2]);
                for ep_ix = 1:size(actv_epochs,1)
                    sig_chunk_ix = [find(hfa_actv.time==actv_epochs(ep_ix,1))...
                        find(hfa_actv.time==actv_epochs(ep_ix,2))];
                    % Find sign of (de)activation
                    if 0<=squeeze(mean(mean(hfa_actv.powspctrm(:,ch_ix,1,sig_chunk_ix(1):sig_chunk_ix(2)),1),4))
                        actv_ep_sign(ep_ix) = 1;
                    else
                        actv_ep_sign(ep_ix) = -1;
                    end
                end
                
                % Code significance and sign: 0=none, -1=deactive, 1=active, 2=both
                if any(actv_ep_sign==1)
                    %                 actv_code{roi_ix}    = [actv_code{roi_ix} 1];
                    actv_g_count(sbj_ix,groi_ix) = actv_g_count(sbj_ix,groi_ix)+1;
                end
                if any(actv_ep_sign==-1)
                    %                 actv_code{roi_ix}    = [actv_code{roi_ix} -1];
                    dact_g_count(sbj_ix,groi_ix) = dact_g_count(sbj_ix,groi_ix)+1;
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
                %             cond_count(roi_ix) = cond_count(roi_ix) + 1;
                if strcmp('stim',event)
                    % Exclude differences after the mean RT for this SBJ
                    if min(mask_chunks(:,1)) < mean_RTs(sbj_ix)
                        cond_g_count(sbj_ix,groi_ix) = cond_g_count(sbj_ix,groi_ix) + 1;
                        both_ix = both_ix + 1;
                    end
                else
                    cond_g_count(sbj_ix,groi_ix) = cond_g_count(sbj_ix,groi_ix) + 1;
                    both_ix = both_ix + 1;
                end
            end
            
            % Check for overlap
            if both_ix==2
                %             both_count(roi_ix) = both_count(roi_ix)+1;
                both_g_count(sbj_ix,groi_ix) = both_g_count(sbj_ix,groi_ix)+1;
            end
        end
    end
    
    clear SBJ SBJ_vars hfa hfa_actv stat einfo actv_ch actv_ch_epochs tmp trial_info
end

%% Plot Percentage of Electrodes Active, Deactive, and Condition Sensitive
% Create and format the plot
fig_name = ['GRP_HFA_ROI_bar_perc_actv_cond_gROI_' event];
figure('Name',fig_name,'units','normalized',...
        'outerposition',[0 0 0.5 0.6],'Visible',fig_vis);
hold on;

% % Create place holder line to initialize axes
% [ax,h1,h2] = plotyy(plt_vars.plt_lim,[find(plot_idx==0,1) find(plot_idx==0,1)],...
%        plt_vars.plt_lim,[find(plot_idx==0,1) find(plot_idx==0,1)]);
% delete(h1);
% delete(h2);

% Compile Data in Plotting Format
bar_vars  = zeros([3 numel(groi_list)]);
for groi_ix = 1:numel(groi_list)
    bar_vars(1,groi_ix) = sum(actv_g_count(:,groi_ix))/sum(elec_g_count(:,groi_ix));
    bar_vars(2,groi_ix) = sum(dact_g_count(:,groi_ix))/sum(elec_g_count(:,groi_ix));
    bar_vars(3,groi_ix) = sum(cond_g_count(:,groi_ix))/sum(elec_g_count(:,groi_ix));
end
scat_vars    = {actv_g_count, dact_g_count, cond_g_count};

% Plot Activations by ROI
b = {};
scat_offsets = linspace(-0.015,0.015,numel(SBJs));
bar_offsets  = linspace(-0.25,0.25,numel(groi_list));   %bar for activation, deactivation, condition
bar_colors   = {'r','b','k'};
for an_ix = 1:3
    % Use "stacked" bars that have empty elements to trick MATLAB into
    % thinking there are multiple elements, which lets me change properties of individual bars
    b{an_ix} = bar(bar_offsets+an_ix,diag(bar_vars(an_ix,:)),0.9,'stacked');
    for groi_ix = 1:numel(groi_list)
        set(b{an_ix}(groi_ix),'FaceColor',groi_colors{groi_ix},'EdgeColor','k');
        % Plot individual subject percentages as scatters on top
        s = scatter(scat_offsets+bar_offsets(groi_ix)+an_ix,...
            scat_vars{an_ix}(:,groi_ix)./elec_g_count(:,groi_ix),50,'k*');
    end
end
if strcmp(event,'stim')
    legend([b{1},s],groi_list{:},'Individuals','Location','northwest');
else
    legend([b{1},s],groi_list{:},'Individuals','Location','northeast');
end

% Plot labels
ax = gca;
% ax.XLabel.String   = 'Time (s)';
% ax.XLabel.FontSize = 14;
ax.XLim    = [0.5 3.5];
ax.XTick   = 1:3;
ax.XTickLabel = {'Activation','Deactivation','Condition Difference'};
ax.XColor  = 'k';

ax.YLabel.String   = '% Electrodes';
ax.YLabel.FontSize = 14;
% ax.YLim            = [0 ymaxs(plot_ix)];
ax.YTick           = ax.YLim(1):0.1:ax.YLim(2);
% ax.YTickLabel      = roi_list;
% ax.YTickLabelRotation = 45;
ax.YColor  = 'k';

ax.Title.String = 'Percentage of Electrodes Showing Significant Effects';
%     ax.Title.String = sprintf('Condition: %.2f; Active: %.2f (+=%.2f;-=%.2f;+-=%.2f); Both: %.2f',...
%         perc_cond,perc_actv,perc_actv_pos,perc_actv_neg,perc_actv_both,perc_both);
ax.Title.FontSize = 16;

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
