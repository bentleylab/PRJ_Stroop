function SBJ10c_HFA_GRP_summary_errbar_perc_GMlim_CSE_RT_ANOVA_ROI(SBJs,stat_id,pipeline_id,an_id,roi_id,...
                                                            atlas_id,gm_thresh,plt_id,plot_out,save_fig,fig_vis,fig_filetype)
% Load HFA analysis results for active, RT correlation, and ANOVA epochs
%   RT correlation: any significance in stat_lim
%   ANOVA factors: any significance in stat_lim, after FDR correction
% OUTPUTS:
%   Bar chart with % active, % deactivated, % RT correlations, % ANOVA factors
% clear all; %close all;
% fig_filetype = 'png';
label_spacer = 0;
groi_label_spacer = '      ';
if ischar(save_fig); save_fig = str2num(save_fig); end
% if isnumeric(actv_win); actv_win = num2str(actv_win); end

%% Data Preparation
% Set up paths
[root_dir, app_dir] = fn_get_root_dir(); ft_dir = [app_dir 'fieldtrip/'];
addpath([root_dir 'PRJ_Stroop/scripts/']);
addpath([root_dir 'PRJ_Stroop/scripts/utils/']);
addpath(ft_dir);
ft_defaults

%% Prep variables
an_vars_cmd = ['run ' root_dir 'PRJ_Stroop/scripts/an_vars/' an_id '_vars.m'];
eval(an_vars_cmd);
plt_vars_cmd = ['run ' root_dir 'PRJ_Stroop/scripts/plt_vars/' plt_id '_vars.m'];
eval(plt_vars_cmd);
stat_vars_cmd = ['run ' root_dir 'PRJ_Stroop/scripts/stat_vars/' stat_id '_vars.m'];
eval(stat_vars_cmd);

% Get condition info
[grp_lab, grp_colors, grp_style] = fn_group_label_styles(model_lab);
% if rt_correlation
[rt_lab, rt_color, rt_style]     = fn_group_label_styles('RT');
% end
conditions = 'CSE';
[cond_lab, cond_colors, cond_style] = fn_condition_label_styles(conditions);
cse_an_id = strrep(an_id,'_sm0_','_sm10_');

% Load all ROI info
[roi_list, roi_colors] = fn_roi_label_styles(roi_id);

% Set up electrode counts
cse_cnt  = zeros([numel(SBJs) numel(roi_list)]);
rt_cnt    = zeros([numel(SBJs) numel(roi_list)]);
grp_cnt   = zeros([numel(SBJs) numel(roi_list) numel(grp_lab)]);
elec_cnt  = zeros([numel(SBJs) numel(roi_list)]);

mean_RTs = NaN(size(SBJs));

%% Load Results
for sbj_ix = 1:numel(SBJs)
    SBJ = SBJs{sbj_ix};
    fprintf('================= Processing: %s =================\n',SBJ);
    % Load variables
    SBJ_vars_cmd = ['run ' root_dir 'PRJ_Stroop/scripts/SBJ_vars/' SBJ '_vars.m'];
    eval(SBJ_vars_cmd);
    
    % Compute mean RT
    if strcmp(event_type,'stim')
        load(strcat(SBJ_vars.dirs.events,SBJ,'_trial_info_final.mat'),'trial_info');
        mean_RTs(sbj_ix) = mean(trial_info.response_time); % converts sec to ms
    end
    % Load data
    load(strcat(SBJ_vars.dirs.proc,SBJ,'_ANOVA_ROI_',stat_id,'_',an_id,'.mat'));
    CSE_filename = strcat(SBJ_vars.dirs.proc,SBJ,'_ROI_',an_id,'_',conditions,'.mat');%cse_an_id
    tmp = load(CSE_filename,'stat'); cse = tmp.stat;
    
    %% Load ROI and GM/WM info
    if strcmp(atlas_id,'Yeo7')
        elec_fname = [SBJ_vars.dirs.recon SBJ '_elec_' pipeline_id '_mni_v_' atlas_id '.mat'];
    else
        elec_fname = [SBJ_vars.dirs.recon SBJ '_elec_' pipeline_id '_pat_' atlas_id '_full.mat'];
    end
    load(elec_fname);
    
    % Sort elecs by stat labels
    cfgs = []; cfgs.channel = stat.label;
    elec = fn_select_elec(cfgs,elec);
    elec.roi = fn_atlas2roi_labels(elec.atlas_lab,atlas_id,roi_id);
    
    % Exclude elecs not in atlas ROIs
    if ~plot_out
        cfgs = []; cfgs.channel = fn_select_elec_lab_match(elec, 'b', atlas_id, []);
        elec = fn_select_elec(cfgs, elec);
        stat = ft_selectdata(cfgs,stat);
        w2 = ft_selectdata(cfgs,w2);
    end
    
    % Get GM probability from tissue labels {'GM','WM','CSF','OUT'}
    if gm_thresh>0
        gm_bin  = elec.tissue_prob(:,1)>gm_thresh;
    else
        gm_bin = ones(size(elec.label));
    end
    
    %% Process parameters
    %!!! is this the best way to do this??? Maybe not...
    sample_rate = (numel(stat.time)-1)/(stat.time(end)-stat.time(1));
    
    % Restrict hfa_actv to stat_lim (should obviate some of below comparisons)
    cfg_lim = [];
    cfg_lim.latency = stat_lim;
    cse      = ft_selectdata(cfg_lim,cse);
    stat     = ft_selectdata(cfg_lim,stat);
    
    % FDR correct pvalues for ANOVA
    qvals = NaN(size(w2.pval));
    for ch_ix = 1:numel(stat.label)
        [~, ~, ~, qvals(:,ch_ix,:)] = fdr_bh(squeeze(w2.pval(:,ch_ix,:)));%,0.05,'pdep','yes');
    end
    
    % Confirm channels and time axis are the same
    %   All the rounding and such is because some stupid rounding errors...
    same_start = round(cse.time(1)*uint8(sample_rate))==round(stat.time(1)*uint8(sample_rate));
    same_end   = round(cse.time(end)*uint8(sample_rate))==round(stat.time(end)*uint8(sample_rate));
    same_numel = size(cse.time,2)==size(stat.time,2);
    if ~same_start || ~same_end || ~same_numel
        error('time axes are not the same across hfa analyses!');
    end
    if ~isempty(setdiff(stat.label,cse.label));
        error('Different electrodes across hfa analyses!');
    end
    
    %% Aggregate results per ROI
    for ch_ix = 1:numel(stat.label)
%         both_ix = 0;
        
        % If elec matches roi_list, get stats
        if any(strcmp(elec.roi{ch_ix},roi_list)) && gm_bin(ch_ix)
            roi_ix = find(strcmp(elec.roi{ch_ix},roi_list));
            elec_cnt(sbj_ix,roi_ix) = elec_cnt(sbj_ix,roi_ix)+1;
            
            % Check if active, get epochs
            if any(cse.mask(ch_ix,1,:))
                cse_cnt(sbj_ix,roi_ix) = cse_cnt(sbj_ix,roi_ix)+1;
%                 % Find significant epoch indices
%                 actv_epochs = actv_ch_epochs{strcmp(stat.label{ch_ix},actv_ch)};
%                 
%                 % Toss late epochs if S-locked
%                 if strcmp('stim',event) %!!! isn't this backwards????
%                     actv_epochs(actv_epochs(:,1)<mean_RTs(sbj_ix),:) = [];
%                 end
                
%                 % Find sign of (de)activation
%                 actv_ep_sign = NaN([1 size(actv_epochs,1)]);
%                 sig_chunk_ix = NaN([1 2]);
%                 for ep_ix = 1:size(actv_epochs,1)
%                     sig_chunk_ix = [find(cse.time==actv_epochs(ep_ix,1))...
%                         find(cse.time==actv_epochs(ep_ix,2))];
%                     % Report sign
%                     if 0<=squeeze(mean(mean(cse.powspctrm(:,ch_ix,1,sig_chunk_ix(1):sig_chunk_ix(2)),1),4))
%                         actv_ep_sign(ep_ix) = 1;
%                     else
%                         actv_ep_sign(ep_ix) = -1;
%                     end
%                 end
                
%                 % Code significance and sign: 0=none, -1=deactive, 1=active, 2=both
%                 if any(actv_ep_sign==1)
%                     %                 actv_code{roi_ix}    = [actv_code{roi_ix} 1];
%                     cI_cnt(sbj_ix,roi_ix) = cI_cnt(sbj_ix,roi_ix)+1;
%                 end
%                 if any(actv_ep_sign==-1)
%                     %                 actv_code{roi_ix}    = [actv_code{roi_ix} -1];
%                     iI_cnt(sbj_ix,roi_ix) = iI_cnt(sbj_ix,roi_ix)+1;
%                 end
%                 both_ix = both_ix + 1;
            end
            
            % Check for RT correlations, get epochs
            if sum(squeeze(stat.mask(ch_ix,1,:)))>0
%                 mask_chunks = fn_find_chunks(squeeze(stat.mask(ch_ix,1,:)));
%                 mask_chunks(squeeze(stat.mask(ch_ix,1,mask_chunks(:,1)))==0,:) = [];
%                 % Convert to time
%                 for ep_ix = 1:size(mask_chunks,1)
%                     mask_chunks(ep_ix,1) = stat.time(mask_chunks(ep_ix,1));
%                     mask_chunks(ep_ix,2) = stat.time(mask_chunks(ep_ix,2));
%                 end
%                 %             cond_count(roi_ix) = cond_count(roi_ix) + 1;
%                 if strcmp('stim',event)
%                     % Exclude differences after the mean RT for this SBJ
%                     if min(mask_chunks(:,1)) < mean_RTs(sbj_ix)
%                         cond_g_count(sbj_ix,groi_ix) = cond_g_count(sbj_ix,groi_ix) + 1;
%                         both_ix = both_ix + 1;
%                     end
%                 else
                rt_cnt(sbj_ix,roi_ix) = rt_cnt(sbj_ix,roi_ix) + 1;
%                 end
            end
            
            % Check for ANOVA group effects
            for grp_ix = 1:numel(grp_lab)
                if any(squeeze(qvals(grp_ix,ch_ix,:))<0.05)
                    grp_cnt(sbj_ix,roi_ix,grp_ix) = grp_cnt(sbj_ix,roi_ix,grp_ix)+1;
                end
            end
        end
    end
    clear SBJ SBJ_vars w2 stat einfo cse tmp trial_info qvals
end

%% Plot Percentage of Electrodes Active, Deactive, and Condition Sensitive
% Create and format the plot
fig_name = ['GRP_HFA_errbar_perc_GMlim' num2str(gm_thresh*100) '_CSE_' stat_id '_' roi_id '_' event_type];
figure('Name',fig_name,'units','normalized',...
        'outerposition',[0 0 0.5 0.6],'Visible',fig_vis);
hold on;

% % Create place holder line to initialize axes
% [ax,h1,h2] = plotyy(plt_vars.plt_lim,[find(plot_idx==0,1) find(plot_idx==0,1)],...
%        plt_vars.plt_lim,[find(plot_idx==0,1) find(plot_idx==0,1)]);
% delete(h1);
% delete(h2);

% Compile Data in Plotting Format
% scat_vars = {cse_cnt, rt_cnt};
bar_data  = zeros([2+numel(grp_lab) numel(roi_list)]);
var_data  = zeros([2+numel(grp_lab) numel(roi_list)]);
for roi_ix = 1:numel(roi_list)
%     bar_vars(1,roi_ix) = sum(cI_cnt(:,roi_ix))/sum(elec_cnt(:,roi_ix));
    bar_data(1,roi_ix) = sum(cse_cnt(:,roi_ix))/sum(elec_cnt(:,roi_ix));
    var_data(1,roi_ix) = std(cse_cnt(:,roi_ix)./elec_cnt(:,roi_ix))/sqrt(numel(SBJs));
    bar_data(2,roi_ix) = sum(rt_cnt(:,roi_ix))/sum(elec_cnt(:,roi_ix));
    var_data(2,roi_ix) = nanstd(rt_cnt(:,roi_ix)./elec_cnt(:,roi_ix))/sqrt(numel(SBJs));
    for grp_ix = 1:numel(grp_lab)
        bar_data(2+grp_ix,roi_ix) = sum(grp_cnt(:,roi_ix,grp_ix))/sum(elec_cnt(:,roi_ix));
        var_data(2+grp_ix,roi_ix) = std(grp_cnt(:,roi_ix,grp_ix)./elec_cnt(:,roi_ix))/sqrt(numel(SBJs));
%         scat_vars{2+grp_ix} = squeeze(grp_cnt(:,:,grp_ix));
    end
end

% Plot Activations by ROI
b = {};
% scat_offsets = linspace(-0.015,0.015,numel(SBJs));
bar_offsets  = linspace(-0.25,0.25,numel(roi_list));   %bar for activation, deactivation, condition
for cond_ix = 1:2+numel(grp_lab)
    % Use "stacked" bars that have empty elements to trick MATLAB into
    % thinking there are multiple elements, which lets me change properties of individual bars
    b{cond_ix} = bar(bar_offsets+cond_ix,diag(bar_data(cond_ix,:)),0.9,'stacked');
    for roi_ix = 1:numel(roi_list)
        set(b{cond_ix}(roi_ix),'FaceColor',roi_colors{roi_ix},'EdgeColor','k');
        % Plot individual subject percentages as scatters on top
        line([bar_offsets(roi_ix) bar_offsets(roi_ix)]+cond_ix,...
            [bar_data(cond_ix,roi_ix)+var_data(cond_ix,roi_ix) bar_data(cond_ix,roi_ix)-var_data(cond_ix,roi_ix)],...
            'Color','k','LineWidth',1.5);
%         has_elecs = squeeze(elec_cnt(:,roi_ix)~=0);
%         s = scatter(scat_offsets(has_elecs)+bar_offsets(roi_ix)+cond_ix,...
%             scat_vars{cond_ix}(has_elecs,roi_ix)./elec_cnt(has_elecs,roi_ix),50,'k*');
    end
end
if strcmp(event_type,'stim')
    legend(b{1},roi_list{:},'Location','northwest');
else
    legend(b{1},roi_list{:},'Location','northeast');
end

% Plot labels
ax = gca;
% ax.XLabel.String   = 'Time (s)';
% ax.XLabel.FontSize = 14;
ax.XLim    = [0.5 2.5+numel(grp_lab)];
ax.XTick   = 1:2+numel(grp_lab);
ax.XColor  = 'k';
ax.XTickLabel = {'CSE','Corr(RT)',grp_lab{:}};

ax.YLabel.String   = 'Proportion of Electrodes';
ax.YLabel.FontSize = 14;
ax.YLim            = [0 0.4];
ax.YTick           = ax.YLim(1):0.1:ax.YLim(2);
% ax.YTickLabel      = roi_list;
% ax.YTickLabelRotation = 45;
ax.YColor  = 'k';

ax.Title.String = 'Proportion of Electrodes Showing Significant Effects';
%     ax.Title.String = sprintf('Condition: %.2f; Active: %.2f (+=%.2f;-=%.2f;+-=%.2f); Both: %.2f',...
%         perc_cond,perc_actv,perc_actv_pos,perc_actv_neg,perc_actv_both,perc_both);
ax.Title.FontSize = 16;

%% Save figure
if save_fig
    fig_dir = [root_dir 'PRJ_Stroop/results/HFA/GRP_summary_errbar_ROI/CSE_'...
        stat_id '_' roi_id '/' an_id '/'];
    if ~exist(fig_dir,'dir')
        [~,~] = mkdir(fig_dir);
    end
    
    fig_filename = [fig_dir fig_name '.' fig_filetype];
    fprintf('Saving %s\n',fig_filename);
    saveas(gcf,fig_filename);
    %eval(['export_fig ' fig_filename]);
end

end
