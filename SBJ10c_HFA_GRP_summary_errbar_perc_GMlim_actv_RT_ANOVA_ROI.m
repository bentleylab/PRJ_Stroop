function SBJ10c_HFA_GRP_summary_errbar_perc_GMlim_actv_RT_ANOVA_ROI(SBJs,proc_id,an_id,stat_id,actv_id,atlas_id,roi_id,...
                                                            gm_thresh,plt_id,plot_out,plot_scat,save_fig,fig_vis,fig_ftype)
% Load HFA analysis results for active, RT correlation, and ANOVA epochs
%   RT correlation: any significance in stat_lim
%   ANOVA factors: any significance in stat_lim, after FDR correction
% OUTPUTS:
%   Bar chart with % active, % deactivated, % RT correlations, % ANOVA factors
% clear all; %close all;
% fig_ftype = 'png';
label_spacer = 0;
groi_label_spacer = '      ';
if ischar(save_fig); save_fig = str2num(save_fig); end
if ischar(plot_scat); plot_scat = str2num(plot_scat); end

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
% plt_vars_cmd = ['run ' root_dir 'PRJ_Stroop/scripts/plt_vars/' plt_id '_vars.m'];
% eval(plt_vars_cmd);

% Get example condition info
load([root_dir 'PRJ_Stroop/data/' SBJs{1} '/04_proc/' SBJs{1} '_smANOVA_ROI_' stat_id '_' an_id '.mat'],'st');
evnt_lab = st.evnt_lab;
[grp_lab, ~, ~] = fn_group_label_styles(st.model_lab);
rt_grp = st.rt_corr;
if rt_grp
    stat_lab = [{'Active','Deactive','corr(RT)'},grp_lab];
else
    stat_lab = [{'Active','Deactive'},grp_lab];
end

% Load all ROI info
[roi_list, roi_colors] = fn_roi_label_styles(roi_id);

% Set up electrode counts
actv_cnt  = zeros([numel(SBJs) numel(roi_list)]);
dact_cnt  = zeros([numel(SBJs) numel(roi_list)]);
% cse_cnt  = zeros([numel(SBJs) numel(roi_list)]);
grp_cnt   = zeros([numel(SBJs) numel(roi_list) numel(grp_lab)]);
elec_cnt  = zeros([numel(SBJs) numel(roi_list)]);
if rt_grp; rt_cnt = zeros([numel(SBJs) numel(roi_list)]); end

if strcmp(evnt_lab,'S'); mean_RTs = NaN(size(SBJs)); end

%% Load Results
for sbj_ix = 1:numel(SBJs)
    SBJ = SBJs{sbj_ix};
    fprintf('================= Processing: %s =================\n',SBJ);
    % Load variables
    SBJ_vars_cmd = ['run ' root_dir 'PRJ_Stroop/scripts/SBJ_vars/' SBJ '_vars.m'];
    eval(SBJ_vars_cmd);
    
    % Load ANOVA
    load(strcat(SBJ_vars.dirs.proc,SBJ,'_smANOVA_ROI_',stat_id,'_',an_id,'.mat'));
    if ~strcmp(st.evnt_lab,evnt_lab); error('mismatching evnt_lab in w2 stats'); end
    
    % Load actv stats
    actv_fname = strcat(SBJ_vars.dirs.proc,SBJ,'_ROI_',an_id,'_',actv_id,'.mat');
    load(actv_fname,'actv');
    tmp = load(actv_fname,'st');
    if ~strcmp(tmp.st.evnt_lab,evnt_lab); error('mismatching evnt_lab in actv stats'); end
    if tmp.st.alpha~=st.alpha; error('mismatching alpha in actv stats'); end
    
    % Load HFA (to get sign of activation)
    hfa_fname = strcat(SBJ_vars.dirs.proc,SBJ,'_ROI_',an_id,'.mat');
    tmp = load(hfa_fname,'hfa'); hfa = tmp.hfa;
    
%     CSE_filename = strcat(SBJ_vars.dirs.proc,SBJ,'_ROI_',cse_an_id,'_',conditions,'.mat');%an_id
%     tmp = load(CSE_filename,'stat'); cse = tmp.stat;
    
    % Compute mean RT
    load(strcat(SBJ_vars.dirs.events,SBJ,'_trial_info_final.mat'),'trial_info');
    if strcmp(evnt_lab,'S')
        mean_RTs(sbj_ix) = mean(trial_info.response_time);
    end
    
    % Restrict hfa to stat_lim (should obviate some of below comparisons)
    cfg_lim = [];
    cfg_lim.latency = st.stat_lim;
    hfa = ft_selectdata(cfg_lim,hfa);
    if rt_grp; stat = ft_selectdata(cfg_lim,stat); end
%     cse      = ft_selectdata(cfg_lim,cse);
    
    % Confirm channels and time axis are the same
    if ~isempty(setdiff(w2.label,hfa.label));
        error('Different electrodes across hfa analyses!');
    end
    
    %% Load ROI and GM/WM info
    if strcmp(atlas_id,'Yeo7')
        elec_fname = [SBJ_vars.dirs.recon SBJ '_elec_' proc_id '_mni_v_' atlas_id '.mat'];
    else
        elec_fname = [SBJ_vars.dirs.recon SBJ '_elec_' proc_id '_pat_' atlas_id '_full.mat'];
    end
    load(elec_fname);
    
    %!!! eventually load this first as my ground truth for elec_cnt
    % Sort elecs by hfa labels
    cfgs = []; cfgs.channel = hfa.label;
    elec = fn_select_elec(cfgs,elec);
    elec.roi = fn_atlas2roi_labels(elec.atlas_lab,atlas_id,roi_id);
    
    % Exclude elecs not in atlas ROIs
    if ~plot_out
        cfgs = []; cfgs.channel = fn_select_elec_lab_match(elec, 'b', atlas_id, roi_id);
        elec = fn_select_elec(cfgs, elec);
%         hfa = ft_selectdata(cfgs,hfa);
%         w2 = ft_selectdata(cfgs,w2);
%         if rt_grp; stat = ft_selectdata(cfgs,stat); end
    end
    
    % Get GM probability from tissue labels {'GM','WM','CSF','OUT'}
    if gm_thresh>0
        gm_bin  = elec.tissue_prob(:,1)>gm_thresh;
    else
        gm_bin = ones(size(elec.label));
    end
    
    %% Aggregate results per ROI
    for ch_ix = 1:numel(elec.label)
        % If elec matches roi_list, get stats
        if any(strcmp(elec.roi{ch_ix},roi_list)) && gm_bin(ch_ix)
            roi_ix = find(strcmp(elec.roi{ch_ix},roi_list));
            elec_cnt(sbj_ix,roi_ix) = elec_cnt(sbj_ix,roi_ix)+1;
            
%             if any(cse.mask(ch_ix,1,:))
%                 cse_cnt(sbj_ix,roi_ix) = cse_cnt(sbj_ix,roi_ix)+1;
            % Check if active, get epochs
            actv_ch_ix = strcmp(actv.label,elec.label{ch_ix});
            if any(squeeze(actv.qval(actv_ch_ix,:))<st.alpha)
                % Find significant epoch indices
                actv_idx = squeeze(actv.qval(actv_ch_ix,:))<st.alpha;
                actv_epochs = fn_find_chunks(actv_idx);
                actv_epochs(actv_idx(actv_epochs(:,1))==0,:) = [];
                
                % Find sign of (de)activation of all epochs
                actv_ep_sign = NaN([1 size(actv_epochs,1)]);
                for ep_ix = 1:size(actv_epochs,1)
                    if mean(actv.avg(actv_ch_ix,actv_epochs(ep_ix,1):actv_epochs(ep_ix,2)),2)>=0
                        actv_ep_sign(ep_ix) = 1;
                    else
                        actv_ep_sign(ep_ix) = -1;
                    end
                end
                
                % Code significance and sign: 0=none, -1=deactive, 1=active, 2=both
                if any(actv_ep_sign==1)
                    actv_cnt(sbj_ix,roi_ix) = actv_cnt(sbj_ix,roi_ix)+1;
                end
                if any(actv_ep_sign==-1)
                    dact_cnt(sbj_ix,roi_ix) = dact_cnt(sbj_ix,roi_ix)+1;
                end
            end
            
            % Check for RT correlations, get epochs
            if rt_grp
                stat_ch_ix = strcmp(stat.label,elec.label{ch_ix});
                if sum(squeeze(stat.mask(stat_ch_ix,1,:)))>0
                    rt_cnt(sbj_ix,roi_ix) = rt_cnt(sbj_ix,roi_ix) + 1;
                end
            end
            
            % Check for ANOVA group effects
            grp_ch_ix = strcmp(w2.label,elec.label{ch_ix});
            for grp_ix = 1:numel(grp_lab)
                if any(squeeze(w2.qval(grp_ix,grp_ch_ix,:))<st.alpha)
                    grp_cnt(sbj_ix,roi_ix,grp_ix) = grp_cnt(sbj_ix,roi_ix,grp_ix)+1;
                end
            end
        end
    end
    clear SBJ SBJ_vars w2 stat elec hfa actv_ch actv_ch_epochs tmp trial_info st %cse
end

%% Plot Percentage of Electrodes Active, Deactive, and Condition Sensitive
if plot_scat
    scat_suffix = '_SBJscat';
else
    scat_suffix = '';
end
% Create and format the plot
fig_name = ['GRP_HFA_errbar_perc_GMlim' num2str(gm_thresh*100) '_' stat_id '_' actv_id '_' roi_id '_' evnt_lab scat_suffix];
figure('Name',fig_name,'units','normalized',...
        'outerposition',[0 0 0.5 0.6],'Visible',fig_vis);
hold on;

% % Create place holder line to initialize axes
% [ax,h1,h2] = plotyy(plt_vars.plt_lim,[find(plot_idx==0,1) find(plot_idx==0,1)],...
%        plt_vars.plt_lim,[find(plot_idx==0,1) find(plot_idx==0,1)]);
% delete(h1);
% delete(h2);

% Compile Data in Plotting Format
scat_vars = {actv_cnt, dact_cnt};%cse_cnt, 
if rt_grp; scat_vars = [scat_vars {rt_cnt}]; end
bar_data  = zeros([numel(stat_lab) numel(roi_list)]);%+1 CSE
var_data  = zeros([numel(stat_lab) numel(roi_list)]);%+1 CSE
for roi_ix = 1:numel(roi_list)
%     bar_vars(1,roi_ix) = sum(cI_cnt(:,roi_ix))/sum(elec_cnt(:,roi_ix));
%     bar_data(1,roi_ix) = sum(cse_cnt(:,roi_ix))/sum(elec_cnt(:,roi_ix));
%     var_data(1,roi_ix) = std(cse_cnt(:,roi_ix)./elec_cnt(:,roi_ix))/sqrt(numel(SBJs));
    bar_data(1,roi_ix) = sum(actv_cnt(:,roi_ix))/sum(elec_cnt(:,roi_ix));
    var_data(1,roi_ix) = nanstd(actv_cnt(:,roi_ix)./elec_cnt(:,roi_ix))/sqrt(numel(SBJs));
    bar_data(2,roi_ix) = sum(dact_cnt(:,roi_ix))/sum(elec_cnt(:,roi_ix));
    var_data(2,roi_ix) = nanstd(dact_cnt(:,roi_ix)./elec_cnt(:,roi_ix))/sqrt(numel(SBJs));
    if rt_grp
        bar_data(3,roi_ix) = sum(rt_cnt(:,roi_ix))/sum(elec_cnt(:,roi_ix));
        var_data(3,roi_ix) = nanstd(rt_cnt(:,roi_ix)./elec_cnt(:,roi_ix))/sqrt(numel(SBJs));
    end
    for grp_ix = 1:numel(grp_lab)
        bar_data(2+rt_grp+grp_ix,roi_ix) = sum(grp_cnt(:,roi_ix,grp_ix))/sum(elec_cnt(:,roi_ix));
        var_data(2+rt_grp+grp_ix,roi_ix) = nanstd(grp_cnt(:,roi_ix,grp_ix)./elec_cnt(:,roi_ix))/sqrt(numel(SBJs));
        scat_vars{2+rt_grp+grp_ix} = squeeze(grp_cnt(:,:,grp_ix));
    end
end

% Plot Activations by ROI
b = {};
scat_offsets = linspace(-0.015,0.015,numel(SBJs));
bar_offsets  = linspace(-0.25,0.25,numel(roi_list));   %bar for activation, deactivation, condition
for cond_ix = 1:2+rt_grp+numel(grp_lab)
    % Use "stacked" bars that have empty elements to trick MATLAB into
    % thinking there are multiple elements, which lets me change properties of individual bars
    b{cond_ix} = bar(bar_offsets+cond_ix,diag(bar_data(cond_ix,:)),0.9,'stacked');
    for roi_ix = 1:numel(roi_list)
        set(b{cond_ix}(roi_ix),'FaceColor',roi_colors{roi_ix},'EdgeColor','k');
        % Add error bar across SBJs
        line([bar_offsets(roi_ix) bar_offsets(roi_ix)]+cond_ix,...
            [bar_data(cond_ix,roi_ix)+var_data(cond_ix,roi_ix) bar_data(cond_ix,roi_ix)-var_data(cond_ix,roi_ix)],...
            'Color','k','LineWidth',1.5);
        % Plot individual subject percentages as scatters on top
        if plot_scat
            has_elecs = squeeze(elec_cnt(:,roi_ix)~=0);
            s = scatter(scat_offsets(has_elecs)+bar_offsets(roi_ix)+cond_ix,...
                scat_vars{cond_ix}(has_elecs,roi_ix)./elec_cnt(has_elecs,roi_ix),50,'k*');
        end
    end
end
% if strcmp(st.evnt_lab,'S')
%     leg_loc = 'northwest';
% else
%     leg_loc = 'northeast';
% end
leg_loc = 'best';%'northeast';
if plot_scat
    legend([b{1},s],roi_list{:},'Individuals','Location',leg_loc);
else
    legend(b{1},roi_list{:},'Location',leg_loc);
end

% Plot labels
ax = gca;
% ax.XLabel.String   = 'Time (s)';
% ax.XLabel.FontSize = 14;
ax.XLim    = [0.5 0.5+numel(stat_lab)];
ax.XTick   = 1:numel(stat_lab);
ax.XColor  = 'k';
ax.XTickLabel = stat_lab;%'CSE'

ax.YLabel.String   = 'Proportion of Electrodes';
ax.YLabel.FontSize = 14;
% ax.YLim            = [0 0.8];%[0 ymaxs(plot_ix)];%
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
    fig_dir = [root_dir 'PRJ_Stroop/results/HFA/GRP_summary_errbar_ROI/' ...
        stat_id '-' actv_id '/' an_id '/' roi_id '/'];
    if ~exist(fig_dir,'dir')
        [~,~] = mkdir(fig_dir);
    end
    
    fig_fname = [fig_dir fig_name '.' fig_ftype];
    fprintf('Saving %s\n',fig_fname);
    saveas(gcf,fig_fname);
end

end
