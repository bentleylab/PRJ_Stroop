function SBJ10c_HFA_GRP_errbar_ROI_stat_comb(SBJs,proc_id,stat_conds,atlas_id,roi_id,...
                                             gm_thresh,plt_id,plot_out,plot_scat,save_fig,fig_vis,fig_ftype)
% Load HFA analysis results for two stat analyses, plot bars with error
% INPUTS:
%   SBJs [cell array str] - subject IDs to plot
%   proc_id [str] - name of analysis pipeline, used to pick elec file
%   stat_conds [cell array] - {{'stat_id1','an_id1','cond1'},...,{'stat_idn','an_idn','condn'}};
%     stat_id_n [str] - ID of statistical analysis
%     an_id [str] - analysis ID for preprocessing, filtering, etc.
%     cond_n [str] - ID of condition/group to overlap with
%       'actv': red for active, blue for deactive, yellow for both
%       NOPE: 'CI': inc vs. con via ft statistics (not run for all patients!)
%       'RT': correlation with RT (red for significant)
%       'CNI': ANOVA of congruence (red for sig)
%       'pCNI': ANOVA of previous trial congruence (red for sig)
%       'PC': ANOVA of proportion congruence (red for sig)
%   atlas_id [str] - {'DK','Dx','Yeo7','Yeo17'}
%   roi_id [str] - ROI grouping by which to color the atlas ROIs
%       'gROI','mgROI','main3' - general ROIs (lobes or broad regions)
%       'ROI','thryROI','LPFC','MPFC','OFC','INS' - specific ROIs (within these larger regions)
%       'Yeo7','Yeo17' - colored by Yeo networks
%       'tissue','tissueC' - colored by tisseu compartment, e.g., GM vs WM vs OUT
%   gm_thresh [float] - threshold of GM % to include electrode (likely = 0)
%   plt_id [str] - ID of the plotting variables
%   plot_out [0/1] - include electrodes that don't have an atlas label or in hemi?
%   plot_scat [0/1] - plot individual SBJ as asterisks scatter on top of bars
%   save_fig [0/1] - save this figure?
%   fig_vis [str] - visible = 'on'/'off'
%   fig_ftype [str] - file extension for figure saving
% OUTPUTS:
%   Bar chart with % active, % deactivated, % RT correlations, % ANOVA factors

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
for st_ix = 1:numel(stat_conds)
    if numel(stat_conds{st_ix})~=3; error('need stat_id, an_id, cond in each stat_cond'); end
end

% Load all ROI info
[roi_list, roi_colors] = fn_roi_label_styles(roi_id);

% Organize IDs
stat_ids = cell(size(stat_conds)); an_ids = cell(size(stat_conds));
cond_ids = cell(size(stat_conds)); ep_labs = cell(size(stat_conds));
for st_ix = 1:numel(stat_conds)
    stat_ids{st_ix} = stat_conds{st_ix}{1};
    an_ids{st_ix}   = stat_conds{st_ix}{2};
    cond_ids{st_ix} = stat_conds{st_ix}{3};
end

%% Load Results
sig_cnt  = zeros([numel(SBJs) numel(stat_ids) numel(roi_list)]);
elec_cnt = zeros([numel(SBJs) numel(roi_list)]);
for sbj_ix = 1:numel(SBJs)
    SBJ = SBJs{sbj_ix};
    fprintf('================= Processing: %s =================\n',SBJ);
    SBJ_vars_cmd = ['run ' root_dir 'PRJ_Stroop/scripts/SBJ_vars/' SBJ '_vars.m'];
    eval(SBJ_vars_cmd);
    
    %% Prepare elec structs
    % Load elec struct
    elec_fname = [SBJ_vars.dirs.recon,SBJ,'_elec_',proc_id,'_pat_',atlas_id,'_full.mat'];
    if exist([elec_fname(1:end-4) '_' roi_id '.mat'],'file')
        elec_fname = [elec_fname(1:end-4) '_' roi_id '.mat'];
    end
    load(elec_fname);
    
    % Match elecs to atlas ROIs
    if ~isfield(elec,'man_adj')
        if strcmp(atlas_id,{'Yeo17'})
            elec.roi = elec.atlas_lab;
        else
            elec.roi = fn_atlas2roi_labels(elec.atlas_lab,atlas_id,roi_id);
        end
    end
    
    % Remove hemi and/or atlas elecs
    if ~plot_out
        roi_elecs = fn_select_elec_lab_match(elec, 'b', atlas_id, roi_id);
    end
    
    % Get GM probability from tissue labels {'GM','WM','CSF','OUT'}
    if gm_thresh>0
        gm_elecs  = elec.label(elec.tissue_prob(:,1)>gm_thresh);
    else
        gm_elecs = elec.label;
    end
    
    cfgs = []; cfgs.channel = intersect(roi_elecs,gm_elecs);
    elec = fn_select_elec(cfgs,elec);
    
    % Count elecs
    roi_idx = zeros(size(elec.label));
    for roi_ix = 1:numel(roi_list)
        roi_idx(strcmp(elec.roi,roi_list{roi_ix})) = roi_ix;
        elec_cnt(sbj_ix,roi_ix) = sum(roi_idx==roi_ix);
    end
    
    %% Load Stats
    if sum(elec_cnt(sbj_ix,:)) > 0
        for st_ix = 1:numel(stat_conds)
            if strcmp(cond_ids{st_ix},'actv')
                load([SBJ_vars.dirs.proc SBJ '_ROI_' an_ids{st_ix} '_' stat_ids{st_ix} '.mat']);
                ep_labs{st_ix} = st.ep_lab;
                for ch_ix = 1:numel(elec.label)
                    stat_ch_ix = strcmp(actv.label,elec.label{ch_ix});
                    if ~any(stat_ch_ix)
                        error([SBJ ' ' elec.label{ch_ix} ' not found in actv!']);
                    end
                    if any(squeeze(actv.qval(stat_ch_ix,:))<st.alpha)
                        sig_cnt(sbj_ix,st_ix,roi_idx(ch_ix)) = sig_cnt(sbj_ix,st_ix,roi_idx(ch_ix))+1;
                    end
                end
            elseif strcmp(cond_ids{st_ix},'RT')
                error('RT not implemented');
            else    % ANOVA
                load([SBJ_vars.dirs.proc SBJ '_smANOVA_ROI_' stat_ids{st_ix} '_' an_ids{st_ix} '.mat']);
                ep_labs{st_ix} = st.ep_lab;
                for ch_ix = 1:numel(elec.label)
                    stat_ch_ix = strcmp(w2.label,elec.label{ch_ix});
                    if ~any(stat_ch_ix)
                        error([SBJ ' ' elec.label{ch_ix} ' not found in w2!']);
                    end
                    if any(squeeze(w2.qval(strcmp(st.groups,cond_ids{st_ix}),stat_ch_ix,:))<st.alpha)
                        sig_cnt(sbj_ix,st_ix,roi_idx(ch_ix)) = sig_cnt(sbj_ix,st_ix,roi_idx(ch_ix))+1;
                    end
                end
            end
            clear w2 st actv stat_ch_ix
        end
    else
        % Print no ROI match
        fprintf(2,'\t%s has no channels in %s ROIs %s\n',SBJ,atlas_id,roi_id);
    end
    clear SBJ SBJ_vars SBJ_vars_cmd elec
end

%% Plot Percentage of Electrodes Active, Deactive, and Condition Sensitive
eval(['run ' root_dir 'PRJ_Stroop/scripts/plt_vars/' plt_id '_vars.m']);
cond_eps = cell(size(stat_conds));
for st_ix = 1:numel(stat_conds)
    cond_eps{st_ix} = [cond_ids{st_ix} '.' ep_labs{st_ix}];
end
if plot_scat; scat_suffix = '_SBJscat'; else scat_suffix = ''; end;
% Create and format the plot
fig_name = ['GRP_HFA_errbar_GM' num2str(gm_thresh*100) '_' atlas_id '_' roi_id '_' strjoin(cond_eps,'-') scat_suffix];
figure('Name',fig_name,'units','normalized',...
        'outerposition',plt.fig_pos,'Visible',fig_vis);
hold on;

% Compile Data in Plotting Format
scat_vars = cell(size(stat_conds));
bar_data  = zeros([numel(stat_conds) numel(roi_list)]);%+1 CSE
var_data  = zeros([numel(stat_conds) numel(roi_list)]);%+1 CSE
for roi_ix = 1:numel(roi_list)
    for st_ix = 1:numel(stat_conds)
        bar_data(st_ix,roi_ix) = sum(sig_cnt(:,st_ix,roi_ix))/sum(elec_cnt(:,roi_ix));
        var_data(st_ix,roi_ix) = nanstd(sig_cnt(:,st_ix,roi_ix)./elec_cnt(:,roi_ix))/sqrt(numel(SBJs));
        scat_vars{st_ix} = squeeze(sig_cnt(:,st_ix,:));
    end
end

% Plot Activations by ROI
b = cell(size(stat_conds));
for st_ix = 1:numel(stat_conds)
    % Use "stacked" bars that have empty elements to trick MATLAB into
    % thinking there are multiple elements, which lets me change properties of individual bars
    b{st_ix} = bar(plt.bar_offsets+st_ix,diag(bar_data(st_ix,:)),0.9,'stacked');
    for roi_ix = 1:numel(roi_list)
        set(b{st_ix}(roi_ix),'FaceColor',roi_colors{roi_ix},'EdgeColor','k');
        % Add error bar across SBJs
        line([plt.bar_offsets(roi_ix) plt.bar_offsets(roi_ix)]+st_ix,...
            [bar_data(st_ix,roi_ix)+var_data(st_ix,roi_ix) bar_data(st_ix,roi_ix)-var_data(st_ix,roi_ix)],...
            'Color',plt.errbar_color,'LineWidth',plt.errbar_width);
        % Plot individual subject percentages as scatters on top
        if plot_scat
            has_elecs = elec_cnt(:,roi_ix)>0;
            s = scatter(plt.scat_offsets(has_elecs)+plt.bar_offsets(roi_ix)+st_ix,...
                scat_vars{st_ix}(has_elecs,roi_ix)./elec_cnt(has_elecs,roi_ix),plt.scat_sz,plt.scat_style);
        end
    end
end
if plot_scat
    legend([b{1},s],roi_list{:},['Patient (n=' num2str(numel(SBJs)) ')'],'Location',plt.legend_loc);
else
    legend(b{1},roi_list{:},'Location',plt.legend_loc);
end

% Plot labels
ax = gca;
ax.XLim       = [0.5 0.5+numel(stat_conds)];
ax.XTick      = 1:numel(stat_conds);
ax.XColor     = plt.ax_color;
ax.XTickLabel = cond_eps;

ax.YLabel.String   = 'Proportion of Electrodes';
ax.YLabel.FontSize = plt.ylab_sz;
ax.YLim            = plt.ylim;%[0 ymaxs(plot_ix)];%
ax.YTick           = ax.YLim(1):plt.y_tick_step:ax.YLim(2);
ax.YColor          = plt.ax_color;

ax.Title.String = 'Proportion of Electrodes Showing Significant Effects';
ax.Title.FontSize = plt.title_sz;

%% Save figure
if save_fig
    fig_dir = [root_dir 'PRJ_Stroop/results/HFA/GRP_errbar_ROI/' strjoin(stat_ids,'-')...
        '/' strjoin(an_ids,'-') '/' strjoin(cond_ids,'-') '/'];
    if ~exist(fig_dir,'dir')
        [~,~] = mkdir(fig_dir);
    end
    
    fig_fname = [fig_dir fig_name '.' fig_ftype];
    fprintf('Saving %s\n',fig_fname);
    saveas(gcf,fig_fname);
end

end
