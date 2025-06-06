function SBJ10c_HFA_GRPavg_onsets_ROI_normRTout_RT_ANOVA(SBJs,stat_id,proc_id,an_id,roi_id,...
                                                    atlas_id,gm_thresh,z_thresh,plt_id,save_fig,fig_vis,fig_ftype)
% Load HFA analysis results for active and condition-differentiating
%   epochs, plot a summary of those time period per electrode
% Normalize all onset times by mean(RT)
% INPUTS:
%   plt_vars.grp_metric [str] - {'avg','mdn','all'}
%       mean/median will compute that metric within each SBJ (variance is across SBJs)
%       all- all electrode onsets are aggregated as if from the same SBJ
% clear all; %close all;
% fig_ftype = 'png';
label_spacer = 0;
groi_label_spacer = '      ';
if ischar(save_fig); save_fig = str2num(save_fig); end
% if isnumeric(actv_win); actv_win = num2str(actv_win); end

%% Data Preparation
% Set up paths
[root_dir, app_dir] = fn_get_root_dir(); ft_dir = [app_dir 'fieldtrip/'];
addpath([root_dir 'PRJ_Stroop/scripts/']);
addpath([root_dir 'PRJ_Stroop/scripts/utils/']);

%% Prep variables
an_vars_cmd = ['run ' root_dir 'PRJ_Stroop/scripts/an_vars/' an_id '_vars.m'];
eval(an_vars_cmd);
plt_vars_cmd = ['run ' root_dir 'PRJ_Stroop/scripts/plt_vars/' plt_id '_vars.m'];
eval(plt_vars_cmd);

% Get condition info
load([root_dir,'PRJ_Stroop/data/',SBJs{1},'/04_proc/',SBJs{1},'_smANOVA_ROI_',stat_id,'_',an_id,'.mat'],'st');

% Get event timing
mean_RTs = zeros(size(SBJs));

% Load all ROI info
[roi_list, roi_colors] = fn_roi_label_styles(roi_id);
if any(strcmp(roi_id,{'mgROI','gROI','main3','lat','deep','gPFC'}))
    roi_field = 'gROI';
else
    roi_field = 'ROI';
end
% if any(strcmp(roi_list,'TMP'))
%     % Turn the yellow to gold for visibility
%     roi_colors{strcmp(roi_list,'TMP')}(2) = 0.8;
% end
if any(strcmp(atlas_id,{'DK','Dx'}))
    view_space = 'pat';
elseif any(strcmp(atlas_id,{'Yeo7','Yeo17'}))
    view_space = 'mni_v';
else
    error(['Unknown atlas_id: ' atlas_id]);
end
if strcmp(plt_vars.grp_metric,'all')
    SBJ_colors = distinguishable_colors(numel(SBJs));
end

%% Load Results
% Set up onset counts
all_onsets  = cell([numel(SBJs) numel(roi_list) numel(st.groups)]);
min_rt      = nan(size(SBJs));
for sbj_ix = 1:numel(SBJs)
    SBJ = SBJs{sbj_ix};
    % Load variables
    SBJ_vars_cmd = ['run ' root_dir 'PRJ_Stroop/scripts/SBJ_vars/' SBJ '_vars.m'];
    eval(SBJ_vars_cmd);
    % Load data
    load(strcat(SBJ_vars.dirs.proc,SBJ,'_smANOVA_ROI_',stat_id,'_',an_id,'.mat'));
    
    % Compute mean RT
    load(strcat(SBJ_vars.dirs.events,SBJ,'_trial_info_final.mat'),'trial_info');
    if ~isempty(strfind(st.ep_lab,'S'))
        % Compute mean RT
        mean_RTs(sbj_ix) = mean(trial_info.response_time);
    end
        
    %% Load ROI and GM/WM info
    elec_fname = [SBJ_vars.dirs.recon SBJ '_elec_' proc_id '_' view_space '_' atlas_id '_final.mat'];
    load(elec_fname);
    
    % Sort elecs by stat labels
    if numel(elec.label)~=numel(w2.label) || ~all(strcmp(elec.label,w2.label))
        error('elec w2 label mismatch');
        cfgs = []; cfgs.channel = stat.label;
        elec = fn_select_elec(cfgs,elec);
    end
    
%     % Get GM probability from tissue labels {'GM','WM','CSF','OUT'}
%     gm_bin  = elec.tissue_prob(:,1)>gm_thresh;
        
    %% Aggregate results per ROI
    for ch_ix = 1:numel(w2.label)
        % If elec matches roi_list and is in GM, get stats
        if any(strcmp(elec.(roi_field){ch_ix},roi_list)) && w2.max_hfa_z(ch_ix)>=z_thresh% && gm_bin(ch_ix)
            roi_ix = find(strcmp(elec.(roi_field){ch_ix},roi_list));
            % Get ANOVA group onsets
            for grp_ix = 1:numel(st.groups)
                if any(squeeze(w2.qval(grp_ix,ch_ix,:))<=st.alpha)
                    sig_onsets = w2.win_lim_s(squeeze(w2.qval(grp_ix,ch_ix,:))<=st.alpha,1);
                    if strcmp(st.evnt_lab,'R')
                        all_onsets{sbj_ix,roi_ix,grp_ix} = ...
                            [all_onsets{sbj_ix,roi_ix,grp_ix} sig_onsets(1)];
                    elseif ~isempty(strfind(st.ep_lab,'S')) && (sig_onsets(1)<mean_RTs(sbj_ix))
                        all_onsets{sbj_ix,roi_ix,grp_ix} = ...
                            [all_onsets{sbj_ix,roi_ix,grp_ix} sig_onsets(1)];
                    else
                        error('st.ep_lab didnt match R or S*');
                    end
                end
            end
            
%             % Get RT correlation onset
%             if sum(squeeze(stat.mask(ch_ix,1,:)))>0
%                 mask_chunks = fn_find_chunks(squeeze(stat.mask(ch_ix,1,:)));
%                 mask_chunks(squeeze(stat.mask(ch_ix,1,mask_chunks(:,1)))==0,:) = [];
%                 % Convert the first onset of significance to time
%                 onset_time = stat.time(mask_chunks(1,1));
%                 % Exclude differences after the mean RT for this SBJ
%                 if strcmp(st.evnt_lab,'R')
%                     all_onsets{sbj_ix,roi_ix,numel(grp_lab)+1} = ...
%                         [all_onsets{sbj_ix,roi_ix,numel(grp_lab)+1} onset_time];
%                 elseif strcmp(st.evnt_lab,'S') && (onset_time<mean_RTs(sbj_ix))
%                     all_onsets{sbj_ix,roi_ix,numel(grp_lab)+1} = ...
%                         [all_onsets{sbj_ix,roi_ix,numel(grp_lab)+1} onset_time];
%                 end
%             end
        end
    end
    
    % Normalize all onset times by mean reaction time
    if ~isempty(strfind(st.evnt_lab,'S'))
        y_label = 'Time (% mean RT)';
        min_rt(sbj_ix) = min(trial_info.response_time)/mean_RTs(sbj_ix);
        for roi_ix = 1:size(all_onsets,2)
            for cond_ix = 1:size(all_onsets,3)
                all_onsets{sbj_ix,roi_ix,cond_ix} = all_onsets{sbj_ix,roi_ix,cond_ix}./mean_RTs(sbj_ix);
            end
        end
    else
        y_label = 'Time (sec)';
    end
    clear SBJ SBJ_vars hfa elec w2 trial_info
end

%% Aggregate/Process onsets per gROI
% Format as struct to fit violinplot
plot_onsets    = cell([numel(st.groups) 1]);
plot_onset_sbj = cell([numel(st.groups) 1]);
good_roi_map   = cell([numel(st.groups) 1]);
for cond_ix = 1:numel(st.groups)
    plot_onsets{cond_ix}  = {};
    plot_onset_sbj{cond_ix}  = {};
    good_roi_map{cond_ix} = zeros(size(roi_list));
    for roi_ix = 1:numel(roi_list)
        if strcmp(plt_vars.grp_metric,'all')
            plot_onsets{cond_ix}{roi_ix} = [all_onsets{:,roi_ix,cond_ix}]';
            plot_onset_sbj{cond_ix}{roi_ix} = [];
        else
            plot_onsets{cond_ix}{roi_ix} = [];
        end
        for sbj_ix = 1:numel(SBJs)
            % Aggregate onsets per ROI within each SBJ
            if strcmp(plt_vars.grp_metric,'all')
                plot_onset_sbj{cond_ix}{roi_ix} = [plot_onset_sbj{cond_ix}{roi_ix};...
                                                    repmat(sbj_ix,size(all_onsets{sbj_ix,roi_ix,cond_ix}))'];
            elseif strcmp(plt_vars.grp_metric,'mdn')
                plot_onsets{cond_ix}{roi_ix} = [plot_onsets{cond_ix}{roi_ix}; nanmedian(all_onsets{sbj_ix,roi_ix,cond_ix})];
            elseif strcmp(plt_vars.grp_metric,'avg')
                plot_onsets{cond_ix}{roi_ix} = [plot_onsets{cond_ix}{roi_ix}; nanmean(all_onsets{sbj_ix,roi_ix,cond_ix})];
            else
                error(['Unknown plt_vars.grp_metric: ' plt_vars.grp_metric]);
            end
            % Report results in text
            %         fprintf('%s , %s: %f (N=%i)\n',SBJs{sbj_ix},groi_list{groi_ix},...
            %             median_onsets(sbj_ix,groi_ix),numel(cond_g_onsets{sbj_ix,groi_ix}));
            %         disp(cond_g_onsets{sbj_ix,groi_ix});
            %         fprintf('\n');
        end
    end
    plot_onsets{cond_ix} = padcat(plot_onsets{cond_ix}{:});
    [~,good_roi_map{cond_ix}] = find(~all(isnan(plot_onsets{cond_ix}),1)==1);
    if strcmp(plt_vars.grp_metric,'all')
        plot_onset_sbj{cond_ix} = padcat(plot_onset_sbj{cond_ix}{:});
    end
end
% save([root_dir 'PRJ_Stroop/data/tmp_onsets_sig_ch.mat'],'-v7.3','sig_ch','plot_onsets','all_onsets');

%% Plot GROI Results
for cond_ix = 1:numel(st.groups)
    % Create and format the plot
    fig_name = ['GRP' plt_vars.grp_metric '_HFA_onsets_' st.groups{cond_ix} '_' roi_id '_' st.ep_lab...
        '_GM' num2str(gm_thresh) '_z' num2str(z_thresh) '_normRTout'];
    figure('Name',fig_name,'units','normalized',...
        'outerposition',[0 0 0.5 0.6],'Visible',fig_vis);
    fprintf('Printing %i onsets across %i groups in %s\n',sum(sum(~isnan(plot_onsets{cond_ix}(:,good_roi_map{cond_ix})))),...
                    numel(good_roi_map{cond_ix}),fig_name);
    
    violins = violinplot(plot_onsets{cond_ix}(:,good_roi_map{cond_ix}),roi_list(good_roi_map{cond_ix}),...
                            'ViolinAlpha',0.3);
                        
    % Adjust plot propeties
    for roi_ix = 1:numel(good_roi_map{cond_ix})
        if strcmp(plt_vars.violin_scat_colors,'SBJ')
            % Change scatter colors to mark SBJ
            violins(roi_ix).ViolinColor = [0.8 0.8 0.8];
            violins(roi_ix).BoxPlot.FaceColor = roi_colors{good_roi_map{cond_ix}(roi_ix)};
            violins(roi_ix).EdgeColor = roi_colors{good_roi_map{cond_ix}(roi_ix)};
            if sum(~isnan(plot_onsets{cond_ix}(:,good_roi_map{cond_ix}(roi_ix))))>1
                scat_colors = zeros([numel(violins(roi_ix).ScatterPlot.XData) 3]);
                for sbj_ix = 1:numel(SBJs)
                    scat_colors(plot_onset_sbj{cond_ix}(:,good_roi_map{cond_ix}(roi_ix))==sbj_ix,:) = repmat(SBJ_colors(sbj_ix,:),...
                        sum(plot_onset_sbj{cond_ix}(:,good_roi_map{cond_ix}(roi_ix))==sbj_ix),1);
                end
                violins(roi_ix).ScatterPlot.MarkerFaceColor = 'flat';   % Necessary for CData to work
                violins(roi_ix).ScatterPlot.MarkerEdgeColor = 'flat';   % Necessary for CData to work
                violins(roi_ix).ScatterPlot.CData = scat_colors;
            else   % violin.MedianColor is plotted over only ScatterPlot point
                violins(roi_ix).MedianColor = ...
                    SBJ_colors(plot_onset_sbj{cond_ix}(1,good_roi_map{cond_ix}(roi_ix)),:);
            end
        else
            % Change violin color to match ROI
            violins(roi_ix).ViolinColor = roi_colors{good_roi_map{cond_ix}(roi_ix)};
        end
    end
    
    % Add label and min RT for perspective
    ax = gca;
    ax.YLabel.String = y_label;
    if ~isempty(strfind(st.evnt_lab,'S'))
        line(xlim,[mean(min_rt) mean(min_rt)]);%,'k--');
    end
    
    if isfield(plt_vars,'mirror') && plt_vars.mirror
        view([90 -90]);
    end
    
    %% Save figure
    if save_fig
        fig_dir = [root_dir 'PRJ_Stroop/results/HFA/GRP_onsets_ROI/'...
            stat_id '/' roi_id '/' an_id '/' plt_id '/'];
        if ~exist(fig_dir,'dir')
            [~] = mkdir(fig_dir);
        end
        
        fig_filename = [fig_dir fig_name '.' fig_ftype];
        fprintf('Saving %s\n',fig_filename);
        saveas(gcf,fig_filename);
        %eval(['export_fig ' fig_filename]);
    end
end

end
