function SBJ10b_ANOVA_plot_SR(SBJ,proc_id,stat_id_s,stat_id_r,an_id_s,an_id_r,...
                              atlas_id,gm_thresh,z_thresh,roi_id,plt_id,save_fig,fig_vis,fig_ftype)
% Plots smANOVA w2 time series per electrode
%   ROI version: one plot per gROI, subplots for ROI
% INPUTS:
%   SBJ [str] - subject ID to plot
%   proc_id [str] - name of analysis pipeline, used to pick elec file
%   stat_id_s [str] - ID of S-locked statistical analysis
%   stat_id_r [str] - ID of R-locked statistical analysis
%   an_id_s [str] - analysis ID for S-locked preprocessing, filtering, etc.
%   an_id_r [str] - analysis ID for S-locked preprocessing, filtering, etc.
%   atlas_id [str] - {'DK','Dx','Yeo7','Yeo17'}
%   gm_thresh [float] - threshold of GM % to include electrode (likely = 0)
%   z_thresh [float] - threshold of HFA z score to include electrode
%   plt_id [str] - ID of the plotting variables
%   save_fig [0/1] - save this figure?
%   fig_vis [str] - visible = 'on'/'off'
%   fig_ftype [str] - file extension for figure saving
if ischar(save_fig); save_fig = str2num(save_fig); end

%% Data Preparation
% Set up paths
[root_dir, app_dir] = fn_get_root_dir(); ft_dir = [app_dir 'fieldtrip/'];
addpath([root_dir 'PRJ_Stroop/scripts/']);
addpath([root_dir 'PRJ_Stroop/scripts/utils/']);
addpath(ft_dir);
ft_defaults

%% Load Results
SBJ_vars_cmd = ['run ' root_dir 'PRJ_Stroop/scripts/SBJ_vars/' SBJ '_vars.m'];
eval(SBJ_vars_cmd);
plt_vars_cmd = ['run ' root_dir 'PRJ_Stroop/scripts/plt_vars/' plt_id '_vars.m'];
eval(plt_vars_cmd);

% Load RTs
load(strcat(SBJ_vars.dirs.events,SBJ,'_trial_info_final.mat'),'trial_info');
srate = trial_info.sample_rate;

% Load smANOVA results
w2 = cell([1 2]); evnt_labs = cell([1 2]);
f_name_s = [SBJ_vars.dirs.proc SBJ '_smANOVA_ROI_' stat_id_s '_' an_id_s '.mat'];
f_name_r = [SBJ_vars.dirs.proc SBJ '_smANOVA_ROI_' stat_id_r '_' an_id_r '.mat'];
tmp = load(f_name_s,'w2'); w2{1} = tmp.w2;
tmp = load(f_name_s,'st'); evnt_labs{1} = tmp.st.evnt_lab;   % handles {'S'} or {'S' 'S'}
tmp = load(f_name_r,'w2'); w2{2} = tmp.w2;
tmp = load(f_name_r,'st'); evnt_labs{2} = tmp.st.evnt_lab;      % assumes {'R'}
clear tmp
load(f_name_s,'st');

% Convert to Percentages
w2{1}.trial = w2{1}.trial.*100;
w2{2}.trial = w2{2}.trial.*100;

evnt_lab = {'S', 'R'};
[cond_lab, cond_colors, cond_style] = fn_condition_label_styles(st.model_lab);
if ~all(strcmp(cond_lab,w2{1}.cond)); error('cond-lab mismatch'); end

% Find significant channels
sig_idx = false([2 numel(w2{1}.cond) numel(w2{1}.label)]);
for sr_ix = 1:2
    for cond_ix = 1:numel(w2{sr_ix}.cond)
        for ch_ix = 1:numel(w2{sr_ix}.label)
            if any(squeeze(w2{sr_ix}.qval(cond_ix,ch_ix,:))<=st.alpha)
                sig_idx(sr_ix,cond_ix,ch_ix) = true;
            end
        end
    end
end

%% Load ROI and GM/WM info
% Load elec file
elec_fname = [SBJ_vars.dirs.recon,SBJ,'_elec_',proc_id,'_pat_',atlas_id,'_final.mat'];
load(elec_fname);
if numel(elec.label)~=numel(w2{1}.label) || ~all(strcmp(elec.label,w2{1}.label))
    error('ch_lab mismatch in elec and w2');
end

% Remove hemi and/or atlas elecs
plot_elecs = fn_select_elec_lab_match(elec, 'b', atlas_id, '');

% Get GM probability from tissue labels {'GM','WM','CSF','OUT'}
if gm_thresh>0
    plot_elecs = intersect(plot_elecs,elec.label(elec.tissue_prob(:,1)>gm_thresh));
end

% Get z_threshold
if z_thresh>0
    sdr_w2 = load([SBJ_vars.dirs.proc SBJ ...
        '_smANOVA_ROI_CNI_PC_S0tmRT_WL1_WS25_D1tRT_R1t5_WL1_WS25' ...
        '_HGm_S2t151_zbtA_sm0_wn100.mat'],'w2');
    plot_elecs = intersect(plot_elecs,elec.label(sdr_w2.w2.max_hfa_z>=z_thresh));
end

% Select elecs (significance, ROI, tissue, z-score)
cfgs = []; cfgs.channel = plot_elecs;
elec = fn_select_elec(cfgs,elec);

%% Plot Results
fig_dir = [root_dir 'PRJ_Stroop/results/HFA/' SBJ '/smANOVA_ts/' stat_id_s '-' stat_id_r '/' ...
    an_id_s '-' an_id_r '/' roi_id '/'];
sig_ln_dir = [fig_dir 'sig_ch/'];
if ~exist(sig_ln_dir,'dir')
    [~] = mkdir(sig_ln_dir);
end

% Find plot limits
max_w2 = max([max(max(max(w2{1}.trial))) max(max(max(w2{2}.trial)))]);
min_w2 = min([min(min(min(w2{1}.trial))) min(min(min(w2{2}.trial)))]);
ylim_fudge = (max_w2-min_w2)*plt_vars.ylim_fudge;
ylims  = [min_w2-ylim_fudge max_w2+ylim_fudge];
yticks = 0:1:ylims(2);

% Create a figure combining statistics for each electrode
for ch_ix = 1:numel(w2{1}.label)
    fig_name = [SBJ '_smANOVA_SR_' w2{1}.cond{cond_ix} '_' w2{1}.label{ch_ix}];
    f = figure('Name',fig_name,'units','normalized',...
        'outerposition',[0 0 1 0.5],'Visible',fig_vis); %twice as wide for the double plot
    ax = gobjects(size(evnt_lab));
    for sr_ix = 1:2
        ax(sr_ix) = subplot(1,2,sr_ix);
        hold on;
        for cond_ix = 1:numel(w2{1}.cond)
            % Plot effect size
            plot(w2{sr_ix}.time,squeeze(w2{sr_ix}.trial(cond_ix,ch_ix,:))',...
                'Color',cond_colors{cond_ix});%,'LineStyle',grp_style{cond_ix});
            
            if sig_idx(sr_ix,cond_ix,ch_ix)
                % Find significant periods
                sig_chunks = fn_find_chunks(squeeze(w2{sr_ix}.qval(cond_ix,ch_ix,:))<=st.alpha);
                sig_chunks(squeeze(w2{sr_ix}.qval(cond_ix,ch_ix,sig_chunks(:,1)))>st.alpha,:) = [];
                if size(sig_chunks,1)==1 && sig_chunks(1,1)==sig_chunks(1,2)
                    % Handle case of single significant window
                    scatter(w2{sr_ix}.time(sig_chunks(1,1)),squeeze(w2{sr_ix}.trial(cond_ix,ch_ix,sig_chunks(1,1))),...
                        50,cond_colors{cond_ix},'o','filled');
                else
                    for sig_ix = 1:size(sig_chunks,1)
                        line(w2{sr_ix}.time(sig_chunks(sig_ix,1):sig_chunks(sig_ix,2)),...
                            squeeze(w2{sr_ix}.trial(cond_ix,ch_ix,sig_chunks(sig_ix,1):sig_chunks(sig_ix,2))),...
                            'Color',cond_colors{cond_ix},...'LineStyle',plt_vars.sig_style,
                            'LineWidth',plt_vars.sig_width);
                    end
                end
            end
        end
        
        % Plot event
        if strcmp(evnt_lab{sr_ix},'S')
            x_tick_lab = plt_vars.plt_lim_S(1):plt_vars.x_step_sz:plt_vars.plt_lim_S(2);
        else
            x_tick_lab = plt_vars.plt_lim_R(1):plt_vars.x_step_sz:plt_vars.plt_lim_R(2);
            % Plot Response Marker
            event_line = line([0 0],ylims,'LineWidth',plt_vars.evnt_width,...
                'Color',plt_vars.evnt_color, 'LineStyle',plt_vars.evnt_style);
        end
        
        % Plotting parameters
        ax(sr_ix).Title.String  = [w2{sr_ix}.label{ch_ix} ' (' elec.ROI{ch_ix} '): ' evnt_lab{sr_ix}];
        ax(sr_ix).Box           = 'off';
        ax(sr_ix).YLim          = ylims;
        ax(sr_ix).YTick         = yticks;
        ax(sr_ix).YLabel.String = '% Variance Explained';
        ax(sr_ix).XLim          = [x_tick_lab(1) x_tick_lab(end)];%[0,size(stat{sr_ix}.time,2)];
        ax(sr_ix).XTick         = x_tick_lab;%0:plt_vars.x_step_sz*srate:size(stat{sr_ix}.time,2);
        %     ax(sr_ix).XTickLabel    = x_tick_lab;
        ax(sr_ix).XLabel.String = 'Time (s)';
        %             ax(2).YColor        = 'k';
        %             legend(main_lines,lgd_lab{:},'Location',plt_vars.legend_loc);
    end
    
    % Save figure
    if save_fig
        fig_fname = [fig_dir fig_name '.' fig_ftype];
        fprintf('Saving %s\n',fig_fname);
        saveas(f,fig_fname);
    end
end


%% Save out list of channels with significant differences
% sig_report_filename = [fig_dir 'sig_ch_list.txt'];
% sig_report = fopen(sig_report_filename,'a');
% fprintf(sig_report,'%s: %s for %s - %s\n',SBJ,stat_id,an_id_s,an_id_r);
% for grp_ix = 1:numel(grp_lab)+1
%     for sr_ix = 1:numel(event_lab)
%         fprintf(sig_report,'%s, %s :\n',lgd_lab{grp_ix},event_lab{sr_ix});
%         if ~isempty(sig_ch{grp_ix,sr_ix})
%             fprintf(sig_report,'\t%s\n',strjoin(sig_ch{grp_ix,sr_ix},','));
%         else
%             fprintf(sig_report,'\n');
%         end
%     end
% end
% fclose(sig_report);

end


%% ========================================================
%                 % Plot var_exp
%                 [axs,grp_lines,rt_line] = plotyy(...
%                     repmat(win_center{sr_ix},size(grp_lab)),squeeze(w2{sr_ix}.trial(:,ch_ix,:))',...
%                     1:numel(stat{sr_ix}.time),squeeze(stat{sr_ix}.rho(ch_ix,:,:))');
%                 for grp_ix = 1:numel(grp_lab)
%                     set(grp_lines(grp_ix),'Color',[grp_colors{grp_ix}],'LineStyle',grp_style{grp_ix});
%                 end
%                 set(rt_line,'Color',rt_color{1},'LineStyle',rt_style{1});
%                 main_lines = [main_lines grp_lines' rt_line];
%                 lgd_lab = {grp_lab{:} 'corr(RT)'};

%                     % RT correlation significant time periods
%                     ylims = ylims2;
%                     yticks = yticks2;
%                     ylab = 'Correlation with RT';
%                     if sum(stat{sr_ix}.mask(ch_ix,:))>0 && strcmp(roi_lab{ch_ix},roi_list{roi_ix})
%                         any_plot = 1;
%                         %                         % Track channels with significant var_exp
%                         %                         if ~iscell(sig_ch{cond_ix,sr_ix})
%                         %                             sig_ch{cond_ix,sr_ix} = {};
%                         %                         end
%                         %                         sig_ch{cond_ix,sr_ix} = {sig_ch{cond_ix,sr_ix}{:} stat{sr_ix}.label{ch_ix}};
%                         plot(1:numel(stat{sr_ix}.time),squeeze(stat{sr_ix}.rho(ch_ix,:,:))',...
%                             'Color',[roi_colors{roi_ix}],'LineStyle','-');
%                         
%                         sig_chunks = fn_find_chunks(squeeze(stat{sr_ix}.mask(ch_ix,:,:)));
%                         sig_chunks(squeeze(stat{sr_ix}.mask(ch_ix,:,sig_chunks(:,1)))==0,:) = [];
%                         sig_times = find(squeeze(stat{sr_ix}.mask(ch_ix,:,:)==1));
%                         if strcmp(plt_vars.sig_type,'bold')
%                             for sig_ix = 1:size(sig_chunks,1)
%                                 sig_times = sig_chunks(sig_ix,1):sig_chunks(sig_ix,2);
%                                 line(sig_times,squeeze(stat{sr_ix}.rho(ch_ix,:,sig_times)),...
%                                     'Color',roi_colors{roi_ix},'LineStyle',plt_vars.sig_style,...
%                                     'LineWidth',plt_vars.sig_width);
%                             end
%                         elseif strcmp(plt_vars.sig_type,'scatter')
%                             scatter(sig_times,repmat(y_sig(cond_ix),size(sig_times)),...
%                                 plt_vars.sig_scat_size2,roi_colors{roi_ix},plt_vars.sig_scat_mrkr2);
%                         end
%                         fprintf('%s RT (%s) -- SIGNIFICANT CLUSTERS FOUND, plotting with significance shading...\n',...
%                             stat{sr_ix}.label{ch_ix},event_lab{sr_ix});%,size(sig_chunks,1));
% %                     else
% %                         fprintf('%s RT (%s) -- NO SIGNIFICANT CLUSTERS FOUND, plotting without significance shading...\n',...
% %                             stat{sr_ix}.label{ch_ix},event_lab{sr_ix});
%                     end
