function SBJ10b_ANOVA_plot_SR_RTcorr_gROIcomb(SBJ,proc_id,stat_id_s,stat_id_r,an_id_s,an_id_r,...
                                        atlas_id,gm_thresh,z_thresh,plot_nsig,plt_id,save_fig,fig_vis,fig_ftype)
% Plots smANOVA w2 time series by ROI groupings
%   gROI version: one plot; subplots per gROI (lines still colored by ROI)
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
%   plot_nsig [0/1] - plot electrodes with no significant epochs?
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

% Load sample rate
load(strcat(SBJ_vars.dirs.events,SBJ,'_trial_info_final.mat'),'trial_info');
srate = trial_info.sample_rate;

% Load smANOVA results
w2 = cell([1 2]); evnt_labs = cell([1 2]);
f_name_s = [SBJ_vars.dirs.proc SBJ '_smANOVA_ROI_' stat_id_s '_' an_id_s '.mat'];
f_name_r = [SBJ_vars.dirs.proc SBJ '_smANOVA_ROI_' stat_id_r '_' an_id_r '.mat'];
tmp = load(f_name_s,'w2'); w2{1} = tmp.w2;
tmp = load(f_name_s,'st'); evnt_labs{1} = tmp.st.evnt_lab{1};   % handles {'S'} or {'S' 'S'}
tmp = load(f_name_r,'w2'); w2{2} = tmp.w2;
tmp = load(f_name_r,'st'); evnt_labs{2} = tmp.st.evnt_lab;      % assumes {'R'}
clear tmp
load(f_name_s,'st');

% Convert to Percentages
w2{1}.trial = w2{1}.trial.*100;
w2{2}.trial = w2{2}.trial.*100;

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
if ~plot_nsig
    plot_elecs = w2{1}.label(squeeze(any(any(sig_idx,1),2)));
    sig_suffix = '_sig';
else
    plot_elecs = w2{1}.label;
    sig_suffix = '';
end

%% Load ROI and GM/WM info
% Get list of ROIs and colors
[all_roi_list, all_roi_colors] = fn_roi_label_styles('all');

% Load elec file
elec_fname = [SBJ_vars.dirs.recon,SBJ,'_elec_',proc_id,'_pat_',atlas_id,'_final.mat'];
load(elec_fname);

% Remove hemi and/or atlas elecs
plot_elecs = intersect(plot_elecs,fn_select_elec_lab_match(elec, 'b', atlas_id, ''));

% Get GM probability from tissue labels {'GM','WM','CSF','OUT'}
if gm_thresh>0
    plot_elecs  = intersect(plot_elecs,elec.label(elec.tissue_prob(:,1)>gm_thresh));
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

% Add ROI plotting color
elec.color  = cell(size(elec.label));
for e_ix = 1:numel(elec.label)
    elec.color(e_ix) = all_roi_colors(strcmp(elec.ROI{e_ix},all_roi_list));
end

%% Plot Results
fig_dir = [root_dir 'PRJ_Stroop/results/HFA/' SBJ '/smANOVA_ts/' stat_id_s...
    '-' stat_id_r '/' an_id_s '-' an_id_r '/gROIcomb/'];
if ~exist(fig_dir,'dir')
    [~] = mkdir(fig_dir);
end

% Find plot limits
max_w2 = max([max(max(max(w2{1}.trial))) max(max(max(w2{2}.trial)))]);
min_w2 = min([min(min(min(w2{1}.trial))) min(min(min(w2{2}.trial)))]);
ylim_fudge = (max_w2-min_w2)*plt_vars.ylim_fudge;
ylims  = [min_w2-ylim_fudge max_w2+ylim_fudge];
yticks = 0:1:ylims(2);

% Create a figure for each condition and each gROI
for cond_ix = 1:numel(w2{1}.cond)
    fig_name = [SBJ '_smANOVA_SR_' w2{1}.cond{cond_ix} '_' atlas_id '_gROI' sig_suffix];
    f = figure('Name',fig_name,'units','normalized',...
        'outerposition',[0 0 1 1],'Visible',fig_vis); %twice as wide for the double plot
    any_plot = 0;
    groi_list = unique(elec.gROI);
    % Subplots for each ROI and event-locked epoch
    for groi_ix = 1:numel(groi_list)
        roi_list = unique(elec.ROI(strcmp(elec.gROI,groi_list{groi_ix})));
        for sr_ix = 1:2
            subplot_ix = fn_rowcol2subplot_ix(numel(groi_list),2,groi_ix,sr_ix);
            ax = subplot(numel(groi_list),2,subplot_ix);
            hold on;
            roi_lines     = gobjects(size(roi_list));
            roi_sig_count = zeros(size(roi_list));
            ch_list = find(strcmp(elec.gROI,groi_list{groi_ix}));
            for e_ix = 1:numel(ch_list)
                ch_ix = find(strcmp(w2{sr_ix}.label,elec.label{ch_list(e_ix)}));
                % Plot significant time periods
                any_plot = 1;
                roi_line = plot(w2{sr_ix}.time,squeeze(w2{sr_ix}.trial(cond_ix,ch_ix,:))',...
                    'Color',[elec.color{ch_list(e_ix)}]);%,'LineStyle',grp_style{cond_ix});
                
                % Track number of elecs per ROI and legend
                roi_ix = strcmp(roi_list,elec.ROI{ch_list(e_ix)});
                roi_lines(roi_ix) = roi_line;
                
                % Plot significant epochs in bold
                if sig_idx(sr_ix,cond_ix,ch_ix)
                    roi_sig_count(roi_ix) = roi_sig_count(roi_ix) + 1;
                    sig_chunks = fn_find_chunks(squeeze(w2{sr_ix}.qval(cond_ix,ch_ix,:))<=st.alpha);
                    sig_chunks(squeeze(w2{sr_ix}.qval(cond_ix,ch_ix,sig_chunks(:,1)))>st.alpha,:) = [];
                    if size(sig_chunks,1)==1 && sig_chunks(1,1)==sig_chunks(1,2)
                        % Handle case of single significant window
                        scatter(w2{sr_ix}.time(sig_chunks(1,1)),squeeze(w2{sr_ix}.trial(cond_ix,ch_ix,sig_chunks(1,1))),...
                            50,elec.color{ch_list(e_ix)},'o','filled');
                    else
                        for sig_ix = 1:size(sig_chunks,1)
                            line(w2{sr_ix}.time(sig_chunks(sig_ix,1):sig_chunks(sig_ix,2)),...
                                squeeze(w2{sr_ix}.trial(cond_ix,ch_ix,sig_chunks(sig_ix,1):sig_chunks(sig_ix,2))),...
                                'Color',[elec.color{ch_list(e_ix)}],...'LineStyle',plt_vars.sig_style,
                                'LineWidth',plt_vars.sig_width);
                        end
                    end
                end
            end
            
            % Plot event
            if strcmp(evnt_labs{sr_ix},'S')
                x_tick_lab = plt_vars.plt_lim_S(1):plt_vars.x_step_sz:plt_vars.plt_lim_S(2);
                legend_loc = 'northeast';
            else
                x_tick_lab = plt_vars.plt_lim_R(1):plt_vars.x_step_sz:plt_vars.plt_lim_R(2);
                legend_loc = 'northwest';
                % Plot Response Marker
                event_line = line([0 0],ylims, 'LineWidth',plt_vars.evnt_width,...
                    'Color',plt_vars.evnt_color, 'LineStyle',plt_vars.evnt_style);
            end
            
            % Plotting parameters
            ax.Title.String  = [groi_list{groi_ix} ' (n sig=' num2str(sum(roi_sig_count)) '; ' strjoin(roi_list,',') '): ' evnt_labs{sr_ix}];
            ax.Box           = 'off';
            ax.YLim          = ylims;
            ax.YTick         = yticks;
            ax.YLabel.String = 'Omega^2';
            ax.XLim          = eval(['plt_vars.plt_lim_' evnt_labs{sr_ix}]);%[min(w2{sr_ix}.win_lim_s(:,1)) max(w2{sr_ix}.win_lim_s(:,2))];
            %                 ax.XTick         = 0:plt_vars.x_step_sz*srate:size(w2{sr_ix}.time,2);
            ax.XTick         = x_tick_lab;
            ax.XLabel.String = 'Time (s)';
            
            % Legend
            roi_legend = cell(size(roi_list));
            for roi_ix = 1:numel(roi_list)
                roi_legend{roi_ix} = [roi_list{roi_ix} ' (n sig=' num2str(roi_sig_count(roi_ix)) ')'];
            end
            legend(roi_lines,roi_legend,'Location',legend_loc);
        end
        
    end
    
    % Save figure
    if any_plot
        if save_fig
            fig_filename = [fig_dir fig_name '.' fig_ftype];
            fprintf('Saving %s\n',fig_filename);
            saveas(gcf,fig_filename);
            %eval(['export_fig ' fig_filename]);
        end
    else
        close(f);
    end
end

end
