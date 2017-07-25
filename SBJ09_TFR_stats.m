function SBJ09_TFR_stats(SBJ,conditions,pipeline_id,an_id,save_data,n_boots,save_stats,plot_vis,save_plot)
% Calculates time frequency representation, computes cluster-based statistics, and plots the results
% INPUTS:
%   SBJ [str] - dataset to be processed
%   conditions [str] - experimental conditions to be compared
%   pipeline_id [str] - which processed data pipeline to get the data
%   an_id [str] - which analysis variabels to use
%   save_data [0/1] - flag to save TFR or not
%   n_boots [int] - iterations of stats permutations
%   save_stats [0/1] - flag to save stats or not
%   plot_vis ['on'/'off'] - make the plot visible or not
%   save_plot [0/1] - save the plots or not

% clear all; %close all;
% Set up paths
addpath('/home/knight/hoycw/PRJ_Stroop/scripts/');
addpath('/home/knight/hoycw/PRJ_Stroop/scripts/utils/');
addpath('/home/knight/hoycw/Apps/fieldtrip/');
ft_defaults

fig_filetype = 'png';

%% Data Preparation
% Load Data
eval(['run /home/knight/hoycw/PRJ_Stroop/scripts/SBJ_vars/' SBJ '_vars.m']);
eval(['run /home/knight/hoycw/PRJ_Stroop/scripts/proc_vars/' pipeline_id '_proc_vars.m']);
eval(['run /home/knight/hoycw/PRJ_Stroop/scripts/an_vars/' an_id '_vars.m']);

load(strcat(SBJ_vars.dirs.preproc,SBJ,'_preproc_',pipeline_id,'.mat'));
load(strcat(SBJ_vars.dirs.events,SBJ,'_trial_info_final.mat'));

% Select Conditions of Interest
[cond_lab, cond_colors, cond_style] = fn_condition_label_styles(conditions);
cond_idx = false([length(cond_lab) length(trial_info.trial_n)]);
for cond_ix = 1:length(cond_lab)
    % Get binary condition index
    cond_idx(cond_ix,:) = logical(fn_condition_index(cond_lab{cond_ix},...
        trial_info.condition_n));
end

%% Select Channel(s)
cfgs = [];
cfgs.channel = SBJ_vars.ch_lab.ROI;
roi = ft_selectdata(cfgs,data);

%% Cut into Trials
if strcmp(event_type,'stim')
    events = trial_info.word_onset;
elseif strcmp(event_type,'resp')
    events = trial_info.resp_onset;
else
    error(stract('ERROR: unknown event_type ',event_type));
end
roi_trl = fn_ft_cut_trials_equal_len(roi,events,trial_info.condition_n',trial_lim_s*roi.fsample);

%% Compute TFRs
% all cfg_tfr options are specified in the an_vars
tfr = {};
n_trials = zeros([1 numel(cond_lab)]);
for cond_ix = numel(cond_lab)
    fprintf('===================================================\n');
    fprintf('------------- TFR Calculations for %s ----------\n',cond_lab{cond_ix});
    fprintf('===================================================\n');
    cfg_tfr.trials = find(cond_idx(cond_ix,:)==1);
    tfr{cond_ix}   = ft_freqanalysis(cfg_tfr, roi_trl);
    
    % Grab n_trials for design matrix
    n_trials(cond_ix) = size(tfr{cond_ix}.trialinfo,1);
end

% %% Baseline Correction --> not doing this since I'm comparing conditions directly, not to baseline
% cfg = [];
% cfg.demean         = demean_yn;
% cfg.baselinewindow = bsln_lim;
% roi_trl = ft_preprocessing(cfg,roi_trl);

%% Run Statistics
% Create design matrix
design = zeros(2,sum(n_trials));
for an_ix = 1:numel(cond_lab)
    if an_ix==1
        design(1,1:n_trials(an_ix)) = an_ix;                                % Conditions (Independent Variable)
        design(2,1:n_trials(an_ix)) = 1:n_trials(an_ix);                    % Trial Numbers
    else
        design(1,sum(n_trials(1:an_ix-1))+1:sum(n_trials(1:an_ix)))= an_ix; % Conditions (Independent Variable)
        design(2,sum(n_trials(1:an_ix-1))+1:sum(n_trials(1:an_ix)))= 1:n_trials(an_ix);
    end
end

% Calculate statistics
cfg_stat = [];
cfg_stat.latency          = stat_lim;
cfg_stat.channel          = {'RSM*','LSM*'};%'all';
cfg_stat.parameter        = 'powspctrm';    %for now, assume power
cfg_stat.method           = 'montecarlo';
cfg_stat.statistic        = 'ft_statfun_indepsamplesT';
cfg_stat.correctm         = 'cluster';
cfg_stat.clusteralpha     = 0.05;   %threshold for a single comparison (time point) to be included in the clust
cfg_stat.clusterstatistic = 'maxsum';
cfg_stat.clustertail      = 0;
cfg_stat.tail             = 0; %two sided
cfg_stat.correcttail      = 'alpha'; %correct the .alpha for two-tailed test (/2)
cfg_stat.alpha            = 0.05;
cfg_stat.numrandomization = n_boots;
cfg_stat.design           = design;
cfg_stat.neighbours       = [];%neighbors;
% cfg_stat.minnbchan        = 0;
cfg_stat.ivar             = 1;  %row of design matrix containing independent variable
% cfg_stat.uvar             = 2;  %row of design matrix containing dependent variable, not needed for indepsamp
[stat] = ft_freqstatistics(cfg_stat, tfr{:});

% Save Results
if save_data
    data_out_filename = strcat(SBJ_vars.dirs.proc,SBJ,'_TFR_stats_',conditions,...
                                '_ROI_',pipeline_id,'_',an_id,'.mat');
    fprintf('Saving %s\n',data_out_filename);
    save(data_out_filename, '-v7.3','stat');
end
%% Plot Results
% stat_full = stat;
% roi_hfa_full = roi_hfa;
% Obtain the descriptive stats
tfr_avg = {};
cfg_avg = [];
cfg_avg.channel = {'RSM*','LSM*'};
for cond_ix = 1:numel(cond_lab)
    tfr_avg{cond_ix} = ft_freqdescriptive(cfg_avg,tfr{cond_ix});
end

% Take the difference
tfr_diff = {};
for cond_ix = 2:numel(cond_lab)
    tfr_diff{cond_ix-1} = tfr_avg{cond_ix}.powspctrm - tfr_avg{cond_ix-1}.powspctrm;
end
stat.raweffect = tfr_diff{1};

%% Cluster plots
cfg_clust = [];
cfg_clust.alpha  = 0.025;
cfg_clust.parameter = 'raweffect';
cfg_clust.zlim   = [-1e-27 1e-27];
cfg_clust.layout = layout;
ft_clusterplot(cfg_clust, stat);

%% Single Plots
cfg_s = [];
cfg_s.latency = [-0.5 2];
cfg_s.channel = 'RAC4-5';
cfg_lay = [];
cfg_lay.layout = 'vertical';
for cond_ix = 1:length(cond_lab)
    tfr_plot = ft_selectdata(cfg_s, tfr{cond_ix});
    layout = ft_prepare_layout(cfg_lay,tfr_plot);
%     subplot(length(tfr_plot.label),1,cond_ix);
    out_file = [SBJ '_TFR_' cond_lab{cond_ix} '_' event_type '_ROI_' pipeline_id '_' an_id];
%         '_Bob.ep' epoch_id '_bsln' bsln_id env_id smooth_id sig_id y_scale_id];
%     if n_chan > 2
    fig_height = 1;
    vis_fig='on';
%     else
%         fig_height = n_chan/3;
%     end
    figure('Name',out_file);%,'units','normalized',...
%         'outerposition',[0 0 1 fig_height],'Visible',vis_fig);
    cfg = [];
    cfg.trials = 'all';%logical(fn_condition_index(cond_lab{cond_ix}, cond_n));
    cfg.latency = [-0.5 2];
    cfg.baseline = bsln_lim; %should be in sec
    cfg.baselinetype = 'db';
    cfg.channel = 'all';
    cfg.colorbar = 'yes';
    cfg.zlim = 'maxabs'; %[-3 3];%
    cfg.parameter = 'powspctrm';
%     cfg.layout = layout;
    % cfg.maskparameter = 'trial';
    % cfg.maskstyle = 'outline';
    cfg.showlabels   = 'yes';
    ft_singleplotTFR(cfg, tfr_plot);
end

%% cfg_s = [];
% cfg_s.latency = [-0.5 2];
% tfr_plot = ft_selectdata(cfg_s, tfr);
% figure;
% for ch_ix = 1:length(tfr.label)
%     subplot(2,1,ch_ix);
%     cfg = [];
%     cfg.baseline = [-.2 0];%bsln_lim; %should be in sec
%     cfg.baselinetype = 'db';
%     cfg.channel = tfr.label{ch_ix};
%     cfg.colorbar = 'yes';
%     cfg.zlim = 'maxabs'; %[-5 5];
%     cfg.parameter = 'powspctrm';
%     % cfg.maskparameter = 'trial';
%     % cfg.maskstyle = 'outline';
%     cfg.showlabels   = 'yes';
%     ft_singleplotTFR(cfg, tfr_plot);
% end


% Create a figure for each ROI
for roi_ix = 1:numel(SBJ_vars.ch_lab.ROI)
    % Select data to plot this ROI
    cfgs = [];
    cfgs.channel = SBJ_vars.ch_lab.ROI{roi_ix};
    stat = ft_selectdata(cfgs,stat_full);
    for an_ix = 1:numel(cond_lab)
        roi_hfa{an_ix} = ft_selectdata(cfgs,roi_hfa_full{an_ix});
    end
    
    % Plot parameters
    roi_name = SBJ_vars.ch_lab.ROI{roi_ix}; if strcmp(roi_name(end),'*');roi_name=roi_name(1:end-1);end
    fig_name = [SBJ '_HFA_stat_' conditions '_' roi_name '_' event_type];
    [plot_rc,~] = fn_num_subplots(numel(stat.label));
    if plot_rc(1)>1; fig_height=1; else fig_height=0.33; end;
    
    figure('Name',fig_name,'units','normalized',...
        'outerposition',[0 0 1 fig_height],'Visible',vis_fig);
    plot_info.fig        = gcf;
    plot_info.x_step     = 0.25*roi.fsample;
    plot_info.x_lab      = trial_lim_s(1):0.25:trial_lim_s(2);
    plot_info.legend_loc = 'southeast';
    plot_info.sig_alpha  = 0.2;
    plot_info.sig_color  = [0.5 0.5 0.5];
    % Stimulus plotting params
    event_info.time      = -trial_lim_s(1)*roi.fsample;
    event_info.name      = 'stim';
    event_info.width     = 2;
    event_info.color     = 'k';
    event_info.style     = '--';
    % Condition plotting params
    cond_info.name       = cond_lab;
    cond_info.style      = cond_style;
    cond_info.color      = cond_colors;
    cond_info.alpha      = [0.5 0.5];
    
    % Plot each channel within this ROI
    for ch_ix = 1:numel(stat.label)
        subplot(plot_rc(1),plot_rc(2),ch_ix);
        plot_info.ax         = gca;
        plot_info.title      = stat.label{ch_ix};
        if ch_ix==1; plot_info.legend=1; else plot_info.legend=0; end;
        
        % Compute means and variance
        means = NaN([numel(cond_lab) size(roi_hfa{1}.avg,2)]);
        var = NaN([numel(cond_lab) size(roi_hfa{1}.avg,2)]);
        for an_ix = 1:numel(cond_lab)
            means(an_ix,:) = roi_hfa{an_ix}.avg(ch_ix,:);
            var(an_ix,:) = squeeze(std(roi_hfa{an_ix}.trial(:,ch_ix,:),[],1)./sqrt(size(roi_hfa{an_ix}.trial,1)))';
        end
        % Find significant time periods
        if sum(stat.mask(ch_ix,:))>0
            mask_chunks = fn_find_chunks(stat.mask(ch_ix,:));
            sig_chunks = mask_chunks;
            sig_chunks(stat.mask(ch_ix,sig_chunks(:,1))==0,:) = [];
            % If stat and roi_erp aren't on same time axis, adjust sig_chunk indices
            if (size(stat.time,2)~=size(roi_hfa{1}.time,2)) || (sum(stat.time==roi_hfa{1}.time)~=numel(stat.time))
                for chunk_ix = 1:size(sig_chunks,1)
                    sig_chunks(chunk_ix,1) = find(roi_hfa{1}.time==stat.time(sig_chunks(chunk_ix,1)));
                    sig_chunks(chunk_ix,2) = find(roi_hfa{1}.time==stat.time(sig_chunks(chunk_ix,2)));
                end
            end
            fprintf('%i SIGNIFICANT CLUSTERS FOUND, plotting with significance shading...\n',size(sig_chunks,1));
            fn_plot_ts_error_bar_sig(plot_info,means,var,sig_chunks,event_info,cond_info);
        else
            fprintf('NO SIGNIFICANT CLUSTERS FOUND, plotting without significance shading...\n');
            fn_plot_ts_error_bar(plot_info,means,var,event_info,cond_info);
        end
    end
    clear stat roi_erp
    
    % Save figure
    if save_fig
        fig_dir = ['/home/knight/hoycw/PRJ_Stroop/results/HFA/' SBJ '/' conditions '/'];
        if ~exist(fig_dir,'dir')
            mkdir(fig_dir);
        end
        fig_filename = [fig_dir fig_name '.' fig_filetype];
        fprintf('Saving %s\n',fig_filename);
        saveas(gcf,fig_filename);
        %eval(['export_fig ' fig_filename]);
    end
end
%%
% % Plot ERPs
% roi_erp_con.mask = stat.mask;
% roi_erp_inc.mask = stat.mask;
% cfgp = [];
% cfgp.showlabels = 'yes';
% cfgp.parameter = 'avg';
% % cfgp.layout = 'ordered';
% cfgp.maskparameter = 'mask';
% ft_singleplotER(cfgp, roi_erp_con, roi_erp_inc);%stat); %roi_erp_con, roi_erp_inc, 

end