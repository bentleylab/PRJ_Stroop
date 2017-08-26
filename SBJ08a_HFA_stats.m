function SBJ08a_HFA_stats(SBJ,conditions,pipeline_id,an_id)
% Calculates high frequency activity, computes cluster-based statistics, and plots the results
% clear all; %close all;

%% Data Preparation
% Set up paths
addpath('/home/knight/hoycw/PRJ_Stroop/scripts/');
addpath('/home/knight/hoycw/PRJ_Stroop/scripts/utils/');
addpath('/home/knight/hoycw/Apps/fieldtrip/');
ft_defaults

SBJ_vars_cmd = ['run /home/knight/hoycw/PRJ_Stroop/scripts/SBJ_vars/' SBJ '_vars.m'];
eval(SBJ_vars_cmd);
an_vars_cmd = ['run /home/knight/hoycw/PRJ_Stroop/scripts/an_vars/' an_id '_vars.m'];
eval(an_vars_cmd);
% eval(['run /home/knight/hoycw/PRJ_Stroop/scripts/proc_vars/' pipeline_id '_proc_vars.m']);

% Load Data
load(strcat(SBJ_vars.dirs.preproc,SBJ,'_preproc_',pipeline_id,'.mat'));
load(strcat(SBJ_vars.dirs.events,SBJ,'_trial_info_final.mat'));

% Select Conditions of Interest
[cond_lab, ~, ~] = fn_condition_label_styles(conditions);
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
% Data
if strcmp(event_type,'stim')
    events = trial_info.word_onset;
elseif strcmp(event_type,'resp')
    events = trial_info.resp_onset;
else
    error(strcat('ERROR: unknown event_type ',event_type));
end
roi_trl = fn_ft_cut_trials_equal_len(roi,events,trial_info.condition_n',trial_lim_s*roi.fsample);

% Extract Baseline if it is timelocked to a different event
if ~strcmp(event_type,bsln_evnt)
    if strcmp(bsln_evnt,'stim')
        bsln_evnts = trial_info.word_onset;
    elseif strcmp(bsln_evnt,'resp')
        bsln_evnts = trial_info.resp_onset;
    else
        error(strcat('ERROR: unknown bsln_evnt ',bsln_evnt));
    end
    % Account for the loss of 1/2 TFR window
    bsln_lim = [bsln_lim(1)-cfg_hfa.t_ftimwin(1)/2 bsln_lim(2)+cfg_hfa.t_ftimwin(1)/2];
    bsln_trl = fn_ft_cut_trials_equal_len(roi,bsln_evnts,trial_info.condition_n',bsln_lim*roi.fsample);
end

%% Compute HFA
roi_hfa = {};
bsln_hfa = {};
n_trials = zeros([1 numel(cond_lab)]);
for cond_ix = 1:numel(cond_lab)
    fprintf('===================================================\n');
    fprintf('------------- TFR Calculations for %s ----------\n',cond_lab{cond_ix});
    fprintf('===================================================\n');
    if strcmp(HFA_type,'multiband')
        cfg_hfa.trials = find(cond_idx(cond_ix,:)==1);
        roi_hfa{cond_ix}   = ft_freqanalysis(cfg_hfa, roi_trl);
        if ~strcmp(event_type,bsln_evnt)
            bsln_hfa{cond_ix}   = ft_freqanalysis(cfg_hfa, bsln_trl);
        end
        
        % Grab n_trials for design matrix
        n_trials(cond_ix) = size(roi_hfa{cond_ix}.trialinfo,1);
    elseif strcmp(HFA_type,'broadband')
        error('Stop using broadband and use multitapers you dummy!');
        %         % Filter to HFA band
        %         cfgpp = [];
        %         cfgpp.hpfilter  = 'yes';
        %         cfgpp.hpfreq    = 70;
        %         cfgpp.lpfilter  = 'yes';
        %         cfgpp.lpfreq    = 150;
        %         roi = ft_preprocessing(cfgpp,roi);
        %         % Hilbert method to extract power
    else
        error('Unknown HFA_type provided');
    end
end

%% Baseline Correction
fprintf('===================================================\n');
fprintf('------------- Baseline Correction for %s ----------\n',cond_lab{cond_ix});
fprintf('===================================================\n');
for cond_ix = 1:numel(cond_lab)
    switch bsln_type
        case {'zscore', 'demean', 'my_relchange'}
            roi_hfa{cond_ix} = fn_bsln_ft_tfr(roi_hfa{cond_ix},bsln_hfa{cond_ix},bsln_type);
        case 'relchange'
            cfgbsln = [];
            cfgbsln.baseline     = bsln_lim;
            cfgbsln.baselinetype = bsln_type;
            cfgbsln.parameter    = 'powspctrm';
            roi_hfa{cond_ix} = ft_freqbaseline(cfgbsln,roi_hfa{cond_ix});
        otherwise
            error(['No baseline implemented for bsln_type: ' bsln_type]);
    end
end

% %% Smooth Power Time Series
% if smooth_pow_ts
%     error('Nope, you dont really want to smooth this before doing stats!');
%     fprintf('===================================================\n');
%     fprintf('------------- Smoothing Power for %s ----------\n',cond_lab{cond_ix});
%     fprintf('===================================================\n');
%     for cond_ix = 1:numel(cond_lab)
%         for ch_ix = 1:size(roi_hfa{1}.label)
%             for f_ix = 1:size(roi_hfa{1}.freq)
%                 if strcmp(lp_yn,'yes') && strcmp(hp_yn,'yes')
%                     roi_hfa{cond_ix}.powspctrm(:,ch_ix,f_ix,:) = fn_EEGlab_bandpass(...
%                         roi_hfa{cond_ix}.powspctrm(:,ch_ix,f_ix,:), roi.fsample, hp_freq, lp_freq);
%                 elseif strcmp(lp_yn,'yes')
%                     roi_hfa{cond_ix}.powspctrm(:,ch_ix,f_ix,:) = fn_EEGlab_lowpass(...
%                         roi_hfa{cond_ix}.powspctrm(:,ch_ix,f_ix,:), roi.fsample, lp_freq);
%                 elseif strcmp(hp_yn,'yes')
%                     error('Why are you only high passing?');
%                 else
%                     error('Why did you say yes smooth but no to both low and high pass?');
%                 end
%             end
%         end
%     end
% end

%% Merge multiple bands
if strcmp(HFA_type,'multiband')
    for cond_ix = 1:numel(cond_lab)
        roi_hfa{cond_ix}.freq = mean(roi_hfa{cond_ix}.freq);
        roi_hfa{cond_ix}.powspctrm = mean(roi_hfa{cond_ix}.powspctrm,3);
    end
end

%% Run Statistics
fprintf('===================================================\n');
fprintf('--------------------- Statistics ------------------\n');
fprintf('===================================================\n');
% Create design matrix
% design = zeros(2,size(roi_erp_con.trial,1) + size(roi_erp_inc.trial,1));
% % Conditions (Independent Variable)
% design(1,1:size(roi_erp_con.trial,1)) = 1;
% design(1,(size(roi_erp_con.trial,1)+1):(size(roi_erp_con.trial,1) + size(roi_erp_inc.trial,1)))= 2;
% % Trial Numbers
% design(2,1:size(roi_erp_con.trial,1)) = 1;
% design(2,(size(roi_erp_con.trial,1)+1):(size(roi_erp_con.trial,1) + size(roi_erp_inc.trial,1)))= 2;
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

% Prepare neighbors layout
% cfgn = [];
% cfgn.method  = 'distance';
% cfgn.layout  = 'ordered';
% cfgn.channel = elecs;
% neighbors    = ft_prepare_neighbours(cfgn,roi_erp_allch{1});
% for ch_ix = 1:numel(roi_erp{1}.label)
%     neighbors(ch_ix).label = roi_erp{1}.label{ch_ix};
%     neighbors(ch_ix).neighblabel = {};
% end

% Calculate statistics
cfg_stat.design           = design;
[stat] = ft_freqstatistics(cfg_stat, roi_hfa{:});

%% Save Results
data_out_filename = strcat(SBJ_vars.dirs.proc,SBJ,'_',conditions,'_ROI_',an_id,'.mat');
fprintf('===================================================\n');
fprintf('--- Saving %s ------------------\n',data_out_filename);
fprintf('===================================================\n');
save(data_out_filename,'roi_hfa','stat');

% %% Plot Results
% stat_full = stat;
% roi_hfa_full = roi_hfa;
% 
% % Create a figure for each ROI
% for roi_ix = 1:numel(proc_vars.ch_lab.ROI)
%     % Select data to plot this ROI
%     cfgs = [];
%     cfgs.channel = proc_vars.ch_lab.ROI{roi_ix};
%     stat = ft_selectdata(cfgs,stat_full);
%     for an_ix = 1:numel(cond_lab)
%         roi_hfa{an_ix} = ft_selectdata(cfgs,roi_hfa_full{an_ix});
%     end
%     
%     % Plot parameters
%     roi_name = proc_vars.ch_lab.ROI{roi_ix}; if strcmp(roi_name(end),'*');roi_name=roi_name(1:end-1);end
%     fig_name = [SBJ '_HFA_stat_' conditions '_' roi_name '_' event_type];
%     [plot_rc,~] = fn_num_subplots(numel(stat.label));
%     if plot_rc(1)>1; fig_height=1; else fig_height=0.33; end;
%     
%     figure('Name',fig_name,'units','normalized',...
%         'outerposition',[0 0 1 fig_height],'Visible',vis_fig);
%     plot_info.fig        = gcf;
%     plot_info.x_step     = 0.25*roi.fsample;
%     plot_info.x_lab      = trial_lim_s(1):0.25:trial_lim_s(2);
%     plot_info.legend_loc = 'southeast';
%     plot_info.sig_alpha  = 0.2;
%     plot_info.sig_color  = [0.5 0.5 0.5];
%     % Stimulus plotting params
%     event_info.time      = -trial_lim_s(1)*roi.fsample;
%     event_info.name      = 'stim';
%     event_info.width     = 2;
%     event_info.color     = 'k';
%     event_info.style     = '--';
%     % Condition plotting params
%     cond_info.name       = cond_lab;
%     cond_info.style      = cond_style;
%     cond_info.color      = cond_colors;
%     cond_info.alpha      = [0.5 0.5];
%     
%     % Plot each channel within this ROI
%     for ch_ix = 1:numel(stat.label)
%         subplot(plot_rc(1),plot_rc(2),ch_ix);
%         plot_info.ax         = gca;
%         plot_info.title      = stat.label{ch_ix};
%         if ch_ix==1; plot_info.legend=1; else plot_info.legend=0; end;
%         
%         % Compute means and variance
%         means = NaN([numel(cond_lab) size(roi_hfa{1}.avg,2)]);
%         var = NaN([numel(cond_lab) size(roi_hfa{1}.avg,2)]);
%         for an_ix = 1:numel(cond_lab)
%             means(an_ix,:) = roi_hfa{an_ix}.avg(ch_ix,:);
%             var(an_ix,:) = squeeze(std(roi_hfa{an_ix}.trial(:,ch_ix,:),[],1)./sqrt(size(roi_hfa{an_ix}.trial,1)))';
%         end
%         % Find significant time periods
%         if sum(stat.mask(ch_ix,:))>0
%             mask_chunks = fn_find_chunks(stat.mask(ch_ix,:));
%             sig_chunks = mask_chunks;
%             sig_chunks(stat.mask(ch_ix,sig_chunks(:,1))==0,:) = [];
%             % If stat and roi_erp aren't on same time axis, adjust sig_chunk indices
%             if (size(stat.time,2)~=size(roi_hfa{1}.time,2)) || (sum(stat.time==roi_hfa{1}.time)~=numel(stat.time))
%                 for chunk_ix = 1:size(sig_chunks,1)
%                     sig_chunks(chunk_ix,1) = find(roi_hfa{1}.time==stat.time(sig_chunks(chunk_ix,1)));
%                     sig_chunks(chunk_ix,2) = find(roi_hfa{1}.time==stat.time(sig_chunks(chunk_ix,2)));
%                 end
%             end
%             fprintf('%i SIGNIFICANT CLUSTERS FOUND, plotting with significance shading...\n',size(sig_chunks,1));
%             fn_plot_ts_error_bar_sig(plot_info,means,var,sig_chunks,event_info,cond_info);
%         else
%             fprintf('NO SIGNIFICANT CLUSTERS FOUND, plotting without significance shading...\n');
%             fn_plot_ts_error_bar(plot_info,means,var,event_info,cond_info);
%         end
%     end
%     clear stat roi_erp
%     
%     % Save figure
%     if save_fig
%         fig_dir = ['/home/knight/hoycw/PRJ_Stroop/results/HFA/' SBJ '/' conditions '/'];
%         if ~exist(fig_dir,'dir')
%             mkdir(fig_dir);
%         end
%         fig_filename = [fig_dir fig_name '.' fig_filetype];
%         fprintf('Saving %s\n',fig_filename);
%         saveas(gcf,fig_filename);
%         %eval(['export_fig ' fig_filename]);
%     end
% end
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
