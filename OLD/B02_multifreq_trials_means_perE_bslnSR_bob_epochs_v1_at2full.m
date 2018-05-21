%% Stroop analysis - High Gamma Power on Single Trials
clear all; %close all
addpath(genpath('/home/knight/hoycw/PRJ_Stroop/scripts/'));
addpath(genpath('/home/knight/hoycw/Apps/export_fig-master/'));
% addpath(genpath('/home/knight/hoycw/Apps/EEGlab/eeglab12_0_2_6b/functions/sigprocfunc/'));

% Analysis Parameters
SBJs       = {'IR35'};%'IR21','IR31','IR32','IR35','IR39'};%'IR39'};%};%
data_ids   = {{'LAC_ft_KLA'}};%{'RC_WM'},{'LAC_WM'},{'IH_CA'},{'LAC_WM'},{'RAC_WM'}};%,'RAC_BP_data'}};%{'RAC_WM'}};
cond_name  = 'CI';              % conditions to average
analyses   = {'ERP','theta','HG'};
%   CI   = only con and inc trial types for all blocks
%   pcon = proportion congruency (con and inc for both mcon and minc)
%   !!!conseq = congruency sequence effects (cC, cI, iC, iI)
epoch_lim = [-200 2000];           % ID of data used for trial rejection (these epochs are NOT used)
buff_lim  = [1000, 1000];            % buffers are for cutting time series and then plotting

% Analysis parameters
HG_type      = 'wideband';                     % only for HG: 'wideband', 'multiband'
event        = 'stim';                          % 'stim'/'resp': event to lock trials
bsln_type    = {'demean','zscore','zscore'};    % 'zscore', 'demean', 'none'
bsln_event   = 's';                    % 's'/'r': event to lock baselining
bsln_lim     = [250, -50];             % ms before and after event
smooth_it    = 1;                      % smooth the result before averaging (0/1)
smooth_freq  = 10;
trial_type   = 'datapad';               % 'nanpad', 'datapad'

% Plotting parameters
save_fig      = 0;
vis_fig       = 'on';                % 'on' or 'off' to determine if plot is shown
sem_alpha     = 0.5;                  % transparency of sem shading (0:1)
clim_perc     = [5 95];              % colormap percentile limits
plot_lim      = [500 500];            % ms to plot before the event and **after RT**
x_step        = 250;                  % step of x tick marks
event_ln_width= 2;
fig_type      = 'png';

%% Process parameters
% Condition Parameters
[cond_lab, cond_colors, cond_style] = fn_condition_label_styles(cond_name);

% Analysis Parameters
all_analysis_id = '';
all_bsln_type = '';
for an_ix = 1:length(analyses)
    [analysis_id{an_ix}, env_it{an_ix}, filt_lim{an_ix}] =...
        fn_B02_analysis_params(analyses{an_ix},HG_type);
    all_analysis_id = strcat(all_analysis_id, analysis_id{an_ix}(1));
    all_bsln_type = strcat(all_bsln_type, bsln_type{an_ix}(1));
end

% Filtering
if smooth_it == 0
    smooth_id = '';
elseif smooth_it == 1
    smooth_id = strcat('_sm',num2str(smooth_freq));
else
    error('smooth_it not in [0,1]');
end

% Epoching/Averaging
epoch_id = [num2str(floor(epoch_lim(1)),'%+d') '.' num2str(floor(epoch_lim(2)),'%+d')];
bsln_id  = ['_' all_bsln_type '.' bsln_event num2str(bsln_lim(1)) '.' num2str(bsln_lim(2))];
if strcmp(trial_type,'nanpad')
    trl_id = 'nan';
elseif strcmp(trial_type,'datapad')
    trl_id = 'data';
else
    error(strcat('Unknown trial averaging type: ',trial_type));
end
if length(bsln_type)~=length(analyses)
    error('Mismatched number of bsln_type and analyses');
end

% Plotting (Y axis)
if (plot_lim(1)>buff_lim(1)) || (plot_lim(2)>buff_lim(2))
    error('plot_lim exceed buff_lim')
end
plot_lim_id = strcat('_xlim',num2str(plot_lim(1)),'.',num2str(plot_lim(2)));
y_scale_id = strcat('_c',num2str(clim_perc(1)),'.',num2str(clim_perc(2)));

%% Analysis Loop
for SBJ_ix = 1:length(SBJs)
    for data_id_ix = 1:length(data_ids{SBJ_ix})
        SBJ = SBJs{SBJ_ix};
        data_id = strcat(SBJ,'_',data_ids{SBJ_ix}{data_id_ix});
        SBJ_dir = strcat('/home/knight/hoycw/PRJ_Stroop/data/',SBJ,'/');
        fig_dir  = strcat('/home/knight/hoycw/PRJ_Stroop/results/multifreq/',SBJ,'/',cond_name,'/');
        if ~exist(fig_dir,'dir')
            mkdir(fig_dir);
        end
        
        %% Load data
        load(strcat(SBJ_dir,'04_proc/',data_id,'.mat'),'data_ecog','header_ecog');
        load(strcat(SBJ_dir,'03_events/',SBJ,'_trial_info_full.mat'));
%         trial_info.resp_onset = manual.trial_info.resp_onset;
%         trial_info.response_time = manual.trial_info.response_time;
        %!!! take out this work around!
        if strcmp(data_id,'IR35_LAC_WM') % Remove LAC4
            data_ecog = data_ecog(1:3,:);
            header_ecog.channel_labels = header_ecog.channel_labels(1:3);
        end
        % need only ok_epochs index, the rest is all in 04_proc file (variables seem same as of 10/10/16...)
        if strcmp(SBJ,'IR32')
            bob_bad_epochs = load(strcat(SBJ_dir,'03_events/',SBJ,'_bob_bad_epochs.mat'),'bad_epochs');
        else
            bob_bad_epochs = load(strcat(SBJ_dir,'03_events/',SBJ,'_bob_bad_epochs.mat'),'bad_epochs');
        end
        bob_bad_epochs = bob_bad_epochs.bad_epochs;
        if strcmp(data_id,'IR39_RAC_BP_data') || strcmp(data_id,'IR39_RAC_ft_KLA')  %!!! SHITTY WORK AROUND!!!
            KLA_ok_epochs = load(strcat(SBJ_dir,'06_epochs/IR39_RAC_WM_',epoch_id,'.mat'),'ok_epochs');
        else
            KLA_ok_epochs = load(strcat(SBJ_dir,'06_epochs/',data_id,'_',epoch_id,'.mat'),'ok_epochs');
        end
        KLA_ok_epochs = KLA_ok_epochs.ok_epochs;
        
        % Correct all info to only good trials
        if header_ecog.sample_rate~=1000
            error('\n!!!WARNING: Sampling rate is not 1 kHz, timing will be off!!\n')
        end
        % Toss epochs that overlap with bad_epochs from Bob
        bad_samples = [];
        time_epoch_idx = fn_epoch_cuts_datapad(1:size(data_ecog,2),trial_info.word_onset,...
            trial_info.resp_onset,buff_lim);
        for bad_epoch_ix = 1:size(bob_bad_epochs,1)
            bad_samples = [bad_samples bob_bad_epochs(bad_epoch_ix,1):bob_bad_epochs(bad_epoch_ix,2)];
        end
        bob_ok_epochs = [];
        for epoch_ix = 1:size(time_epoch_idx,1)
            if isempty(intersect(time_epoch_idx(epoch_ix,:),bad_samples));
                bob_ok_epochs = [bob_ok_epochs epoch_ix];
            end
        end
        ok_epochs = intersect(bob_ok_epochs,KLA_ok_epochs);
        fprintf('Num Bob ok epochs: %i\n',length(bob_ok_epochs));
        fprintf('Num KLA ok epochs: %i\n',length(ok_epochs));
        fprintf('Overlap          : %i\n',length(ok_epochs));
        
        % Toss trials with no RT
        no_RT_epochs = find(isnan(trial_info.resp_onset)==1);
        good_epochs = setdiff(ok_epochs,no_RT_epochs);
        fprintf('Num trials w/o RT: %i\n',sum(no_RT_epochs));
        fprintf('FINAL NUM GOOD TRIALS: %i\n',length(good_epochs));
        % con = 1-3, neu = 4-6, inc = 7-9
        % within those: same, mic, mcon
        % trial_type = NaN(size(trial_info.condition_n));
        word_onset = trial_info.word_onset(good_epochs);
        resp_onset = trial_info.resp_onset(good_epochs);
        RTs        = round(1000*trial_info.response_time(good_epochs)); % converts sec to ms
        [RTs_sorted,RTs_sort_idx] = sort(RTs);
        max_RT     = max(RTs);
        
        % Conditional Behavior
        cond_idx = false([length(cond_lab) length(good_epochs)]);
        RTs_sort_cond = {};
        for cond_ix = 1:length(cond_lab)
            % Get binary condition index
            cond_idx(cond_ix,:) = logical(fn_condition_index(cond_lab{cond_ix},...
                trial_info.condition_n(good_epochs)));
            RTs_sort_cond{cond_ix} = RTs_sorted(cond_idx(cond_ix,:));
        end
            
        % Trial Parameters
        n_chan = size(data_ecog,1);
        if strcmp(event,'stim')
            lock_events = word_onset;
            trial_len = buff_lim(1)+max_RT+buff_lim(2)+1;
        else
            lock_events = resp_onset;
            trial_len = buff_lim(1)+buff_lim(2)+1;
        end
        if strcmp(bsln_event,'s')
            bsln_events = word_onset;
        else
            bsln_events = resp_onset;
        end
        
        %% Analysis Computations
        % Pre-Allocate trials
        trials     = NaN([length(analyses) n_chan length(good_epochs) trial_len]);
        trial_mean = NaN([length(analyses) n_chan length(cond_lab) trial_len]);
        trial_sem  = NaN([length(analyses) n_chan length(cond_lab) trial_len]);
        minmax     = NaN([length(analyses) n_chan length(cond_lab) 2]);   %separate ylim per chan
        for an_ix = 1:length(analyses)
            if ~isempty(filt_lim{an_ix})
                % Filter data
                data_filt = NaN([n_chan length(filt_lim{an_ix}) size(data_ecog,2)]);
                for ch_ix = 1:n_chan
                    for cf_ix = 1:length(filt_lim{an_ix})
                        % Filter data
                        data_filt(ch_ix,cf_ix,:) = fn_EEGlab_bandpass(data_ecog(ch_ix,:),...
                            header_ecog.sample_rate, filt_lim{an_ix}{cf_ix}(1),...
                            filt_lim{an_ix}{cf_ix}(2));
                        % Envelope
                        if env_it{an_ix} ==1
                            data_filt(ch_ix,cf_ix,:) = abs(hilbert(data_filt(ch_ix,cf_ix,:)));
                        end
                    end
                end
                % Average to get single HG time series
                data_proc = mean(data_filt,2);
            else
                % ERP case, no filtering
                data_proc = data_ecog;
            end
            
            %% Format trials
            for ch_ix = 1:size(data_proc,1)
                % Smooth if needed
                if smooth_it == 1
                    data_proc(ch_ix,:) = fn_EEGlab_lowpass(data_proc(ch_ix,:),...
                        header_ecog.sample_rate, smooth_freq);
                end
                
                % Cut trials and baseline all at once
                trials(an_ix,ch_ix,:,:) = fn_epoch_bsln_combo(data_proc(ch_ix,:), lock_events,...
                    resp_onset, buff_lim, trial_type, bsln_events, bsln_lim, bsln_type{an_ix});
            end
            
            %% ERPs by condition
            for cond_ix = 1:length(cond_lab)
                % Average trials matching condition index
                trial_mean(an_ix,:,cond_ix,:) = squeeze(nanmean(trials(an_ix,:,cond_idx(cond_ix,:),:),3));
                % Standarad Error of Mean
                trial_sem(an_ix,:,cond_ix,:) = squeeze(nanstd(trials(an_ix,:,cond_idx(cond_ix,:),:),...
                                                0,3))./sqrt(sum(cond_idx(cond_ix,:)));
                % Plot Limits per electrode for mean
                minmax(an_ix,:,cond_ix,1) = min(trial_mean(an_ix,:,:),[],3);
                minmax(an_ix,:,cond_ix,2) = max(trial_mean(an_ix,:,:),[],3);
            end
            
            
        end
        
        %% Plot ERPs
        % Plotting parameters
        win_on  = buff_lim(1)-plot_lim(1)+1;
        if strcmp(event,'stim')
            win_off = buff_lim(1)+max_RT+plot_lim(2)+1;         % Capture even longest trials
        else
            win_off = buff_lim(1)+plot_lim(2)+1;
        end
        x_lab   = -plot_lim(1):x_step:max_RT+plot_lim(2);
        ylims = zeros([1 2]);
        clims = zeros([1 2]);
        
        for ch_ix = 1:n_chan
            % Plot Condition ERPs
            out_file = [data_id '_' header_ecog.channel_labels{ch_ix} '_' all_analysis_id...
                '_trials_' event '_' cond_name '_Bob.ep' epoch_id '_' trl_id...
                '_bsln' bsln_id smooth_id plot_lim_id y_scale_id];
            fig_height = length(analyses)/3;
            figure('Name',out_file,'units','normalized',...
                'outerposition',[0 0 1 fig_height],'Visible',vis_fig);
            for an_ix = 1:length(analyses)
                ylims(1) = min(minmax(an_ix,ch_ix,:,1),[],3);
                ylims(2) = max(minmax(an_ix,ch_ix,:,2),[],3);
                clims(1) = prctile(reshape(trials(an_ix,ch_ix,:,:),...
                    [1 size(trials,3)*size(trials,4)]),clim_perc(1));
                clims(2) = prctile(reshape(trials(an_ix,ch_ix,:,:),...
                    [1 size(trials,3)*size(trials,4)]),clim_perc(2));
                
                % Plot Event Locked Average
                subplot(length(analyses),length(cond_lab)+1,...
                    an_ix*(length(cond_lab)+1)-length(cond_lab));
                hold on;
                ebars = {};
                main_lines = [];
                for cond_ix = 1:length(cond_lab)
                    ebar{cond_ix} = shadedErrorBar([],...
                        trial_mean(an_ix,ch_ix,cond_ix,win_on:win_off),...
                        trial_sem(an_ix,ch_ix,cond_ix,win_on:win_off),...
                        {'Color',[cond_colors{cond_ix}],...
                        'LineStyle',cond_style{cond_ix}},sem_alpha);
                    cond_lab_legend{cond_ix} = [cond_lab{cond_ix} '-' num2str(sum(cond_idx(cond_ix,:)))];
                    main_lines = [main_lines ebar{cond_ix}.mainLine];
                end
                ylim(ylims);
                event_line = line([plot_lim(1) plot_lim(1)],ylim,...
                    'LineWidth',event_ln_width,'Color','k');
                % Axes and Labels
                ax = gca;
                ax.XLim = [0,win_off];
                ax.XTick = 0:x_step:win_off;
                ax.XTickLabel = x_lab;
                main_lines = [main_lines event_line];
                title(strcat(header_ecog.channel_labels(ch_ix), ': ', analysis_id{an_ix},...
                    ' (', event,'-locked)'));
                if an_ix==1
                    legend(main_lines,cond_lab_legend{:},event,'Location','southeast');%,'inc-con'
                end
                
                % Plot Single Trials Per Condition
                for cond_ix = 1:length(cond_lab)
                    subplot(length(analyses),length(cond_lab)+1,...
                        an_ix*(length(cond_lab)+1)-length(cond_lab)+cond_ix); hold on;
                    
                    imagesc(squeeze(trials(an_ix,ch_ix,RTs_sort_idx(cond_idx(cond_ix,:)),win_on:win_off)));
                    set(gca,'YDir','normal');
                    if strcmp(event,'stim')
                        scat = scatter(RTs_sort_cond{cond_ix}+plot_lim(1),1:sum(cond_idx(cond_ix,:)),...
                            'MarkerFaceColor',[cond_colors{cond_ix}],...
                            'MarkerEdgeColor','k');
                    end
                    ylim([1 sum(cond_idx(cond_ix,:))]);
                    event_line = line([plot_lim(1) plot_lim(1)],ylim,...
                        'LineWidth',event_ln_width,'Color','k');
                    % Axes and Labels
                    cbar = colorbar;
                    caxis(clims);
                    ax = gca;
                    ax.XLim = [0,win_off];
                    ax.XTick = 0:x_step:win_off;
                    ax.XTickLabel = x_lab;
                    
                    title(strcat(header_ecog.channel_labels(ch_ix), ':', analysis_id{an_ix},...
                        ',',cond_lab{cond_ix},' (Single Trial, ',event,'-locked)'));
                end
            end
            
            % Save plots
            out_filename = [fig_dir out_file '.' fig_type];
            if save_fig ==1
                fprintf('Saving %s\n',out_filename);
                eval(['export_fig ' out_filename]);
            end
            % saveas(gcf,out_filename);
        end
        
        clear data_ecog header_ecog trial_info SBJ data_id good_epochs ok_epochs no_RT_epochs word_onset...
            resp_onset RTs max_RT data_filt data_proc...
            trials cond_idx an_ix ch_ix trial_mean trial_sem out_file...
            minmax y_bound
    end
end
