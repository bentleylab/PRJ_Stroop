%% Stroop analysis - High Gamma Power
clear all; %close all
addpath(genpath('/home/knight/hoycw/PRJ_Stroop/scripts/'));
addpath(genpath('/home/knight/hoycw/Apps/export_fig-master/'));
addpath(genpath('/home/knight/hoycw/Apps/EEGlab/eeglab12_0_2_6b/functions/sigprocfunc/'));

% Analysis Parameters
SBJs       = {'IR35'};%'IR21','IR31','IR35','IR39'};%'IR32','IR39'};%};%
data_ids   = {{'LAC_WM'}};%{'RC_WM'},{'LAC_WM'},{'LAC_WM'},{'RAC_WM','RAC_BP_data'}};%{'IH_CA'},{'RAC_WM'}};%,};%
cond_name  = 'pcon';              % conditions to average
%   CNI  = basic trial types (con, neu, inc) for all blocks
%   CI   = only con and inc trial types for all blocks
%   pcon = proportion congruency (con and inc for both mcon and minc)
%   !!!conseq = congruency sequence effects (cC, cI, iC, iI)
epoch_lim  = [-200 2000];           % ID of data used for trial rejection (these epochs are NOT used)
buff_lim   = [1000, 1000];            % buffers are for cutting time series and then plotting
sig_test   = 0;                      % run significance testing via permutations across conditions
sig_lim_s  = [0 500];               % window before stim and **after STIM** for significance testing
sig_lim_r  = [500 0];               % window before and after RESP for significance testing
sig_len    = 50;                    % window size for significance testing

% Analysis parameters
phs_fband    = 'theta';                  % 'ERP', 'HG', 'theta'
amp_fband    = 'HG';                  % 'ERP', 'HG', 'theta'
HG_type      = 'wideband';            % only for HG: 'wideband', 'multiband'
n_boots      = 1000;                  % n times to shuffle PLVs for significance testing

% Don't mess with these:
smooth_it    = 0;%!!!keep 0!                      % smooth the result before averaging (0/1)
smooth_freq  = 10;
% hp_freq      = 0.5;                    % high pass filter ONLY for ERPs with baseline='none'
trial_type   = 'datapad';               % 'nanpad', 'datapad'

% Plotting parameters
save_fig      = 0;
vis_fig       = 'on';                % 'on' or 'off' to determine if plot is shown
% plot_sem      = 1;                    % plot the standard error of the mean traces? 0/1
sem_alpha     = 0.5;                  % transparency of sem shading (0:1)
y_symmetry    = 0;                    % symmetric y-axis? 0/1
y_common      = 1;                    % plots share y_lim? 0=all free;1=S/R common;2=all common
plot_lim_s    = [250 750];            % time window to plot around stimulus
plot_lim_r    = [500 500];            % time window to plot around response
x_step        = 250;                  % step of x tick marks
% plot_prop_con = 0;
% manual_cond_colors = {'c', 'k', 'r'}; %condition colors if not auto assigned below
% cond_lines    = {'-', '-', '--','--'};    % colors for cond_lab plotting
% block_colors  = {'r', 'g', 'b'};    % colors for [mcon, same, minc]
% prop_con_lab  = {'con-mcon', 'con-same', 'con-minc', 'neu-mcon', 'neu-same',...
%     'neu-minc', 'inc-mcon', 'inc-same', 'inc-minc'};
fig_type      = 'eps';

%% Process parameters
% Condition Types
[cond_lab, cond_colors, cond_style] = fn_condition_label_styles(cond_name);
if (length(cond_lab)~=2) && (sig_test)
    error('Significance testing loop only written for two conditions! (e.g., CI)');
end

% High Gamma Extraction
[phs_fband_id, ~, phs_filt_lim] = fn_B02_analysis_params(phs_fband,HG_type);
[amp_fband_id, ~, amp_filt_lim] = fn_B02_analysis_params(amp_fband,HG_type);
if length(phs_filt_lim)>1
    error('Phase data can only have one center frequency!')
end

% Filtering
if smooth_it == 0
    smooth_id = '';
elseif smooth_it == 1
    smooth_id = strcat('_sm',num2str(smooth_freq));
% elseif smooth_it == 2
%     smooth_id = ['_sm.avgs.' num2str(smooth_freq)];
else
    error('smooth_it not in [0,1]');
end

% Epoching/Averaging
epoch_id = [num2str(floor(epoch_lim(1)),'%+d') '.' num2str(floor(epoch_lim(2)),'%+d')];
if sig_test
    sig_id   = ['_sigS' num2str(sig_lim_s(1)) '.' num2str(sig_lim_s(2))...
        '_R' num2str(sig_lim_r(1)) '.' num2str(sig_lim_r(2))];
else
    sig_id = '';
end
if strcmp(trial_type,'nanpad')
    trl_id = 'nan';
elseif strcmp(trial_type,'datapad')
    trl_id = 'data';
else
    error(strcat('Unknown trial averaging type: ',trial_type));
end

% Plotting (Y axis)
y_scale_id = strcat('_yc',num2str(y_common),'s',num2str(y_symmetry));

%% Analysis Loop
for SBJ_ix = 1:length(SBJs)
    for data_id_ix = 1:length(data_ids{SBJ_ix})
        SBJ = SBJs{SBJ_ix};
        data_id = strcat(SBJ,'_',data_ids{SBJ_ix}{data_id_ix});
        SBJ_dir = strcat('/home/knight/hoycw/PRJ_Stroop/data/',SBJ,'/');
        fig_dir  = strcat('/home/knight/hoycw/PRJ_Stroop/results/PAC/',SBJ,'/',cond_name,'/');
        if ~exist(fig_dir,'dir')
            mkdir(fig_dir);
        end
        
        %% Load data
        load(strcat(SBJ_dir,'04_proc/',data_id,'.mat'),'data_ecog','header_ecog');
        load(strcat(SBJ_dir,'03_events/',data_id,'_trial_info.mat'));
%         trial_info.resp_onset = manual.trial_info.resp_onset;
%         trial_info.response_time = manual.trial_info.response_time;
        %!!! take out this work around!
        if strcmp(data_id,'IR35_LAC_WM') % Remove LAC4
            data_ecog = data_ecog(1:3,:);
            header_ecog.channel_labels = header_ecog.channel_labels(1:3);
        end
        % need only ok_epochs index, the rest is all in 04_proc file (variables seem same as of 10/10/16...)
        if strcmp(SBJ,'IR32')
            bob_bad_epochs = load(strcat(SBJ_dir,'02_preproc/',SBJ,'_bad_epochs.mat'),'bad_epochs');
        else
            bob_bad_epochs = load(strcat(SBJ_dir,'02_preproc/',SBJ,'_bad_epochs.mat'),'bad_epochs');
        end
        bob_bad_epochs = bob_bad_epochs.bad_epochs;
        if strcmp(data_id,'IR39_RAC_BP_data')   %!!! SHITTY WORK AROUND!!!
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
        max_RT     = max(RTs);
        
        % Conditional Behavior
        cond_idx = false([length(cond_lab) length(good_epochs)]);
        for cond_ix = 1:length(cond_lab)
            % Get binary condition index
            cond_idx(cond_ix,:) = logical(fn_condition_index(cond_lab{cond_ix},...
                trial_info.condition_n(good_epochs)));
        end
        % Plot/Analysis Parameters
        n_chan = size(data_ecog,1);
        trial_len_s = buff_lim(1)+max_RT+buff_lim(2)+1;
        trial_len_r = buff_lim(1)+buff_lim(2)+1;
        
        %% Phase-Amplitude Coupling Calculation
        % Filter phase data
        p_phs_data = NaN(size(data_ecog));
        a_phs_data = NaN(size(data_ecog));
        for ch_ix = 1:n_chan
            % Filter phase data, extract phase time series
            p_phs_data(ch_ix,:) = fn_EEGlab_bandpass(data_ecog(ch_ix,:),...
                header_ecog.sample_rate,phs_filt_lim{1}(1),phs_filt_lim{1}(2));
            p_phs_data(ch_ix,:) = angle(hilbert(p_phs_data(ch_ix,:)));    %ok not to flip b/c single time series
        
            % Filter amplitude data
            amp_tmp = NaN([length(amp_filt_lim) size(data_ecog,2)]);
            for cf_ix = 1:length(amp_filt_lim)
                % Filter data
                amp_tmp(cf_ix,:) = fn_EEGlab_bandpass(data_ecog(ch_ix,:),...
                    header_ecog.sample_rate,amp_filt_lim{cf_ix}(1),amp_filt_lim{cf_ix}(2));
                % Envelope
                amp_tmp(cf_ix,:) = abs(hilbert(amp_tmp(cf_ix,:)));
            end
            % Average over all envelopes (e.g., multiband if necessary)
            a_phs_data(ch_ix,:) = mean(amp_tmp,1);
            % Filter amplitude data to phase frequency, extract phase
            a_phs_data(ch_ix,:) = fn_EEGlab_bandpass(a_phs_data(ch_ix,:),...
                header_ecog.sample_rate,phs_filt_lim{cf_ix}(1),phs_filt_lim{cf_ix}(2));
            a_phs_data(ch_ix,:) = angle(hilbert(a_phs_data(ch_ix,:)));
        end
        phs_diff = p_phs_data-a_phs_data;
%         % Smooth
%         if smooth_it == 1
%             for ch_ix = 1:size(data_ecog,1)
%                 data_hg(ch_ix,:) = eegfilt(data_hg(ch_ix,:), header_ecog.sample_rate, [], smooth_freq);
%             end
%         end
        
        %% Format trials
        % Cut Trials
        trials_stim = NaN([size(phs_diff,1),length(good_epochs),trial_len_s]);
        trials_resp = NaN([size(phs_diff,1),length(good_epochs),trial_len_r]);
        for ch_ix = 1:n_chan
            % Cut trials and baseline all at once
            trials_stim(ch_ix,:,:) = fn_epoch_cuts_datapad(phs_diff(ch_ix,:),word_onset,resp_onset,buff_lim);
            trials_resp(ch_ix,:,:) = fn_epoch_cuts_datapad(phs_diff(ch_ix,:),resp_onset,resp_onset,buff_lim);
        end
        
%         %% Calculate PLV time series (PAC)
%         PLV_time = fn_PAC_PLVenvPh_ePair_time(phase1, phase2, starts, ends)
%         PLV_time = abs(mean(exp(1i*angle_diff),1));
        %% ERPs by condition
        % separate ylim per chan
        mean_stim = NaN([n_chan length(cond_lab) trial_len_s]);
        mean_resp = NaN([n_chan length(cond_lab) trial_len_r]);
        minmax    = NaN([n_chan length(cond_lab) 2 2]);    %chan, cond, s/r, low/hi
        for cn_ix = 1:length(cond_lab)
            % Average trials matching condition index
            mean_stim(:,cn_ix,:) = squeeze(abs(nanmean(exp(1i*trials_stim(:,cond_idx(cn_ix,:),:)),2)));
            mean_resp(:,cn_ix,:) = squeeze(abs(nanmean(exp(1i*trials_resp(:,cond_idx(cn_ix,:),:)),2)));
            
            % Get bounds across channels for plotting on equal axes
%             if y_common
%                 ch_con_min = min([min(min(eval(['mean_stim.' cond_lab{ix}]),[],2))...
%                     min(min(eval(['mean_resp.' cond_lab{ix}])))]);
%                 ch_con_max = max([max(max(eval(['mean_stim.' cond_lab{ix}])))...
%                     max(max(eval(['mean_resp.' cond_lab{ix}])))]);
%             else
            minmax(:,cn_ix,:,1) = [min(mean_stim(:,cn_ix,:),[],3),...
                min(mean_resp(:,cn_ix,:),[],3)];
            minmax(:,cn_ix,:,2) = [max(mean_stim(:,cn_ix,:),[],3),...
                max(mean_resp(:,cn_ix,:),[],3)];
%             !!! now fix thisif (abs(ch_con_min) > minmax) || (ch_con_max > minmax)
%                 minmax = max([abs(ch_con_min) ch_con_max]);
%             end
            
        end
        
        %% Calculate pval difference between these
        if sig_test
            n_wins = int(sig_len
            for wn_ix = 1:size(sig_wins)
                
                plv_pval = NaN([n_chan 2]);             % channel, s/r
                pval_win_s = [buff_lim(1)-sig_lim_s(1):buff_lim(1)+sig_lim_s(2)];
                pval_win_r = [buff_lim(1)-sig_lim_r(1):buff_lim(1)+sig_lim_r(2)];
                for ch_ix = 1:n_chan
                    % Real PLV difference
                    plv_diff_s = mean(mean_stim(ch_ix,1,pval_win_s))-mean(mean_stim(ch_ix,2,pval_win_s));
                    plv_diff_r = mean(mean_resp(ch_ix,1,pval_win_r))-mean(mean_resp(ch_ix,2,pval_win_r));
                    % Shuffle Conditions to get null difference distribution & pval
                    [plv_pval(ch_ix,1),plv_surr_hist_s] = fn_PLV_time_surr_epoch_diff_rand_PLV_ts(plv_diff_s,...
                        squeeze(mean_stim(ch_ix,1,pval_win_s))',squeeze(mean_stim(ch_ix,2,pval_win_s))',n_boots);
                    [plv_pval(ch_ix,2),plv_surr_hist_r] = fn_PLV_time_surr_epoch_diff_rand_PLV_ts(plv_diff_r,...
                        squeeze(mean_resp(ch_ix,1,pval_win_r))',squeeze(mean_resp(ch_ix,2,pval_win_r))',n_boots);
                end
            end
        end
        
        %% Plot ERPs
        % Plotting parameters
        x_lab_stim   = -plot_lim_s(1):x_step:plot_lim_s(2);
        x_lab_resp   = -plot_lim_r(1):x_step:plot_lim_r(2);
        stim_win_on  = buff_lim(1)-plot_lim_s(1)+1;
        resp_win_on  = buff_lim(1)-plot_lim_r(1)+1;
        stim_win_off = buff_lim(1)+plot_lim_s(2)+1;
        resp_win_off = buff_lim(1)+plot_lim_r(2)+1;
        
        ylims = zeros([size(trials_stim,1) 2 2]);
        ylims(:,:,1) = min(minmax(:,:,:,1),[],2);
        ylims(:,:,2) = max(minmax(:,:,:,2),[],2);
        if y_common==1
            min_SR = min(ylims(:,:,1),[],2);
            ylims(:,:,1) = horzcat(min_SR,min_SR);
            max_SR = max(ylims(:,:,2),[],2);
            ylims(:,:,2) = horzcat(max_SR,max_SR);
        elseif y_common==2
            min_all = min(min(ylims(:,:,1),[],2));
            ylims(:,:,1) = repmat(min_all,size(ylims,1),size(ylims,2));
            max_all = max(max(ylims(:,:,2),[],2));
            ylims(:,:,2) = repmat(max_all,size(ylims,1),size(ylims,2));          
        end
%         if y_symmetry==1
%         end
%         y_bound      = ceil(minmax/y_bound_chunk)*y_bound_chunk;
        
        % Plot Condition ERPs
        out_file = [data_id '_PAC' phs_fband_id '.' amp_fband_id '_' cond_name...
            '_Bob.ep' epoch_id sig_id '_' trl_id smooth_id y_scale_id];
        if n_chan > 2
            fig_height = 1;
        else
            fig_height = n_chan/3;
        end
        fig = figure('Name',out_file,'units','normalized',...
            'outerposition',[0 0 1 fig_height],'Visible',vis_fig);
        for ch_ix = 1:n_chan
            % Plot Stim Locked
            subplot(n_chan,2,ch_ix*2-1); hold on;
%             ebar = {};
            for cn_ix = 1:length(cond_lab)
%                 if n_chan==1
% %                     eval(['ebar{cn_ix} = shadedErrorBar([],'...
% %                         'mean_stim.' cond_lab{cn_ix} '(stim_win_on:stim_win_off),'...
% %                         'sem_stim.' cond_lab{cn_ix} '(stim_win_on:stim_win_off),'...
% %                         '{''Color'',[cond_colors{cn_ix}],''LineStyle'',cond_style{cn_ix}},sem_alpha);']);
%                     plot(mean_stim(ch_ix,cn_ix,stim_win_on:stim_win_off),'Color',cond_colors{cn_ix});
%                 else
%                     eval(['ebar{cn_ix} = shadedErrorBar([],'...
%                         'mean_stim.' cond_lab{cn_ix} '(ch,stim_win_on:stim_win_off),'...
%                         'sem_stim.' cond_lab{cn_ix} '(ch,stim_win_on:stim_win_off),'...
%                         '{''Color'',[cond_colors{cn_ix}],''LineStyle'',cond_style{cn_ix}},sem_alpha);']);
                plot(squeeze(mean_stim(ch_ix,cn_ix,stim_win_on:stim_win_off)),...
                    'Color',[cond_colors{cn_ix}],'LineStyle',cond_style{cn_ix});
%                 end
            end
            %     plot(mean_diff_i_c_stim(ch,:),'g');
            ylim([ylims(ch_ix,1,1) ylims(ch_ix,1,2)]);
            stim_line = line([plot_lim_s(1) plot_lim_s(1)],ylim);%, 'k--');
            rt_line = line([plot_lim_s(1)+nanmean(RTs), plot_lim_s(1)+nanmean(RTs)],ylim,'LineStyle','--');
            % Axes and Labels
            ax = gca;
            ax.XLim = [0,sum(plot_lim_s)];
            ax.XTick = 0:x_step:sum(plot_lim_s);
            ax.XTickLabel = x_lab_stim;
            if sig_test
                title_str = strcat('Stim Locked: ',header_ecog.channel_labels(ch_ix),...
                '-p=',num2str(plv_pval(ch_ix,1)));
            else
                title_str = strcat('Stim Locked: ',header_ecog.channel_labels(ch_ix));
            end
            title(title_str);
            for cn_ix = 1:length(cond_lab)
                cond_lab_legend{cn_ix} = [cond_lab{cn_ix} '-' num2str(sum(cond_idx(cn_ix,:)))];
            end
            if ch_ix==1
                legend(cond_lab_legend{:},'stim','mean RT','Location','southeast');%,'inc-con'
            end
            
            % Plot RT locked
            subplot(n_chan,2,ch_ix*2); hold on;
%             ebar = {};
            for cn_ix = 1:length(cond_lab)
%                 if n_chan==1
%                     eval(['ebar{cn_ix} = shadedErrorBar([],'...
%                         'mean_resp.' cond_lab{cn_ix} '(resp_win_on:resp_win_off),'...
%                         'sem_resp.' cond_lab{cn_ix} '(resp_win_on:resp_win_off),'...
%                         '{''Color'',[cond_colors{cn_ix}],''LineStyle'',cond_style{cn_ix}},sem_alpha);']);
%                     eval(['plot(mean_resp.' cond_lab{cn_ix} '(resp_win_on:resp_win_off),''Color'',cond_colors{cn_ix});']);
%                 else
%                     eval(['ebar{cn_ix} = shadedErrorBar([],'...
%                         'mean_resp.' cond_lab{cn_ix} '(ch,resp_win_on:resp_win_off),'...
%                         'sem_resp.' cond_lab{cn_ix} '(ch,resp_win_on:resp_win_off),'...
%                         '{''Color'',[cond_colors{cn_ix}],''LineStyle'',cond_style{cn_ix}},sem_alpha);']);
                plot(squeeze(mean_resp(ch_ix,cn_ix,resp_win_on:resp_win_off)),...
                    'Color',[cond_colors{cn_ix}],'LineStyle',cond_style{cn_ix});
%                 end
            end
            ylim([ylims(ch_ix,2,1) ylims(ch_ix,2,2)]);
            rt_line_r = line([plot_lim_r(1) plot_lim_r(1)],ylim);%, 'k--');
            ax = gca;
            ax.XLim = [0,sum(plot_lim_r)];
            ax.XTick = 0:x_step:sum(plot_lim_r);
            ax.XTickLabel = x_lab_resp;
            if sig_test
                title_str = strcat('Resp Locked: ',header_ecog.channel_labels(ch_ix),...
                '-p=',num2str(plv_pval(ch_ix,2)));
            else
                title_str = strcat('Resp Locked: ',header_ecog.channel_labels(ch_ix));
            end
            title(title_str);
            if ch_ix==1
                legend(cond_lab_legend{:},'resp','Location','southeast');%,'inc-con'
            end
        end
        
        % Save plots
        out_filename = [fig_dir out_file '.' fig_type];
        if save_fig ==1
            fprintf('Saving %s\n',out_filename);
            eval(['export_fig ' out_filename]);
        end
        % saveas(gcf,out_filename);
        
        clear data_ecog a_phs_data phs_diff header_ecog trial_info SBJ data_id...
            amp_tmp n_chan good_epochs ok_epochs no_RT_epochs word_onset...
            resp_onset RTs max_RT avg_ch_ix avg_ch_lab data_filt data_hg...
            trials_stim trials_resp cond_idx mean_stim mean_resp out_file...
            minmax ylims y_bound
    end
end
