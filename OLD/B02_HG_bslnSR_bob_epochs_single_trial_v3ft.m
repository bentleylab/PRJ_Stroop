%% Stroop analysis - High Gamma Power on Single Trials
clear all; %close all
addpath(genpath('/home/knight/hoycw/PRJ_Stroop/scripts/'));
addpath(genpath('/home/knight/hoycw/Apps/export_fig-master/'));
addpath(genpath('/home/knight/hoycw/Apps/EEGlab/eeglab12_0_2_6b/functions/sigprocfunc/'));

% Analysis Parameters
SBJs       = {'IR39'};%'IR21','IR31','IR35','IR39'};%'IR32','IR39'};%};%
data_ids   = {{'RAM*','RIN*','RAC*','LOF*','ROF*'}};
%{'RC_WM'},{'LAC_WM'},{'LAC_WM'},{'RAC_WM'}};%{'IH_CA'},{'RAC_WM'}};%,};%
cond_lab   = {'CI'};              % conditions to average
%   CI   = only con and inc trial types for all blocks
%   pcon = proportion congruency (con and inc for both mcon and minc)
%   !!!conseq = congruency sequence effects (cC, cI, iC, iI)
epoch_lim = [-200 2000];           % ID of data used for trial rejection (these epochs are NOT used)
buff_lim  = [500, 500];            % buffers are for cutting time series and then plotting

% Analysis parameters
analysis     = 'HG';                  % 'ERP', 'HG', 'theta'
filt_type    = 'wideband';            % only for HG: 'wideband', 'multiband'
event        = 'stim';                 % 'stim'/'resp': event to lock trials
bsln_type    = 'zscore';               % 'zscore', 'demean', 'none'
bsln_event   = 's';                    % 's'/'r': event to lock baselining
bsln_lim     = [400, -200];             % ms before and after event
smooth_it    = 1;                      % smooth the result before averaging (0/1)
smooth_freq  = 10;
trial_type   = 'datapad';               % 'nanpad', 'datapad'
% avg_channels = 0;                      % average across channels

% Plotting parameters
save_fig      = 1;
vis_fig       = 'off';                % 'on' or 'off' to determine if plot is shown
sem_alpha     = 0.5;                  % transparency of sem shading (0:1)
clim_perc     = [5 95];              % colormap percentile limits
plot_lim      = [500 500];            % ms to plot before the event and **after RT**
x_step        = 250;                  % step of x tick marks
fig_type      = 'png';
% plot_prop_con = 0;
% manual_cond_colors = {'c', 'k', 'r'}; %condition colors if not auto assigned below
% cond_lines    = {'-', '-', '--','--'};    % colors for cond_lab plotting
% block_colors  = {'r', 'g', 'b'};    % colors for [mcon, same, minc]
% prop_con_lab  = {'con-mcon', 'con-same', 'con-minc', 'neu-mcon', 'neu-same',...
%     'neu-minc', 'inc-mcon', 'inc-same', 'inc-minc'};

%% Process parameters
% Condition Types
if length(cond_lab) == 1
    switch cond_lab{1}
        case 'pcon'
            cond_id = 'pcon';
            cond_lab = {'con_mcon', 'con_minc', 'inc_mcon', 'inc_minc'};
            cond_colors = {[0 0 0], [0 0 0], [228,26,28]./256, [228,26,28]./256};    % colors for cond_lab plotting
            cond_style = {'-', '--', '-','--'};    % colors for cond_lab plotting
        case 'CNI'
            cond_id = 'CNI';
            cond_lab = {'con', 'neu', 'inc'};
            cond_colors = {[55,126,184]./256, [0 0 0], [228,26,28]./256};
            cond_style = {'-', '-', '-'};    % colors for cond_lab plotting
        case 'CI'
            cond_id = 'CI';
            cond_lab = {'con', 'inc'};
            cond_colors = {[55,126,184]./256, [228,26,28]./256};
            cond_style = {'-', '-'};    % colors for cond_lab plotting
        case 'conseq'
            cond_id = 'conseq';
        otherwise
            error(strcat('Only one, unrecognized condition offered: ',cond_lab{:}));
    end
else
    cond_id = 'cst';
    for c_ix = 1:length(cond_lab)
        cond_id = [cond_id fn_convert_condition_lab2num(cond_lab{c_ix})];
    end
    cond_colors = manual_cond_colors;
    if length(cond_lab)~=length(cond_colors)
        error('Mismatched condition labels and colors');
    end
end

% High Gamma Extraction
if strcmp(analysis,'HG')
    env_it = 1;
    env_id = '_env';
    if strcmp(filt_type,'wideband')
        analysis_id = 'HGw';
        filt_cf = [110];
        filt_bw = 80;
    elseif strcmp(filt_type,'multiband')
        analysis_id = 'HGm';
        filt_cf = [70:10:150];
        filt_bw = 10;
    else
        error(strcmp('Unknown HG filter type:',filt_type));
    end
    y_bound_chunk = 0.02;       % extend ylim by chunks of this size
elseif strcmp(analysis,'theta')
    env_it = 1;
    env_id = '_env';
    analysis_id = 'theta';
    filt_cf = [6];
    filt_bw = 4;
    y_bound_chunk = 0.5;
elseif strcmp(analysis,'ERP')
    env_it = 0;
    env_id = '';
    analysis_id = 'ERP';
    filt_cf = [];
    y_bound_chunk = 5;
else
    error(strcat('Unknown analysis selection: ',analysis));
end
for f_ix = 1:length(filt_cf)
    filt_lim{f_ix} = [filt_cf(f_ix)-filt_bw/2 filt_cf(f_ix)+filt_bw/2];
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
bsln_id  = ['_' bsln_type(1) bsln_event num2str(bsln_lim(1)) '.' num2str(bsln_lim(2))];
if strcmp(trial_type,'nanpad')
    trl_id = 'nan';
elseif strcmp(trial_type,'datapad')
    trl_id = 'data';
else
    error(strcat('Unknown trial averaging type: ',trial_type));
end
% if avg_channels == 1
%     grand_id = '_ch.avg';
% else
grand_id = '';
% end

% Plotting (Y axis)
if (plot_lim(1)>buff_lim(1)) || (plot_lim(2)>buff_lim(2))
    error('plot_lim exceed buff_lim')
end
y_scale_id = strcat('_c',num2str(clim_perc(1)),'.',num2str(clim_perc(2)));

%% Analysis Loop
for SBJ_ix = 1:length(SBJs)
    for data_id_ix = 1:length(data_ids{SBJ_ix})
        SBJ = SBJs{SBJ_ix};
        data_id = strcat(SBJ,'_ft_preproc');
        SBJ_dir = strcat('/home/knight/hoycw/PRJ_Stroop/data/',SBJ,'/');
        if strcmp(analysis,'ERP')
            fig_dir  = strcat('/home/knight/hoycw/PRJ_Stroop/results/ERPs/',SBJ,'/',cond_id,'/');
        else
            fig_dir  = strcat('/home/knight/hoycw/PRJ_Stroop/results/power/',SBJ,'/',cond_id,'/');
        end
        if ~exist(fig_dir,'dir')
            mkdir(fig_dir);
        end
        
        %% Load data
        cfg = [];
        cfg.channel = data_ids{SBJ_ix}{data_id_ix};
        cfg.datafile = strcat(SBJ_dir,'02_preproc/',data_id,'.mat');
        data_ft = ft_preprocessing(cfg);
        data_ecog = data_ft.trial{1};
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
        if strcmp(data_id,'IR39_RAC_BP_data') || strcmp(data_id,'IR39_ft_preproc')%!!! SHITTY WORK AROUND!!!
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
        [RTs_sorted,RT_sort_idx] = sort(RTs);
        max_RT     = max(RTs);
        
%         if avg_channels == 1
%             avg_ch_ix = fn_ch_to_combine_ix(data_id);
%             if isempty(avg_ch_ix)
%                 avg_ch_ix = 1:length(header_ecog.channel_labels);
%             end
%             avg_ch_lab = strcat('{',header_ecog.channel_labels{avg_ch_ix(1)});
%             if length(avg_ch_ix) > 1
%                 for ch_lab_ix = 2:length(avg_ch_ix)
%                     avg_ch_lab = [avg_ch_lab ',' header_ecog.channel_labels{avg_ch_ix(ch_lab_ix)}];
%                 end
%             end
%             avg_ch_lab = strcat(avg_ch_lab,'}');
%         end
               
        %% High Gamma Extraction
        if strcmp(analysis,'HG') || strcmp(analysis,'theta')
            % Filter data
            data_filt = NaN([size(data_ecog,1) length(filt_cf) size(data_ecog,2)]);
            for ch_ix = 1:size(data_ecog,1)
                for cf_ix = 1:length(filt_cf)
                    % Filter data
                    data_filt(ch_ix,cf_ix,:) = eegfilt(data_ecog(ch_ix,:), ...
                        header_ecog.sample_rate, filt_lim{cf_ix}(1), []);
                    data_filt(ch_ix,cf_ix,:) = eegfilt(data_filt(ch_ix,cf_ix,:), ...
                        header_ecog.sample_rate, [], filt_lim{cf_ix}(2));
                    % Envelope
                    if env_it ==1
                        data_filt(ch_ix,cf_ix,:) = abs(hilbert(data_filt(ch_ix,cf_ix,:)));
                    end
                    % Baseline?
                end
            end
            % Average to get single HG time series
            data_hg = mean(data_filt,2);
        else
            % ERP case, no filtering
            data_hg = data_ecog;
        end
        
        % Smooth
        if smooth_it == 1
            for ch_ix = 1:size(data_ecog,1)
                data_hg(ch_ix,:) = eegfilt(data_hg(ch_ix,:), header_ecog.sample_rate, [], smooth_freq);
            end
        end
        
        %% Format trials
        % Cut Trials
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
        trials = NaN([size(data_hg,1),length(good_epochs),trial_len]);
        % Cut trials and baseline all at once
        for ch = 1:size(data_hg,1)
            trials(ch,:,:) = fn_epoch_bsln_combo(data_hg(ch,:),lock_events,resp_onset,buff_lim,trial_type,...
                                                        bsln_events,bsln_lim,bsln_type);
        end
        
        %% ERPs by condition
        cond_idx = {};
        RTs_sort_cond = {};
        trial_mean = {};
        trial_sem = {};
        minmax = zeros([size(trials,1) length(cond_lab) 2]);    %separate ylim per chan
        for ix = 1:length(cond_lab)
            % Get binary condition index
            cond_idx{ix} = logical(fn_condition_index(cond_lab{ix}, trial_info.condition_n(good_epochs)));
            RTs_sort_cond{ix} = RTs_sorted(cond_idx{ix});
            
            % Average trials matching condition index
            trial_mean{ix} = squeeze(nanmean(trials(:,cond_idx{ix}==1,:),2));
            % Standarad Error of Mean
            trial_sem{ix} = squeeze(nanstd(trials(:,cond_idx{ix}==1,:),0,2))./sqrt(sum(cond_idx{ix}));
            % Plot Limits per electrode for mean
            minmax(:,ix,1) = min(trial_mean{ix},[],2);
            minmax(:,ix,2) = max(trial_mean{ix},[],2);
        end
        
        %% Plot ERPs
        % Plotting parameters
        n_chan = size(data_ecog,1);
        win_on  = buff_lim(1)-plot_lim(1)+1;
        if strcmp(event,'stim')
            win_off = buff_lim(1)+max_RT+plot_lim(2)+1;         % Capture even longest trials
        else
            win_off = buff_lim(1)+plot_lim(2)+1;
        end
        x_lab   = -plot_lim(1):x_step:max_RT+plot_lim(2);
        ylims = zeros([size(trials,1) 2]);
        ylims(:,1) = min(minmax(:,:,1),[],2);
        ylims(:,2) = max(minmax(:,:,2),[],2);
        clims   = zeros([n_chan 2]);    %separate clim per chan
        for ch = 1:size(trials,1)
            clims(ch,1) = prctile(reshape(trials(ch,:,:),[1 size(trials,2)*size(trials,3)]),clim_perc(1));
            clims(ch,2) = prctile(reshape(trials(ch,:,:),[1 size(trials,2)*size(trials,3)]),clim_perc(2));
        end
        
        % Plot Condition ERPs
        out_file = [data_id '_' analysis_id '_trials_' event '_' cond_id grand_id...
            '_Bob.ep' epoch_id '_' trl_id '_bsln' bsln_id env_id smooth_id y_scale_id];
        if n_chan > 2
            fig_height = 1;
        else
            fig_height = n_chan/3;
        end
        figure('Name',out_file,'units','normalized',...
            'outerposition',[0 0 1 fig_height],'Visible',vis_fig);
        for ch = 1:n_chan
            % Plot Event Locked Average
            subplot(n_chan,length(cond_lab)+1,ch*(length(cond_lab)+1)-length(cond_lab)); hold on;
            ebars = {};
            main_lines = [];
            for cond_ix = 1:length(cond_lab)
                if n_chan==1
                    ebar{cond_ix} = shadedErrorBar([],...
                        trial_mean{cond_ix}(win_on:win_off),...
                        trial_sem{cond_ix}(win_on:win_off),...
                        {'Color',[cond_colors{cond_ix}],...
                        'LineStyle',cond_style{cond_ix}},sem_alpha);
                else
                    ebar{cond_ix} = shadedErrorBar([],...
                        trial_mean{cond_ix}(ch,win_on:win_off),...
                        trial_sem{cond_ix}(ch,win_on:win_off),...
                        {'Color',[cond_colors{cond_ix}],...
                        'LineStyle',cond_style{cond_ix}},sem_alpha);
                end
                cond_lab_legend{cond_ix} = [cond_lab{cond_ix} '-' num2str(sum(cond_idx{cond_ix}))];
                main_lines = [main_lines ebar{cond_ix}.mainLine];
            end
            ylim(ylims(ch,:));
            event_line = line([plot_lim(1) plot_lim(1)],ylim);%, 'k--');
            % Axes and Labels
            ax = gca;
            ax.XLim = [0,win_off];
            ax.XTick = 0:x_step:win_off;
            ax.XTickLabel = x_lab;
            main_lines = [main_lines event_line];
%             if avg_channels == 1
%                 title(strcat('Mean',avg_ch_lab, ': ', event,'-locked ',analysis_id));
%             else
            title(strcat(header_ecog.channel_labels(ch), ': ', event,'-locked ',analysis_id));
%             end
            if ch==1
                legend(main_lines,cond_lab_legend{:},event,'Location','southeast');%,'inc-con'
            end
            
            % Plot Single Trials Per Condition
            for cond_ix = 1:length(cond_lab)
                subplot(n_chan,length(cond_lab)+1,...
                    ch*(length(cond_lab)+1)-length(cond_lab)+cond_ix); hold on;
                
                if n_chan==1
                    imagesc(squeeze(trials(RT_sort_idx(cond_idx{cond_ix}),win_on:win_off)));
                else
                    imagesc(squeeze(trials(ch,RT_sort_idx(cond_idx{cond_ix}),win_on:win_off)));
                end
                set(gca,'YDir','normal');
                %!!! check y order by making a line all 100s
                if strcmp(event,'stim')
                    scat = scatter(RTs_sort_cond{cond_ix}+plot_lim(1),1:sum(cond_idx{cond_ix}),...
                        'MarkerFaceColor',[cond_colors{cond_ix}],...
                        'MarkerEdgeColor','k');
                end
                ylim([1 sum(cond_idx{cond_ix})]);
                event_line = line([plot_lim(1) plot_lim(1)],ylim);%, 'k--');
                % Axes and Labels
                cbar = colorbar;
                caxis(clims(ch,:));
                ax = gca;
                ax.XLim = [0,win_off];
                ax.XTick = 0:x_step:win_off;
                ax.XTickLabel = x_lab;
                
%                 if avg_channels == 1
%                     title(strcat('Mean',avg_ch_lab, ': Single Trial ',...
%                         event,'-locked ',analysis_id));
%                 else
                title(strcat(header_ecog.channel_labels(ch), ' ', cond_lab{cond_ix}, ': Single Trial ',...
                    event,'-locked ',analysis_id));
%                 end
            end
        end
        
        % Save plots
        out_filename = [fig_dir out_file '.' fig_type];
        if save_fig ==1
            fprintf('Saving %s\n',out_filename);
            eval(['export_fig ' out_filename]);
        end
        % saveas(gcf,out_filename);
        
        clear data_ecog header_ecog trial_info SBJ data_id good_epochs ok_epochs no_RT_epochs word_onset...
            resp_onset RTs max_RT avg_ch_ix avg_ch_lab data_filt data_hg...
            trials cond_idx trial_mean trial_sem out_file...
            minmax y_bound
    end
end
