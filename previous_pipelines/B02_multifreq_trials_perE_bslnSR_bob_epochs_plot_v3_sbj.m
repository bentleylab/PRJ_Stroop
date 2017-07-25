%% Stroop analysis - High Gamma Power on Single Trials
% v3 just does the filtering and saves out (not even cutting trials)
clear all; %close all
addpath(genpath('/home/knight/hoycw/PRJ_Stroop/scripts/'));
addpath(genpath('/home/knight/hoycw/Apps/export_fig-master/'));
addpath('/home/knight/hoycw/Apps/fieldtrip/');
ft_defaults
% addpath(genpath('/home/knight/hoycw/Apps/EEGlab/eeglab12_0_2_6b/functions/sigprocfunc/'));

% Analysis Parameters
data_info = {{'IR21','ft_preproc',{'RC*','ROF*','LIN*'}},...
             {'IR31','ft_preproc',{'RAC*','LIN*','LAC*'}},...
             {'IR35','ft_preproc',{'LAC*','RIN*','LPC*','LIN*'}},...
             {'IR39','ft_preproc',{'RAC*','ROF*','RIN*','LOF*','LIN*','RAM7-RAM8'}}};
%              {'IR32','IH_CA',{''}},...
cond_name = 'CI';              % conditions to average
analysis  = 'ERP';%,'HG'};%'theta','beta'
plot_rois = {'INS','MCC','MFG'};
events    = {'stim','resp'};
%   CI   = only con and inc trial types for all blocks
%   pcon = proportion congruency (con and inc for both mcon and minc)
%   !!!conseq = congruency sequence effects (cC, cI, iC, iI)
epoch_lim = [-200 2000];           % ID of data used for trial rejection (these epochs are NOT used)
buff_lim  = [1000, 1000];            % buffers are for cutting time series and then plotting
sig_lim   = [500, 500];

% Analysis parameters
HG_type      = 'wideband';                     % only for HG: 'wideband', 'multiband'
event        = 'SvsR';                          % 'stim'/'resp': event to lock trials
bsln_type    = 'demean';%'zscore'};    % 'zscore', 'demean', 'none'
bsln_event   = 's';                    % 's'/'r': event to lock baselining
bsln_lim     = [250, -50];             % ms before and after event
smooth_it    = 1;                      % smooth the result before averaging (0/1)
smooth_freq  = 10;
trial_type   = 'datapad';               % 'nanpad', 'datapad'
sig_win_len  = 50;                      % Length of significance testing window in ms
sig_win_step = 25;                      % Step size to slide teh sig window forward
sig_clust_len= 100;                     % length in ms over which consecutive windows must be significant
sig_nboots   = 1000;                    % number of iterations for permutation testing

% Plotting parameters
save_fig      = 1;
vis_fig       = 'on';                % 'on' or 'off' to determine if plot is shown
sem_alpha     = 0.5;                  % transparency of sem shading (0:1)
clim_perc     = [5 95];              % colormap percentile limits
plot_lim      = [500 500];            % ms to plot before the event and **after RT**
stim_plot_max = 1500;
x_step        = 250;                  % step of x tick marks
event_ln_width= 2;
fig_type      = 'png';

%% Process parameters
% Condition Parameters
[cond_lab, cond_colors, cond_style] = fn_condition_label_styles(cond_name);

% % Analysis Parameters
% all_analysis_id = '';
% all_bsln_type = '';
% for an_ix = 1:length(analysis)
%     [analysis_id{an_ix}, env_it{an_ix}, filt_lim{an_ix}] =...
%         fn_B02_analysis_params(analysis{an_ix},HG_type);
%     all_analysis_id = strcat(all_analysis_id, analysis_id{an_ix}(1));
%     all_bsln_type = strcat(all_bsln_type, bsln_type{an_ix}(1));
% end

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
bsln_id  = ['_' bsln_type '.' bsln_event num2str(bsln_lim(1)) '.' num2str(bsln_lim(2))];
if strcmp(trial_type,'nanpad')
    trl_id = 'nan';
elseif strcmp(trial_type,'datapad')
    trl_id = 'data';
else
    error(strcat('Unknown trial averaging type: ',trial_type));
end
% if length(bsln_type)~=length(analysis)
%     error('Mismatched number of bsln_type and analyses');
% end
sig_id = ['sig' num2str(sig_lim(1)) '.' num2str(sig_lim(2))];

% Plotting (Y axis)
if (plot_lim(1)>buff_lim(1)) || (plot_lim(2)>buff_lim(2))
    error('plot_lim exceed buff_lim')
end
plot_lim_id = strcat('_xlim',num2str(plot_lim(1)),'.',num2str(plot_lim(2)));
y_scale_id = strcat('_c',num2str(clim_perc(1)),'.',num2str(clim_perc(2)));

roi_list = {};

%% Analysis Loop
for data_ix = 1:length(data_info)
    SBJ = data_info{data_ix}{1};
    data_id = strcat(SBJ,'_',data_info{data_ix}{2});
    ch_id = data_info{data_ix}{3}; 
    SBJ_dir = strcat('/home/knight/hoycw/PRJ_Stroop/data/',SBJ,'/');
    fprintf('=====================================================\n');
    fprintf('===================== %s ========================\n',SBJ);
    
    %% Load data
    % header_ecog should be the same for all of them
    load(strcat(SBJ_dir,'04_proc/',data_id,'_',analysis,smooth_id),'header_ecog');
    if strcmp(SBJ,'IR32')
        load(strcat(SBJ_dir,'03_events/',data_id,'_trial_info.mat'));
    else
        load(strcat(SBJ_dir,'03_events/',SBJ,'_trial_info_full.mat'));
    end
    
    %% ROI info
    roi_filename = strcat(SBJ_dir,'04_proc/',SBJ,'_sfn_elec_ROIs.csv');
    einfo{data_ix} = textread(roi_filename,'%s','delimiter',',');
    einfo{data_ix} = reshape(einfo{data_ix},[3, length(header_ecog.channel_labels)]); % elec_name,ROI,tissue(GG,GW,GW,WW,BS)
    gm_elec_idx = logical(zeros([1 size(einfo{data_ix},2)]));
    for e_ix = 1:size(einfo{data_ix},2)
        if ~isempty(strfind(einfo{data_ix}{3,e_ix},'G'))
            gm_elec_idx(e_ix) = true;
        end
    end
    roi_list{data_ix} = unique(einfo{data_ix}(2,gm_elec_idx));
    
    %% Process RTs
    load(strcat(SBJ_dir,'04_proc/',data_id,'_good_epochs.mat')); %loads ok_epochs
    fprintf('Num ok_epochs : %i\n',length(ok_epochs));
    
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
    
    %% Load Data
    load(strcat(SBJ_dir,'04_proc/',data_id,'_',analysis,smooth_id),'data_proc');
    % Trial Parameters
    n_chan = length(header_ecog.channel_labels);
    trial_mean{data_ix} = {};
    minmax{data_ix}     = {};
    for evnt_ix = 1:2
        event = events{evnt_ix};
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
        %     %         !!! assumes 1kHz sample rate, also only 2 conditions!
        %     sig_wins = [buff_lim(1)-sig_lim(1):sig_win_step:trial_len-buff_lim(2)+sig_lim(2)-sig_win_len]+1;
        %     sig_wins = vertcat(sig_wins,sig_wins+sig_win_len);
        
        % Pre-Allocate trials
        trials     = NaN([n_chan length(good_epochs) trial_len]);
        trial_mean{data_ix}{evnt_ix} = NaN([n_chan length(cond_lab) trial_len]);
        %     trial_sem{data_ix}  = NaN([n_chan length(cond_lab) trial_len]);
        minmax{data_ix}{evnt_ix}     = NaN([n_chan length(cond_lab) 2]);   %separate ylim per chan
        %     sig_pvals{data_ix}  = NaN([n_chan size(sig_wins,2)]);
        %     sig_times{data_ix}  = {};
        %     for an_ix = 1:length(analysis)
        
        %% Format trials
        for ch_ix = 1:size(data_proc,1)
            % Cut trials and baseline all at once
            trials(ch_ix,:,:) = fn_epoch_bsln_combo(data_proc(ch_ix,:), lock_events,...
                resp_onset, buff_lim, trial_type, bsln_events, bsln_lim, bsln_type);
        end
        
        %% ERPs by condition
        for cond_ix = 1:length(cond_lab)
            % Average trials matching condition index
            trial_mean{data_ix}{evnt_ix}(:,cond_ix,:) = squeeze(nanmean(trials(:,cond_idx(cond_ix,:),:),2));
            %             % Standarad Error of Mean
            %             trial_sem{data_ix}(:,cond_ix,:) = squeeze(nanstd(trials{data_ix}(:,cond_idx(cond_ix,:),:),...
            %                 0,3))./sqrt(sum(cond_idx(cond_ix,:)));
            % Plot Limits per electrode for mean
            minmax{data_ix}{evnt_ix}(:,cond_ix,1) = min(trial_mean{data_ix}{evnt_ix}(:,cond_ix,:),[],3);
            minmax{data_ix}{evnt_ix}(:,cond_ix,2) = max(trial_mean{data_ix}{evnt_ix}(:,cond_ix,:),[],3);
        end
    end
    clear data_proc header_ecog trial_info SBJ data_id good_epochs ok_epochs no_RT_epochs word_onset...
        resp_onset RTs max_RT cond_idx ch_ix trials
end
    
%% Plot ERPs
for data_ix = 1:length(data_info)
    SBJ = data_info{data_ix}{1};
    
    fig_dir  = strcat('/home/knight/hoycw/PRJ_Stroop/results/multifreq/',SBJ,'/',cond_name,'/',analysis,'/');
    if ~exist(fig_dir,'dir')
        mkdir(fig_dir);
    end
    
    % Plot Condition ERPs
    fig_file = [SBJ '_' analysis '_' strjoin(plot_rois,'.') '_'...
        cond_name '_' event '_Bob.ep' epoch_id '_' trl_id...
        '_bsln' bsln_id smooth_id plot_lim_id y_scale_id];
    fig_height = 1;%length(analysis)/3;
    figure('Name',fig_file,'units','normalized',...
        'outerposition',[0 0 1 fig_height],'Visible',vis_fig);
    
    for roi_ix = 1:length(plot_rois)
        %         ylims = zeros([1 2]);
        %         clims = zeros([1 2]);
        %         ylims(1) = min(minmax(ch_ix,:,1),[],3);
        %         ylims(2) = max(minmax(ch_ix,:,2),[],3);
        %         clims(1) = prctile(reshape(trials(ch_ix,:,:),...
        %             [1 size(trials,3)*size(trials,4)]),clim_perc(1));
        %         clims(2) = prctile(reshape(trials(an_ix,ch_ix,:,:),...
        %             [1 size(trials,3)*size(trials,4)]),clim_perc(2));
        for evnt_ix = 1:length(events)
            event = events{evnt_ix};
            
            % Plotting parameters
            win_on  = buff_lim(1)-plot_lim(1)+1;
            if strcmp(event,'stim')
                win_off = buff_lim(1)+stim_plot_max+1;         % Capture even longest trials
            else
                win_off = buff_lim(1)+plot_lim(2)+1;
            end
            x_lab   = -plot_lim(1):x_step:win_off-buff_lim(1);
            
            % Plot Event Locked Average
            subplot(length(plot_rois),length(events),...
                fn_RC2subplot_ix(length(plot_rois),length(events),roi_ix,evnt_ix));
            hold on;
            ch_cnt = 0;
            for ch_ix = 1:size(trial_mean{data_ix}{evnt_ix},1)
                % ROI selection
                if (strcmp(einfo{data_ix}{2,ch_ix},plot_rois{roi_ix}))% && (gm_elec_idx{data_ix}(ch_ix)==1)
                    ch_cnt = ch_cnt+1;
                    for cond_ix = 1:length(cond_lab)
                        plot(squeeze(trial_mean{data_ix}{evnt_ix}(ch_ix,cond_ix,win_on:win_off)),...
                            'Color',[cond_colors{cond_ix}]);%,'LineStyle',cond_style{cond_ix}},sem_alpha);
                    end
                end
            end
            %         ylim(ylims);
            event_line = line([plot_lim(1) plot_lim(1)],ylim,...
                'LineWidth',event_ln_width,'Color','k');
            %         % Significance marker
            %         scatter(sig_times{an_ix}{ch_ix}-plot_lim(1),zeros([1 length(sig_times{an_ix}{ch_ix})]),...
            %             'Marker','.','LineWidth',0.1,'MarkerEdgeColor','k');
            % Axes and Labels
            ax = gca;
            ax.XLim = [0,win_off-win_on];
            ax.XTick = 0:x_step:win_off;
            ax.XTickLabel = x_lab;
            %         main_lines = [main_lines event_line];
            title(strcat(plot_rois{roi_ix},' (', event,'-locked): ',num2str(ch_cnt),'ch'));
            %         if an_ix==1
            %             legend(main_lines,cond_lab_legend{:},event,'Location','southeast');%,'inc-con'
            %         end
        end
    end
    
    % Save plots
    out_filename = [fig_dir fig_file '.' fig_type];
    if save_fig ==1
        fprintf('Saving %s\n',out_filename);
        eval(['export_fig ' out_filename]);
    end
    % saveas(gcf,out_filename);
end




% % Plot Single Trials Per Condition
% for cond_ix = 1:length(cond_lab)
%     subplot(length(analysis),length(cond_lab)+1,...
%         an_ix*(length(cond_lab)+1)-length(cond_lab)+cond_ix); hold on;
%     
%     imagesc(squeeze(trials(an_ix,ch_ix,RTs_sort_idx(cond_idx(cond_ix,:)),win_on:win_off)));
%     set(gca,'YDir','normal');
%     if strcmp(event,'stim')
%         scat = scatter(RTs_sort_cond{cond_ix}+plot_lim(1),1:sum(cond_idx(cond_ix,:)),...
%             'MarkerFaceColor',[cond_colors{cond_ix}],...
%             'MarkerEdgeColor','k');
%     end
%     ylim([1 sum(cond_idx(cond_ix,:))]);
%     event_line = line([plot_lim(1) plot_lim(1)],ylim,...
%         'LineWidth',event_ln_width,'Color','k');
%     % Axes and Labels
%     cbar = colorbar;
%     caxis(clims);
%     ax = gca;
%     ax.XLim = [0,win_off-win_on];
%     ax.XTick = 0:x_step:win_off;
%     ax.XTickLabel = x_lab;
%     
%     title(strcat(header_ecog.channel_labels(ch_ix), ':', analysis_id{an_ix},...
%         ',',cond_lab{cond_ix},' (Single Trial, ',event,'-locked)'));
% end

%         %% Sliding Window Statistics       
%         sig_times{data_ix}{an_ix} = {};
%         for ch_ix = 1:n_chan
%             fprintf('=========================Significance Testing: %s %s %s=========================\n',...
%                 SBJ,header_ecog.channel_labels{ch_ix},analyses{an_ix})
%             sig_times{an_ix}{ch_ix} = [];
%             for sig_win_ix = 1:size(sig_wins,2)
%                 % Permute trial means within each window
%                 sig_win = sig_wins(1,sig_win_ix):sig_wins(2,sig_win_ix);
%                 [sig_pvals(an_ix,ch_ix,sig_win_ix),~] = fn_sig_permute_labels_mean_diff(...
%                     squeeze(mean(trials(an_ix,ch_ix,cond_idx(1,:),sig_win),4)),...
%                     squeeze(mean(trials(an_ix,ch_ix,cond_idx(2,:),sig_win),4)),sig_nboots);
%             end
%             % Adjust for multiple comparisons- must use BHFDR method
%             % because Storey's method falls apart with <1000 p values: 
%             %https://www.mathworks.com/matlabcentral/answers/144113-mafdr-interpreting-q-values-vs-bhfdr-adjusted-p-values
%             sig_qvals = mafdr(squeeze(sig_pvals(an_ix,ch_ix,:)),'BHFDR',true);
%             sig_ts = logical(sig_qvals<0.05);
%             % calculate number of consecutive windows to meet the cluster length
%             if sum(sig_ts)>0
%                 chunk_edges = find(diff(sig_ts));
%                 for chunk_ix = 0:length(chunk_edges)
%                     if chunk_ix==0
%                         chunk_start = 1;
%                     else
%                         chunk_start = chunk_edges(chunk_ix)+1;
%                     end
%                     if chunk_ix==length(chunk_edges)
%                         chunk_end = size(sig_wins,2);
%                     else
%                         chunk_end = chunk_edges(chunk_ix+1);
%                     end
%                     fprintf('chunk %i (sig=%i): %i-%i\n',chunk_ix,sig_ts(chunk_start),sig_wins(1,chunk_start),sig_wins(2,chunk_end));
%                     if (sig_ts(chunk_start)==1) && (sig_wins(2,chunk_end)-sig_wins(1,chunk_start)>=sig_clust_len)
%                         sig_times{an_ix}{ch_ix} = [sig_times{an_ix}{ch_ix} sig_wins(1,chunk_start):sig_wins(2,chunk_end)];
%                     end
%                 end
%             end
% %             figure;hold on;
% %             plot(squeeze(sig_pvals(an_ix,ch_ix,:)),'b');
% %             plot(squeeze(sig_qvals),'r');
% %             scatter(find(sig_ts==1),ones([1 length(find(sig_ts==1))])*0.5);
% %             title(strcat(SBJ,'_',header_ecog.channel_labels{ch_ix},'_',analyses{an_ix}));
%             clear sig_qvals sig_ts
%         end
%     end
%     sig_times_file = [SBJ_dir '04_proc/' data_id '_' all_analysis_id...
%         '_' event '_' cond_name '_Bob.ep' epoch_id '_' trl_id '_bsln' bsln_id smooth_id '_' sig_id '_times.mat'];
%     save(sig_times_file,'sig_times','sig_wins','sig_pvals','header_ecog');%'ch_ids',