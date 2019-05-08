function SBJ08d_HFA_ROL(SBJ, an_id, rol_id, plot_qa, plot_stack, fig_vis, fig_ftype)
% Find single trial onset latencies of HFA
%   Adapted from code by Eleonora Bartoli and Brett Foster
% INPUTS:
%   SBJ
%   proc_id
%   an_id
%   plot_qa [0/1] - plot a random subset of trials for quality assurance
%   plot_stack [0/1] - plot single trial power stacks with ROL and RT marked
%   fig_vis
%   fig_ftype
% OUTPUTS:
%   ____ [float mat] - matrix size(n_trl,2) of onsets per trial using
%       (regression, derivative) methods

%% Set up paths
if exist('/home/knight/hoycw/','dir');root_dir='/home/knight/hoycw/';ft_dir=[root_dir 'Apps/fieldtrip/'];
else root_dir='/Volumes/hoycw_clust/';ft_dir='/Users/colinhoy/Code/Apps/fieldtrip/';end
addpath([root_dir 'PRJ_Stroop/scripts/']);
addpath([root_dir 'PRJ_Stroop/scripts/utils/']);
addpath(ft_dir);
ft_defaults

%% Data Preparation
SBJ_vars_cmd = ['run ' root_dir 'PRJ_Stroop/scripts/SBJ_vars/' SBJ '_vars.m'];
eval(SBJ_vars_cmd);
an_vars_cmd = ['run ' root_dir 'PRJ_Stroop/scripts/an_vars/' an_id '_vars.m'];
eval(an_vars_cmd);
rol_vars_cmd = ['run ' root_dir 'PRJ_Stroop/scripts/an_vars/' rol_id '_vars.m'];
eval(rol_vars_cmd);

% Load Data
hfa_fname = strcat(SBJ_vars.dirs.proc,SBJ,'_ROI_',an_id,'.mat');
load(hfa_fname);
load(strcat(SBJ_vars.dirs.events,SBJ,'_trial_info_final.mat'));

% Sort RTs
[~,rt_sort_idx] = sort(trial_info.response_time);
rt_sort_colors = parula(numel(trial_info.response_time));

% Output prep
qa_main_dir = [root_dir 'PRJ_Stroop/results/HFA/' SBJ '/ROL/' an_id '/QA/'];
stack_dir = [root_dir 'PRJ_Stroop/results/HFA/' SBJ '/ROL/' an_id '/stack/'];
if ~isdir(stack_dir)
    mkdir(stack_dir);
end

%% Select time range for onsets
cfgs = [];
cfgs.latency = [rol_trl_lim_s(1) max(trial_info.response_time)+rol_trl_lim_s(2)];
hfa = ft_selectdata(cfgs, hfa);

% Remove post

%% Compute parameters for ROLs
% Power threshold per channel over entire trial
thresh = quantile(reshape(hfa.powspctrm,[numel(hfa.label) size(hfa.powspctrm,1)*size(hfa.powspctrm,4)]),quant_thresh,2);

% Time threshold
%!!! fix this naming for the different windows!
sample_rate = (numel(hfa.time)-1)/(hfa.time(end)-hfa.time(1));
min_actv_dp = round(min_actv_s*sample_rate);

% % Create pseudo stimulus distribution
% %   correlation with all zeros has no variance and is undefined (divide by 0 = NaN)
% %   therefore, correlation with "stimulus" will be gaussian around 0 but with small variance
% pseudo_stim_times = normrnd(0,stim_var,size(trial_info.response_time));

%% Response Onset Latency Analyses
% response onset latencies (in sec) per channel and trial for [derivative, linear regression] methods 
rol_times  = nan([numel(hfa.label) size(hfa.powspctrm,1) numel(rol_lab)]);   % ROL in sec [D/L]
actv_len   = nan([numel(hfa.label) size(hfa.powspctrm,1)]);     % # dp above threshold in ROL window
actv_amp   = nan([numel(hfa.label) size(hfa.powspctrm,1)]);     % mean power above threshold in ROL window
rol_win_sz = nan([numel(hfa.label) size(hfa.powspctrm,1)]);     % ROL window size in dp
evnt_lock  = zeros([numel(hfa.label) numel(rol_lab)]);          % event that best explains channel's ROL
plt_trls   = randi(size(hfa.powspctrm,1),round(trl_plt_perc*size(hfa.powspctrm,1)),1);
for ch_ix = 1:10:numel(hfa.label)
    fprintf('-------------------------- %s (%i/%i) --------------------------\n',hfa.label{ch_ix},ch_ix,numel(hfa.label));
    %% Smooth data
    ch_data_orig = squeeze(hfa.powspctrm(:,ch_ix,:,:));
    ch_data = sgolayfilt(ch_data_orig',sgfilt_ord,sgfilt_win)';     % Filters columns, so flip (2x)
    
    %% Compute ROLs
    above_lens = [];
    above_amps = [];
    above_kept = [];
    above_mat  = zeros(size(ch_data));
    n_ep_above = zeros([size(hfa.powspctrm,1) 1]);
    for trl_ix = 1:size(hfa.powspctrm,1)
        trl_data = ch_data(trl_ix,:);
        % Plot this trial if in random subset
        if any(plt_trls==trl_ix) && plot_qa
            qa_dir = [qa_main_dir hfa.label{ch_ix} '/'];
            if ~isdir(qa_dir)
                mkdir(qa_dir);
            end
            fig_name = [SBJ '_' hfa.label{ch_ix} '_t' num2str(trl_ix)];
            figure('Name',fig_name,'Visible',fig_vis);%'off');    % always 'off' because too many plots
            subplot(3,1,1); hold on;
            title('Full Trial');
            plot(hfa.time,ch_data_orig(trl_ix,:),'k');
            plot(hfa.time,trl_data,'r');
            line([hfa.time(1) hfa.time(end)],[thresh(ch_ix) thresh(ch_ix)],...
                 'Color','b','LineStyle','--');
            ax = gca;
            ax.XLabel.String = 'Time (s)';
            ax.YLabel.String = 'HFA z-score';
            ax.XLim = [hfa.time(1) hfa.time(end)];
        end
        
        % Find times above threshold
%         above_th  = find(trl_data>thresh(trl_ix));
%!!! convert to (ep,2) column vectors instead of 2 variables!
        beg_above = find(diff([0 trl_data>thresh(ch_ix) 0])==1);
        end_above = find(diff([0 trl_data>thresh(ch_ix) 0])==-1)-1;
        len_above = end_above-beg_above;
        
        % Remove epochs shorter than win_len
        beg_above = beg_above(len_above>min_actv_dp);
        end_above = end_above(len_above>min_actv_dp);
        
        %% Compute ROLs for windows above threshold for at least win_len
        n_ep_above(trl_ix) = numel(beg_above);
        if any(beg_above)
            % Take epoch with highest overall power
            if numel(beg_above)>1
                ep_mean = zeros(size(beg_above));
                ep_lens = zeros(size(beg_above));
                ep_kept = zeros(size(beg_above));
                for ep_ix = 1:numel(beg_above)
                    ep_mean(ep_ix) = mean(trl_data(beg_above(ep_ix):end_above(ep_ix)));
                    ep_lens(ep_ix) = end_above(ep_ix)-beg_above(ep_ix)+1;
                    
                    above_mat(trl_ix,beg_above(ep_ix):end_above(ep_ix)) = ep_mean(ep_ix);
                end
                
                % Select biggest epoch
                [~, max_ix] = max(ep_mean);
                above_idx = beg_above(max_ix):end_above(max_ix);
                ep_kept(max_ix) = 1;
            else
                above_idx = beg_above:end_above;
                
                ep_mean = mean(trl_data(above_idx));
                ep_lens = numel(above_idx);
                ep_kept = 1;
                above_mat(trl_ix,above_idx) = ep_mean;
            end
            
            % Keep stats on all epochs above thresh
            above_amps = [above_amps ep_mean];
            above_lens = [above_lens ep_lens];
            above_kept = [above_kept ep_kept];
            
            % Define start and end of onset search epoch
            rol_lim = round([above_idx(1)+rol_lim_s(1)*sample_rate ...
                             above_idx(1)+rol_lim_s(2)*sample_rate]);
            
            % Check search epoch is within trial limits
            if rol_lim(1)<=0
                rol_lim(1) = 1;
            end
            if rol_lim(2)>numel(trl_data)
                rol_lim(2) = numel(trl_data);
            end
            
            % Grab search epoch indices and data
            search_ep_ix = rol_lim(1):rol_lim(2);
            search_ep = trl_data(search_ep_ix);
            
            % Record ROL window properties
            rol_win_sz(ch_ix,trl_ix) = rol_lim(2)-rol_lim(1);
            actv_len(ch_ix,trl_ix)   = numel(above_idx);
            actv_amp(ch_ix,trl_ix)   = mean(trl_data(above_idx));
%             if rol_win_sz(ch_ix,trl_ix) < round((rol_lim_s(2)-rol_lim_s(1))*sample_rate)
%                 fprintf('\tWarning: trial %i rol window is %i dp long (not %i)\n',...
%                         trl_ix,rol_win_sz(ch_ix,trl_ix),round((rol_lim_s(2)-rol_lim_s(1)))*sample_rate);
%             end
            
            % Plot search epoch for ROL
            if any(plt_trls==trl_ix) && plot_qa
                ylims = get(gca,'ylim');
                patch([hfa.time(rol_lim(1)) hfa.time(rol_lim(1)) hfa.time(rol_lim(2)) hfa.time(rol_lim(2))],...
                    [ylims(1) ylims(2) ylims(2) ylims(1)],...
                    [0.5 0.5 0.5], 'FaceAlpha', 0.4);
            end
            
            %% DERIVATIVE METHOD: Median of steepest slope and slope onset
            % Find steepest slope (maximum of derivative)
            deriv = diff(search_ep,1);
            [~, max_deriv_ix] = max(deriv);
            
            % Find start of upward slope (preceeding local min/zero crossing)
            local_mins = find(diff(sign(deriv),1));
            early_local_mins = local_mins(local_mins<max_deriv_ix);
            if ~isempty(early_local_mins)
                slope_beg_ix = early_local_mins(end);
            else
                % If no local min in rol_lim, take first point of window
                slope_beg_ix = 1;
            end
            
            % Define ROL as median of [slope_start max(slope)]
            rol_deriv_ix = round(median([slope_beg_ix max_deriv_ix]));
            rol_times(ch_ix,trl_ix,1) = hfa.time(search_ep_ix(rol_deriv_ix));
            
            %% Plot ROL checks
            if any(plt_trls==trl_ix) && plot_qa
                subplot(3,1,2); hold on;
                title('Derivative Method');
                
                % Plot search epoch
                plot(hfa.time(search_ep_ix),search_ep,'k','LineWidth',2);
                
%                 % Plot first order derivative (scale derivative to make it easily visible)
%                 d = plot(hfa.time(search_ep_ix(1:end-1)),diff(search_ep,1)*deriv_scale,'g');
                
                % Plot local min and steepest slope
                s = line([hfa.time(search_ep_ix(slope_beg_ix)) hfa.time(search_ep_ix(slope_beg_ix))],ylim,...
                     'LineStyle','--','Color','b');
                m = line([hfa.time(search_ep_ix(max_deriv_ix)) hfa.time(search_ep_ix(max_deriv_ix))],ylim,...
                     'LineStyle','--','Color','b');
                % Plot estimated ROL
                r = line([hfa.time(search_ep_ix(rol_deriv_ix)) hfa.time(search_ep_ix(rol_deriv_ix))],...
                         get(gca,'ylim'), 'Color','r');
                legend([s r m],{'local min','ROL','steepest'},'Location','best');
                
                % Set up plot for linear method
                subplot(3,1,3); hold on;
                title('Linear Fit Method');
                plot(hfa.time(search_ep_ix),search_ep,'k','LineWidth',2);
            end
            
            %% LINEAR REGRESSION METHOD:
            % Set up sliding window parameters
            win_len  = round(reg_win_len_s*sample_rate);
            win_step = round(reg_win_stp_s*sample_rate);
            win_lim  = zeros([ceil((numel(search_ep)-win_len-1)/win_step) 2]);           
%             win_lim = fn_sliding_window_lim(search_ep,round(reg_win_len_s*sample_rate),...
%                                                 round(reg_win_stp_s*sample_rate));
                        
            % Linear regression in each window
            lm_fit = zeros(size(win_lim,1),5);
            for win_ix = 1:size(win_lim,1)
                % Set window parameters
                if win_ix==1
                    win_lim(win_ix,:) = [1 win_len];
                else
                    win_lim(win_ix,:) = [win_lim(win_ix-1,1)+win_step win_lim(win_ix-1,2)+win_step];
                end
                
                % coeff(1) = slope; coeff(2) = intercept
                lm_coeff   = polyfit(hfa.time(search_ep_ix(win_lim(win_ix,1):win_lim(win_ix,2))),...
                                  search_ep(win_lim(win_ix,1):win_lim(win_ix,2)),1);
                model   = lm_coeff(1)*hfa.time(search_ep_ix(win_lim(win_ix,1):win_lim(win_ix,2)))+lm_coeff(2);
                abs_err = abs(search_ep(win_lim(win_ix,1):win_lim(win_ix,2))-model);
                
                % Outputs: slope, intercept, max error-residual, win onset, win offset
                lm_fit(win_ix,:) = [lm_coeff max(abs_err) hfa.time(search_ep_ix(win_lim(win_ix,1))) ...
                                    hfa.time(search_ep_ix(win_lim(win_ix,2)))];
                % Plot the fits
                if any(plt_trls==trl_ix) && plot_qa
                    plot(hfa.time(search_ep_ix(win_lim(win_ix,1):win_lim(win_ix,2))),model,'Color','b');
                end
            end
            
            % Sort the slopes from biggest to smallest
            [sorted_slopes, sorted_slope_idx] = sort(lm_fit(:,1), 'descend');
            
            % Remove negative slopes
            sorted_slope_idx = sorted_slope_idx(sorted_slopes>0);
            sorted_slopes = sorted_slopes(sorted_slopes>0);
            
            % Take biggest slopes
            if numel(sorted_slope_idx)>=n_big_slopes
                big_slope_ix = sorted_slope_idx(1:n_big_slopes);
            else
                big_slope_ix = sorted_slope_idx;
            end
            
            % Sort residual errors of big slopes from smallest to biggest
            [~, sorted_err_ix] = sort(lm_fit(big_slope_ix,3), 'ascend');
            %get smallest error data, if any
            
            % Take the big slope with smallest error
            if numel(sorted_slope_idx)>0
                lm_result = lm_fit(big_slope_ix(sorted_err_ix(1)),:);
                rol_times(ch_ix,trl_ix,2) = lm_result(4);
                
                % Plot lm result (onset of best window)
                if any(plt_trls==trl_ix) && plot_qa
                    line([lm_result(4) lm_result(4)],get(gca,'ylim'),'Color','r');
                    
                    % Plot both ROLs on original time series
                    subplot(3,1,1);
                    line([rol_times(ch_ix,trl_ix,1) rol_times(ch_ix,trl_ix,1)],...
                         get(gca,'ylim'),'Color','r');
                    line([rol_times(ch_ix,trl_ix,2) rol_times(ch_ix,trl_ix,2)],...
                         get(gca,'ylim'),'Color','r');
                    
                    % Save QA plot
                    saveas(gcf,[qa_dir fig_name '.' fig_ftype]);
                    if strcmp(fig_vis,'off')
                        close(gcf);
                    end
                end
            end
        end
    end
    
    %% Plot overall QA for this channel
    if plot_qa
        fig_name = [SBJ '_' hfa.label{ch_ix} '_QA_summary'];
        figure('Name',fig_name,'Visible',fig_vis,'Units','normalized','OuterPosition',[0 0 0.4 0.4]);
        % ROL window sizes
        subplot(1,3,1); hold on;
        small_win_cnt = sum(rol_win_sz(ch_ix,:)/sample_rate < rol_lim_s(2)-rol_lim_s(1)-0.001);
        title(['Full = ' num2str(rol_lim_s(2)-rol_lim_s(1)) '; # smaller win = ' ...
              num2str(small_win_cnt) ' (' num2str(100*small_win_cnt/sum(~isnan(rol_win_sz(ch_ix,:))),'%.01f') '%)']);
        histogram(rol_win_sz(ch_ix,:)/sample_rate, n_hist_bins);
        xlabel('ROL window size (ms)');
        
        % ROL activation lengths vs. amplitudes
        subplot(1,3,2); hold on;
        title(['Onset Epochs; ' num2str(mean(n_ep_above)) ' ep/trl']);
        % Plot Epochs not analyzed
        above_len_s = above_lens/sample_rate;
        toss_scat = scatter(above_len_s(~above_kept), above_amps(~above_kept), 'd', 'MarkerEdgeColor',[0.5 0.5 0.5]);
        % Plot Epochs analyzed
        kept_len_s = actv_len(ch_ix,:)/sample_rate;
        kept_scat = scatter(kept_len_s, actv_amp(ch_ix,:),'o','MarkerEdgeColor','b');
        toss_mean = scatter(nanmean(above_len_s(~above_kept)), nanmean(above_amps(~above_kept)), ...
                            toss_scat.SizeData*2, 'd', 'MarkerEdgeColor',[0.5 0.5 0.5], 'MarkerFaceColor', 'k');
        kept_mean = scatter(nanmean(kept_len_s),nanmean(actv_amp(ch_ix,:)),...
                            kept_scat.SizeData*2,'o','MarkerEdgeColor','b','MarkerFaceColor','r');
        tmp_xlim = get(gca,'XLim');
        thresh_line = line(xlim,[thresh(ch_ix) thresh(ch_ix)],'Color','r');
        set(gca,'XLim',tmp_xlim);
        xlabel('Time above threshold (ms)');
        ylabel('Mean epoch amplitude');
        legend([toss_mean kept_mean thresh_line],{['Toss (' num2str(sum(~above_kept)) ')'],...
               ['Kept (' num2str(sum(above_kept)) ')'],'Thresh'},'Location','best');
        
        % ROL times vs. amplitudes
        subplot(1,3,3); hold on;
        title('ROL vs. Amplitude');
        % Plot Epochs analyzed
        rol_scat = scatter(rol_times(ch_ix,:,1), actv_amp(ch_ix,:),[],rt_sort_colors,'o');
        rol_mean = scatter(nanmean(rol_times(ch_ix,:,1)), nanmean(actv_amp(ch_ix,:)), ...
                            rol_scat.SizeData*2, 'o', 'MarkerEdgeColor','b','MarkerFaceColor','r');
%         tmp_xlim = get(gca,'XLim');
%         thresh_line = line(xlim,[thresh(ch_ix) thresh(ch_ix)],'Color','r');
%         set(gca,'XLim',tmp_xlim);
        xlabel('ROL (s)');
        ylabel('Mean epoch amplitude');
%         cbar = colorbar;
%         cbar.Label.String = 'RT';
        legend(rol_mean,'mean','Location','best');
%         legend([r_mean kept_mean thresh_line],{['Toss (' num2str(sum(~above_kept)) ')'],...
%                ['Kept (' num2str(sum(above_kept)) ')'],'Thresh'},'Location','best');
        
        % Save figure
        saveas(gcf,[qa_dir fig_name '.' fig_ftype]);
        if strcmp(fig_vis,'off')
            close(gcf);
        end
    end
    
    %% Stimulus vs. Response Modeling of ROLs
    % Flag outlier ROLs
    
    % Regression based method:
    outlier_idx = false([size(rol_times,2) numel(rol_lab)]);
    stats       = cell([numel(rol_lab) numel(evnt_lab)]);
    err_var     = nan([numel(rol_lab) numel(evnt_lab)]);
    err         = nan([numel(rol_lab) numel(evnt_lab)]);
    for rol_ix = 1:numel(rol_lab)
        outlier_idx(:,rol_ix) = abs(rol_times(ch_ix,:,rol_ix)-nanmean(rol_times(ch_ix,:,rol_ix))) > ...
            rol_outlier_std*nanstd(squeeze(rol_times(ch_ix,:,rol_ix)));
        model = ones(size(trial_info.trial_n(~outlier_idx(:,rol_ix))));
        for evnt_ix = 1:numel(evnt_lab)
            if evnt_ix==1
                % Stimulus model
                %   Flip to response locked and try to fit the stimulus times
%                 model = [-trial_info.response_time(rt_sort_idx(~outlier_idx(:,rol_ix))) ones(size(trial_info.trial_n(~outlier_idx(:,rol_ix))))];
                fit_data = squeeze(rol_times(ch_ix,~outlier_idx(:,rol_ix),rol_ix))'-trial_info.response_time(rt_sort_idx(~outlier_idx(:,rol_ix)));
            else
                % Response model
%                 model = [trial_info.response_time(rt_sort_idx(~outlier_idx(:,rol_ix))) ones(size(trial_info.trial_n(~outlier_idx(:,rol_ix))))];
                fit_data = squeeze(rol_times(ch_ix,~outlier_idx(:,rol_ix),rol_ix))';
            end
            % stats: (R^2, Fval, pval, error varaince)
            [~, ~, resid, ~, stats{rol_ix,evnt_ix}] = regress(...
                fit_data, model);
            err(rol_ix,evnt_ix) = nansum(abs(resid));
            err_var(rol_ix,evnt_ix) = stats{rol_ix,evnt_ix}(4);
        end
        
        % Assign event label (index)
        [~,evnt_lock(ch_ix,rol_ix)] = min(err(rol_ix,:));
    end
    if any(diff(evnt_lock(ch_ix,:)))
        warning(['\tWARNING!!! ' hfa.label{ch_ix} ' has different event assignments across ROL methods!']);
    end
    
    % Plotting data to get a feeling...
%     for rol_ix = 1:2
%         outlier_idx(:,rol_ix) = abs(rol_times(ch_ix,:,rol_ix)-nanmean(rol_times(ch_ix,:,rol_ix))) > ...
%             rol_outlier_std*nanstd(squeeze(rol_times(ch_ix,:,rol_ix)));
%         rol_var(rol_ix,1) = nanstd(squeeze(rol_times(ch_ix,~outlier_idx(:,rol_ix),rol_ix)));
%         rol_var(rol_ix,2) = nanstd(squeeze(rol_times(ch_ix,~outlier_idx(:,rol_ix),rol_ix))'-trial_info.response_time(~outlier_idx(:,rol_ix)));
%         
%         [~,sr_ix] = min(rol_var(rol_ix,:));
%         fprintf('%s: %s var([S, R]) = [%f vs. %f]\n',rol_lab{rol_ix}, evnt_lab{sr_ix}, rol_var(rol_ix,1), rol_var(rol_ix,2));
%     end
%     
%     % Plot variance of ROLs vs. RT-ROL
%     figure;
%     for rol_ix = 1:2
%         subplot(1,2,rol_ix); hold on;
%         histogram(rol_times(ch_ix,~outlier_idx(:,rol_ix),rol_ix),n_hist_bins,'FaceColor','b');
%         histogram(rol_times(ch_ix,~outlier_idx(:,rol_ix),rol_ix)'-trial_info.response_time(~outlier_idx(:,rol_ix)),n_hist_bins,'FaceColor','r');
%     end
%     
%     figure;
%     scatter(rol_times(ch_ix,~outlier_idx(rt_sort_idx,rol_ix),rol_ix), ...
%             rol_times(ch_ix,~outlier_idx(rt_sort_idx,rol_ix),rol_ix)'-trial_info.response_time(~outlier_idx(:,rol_ix)),...
%             [],rt_sort_colors(~outlier_idx(rt_sort_idx,rol_ix)));
%     xlabel('ROL time (s)');
%     ylabel('ROL - RT');
    
    %% Plot single trial stacks with ROLs and RTs
    if plot_stack
        fig_name = [SBJ '_' evnt_lab{evnt_lock(ch_ix,1)} '_' hfa.label{ch_ix} '_ROL_stack'];
        figure('Name',fig_name,'units','normalized',...
               'outerposition',[0 0 0.5 0.6],'Visible',fig_vis);
        
        % PLOT ROLS
        subplot(1,2,1); hold on;
        % Get color limits
        clims = [prctile(ch_data(:),clim_perc(1)) prctile(ch_data(:),clim_perc(2))];
        
        % Plot single trial stack
        imagesc(hfa.time,1:numel(trial_info.trial_n),ch_data(rt_sort_idx,:));
        set(gca,'YDir','normal');
        % Plot stimulus onset
        line([0 0],ylim,'Color','k');
        % Plot ROLs
        d_scat = scatter(rol_times(ch_ix,rt_sort_idx,1),1:size(rol_times,2),rol_mrkr_sz,deriv_marker,'MarkerEdgeColor',deriv_color);
        l_scat = scatter(rol_times(ch_ix,rt_sort_idx,2),1:size(rol_times,2),rol_mrkr_sz,lin_marker,'MarkerEdgeColor',lin_color);
        % Plot ROL Means
        d_med = line([nanmedian(rol_times(ch_ix,:,1)) nanmedian(rol_times(ch_ix,:,1))],[1 size(rol_times,2)],...
            'Color',deriv_color,'LineStyle','--','LineWidth',2);
        l_med = line([nanmedian(rol_times(ch_ix,:,2)) nanmedian(rol_times(ch_ix,:,2))],[1 size(rol_times,2)],...
            'Color',lin_color,'LineStyle',':','LineWidth',2);
        % Plot RTs
        rt_scat = scatter(trial_info.response_time(rt_sort_idx),1:numel(trial_info.trial_n),rt_mrkr_sz,...
            rt_marker,'MarkerEdgeColor',rt_color,'MarkerFaceColor',rt_color);
        
        % Axis parameters
        ax = gca;
        ax.Title.String = [num2str(sum(~isnan(rol_times(ch_ix,:,1)))) ' ROLs (' ...
                           num2str(100*sum(~isnan(rol_times(ch_ix,:,1)))/size(rol_times,2),'%.01f'),...
                           '%); Error (S,R): Der = (',num2str(err(1,1),'%.01f'),',',num2str(err(1,2),'%.01f'),...
                           '); Lin = (',num2str(err(2,1),'%.01f'),',',num2str(err(2,2),'%.01f'),')'];
        ax.XLim = [hfa.time(1) hfa.time(end)];
        ax.YLim = [1 numel(trial_info.trial_n)];
        ax.XLabel.String = 'Time (s)';
        ax.YLabel.String = 'Trials';
        colorbar; caxis(clims);
        legend([d_scat l_scat d_med l_med rt_scat],{'Der','Lin','D avg','L avg','RT'},'Location','southeast');
        
        % PLOT ABOVE THRESHOLD EPOCHS
        subplot(1,2,2); hold on;
        % Plot single trial stack
        imagesc(hfa.time,1:numel(trial_info.trial_n),above_mat(rt_sort_idx,:));
        set(gca,'YDir','normal');
        % Plot stimulus onset
        line([0 0],ylim,'Color','k');
        % Plot ROL Means
        d_med = line([nanmedian(rol_times(ch_ix,:,1)) nanmedian(rol_times(ch_ix,:,1))],[1 size(rol_times,2)],...
            'Color',deriv_color,'LineStyle','--','LineWidth',2);
        l_med = line([nanmedian(rol_times(ch_ix,:,2)) nanmedian(rol_times(ch_ix,:,2))],[1 size(rol_times,2)],...
            'Color',lin_color,'LineStyle',':','LineWidth',2);
        % Plot RTs
        rt_scat = scatter(trial_info.response_time(rt_sort_idx),1:numel(trial_info.trial_n),rt_mrkr_sz,...
            rt_marker,'MarkerEdgeColor',rt_color,'MarkerFaceColor',rt_color);
        
        % Axis parameters
        ax = gca;
        ax.Title.String = ['Thresholded Epochs: ' num2str(mean(n_ep_above)) ' ep/trl'];
        ax.XLim = [hfa.time(1) hfa.time(end)];
        ax.YLim = [1 numel(trial_info.trial_n)];
        ax.XLabel.String = 'Time (s)';
        ax.YLabel.String = 'Trials';
        colorbar; caxis(clims);
        legend([d_med l_med rt_scat],{'D avg','L avg','RT'},'Location','southeast');
        
        % Save figure
        saveas(gcf,[stack_dir fig_name '.' fig_ftype]);
        if strcmp(fig_vis,'off')
            close(gcf);
        end
    end
    
end

%% Save ROL output
rol_fname = [SBJ_vars.dirs.proc SBJ '_ROI_' an_id '_' rol_id '.mat'];
save(rol_fname,'-v7.3','rol_times','evnt_lock','thresh','err','err_var','actv_len','actv_amp','rol_win_sz');

end
