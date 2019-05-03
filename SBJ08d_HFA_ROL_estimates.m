function SBJ08d_HFA_ROL(SBJ, proc_id, an_id, rol_id, plot_example)
% Find single trial onset latencies of HFA
%   Adapted from code by Eleonora Bartoli and Brett Foster
% INPUTS:
%   SBJ
%   proc_id
%   an_id
%   plot_example [0/1] - plot example trials? takes forever!
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

% Output prep
fig_dir = [SBJ_vars.dirs.proc 'ROL/' an_id '/'];
if ~isdir(fig_dir)
    mkdir(fig_dir);
end

%% Compute parameters for ROLs
% Power threshold per channel over entire trial
thresh = quantile(reshape(hfa.powspctrm,[numel(hfa.label) size(hfa.powspctrm,1)*size(hfa.powspctrm,4)]),quant_thresh,2);

% Time threshold
%!!! fix this naming for the different windows!
sample_rate = (numel(hfa.time)-1)/(hfa.time(end)-hfa.time(1));
win_dp = round(win_len*sample_rate);

%% Select time range for onsets
rol_lim = [0 max(trial_info.response_time)+0.7];
cfgs = [];
cfgs.latency = rol_lim;
hfa = ft_selectdata(cfgs, hfa);

%% Compute ROLs
rol = zeros([numel(hfa.label) size(hfa.powspctrm,1) 2]);
for ch_ix = 1:numel(hfa.label)
    for trl_ix = 1:size(hfa.powspctrm,1)
        trl_data = squeeze(hfa.powspctrm(trl_ix,ch_ix,:,:))';    % column
        % Find times above threshold
%         above_th  = find(trl_data>thresh(trl_ix));
%!!! convert to (ep,2) column vectors instead of 2 variables!
        beg_above = find(diff([0 trl_data>thresh(trl_ix) 0])==1);
        end_above = find(diff([0 trl_data>thresh(trl_ix) 0])==-1)-1;
        len_above = end_above-beg_above;
        
        % Remove epochs shorter than win_len
        beg_above = beg_above(len_above>win_dp);
        end_above = end_above(len_above>win_dp);
        
        %% Compute ROLs for windows above threshold for at least win_len
        if any(beg_above)
            % Take epoch with largest power change
            if numel(beg_above)>1
                ep_mean = zeros(size(beg_above));
                for ep_ix = 1:numel(beg_above)
                    ep_mean = mean(trl_data(beg_above(ep_ix):end_above(ep_ix)));
                end
                [~, max_ix] = max(ep_mean);
                ep_idx = beg_above(max_ix):end_above(max_ix);
            else
                ep_idx = beg_above:end_above;
            end
            
            % Define start and end of onset search epoch
            rol_lim = round([ep_idx(1)+rol_lim_s(1)*sample_rate ...
                             ep_idx(1)+rol_lim_s(2)*sample_rate]);
            
            % Check search epoch is within trial limits
            if rol_lim(1)<=0
                rol_lim(1) = 1;
            end
            if rol_lim(2)>numel(trl_data)
                rol_lim(2) = numel(trl_data);
            end
            
            % Grab search epoch
            search_ep = trl_data(rol_lim(1):rol_lim(2));
            
            % Prepare figure to plot ROL checks
            figure(trl_ix);
            % Plot search epoch
            plot(hfa.time(rol_lim(1):rol_lim(2)),trl_data(rol_lim(1):rol_lim(2)),'k','LineWidth',2);
            hold on;
            % Plot first order derivative (multiply it by 10 to make it easily visible)
            plot(hfa.time(rol_lim(1):rol_lim(2)-1),diff(trl_data(rol_lim(1):rol_lim(2)),1)*10,'r');
            
            %% DERIVATIVE METHOD
            !!! start again here
            %calculate first order differential: get max (=steeper slope)
            %and find the preceeding zerocrossing (local minimum before peak)
            deriv = diff(trl_data(rol_lim(1):rol_lim(2)),1);
            [value,index_diff]=max(deriv);
            
            zerocross=find(diff(sign(deriv),1));
            zerocross_before_peak=find(zerocross<index_diff);
            if ~isempty(zerocross_before_peak)
                index_diff_zerocross=zerocross(zerocross_before_peak(end));
                index_diff_onset=round(median([index_diff_zerocross index_diff]));
            else
                index_diff_zerocross=1;
                index_diff_onset=round(median([index_diff_zerocross index_diff]));
            end
            
            plot(T_pow_win(index_diff),value*10,'o-r','LineWidth',4)
            plot(T_pow_win(index_diff_zerocross),value*10,'o-r','LineWidth',4)
            plot([T_pow_win(index_diff_onset) T_pow_win(index_diff_onset)] ,get(gca,'ylim'),'r')
            
            %% linear model as in the Neuron Foster paper:
            
            % initialize the start of the time window
            % 50 samples wide, move the window of 10 samples
            % each time
            % Foster 2015: 100ms wins, moved in 20 steps (90% overlap)
            
            T3=51; %changed from 200 to 100 %to 50 samples
            T0=1; % this number will be increased in the cycle below
            
            %define everything in samples: how many sliding steps do you
            %need to cover the T1-T2 window?
            size_tot_win=win_end-win_start; %300samples or less
            size_mov_win=T3; %50samples
            size_movement=10; %10samples
            num_tw=round((size_tot_win-size_mov_win)/size_movement);
            if num_tw==0
                num_tw=1;
            end
            
            %prepare empty array for the results of the lm regression
            lm_fit=zeros(num_tw,5);
            % start cycle on moving windows
            for tw=1:num_tw
                %if window size is not big enough, depending
                %on previous adjustments
                if T0+T3>size(T_pow_win,2)
                    T3=size(T_pow_win,2)-T0;
                end
                % linear fit: first coef is slope, second is
                % intercept
                coeff=polyfit(T_pow_win(T0:T0+T3),max_pow_win(T0:T0+T3),1);
                pow_predict=coeff(1)*T_pow_win(T0:T0+T3)+coeff(2);
                abs_err=abs(max_pow_win(T0:T0+T3)-pow_predict);
                %store slope, intercept, max error-residual,
                %start of window and end of window
                lm_fit(tw,:)=[coeff max(abs_err) time_vec(win_start+T0) time_vec(win_start+T0+T3)];
                %plot the fit
                plot(T_pow_win(T0:T0+T3),pow_predict,'o')
                %move the time window forward of 10 samples-10msec
                T0=T0+size_movement;
            end
            
            %sort the slopes from biggest to smallest
            [sorted_slopes,sorting_slopes]=sort(lm_fit(:,1),'descend');
            
            % avoid getting negative slopes, even if they are big:
            increasing_fit = sign(sorted_slopes)==1;
            sorting_slopes = sorting_slopes(increasing_fit);
            sorted_slopes = sorted_slopes(increasing_fit);
            
            %take five biggest slopes
            if numel(sorting_slopes)>=5
                five_big_slopes=sorting_slopes(1:5);
            else
                five_big_slopes=sorting_slopes(1:numel(sorting_slopes));
            end
            
            %and sort the erros
            [sorted_error,sorting_error]=sort(lm_fit(five_big_slopes,3),'ascend');
            %get smallest error data, if any
            
            if numel(sorting_slopes)>0
                result=lm_fit(five_big_slopes(sorting_error(1)),:);
                
                %plot lm result
                plot([result(4) result(4)],[get(gca,'ylim')],'b')
                
                %% Store regression fit as the first value, diff as the second
                onset(trl_ix,:)=[result(4) T_pow_win(index_diff_onset)];
            else
                onset(trl_ix,:)=[NaN T_pow_win(index_diff_onset)];
            end
            
        else
            % No windows above threshold for win_len
            rol(trl_ix,:) = [NaN NaN];
        end
        
        if plot_and_save == 1
            %trial info
            trial_s=sprintf('%d',trl_ix);
            %plot singletrials spectrogram
            %figure(trialcounter)
            subplot(2,2,2)
            contourf(time_vec,find(freq_range_ind),squeeze(signal(find(freq_range_ind),time_range(1):time_range(2),trl_ix)),40,'linecolor','none');
            hold on; plot([onset(trl_ix,1) onset(trl_ix,1)],get(gca,'ylim'),'-w','Linewidth',2)
            plot([onset(trl_ix,2) onset(trl_ix,2)],get(gca,'ylim'),'-r','Linewidth',2)
            %plot singletrials gamma trace
            subplot(2,1,2)
            plot(time_vec,trialtime_vec(trl_ix,:),'k','LineWidth',2)
            title(['Trial number ' num2str(trl_ix)])
            hold on
            plot(get(gca,'xlim'),[th(trl_ix) th(trl_ix)])
            plot([0 0],get(gca,'ylim'),'-.k','Linewidth',1)
            plot([onset(trl_ix,1) onset(trl_ix,1)],get(gca,'ylim'),'-b','Linewidth',2)
            plot([onset(trl_ix,2) onset(trl_ix,2)],get(gca,'ylim'),'-r','Linewidth',2)
            
            saveas(gcf,[save_dir 'ROL_trial_' trial_s],'png')
            close
        end
        
    end
end

end
