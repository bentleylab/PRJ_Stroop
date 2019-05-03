%% ROL_estimates

% get relative onset latency (from stim onset) of gamma power
% code adapted from HBM paper on stopping, has both the Foster method and
% the derivative method

% Eleonora

%% PREPARE SIGNAL FOR INPUT:
%
% freq = vector of frequency values in Hz (eg freq = [2:1:201]; %Hz
% time = vector of time values in ms (0 is onset od stimulus):
% time=[-1000:1:1000]; msec
% signal: freq by time by trial matrix of power/amplitude values

% OUTUP: matrix of as many rows as trials, 2 columns (1st is regression
% method, second is derivative)
% 
%--------------------------------------------------------------------------
%example to see how it works: change here to get it to read your data
load('/Volumes/hoycw_clust/PRJ_Stroop/data/eleonora_ROL_test_data.mat');
freq=test_data.freq; time=test_data.time;
signal=test_data.spec_matrix;
%--------------------------------------------------------------------------

% function onset = ROL_vis_contrast(signal,freq, time)

% do you want to plot signle trial spectrograms and traces to visualize how
% the rol estimator is perfoming? Note: it will save a file for each
% trial: so if you have many trials this will take a long time, just use it
% in debugging

plot_and_save=1;

% create ROL folder: figures will be saved here:
save_dir = ['~/Desktop/Figures/ROL/'];
if ~isdir(save_dir)
    mkdir(save_dir)
end

% GET POWER within a given range & CROP TIME AROUND STIM ONSET:
% BBG: 
freq_range = [70 150];
%100 before stim on and 700 after:
time_range_msec = [-99 1000]; 

%% get the relevant signals from the signal matrix:
freq_range_ind = freq>=freq_range(1) & freq<=freq_range(2);
time_range = [find(time==time_range_msec(1)) find(time==time_range_msec(2))];
time_vec = time(time_range(1):time_range(2));

timetrial_vec=squeeze(mean(signal(freq_range_ind,time_range(1):time_range(2),:),1));

%% ROL routing: derivative method and regression method, init params

% timetrial_vec is time by trials, while my old code had trial by time, for
% sake of brevity, just transpose it and use the old code:

trialtime_vec=timetrial_vec';

%get the 75th percentile for each baseline-trial: will be our power values threshold
th=quantile(squeeze(mean(signal(freq_range_ind,:,:),1))',.75,2);

% min window length threshold (Foster2015 uses 100ms, my HBM uses 40ms)
time_th = 40; %40 msec = 40 samples, tf data is downsampled to 1kHz

%create empty array to store results: as many rows as trials, 2 columns:
onset=zeros(size(th,1),2);

for trialcounter =1:size(trialtime_vec,1)
    
    %POWER above threshold: logical index
    above_th=find(trialtime_vec(trialcounter,:)>th(trialcounter));
    
    %start and end value of the signal above threshold:
    %find rising and falling edges (using diff,add zeros before after to
    %pad) and check the window length between rising (1) and falling (-1)
    startind_above_th = find(diff([0 trialtime_vec(trialcounter,:)>th(trialcounter) 0])==1);
    endind_above_th = find(diff([0 trialtime_vec(trialcounter,:)>th(trialcounter) 0])==-1)-1;
    time_above_th = abs(startind_above_th-endind_above_th);
    
    %if at least one series of power values is above th for at least time th...
    % get onset values, otherwise, just spit out NaNs
    
    if any(time_above_th>time_th)
        
        %if more than one series is above the th,
        %get only those more than 40ms and get rid of others
        
        if numel(time_above_th)>1
            series_above=time_above_th>time_th;
            if sum(series_above)==1
                above_th = startind_above_th(series_above):endind_above_th(series_above);
                %just one longer than 40ms
            elseif sum(series_above)~=1
                %more than one peak: get the one with the biggest
                %%change
                series_ind=find(series_above);
                pow = zeros(numel(series_ind),1);
                for npeaks=1:numel(series_ind)
                    pow(npeaks) = mean(trialtime_vec(trialcounter,startind_above_th(series_ind(npeaks)):endind_above_th(series_ind(npeaks))));
                end
                [maxval, indmax]=max(pow);
                above_th = startind_above_th(series_ind(indmax)):endind_above_th(series_ind(indmax));
                %above_th=above_th(cumsum(time_above_th(1:series_above(1)-1))+1:end);
            end
        end
        
        %window around peak: 200 before and 100 after (Foster 2015)
        T1=200;
        T2=100;
        
        %define start and end of the window with respect to the first
        %value above th (indices, not actual time values):
        win_start=above_th(1)-T1;
        win_end=above_th(1)+T2;
        
        %check that we are not exceeding the trial window indices
        if win_start<=0
            win_start=1;
        end
        if win_end>size(trialtime_vec,2)
            win_end=size(trialtime_vec,2);
        end
        
        % define a power and time series based on this window
        max_pow_win=trialtime_vec(trialcounter,win_start:win_end);
        T_pow_win=time_vec(win_start:win_end);
        
        %prepare figure
        figure(trialcounter)
%         subplot(2,2,1)
        plot(T_pow_win,max_pow_win,'k','LineWidth',2)
        hold on
        %plot first order differential (multiply it by 10 to make it
        %easily visible) in RED
        plot(time_vec(win_start:win_end-1),diff(max_pow_win,1)*10,'r')
        
        %% DERIVATIVE METHOD
        
        %calculate first order differential: get max (=steeper slope)
        %and find the preceeding zerocrossing (local minimum before peak)
        derivative=diff(max_pow_win,1);
        [value,index_diff]=max(derivative);
        
        zerocross=find(diff(sign(derivative),1));
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
            onset(trialcounter,:)=[result(4) T_pow_win(index_diff_onset)];
        else
            onset(trialcounter,:)=[NaN T_pow_win(index_diff_onset)];
        end
        
    else
        %no significant time points? Return NaNs
        onset(trialcounter,:)=[NaN NaN];
    end
    
    if plot_and_save == 1
        %trial info
        trial_s=sprintf('%d',trialcounter);
        %plot singletrials spectrogram
        %figure(trialcounter)
        subplot(2,2,2)
        contourf(time_vec,find(freq_range_ind),squeeze(signal(find(freq_range_ind),time_range(1):time_range(2),trialcounter)),40,'linecolor','none');
        hold on; plot([onset(trialcounter,1) onset(trialcounter,1)],get(gca,'ylim'),'-w','Linewidth',2)
        plot([onset(trialcounter,2) onset(trialcounter,2)],get(gca,'ylim'),'-r','Linewidth',2)
        %plot singletrials gamma trace
        subplot(2,1,2)
        plot(time_vec,trialtime_vec(trialcounter,:),'k','LineWidth',2)
        title(['Trial number ' num2str(trialcounter)])
        hold on
        plot(get(gca,'xlim'),[th(trialcounter) th(trialcounter)])
        plot([0 0],get(gca,'ylim'),'-.k','Linewidth',1)
        plot([onset(trialcounter,1) onset(trialcounter,1)],get(gca,'ylim'),'-b','Linewidth',2)
        plot([onset(trialcounter,2) onset(trialcounter,2)],get(gca,'ylim'),'-r','Linewidth',2)
        
        saveas(gcf,[save_dir 'ROL_trial_' trial_s],'png')
        close
    end
    
end

             
% end