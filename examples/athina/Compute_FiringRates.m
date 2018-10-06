function [] = Compute_FiringRates(path_spikes, events_file, micro,Timing)

% path_spikes: the folder containing wave_clus output files
% events_file: a file containing an 'events' structure, fieldtrip-style
% micro: a 1x1 struct, with fields .Names: the names of microchannels we
% want to analyse, and .clusters2test: the number of cluster we want to
% assess for each microchannel
% Timing: a 1x1 struct, with fields .base: baseline length (ms) .post:
% poststimulus interval length (ms), .bin_l: length of binning windows
% (ms), .base_bin: length of baseline for binned data.


% Read events, to get trial info:
disp('Reading events...')
load(events_file)

fts = hdr.FirstTimeStamp;

ts = ([event(:).timestamp])/10^3; %timestamps
sample = ([event(:).sample]); %triggers recorded in neuralynx
cond = ([event(:).value]); % experimental conditions

nu_micro = length(micro.Names);
offset = ts(1) - uint64(1000*(Timing.post+Timing.base)); % subtract the offset of first time-stamp
ts = double(ts-offset);

% loop across all wires
for ch = 1:nu_micro
    
    con2check = char(micro.Names(ch));
    spfile = [path_spikes, 'times_', con2check, micro.file_ext, '.mat'];
    
    if exist(spfile) == 2
        
        load(spfile)
        
        spike_timing = cluster_class(:,2); % These are in ms already!!
        
        clusters = sort(unique(cluster_class(:,1)));
        % make a pseudo-time course of spiking activity, for each spike cluster:
        tc = zeros(length(clusters),round(ts(end)+1000));
        for cl = 1:length(clusters)
            cluster_timing = uint64(spike_timing(find(cluster_class(:,1) == clusters(cl)))) -offset;
            tc(cl,round(cluster_timing)) = 1;
        end
        
        trials = zeros(Timing.base+Timing.post,length(event),length(clusters));
        for tt = 1:length(event)
            trials(:,tt,:) = tc(:,round(ts(tt)-Timing.base):round(ts(tt)+Timing.post-1))';
        end
        
        % bin firing rates:
        [trials_binned, trials] = bin_spikes(trials, Timing.bin_l, Timing.base_bin, Timing.base);
        t_v = [-Timing.base/1000:1/(1000/Timing.bin_l):Timing.post/1000];
        t_v(end) = [];
                
        fig = figure;
        set(fig, 'Position', [100,100, 300, 800])
        for cl = 1:length(clusters)
            
            subplot(length(clusters),1,cl)
            bar(t_v, squeeze(mean(trials_binned(:,:,cl),2)), 'FaceColor', [0.5,0.5,0.5],'EdgeColor', [0.5,0.5,0.5] )
        end
        saveas(fig, [path_spikes, 'Histograms_', micro.p_id, '_', con2check, '.png'])
        
        close
        clear trials_binned spikes spike_timing fig cluster_class
    end
end


function [trials_binned,trials] = bin_spikes(trials, win,base_w,base_unbinned)

% downsampling factor:
downs = size(trials,1)/win;
trials_binned = zeros(downs, size(trials,2), size(trials,3));

cou = 1;
for ww = 1:win:size(trials,1)
    trials_binned(cou,:,:) = (sum(trials(ww:ww+win-1,:,:),1))/win;
    cou = cou+1;
end

base = squeeze(mean(trials_binned(1:base_w,:,:),1));
for ww = 1:downs
    trials_binned(ww,:,:) = (squeeze(trials_binned(ww,:,:))-base);
end

% baseline correct also the non binned data:
base = squeeze(mean(trials(1:base_unbinned,:,:),1));
for ww = 1:size(trials,1)
    trials(ww,:,:) = (squeeze(trials(ww,:,:))-base);
end
