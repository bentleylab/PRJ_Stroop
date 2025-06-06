%% Load clustering results and plot stim-locked time histogram of spikes
if exist('/home/knight/hoycw/','dir');root_dir='/home/knight/hoycw/';ft_dir=[root_dir 'Apps/fieldtrip/'];
else root_dir='/Volumes/hoycw_clust/';ft_dir='/Users/colinhoy/Code/Apps/fieldtrip/';end
SBJ = 'IR75';
micro_ch_lab = 'mfoa4';
block_suffix = '';

%% Add paths
addpath([root_dir 'PRJ_Stroop/scripts/']);
addpath([root_dir 'PRJ_Stroop/scripts/utils/']);
addpath(genpath('/Users/colinhoy/Code/Apps/wave_clus/'));
addpath(genpath('/Users/colinhoy/Code/Apps/UR_NLX2MAT_releaseDec2015/'));
addpath(ft_dir);
ft_defaults

%% Load Data
% SBJ_vars
SBJ_vars_cmd = ['run ' root_dir 'PRJ_Stroop/scripts/SBJ_vars/' SBJ '_vars.m'];
eval(SBJ_vars_cmd);

% Trial Info
load([SBJ_vars.dirs.events SBJ '_trial_info_auto',block_suffix,'.mat']);

% Clustering Results
%   spikes [n_spikes, time] - waveforms
%   cluster_class [n_spikes, 2] - cluster class (int), then time in ms
load([SBJ_vars.dirs.preproc 'micro_clusters/semi_auto/times_' micro_ch_lab SBJ_vars.ch_lab.suffix '.mat']);
micro_hdr = ft_read_header([SBJ_vars.dirs.SU 'micro/' micro_ch_lab SBJ_vars.ch_lab.suffix '.ncs']);

%% Process data
% Adjust spike times from 1000 kHz time stamps to sec
cluster_class(:,2) = (cluster_class(:,2)-double(micro_hdr.FirstTimeStamp)/1000)/1000;
% Account for analysis time
if numel(SBJ_vars.analysis_time{b_ix})>1
    error('havent set up processing for multi block concat!');
end
cluster_class(:,2) = cluster_class(:,2) - SBJ_vars.analysis_time{b_ix}{1}(1);

%% Plot PSTH for each cluster across all conditions
units = sort(unique(cluster_class(:,1)));
clus_lab = {};
clus_colors = distinguishable_colors(numel(units));
spike_times = {[numel(units) 1]};
for clus_ix = 1:numel(units)
    spike_times{clus_ix} = cluster_class(cluster_class(:,1)==clus_ix,2);
    clus_lab{clus_ix} = ['Cluster ' num2str(clus_ix)];
end

psth_lim = [-0.25 2];
bin_sz = 0.05;    % histogram bins in s
bins = psth_lim(1):bin_sz:psth_lim(2);
psth = zeros([numel(units) numel(bins)-1]);
figure; hold on;
for clus_ix = 1:numel(units)
    for t_ix = 1:numel(trial_info.trial_n)
        for bin_ix = 1:numel(bins)-1
            if any(spike_times{clus_ix}>=trial_info.word_time(t_ix)+bins(bin_ix) & spike_times{clus_ix}<trial_info.word_time(t_ix)+bins(bin_ix+1))
                psth(clus_ix,bin_ix) = psth(clus_ix,bin_ix)+sum(spike_times{clus_ix}>=trial_info.word_time(t_ix)+bins(bin_ix) &...
                    spike_times{clus_ix}<trial_info.word_time(t_ix)+bins(bin_ix+1));
            end
        end
    end
    plot(bins(1:end-1),psth(clus_ix,:),clus_colors{clus_ix});
end
line([bins(bins==0) bins(bins==0)], ylim, 'Color','k','LineStyle','--');
xlabel('Time (s)');
ylabel('# Spikes');
legend(clus_lab,'Stim');

%% Plot PSTH for each cluster per condition
[cond_lab, cond_colors, ~] = fn_condition_label_styles('CNI');
cond_trials = {};
for cond_ix = 1:numel(cond_lab)
	cond_trials{cond_ix} = find(fn_condition_index([cond_lab{cond_ix}], trial_info.condition_n)==1);
end

psth_cond = zeros([numel(units) numel(cond_lab) numel(bins)-1]);
for clus_ix = 1:numel(units)
    figure; hold on;
    for cond_ix = 1:numel(cond_lab)
        for t_cond_ix = 1:numel(cond_trials{cond_ix})
            t_ix = cond_trials{cond_ix}(t_cond_ix);
            for bin_ix = 1:numel(bins)-1
                if any(spike_times{clus_ix}>=trial_info.word_time(t_ix)+bins(bin_ix) & spike_times{clus_ix}<trial_info.word_time(t_ix)+bins(bin_ix+1))
                    psth_cond(clus_ix,cond_ix,bin_ix) = psth_cond(clus_ix,cond_ix,bin_ix)+sum(spike_times{clus_ix}>=trial_info.word_time(t_ix)+bins(bin_ix) &...
                        spike_times{clus_ix}<trial_info.word_time(t_ix)+bins(bin_ix+1));
                end
            end
        end
        plot(bins(1:end-1),squeeze(psth_cond(clus_ix,cond_ix,:)),'Color',[cond_colors{cond_ix}]);
    end
    line([bins(bins==0) bins(bins==0)], ylim, 'Color','k','LineStyle','--');
    xlabel('Time (s)');
    ylabel('# Spikes');
    title(clus_lab{clus_ix});
    legend(cond_lab,'Stim');
end

