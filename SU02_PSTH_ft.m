%% Fieldtrip-based PSTH analysis
if exist('/home/knight/hoycw/','dir');root_dir='/home/knight/hoycw/';ft_dir=[root_dir 'Apps/fieldtrip/'];
else root_dir='/Volumes/hoycw_clust/';ft_dir='/Users/colinhoy/Code/Apps/fieldtrip/';end
SBJ = 'IR75';
an_id = 'PSTH_S_trl2to150_bn20';
conditions = 'CNI';
pipeline_id = 'SU_nlx';

%% Add paths
addpath([root_dir 'PRJ_Stroop/scripts/']);
addpath([root_dir 'PRJ_Stroop/scripts/utils/']);
addpath(genpath('/Users/colinhoy/Code/Apps/wave_clus/'));
addpath(genpath('/Users/colinhoy/Code/Apps/UR_NLX2MAT_releaseDec2015/'));
addpath(ft_dir);
ft_defaults

%% Processing Variables
an_vars_cmd = ['run ' root_dir 'PRJ_Stroop/scripts/an_vars/' an_id '_vars.m'];
eval(an_vars_cmd);
proc_vars_cmd = ['run ' root_dir 'PRJ_Stroop/scripts/proc_vars/' pipeline_id '_proc_vars.m'];
eval(proc_vars_cmd);
SBJ_vars_cmd = ['run ' root_dir 'PRJ_Stroop/scripts/SBJ_vars/' SBJ '_vars.m'];
eval(SBJ_vars_cmd);
if numel(SBJ_vars.analysis_time)>1
    error('havent set up multi-run processing yet!');
elseif numel(SBJ_vars.analysis_time{1})>1
    error('havent set up processing for multi block concat!');
end

%% Load Trial Info
load([SBJ_vars.dirs.events SBJ '_trial_info_auto.mat']);
% Convert events to seconds
if strcmp(event_lab,'S')
    events = trial_info.word_onset/proc_vars.resample_freq;
elseif strcmp(event_lab,'R')
    events = trial_info.resp_onset/proc_vars.resample_freq;
else
    error(['Bad event_lab: ' event_lab]);
end
% Compile cond_type, RT, trial_n
[cond_lab, cond_colors, ~] = fn_condition_label_styles(conditions);
cond_idx = false([length(cond_lab) length(trial_info.trial_n)]);
for cond_ix = 1:length(cond_lab)
    % Get binary condition index
    cond_idx(cond_ix,:) = logical(fn_condition_index(cond_lab{cond_ix},...
        trial_info.condition_n));
end

[cond_mat,~] = find(cond_idx);
cond_mat = horzcat(cond_mat,trial_info.response_time,[1:numel(cond_mat)]');
cond_mat = sortrows(cond_mat,[1 2]);
% RT_mean = mean(round(1000*trial_info.response_time)); % converts sec to ms
% % Add in the baseline offset to plot correctly
% RT_mean = RT_mean-plt_vars.plt_lim_S(1)*1000;

%% Load the output of semi-automatic clustering with wave_clus
spike = fn_load_wave_clus_ft(SBJ,'');

%% Convert into trials
% Define spike times relative to trial events, sorted by condition
cfg = [];
cfg.trlunit = 'timestamps';
cfg.timestampspersecond = 1;
cfg.trl = [events(cond_mat(:,3))+trial_lim_s(1),...        % start of the trial including baseline
           events(cond_mat(:,3))+trial_lim_s(2),...    % end of trial
           repmat(trial_lim_s(1),length(events),1),...         %time of event relative to start of trial
           trial_info.condition_n(cond_mat(:,3))'];                 % trial type
spike_trl = ft_spike_maketrials(cfg,spike);

%% Get rid of units with insufficient number of spikes
if ~isnumeric(min_spikes)
    if ~strcmp(min_spikes,'n_trials'); error('unknown min_spikes');end
    min_spikes = numel(trial_info.word_onset);
end
too_small = true(size(spike_trl.label));
for u = 1:numel(spike_trl.label)
    if numel(spike_trl.timestamp{u})>=min_spikes
        too_small(u) = false;
    end
end
cfgs = [];cfgs.channel = spike_trl.label(~too_small);
spike_trl = ft_selectdata(cfgs,spike_trl);

%% Characterize ISI Distribution
if plot_ISI
    cfg = [];
    cfg.bins = [0:0.0005:0.1];
    cfg.param = 'coeffvar';     % coefficient of variation (sd/mn of ISIs)
    isih = ft_spike_isi(cfg,spike_trl);
    %   isih.isi - ISI per spike (with respect to previous spike)
    %   isih.avg - average ISI histogram per unit
    %   isih.coeffvar - summary stat of ISI histogram (see Shinomoto et al., 2009 apparently...)
    
    fig_dir = [SBJ_vars.dirs.preproc 'micro_clusters/semi_auto/ISI_dist/' an_id '/'];
    if ~exist(fig_dir,'dir'); mkdir(fig_dir); end
    for u = 1:numel(spike_trl.label)
        cfg              = [];
        cfg.spikechannel = isih.label{u};
        cfg.interpolate  = 5; % interpolate at 5 times the original density
        cfg.window       = 'gausswin'; % use a gaussian window to smooth
        cfg.winlen       = 0.004; % the window by which we smooth has size 4 by 4 ms
        cfg.colormap     = jet(300); % colormap
        cfg.scatter      = 'no'; % do not plot the individual isis per spike as scatters
        fig_name = [SBJ '_' isih.label{u} '_ISI_dist'];
        f = figure('Name',fig_name,'Visible','off');
        ft_spike_plot_isireturn(cfg,isih);
        title(['Total spikes = ' num2str(numel(spike_trl.timestamp{u}))]);
        saveas(gcf,[fig_dir fig_name '.png']);
        close(f);
    end
end

%% Compute PSTH
cfg             = []; 
cfg.binsize     = bin_sz; % in sec or {'scott','sqrt'} to estimate the optimal bin size from the data
cfg.outputunit  = 'rate'; % {'rate','spikecount','proportion'}
cfg.latency     = trial_lim_s;
cfg.vartriallen = 'yes';
cfg.keeptrials  = 'no'; %'yes' for stats later?
% cfg.trials      = find(trial_info.condition_n==cond_n);
psth = ft_spike_psth(cfg,spike_trl);
%   psth.time        = center histogram bin points
%  	psth.fsample     = 1/binsize;
%   psth.avg         = contains average PSTH per unit 
%   psth.trial       = contains PSTH per unit per trial 
%   psth.var         = contains variance of PSTH per unit across trials

%% Plot Raster
fig_dir = [root_dir 'PRJ_Stroop/results/SU/raster/' an_id '/'];
if ~exist(fig_dir,'dir'); mkdir(fig_dir); end

cfg              = [];
cfg.topplotfunc  = 'line'; % plot as a line
cfg.latency      = 'maxperiod';%trial_lim_s;
cfg.errorbars    = 'std'; % plot with the standard deviation
cfg.interactive  = 'no'; % toggle off interactive mode
for u = 1:numel(spike_trl.label)
    cfg.spikechannel = spike_trl.label(u);
    fig_name = [SBJ '_' spike_trl.label{u} '_PSTH_raster_' event_lab];
    figure('Name',fig_name,'Visible','on');
    f = ft_spike_plot_raster(cfg, spike_trl, psth);
    
    % Add RTs
    axes(f.hdl.axRaster);
    for cond_ix = 1:numel(cond_lab)
        idx = cond_mat(:,1)==cond_ix;
        scatter(cond_mat(idx,2)-trial_lim_s(1),find(idx),'.',...
            'MarkerEdgeColor',[cond_colors{cond_ix}]);%,'MarkerEdgeColor');
    end
    
    % Save
    saveas(gcf,[fig_dir fig_name '.png']);
    close(gcf);
end
