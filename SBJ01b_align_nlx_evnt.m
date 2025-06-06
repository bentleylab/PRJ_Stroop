function SBJ01b_align_nlx_evnt(SBJ, proc_id, block_ix, save_it)
% Load, resample, and align photodiode and microphone from neuralynx system to clinical data
% save_it == 0: don't save plots, compare to raw data
%         == 1: save plots and data, compare to import
if block_ix~=1
    error('SBJ01b not ready for multi-block runs yet!');
end

%% Load, preprocess, and save out photodiode
if exist('/home/knight/','dir');root_dir='/home/knight/';app_dir=[root_dir 'hoycw/Apps/'];
else root_dir='/Volumes/hoycw_clust/';app_dir='/Users/colinhoy/Code/Apps/';end

%% Paths
addpath([root_dir 'PRJ_Stroop/scripts/']);
addpath([root_dir 'PRJ_Stroop/scripts/utils/']);
addpath(genpath([app_dir 'wave_clus/']));
addpath(genpath([app_dir 'UR_NLX2MAT_releaseDec2015/']));
addpath([app_dir 'fieldtrip/']);
ft_defaults

%% SBJ vars
eval(['run ' root_dir 'PRJ_Stroop/scripts/SBJ_vars/' SBJ '_vars.m']);
if ~isfield(SBJ_vars.dirs,'nlx')
    error([SBJ ' does not have NLX data, run only SBJ01a!']);
end
if numel(SBJ_vars.analysis_time{block_ix})>1
    error('havent set up processing for multi block concat within a run!');
end

if numel(SBJ_vars.block_name)>1
    block_suffix = ['_' SBJ_vars.block_name{block_ix}];
else
    block_suffix = SBJ_vars.block_name{block_ix};   % should just be ''
end

eval(['run ' root_dir 'PRJ_Stroop/scripts/proc_vars/' proc_id '_vars.m']);
if length(SBJ_vars.ch_lab.nlx_nk_align)~=2, error('why not 2 NLX-NK align channels?'); end

%% Read photodiode, NLX macro, clinical data
% First check for invalid samples
check_neuralynx_validsamples([SBJ_vars.dirs.nlx 'photo/']);

macro_fnames = SBJ_vars.ch_lab.nlx_nk_align;
for m_ix = 1:numel(macro_fnames)
    macro_fnames{m_ix} = [SBJ_vars.dirs.nlx 'macro/' SBJ_vars.ch_lab.nlx_nk_align{m_ix}...
        SBJ_vars.ch_lab.nlx_suffix '.ncs'];
end
try
    % Neuralynx photodiode
    evnt       = ft_read_neuralynx_interp({[SBJ_vars.dirs.nlx 'photo/' ...
        SBJ_vars.ch_lab.photod{block_ix} SBJ_vars.ch_lab.nlx_suffix '.ncs']});
    % Neuralynx mic
    mic       = ft_read_neuralynx_interp({[SBJ_vars.dirs.nlx 'mic/' ...
        SBJ_vars.ch_lab.mic{block_ix} SBJ_vars.ch_lab.nlx_suffix '.ncs']});
    % Neuralynx macro channel
    macro = ft_read_neuralynx_interp(macro_fnames);
catch ME
    if strcmp(ME.message, 'Sample points must be unique.')
        fprintf(2,'ft_read_neuralynx_interp error, trying ft_preprocessing...');
        cfg_load = [];
        cfg_load.dataset = [SBJ_vars.dirs.nlx 'photo/' ...
            SBJ_vars.ch_lab.photod{block_ix} SBJ_vars.ch_lab.nlx_suffix '.ncs'];
        evnt = ft_preprocessing(cfg_load);
        % Neuralynx mic
        cfg_load.dataset = [SBJ_vars.dirs.nlx 'mic/' ...
            SBJ_vars.ch_lab.mic{block_ix} SBJ_vars.ch_lab.nlx_suffix '.ncs'];
        mic       = ft_preprocessing(cfg_load);
        % Neuralynx macro channel
        cfg_load.dataset = macro_fnames{1};
        macro1 = ft_preprocessing(cfg_load);
        cfg_load.dataset = macro_fnames{2};
        macro2 = ft_preprocessing(cfg_load);
        macro = ft_appenddata([],macro1,macro2);
    else % Re-throw the error if it's not the specific one you're looking for
        rethrow(ME);
    end
end
macro.label = SBJ_vars.ch_lab.nlx_nk_align;

% Deal with multi-block within a run (error for now...)
% if numel(SBJ_vars.analysis_time{block_ix})>1
%     data_blocks = {};
%     for block_ix = 1:numel(SBJ_vars.block_name)
%         load([SBJ_vars.dirs.import SBJ '_' srate_str 'hz_' SBJ_vars.block_name{block_ix} '.mat']);
%         data_blocks{block_ix} = data;
%     end
%     data = fn_concat_blocks(data_blocks);
% end

% Nihon Kohden clinical channel
if ~save_it
    clin_fname = SBJ_vars.dirs.raw_filename{block_ix};
else
    if any(SBJ_vars.low_srate)
        srate_str = num2str(SBJ_vars.low_srate(block_ix));
    else
        srate_str = num2str(proc.resample_freq);
    end
    clin_fname = [SBJ_vars.dirs.import SBJ '_' srate_str 'hz' block_suffix '.mat'];
end
load(clin_fname);

cfgs         = [];
cfgs.channel = SBJ_vars.ch_lab.nlx_nk_align;
clin         = ft_selectdata(cfgs,data);
% clin_orig    = clin;

%% Preprocess
% Cut NLX to nlx_analysis_time
if isfield(SBJ_vars,'nlx_analysis_time')
    if numel(SBJ_vars.nlx_analysis_time{block_ix})>1
        error('not built to handle multiple cuts within a single nlx run!');
    end
    cfgs = [];
    cfgs.latency = SBJ_vars.nlx_analysis_time{block_ix}{1};
    evnt  = ft_selectdata(cfgs, evnt);
    mic   = ft_selectdata(cfgs, mic);
    macro = ft_selectdata(cfgs, macro);
end
if any(isnan(evnt.trial{1}(:))) || any(isnan(mic.trial{1}(:))) || any(isnan(macro.trial{1}(:)))
    error('NaNs detected in NLX data, check for discontinuities!');
end

% Inversion on NLX data
if SBJ_vars.nlx_macro_inverted
    macro.trial{1} = macro.trial{1}*-1;
end
if SBJ_vars.photo_inverted
    evnt.trial{1} = evnt.trial{1}*-1;
    mic.trial{1}  = mic.trial{1}*-1;
end

% Set aside full sample rate mic for listening
mic_data = mic.trial{1};
mic_full_srate = mic.fsample;
    
% Bipolar rereference
if numel(SBJ_vars.ch_lab.nlx_nk_align)>1
    cfg = [];
    cfg.channel = macro.label;
    cfg.montage.labelold = cfg.channel;
    [cfg.montage.labelnew, cfg.montage.tra, ~] = fn_create_ref_scheme_bipolar(SBJ_vars.ch_lab.nlx_nk_align);
    macro = ft_preprocessing(cfg,macro);
    
    cfg = [];
    cfg.channel = clin.label;
    cfg.montage.labelold = cfg.channel;
    [cfg.montage.labelnew, cfg.montage.tra, ~] = fn_create_ref_scheme_bipolar(SBJ_vars.ch_lab.nlx_nk_align);
    clin = ft_preprocessing(cfg,clin);
end

% % Cut clincial macro to analysis time
% if numel(SBJ_vars.analysis_time{block_ix})>1
%     error('havent set up processing for multi block concat!');
% end
cfgs = []; cfgs.latency = SBJ_vars.analysis_time{block_ix}{1};
clin = ft_selectdata(cfgs,clin);
clin.time{1} = clin.time{1}-SBJ_vars.analysis_time{block_ix}{1}(1);

% Match sampling rates
if macro.fsample > clin.fsample
    fprintf('downsampling Neuralynx from %d to %d Hz\n', macro.fsample, clin.fsample)
    cfgr = [];
    cfgr.demean     = 'yes';
    cfgr.resamplefs = clin.fsample;
    macro = ft_resampledata(cfgr, macro);
elseif macro.fsample < clin.fsample
    error('Clinical data has higher sampling rate than NLX macro, whats going on?');
end
if evnt.fsample > clin.fsample
    fprintf('downsampling Neuralynx photodiode from %d to %d Hz\n', evnt.fsample, clin.fsample)
    cfgr = [];
    cfgr.demean     = 'yes';
    cfgr.resamplefs = clin.fsample;
    evnt = ft_resampledata(cfgr, evnt);
end
if mic.fsample > clin.fsample
    fprintf('downsampling Neuralynx photodiode from %d to %d Hz\n', mic.fsample, clin.fsample)
    cfgr = [];
    cfgr.demean     = 'yes';
    cfgr.resamplefs = clin.fsample;
    mic = ft_resampledata(cfgr, mic);
end

% % Remove extreme values
% clin_thresh  = proc.nlx_nk_align_std_thresh*std(clin.trial{1});
% macro_thresh = proc.nlx_nk_align_std_thresh*std(macro.trial{1});
% clin.trial{1}((clin.trial{1}>median(clin.trial{1})+clin_thresh)|(clin.trial{1}<median(clin.trial{1})-clin_thresh)) = median(clin.trial{1});
% macro.trial{1}((macro.trial{1}>median(macro.trial{1})+macro_thresh)|(macro.trial{1}<median(macro.trial{1})-macro_thresh)) = median(macro.trial{1});

%% Compare PSDs
fn_plot_PSD_1by1_compare(clin.trial{1},macro.trial{1},clin.label,macro.label,...
    clin.fsample,'clinical','macro');
if save_it
    saveas(gcf,[SBJ_vars.dirs.import SBJ '_nlx_nk_PSD_compare_' macro.label{1} block_suffix '.png']);
end

%% Compute cross correlation at varying time lags
if numel(clin.trial{1}) <= numel(macro.trial{1})
    warning('Clinical data is smaller than macro data, recut clinical block!');
end

% Arjen's way: synchronize nihon kohden and neuralynx timeseries
%   Find cross-variance and lags by shifting macro along clin
[covar, lags] = xcov(clin.trial{1}', macro.trial{1}');
%   Find the peak (previously: get sharp peak; remove very long trends with
%       smooth (moving average) to find big spike in covariance)
[~, idx]      = max(covar);% - smooth(covar, clin.fsample*10));   
if isfield(SBJ_vars,'nlx_nk_align_force')
    algorithm_idx = idx;
    idx = SBJ_vars.nlx_nk_align_force;
end

%% Plot match
figure; subplot(2,1,1);
hold on; plot(lags,covar);
hold on; plot(lags(idx),covar(idx),'k*');
if exist('algorithm_idx','var')
    plot(lags(algorithm_idx),covar(algorithm_idx),'r*');
    legend('corrected','alogrithm');
else
    legend('algorithm');
end
ylabel('correlation');
xlabel('lag');
subplot(2,1,2);
t = 1:numel(clin.time{1}); 
hold on; plot(t, zscore(clin.trial{1}));
t2 = lags(idx):lags(idx)+numel(macro.trial{1})-1;
hold on; plot(t2, zscore(macro.trial{1})+10);
t3 = lags(idx):lags(idx)+numel(evnt.time{1})-1; % ignore any offset between photo and chan
hold on; plot(t3, zscore(evnt.trial{1})+20);
legend('NK', 'NLX', 'NLX photo');
title(macro.label{1});

%% Save figure and data
if save_it
    saveas(gcf,[SBJ_vars.dirs.import SBJ '_nlx_nk_macro_alignment_' macro.label{1} '.fig']);
    % print([subj(1).datadir 'datafiles/sync_nk-nl_' subjectm(9:end) '_' num2str(tcgver) '_' num2str(d)], '-dpdf');
    
    %% Create photodiode channel matched to clinical data
    evnt_nlx   = evnt;
    mic_nlx    = mic;
    evnt       = clin;
    evnt.label = {SBJ_vars.ch_lab.photod{1}, SBJ_vars.ch_lab.mic{1}};
    % Create dummy time series of the median of the photodiode
    evnt.trial{1} = [ones(1,numel(clin.trial{1})).*median(evnt.trial{1});
                     ones(1,numel(clin.trial{1})).*median(mic.trial{1})];
    % Add in photodiode data for segments when NLX overlaps
    evnt.trial{1}(1,t3(t3>0 & t3<numel(evnt.trial{1}))) = evnt_nlx.trial{1}(t3>0 & t3<numel(evnt.trial{1}));
    evnt.trial{1}(2,t3(t3>0 & t3<numel(evnt.trial{1}))) = mic_nlx.trial{1}(t3>0 & t3<numel(evnt.trial{1}));
    
    %% Save out mic.wav for listening
    % Rescale to prevent clipping, add 0.05 fudge factor
    mic_data_rescale = mic_data./(max(abs(mic_data))+0.05);
    mic_data_fname = strcat(SBJ_vars.dirs.import,SBJ,'_mic_recording',block_suffix,'.wav');
    fprintf('Saving %s\n',mic_data_fname);
    audiowrite(mic_data_fname,mic_data_rescale,mic_full_srate);
    
    %% Save data out
    evnt_out_fname = strcat(SBJ_vars.dirs.import,SBJ,'_evnt',block_suffix,'.mat');
    fprintf('Saving %s\n',evnt_out_fname);
    save(evnt_out_fname, '-v7.3', 'evnt');
end

end
