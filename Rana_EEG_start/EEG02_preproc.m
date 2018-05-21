%% EEG02_Preprocessing
SBJ = 'IR57';
pipeline_id = 'eeg_ft';

% Parameters
psd_filt     = 'yes';          % plot psds after filtering?
psd_reref    = 'yes';          % plot psds after rereferencing?
psd_fig_type = 'jpg';

%% Add paths
addpath('/home/knight/hoycw/PRJ_Stroop/scripts/');
addpath('/home/knight/hoycw/PRJ_Stroop/scripts/utils/');
addpath('/home/knight/hoycw/Apps/fieldtrip/');
ft_defaults

%% Load the data
eval(['run /home/knight/hoycw/PRJ_Stroop/scripts/SBJ_vars/' SBJ '_vars.m']);
eval(['run /home/knight/hoycw/PRJ_Stroop/scripts/proc_vars/' pipeline_id '_proc_vars.m']);

eeg_filename = strcat(SBJ_vars.dirs.import,SBJ,'_eeg_',num2str(proc_vars.resample_freq),'hz.mat');
load(eeg_filename);
eog_filename = strcat(SBJ_vars.dirs.import,SBJ,'_eog_',num2str(proc_vars.resample_freq),'hz.mat');
load(eog_filename);%data = ft_preprocessing(cfg);

bad_filename = [SBJ_vars.dirs.events SBJ '_eeg_bad_epochs_preclean.mat'];
load(bad_filename);

%% Filtering
fprintf('============== High and Low Pass Filtering %s ==============\n',SBJ);
bs_freq_lim = NaN([length(SBJ_vars.notch_freqs) 2]);
for f_ix = 1:length(SBJ_vars.notch_freqs)
    bs_freq_lim(f_ix,:) = fn_freq_lim_from_CFBW(SBJ_vars.notch_freqs(f_ix),SBJ_vars.bs_width);
end
cfg          = [];
cfg.bsfilter = 'yes';
cfg.bsfreq   = bs_freq_lim;
cfg.demean   = proc_vars.demean_yn;
cfg.lpfilter = proc_vars.lp_yn;
cfg.lpfreq   = proc_vars.lp_freq;
cfg.hpfilter = proc_vars.hp_yn;
cfg.hpfreq   = proc_vars.hp_freq;
if isfield(proc_vars,'hp_order')
    cfg.hpfiltord = proc_vars.hp_order;
end
eeg_filt = ft_preprocessing(cfg,eeg);
% data_bp = eeg_filt;

if strcmp(psd_filt,'yes')
    psd_dir = strcat(SBJ_vars.dirs.preproc,'PSDs/eeg_filt/');
    if ~exist(psd_dir,'dir')
        mkdir(psd_dir);
    end
    fn_plot_PSD_1by1_compare_save(eeg.trial{1},eeg_filt.trial{1},eeg.label,eeg_filt.label,...
        eeg_filt.fsample,strcat(psd_dir,SBJ,'_eeg_PSD_filt'),'orig','filt',psd_fig_type);
end


% %% Re-referencing
% cfg = [];
% cfg.reref      = 'yes';
% cfg.refchannel = setdiff(eeg.label,SBJ_vars.ch_lab.eeg_bad);
% cfg.refmethod  = 'avg';
% eeg_filt_car = ft_preprocessing(cfg, eeg_filt);
% 
% if strcmp(psd_car,'yes')
%     psd_dir = strcat(SBJ_vars.dirs.preproc,'PSDs/eeg/');
%     if ~exist(psd_dir,'dir')
%         mkdir(psd_dir);
%     end
%     fn_plot_PSD_1by1_compare_save(eeg_car.trial{1},eeg_car_filt.trial{1},eeg_car.label,eeg_car_filt.label,...
%         eeg_car_filt.fsample,strcat(psd_dir,SBJ,'_eeg_PSD_car_filt'),'car','car_filt',psd_fig_type);
% end

%% CZ laplacian
% if isempty(strmatch('C3',eeg_filt.label)) || isempty(strmatch('C4',eeg_filt.label))
%     error('ERROR: No C3 or C4 detected! You need a different refrencing strategy...');
% end
cfgs = [];
cfgs.channel = SBJ_vars.ch_lab.CZ_lap_ref;
cz_neigh = ft_selectdata(cfgs,eeg_filt);
cz_neigh_car = mean(cz_neigh.trial{1},1);
% cfgm = [];
% cfgm.parameter = 'trial';
% cfgm.operation = '(x1+x2)/2';
% cz_neigh_car = ft_math(cfgm,cz_neigh);
% 
cfgs.channel = 'CZ';
cz_lap = ft_selectdata(cfgs,eeg_filt);
% cfg = [];
% cfg.operation = 'x1';
% cz_lap = ft_math(cfg,cz,;
cz_lap.trial{1} = cz_lap.trial{1}-cz_neigh_car;

%% Plot Comparisons
sec = 60;

figure;
title('CZ, Neighbor CAR, CZ LaPlacian (all after filtering)');
hold on;
plot(cz_lap.time{1}(1:sec*cz_lap.fsample),cz_lap.trial{1}(1:sec*cz_lap.fsample)+400,'k');
plot(cz.time{1}(1:sec*cz_lap.fsample),cz.trial{1}(1:sec*cz_lap.fsample)+200,'b');
plot(cz_lap.time{1}(1:sec*cz_lap.fsample),cz_neigh_car(1:sec*cz_lap.fsample),'r');
for ch_ix = 1:numel(cz_neigh.label)
    plot(cz_lap.time{1}(1:sec*cz_lap.fsample),cz_neigh.trial{1}(ch_ix,1:sec*cz_lap.fsample)-200*ch_ix,'g');
end
legend('CZ-LP','CZ','CZ-neigh-CAR');

if strcmp(psd_reref,'yes')
    psd_dir = strcat(SBJ_vars.dirs.preproc,'PSDs/eeg/');
    if ~exist(psd_dir,'dir')
        mkdir(psd_dir);
    end
    fn_plot_PSD_1by1_compare_save(cz.trial{1},cz_lap.trial{1},cz.label,cz_lap.label,...
        cz.fsample,strcat(psd_dir,SBJ,'_cz_PSD_filt_lp'),'filt','filt-lp',psd_fig_type);
end

%% EOG processing
% Filter
eog_orig = eog;
cfg          = [];
cfg.bsfilter = 'yes';
cfg.bsfreq   = bs_freq_lim;
cfg.demean   = proc_vars.demean_yn;
cfg.lpfilter = proc_vars.lp_yn;
cfg.lpfreq   = proc_vars.lp_freq;
cfg.hpfilter = proc_vars.hp_yn;
cfg.hpfreq   = proc_vars.hp_freq;
if isfield(proc_vars,'hp_order')
    cfg.hpfiltord = proc_vars.hp_order;
end
eog = ft_preprocessing(cfg,eog);

% lul = eog_filt.trial{1}(1,:)-eog_filt.trial{1}(2,:);
% rul = eog_filt.trial{1}(3,:)-eog_filt.trial{1}(4,:);
% url = eog_filt.trial{1}(3,:)-eog_filt.trial{1}(1,:);
% lrl = eog_filt.trial{1}(4,:)-eog_filt.trial{1}(2,:);
% eog_filt_ref = eog_filt;
% eog_filt_ref.trial{1} = vertcat(lul,rul,url,lrl);
% eog_filt_ref.label = {'lul','rul','url','lrl'};

% % Plot Comparisons
% figure;
% hold on;
% data = eog_filt_ref;
% title('eog-filt-ref');
% spacer = 200;
% sec = 500;
% for ch_ix = 1:numel(data.label)
%     if ~any(strcmp(data.label{ch_ix},SBJ_vars.ch_lab.eeg_bad))
%         plot(data.time{1}(1:sec*data.fsample),data.trial{1}(ch_ix,1:sec*data.fsample)+spacer*ch_ix);
%     end
%     legend(data.label);
% end
% 
% error('fix to be eog vs eog_orig');
% if strcmp(psd_reref,'yes')
%     psd_dir = strcat(SBJ_vars.dirs.preproc,'PSDs/eeg/');
%     if ~exist(psd_dir,'dir')
%         mkdir(psd_dir);
%     end
%     fn_plot_PSD_1by1_compare_save(eog.trial{1},eog_filt.trial{1},eog.label,eog_filt.label,...
%         eog.fsample,strcat(psd_dir,SBJ,'_eog_PSD_filt'),'eog','eog_filt',psd_fig_type);
% end

%% Save the data
eeg_out_filename = strcat(SBJ_vars.dirs.preproc,SBJ,'_cz_lap.mat');
save(eeg_out_filename, '-v7.3', 'cz');
eog_out_filename = strcat(SBJ_vars.dirs.preproc,SBJ,'_eog_filt.mat');
save(eog_out_filename, '-v7.3', 'eog');
