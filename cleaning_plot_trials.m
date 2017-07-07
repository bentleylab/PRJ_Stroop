%% Segment and plot the preclean (preproc) data
addpath('/home/knight/hoycw/Apps/fieldtrip/'); % read .besa
ft_defaults
clear all
%%
SBJ = 'IR54';
PRJ = 'Stroop';
SBJ_dir = strcat('/home/knight/hoycw/PRJ_',PRJ,'/data/',SBJ,'/');
preclean_filename = fullfile(SBJ_dir,'02_preproc',[SBJ '_preclean_B1_ft.mat']);

%% Load Data
load(preclean_filename);
% cfg = [];
% cfg.dataset = preclean_filename;
% cfg.continuous = 'yes';
% cfg.channel = 'all';
% data_preclean = ft_preprocessing(cfg);   % just load the data, don't process it


%% Cut into trials
% cfg = [];
% % cfg.dataset = strcat(SBJ_dir,'04_proc/',data_id,'_data_only.mat');
% % Define trials as my whole period of interest (buffer to buffer, no offset)
% cfg.trl = [word_onset-buff_lim(1), ...        % trial onset
%     resp_onset+buff_lim(2),...    % trial offset
%     repmat(0,length(ok_epochs),1),...         %baseline before each trial
%     cond_n];                 % trial type
% % cfg.continuous = 'yes';
% trials = ft_redefinetrial(cfg, data_preclean);


%% PLot
cfg_plot = [];
cfg_plot.viewmode = 'vertical'; %trace per line
cfg_plot.continuous = 'yes';    % not trials
cfg_plot.plotlabels = 'yes';    %to know what to cut
cfg_plot.blocksize = 10;        % good for bob, can adjust
plot_out = ft_databrowser(cfg_plot, data_notch); %pulls up your plotted data
