%% Plot preprocessed data to see if it's really clean
%   1) Load preproc data (downsample, highpass, trim, notch)
%   2) Plot data to identify bad channels and epochs
clear all; %close all;

addpath(genpath('/home/knight/hoycw/PRJ_Stroop/scripts/_TOOLBOXES/'));
addpath('/home/knight/hoycw/Apps/fieldtrip/');
ft_defaults

SBJ = 'IR31';
SBJ_dir = fullfile('/home/knight/hoycw/PRJ_Stroop/data/',SBJ,'/');
preproc_dir = fullfile(SBJ_dir,'02_preproc/');
preclean_filename = strcat(preproc_dir,SBJ,'_ft_preproc.mat');

plot_psd = '1by1';

%% Load the data
fprintf('============== Loading %s ==============\n',SBJ);
load(preclean_filename);
load(strcat(SBJ_dir,'03_events/',SBJ,'_bob_bad_epochs.mat'));

%% Plot data: check good channels, refs, artifacts, etc.
cfg_plot = [];
cfg_plot.viewmode   = 'vertical'; %trace per line
cfg_plot.continuous = 'yes';    % not trials
cfg_plot.plotlabels = 'yes';    %to know what to cut
cfg_plot.blocksize  = 10;        % good for bob, can adjust
cfg_plot.ylim       = [-100 100];
cfg_plot.artfctdef.visual.artifact = bad_epochs;

% need to load these variables
% cfg_plot.event = [word_onset./1000, resp_onset./1000];%
% cfg_plot.plotevents = 'yes';
% cfg_plot.ploteventlabels = 'trial=1';
% plot_ch = data_notch.label;
% cfg_plot.artfctdef.visual.artifact = bad_epochs;

plot_out = ft_databrowser(cfg_plot, data); %pulls up your plotted data

%% 
if strcmp(plot_psd,'all')
    fn_plot_PSD_all(data.trial{1},data.label,data.fsample,1);
elseif strcmp(plot_psd,'1by1')
    fn_plot_PSD_1by1(data.trial{1},data.label,data.fsample);
end

%%
% bad_epochs = plot_out.artfctdef.visual.artifact;
% save(strcat('/home/knight/hoycw/PRJ_Stroop/data/',SBJ,'/02_preproc/',SBJ,'_bad_epochs.mat'),'bad_epochs');

