addpath('/home/knight/hoycw/PRJ_Stroop/scripts/');
addpath('/home/knight/hoycw/PRJ_Stroop/scripts/utils/');
addpath('/home/knight/hoycw/Apps/fieldtrip/');
ft_defaults

%%
SBJ = 'IR57';
eval(['run /home/knight/hoycw/PRJ_Stroop/scripts/SBJ_vars/' SBJ '_vars.m']);

cfg = [];
cfg.dataset = SBJ_vars.dirs.raw_filename;
cfg.continuous = 'yes';
cfg.channel = 'all';
data = ft_preprocessing(cfg);

load([SBJ_vars.dirs.events,SBJ,'_bob_bad_epochs_preclean.mat']);
% load([SBJ_vars.dirs.import,SBJ,'_eeg_1000hz.mat']);
load('~/PRJ_Stroop/scripts/utils/cfg_plot.mat');
cfg_plot.artfctdef.visual.artifact = bad_epochs;

% cfg = [];
% both = ft_appenddata(cfg,data,eeg);

%%
out = ft_databrowser(cfg_plot,data);

%%
eeg_bad_epochs = out.artfctdef.visual.artifact;
bad_filename = [SBJ_vars.dirs.events SBJ '_eeg_bad_epochs_preclean.mat'];
save(bad_filename, 'eeg_bad_epochs');