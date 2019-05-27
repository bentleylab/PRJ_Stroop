if exist('/home/knight/hoycw/','dir');root_dir='/home/knight/hoycw/';ft_dir=[root_dir 'Apps/fieldtrip/'];
else root_dir='/Volumes/hoycw_clust/';ft_dir='/Users/colinhoy/Code/Apps/fieldtrip/';end

%%
SBJ = 'CP24';
eval(['run ' root_dir 'PRJ_Stroop/scripts/SBJ_vars/' SBJ '_vars.m']);

load([SBJ_vars.dirs.import SBJ '_1000hz_R1.mat']);

cfg = [];
cfg.method = 'runica';
ica = ft_componentanalysis(cfg,data); %took ~7-8 min

%%
% cfg = [];
% cfg.layout   = 'biosemi64.lay';
% cfg.channel  = 'all';
% cfg.viewmode = 'component';
load([root_dir 'PRJ_Stroop/scripts/utils/cfg_plot.mat']);
ft_databrowser(cfg_plot, ica);
ft_databrowser(cfg_plot, data);

%%
cfg = [];
cfg.component = [2 3];
cfg.demean = 'no';
clean_data = ft_rejectcomponent(cfg, ica);
% seems to work okay! at the very least decreases amplitude, which will
% help a lot with single trial HFA effects