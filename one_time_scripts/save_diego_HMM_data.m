%% Export Data for Diego Decoding + HMM

if exist('/home/knight/hoycw/','dir');root_dir='/home/knight/hoycw/';ft_dir=[root_dir 'Apps/fieldtrip/'];
else root_dir='/Volumes/hoycw_clust/';ft_dir='/Users/colinhoy/Code/Apps/fieldtrip/';end
%%
addpath([root_dir 'PRJ_Stroop/scripts/']);
addpath([root_dir 'PRJ_Stroop/scripts/utils/']);
addpath(ft_dir);
ft_defaults

%%
SBJ = 'IR35';
proc_id = 'main_ft';
eval(['run ' root_dir 'PRJ_Stroop/scripts/SBJ_vars/' SBJ '_vars.m']);

%% Load Data
an_id = 'HGm_S2t151_zbtA_sm0_wn100';
load(strcat(SBJ_vars.dirs.events,SBJ,'_trial_info_final.mat'));
load(strcat(SBJ_vars.dirs.proc,SBJ,'_ROI_',an_id,'.mat'));
load(strcat(SBJ_vars.dirs.recon,SBJ,'_elec_main_ft_pat_Dx_final.mat'));

%%
good_chan = {'LAC2-3','LAC5-6','LPC5-6','LOF2-3','RIN3-4','RIN5-6'};
ok_chan = {'LIN4-5','LAC7-8','ROF1-2','LTH6-7','ROF5-6'};

cfg = [];
cfg.channel = good_chan;
good_data = ft_selectdata(cfg,hfa);
good_elec = fn_select_elec(cfg,elec);
cfg.channel = ok_chan;
ok_data = ft_selectdata(cfg,hfa);
ok_elec = fn_select_elec(cfg,elec);

%%
diego_dir = [SBJ_vars.dirs.proc 'diego_example_data/'];
[~] = mkdir(diego_dir);
good_fname = [diego_dir SBJ '_good_hfa.mat'];
save(good_fname,'-v7.3','an_id','good_data','good_elec');
ok_fname = [diego_dir SBJ '_ok_hfa.mat'];
save(ok_fname,'-v7.3','an_id','ok_data','ok_elec');
