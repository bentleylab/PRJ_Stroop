restoredefaultpath;
root_dir = '/Users/anaskhan/Desktop/';
ft_dir = [root_dir 'Apps/fieldtrip/'];

addpath(ft_dir);
ft_defaults;

addpath([root_dir 'PRJ_Stroop/scripts']);
addpath([root_dir 'PRJ_Stroop/scripts/utils']);
%%
proc_id = 'main_ft';
SBJ = 'IR82';
SBJ_dir = ['/Users/anaskhan/Desktop/PRJ_Stroop/data/' SBJ '/05_recon/'];
eval(['run ' root_dir 'PRJ_Stroop/scripts/SBJ_vars/' SBJ '_vars.m']);
%%

load(strcat(SBJ_vars.dirs.events,SBJ,'_trial_info_final.mat'));
%%
atlas_id = 'Dx';
load([SBJ_vars.dirs.recon,SBJ,'_elec_',proc_id,'_pat_',atlas_id,'_final.mat']);

% fn_view_recon(SBJ, '', 'ortho', 'pat', '', 1, 'b', 1)
%%
save([SBJ '_elec_main_ft_pat_Dx_final_paper.mat'],'elec')