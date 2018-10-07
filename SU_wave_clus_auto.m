function SU_wave_clus_auto(SBJ)
%% Run wave_clus automated spike sorting
% Paths
if exist('/home/knight/hoycw/','dir');root_dir='/home/knight/hoycw/';app_dir=[root_dir 'Apps/'];
else root_dir='/Volumes/hoycw_clust/';app_dir='/Users/colinhoy/Code/Apps/';end
addpath(genpath([app_dir 'wave_clus/']));

eval(['run ' root_dir 'PRJ_Stroop/scripts/SBJ_vars/' SBJ '_vars.m']);

cd([SBJ_vars.dirs.SU 'micro/']);
Get_spikes('all');
Do_clustering('all');

end