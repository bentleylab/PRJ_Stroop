function SU_wave_clus_auto(SBJ,path)
%% Run wave_clus automated spike sorting
% Paths
if exist('/home/knight/hoycw/','dir');root_dir='/home/knight/hoycw/';app_dir=[root_dir 'Apps/'];
else root_dir='/Volumes/hoycw_clust/';app_dir='/Users/colinhoy/Code/Apps/';end
addpath(genpath([app_dir 'wave_clus/']));

%% Do clustering
eval(['run ' root_dir 'PRJ_Stroop/scripts/SBJ_vars/' SBJ '_vars.m']);
if isempty(path) % .ncs files (default)
    cd([SBJ_vars.dirs.SU 'micro/']);
    Get_spikes('all');
    Do_clustering('all');
else            % .mat files
    cd(path);
    load('nse_wave_clus_param.mat');
    % Get list of spike files
    mat_fnames = dir([path '/*mat']);
    mat_fnames(strcmp({mat_fnames.name},'nse_wave_clus_param.mat')) = [];
    Get_spikes({mat_fnames.name},'par',par);
    Do_clustering('all');
end

end