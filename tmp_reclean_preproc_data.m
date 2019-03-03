SBJ = 'IR21';
pipeline_id = 'main_ft';

if exist('/home/knight/hoycw/','dir');root_dir='/home/knight/hoycw/';ft_dir=[root_dir 'Apps/fieldtrip/'];
else root_dir='/Volumes/hoycw_clust/';ft_dir='/Users/colinhoy/Code/Apps/fieldtrip/';end

SBJ_vars_cmd = ['run ' root_dir 'PRJ_Stroop/scripts/SBJ_vars/' SBJ '_vars.m'];
eval(SBJ_vars_cmd);

%%
b_ix = 1;
if numel(SBJ_vars.raw_file)==1 || isfield(SBJ_vars.dirs,'nlx')
    block_suffix = '';
else
    block_suffix = strcat('_',SBJ_vars.block_name{b_ix});
end

SBJ00b_view_preclean(SBJ,b_ix,0,'reorder',{},'bad_epochs','load');
% load(strcat(SBJ_vars.dirs.preproc,SBJ,'_preclean',block_suffix,'.mat'));
% raw = data;
% 
% bad_preclean = load(strcat(SBJ_vars.dirs.events,SBJ,'_bob_bad_epochs_preclean',block_suffix,'.mat'));

%%
load(strcat(SBJ_vars.dirs.preproc,SBJ,'_preproc_',pipeline_id,'.mat'));
preclean_ep_at = fn_compile_epochs_full2at(SBJ,pipeline_id);

%% Plot data with bad_epochs highlighted
load(strcat(root_dir,'emodim/scripts/utils/cfg_plot.mat'));

cfgpp = cfg_plot;
cfgpp.artfctdef.visual.artifact = preclean_ep_at;
% if isfield(data,'sampleinfo')
%     data = rmfield(data,'sampleinfo');
% end
out = ft_databrowser(cfgpp,data);

% bad_at = fn_convert_epochs_full2at(bad_preclean.bad_epochs,SBJ_vars.analysis_time{b_ix},...
%                                             strcat(SBJ_vars.dirs.preproc,SBJ,'_preclean',block_suffix,'.mat'),1);

%% View correlations of data
fn_view_correlations(data);
