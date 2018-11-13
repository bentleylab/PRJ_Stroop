function SBJ13c_FOOOF_peak_params(SBJ,an_id)
%% Load FOOOF results and plot them

%% Check which root directory
if exist('/home/knight/hoycw/','dir');root_dir='/home/knight/hoycw/';ft_dir=[root_dir 'Apps/fieldtrip/'];
else root_dir='/Volumes/hoycw_clust/';ft_dir='/Users/colinhoy/Code/Apps/fieldtrip/';end

%% Set Up Directories
addpath([root_dir 'PRJ_Stroop/scripts/']);
addpath([root_dir 'PRJ_Stroop/scripts/utils/']);
addpath(ft_dir);
ft_defaults

%% Data Preparation
% Load Data
eval(['run ' root_dir 'PRJ_Stroop/scripts/SBJ_vars/' SBJ '_vars.m']);
eval(['run ' root_dir 'PRJ_Stroop/scripts/an_vars/' an_id '_vars.m']);

% Load ROI and GM/WM info
[roi_list, roi_colors] = fn_roi_label_styles(roi_id);
if strcmp(atlas_id,'Yeo7') || strcmp(atlas_id,'Yeo17')
    elec_space = 'mni_v';
else
    elec_space = 'pat';
end
elec_fname = [SBJ_vars.dirs.recon SBJ '_elec_' pipeline_id '_' elec_space '_' atlas_id '.mat'];
load(elec_fname);

%% Load FOOOF outputs
fooof_dir = [SBJ_vars.dirs.preproc an_id '_results/'];
peak_params = cell(numel(elec.label),1);
bkgd_params = cell(numel(elec.label),1);
for e = 1:numel(elec.label)
    % Loads:
    %   background_params: [offset, slope]
    %   peak_params: [CF, Amp, BW] with 1 row per peak
    %   r_squared: float of model fit
    %   error: float of root mean squared error
    %   gaussian_params: [n_peaks x 3]
    tmp = load([fooof_dir SBJ '_' an_id '_fooof_' elec.label{e}]);
    peak_params{e} = tmp.peak_params;
    bkgd_params{e} = tmp.background_params;
    
    elec.roi{e} = fn_atlas2roi_labels(elec.label(e),atlas_id,roi_id);
end

%% 


end