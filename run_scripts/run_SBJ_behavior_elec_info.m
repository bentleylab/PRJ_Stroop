%% Run scripts to get data for Thesis Table 1
if exist('/home/knight/hoycw/','dir');root_dir='/home/knight/hoycw/';ft_dir=[root_dir 'Apps/fieldtrip/'];
else root_dir='/Volumes/hoycw_clust/'; ft_dir='/Users/colinhoy/Code/Apps/fieldtrip/';end

%%
addpath([root_dir 'PRJ_Stroop/scripts/']);
addpath([root_dir 'PRJ_Stroop/scripts/utils/']);
addpath(ft_dir);
ft_defaults

%%
SBJs = {'CP24','IR21','IR26','IR31','IR32','IR35','IR39','IR41',...
        'IR52','IR54','IR57','IR61','IR65','IR67','IR68','IR72','IR74'};

proc_id  = 'main_ft';
roi_id   = 'gROI';
atlas_id = 'Dx';
[roi_list, ~] = fn_roi_label_styles(roi_id);
if any(strcmp(roi_id,{'mgROI','gROI','main3','lat','deep','gPFC'}))
    roi_field = 'gROI';
else
    roi_field = 'ROI';
end

%% Get Behavioral and Electrode Information
roi_cnt = zeros([numel(roi_list) numel(SBJs)]);
n_trials = zeros(size(SBJs));
n_var_rej = nan(size(SBJs));
for s = 1:numel(SBJs)
    SBJ = SBJs{s};
    fprintf('%s:\n',SBJ);
    eval(['run ' root_dir 'PRJ_Stroop/scripts/SBJ_vars/' SBJ '_vars.m']);
    load(strcat(SBJ_vars.dirs.events,SBJ,'_trial_info_final.mat'),'trial_info');
    elec_fname = [SBJ_vars.dirs.recon SBJ '_elec_main_ft_pat_Dx_final.mat'];
    load(elec_fname);
    
    fprintf('\tsrate = %.1f\n',trial_info.sample_rate);
    
    % Number of trials
    n_trials(s) = numel(trial_info.trial_n);
    fprintf('\tn_trials = %d\n',n_trials(s));
    
    % RT
    fprintf('\tmean_rt = %.3f\n',mean(trial_info.response_time));
    
    % Errors
    fprintf('\tn_errors = %d\n',numel(trial_info.bad_trials.error));
    
    % Variance rejection
    n_var_rej(s) = numel(trial_info.bad_trials.variance);
    fprintf('\tn_var_rej = %d\n',numel(trial_info.bad_trials.variance));
    
    % Electrodes
    for roi_ix = 1:numel(roi_list)
        roi_cnt(roi_ix,s) = sum(strcmp(elec.(roi_field),roi_list{roi_ix}));
        fprintf('\t%s count = %d\n',roi_list{roi_ix},roi_cnt(roi_ix,s));
    end
    
    clear SBJ SBJ_vars trial_info elec
end

fprintf('n_trials mean +/- SD = %.1f +/- %.1f\n',mean(n_trials),std(n_trials));
fprintf('n_trials min = %d; max = %d\n',min(n_trials),max(n_trials));
fprintf('mean n_var_rej = %.2f\n',mean(n_var_rej));
fprintf('std n_var_rej = %.2f\n',std(n_var_rej));