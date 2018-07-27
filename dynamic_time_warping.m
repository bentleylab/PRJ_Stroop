%% Dynamic time warping testing

%% Set Up
if exist('/home/knight/hoycw/','dir');root_dir='/home/knight/hoycw/';ft_dir=[root_dir 'Apps/fieldtrip/'];
else root_dir='/Volumes/hoycw_clust/';ft_dir='/Users/colinhoy/Code/Apps/fieldtrip/';end

%% Paths
addpath([root_dir 'PRJ_Stroop/scripts/']);
addpath([root_dir 'PRJ_Stroop/scripts/utils/']);
addpath(genpath([ft_dir(1:end-10) 'dPCA/']));

%% Parameters
SBJ = 'IR35';
an_id = 'HGm_S_zbtS_trl2to251_sm0_wn100_stat0';
conditions = 'CNI';

%% Load Data
eval(['run ' root_dir 'PRJ_Stroop/scripts/SBJ_vars/' SBJ '_vars.m']);

% Data
load(strcat(SBJ_vars.dirs.events,SBJ,'_trial_info_final.mat'),'trial_info');
load([SBJ_vars.dirs.proc SBJ '_HFA_ROI_' an_id '.mat']);

% Design matrix
[cond_lab, cond_colors, ~] = fn_condition_label_styles(conditions);
% cond_idx = false([length(cond_lab) length(trial_info.trial_n)]);
design_filename = [SBJ_vars.dirs.proc SBJ '_' conditions '_design_mtx.mat'];
if exist(design_filename)
    load(design_filename);
else
    design = zeros([length(trial_info.trial_n) length(cond_lab)+1]);
    for cond_ix = 1:length(cond_lab)
        % Get binary condition index
        design(:,cond_ix) = fn_condition_index(cond_lab{cond_ix},...
            trial_info.condition_n);
    end
    design(:,numel(cond_lab)+1) = trial_info.response_time;
    save(design_filename,'-v7.3','design');
end
design = logical(design(:,1:numel(cond_lab)));   % lose RTs for now...

%% Dynamic Time Warping
% NEED SOMETHING THAT WORKS WITH OLDER MATLAB VERSIONS! this is 2016a+
distances = zeros([numel(hfa.label) numel(trial_info.trial_n)]);
x_steps   = zeros([numel(hfa.label) numel(trial_info.trial_n) numel(hfa.time)]);
y_steps   = zeros([numel(hfa.label) numel(trial_info.trial_n) numel(hfa.time)]);
for ch_ix = 1
    for t_ix = 2:numel(trial_info.trial_n)
%         [distances(ch_ix,t_ix), x_steps(ch_ix,t_ix,:), y_steps(ch_ix,t_ix,:)] = ...
%             dtw(squeeze(hfa.powspctrm(1,ch_ix,1,:)), squeeze(hfa.powspctrm(t_ix,ch_ix,1,:)));
    end
end
