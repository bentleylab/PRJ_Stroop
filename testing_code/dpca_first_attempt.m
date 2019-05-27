%% Demixed PCA attempt 1
% trialNum: N x S x D
% firingRates: N x S x D x T x maxTrialNum
% firingRatesAverage: N x S x D x T
%
% N is the number of neurons
% S is the number of stimuli conditions (F1 frequencies in Romo's task)
% D is the number of decisions (D=2)
% T is the number of time-points (note that all the trials should have the
% same length in time!)
%
% trialNum -- number of trials for each neuron in each S,D condition (is
% usually different for different conditions and different sessions)
%
% firingRates -- all single-trial data together, massive array. Here
% maxTrialNum is the maximum value in trialNum. E.g. if the number of
% trials per condition varied between 1 and 20, then maxTrialNum = 20. For
% the neurons and conditions with less trials, fill remaining entries in
% firingRates with zeros or nans.
%
% firingRatesAverage -- average of firingRates over trials (5th dimension).

%% Set Up
if exist('/home/knight/hoycw/','dir');root_dir='/home/knight/hoycw/';ft_dir=[root_dir 'Apps/fieldtrip/'];
else root_dir='/Volumes/hoycw_clust/';ft_dir='/Users/colinhoy/Code/Apps/fieldtrip/';end
addpath([root_dir 'PRJ_Stroop/scripts/']);
addpath([root_dir 'PRJ_Stroop/scripts/utils/']);
addpath(genpath([ft_dir(1:end-10) 'dPCA/']));

%% Parameters
SBJ = 'IR32';
an_id = 'HGm_R_zbtS_trl5to101_sm0_wn100_stat5to1';%'HGm_S_zbtS_trl2to151_sm0_wn100_stat15';%
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

%% dPCA set up
% Create data matrices
N = numel(hfa.label);  % number of neurons
T = size(hfa.powspctrm,4);  % number of time points
S = 3;                      % number of stimuli
E = max(sum(design,1));     % maximal number of trial repetitions
% Not using D because I don't have errors, color response, etc.
% D = 2;          % number of decisions

data = NaN([N S T E]);
for cond_ix = 1:numel(cond_lab)
    for ch_ix = 1:numel(hfa.label)
        data(ch_ix,cond_ix,:,1:sum(design(:,cond_ix))) = squeeze(hfa.powspctrm(design(:,cond_ix),ch_ix,:))';
    end
end
data_avg = nanmean(data,4);

trial_n = repmat(sum(design,1), [N 1]);

% Define Parameter Groupings
%   parameters:
%       1 - condition
%       2 - time
%       [1 2] - condition/time interaction
%   groupings:
%       {1, [1 2]} - condition grouped with time
%       {2}        - condition independent
param_grps  = {{1, [1 2]}, {2}};
marg_names  = {'Condition', 'Condition-Independent'};
marg_colors = [23 100 171; 187 20 25]/256;  %; 150 150 150; 114 97 171

events = find(hfa.time==0);

% check consistency between trialNum and firingRates
for ch_ix = 1:size(data,1)
    for cond_ix = 1:size(data,2)
        assert(isempty(find(isnan(data(ch_ix,cond_ix,:,1:sum(design(:,cond_ix)))), 1)), 'Something is wrong!')
    end
end

%% Step 1: PCA of the dataset
X = data_avg(:,:);
X = bsxfun(@minus, X, mean(X,2));

[W,~,~] = svd(X, 'econ');
W = W(:,1:20);

% minimal plotting
% dpca_plot(data_avg, W, W, @dpca_plot_default);

% computing explained variance
explVar = dpca_explainedVariance(data_avg, W, W, ...
    'combinedParams', param_grps);

% a bit more informative plotting
dpca_plot(data_avg, W, W, @dpca_plot_default, ...
    'explainedVar', explVar, ...
    'time', hfa.time,                        ...
    'timeEvents', events,               ...
    'marginalizationNames', marg_names, ...
    'marginalizationColours', marg_colors);

%% Step 2: PCA in each marginalization separately

dpca_perMarginalization(data_avg, @dpca_plot_default, ...
   'combinedParams', param_grps);

%% Step 3: dPCA without regularization and ignoring noise covariance

% This is the core function.
% W is the decoder, V is the encoder (ordered by explained variance),
% whichMarg is an array that tells you which component comes from which
% marginalization

tic
[W,V,whichMarg] = dpca(data_avg, 20, ...
    'combinedParams', param_grps);
toc

explVar = dpca_explainedVariance(data_avg, W, V, ...
    'combinedParams', param_grps);

dpca_plot(data_avg, W, V, @dpca_plot_default, ...
    'explainedVar', explVar, ...
    'marginalizationNames', marg_names, ...
    'marginalizationColours', marg_colors, ...
    'whichMarg', whichMarg,                 ...
    'time', hfa.time,                        ...
    'timeEvents', events,               ...
    'legendSubplot', 16);
%     'timeMarginalization', 3, ...

%% Step 4: dPCA with regularization

% This function takes some minutes to run. It will save the computations 
% in a .mat file with a given name. Once computed, you can simply load 
% lambdas out of this file:
%   load('tmp_optimalLambdas.mat', 'optimalLambda')

% Please note that this now includes noise covariance matrix Cnoise which
% tends to provide substantial regularization by itself (even with lambda set
% to zero).

optimalLambda = dpca_optimizeLambda(data_avg, data, trial_n, ...
    'combinedParams', param_grps, ...
    'simultaneous', true, ...
    'numRep', 3, ...  % increase this number to ~10 for better accuracy
    'filename', [SBJ_vars.dirs.proc 'tmp_dpca_' an_id '_optimalLambdas.mat']);

Cnoise = dpca_getNoiseCovariance(data_avg, ...
    data, trial_n, 'simultaneous', true);

[W,V,whichMarg] = dpca(data_avg, 20, ...
    'combinedParams', param_grps, ...
    'lambda', optimalLambda, ...
    'Cnoise', Cnoise);

explVar = dpca_explainedVariance(data_avg, W, V, ...
    'combinedParams', param_grps);

dpca_plot(data_avg, W, V, @dpca_plot_default, ...
    'explainedVar', explVar, ...
    'marginalizationNames', marg_names, ...
    'marginalizationColours', marg_colors, ...
    'whichMarg', whichMarg,                 ...
    'time', hfa.time,                        ...
    'timeEvents', events,               ...
    'legendSubplot', 16);
%     'timeMarginalization', 3,           ...

%% Optional: estimating "signal variance"

% explVar = dpca_explainedVariance(data_avg, W, V, ...
%     'combinedParams', param_grps, ...
%     'Cnoise', Cnoise, 'numOfTrials', trial_n);
% 
% % Note how the pie chart changes relative to the previous figure.
% % That is because it is displaying percentages of (estimated) signal PSTH
% % variances, not total PSTH variances. See paper for more details.
% 
% dpca_plot(firingRatesAverage, W, V, @dpca_plot_default, ...
%     'explainedVar', explVar, ...
%     'marginalizationNames', margNames, ...
%     'marginalizationColours', margColours, ...
%     'whichMarg', whichMarg,                 ...
%     'time', time,                        ...
%     'timeEvents', timeEvents,               ...
%     'timeMarginalization', 3,           ...
%     'legendSubplot', 16);

