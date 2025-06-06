function SBJ10a_corrRT_regressRT_ANOVA(SBJ,an_id,stat_id)
%% Run ANOVA with potential RT regression before
% Set up paths
addpath('/home/knight/hoycw/PRJ_Stroop/scripts/');
addpath('/home/knight/hoycw/PRJ_Stroop/scripts/utils/');
addpath('/home/knight/hoycw/Apps/fieldtrip/');
ft_defaults

rng('shuffle'); % seed randi with time

%% Load Data
eval(['run /home/knight/hoycw/PRJ_Stroop/scripts/SBJ_vars/' SBJ '_vars.m']);
eval(['run /home/knight/hoycw/PRJ_Stroop/scripts/an_vars/' an_id '_vars.m']);
eval(['run /home/knight/hoycw/PRJ_Stroop/scripts/stat_vars/' stat_id '_vars.m']);

load(strcat(SBJ_vars.dirs.events,SBJ,'_trial_info_final.mat'));
load(strcat(SBJ_vars.dirs.proc,SBJ,'_HFA_ROI_',an_id,'.mat'));

% Select data in stat window
cfg_trim = [];
cfg_trim.latency = [stat_lim(1) stat_lim(2)+0.001]; %add a data point to get a full window over the tail end
hfa = ft_selectdata(cfg_trim,hfa);

%% Build Design Matrix
design = zeros([numel(trial_info.trial_n) numel(groups)]);
for grp_ix = 1:numel(groups)
    for level_ix = 1:numel(levels)
        design(:,grp_ix) = design(:,grp_ix) + ...
            level_ix*fn_condition_index(levels{grp_ix}{level_ix}, trial_info.condition_n)';
    end
end

%% Run correlations with RT
if rt_correlation
    cfg_rt.design           = zscore(trial_info.response_time);
    cfg_rt.ivar             = 1;
    stat = ft_freqstatistics(cfg_rt, hfa);
end
% OUTPUT:
%   .rho   = correlation coefficients
%   .stat  = t values
%   .prob  = p values

%% Regress off Reaction Time
if regress_rt
    cfg_conf = [];
    cfg_conf.model    = 'yes';
    cfg_conf.confound = zscore(trial_info.response_time);
    hfa = ft_regressconfound(cfg_conf, hfa);
end
% OUTPUT:
%   .beta  = weights
%   .model = confounds * weights = X * X\Y
%   .powspctrm = Yclean = Y - X * X\Y (i.e., the residuals after subtracting the model off
% NOTE: .stat and .prob output if .statistics = 'yes', but DO NOT TRUST THEM!!! (e.g., p values of 2...)

%% Run ANOVA
% Sliding window parameters
win_lim    = fn_sliding_window_lim(squeeze(hfa.powspctrm(1,1,1,:)),win_len,win_step);
win_center = round(mean(win_lim,2));

% Create structure for w2 in fieldtrip style
w2.cond   = groups;
w2.time   = hfa.time(win_center);
w2.label  = hfa.label;
w2.dimord = 'rpt_chan_time';
w2.trial  = zeros([numel(w2.cond) length(hfa.label) length(w2.time)]);
w2.boot   = zeros([numel(w2.cond) length(hfa.label) length(w2.time) n_boots]);
w2.pval   = w2.trial;

% Compute ANOVA and Explained Variance
for ch_ix = 1:length(hfa.label)
    fprintf('============================ %s (%i / %i) ============================\n',...
        hfa.label{ch_ix},ch_ix,numel(hfa.label));
    fprintf('   time = ');
    for t_ix = 1:numel(win_center);
        fprintf('%.3f, ',hfa.time(win_center(t_ix)));
        % Compute ANOVA (real 1st iteration, then bootstrap with randomized design matrix
        for boot_ix = 1:n_boots+1
            % Randomize design matrix after first real model
            grp_col = {};
            if boot_ix == 1 % Real Model
                for grp_ix = 1:numel(groups)
                    grp_col = {grp_col{:} design(:,grp_ix)};
                end
            else            % Randomized Model
                for grp_ix = 1:numel(groups)
                    grp_col = {grp_col{:} design(randi(size(design,1),size(design(:,grp_ix))),grp_ix)};
                end
            end
            [p table] = anovan(squeeze(nanmean(hfa.powspctrm(:,ch_ix,1,win_lim(t_ix,1):win_lim(t_ix,2)),4)), ...
                grp_col, ...% design(:,2)}, ...
                'model', 'linear', 'sstype', 2, ...% 'continuous', strmatch('RT',w2.cond),
                'varnames', w2.cond, 'display', 'off');
            
            % Calculate w2 (debiased effect size; multiply with 100 to get PEV)
            %   table(:,2) = sum of squares
            %   table(:,3) = dof
            %   table(:,5) = mean square error
            %   table(:,6) = F statistic for main effects, only for rows = groups
            %   table(:,7) = p value for main effects, only for rows = groups
            %   rows: 1 = labels, 2:2+n_cond = groups, end-1 = error, end = total
            %   w2 = (SSbw - df*MSE) / (SStot + MSE)
            mse = table{numel(w2.cond)+2,5};
            for factor_ix = 1:numel(w2.cond)
                factor_row = strmatch(w2.cond{factor_ix},table(:,1));
                if boot_ix == 1
                    w2.trial(factor_ix,ch_ix,t_ix) = (table{factor_row,2} - (table{factor_row,3} * mse))/...
                        (table{end,2} + mse);
                else
                    w2.boot(factor_ix,ch_ix,t_ix,boot_ix-1) = (table{factor_row,2} - (table{factor_row,3} * mse))/...
                        (table{end,2} + mse);
                end
            end
            clear p table
        end
        
        % Compute Significance accross bootstraps
        for factor_ix = 1:numel(w2.cond)
            w2_false_pos = find(w2.boot(factor_ix,ch_ix,t_ix,:) > w2.trial(factor_ix,ch_ix,t_ix));
            w2.pval(factor_ix,ch_ix,t_ix) = numel(w2_false_pos)/n_boots;
        end
    end
    fprintf('\n');
end

%% Save Results
f_name = [SBJ_vars.dirs.proc SBJ '_ANOVA_ROI_' stat_id '_' an_id '.mat'];
if rt_correlation
    save(f_name,'-v7.3','hfa','w2','stat');
else
    save(f_name,'-v7.3','hfa','w2');
end

end
