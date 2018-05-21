%% Compile Task Design Matrix (with RT)
clear all
SBJ = 'IR54';
pipeline_id = 'main_ft';
an_id = 'HGm_S_zbtS_trl2to15_sm10_wn30_stat15';
% ch_ix = 

% Load Data
SBJ_vars_cmd = ['run /home/knight/hoycw/PRJ_Stroop/scripts/SBJ_vars/' SBJ '_vars.m'];
eval(SBJ_vars_cmd);
an_vars_cmd = ['run /home/knight/hoycw/PRJ_Stroop/scripts/an_vars/' an_id '_vars.m'];
eval(an_vars_cmd);
% eval(['run /home/knight/hoycw/PRJ_Stroop/scripts/proc_vars/' pipeline_id '_proc_vars.m']);

load(strcat(SBJ_vars.dirs.events,SBJ,'_trial_info_final.mat'));
load(strcat(SBJ_vars.dirs.proc,SBJ,'_HFA_ROI_',an_id,'.mat'));

%% Build Regressors
% [cond_lab_CNI, ~, ~] = fn_condition_label_styles('CNI');
% [cond_lab_pcon, ~, ~] = fn_condition_label_styles('pcon');
% cond_lab = horzcat(cond_lab_CNI,cond_lab_pcon,'RT');
% 
% design = zeros([numel(trial_info.trial_n) numel(cond_lab)]);    %+1 not needed b/c fitglm adds it
% % Trial Type
% design(:,1) = fn_condition_index('con', trial_info.condition_n);
% design(:,2) = fn_condition_index('neu', trial_info.condition_n);
% design(:,3) = fn_condition_index('inc', trial_info.condition_n);
% % Block Type
% design(:,4) = fn_condition_index('mcon', trial_info.condition_n);
% design(:,5) = fn_condition_index('same', trial_info.condition_n);
% design(:,6) = fn_condition_index('minc', trial_info.condition_n);
% % Reaction Time (Confound)
% design(:,7) = trial_info.response_time;
% % Intercept - automatically added by fitglm
% % design(:,8) = ones([1 numel(trial_info.trial_n)]);
% 

[cond_lab, ~, ~] = fn_condition_label_styles('CI');
cond_lab = horzcat(cond_lab,'RT');

design = zeros([numel(trial_info.trial_n) numel(cond_lab)]);    %+1 not needed b/c fitglm adds it
% Trial Type
design(:,1) = fn_condition_index('con', trial_info.condition_n);
% design(:,2) = fn_condition_index('neu', trial_info.condition_n);
design(:,2) = fn_condition_index('inc', trial_info.condition_n);
% Reaction Time (Confound)
design(:,3) = zscore(trial_info.response_time);

%% Run ANOVA
win_len    = 200;
win_step   = 50;
% Select data in stat window
cfgs = [];
cfgs.latency = [stat_lim(1) stat_lim(2)+win_step/1000]; %add a data point to get a full window over the tail end
% cfgs.channel = 'RAC6-7';
hfa = ft_selectdata(cfgs,hfa);

win_lim    = fn_sliding_window_lim(squeeze(hfa.powspctrm(1,1,1,:)),win_len,win_step);
win_center = round(mean(win_lim,2));

% Compute ANOVA and Explained Variance
for ch_ix = 1:length(hfa.label)
    for win_ix = 1:numel(win_center);
        % Run GLM
        model = fitglm(design, squeeze(nanmean(hfa.powspctrm(:,ch_ix,1,win_lim(win_ix,1):win_lim(win_ix,2)),4)), ...
            'linear', 'CategoricalVars',[1:numel(cond_lab)-1], 'VarNames',{cond_lab{:} 'HFA'});
        
        
        % Calculate w2 (debiased effect size; multiply with 100 to get PEV)
        %   w2 = (SSbw - df*MSE) / (SStot + MSE)
        %   table(:,2) = sum of squares
        %   table(:,3) = dof
        %   table(:,5) = mean square error
        %   table(:,6) = F statistic for main effects, only for rows = factors
        %   table(:,7) = p value for main effects, only for rows = factors
        %   rows: 1 = labels, 2:2+n_cond = factors, end-1 = error, end = total
%         mse = table{numel(w2.cond)+2,5};
%         for factor_ix = 1:numel(w2.cond)
%             factor_row = strmatch(w2.cond{factor_ix},table(:,1));
%             w2.trial(factor_ix,ch_ix,win_ix) = (table{factor_row,2} - (table{factor_row,3} * mse))/...
%                                             (table{end,2} + mse);
%         end
%         w2.trial(2,ch_ix,t) = (table{3,2}-(table{3,3}*table{numel(w2.cond)+2,5}))/...
%                                 (table{end,2}+table{numel(w2.cond)+2,5});
%         w2.trial(3,ch_ix,t) = (table{4,2}-(table{4,3}*table{numel(w2.cond)+2,5}))/...
%                                 (table{end,2}+table{numel(w2.cond)+2,5});
        
    end
end
