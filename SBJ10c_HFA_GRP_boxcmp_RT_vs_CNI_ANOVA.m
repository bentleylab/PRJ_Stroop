function SBJ10c_HFA_GRP_boxcmp_RT_vs_CNI_ANOVA(SBJs,stat_id,ci_an_id,an_id,save_out)
% Load HFA analysis results for original CI HFA test to compare to RT correlation+ANOVA results
% Original analysis:
%   'CI'- original clusterbased permutation corrected via FT results
% Compared to new analyses:
%   RT correlation: correlation with RT and cluster-based perm. corrected via FT
%   CNI ANOVA factor: FDR corrected sliding window ANOVA after regressing off RT as a confound
% OUTPUTS:
%   number of electrodes sig vs. non-sig across the two analyses,
%   separately for orig vs. RT and orig vs. CNI
% clear all; %close all;
if ischar(save_out); save_out = str2num(save_out); end

%% Data Preparation
% Set up paths
addpath('/home/knight/hoycw/PRJ_Stroop/scripts/');
addpath('/home/knight/hoycw/PRJ_Stroop/scripts/utils/');
addpath('/home/knight/hoycw/Apps/fieldtrip/');
ft_defaults

%% Prep variables
conditions = 'CI';
rt_lab     = 'RT';
cni_lab    = {'CNI'};

an_vars_cmd = ['run /home/knight/hoycw/PRJ_Stroop/scripts/an_vars/' an_id '_vars.m'];
eval(an_vars_cmd);
stat_vars_cmd = ['run /home/knight/hoycw/PRJ_Stroop/scripts/stat_vars/' stat_id '_vars.m'];
eval(stat_vars_cmd);

[grp_lab, ~, ~] = fn_group_label_styles(model_lab);
cni_grp_ix = strmatch(cni_lab,grp_lab);

% Set up electrode counts
ci_cnt    = cell(size(SBJs));
rt_cnt    = cell(size(SBJs));
cni_cnt   = cell(size(SBJs));
elec_cnt  = zeros([1 numel(SBJs)]);

%% Load Results
for sbj_ix = 1:numel(SBJs)
    SBJ = SBJs{sbj_ix};
    fprintf('================= Processing: %s =================\n',SBJ);
    % Load variables
    SBJ_vars_cmd = ['run /home/knight/hoycw/PRJ_Stroop/scripts/SBJ_vars/' SBJ '_vars.m'];
    eval(SBJ_vars_cmd);
    
    % Load data
    tmp = load(strcat(SBJ_vars.dirs.proc,SBJ,'_',conditions,'_ROI_',ci_an_id,'.mat'),'stat');
    ci_stat = tmp.stat;
    load(strcat(SBJ_vars.dirs.proc,SBJ,'_ANOVA_ROI_',stat_id,'_',an_id,'.mat'),'stat','w2');
    
    ci_cnt{sbj_ix}  = zeros([1 numel(stat.label)]);
    rt_cnt{sbj_ix}  = zeros([1 numel(stat.label)]);
    cni_cnt{sbj_ix} = zeros([1 numel(stat.label)]);
    %% Process parameters
    %!!! is this the best way to do this??? Maybe not...
    sample_rate = (numel(stat.time)-1)/(stat.time(end)-stat.time(1));
    
    % FDR correct pvalues for ANOVA
    qvals = NaN(size(w2.pval));
    for ch_ix = 1:numel(stat.label)
        [~, ~, ~, qvals(:,ch_ix,:)] = fdr_bh(squeeze(w2.pval(:,ch_ix,:)));%,0.05,'pdep','yes');
    end
    
    % clip the end of stat_lim
    cfgs = [];
    cfgs.latency = stat_lim;
    stat = ft_selectdata(cfgs,stat);
    
    % Confirm channels and time axis are the same
    %   All the rounding and such is because some stupid rounding errors...
    same_start = round(ci_stat.time(1)*uint8(sample_rate))==round(stat.time(1)*uint8(sample_rate));
    same_end   = round(ci_stat.time(end)*uint8(sample_rate))==round(stat.time(end)*uint8(sample_rate));
    same_numel = size(ci_stat.time,2)==size(stat.time,2);
    if ~same_start || ~same_end || ~same_numel
        error('time axes are not the same across hfa analyses!');
    end
    if ~isempty(setdiff(stat.label,ci_stat.label)) || ~isempty(setdiff(ci_stat.label,stat.label))
        error('Different electrodes across hfa analyses!');
    end
    
    %% Aggregate results per ROI
    for ch_ix = 1:numel(stat.label)
        elec_cnt(sbj_ix) = elec_cnt(sbj_ix)+1;
        
        % Check for condition differences, get epochs
        if sum(squeeze(ci_stat.mask(ch_ix,1,:)))>0
            ci_cnt{sbj_ix}(ch_ix) = 1;
        end
        
        % Check for RT correlations, get epochs
        if sum(squeeze(stat.mask(ch_ix,1,:)))>0
            rt_cnt{sbj_ix}(ch_ix) = 1;
        end
        
        % Check for ANOVA group effects
        if any(squeeze(qvals(cni_grp_ix,ch_ix,:))<0.05)
            cni_cnt{sbj_ix}(ch_ix) = 1;
        end
    end
    clear SBJ SBJ_vars w2 ci_stat stat tmp trial_info qvals
end

%% Print positive and negative findings across analyses
% CI vs. RT
tp = 0; fp = 0; fn = 0; tn = 0;
op = 0; on = 0; np = 0; nn = 0;
for sbj_ix = 1:numel(SBJs)
    tp = tp + sum(ci_cnt{sbj_ix}==1 & rt_cnt{sbj_ix}==1);
    fp = fp + sum(ci_cnt{sbj_ix}==1 & rt_cnt{sbj_ix}==0);
    fn = fn + sum(ci_cnt{sbj_ix}==0 & rt_cnt{sbj_ix}==1);
    tn = tn + sum(ci_cnt{sbj_ix}==0 & rt_cnt{sbj_ix}==0);
    
    op = op + sum(ci_cnt{sbj_ix}==1);
    on = on + sum(ci_cnt{sbj_ix}==0);
    np = np + sum(rt_cnt{sbj_ix}==1);
    nn = nn + sum(rt_cnt{sbj_ix}==0);
end
total = sum(elec_cnt);
fprintf('============== CI vs. RT ==============\n');
fprintf('\t\tRT+\tRT-\tTotal#\n');
fprintf('\tCI+\t%i\t%i\t%i\n',tp,fp,op);
fprintf('\tCI-\t%i\t%i\t%i\n',fn,tn,on);
fprintf('\t\t%i\t%i\t%i\n',np,nn,total);
fprintf('\n\n');

% CI vs CNI
tp = 0; fp = 0; fn = 0; tn = 0;
op = 0; on = 0; np = 0; nn = 0;
for sbj_ix = 1:numel(SBJs)
    tp = tp + sum(ci_cnt{sbj_ix}==1 & cni_cnt{sbj_ix}==1);
    fp = fp + sum(ci_cnt{sbj_ix}==1 & cni_cnt{sbj_ix}==0);
    fn = fn + sum(ci_cnt{sbj_ix}==0 & cni_cnt{sbj_ix}==1);
    tn = tn + sum(ci_cnt{sbj_ix}==0 & cni_cnt{sbj_ix}==0);
    
    op = op + sum(ci_cnt{sbj_ix}==1);
    on = on + sum(ci_cnt{sbj_ix}==0);
    np = np + sum(cni_cnt{sbj_ix}==1);
    nn = nn + sum(cni_cnt{sbj_ix}==0);
end
total = sum(elec_cnt);
fprintf('============== CI vs. CNI ==============\n');
fprintf('\t\tCNI+\tCNI-\tTotal#\n');
fprintf('\tCI+\t%i\t%i\t%i\n',tp,fp,op);
fprintf('\tCI-\t%i\t%i\t%i\n',fn,tn,on);
fprintf('\t\t%i\t%i\t%i\n',np,nn,total);

% CNI vs RT
tp = 0; fp = 0; fn = 0; tn = 0;
op = 0; on = 0; np = 0; nn = 0;
for sbj_ix = 1:numel(SBJs)
    tp = tp + sum(rt_cnt{sbj_ix}==1 & cni_cnt{sbj_ix}==1);
    fp = fp + sum(rt_cnt{sbj_ix}==1 & cni_cnt{sbj_ix}==0);
    fn = fn + sum(rt_cnt{sbj_ix}==0 & cni_cnt{sbj_ix}==1);
    tn = tn + sum(rt_cnt{sbj_ix}==0 & cni_cnt{sbj_ix}==0);
    
    op = op + sum(rt_cnt{sbj_ix}==1);
    on = on + sum(rt_cnt{sbj_ix}==0);
    np = np + sum(cni_cnt{sbj_ix}==1);
    nn = nn + sum(cni_cnt{sbj_ix}==0);
end
total = sum(elec_cnt);
fprintf('============== RT vs. CNI ==============\n');
fprintf('\t\tCNI+\tCNI-\tTotal#\n');
fprintf('\tRT+\t%i\t%i\t%i\n',tp,fp,op);
fprintf('\tRT-\t%i\t%i\t%i\n',fn,tn,on);
fprintf('\t\t%i\t%i\t%i\n',np,nn,total);

%% Save figure
if save_out
    out_dir = ['/home/knight/hoycw/PRJ_Stroop/results/HFA/'];    
    filename = [out_dir 'GRP_HFA_boxcmp_CI_vs_' stat_id '_' an_id '.mat'];
    fprintf('Saving %s\n',filename);
    save(filename,'ci_cnt','cni_cnt','rt_cnt','elec_cnt');
end

end
