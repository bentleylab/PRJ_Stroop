function SBJ10c_HFA_GRP_boxcmp_RT_ANOVA_an_ids(SBJs,stat_id,an_id1,an_id2,save_out)
% Load HFA analysis results for two an_ids to compare RT correlation+ANOVA results
% Analyses compared:
%   an_id1/2: two analysis pipelines (e.g., versions of HFA) to compare
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
rt_lab     = 'RT';
cni_lab    = {'CNI'};

an1_vars_cmd = ['run /home/knight/hoycw/PRJ_Stroop/scripts/an_vars/' an_id1 '_vars.m'];
eval(an1_vars_cmd);
an2_vars_cmd = ['run /home/knight/hoycw/PRJ_Stroop/scripts/an_vars/' an_id2 '_vars.m'];
eval(an2_vars_cmd);
stat_vars_cmd = ['run /home/knight/hoycw/PRJ_Stroop/scripts/stat_vars/' stat_id '_vars.m'];
eval(stat_vars_cmd);

[grp_lab, ~, ~] = fn_group_label_styles(model_lab);
cni_grp_ix = strmatch(cni_lab,grp_lab);

% Set up electrode counts
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
    stat = {}; w2 = {};
    tmp = load(strcat(SBJ_vars.dirs.proc,SBJ,'_ANOVA_ROI_',stat_id,'_',an_id1,'.mat'));
    stat{1} = tmp.stat; w2{1} = tmp.w2;
    tmp = load(strcat(SBJ_vars.dirs.proc,SBJ,'_ANOVA_ROI_',stat_id,'_',an_id2,'.mat'));
    stat{2} = tmp.stat; w2{2} = tmp.w2;
    
    rt_cnt{sbj_ix}  = zeros([2 numel(stat{1}.label)]);
    cni_cnt{sbj_ix} = zeros([2 numel(stat{1}.label)]);
    %% Process parameters
    %!!! is this the best way to do this??? Maybe not...
    sample_rate = (numel(stat{1}.time)-1)/(stat{1}.time(end)-stat{1}.time(1));
    
    % FDR correct pvalues for ANOVA
    qvals = cell([2 1]);
    for an_ix = 1:2
        qvals{an_ix} = NaN(size(w2{an_ix}.pval));
        for ch_ix = 1:numel(stat{an_ix}.label)
            [~, ~, ~, qvals{an_ix}(:,ch_ix,:)] = fdr_bh(squeeze(w2{an_ix}.pval(:,ch_ix,:)));%,0.05,'pdep','yes');
        end
    end
    
    % clip the end of stat_lim
    cfgs = [];
    cfgs.latency = stat_lim;
    stat{1} = ft_selectdata(cfgs,stat{1});
    stat{2} = ft_selectdata(cfgs,stat{2});

    % Confirm channels and time axis are the same
    %   All the rounding and such is because some stupid rounding errors...
    same_start = round(stat{1}.time(1)*uint8(sample_rate))==round(stat{2}.time(1)*uint8(sample_rate));
    same_end   = round(stat{1}.time(end)*uint8(sample_rate))==round(stat{2}.time(end)*uint8(sample_rate));
    same_numel = size(stat{1}.time,2)==size(stat{2}.time,2);
    if ~same_start || ~same_end || ~same_numel
        error('time axes are not the same across hfa analyses!');
    end
    if ~isempty(setdiff(stat{1}.label,stat{2}.label)) || ~isempty(setdiff(stat{2}.label,stat{1}.label))
        error('Different electrodes across hfa analyses!');
    end
    
    %% Aggregate results per ROI
    for ch_ix = 1:numel(stat{1}.label)
        elec_cnt(sbj_ix) = elec_cnt(sbj_ix)+1;
        
        for an_ix = 1:2
            % Check for RT correlations
            if sum(squeeze(stat{an_ix}.mask(ch_ix,1,:)))>0
                rt_cnt{sbj_ix}(an_ix,ch_ix) = 1;
            end
            
            % Check for ANOVA group effects
            if any(squeeze(qvals{an_ix}(cni_grp_ix,ch_ix,:))<0.05)
                cni_cnt{sbj_ix}(an_ix,ch_ix) = 1;
            end
        end
    end
    clear SBJ SBJ_vars w2 stat tmp trial_info qvals
end

%% Print positive and negative findings across analyses
% RT1 vs. RT2
tp = 0; fp = 0; fn = 0; tn = 0;
op = 0; on = 0; np = 0; nn = 0;
for sbj_ix = 1:numel(SBJs)
    tp = tp + sum(rt_cnt{sbj_ix}(1,:)==1 & rt_cnt{sbj_ix}(2,:)==1);
    fp = fp + sum(rt_cnt{sbj_ix}(1,:)==1 & rt_cnt{sbj_ix}(2,:)==0);
    fn = fn + sum(rt_cnt{sbj_ix}(1,:)==0 & rt_cnt{sbj_ix}(2,:)==1);
    tn = tn + sum(rt_cnt{sbj_ix}(1,:)==0 & rt_cnt{sbj_ix}(2,:)==0);
    
    op = op + sum(rt_cnt{sbj_ix}(1,:)==1);
    on = on + sum(rt_cnt{sbj_ix}(1,:)==0);
    np = np + sum(rt_cnt{sbj_ix}(2,:)==1);
    nn = nn + sum(rt_cnt{sbj_ix}(2,:)==0);
end
total = sum(elec_cnt);
fprintf('============== RTs: ==============\n');
fprintf('an_id1: %s\n',an_id1);
fprintf('an_id2: %s\n',an_id2);
fprintf('\t\tid2+\tid2-\tTotal#\n');
fprintf('\tid1+\t%i\t%i\t%i\n',tp,fp,op);
fprintf('\tid1-\t%i\t%i\t%i\n',fn,tn,on);
fprintf('\t\t%i\t%i\t%i\n',np,nn,total);
fprintf('\n\n');

% CNI1 vs. CNI2
tp = 0; fp = 0; fn = 0; tn = 0;
op = 0; on = 0; np = 0; nn = 0;
for sbj_ix = 1:numel(SBJs)
    tp = tp + sum(cni_cnt{sbj_ix}(1,:)==1 & cni_cnt{sbj_ix}(2,:)==1);
    fp = fp + sum(cni_cnt{sbj_ix}(1,:)==1 & cni_cnt{sbj_ix}(2,:)==0);
    fn = fn + sum(cni_cnt{sbj_ix}(1,:)==0 & cni_cnt{sbj_ix}(2,:)==1);
    tn = tn + sum(cni_cnt{sbj_ix}(1,:)==0 & cni_cnt{sbj_ix}(2,:)==0);
    
    op = op + sum(cni_cnt{sbj_ix}(1,:)==1);
    on = on + sum(cni_cnt{sbj_ix}(1,:)==0);
    np = np + sum(cni_cnt{sbj_ix}(2,:)==1);
    nn = nn + sum(cni_cnt{sbj_ix}(2,:)==0);
end
total = sum(elec_cnt);
fprintf('============== CNI ==============\n');
fprintf('an_id1: %s\n',an_id1);
fprintf('an_id2: %s\n',an_id2);
fprintf('\t\tid2+\tid2-\tTotal#\n');
fprintf('\tid1+\t%i\t%i\t%i\n',tp,fp,op);
fprintf('\tid1-\t%i\t%i\t%i\n',fn,tn,on);
fprintf('\t\t%i\t%i\t%i\n',np,nn,total);

%% Save figure
if save_out
    out_dir = ['/home/knight/hoycw/PRJ_Stroop/results/HFA/'];    
    filename = [out_dir 'GRP_HFA_boxcmp_' stat_id '_' an_id1 '_vs_' an_id2 '.mat'];
    fprintf('Saving %s\n',filename);
    save(filename,'cni_cnt','rt_cnt','elec_cnt');
end

end
