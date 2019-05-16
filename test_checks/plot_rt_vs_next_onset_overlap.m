% original script: SBJs = {'IR21','IR27','IR31','IR32','IR35','IR39','IR41','IR52','IR54','IR57','IR61'};
SBJs = {'IR21','IR31','IR32','IR35','IR39','IR41','IR52','IR54','IR57','IR61','IR65','IR68','IR72','IR74'};%ALMOST: 'CP24','IR26',  %NEVER: 'IR27','IR37','IR48',
overlap_thresh = 1.25;

if exist('/home/knight/hoycw/','dir');root_dir='/home/knight/hoycw/';ft_dir=[root_dir 'Apps/fieldtrip/'];
else root_dir='/Volumes/hoycw_clust/';ft_dir='/Users/colinhoy/Code/Apps/fieldtrip/';end

%% Set up paths
addpath([root_dir 'PRJ_Stroop/scripts/']);
addpath([root_dir 'PRJ_Stroop/scripts/utils/']);
addpath(ft_dir);
ft_defaults

%%
for sbj_ix = 1:numel(SBJs)
    SBJ = SBJs{sbj_ix};
    % Load variables
    SBJ_vars_cmd = ['run ' root_dir 'PRJ_Stroop/scripts/SBJ_vars/' SBJ '_vars.m'];
    eval(SBJ_vars_cmd);
    
    % Load trial_info
    load([SBJ_vars.dirs.events SBJ '_trial_info_final.mat']);
    
    % Get difference between response time and next trial onset
    diffs = [];
    for t_ix = 1:numel(trial_info.trial_n)-1    % don't do the last one
        if trial_info.trial_n(t_ix+1)==trial_info.trial_n(t_ix)+1   % don't count skipped trials
            diffs = [diffs (trial_info.word_onset(t_ix+1)-trial_info.resp_onset(t_ix))/trial_info.sample_rate];
        end
    end
%     if any(diffs<overlap_thresh) % given overlap_thresh post-R analysis period, how many are bad?
        fprintf('%s trials <1s post-R:\t%i / %i\n',SBJ,sum(diffs<overlap_thresh),numel(trial_info.trial_n));
        figure;
        histogram(diffs,linspace(-0.3,3,20));
        title([SBJ ' RT to next onset histogram']);
%     end
end