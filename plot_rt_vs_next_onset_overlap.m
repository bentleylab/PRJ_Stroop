SBJs = {'IR21','IR27','IR31','IR32','IR35','IR39','IR41','IR52','IR54','IR57','IR61'};

for sbj_ix = 1:numel(SBJs)
    SBJ = SBJs{sbj_ix};
    % Load variables
    SBJ_vars_cmd = ['run /home/knight/hoycw/PRJ_Stroop/scripts/SBJ_vars/' SBJ '_vars.m'];
    eval(SBJ_vars_cmd);
    
    % Load trial_info
    load([SBJ_vars.dirs.events SBJ '_trial_info_final.mat']);
    
    % Get difference between response time and next trial onset
    diffs = [];
    for t_ix = 1:numel(trial_info.trial_n)-1    % don't do the last one
        if trial_info.trial_n(t_ix+1)==trial_info.trial_n(t_ix)+1
            diffs = [diffs (trial_info.word_onset(t_ix+1)-trial_info.resp_onset(t_ix))/trial_info.sample_rate];
        end
    end
    figure;
    histogram(diffs,linspace(-0.3,3,20));
    title([SBJ ' RT to next onset histogram']);
end