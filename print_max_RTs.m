SBJs = {'IR21','IR31','IR32','IR35','IR39','IR41','IR54','IR57'};
% max_RTs = zeros(size(SBJs));
% figure;
str = 'not yet';
for s = 1:numel(SBJs)
    tmp = load(['/home/knight/hoycw/PRJ_Stroop/data/' SBJs{s} '/03_events/' SBJs{s} '_trial_info_manual.mat']);
%     ITI = diff(trial_info.word_onset)/1000;
%     max_RTs(s) = min(ITI);
% %     max_RTs(s) = max(trial_info.response_time);
    if ~isempty(find(tmp.trial_info.response_time<0))
        str = 'yes!';
    else
        str = 'no';
    end
    fprintf('%s = %s\n',SBJs{s},str);
%     histogram(trial_info.response_time);
% %     set(gca,'xlim',[0 4]);
%     title(SBJs{s});
%     pause;
    clear tmp
end
% fprintf('Maximum RT = %f\n',max(max_RTs));