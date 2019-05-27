SBJs = {'IR21','IR31','IR32','IR35','IR39','IR41','IR54','IR57'};
% max_RTs = zeros(size(SBJs));
% figure;
str = 'not yet';
for s = 1:numel(SBJs)
    tmp = load(['/home/knight/hoycw/PRJ_Stroop/data/' SBJs{s} '/03_events/' SBJs{s} '_trial_info_final.mat']);
%     max_RTs(s) = min(ITI);
    max_RTs(s) = max(tmp.trial_info.response_time);
    fprintf('%s = %f\n',SBJs{s},max_RTs(s));
    % Print RTs below 0.3s
%     rt_over2 = find(tmp.trial_info.response_time<0.3);
%     fprintf('%s = %i\n',SBJs{s},numel(rt_over2));
%     if ~isempty(rt_over2)
%         tmp.trial_info.response_time(rt_over2)
%     end
%     if ~isempty(find(tmp.trial_info.response_time<0))
%         str = 'yes!';
%     else
%         str = 'no';
%     end
%     histogram(trial_info.response_time);
% %     set(gca,'xlim',[0 4]);
%     title(SBJs{s});
%     pause;
    clear tmp
end
% fprintf('Maximum RT = %f\n',max(max_RTs));