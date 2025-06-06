[root_dir,~] = fn_get_root_dir();
SBJs = {'CP24','IR21','IR26','IR31','IR32','IR35','IR39','IR41',...
        'IR52','IR54','IR57','IR61','IR65','IR67','IR68','IR72','IR74'};
% max_RTs = zeros(size(SBJs));
% figure;
min_rt_thresh  = 0.35;
min_rt_thresh2 = 0.3;
str = 'not yet';
for s = 1:numel(SBJs)
    load([root_dir 'PRJ_Stroop/data/' SBJs{s} '/03_events/' SBJs{s} '_trial_info_final.mat']);
%     max_RTs(s) = min(ITI);
%     max_RTs(s) = max(tmp.trial_info.response_time);
    fprintf('%s = %.3f - %.3f\n',SBJs{s},min(trial_info.response_time),max(trial_info.response_time));
    
    % Check that ti.response_time is accurate
    %   IR67 is off by +/-1ms, likely due to sampling rate = 500 Hz
    if any((trial_info.response_time-(trial_info.resp_onset-trial_info.word_onset)./trial_info.sample_rate)>0.0011)
        error(['bad rts in ' SBJs{s}]);
    end
    
    % Print RTs below 0.4s
    short_rts = find(trial_info.response_time<min_rt_thresh);
    if ~isempty(short_rts)
        fprintf(2,'\t%s < %.2f = %i\n',SBJs{s},min_rt_thresh,numel(short_rts));
        if any(trial_info.response_time<min_rt_thresh2)
            fprintf(2,'\t%s < %.2f = %i\n',SBJs{s},min_rt_thresh2,sum(trial_info.response_time<min_rt_thresh2));
        end
%         tmp.trial_info.response_time(short_rts)
    end
%     if ~isempty(find(tmp.trial_info.response_time<0))
%         str = 'yes!';
%     else
%         str = 'no';
%     end
%     histogram(trial_info.response_time);
% %     set(gca,'xlim',[0 4]);
%     title(SBJs{s});
%     pause;
    clear trial_info
end
% fprintf('Maximum RT = %f\n',max(max_RTs));