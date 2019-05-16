itis = cell(size(SBJs));
for s = 1:numel(SBJs)
    SBJ = SBJs{s};
    load([root_dir 'PRJ_Stroop/data/' SBJ '/03_events/' SBJ '_trial_info_final.mat']);
    block_edges = [diff(trial_info.block_n)];
    trial_skips = [diff(trial_info.trial_n)];
    itis{s} = diff(trial_info.word_onset)/trial_info.sample_rate;
    good_idx = ones(size(itis{s}));
    good_idx(block_edges~=0) = 0;
    good_idx(trial_skips~=1) = 0;
    itis{s}(~good_idx) = [];
    if any(itis{s}>3) || any(itis{s}<2.5)
        fprintf(2,'%s has bad ITIs!\n',SBJ);
    end
    figure;
    histogram(itis{s},100);
    title(SBJ);
    saveas(gcf,[root_dir 'PRJ_Stroop/results/ITI/' SBJ '_ITIs.png']);
    clear trial_info
end

figure;
histogram(vertcat(itis{:}),100);
title('all');
saveas(gcf,[root_dir 'PRJ_Stroop/results/ITI/all_SBJ_ITIs.png']);