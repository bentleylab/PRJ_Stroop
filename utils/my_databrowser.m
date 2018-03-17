figure;
hold on;
data = eog_filt_c34;
title('eog-filt-c34');
spacer = 200;
sec = 500;
for ch_ix = 1:numel(data.label)
    if ~any(strcmp(data.label{ch_ix},SBJ_vars.ch_lab.eeg_bad))
        plot(data.time{1}(1:sec*data.fsample),data.trial{1}(ch_ix,1:sec*data.fsample)+spacer*ch_ix);
    end
    legend(data.label);
end