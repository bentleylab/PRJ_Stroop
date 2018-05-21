%% Testing CAR methods in fieldtrip
% if all channels are included
raw = probe_data.trial{1};
avg = mean(raw,1); % assumes no ref_excluded
for row_ix = 1:size(raw,1); reref(row_ix,:) = raw(row_ix,:)-avg;end

% Compare to un-re-referenced signal
figure;
for row_ix = 1:size(raw,1);
    plot(avg(1:500),'k--');
    hold on;
    plot(raw(row_ix,1:500),'b');
    plot(reref(row_ix,1:500),'r');
    pause;hold off;
end

% Compare to fieldtrip
figure;
ts_len = 1000;
for row_ix = 1:size(reref,1);
    plot(data_reref{d}.trial{1}(row_ix,1:ts_len),'k');
    hold on;
%     plot(raw(row_ix,1:500),'b');
    plot(reref(row_ix,1:ts_len),'r');
    pause;hold off;
end


%% If some channels should be excluded...

% ch1 = probe_data.trial{1}(1,:);
total_signal = zeros(size(probe_data.trial{1}(1,:)));
n_signals = 0;
for ch = 1:size(probe_data.trial{1},1)
    if isempty(strmatch(probe_data.label{ch},SBJ_vars.ref_exclude))
        total_signal = total_signal + probe_data.trial{1}(ch,:);
        n_signals = n_signals+1;
    end
end
CAR = total_signal/n_signals;
reref = probe_data.trial{1};
for ch = 1:size(probe_data.trial{1},1)
    reref(ch,:) = reref(ch,:)-CAR;
end

figure;
ts_len = 1000;
for row_ix = 1:size(reref,1);
    plot(data_reref{d}.trial{1}(row_ix,1:ts_len),'k');
    hold on;
%     plot(raw(row_ix,1:500),'b');
    plot(reref(row_ix,1:ts_len),'r');
    pause;hold off;
end
